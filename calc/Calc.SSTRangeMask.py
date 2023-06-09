"""
Calc.SSTRangeMask.py
==========================
Calculate and output the range mask that can be used in pacemaker simulation.
"""

import numpy as np
import xarray as xr
import scipy.interpolate as sciinterp
import skimage.filters as skifilters
import matplotlib.pyplot as plt
import matplotlib as mpl
import cartopy.crs as ccrs
import pop_tools as popt
import os
import sys
sys.path.append('../')
import tools.tool_GetClimData as tool_GetCD
import tools.tool_ClimData_Preprocessing as tool_CP

def Create_SSTRangeMask_IO(Lat, Lon, LandFraction):

	"""
	Create the range mask for Indian Ocean (IO) region
	==========================
	Input:

		Lat (numpy array): latitude

		Lon (numpy array): longitude

		LandFraction (numpy array): land fraction

	Output:

		SSTRangeMask (numpy array): the range mask
	==========================
	"""

	if (LandFraction.shape != (Lat.size, Lon.size)):

		raise ValueError('The shape of LandFraction must equal to (latitude, longitude)')
	
	# Fill in a rectangle region
	SSTRangeMask = np.where((Lat[:, None]>=-20)&(Lat[:, None]<=20)&(Lon[None, :]>=30)&(Lon[None, :]<=120), True, False)

	# Exclude the north-east corner of array
	SSTRangeMask = np.where((Lat[:, None]>=-8)&(Lon[None, :]>=105), False, SSTRangeMask)
	SSTRangeMask = np.where((Lat[:, None]>=3)&(Lon[None, :]>=98), False, SSTRangeMask)
	SSTRangeMask = np.where((Lat[:, None]>=-1)&(Lon[None, :]>=101), False, SSTRangeMask)

	# Land mask
	LandMask   = tool_CP.Get_LandMask(LandFraction, MaskType='Ocean')
	SSTRangeMask = np.where(LandMask, SSTRangeMask, False)
	
	return SSTRangeMask

def Create_SSTRangeMask_WP(Lat, Lon, LandFraction):

	"""
	Create the range mask for Western Pacific (WP) region
	==========================
	Input:
		
		Lat (numpy array): latitude
		
		Lon (numpy array): longitude

		LandFraction (numpy array): land fraction

	Output:

		SSTRangeMask (numpy array): the range mask
	==========================
	"""

	if (LandFraction.shape != (Lat.size, Lon.size)):

		raise ValueError('The shape of LandFraction must equal to (latitude, longitude)')
	
	# Fill in a rectangle region
	SSTRangeMask = np.where((Lat[:, None]>=-20)&(Lat[:, None]<=20)&(Lon[None, :]>=100)&(Lon[None, :]<=170), True, False)

	# Exclude the north-east corner of array
	SSTRangeMask = np.where((Lat[:, None]<=-10)&(Lon[None, :]<=120), False, SSTRangeMask)
	SSTRangeMask = np.where((Lat[:, None]<=-8)&(Lon[None, :]<=114), False, SSTRangeMask)
	SSTRangeMask = np.where((Lat[:, None]<=-6)&(Lon[None, :]<=108), False, SSTRangeMask)
	SSTRangeMask = np.where((Lat[:, None]<=-5)&(Lon[None, :]<=104), False, SSTRangeMask)
	SSTRangeMask = np.where((Lat[:, None]<=-2)&(Lon[None, :]<=103), False, SSTRangeMask)

	# Land mask
	LandMask   = tool_CP.Get_LandMask(LandFraction, MaskType='Ocean')
	SSTRangeMask = np.where(LandMask, SSTRangeMask, False)
	
	return SSTRangeMask

def Calc_Data_Smoothing(Data, Gaussian_sigma=1):

	Data = skifilters.gaussian(Data, sigma=Gaussian_sigma, mode='wrap', preserve_range=True)

	return Data

def Output_SSTRangeMask(Range, SSTRangeMask, File_Path='../src/SSTRangeMask/', Bool_Output_POP=True):

	"""
	Output the range mask.
	==========================
	Input:

		Range (str): the name of the range. This would be contained in the output file name

		SSTRangeMask (numpy array): range mask

		File_Path (str): optional. The path for the file to be output

	Output:

		None
	==========================
	"""

	# Set output grid: 0.9x1.25
	LatLon_Grid     = tool_GetCD.Get_Grid('0.9x1.25', Bool_Meshgrid=False)
	LatLon_Grid_Lat = LatLon_Grid['Lat']
	LatLon_Grid_Lon = LatLon_Grid['Lon']

	# Check whether the file path has existed
	if (not os.path.exists(File_Path)): os.makedirs(File_Path)

	# Write to nc file
	xr.Dataset( \
		data_vars={\
			'SSTRangeMask': xr.DataArray(data=SSTRangeMask, dims=['Lat', 'Lon'], attrs={'description': 'range mask (1=True, 0=False)'}), \
		}, \
		coords={\
			'Lat': LatLon_Grid_Lat, \
			'Lon': LatLon_Grid_Lon, \
		}, \
		attrs={\
			'title': 'Range Mask: {Range}'.format(Range=Range), \
		} \
	).to_netcdf('{File_Path}/SSTRangeMask_{Range}.nc'.format(File_Path=File_Path, Range=Range))

	# Output to POP
	if (Bool_Output_POP):

		# Get pop grid
		POP_Grid     = popt.get_grid('POP_gx1v7')
		POP_Grid_Lat = POP_Grid.TLAT.values
		POP_Grid_Lon = POP_Grid.TLONG.values

		# Get 0.9x1.25 grid
		LatLon_Grid = tool_GetCD.Get_Grid('0.9x1.25')
		LatLon_Grid_Lat = LatLon_Grid['Lat']
		LatLon_Grid_Lon = LatLon_Grid['Lon']

		# Interpolation
		Interp_Linear       = sciinterp.LinearNDInterpolator((LatLon_Grid_Lat.flatten(), LatLon_Grid_Lon.flatten()), SSTRangeMask.flatten())
		Interp_Nearest      = sciinterp.NearestNDInterpolator((LatLon_Grid_Lat.flatten(), LatLon_Grid_Lon.flatten()), SSTRangeMask.flatten())
		Data_Interp_Linear  = Interp_Linear(np.concatenate((POP_Grid_Lat.flatten()[:, None], POP_Grid_Lon.flatten()[:, None]), axis=1))
		Data_Interp_Nearest = Interp_Nearest(np.concatenate((POP_Grid_Lat.flatten()[:, None], POP_Grid_Lon.flatten()[:, None]), axis=1))
		SSTRangeMask_POP    = np.where(np.isnan(Data_Interp_Linear), Data_Interp_Nearest, Data_Interp_Linear).reshape(POP_Grid_Lat.shape)

		# Write to nc file
		xr.Dataset( \
			data_vars={\
				'SSTRangeMask': xr.DataArray(data=SSTRangeMask_POP, coords={'TLAT': (['y', 'x'], POP_Grid_Lat), 'TLONG': (['y', 'x'], POP_Grid_Lon)}, dims=['y', 'x']), \
			}, \
			attrs={\
				'title': 'Range Mask: {Range}'.format(Range=Range), \
			} \
		).to_netcdf('{File_Path}/SSTRangeMask_{Range}.gx1v7.nc'.format(File_Path=File_Path, Range=Range))

	return

def Plot_SSTRangeMask(Plot_Data, Plot_Config):

	# Create figure object
	fig, ax = plt.subplots(figsize=(5, 3), dpi=300, subplot_kw={'projection': ccrs.PlateCarree()})

	# Plot
	Pcolormesh_1 = ax.pcolormesh(\
		Plot_Data['Lon'], Plot_Data['Lat'], Plot_Data['SSTRangeMask'], \
		cmap=mpl.colormaps['Greens'], vmin=0, vmax=1, \
		transform=ccrs.PlateCarree(), \
	)
	
	# Configuration
	ax.set_title('Map: SST Region Mask ({Range})'.format(Range=Plot_Config['Plot_Range']))
	gl = ax.gridlines(draw_labels=True)
	gl.ylocator = mpl.ticker.FixedLocator([-80, -60, -40, -20, 0, 20, 40, 60, 80])
	gl.top_labels = False
	gl.left_labels = False
	ax.coastlines(resolution='110m')

	# Save figure
	Figure_Path = '../output/Output_Figure/Plot.Map.SSTRegionMask/'
	Figure_Name = 'Plot.Map.SSTRegionMask.{Range}.png'.format(Range=Plot_Config['Plot_Range'])
	if not (os.path.exists(Figure_Path)): os.makedirs(Figure_Path)
	plt.savefig(Figure_Path + Figure_Name)

	return

if (__name__ == '__main__'):

	"""
	When this file is called, calculate the range masks over IO and WP regions and write to new nc files.
	"""

	RefInfo = tool_CP.Get_RefInfo()

	# Calculate range mask for IO, WP
	SSTRangeMask_IO = Create_SSTRangeMask_IO(RefInfo['Lat'], RefInfo['Lon'], RefInfo['LandFraction'])
	SSTRangeMask_WP = Create_SSTRangeMask_WP(RefInfo['Lat'], RefInfo['Lon'], RefInfo['LandFraction'])
	SSTRangeMask_IO = np.where(SSTRangeMask_IO, 1, 0)
	SSTRangeMask_WP = np.where(SSTRangeMask_WP, 1, 0)

	# Smooth the range mask
	SSTRangeMask_IO = Calc_Data_Smoothing(SSTRangeMask_IO, 1.5)
	SSTRangeMask_WP = Calc_Data_Smoothing(SSTRangeMask_WP, 1.5)

	# Output to nc file
	Output_SSTRangeMask('IO', SSTRangeMask_IO)
	Output_SSTRangeMask('WP', SSTRangeMask_WP)
	Output_SSTRangeMask('IOWP', np.maximum(SSTRangeMask_IO, SSTRangeMask_WP))

	# Plot
	Plot_SSTRangeMask(\
		{'SSTRangeMask': SSTRangeMask_IO, 'Lat': RefInfo['Lat'], 'Lon': RefInfo['Lon']}, \
		{'Plot_Range': 'IO'}, \
	)
	Plot_SSTRangeMask(\
		{'SSTRangeMask': SSTRangeMask_WP, 'Lat': RefInfo['Lat'], 'Lon': RefInfo['Lon']}, \
		{'Plot_Range': 'WP'}, \
	)
	Plot_SSTRangeMask(\
		{'SSTRangeMask': np.maximum(SSTRangeMask_IO, SSTRangeMask_WP), 'Lat': RefInfo['Lat'], 'Lon': RefInfo['Lon']}, \
		{'Plot_Range': 'IOWP'}, \
	)