"""
tool_SSTRangeMask.py
==========================
Calculate and output the range mask that can be used in pacemaker simulation.
"""

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib as mpl
import cartopy.crs as ccrs
import os
import tool_ClimData_Preprocessing as tool_CP

def Create_SSTRangeMask_IO(Lat, Lon, LandFraction):

	"""
	Create the range mask for Indian Ocean (IO) region
	==========================
	Argument:
		
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
	Argument:
		
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

def Output_SSTRangeMask(Range, SSTRangeMask, FilePath='../src/SSTRangeMask/'):

	"""
	Output the range mask
	==========================
	Argument:
		
		Range (str): the name of the range. This would be contained in the output file name
		
		SSTRangeMask (numpy array): range mask

		FilePath (str): the path for the file to be output

	Output:

		None
	==========================
	"""

	# Check whether the file path has existed
	if (not os.path.exists(FilePath)): os.makedirs(FilePath)

	# Write to nc file
	# Create new file
	ncFile       = nc.Dataset('{FilePath}/SSTRangeMask_{Range}.nc'.format(FilePath=FilePath, Range=Range), mode='w', format='NETCDF4')
	ncFile.title = 'Range Mask: {Range}'.format(Range=Range)
	
	# Create new dimensions
	Dim_Lat = ncFile.createDimension('Lat', SSTRangeMask.shape[0])
	Dim_Lon = ncFile.createDimension('Lon', SSTRangeMask.shape[1])

	# Create new variables
	Var_SSTRangeMask = ncFile.createVariable('SSTRangeMask', np.float32, ('Lat', 'Lon'))
	Var_SSTRangeMask.long_name = 'range mask (1=True, 0=False)'
	Var_SSTRangeMask[:] = SSTRangeMask

	ncFile.close()

	return

def Plot_SSTRangeMask(Plot_Data, Plot_Config):

	Figure_Path = '../output/Output_Figure/src_Map_SSTRegionMask/'
	Figure_Name = 'src_Map_SSTRegionMask_{Range}.png'.format(Range=Plot_Config['Plot_Range'])

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

	# Output to nc file
	Output_SSTRangeMask('IO', SSTRangeMask_IO)
	Output_SSTRangeMask('WP', SSTRangeMask_WP)
	Output_SSTRangeMask('IOWP', SSTRangeMask_IO | SSTRangeMask_WP)

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
		{'SSTRangeMask': SSTRangeMask_IO | SSTRangeMask_WP, 'Lat': RefInfo['Lat'], 'Lon': RefInfo['Lon']}, \
		{'Plot_Range': 'IOWP'}, \
	)