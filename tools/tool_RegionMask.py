"""
tool_RegionMask.py
==========================
Calculate and output the regional mask that can be used in pacemaker simulation.
"""

import numpy as np
import netCDF4 as nc
import os
import tool_ClimData_Preprocessing as tool_CP

def Create_RegionMask_IO(Lat, Lon, LandFraction):

	"""
	Create the regional mask for Indian Ocean (IO) region
	==========================
	Argument:
		
		Lat (numpy array): latitude
		
		Lon (numpy array): longitude

		LandFraction (numpy array): land fraction

	Output:

		RegionMask (numpy array): the regional mask
	==========================
	"""

	if (LandFraction.shape != (Lat.size, Lon.size)):

		raise ValueError('The shape of LandFraction must equal to (latitude, longitude)')
	
	# Fill in a rectangle region
	RegionMask = np.where((Lat[:, None]>=-30)&(Lat[:, None]<=30)&(Lon[None, :]>=30)&(Lon[None, :]<=120), True, False)

	# Exclude the north-east corner of array
	RegionMask = np.where((Lat[:, None]>=-8)&(Lon[None, :]>=105), False, RegionMask)
	RegionMask = np.where((Lat[:, None]>=3)&(Lon[None, :]>=98), False, RegionMask)
	RegionMask = np.where((Lat[:, None]>=-1)&(Lon[None, :]>=101), False, RegionMask)

	# Land mask
	LandMask   = tool_CP.Get_LandMask(LandFraction, MaskType='Ocean')
	RegionMask = np.where(LandMask, RegionMask, False)
	
	return RegionMask

def Create_RegionMask_WP(Lat, Lon, LandFraction):

	"""
	Create the regional mask for Western Pacific (WP) region
	==========================
	Argument:
		
		Lat (numpy array): latitude
		
		Lon (numpy array): longitude

		LandFraction (numpy array): land fraction

	Output:

		RegionMask (numpy array): the regional mask
	==========================
	"""

	if (LandFraction.shape != (Lat.size, Lon.size)):

		raise ValueError('The shape of LandFraction must equal to (latitude, longitude)')
	
	# Fill in a rectangle region
	RegionMask = np.where((Lat[:, None]>=-20)&(Lat[:, None]<=20)&(Lon[None, :]>=100)&(Lon[None, :]<=170), True, False)

	# Exclude the north-east corner of array
	RegionMask = np.where((Lat[:, None]<=-10)&(Lon[None, :]<=120), False, RegionMask)
	RegionMask = np.where((Lat[:, None]<=-8)&(Lon[None, :]<=114), False, RegionMask)
	RegionMask = np.where((Lat[:, None]<=-6)&(Lon[None, :]<=108), False, RegionMask)
	RegionMask = np.where((Lat[:, None]<=-5)&(Lon[None, :]<=104), False, RegionMask)
	RegionMask = np.where((Lat[:, None]<=-2)&(Lon[None, :]<=103), False, RegionMask)

	# Land mask
	LandMask   = tool_CP.Get_LandMask(LandFraction, MaskType='Ocean')
	RegionMask = np.where(LandMask, RegionMask, False)
	
	return RegionMask

def Output_RegionMask(Region, RegionMask, FilePath='../src/RegionMask/'):

	"""
	Output the regional mask
	==========================
	Argument:
		
		Region (str): the name of the region. This would be contained in the output file name
		
		RegionMask (numpy array): regional mask

		FilePath (str): the path for the file to be output

	Output:

		None
	==========================
	"""

	# Check whether the file path has existed
	if (not os.path.exists(FilePath)): os.makedirs(FilePath)

	# Write to nc file
	# Create new file
	ncFile       = nc.Dataset('{FilePath}/RegionMask_{Region}.nc'.format(FilePath=FilePath, Region=Region), mode='w', format='NETCDF4')
	ncFile.title = 'Regional Mask: {Region}'.format(Region=Region)
	
	# Create new dimensions
	Dim_Lat = ncFile.createDimension('Lat', RegionMask.shape[0])
	Dim_Lon = ncFile.createDimension('Lon', RegionMask.shape[1])

	# Create new variables
	Var_RegionMask = ncFile.createVariable('RegionMask', np.float32, ('Lat', 'Lon'))
	Var_RegionMask.long_name = 'regional mask'
	Var_RegionMask[:] = RegionMask

	ncFile.close()

	return

if (__name__ == '__main__'):

	RefInfo = tool_CP.Get_RefInfo()

	RegionMask_IO = Create_RegionMask_IO(RefInfo['Lat'], RefInfo['Lon'], RefInfo['LandFraction'])
	Output_RegionMask('IO', RegionMask_IO)

	RegionMask_WP = Create_RegionMask_WP(RefInfo['Lat'], RefInfo['Lon'], RefInfo['LandFraction'])
	Output_RegionMask('WP', RegionMask_WP)