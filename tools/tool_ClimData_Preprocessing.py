"""
tool_ClimData_Preprocessing.py
==========================
The tools used to process the climate dataset.
"""

import numpy as np
import netCDF4 as nc
import json

Config = json.load(open('../config.json'))

def Get_RefData(\
		RefDataset='CESM2_FHIST_BGC_Historical', \
		FileName='f.e21.FHIST_BGC.f09_f09.historical.ersstv5.goga.ens01.cam.h0.LANDFRAC.188001-201412.nc', \
		Var=None, \
	):
	
	"""
	Connect to the reference dataset. The path of required dataset is provided in config.json
	==========================
	Argument:

		RefDataset (str): the name of the reference dataset
		
		FileName (str): the file name of the connected reference dataset
		
		Var (str): the variable to be read

	Output:

		Data (numpy array): the data in the nc file
	==========================
	"""

	# Check whether argument Var is empty
	if (Var is None):

		raise ValueError('The argument Var must not be empty.')

	# Connect to nc file
	FilePath = str(Config['DataPath']['Ref_{RefDataset}'.format(RefDataset=RefDataset)]) + '/{FileName}'.format(FileName=FileName)
	Data = nc.Dataset(FilePath).variables[Var][:].squeeze().filled(np.nan)

	return Data

def Get_RefInfo():
	
	"""
	Get the reference information from the reference dataset. The path of required reference dataset is provided in config.json
	==========================
	Output:

		RefInfo (dict): the reference information
	==========================
	"""

	RefInfo = {\
		'Lat': Get_RefData(Var='lat'), \
		'Lon': Get_RefData(Var='lon'), \
		'LandFraction': Get_RefData(Var='LANDFRAC'), \
	}
	
	return RefInfo

def Get_RangeBoundary(Range):

	"""
	Get the range boundary
	==========================
	Argument:

		Range (str): the range to get boundary

	Output:

		RangeBoundary (tuple): (minimum latitude, maximum latitude, minimum longitude, maximum longitude), respectively
	==========================
	"""

	if (Range == 'MC_Analysis'):
		
		return -10, 10, 90, 140
	
	elif (Range == 'Borneo'):

		return -5, 8, 108, 120

	else:
		
		return None

def Get_RangeMask(Range, Lat=None, Lon=None):

	"""
	Get the range mask
	==========================
	Argument:

		Range (str): the range to generate range mask

		Lat (numpy array): latitude. If there's no value provided, the latitude in the reference dataset would be used

		Lon (numpy array): longitude. If there's no value provided, the longitude in the reference dataset would be used
	
	Output:

		RangeMask (numpy array): range mask
	==========================
	"""

	# Check whether the arguments Lat and Lon are given
	if ((Lat is None) or (Lon is None)):
		
		RefInfo = Get_RefInfo()
		Lat     = RefInfo['Lat']
		Lon     = RefInfo['Lon']
	
	# Get the range boundary
	Lat_Min, Lat_Max, Lon_Min, Lon_Max = Get_RangeBoundary(Range)

	# Create the range mask
	RangeMask = np.where((Lat[:, None]>=Lat_Min)&(Lat[:, None]<=Lat_Max)&(Lon[None, :]>=Lon_Min)&(Lon[None, :]<=Lon_Max), True, False)

	return RangeMask

def Get_LandMask(LandFraction=None, MaskType='Land'):
	
	"""
	Get land mask from land fraction data.
	==========================
	Argument:

		LandFraction (numpy array): the land fraction data. If there's no value provided, the land fraction in the reference dataset would be used

		MaskType (str): 'Land' for land region (default); 'Ocean' for ocean region

	Output:

		LandMask (numpy array): the mask of land (or ocean)
	==========================
	"""
	
	# Check if the argument LandFraction is empty
	if (LandFraction is None):

		LandFraction = Get_RefInfo()['LandFraction']

	# Determine what tpye of mask should be return
	if (MaskType == 'Land'):

		LandMask = np.where(LandFraction>=0.5, True, False)

	elif (MaskType == 'Ocean'):

		LandMask = np.where(LandFraction<0.5, True, False)
	
	else:

		raise ValueError('The argument MaskType must be Land or Ocean')

	return LandMask

def Calc_SeasonalCycle(Data, Time_Axis=0, Data_Frequency='Monthly'):

	"""
	Calculate the seasonal cycle of the data
	==========================
	Argument:

		Data (numpy array)

		Time_Axis (int): optional. The axis number of time dimension. Default is 0 (the first dimension)

		Data_Frequency (str): optional: 'Monthly' or 'Daily'. The frequency of time dimension of the data
	
	Output:

		SeasonalCycle (numpy array)
	==========================
	"""

	# Check whether the dimensions are correct
	if (not isinstance(Time_Axis, int)):

		raise ValueError('The argument Time_Axis must not integer.')
	
	# Set the period of the seasonal cycle
	if (Data_Frequency == 'Monthly'):

		Time_Period = 12

	if (Data_Frequency == 'Daily'):

		Time_Period = 365
	
	# Check whether the time dimension is a multiple of 12 (or 365) if Data_Frequency is monthly (or daily)
	if (Data.shape[Time_Axis] % Time_Period != 0) & (Data_Frequency == 'Monthly'):

		raise ValueError('If the argument Data_Frequency is monthly, the time dimension of Data must be a multiple of 12. If the argument Data_Frequency is daily, the time dimension of Data must be a multiple of 365.')

	# Calculate seasonal cycle
	SeasonalCycle = np.nanmean(Data.reshape(tuple([*Data.shape[:Time_Axis], Data.shape[Time_Axis]//Time_Period, Time_Period, *Data.shape[Time_Axis+1:]])), axis=Time_Axis)

	return SeasonalCycle