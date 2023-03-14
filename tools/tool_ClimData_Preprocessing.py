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
	FilePath = str(Config['Data_Path']['Ref_{RefDataset}'.format(RefDataset=RefDataset)]) + '/{FileName}'.format(FileName=FileName)
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
	
	elif (Range == 'Borneo_Analysis'):

		return -5, 8, 108, 120

	elif (Range == 'NH_Analysis'):

		return 0, 90, 0, 360

	elif (Range == 'SH_Analysis'):

		return -90, 0, 0, 360

	else:
		
		raise ValueError('The argument Range is not available.')

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
	
	# Check whether the argument Range is SST range
	if (Range in ['IO_Analysis', 'WP_Analysis', 'IOWP_Analysis']):
		
		RangeMask = nc.Dataset('../src/SSTRangeMask_{Range}.nc'.format(Range=Range.split('_')[0])).variables['SSTRangeMask'][:]
		RangeMask = np.where(RangeMask==1, True, False)
		
	else:

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

def Calc_SpatialAverage(\
		Data, \
		Spatial_Axis=(-2, -1), \
		LatWeighted=True, \
		LandMask=None, \
		RangeMask=None, \
	):

	"""
	Calculate the spatial average of the data
	==========================
	Argument:

		Data (numpy array)

		Spatial_Axis (tuple of int): optional. The axis numbers of spatial (latitude and longtitude) dimension. Default is (-2, -1) (the last two dimensions)

		LatWeighted (boolean): Default is True. If the argument is True, the spatial average would be calculated considering latitude weighted

		LandMask (str): Default is None. If the argument is given, the land mask would be used to filter data

		RangeMask (str): Default is None. If the argument is given, the range mask would be used to filter data

	Output:

		SpatialAverage (numpy array)
	==========================
	"""

	# Check whether the dimensions are correct
	if not (isinstance(Spatial_Axis, tuple)):

		raise ValueError('The argument Spatial_Axis must be a tuple.')
	
	if (len(Spatial_Axis) != 2):

		raise ValueError('The length of the argument Spatial_Axis must be 2.')

	# Create latitude weights
	Weight_Lat = np.cos(np.deg2rad(Get_RefInfo()['Lat']))[:, None]
	
	# Create land mask weights
	if (LandMask is None):
		
		Weight_LandMask = np.full((Data.shape[Spatial_Axis[0]], Data.shape[Spatial_Axis[1]]), 1)

	else:

		Weight_LandMask = np.where(Get_LandMask(MaskType=LandMask), 1, 0)
	
	# Create range mask weights
	if (RangeMask is None):
		
		Weight_RangeMask = np.full((Data.shape[Spatial_Axis[0]], Data.shape[Spatial_Axis[1]]), 1)

	else:

		Weight_RangeMask = np.where(Get_RangeMask(RangeMask), 1, 0)
	
	# Calculate spatial average
	if (LatWeighted):

		Weights = np.broadcast_to((Weight_Lat*Weight_LandMask*Weight_RangeMask)[None, :], Data.shape)
	
	else:

		Weights = np.broadcast_to((Weight_LandMask*Weight_RangeMask)[None, :], Data.shape)
	
	SpatialAverage = np.ma.average(np.ma.MaskedArray(Data, mask=np.isnan(Data)), axis=Spatial_Axis, weights=Weights)
		
	return SpatialAverage

def Calc_SeasonalCycle(Data, Time_Axis=0, Data_Frequency='Monthly'):

	"""
	Calculate the seasonal cycle of the data
	==========================
	Argument:

		Data (numpy array)

		Time_Axis (int): optional. The axis number of time dimension. Default is 0 (the first dimension)

		Data_Frequency (str): optional. 'Monthly' or 'Daily'. The frequency of time dimension of the data
	
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

		raise ValueError(\
			'If the argument Data_Frequency is monthly, '\
			'the time dimension of Data must be a multiple of 12. '\
			'If the argument Data_Frequency is daily, '\
			'the time dimension of Data must be a multiple of 365.'\
		)

	# Calculate seasonal cycle
	SeasonalCycle = np.nanmean(Data.reshape(tuple([*Data.shape[:Time_Axis], Data.shape[Time_Axis]//Time_Period, Time_Period, *Data.shape[Time_Axis+1:]])), axis=Time_Axis)

	return SeasonalCycle

def Calc_AnnualMean(Data, Time_Axis=0, Data_Frequency='Monthly'):

	"""
	Calculate the annual mean of the data
	==========================
	Argument:

		Data (numpy array)

		Time_Axis (int): optional. The axis number of time dimension. Default is 0 (the first dimension)

		Data_Frequency (str): optional. 'Monthly' or 'Daily'. The frequency of time dimension of the data
	
	Output:

		AnnualMean (numpy array)
	==========================
	"""

	# Check whether the dimensions are correct
	if (not isinstance(Time_Axis, int)):

		raise ValueError('The argument Time_Axis must not integer.')
	
	# Set the period of the annual mean
	if (Data_Frequency == 'Monthly'):

		Time_Period = 12

	if (Data_Frequency == 'Daily'):

		Time_Period = 365
	
	# Check whether the time dimension is a multiple of 12 (or 365) if Data_Frequency is monthly (or daily)
	if (Data.shape[Time_Axis] % Time_Period != 0) & (Data_Frequency == 'Monthly'):

		raise ValueError(\
			'If the argument Data_Frequency is monthly, '\
			'the time dimension of Data must be a multiple of 12. '\
			'If the argument Data_Frequency is daily, '\
			'the time dimension of Data must be a multiple of 365.'\
		)

	# Calculate annual mean
	AnnualMean = np.nanmean(Data.reshape(tuple([*Data.shape[:Time_Axis], Data.shape[Time_Axis]//Time_Period, Time_Period, *Data.shape[Time_Axis+1:]])), axis=Time_Axis+1)

	return AnnualMean