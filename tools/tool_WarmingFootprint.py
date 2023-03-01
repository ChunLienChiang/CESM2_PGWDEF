"""
tool_WarmingFootprint.py
==========================
Calculate and output the warming footprint for TS over the ocean region.
"""

import numpy as np
import netCDF4 as nc
import os
import statsmodels.api as sm
import tool_ClimData_Preprocessing as tool_CP

def Create_WarmingFootprint(Data_TS, Data_Time=None, Data_TimeRange=None):

	"""
	Create the warming footprint
	==========================
	Argument:
		
		Data_TS (numpy array): 3-d array of TS data
		
		Data_Time (numpy array): optional, 1-d array of time

	Output:

		WarmingFootprint (numpy array): 2-d array of the warming footprint
	==========================
	"""

	# Check whether the dimensions are all correct
	if (np.ndim(Data_TS) != 3):

		raise ValueError('The number of the argument Data_TS must be 3.')

	# Check whether the argument Data_Time is given
	if (not Data_Time is None):

		if (np.ndim(Data_Time) != 1):

			raise ValueError('The number of the argument Data_Time must be 1.')

		elif (Data_Time.shape[0] != Data_TS.shape[0]):

			raise ValueError('The length of the first dimension of the argument Data_TS must equal to the length of the argument Data_Time.')
	
	else:

		Data_Time = np.arange(Data_TS.shape[0])
	
	# Check whether the argument Data_TimeRange is correct
	if (not Data_TimeRange is None):

		if (not isinstance(Data_TimeRange, list)):

			raise ValueError('The argument Data_TimeRange must be a list.')

		if not (isinstance(Data_TimeRange[0], str) and isinstance(Data_TimeRange[1], str)):

			raise ValueError('The items in the argument Data_TimeRange must be strings.')

		# Crop data along time dimension
		Data_TS   = Data_TS[(Data_Time>=int(Data_TimeRange[0]))&(Data_Time<=int(Data_TimeRange[1])), ...]
		Data_Time = Data_Time[(Data_Time>=int(Data_TimeRange[0]))&(Data_Time<=int(Data_TimeRange[1])), ...]
	
	# Create new arrays to save regression coefficients and p-values
	Footprint        = np.full(Data_TS.shape[1:], np.nan)
	Footprint_pvalue = np.full(Data_TS.shape[1:], np.nan)

	# Get land mask
	LandMask         = tool_CP.Get_LandMask(MaskType='Ocean')

	# Calculate the linear regression of each grid point
	for ind_Lat in np.arange(Data_TS.shape[-2]):

		for ind_Lon in np.arange(Data_TS.shape[-1]):
			
			# Check whether the grid point is in the ocean region
			if (LandMask[ind_Lat, ind_Lon] is False): continue

			# OLS
			Model_Results                      = sm.OLS(Data_TS[:, ind_Lat, ind_Lon], sm.add_constant(Data_Time)).fit()
			Footprint[ind_Lat, ind_Lon]        = Model_Results.params[1]
			Footprint_pvalue[ind_Lat, ind_Lon] = Model_Results.pvalues[1]

	return Footprint, Footprint_pvalue

def Output_WarmingFootprint(Footprint, Footprint_pvalue, FilePath='../src/WarmingFootprint/'):

	"""
	Output the warming footprint
	==========================
	Argument:
		
		Footprint (numpy array): warming footprint
		
		Footprint_pvalue (numpy array): p-values of warming footprint

		FilePath (str): the path for the file to be output

	Output:

		None
	==========================
	"""

	# Check whether the file path has existed
	if (not os.path.exists(FilePath)): os.makedirs(FilePath)

	# Write to nc file
	# Create new file
	ncFile       = nc.Dataset('{FilePath}/WarmingFootprint.nc'.format(FilePath=FilePath), mode='w', format='NETCDF4')
	ncFile.title = 'Warming footprint'
	
	# Create new dimensions
	Dim_Lat = ncFile.createDimension('Lat', Footprint.shape[0])
	Dim_Lon = ncFile.createDimension('Lon', Footprint.shape[1])

	# Create new variables
	Var_Footprint = ncFile.createVariable('WarmingFootprint', np.float32, ('Lat', 'Lon'))
	Var_Footprint.long_name = 'warming footprint (the linear trend of TS)'
	Var_Footprint.units = 'K per month'
	Var_Footprint[:] = Footprint

	Var_Footprint_pvalue = ncFile.createVariable('WarmingFootprint_pvalue', np.float32, ('Lat', 'Lon'))
	Var_Footprint_pvalue.long_name = 'pvalues of warming footprint'
	Var_Footprint_pvalue[:] = Footprint_pvalue

	ncFile.close()

	return

if (__name__ == '__main__'):

	"""
	When this file is called, calculate the warming footprint and write to new nc files.
	"""

	Data_TimeRange = [\
		'19600101', \
		'20141231', \
	]

	# Get TS data from reference dataset (FHIST)
	Data_TS = tool_CP.Get_RefData(\
		FileName='f.e21.FHIST_BGC.f09_f09.historical.ersstv5.goga.ens01.cam.h0.TS.188001-201412.nc', \
		Var='TS', \
	)
	Data_Date = tool_CP.Get_RefData(\
		FileName='f.e21.FHIST_BGC.f09_f09.historical.ersstv5.goga.ens01.cam.h0.TS.188001-201412.nc', \
		Var='date', \
	)

	# Calculate warming footprint
	Footprint, Footprint_pvalue = Create_WarmingFootprint(Data_TS, Data_Date, Data_TimeRange)

	# Output to nc file
	Output_WarmingFootprint(Footprint, Footprint_pvalue)