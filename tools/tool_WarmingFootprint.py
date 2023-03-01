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

def Create_WarmingFootprint(Data_TS):

	"""
	Create the warming footprint
	==========================
	Argument:
		
		Data_TS (numpy array): 3-d array of TS data, monthly data
		
	Output:

		WarmingFootprint (numpy array): 2-d array of the warming footprint
	==========================
	"""

	# Check whether the dimensions are all correct
	if (np.ndim(Data_TS) != 3):

		raise ValueError('The number of the argument Data_TS must be 3.')
	
	if (Data_TS.shape[0] % 12 != 0):

		raise ValueError('The length of the first dimension of the argument Data_TS must be a multiple of 12.')

	# Create new arrays to save regression coefficients and p-values
	Footprint        = np.full((12, *Data_TS.shape[1:]), np.nan)
	Footprint_pvalue = np.full((12, *Data_TS.shape[1:]), np.nan)

	# Get land mask
	LandMask         = tool_CP.Get_LandMask(MaskType='Ocean')

	# Calculate the linear regression of each grid point
	for ind_Month in np.arange(12):

		for ind_Lat in np.arange(Data_TS.shape[-2]):

			for ind_Lon in np.arange(Data_TS.shape[-1]):
				
				# Check whether the grid point is in the ocean region
				if (LandMask[ind_Lat, ind_Lon] is False): continue

				# OLS
				Model_Results                                 = sm.OLS(Data_TS[ind_Month::12, ind_Lat, ind_Lon], sm.add_constant(np.arange(Data_TS.shape[0]//12))).fit()
				Footprint[ind_Month, ind_Lat, ind_Lon]        = Model_Results.params[1]
				Footprint_pvalue[ind_Month, ind_Lat, ind_Lon] = Model_Results.pvalues[1]

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
	Dim_Month = ncFile.createDimension('Month', Footprint.shape[0])
	Dim_Lat   = ncFile.createDimension('Lat', Footprint.shape[-2])
	Dim_Lon   = ncFile.createDimension('Lon', Footprint.shape[-1])

	# Create new variables
	Var_Footprint = ncFile.createVariable('WarmingFootprint', np.float32, ('Month', 'Lat', 'Lon'))
	Var_Footprint.long_name = 'warming footprint (the linear trend of TS)'
	Var_Footprint.units = 'K per year'
	Var_Footprint[:] = Footprint

	Var_Footprint_pvalue = ncFile.createVariable('WarmingFootprint_pvalue', np.float32, ('Month', 'Lat', 'Lon'))
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
	Data_Time = tool_CP.Get_RefData(\
		FileName='f.e21.FHIST_BGC.f09_f09.historical.ersstv5.goga.ens01.cam.h0.TS.188001-201412.nc', \
		Var='date', \
	)
	
	Data_TS   = Data_TS[(Data_Time>=int(Data_TimeRange[0]))&(Data_Time<=int(Data_TimeRange[1])), ...]
	
	# Calculate warming footprint
	Footprint, Footprint_pvalue = Create_WarmingFootprint(Data_TS)

	# Output to nc file
	Output_WarmingFootprint(Footprint, Footprint_pvalue)