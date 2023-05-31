"""
Calc.LUCForcing.py
==========================
Calculate and output the LUC (land-use change) forcing for PFT in the MC region.
"""

import numpy as np
import xarray as xr
import json
import os
import shutil
import sys
sys.path.append('../')
import tools.tool_ClimData_Preprocessing as tool_CP

Config = json.load(open('../config.json'))

def Create_LUCForcing(\
		File_flanduse, \
		Year_Start=1960, \
		Forcing_Range='MC_Analysis', \
	):

	"""
	Calculate the LUC forcing in the MC region
	==========================
	Argument:
		
		File_flanduse (str): the file path + name of flanduse_timeseries
		
		Year_Start (int): the starting year of detrend. Default is 1960

		Forcing_Range (str): the range of forcing. Default is 'MC_Analysis'

	Output:

		flanduse_fixLUC (numpy array): the PFT data that PFTs have the same value as the Year_Start
	==========================
	"""
	
	# Get the range mask
	Mask_Range      = tool_CP.Get_RangeMask(Forcing_Range)

	# Connect to nc file
	ncFile_flanduse = xr.open_dataset(File_flanduse)
	PFT_flanduse    = ncFile_flanduse['PCT_NAT_PFT'].values
	Year_flanduse   = ncFile_flanduse['time'].values

	# Set the value after the Year_Start to be the same as the Year_Start
	flanduse_fixLUC = PFT_flanduse.copy()
	flanduse_fixLUC = np.where((Year_flanduse >= Year_Start)[:, None, None, None] & Mask_Range[None, None, ...], PFT_flanduse[Year_Start==Year_flanduse, ...], flanduse_fixLUC)

	return flanduse_fixLUC

def Output_LUCForcing(\
		FileName_flanduse_CTL, \
		flanduse_fixLUC, \
		FilePath='../src/LUCForcing/', \
	):

	"""
	Output the LUC Forcing
	==========================
	Argument:
		
		FileName_flanduse_Ref (str): the file path + name of flanduse_timeseries

		flanduse_Detrend (numpy array): the detrended flanduse data

		FilePath (str): the path for the file to be output

	Output:

		None
	==========================
	"""

	# Check whether the file path has existed
	if (not os.path.exists(FilePath)): os.makedirs(FilePath)

	# Set new file names
	FileName_CTL    = '{FilePath}{FileName}_CESM2_PGWDEF_CTL.nc'.format(FilePath=FilePath, FileName='_'.join(FileName_flanduse_CTL.split('/')[-1].split('_')[:2]))
	FileName_fixLUC = '{FilePath}{FileName}_CESM2_PGWDEF_fixLUC.nc'.format(FilePath=FilePath, FileName='_'.join(FileName_flanduse_CTL.split('/')[-1].split('_')[:2]))
	
	# Copy files
	shutil.copy(FileName_flanduse_CTL, FileName_CTL)
	shutil.copy(FileName_flanduse_CTL, FileName_fixLUC)

	# Write to nc file: Modify the values of PCT_NAT_PFT in new FileName_fixLUC nc file
	ncFile_flanduse_fixLUC = xr.open_dataset(FileName_fixLUC)
	ncFile_flanduse_fixLUC['PCT_NAT_PFT'][:] = flanduse_fixLUC
	ncFile_flanduse_fixLUC.to_netcdf(FileName_fixLUC)

	return

if (__name__ == '__main__'):

	# The file name of flanduse and fsurdat
	CESM_Inputdata_Path = Config['Data_Path']['CESM_Inputdata']
	FileName_flanduse   = CESM_Inputdata_Path + 'lnd/clm2/surfdata_map/release-clm5.0.18/landuse.timeseries_0.9x1.25_SSP5-8.5_78pfts_CMIP6_simyr1850-2100_c190214.nc'

	# Calculate LUC forcing
	flanduse_fixLUC = Create_LUCForcing(FileName_flanduse)

	# Output to nc file
	Output_LUCForcing(FileName_flanduse, flanduse_fixLUC)