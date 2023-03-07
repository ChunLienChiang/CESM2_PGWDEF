"""
tool_LUCForcing.py
==========================
Calculate and output the LUC (land-use change) forcing for PFT in the MC region.
"""

import numpy as np
import netCDF4 as nc
import json
import os
import shutil
import tool_ClimData_Preprocessing as tool_CP

Config = json.load(open('../config.json'))

def Create_LUCForcing(\
		File_flanduse, \
		File_fsurdat, \
		Year_CTL=1960, \
		Year_DEF=2020, \
		Forcing_Range='MC_Analysis', \
	):

	"""
	Calculate the LUC forcing in the MC region
	==========================
	Argument:
		
		File_flanduse (str): the file path + name of flanduse_timeseries
		
		File_fsurdat (str): the file path + name of fsurdat

		Year_CTL (str): the year of the CTL state. Default is 1960

		Year_DEF (str): the year of the DEF state. Default is 2020

	Output:

		PFT_fsurdat_CTL (numpy array): the PFT data in the CTL state

		PFT_fsurdat_DEF (numpy array): the PFT data in the DEF state
	==========================
	"""
	
	# Get the range mask
	RangeMask = tool_CP.Get_RangeMask(Forcing_Range)

	# Connect to nc file
	ncFile_flanduse = nc.Dataset(File_flanduse)
	ncFile_fsurdat  = nc.Dataset(File_fsurdat)

	# Read PFT
	PFT_flanduse    = ncFile_flanduse.variables['PCT_NAT_PFT'][:]
	PFT_fsurdat     = ncFile_fsurdat.variables['PCT_NAT_PFT'][:]
	
	# Read year
	Year_flanduse   = ncFile_flanduse.variables['time'][:]

	# Read LANDFRAC_PFT
	LANDFRAC_PFT_flanduse = ncFile_flanduse.variables['LANDFRAC_PFT'][:]

	# Extract the PFT at the years of CTL and DEF states
	PFT_flanduse_Ref = PFT_flanduse[0, ...]
	PFT_flanduse_CTL = PFT_flanduse[(Year_flanduse==Year_CTL), ...].squeeze()
	PFT_flanduse_DEF = PFT_flanduse[(Year_flanduse==Year_DEF), ...].squeeze()

	# Create fsurdat for the years of CTL and DEF states
	PFT_fsurdat_CTL  = \
		PFT_fsurdat + \
		(PFT_flanduse_CTL - PFT_flanduse_Ref)
	PFT_fsurdat_DEF  = \
		PFT_fsurdat + \
		(PFT_flanduse_CTL - PFT_flanduse_Ref) + \
		(PFT_flanduse_DEF - PFT_flanduse_CTL) * np.where(RangeMask, 1, 0)

	""" ------------------------------------------------------------------------------------------ """
	# Try to print
	PFT_fsurdat_CTL_SptAvg = np.nanmean(np.where(RangeMask&(LANDFRAC_PFT_flanduse>=1), PFT_fsurdat_CTL, np.nan), axis=(-2, -1))
	PFT_fsurdat_DEF_SptAvg = np.nanmean(np.where(RangeMask&(LANDFRAC_PFT_flanduse>=1), PFT_fsurdat_DEF, np.nan), axis=(-2, -1))
	ChangeRatio = np.round(100*(PFT_fsurdat_DEF_SptAvg-PFT_fsurdat_CTL_SptAvg)/PFT_fsurdat_CTL_SptAvg, 3)
	np.set_printoptions(suppress = True)
	print(np.round(PFT_fsurdat_CTL_SptAvg, 4))
	print(np.round(PFT_fsurdat_DEF_SptAvg, 4))

	print('0  Bare Ground:', ChangeRatio[0], '%')
	print('1  Needleleaf evergreen tree - temperate:', ChangeRatio[1], '%')
	print('4  Broadleaf evergreen tree - tropical:', ChangeRatio[4], '%')
	print('6  Broadleaf deciduous tree - tropical:', ChangeRatio[6], '%')
	print('10 Broadleaf deciduous shrub - temperate:', ChangeRatio[10], '%')
	print('13 C3 grass:', ChangeRatio[13], '%')
	print('14 C4 grass:', ChangeRatio[14], '%')

	print(np.sum(PFT_fsurdat_CTL_SptAvg))
	print(np.sum(PFT_fsurdat_DEF_SptAvg))
	""" ------------------------------------------------------------------------------------------ """

	return PFT_fsurdat_CTL, PFT_fsurdat_DEF

def Output_LUCForcing(\
		PFT_fsurdat_CTL, \
		PFT_fsurdat_DEF, \
		Ref_fsurdat=Config['CESM_Inputdata']+'lnd/clm2/surfdata_map/release-clm5.0.18/surfdata_0.9x1.25_hist_78pfts_CMIP6_simyr1850_c190214.nc', \
		FilePath='../src/LUCForcing/', \
	):

	"""
	Output the LUC Forcing
	==========================
	Argument:
		
		PFT_fsurdat_CTL (numpy array): the PFT data in the CTL state

		PFT_fsurdat_DEF (numpy array): the PFT data in the DEF state

		Ref_fsurdat (str): the path of reference fsurdat file

		FilePath (str): the path for the file to be output

	Output:

		None
	==========================
	"""

	# Check whether the file path has existed
	if (not os.path.exists(FilePath)): os.makedirs(FilePath)

	# Set new file names
	FileName_CTL = '{FilePath}{FileName}_CESM2_PGWDEF_CTL.nc'.format(FilePath=FilePath, FileName='_'.join(Ref_fsurdat.split('/')[-1].split('_')[:2]))
	FileName_DEF = '{FilePath}{FileName}_CESM2_PGWDEF_DEF.nc'.format(FilePath=FilePath, FileName='_'.join(Ref_fsurdat.split('/')[-1].split('_')[:2]))
	
	# Copy files
	shutil.copy(Ref_fsurdat, FileName_CTL)
	shutil.copy(Ref_fsurdat, FileName_DEF)

	# Write to nc file
	# Modify the values of PCT_NAT_PFT in new nc file
	# CTL
	ncFile = nc.Dataset(FileName_CTL, mode='r+', format='NETCDF4')
	ncFile.variables['PCT_NAT_PFT'][:] = PFT_fsurdat_CTL
	ncFile.close()

	# DEF
	ncFile = nc.Dataset(FileName_DEF, mode='r+', format='NETCDF4')
	ncFile.variables['PCT_NAT_PFT'][:] = PFT_fsurdat_DEF
	ncFile.close()

	return

if (__name__ == '__main__'):

	"""
	When this file is called, calculate the LUC forcing and write to new nc files.
	"""

	# The file name of flanduse and fsurdat
	CESM_Inputdata_Path = Config['DataPath']['CESM_Inputdata']
	FileName_flanduse   = CESM_Inputdata_Path + 'lnd/clm2/surfdata_map/release-clm5.0.18/landuse.timeseries_0.9x1.25_SSP5-8.5_78pfts_CMIP6_simyr1850-2100_c190214.nc'
	FileName_fsurdat    = CESM_Inputdata_Path + 'lnd/clm2/surfdata_map/release-clm5.0.18/surfdata_0.9x1.25_hist_78pfts_CMIP6_simyr1850_c190214.nc'

	# Calculate LUC forcing
	PFT_fsurdat_CTL, PFT_fsurdat_DEF = Create_LUCForcing(FileName_flanduse, FileName_fsurdat)

	# Output to nc file
	Output_LUCForcing(PFT_fsurdat_CTL, PFT_fsurdat_DEF)
