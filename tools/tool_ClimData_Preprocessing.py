"""
tool_ClimData_Preprocessing.py
==========================
The tools used to process the climate dataset.
"""

import numpy as np
import netCDF4 as nc
import json

Config = json.load(open('../config.json'))

def Connect_RefDataset(\
		RefDataset='CESM2_BHIST_Historical', \
		FileName='b.e21.BHIST.f09_g17.CMIP6-historical.003.cam.h0.LANDFRAC.185001-201412.nc', \
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
		'Lat': Connect_RefDataset(Var='lat'), \
		'Lon': Connect_RefDataset(Var='lon'), \
		'LandFraction': Connect_RefDataset(Var='LANDFRAC'), \
	}
	
	return RefInfo

def Get_LandMask(LandFraction, MaskType='Land'):
	
	"""
	Get land mask from land fraction data.
	==========================
	Argument:

		LandFraction (numpy array): the land fraction data

		MaskType (str): 'Land' for land region (default); 'Ocean' for ocean region

	Output:

		LandMask (numpy array): the mask of land (or ocean)
	==========================
	"""

	if (MaskType == 'Land'):

		LandMask = np.where(LandFraction>=0.5, True, False)

	elif (MaskType == 'Ocean'):

		LandMask = np.where(LandFraction<0.5, True, False)
	
	else:

		raise ValueError('The argument MaskType must be Land or Ocean')

	return LandMask