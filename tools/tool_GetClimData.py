"""
tool_GetClimData.py
==========================
The tools obtaining climate raw data.
"""

import numpy as np
import netCDF4 as nc
import os
import json

Config = json.load(open('../config.json'))

def Get(Exp, Var, Component='atm', Data_TimeRange=None, Data_Frequency='Monthly'):
	
	"""
	Get climate raw data from disk
	==========================
	Argument:

		Exp (str): the experiment name to be obtained

		Var (str): the variable name to be obtained

		Component (str): optional. The component name where the variable is output. Default is 'atm'.

		Data_TimeRange (list of int): optional. The list containing start and end year of climate raw data. If there's no value provided, the default value would be used from config.json

		Data_Frequency (str): optional. 'Monthly', 'Daily', or 'Hourly'. Default is 'Monthly'. The frequency of time dimension of the data
	
	Output:

		Data (numpy array)
	==========================
	"""

	# Check whether the argument Exp is available
	if not (Exp.split('_')[2] in Config['Data_Name']['Main']):

		raise ValueError('Given experiment does not exist.')
	
	# Set time range from the argument or from config.json
	if (Data_TimeRange is None):
	
		Data_TimeRange = Config['Data_TimeRange'][Exp.split('_')[2]]

	if not (isinstance(Data_TimeRange, list)):
		
		raise ValueError('The argument Data_TimeRange must be a list.')

	if (len(Data_TimeRange) != 2):
		
		raise ValueError('The length of the argument Data_TimeRange must be 2.')

	# Set time range as a numpy array
	Data_TimeRange = np.arange(Data_TimeRange[0], Data_TimeRange[1]+1)

	# Set data information
	Data_Path = Config['Data_Path']['Main'] + Exp + '/' + Component + '/hist/'

	# Get the file list
	if (Data_Frequency == 'Monthly'):
		
		File_List = [i for i in os.listdir(Data_Path) if (i.split('.')[-3]=='h0') and (int(i.split('.')[-2].split('-')[0]) in Data_TimeRange)]
	
	if (Data_Frequency == 'Daily'):
		
		File_List = [i for i in os.listdir(Data_Path) if (i.split('.')[-3]=='h1') and (int(i.split('.')[-2].split('-')[0]) in Data_TimeRange)]
	
	if (Data_Frequency == 'Daily'):
		
		File_List = [i for i in os.listdir(Data_Path) if (i.split('.')[-3]=='h2') and (int(i.split('.')[-2].split('-')[0]) in Data_TimeRange)]

	# Get climate data
	Data = nc.MFDataset([Data_Path + File for File in File_List]).variables[Var][:]

	return Data