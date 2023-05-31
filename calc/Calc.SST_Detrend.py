"""
Calc.SST_Detrend.py
===============================
1. Remove the linear trend of SST between 1960-2020 and save the detrended SST as a new dataset.
2. Plot the line chart of SST and detrended SST.
3. Output the detrended SST as a new dataset (under the new directory).
4. Output the regrided detrended/original SST as a new nc file (under the new directory).
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr
import json
import os
import datetime as dt

Config = json.load(open('../config.json'))

def Get_Data():

	"""
	Get ModifiedSST_ERSST_gx1v7 data
	=========================================
	Output:
		Data: ModifiedSST_ERSST_gx1v7 data
	"""

	# Print message
	print('Get ModifiedSST_ERSST_gx1v7 data')

	# Path of ModifiedSST_ERSST_gx1v7 data
	Data_Path = Config['Data_Path']['SST_Pacemaker'] + Config['Data_Name']['SST_Pacemaker']['CTL'] + '/'

	# Read all nc files that ends with number under the path
	File_List = [i for i in os.listdir(Data_Path) if i.split('.')[-1].isnumeric()]
	File_List.sort()

	# Create new lists to store data
	Data            = []
	Data_OriginTime = []
	Data_NewTime    = []

	for i_File in File_List:
		
		# Read data: SST
		Data_New      = xr.open_dataset(Data_Path + i_File)

		# Read data: Time from file name
		Data_NewTime_New = '-'.join(i_File.split('.')[-3:-1])
		Data_NewTime_New = pd.Timestamp(dt.datetime.strptime(Data_NewTime_New.split('-')[0] + '-' + Data_NewTime_New.split('-')[-1], '%Y-%j'))

		# Append data
		Data.append(Data_New)
		Data_OriginTime.append('.'.join(i_File.split('.')[-3:]))
		Data_NewTime.append(Data_NewTime_New)
	
	# Write new nc file with combined Data and Data_NewTime
	Data = xr.concat(Data, dim = 'time')
	Data = Data.assign_coords(time = Data_NewTime)

	# Add OriginTime as a new variable "filename"
	Data['filename'] = (['time'],  Data_OriginTime)

	# Output concatenated data to new nc file
	if (os.path.exists(Data_Path + 'ModifiedSST_ERSST_gx1v7.nc.All')): os.remove(Data_Path + 'ModifiedSST_ERSST_gx1v7.nc.All')
	Data.to_netcdf(Data_Path + 'ModifiedSST_ERSST_gx1v7.nc.All')

	return Data

def Remove_Trend(Data):

	"""
	Remove the linear trend of SST between 1960-2020 and shift the data to make the first data after 1960-01-01 to be equal to the first data after 1960-01-01 in original SST data.
	====================================================================
	Input:
		Data: (xarray dataset) ModifiedSST_ERSST_gx1v7 data.
	Output:
		Data_Modified_SST: (xarray dataset) Detrended SST.
	"""

	# Print message
	print('Remove the trend')

	# Extract the data between 1960-2020
	Data_Modified      = Data.copy(deep=True)
	Data_Modified      = Data_Modified.sel(time = slice('1960-01-01', '2020-01-31'))
	Data_Modified_Time = Data_Modified['time'].values
	Data_Modified_SST  = Data_Modified['SST'].values

	# Extract time and convert it to "the days since 1960-01-01"
	Data_Modified_Time = np.array([(i - pd.Timestamp('1960-01-01')).days for i in Data_Modified_Time])

	# Regress each grid point of SST on Data_Modified_Time
	for i_Lat in range(Data_Modified_SST.shape[1]):

		for i_Lon in range(Data_Modified_SST.shape[2]):

			# Regress SST on Data_Modified_Time
			Reg = np.polyfit(Data_Modified_Time, Data_Modified_SST[:, i_Lat, i_Lon], 1)

			# Remove the linear trend
			Data_Modified_SST[:, i_Lat, i_Lon] = Data_Modified_SST[:, i_Lat, i_Lon] - (Reg[0] * Data_Modified_Time + Reg[1])

	# Shift the Data_Modified_SST to make the first data after 1960-01-01 to be equal to the first data after 1960-01-01 in original SST data
	Data_Modified_SST = \
		Data_Modified_SST + \
		(Data['SST'].sel(time = slice('1960-01-01', '1965-01-01')).values[0, ...] - Data_Modified_SST[0, ...])[None, ...]

	# Concatenate the original SST data before 1960-01-01 and the detrended SST data after 1960-01-01
	Data_Modified_SST = np.concatenate([Data['SST'].sel(time = slice(None, '1960-01-01')).values, Data_Modified_SST], axis = 0)

	# Convert Data_Modified_SST to xarray dataset
	Data_Modified_SST = xr.DataArray(Data_Modified_SST, coords = Data['SST'].coords, dims = Data['SST'].dims)
	Data_Modified_SST = Data_Modified_SST.to_dataset(name = 'SST')

	return Data_Modified_SST

def Plot_Linechart(Plot_Data, Plot_Config):

	# Print message
	print('Plot line chart')

	# Create figure
	fig, ax = plt.subplots(figsize=(5, 5), dpi=300)

	# Plot: SST
	ax.plot(Plot_Data['Time'], Plot_Data['SST'], color='black', linewidth=1.0, label='SST', zorder=2)
	ax.plot(Plot_Data['Time'], Plot_Data['SST_Detrend'], color='red', linewidth=1.0, label='SST_Detrend', zorder=1)

	# Plot configuration
	ax.set_xlim([Plot_Data['Time'][0]-1, Plot_Data['Time'][-1]+1])
	ax.set_xlabel('Year')
	ax.set_ylabel('SST (degC)')
	ax.legend(loc='upper left', fontsize=10)

	# Save figure
	plt.tight_layout()
	Output_Path = '../output/Output_Figure/Plot.Lineplot.SST_Detrend/'
	Output_Name = 'Plot.Lineplot.SST_Detrend.png'
	if not os.path.exists(Output_Path): os.makedirs(Output_Path)
	plt.savefig(Output_Path + Output_Name, bbox_inches='tight')

	return

def Output_Data(Data, File_List):

	# Print message
	print('Output data')

	# Set output path
	Output_Path = '../src/SST/ModifiedSST_ERSST_gx1v7_Detrend/'
	if not os.path.exists(Output_Path): os.makedirs(Output_Path)

	# Output data
	File_List = ['ModifiedSST_ERSST_gx1v7.nc.' + i for i in File_List]
	
	for ind_File in range(len(File_List)):

		# Output data
		xr.DataArray(
			data=Data['SST'].values[ind_File, ...],
			dims=['nlat', 'nlon'],
		).to_dataset(name='SST').to_netcdf(Output_Path + File_List[ind_File])
	
	# Output concatenated data to new nc file
	if (os.path.exists(Output_Path + 'ModifiedSST_ERSST_gx1v7.nc.All')): os.remove(Output_Path + 'ModifiedSST_ERSST_gx1v7.nc.All')
	xr.DataArray(
		data=Data['SST'].values,
		dims=['time', 'nlat', 'nlon'],
		coords={'time': Data['time'].values},
	).to_dataset(name='SST').to_netcdf(Output_Path + 'ModifiedSST_ERSST_gx1v7.nc.All')

	return

if (__name__ == '__main__'):

	# Read data: SST
	Data = Get_Data()

	# Remove the linear trend of SST between 1960-2020
	Data_Detrend = Remove_Trend(Data)

	# Output data
	Output_Data(Data_Detrend, Data.filename.values)

	# ================================================================
	# Plot two line on the same line chart
	Plot_Data = {
		'Time'       : np.arange(1880, 2020), \
		'SST'        : np.nanmean(np.nanmean(Data['SST'].values, axis=(1,2))[1:-1].reshape((-1,12)), axis=-1), \
		'SST_Detrend': np.nanmean(np.nanmean(Data_Detrend['SST'].values, axis=(1,2))[1:-1].reshape((-1,12)), axis=-1), \
	}

	Plot_Config = {}

	Plot_Linechart(Plot_Data, Plot_Config)