"""
tool_Regrid.py
==============================
Convert pop grid to 0.9x1.25 grid using scipy.interpolate.griddata
"""

import numpy as np
import xarray as xr
import datetime as dt
import scipy.interpolate as sciinterp
import multiprocessing as mp
import json
import pop_tools as popt
import sys
sys.path.append('../')
from tools import tool_GetClimData as tool_GetCD

Config = json.load(open('../config.json'))

def _process_Interp(Data, Time):

	# Print time
	print('Calculate interpolation: ' + str(Time).split('T')[0])

	# Create interpolator
	Interp_Linear  = sciinterp.LinearNDInterpolator((POP_Grid_Lat.flatten(), POP_Grid_Lon.flatten()), Data.flatten())
	Interp_Nearest = sciinterp.NearestNDInterpolator((POP_Grid_Lat.flatten(), POP_Grid_Lon.flatten()), Data.flatten())

	# Interpolation
	Data_Interp_Linear = Interp_Linear(LatLon_Grid_Lat, LatLon_Grid_Lon)
	Data_Interp_Nearest = Interp_Nearest(LatLon_Grid_Lat, LatLon_Grid_Lon)
	Data_Interp = np.where(np.isnan(Data_Interp_Linear), Data_Interp_Nearest, Data_Interp_Linear)

	return Data_Interp

if (__name__ == '__main__'):

	# Set input grid: pop gx1v7
	POP_Grid = popt.get_grid('POP_gx1v7')
	POP_Grid_Lat = POP_Grid.TLAT.values
	POP_Grid_Lon = POP_Grid.TLONG.values

	# Set output grid: 0.9x1.25
	LatLon_Grid = tool_GetCD.Get_Grid('0.9x1.25')
	LatLon_Grid_Lat = LatLon_Grid['Lat']
	LatLon_Grid_Lon = LatLon_Grid['Lon']

	for i_Dataset in ['CTL', 'fixSST']:
	   
		# Read data
		Data = tool_GetCD.Get_SST_Pacemaker(i_Dataset)
		Data_Time = tool_GetCD.Get_SST_Pacemaker(i_Dataset, Var='time')

		# Create output array
		Data_Interp = np.full((Data.shape[0], *LatLon_Grid_Lat.shape), np.nan)

		Pool_Results = mp.Pool(30).starmap(\
			_process_Interp, \
			zip([np.squeeze(i) for i in np.split(Data, Data.shape[0], axis=0)], Data_Time), \
		)
		Data_Interp = np.concatenate([i[None, ...] for i in Pool_Results], axis=0)

		# Save data
		xr.Dataset({\
			'SST': (['Time', 'Lat', 'Lon'], Data_Interp),\
			'Lat': (['Lat'], LatLon_Grid_Lat[:, 0]),\
			'Lon': (['Lon'], LatLon_Grid_Lon[0, :]),\
			'Time': (['Time'], Data_Time)\
		}).to_netcdf(Config['Data_Path']['SST_Pacemaker'] + Config['Data_Name']['SST_Pacemaker'][i_Dataset] + '/ModifiedSST_ERSST_gx1v7.nc.All_0.9x1.25')