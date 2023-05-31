"""
Plot.Map.PFT.py
=============================
Output the figure of the map of PFTs[4].
"""

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import xarray as xr
import os
import sys
sys.path.append('../')
from tools import tool_ClimData_Preprocessing as tool_CDPrep

def Plot_Map(Plot_Data, Plot_Config):
	
	# Create figure object
	fig, ax = plt.subplots(figsize=(5, 3), dpi=300, subplot_kw={'projection': ccrs.PlateCarree()})

	# Plot: PFT
	PFT_pcolormesh = ax.pcolormesh(Plot_Data['Lon'], Plot_Data['Lat'], Plot_Data['PFT'], cmap='PiYG', vmin=-3, vmax=3, transform=ccrs.PlateCarree())
	fig.colorbar(PFT_pcolormesh, label='PFT change (%)', orientation='vertical', shrink=1, fraction=0.05, pad=0.05)

	# Plot: rectangle for MC
	Lat_Min, Lat_Max, Lon_Min, Lon_Max = tool_CDPrep.Get_RangeBoundary('MC_Analysis')
	ax.add_patch(plt.Rectangle(xy=(Lon_Min, Lat_Min), width=(Lon_Max - Lon_Min), height=(Lat_Max - Lat_Min), fill=False, transform=ccrs.PlateCarree(), color='Black', linewidth=1.5))

	# Configuration
	Lat_Min, Lat_Max, Lon_Min, Lon_Max = tool_CDPrep.Get_RangeBoundary('MC_Map')
	ax.set_extent([Lon_Min, Lon_Max, Lat_Min, Lat_Max], crs=ccrs.PlateCarree())
	ax.coastlines(resolution='10m')
	ax.set_title('Map: PFT {}'.format(Plot_Config['n_PFT']), pad=12)

	# Save figure
	plt.tight_layout()
	Figure_Path = '../output/Output_Figure/Plot.Map.PFT/'
	Figure_Name = 'Plot.Map.PFT.{}.png'.format(Plot_Config['n_PFT'])
	if not (os.path.exists(Figure_Path)): os.makedirs(Figure_Path)
	plt.savefig(Figure_Path + Figure_Name)

	return

if (__name__ == '__main__'):

	# Get data
	Data_PFT_CTL    = xr.open_dataset('../src/LUCForcing/landuse.timeseries_0.9x1.25_CESM2_PGWDEF_CTL.nc')['PCT_NAT_PFT'].sel(time=slice('2020', '2020')).values
	Data_PFT_fixLUC = xr.open_dataset('../src/LUCForcing/landuse.timeseries_0.9x1.25_CESM2_PGWDEF_fixLUC.nc')['PCT_NAT_PFT'].sel(time=slice('2020', '2020')).values

	# Calculate spatial average
	Data_PFT_CTL    = np.nanmean(Data_PFT_CTL, axis=0)
	Data_PFT_fixLUC = np.nanmean(Data_PFT_fixLUC, axis=0)

	# Get reference information
	RefInfo         = tool_CDPrep.Get_RefInfo()

	# Plot
	for i_PFT in range(0, 15):

		Plot_Data = {\
			'PFT'  : Data_PFT_CTL[i_PFT, ...] - Data_PFT_fixLUC[i_PFT, ...], \
			'Lat'  : RefInfo['Lat'], \
			'Lon'  : RefInfo['Lon'], \
		}

		Plot_Config = {\
			'n_PFT': i_PFT, \
		}

		Plot_Map(Plot_Data, Plot_Config)