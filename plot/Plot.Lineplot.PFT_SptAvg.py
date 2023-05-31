"""
Plot.Lineplot.PFT_SptAvg.py
==========================
Output the figure of the time series of PFT spatial average.
"""

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import os
import sys
sys.path.append('../')
from tools import tool_ClimData_Preprocessing as tool_CDPrep

def Plot_Lineplot(Plot_Data, Plot_Config):

	# Create figure object
	fig, ax = plt.subplots(figsize=(5, 5), dpi=300)

	# Plot
	x_array = np.arange(1880, 2021)
	ax.plot(x_array, Plot_Data['PFT_CTL'], color='Black', zorder=2, label='CTL')
	ax.plot(x_array, Plot_Data['PFT_fixLUC'], color='Red', zorder=1, label='fixLUC')

	# Configuration
	ax.set_xlabel('Year')
	ax.set_xlim([1880, 2020])
	ax.set_ylabel('PFT spatial average (%)')
	ax.set_title('Lineplot: PFT spatial average at MC', pad=12)
	ax.legend(loc='upper right')

	# Save figure
	plt.tight_layout()
	Figure_Path = '../output/Output_Figure/Plot.Lineplot.PFT_SptAvg/'
	Figure_Name = 'Plot.Lineplot.PFT_SptAvg.{}.png'.format(Plot_Config['n_PFT'])
	if not (os.path.exists(Figure_Path)): os.makedirs(Figure_Path)
	plt.savefig(Figure_Path + Figure_Name)

	return

if (__name__ == '__main__'):

	# Get data
	Data_PFT_CTL    = xr.open_dataset('../src/LUCForcing/landuse.timeseries_0.9x1.25_CESM2_PGWDEF_CTL.nc')['PCT_NAT_PFT'].sel(time=slice('1880', '2020')).values
	Data_PFT_fixLUC = xr.open_dataset('../src/LUCForcing/landuse.timeseries_0.9x1.25_CESM2_PGWDEF_fixLUC.nc')['PCT_NAT_PFT'].sel(time=slice('1880', '2020')).values
	
	# Calculate spatial average
	Data_PFT_CTL    = tool_CDPrep.Calc_SpatialAverage(Data_PFT_CTL, LandMask='Land', RangeMask='MC_Analysis')
	Data_PFT_fixLUC = tool_CDPrep.Calc_SpatialAverage(Data_PFT_fixLUC, LandMask='Land', RangeMask='MC_Analysis')

	# Plot
	for i_PFT in range(0, 15):

		Plot_Data = {\
			'PFT_CTL'   : Data_PFT_CTL[:, i_PFT], \
			'PFT_fixLUC': Data_PFT_fixLUC[:, i_PFT], \
		}

		Plot_Config = {\
			'n_PFT'     : i_PFT, \
		}

		Plot_Lineplot(Plot_Data, Plot_Config)