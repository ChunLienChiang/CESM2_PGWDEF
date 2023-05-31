"""
Plot.Lineplot.TS_SptAvg.py
==========================
Output the figure of the time series of TS spatial average.
"""

import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import os
import sys
sys.path.append('../')
from tools import tool_GetClimData as tool_GetCD
from tools import tool_ClimData_Preprocessing as tool_CDPrep

def Plot_Lineplot(Plot_Data, Plot_Config):

	# Create figure object
	fig, ax = plt.subplots(figsize=(5, 5), dpi=300)

	# Plot
	x_array = np.arange(1880, 2020)
	ax.plot(x_array, Plot_Data['SST_CTL'], color='Black', zorder=2, label='CTL')
	ax.plot(x_array, Plot_Data['SST_fixSST'], color='Red', zorder=1, label='fixSST')

	# Configuration
	ax.set_xlabel('Year')
	ax.set_ylabel('TS (K)')
	ax.set_title('Lineplot: SST spatial average at IO&WP', pad=12)
	ax.legend(loc='upper right')

	# Save figure
	Figure_Path = '../output/Output_Figure/Plot.Lineplot.TS_SptAvg/'
	Figure_Name = 'Plot.Lineplot.TS_SptAvg.png'
	if not (os.path.exists(Figure_Path)): os.makedirs(Figure_Path)
	plt.savefig(Figure_Path + Figure_Name)

	return

if (__name__ == '__main__'):

	# Get data
	Data_SST_CTL    = tool_GetCD.Get_SST_Pacemaker('CTL', Grid='0.9x1.25')
	Data_SST_fixSST = tool_GetCD.Get_SST_Pacemaker('fixSST', Grid='0.9x1.25')
	
	# Calculate annual mean
	Data_SST_CTL    = tool_CDPrep.Calc_AnnualMean(Data_SST_CTL[1:-1, ...])
	Data_SST_fixSST = tool_CDPrep.Calc_AnnualMean(Data_SST_fixSST[1:-1, ...])

	# Calculate spatial average
	Data_SST_CTL    = tool_CDPrep.Calc_SpatialAverage(Data_SST_CTL, LandMask='Ocean', RangeMask='IOWP_Analysis')
	Data_SST_fixSST = tool_CDPrep.Calc_SpatialAverage(Data_SST_fixSST, LandMask='Ocean', RangeMask='IOWP_Analysis')

	# Plot
	Plot_Data = {\
		'SST_CTL'   : Data_SST_CTL, \
		'SST_fixSST': Data_SST_fixSST, \
	}

	Plot_Config = {\
		\
	}

	Plot_Lineplot(Plot_Data, Plot_Config)