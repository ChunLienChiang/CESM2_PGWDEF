"""
Plot_Lineplot_TS_SptAvg.py
==========================
The script to output the figure of the time series of global TS spatial average.
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

    Figure_Path = '../output/Output_Figure/Plot_Lineplot_TS_SpaAvg/'
    Figure_Name = 'Plot_Lineplot_TS_SpaAvg.png'

    # Create figure object
    fig, ax = plt.subplots(figsize=(5, 5), dpi=150)

    # Plot
    x_array = np.arange(Plot_Data['TS_Global_SptAvg'].size) + 1
    ax.plot(x_array, Plot_Data['TS_Global_SptAvg'], color='Black')
    ax.plot(x_array, Plot_Data['TS_NHL_SptAvg'], color='Red')
    ax.plot(x_array, Plot_Data['TS_SHL_SptAvg'], color='Blue')
    ax.plot(x_array, Plot_Data['TS_NHO_SptAvg'], color='Red', linestyle='dashed')
    ax.plot(x_array, Plot_Data['TS_SHO_SptAvg'], color='Blue', linestyle='dashed')

    # Configuration
    ax.set_xlabel('Year')
    ax.set_ylabel('TS (K)')
    ax.set_title('Lineplot: TS spatial average', pad=12)
    ax.text(0.99, 0, dt.datetime.today().strftime('%Y/%m/%d %H:%M'), ha='right', va='bottom', fontsize=7, transform=ax.transAxes)

    # Save figure
    if not (os.path.exists(Figure_Path)): os.makedirs(Figure_Path)
    plt.savefig(Figure_Path + Figure_Name)

    return

if (__name__ == '__main__'):

    # Get cliamte raw data: TS
    Data_TS = tool_GetCD.Get('CESM2_PGWDEF_CTL_e0', 'TS')

    # Calculate annual mean
    Data_TS = tool_CDPrep.Calc_AnnualMean(Data_TS)
    
    # Plot
    Plot_Data = {\
        'TS_Global_SptAvg': tool_CDPrep.Calc_SpatialAverage(Data_TS), \
        'TS_NHL_SptAvg': tool_CDPrep.Calc_SpatialAverage(Data_TS, LandMask='Land', RangeMask='NH_Analysis'), \
        'TS_SHL_SptAvg': tool_CDPrep.Calc_SpatialAverage(Data_TS, LandMask='Land', RangeMask='SH_Analysis'), \
        'TS_NHO_SptAvg': tool_CDPrep.Calc_SpatialAverage(Data_TS, LandMask='Ocean', RangeMask='NH_Analysis'), \
        'TS_SHO_SptAvg': tool_CDPrep.Calc_SpatialAverage(Data_TS, LandMask='Ocean', RangeMask='SH_Analysis'), \
    }
    Plot_Config = {\
        \
    }
    Plot_Lineplot(Plot_Data, Plot_Config)