"""
Plot_Lineplot_TS_SptAvg.py
==========================
The script to output the figure of the time series of global TS spatial average.
"""

import numpy as np
import matplotlib.pyplot as plt
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
    plt.plot(np.arange(Plot_Data['TS_SptAvg'].size), Plot_Data['TS_SptAvg'])

    # Configuration

    # Save figure
    if not (os.path.exists(Figure_Path)): os.makedirs(Figure_Path)
    plt.savefig(Figure_Path + Figure_Name)

    return

if (__name__ == '__main__'):

    # Get cliamte raw data: TS
    Data_TS = tool_GetCD.Get('CESM2_PGWDEF_CTL_e0', 'TS')
    
    # Calculate global spatial average
    Data_TS = tool_CDPrep.Calc_SpatialAverage(Data_TS)

    # Plot
    Plot_Data = {\
        'TS_SptAvg': Data_TS\
    }
    Plot_Config = {\
        \
    }
    Plot_Lineplot(Plot_Data, Plot_Config)