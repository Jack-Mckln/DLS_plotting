# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 11:25:31 2024

@author: jackm
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from pathlib import Path
import seaborn as sns


def read_RH(RH_file, delimiter = ','):
    RH_csv = np.loadtxt(RH_file, delimiter = delimiter) #importing the data from the CSV file
    RH_T = RH_csv.T # transposing the data to have all time and all RH in one list
    
    time_vals = RH_T[0]
    RH_vals = RH_T[1]
    
    return time_vals, RH_vals

RH_file = Path(r'C:\Users\jackm\OneDrive - University of Bath\PhD\Acoustic levtitation beamtimes\Diamond beamtime February 2024\Diamond RH data\Drop5_DPPC_POPC_vesicles_NaCl7-5mg-ml23_02_2024-01-57.txt')

RH_t, RH = read_RH(RH_file, delimiter = ',')

#making figure and subplots
fig, (ax1, ax2) = plt.subplots(1 , 2, sharey=True, width_ratios= [1,.2])

#RH = [[1,2,3],[4,5,6]]

ax1.plot(RH, RH_t)
ax2.plot(RH_t, RH)

plt.show()


# Generate sample data
data1 = np.random.rand(30, 30)  # Sample data for left plot (Seaborn heatmap)
x2 = np.linspace(0, 10, 1300)    # Sample x-values for right plot (line plot)
y2 = np.sin(x2)                  # Sample y-values for right plot

# Create figure and subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6), gridspec_kw={'width_ratios': [1, 5]})

# Plot left subplot (Seaborn heatmap)
sns.heatmap(data1, ax=ax1)
ax1.set_title('Seaborn Heatmap')

# Plot right subplot (line plot)
ax2.plot(x2, y2, color='blue')
ax2.set_title('Line Plot')

# Set y-limits to match the range of both datasets
y_min = min(data1.min(), y2.min())
y_max = max(data1.max(), y2.max())
ax1.set_ylim(y_min, y_max)
ax2.set_ylim(y_min, y_max)

# Show plot
plt.tight_layout()
plt.show()