# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 15:09:28 2024

@author: jackm
"""
import numpy as np
import h5py
import numpy as np 
import os
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import re
import ciso8601 #for parsing time strings quickly
import datetime # for unix epoch conversion
import socket #


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
