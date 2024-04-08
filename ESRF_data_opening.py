# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 16:08:41 2024

@author: jackm
"""
import h5py
import numpy as np 
import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib as mpl
import seaborn as sns
import glob
import math
from pathlib import Path
import re


### Setting up the file paths ###

mother_directory=Path(r"C:\Users\jackm\OneDrive - University of Bath\PhD\Acoustic levtitation beamtimes\Data ESRF February 2024\X-Ray data") #specify the directory as a path object
DLS_mother_directory = Path(r"C:\Users\jackm\OneDrive - University of Bath\PhD\Acoustic levtitation beamtimes\Diamond beamtime February 2024\Diamond X-ray data\processed")


droplet_directory = mother_directory / 'Drop16_POPC_vesicles_H2O_0-75perc_Hamilton'
DLS_drop_dir = DLS_mother_directory / 'Drop16_2-1_DPPC-POPC_7-5mgml_NaCl7-5mgml_ves_RHramp_3'
#background_name="i22-636979_saxs_Transmission_IvsQ_processed.nxs"

scan_no = '00001' 
DLS_scan_no = '736088'

#DLS_SAXS_str = f'i22-{}_saxs_Transmission_IvsQ_processed'
#DLS_WAXS_str  = f'i22-{}_waxs_Transmission_IvsQ_processed'
SAXS_file = f"i22-{DLS_scan_no}_saxs_Transmission_IvsQ_processed.nxs"
WAXS_file = f"drop16_POPC_veiscles_H2O_waxs_{scan_no}_00_ave.h5"

SAXS_path = DLS_drop_dir / SAXS_file     #joining the file name to the directory path
WAXS_path = droplet_directory / WAXS_file
#background_path = droplet_directory / background_name # these two notations are equivalent

file_list=sorted(droplet_directory.glob("*ave.h5")) #this populates a list with al the names of files 
                                                                                #in the directory that end like this
                                                                                # and sorts them




### Function to open and read the correct data ###
def read_nxs(filename, source): 
    '''
    This will return a tuple in the form: (q_values, data).
    
    The data will have a shape dependent on if multiple SAXS patterns are stored 
    in the .nxs file. 
    
    For example, data from different positions along a capillary will be in the 
    form [x,y,z]:
        x = frame number (multiple frames if .nxs is a scan across a capillary/droplet etc.)
        y = intensity 
        z = q axis with a length of the number of q values 
    
    if I want to access the first frame of a scan, I would type:
        data[0,:,:]
        i.e. all of the intensities (y) at every q (z) for the first (0th) frame 
    
    check the shape of your data by using data.shape
    
    If it's just a single SAXS image, "data" will just be the intensity value at each q value
    
    '''
    
    with h5py.File(filename,'r') as hdf:
       # base_items = list(hdf.items()) # See what base items there are in the nxs
       # print('items in base directory:', base_items
        if source == 'ESRF':
            data = np.array(hdf.get('entry_0000/PyFAI/result_ave/data'))
            q_values = np.array(hdf.get('entry_0000/PyFAI/result_ave/q'))
            
        elif source == 'DLS':
            data = np.array(hdf.get('processed/result/data'))
            q_values = np.array(hdf.get('processed/result/q'))
            
        else:
            print("Error: Source must be either 'DLS' or 'ESRF'")
        #pro = hdf.get('processed') # Access the processed group
        #pro_items = list(pro.items())
        #print('Items in Pro:', pro_items)
    
        #res = pro.get('/processed/result') # Access the 'result' subgroup using its path in the file
        #res_items = list(res.items())
        #print(res_items)
    
        #data = np.array(res.get('data')) # Make a numpy array out of the 'data' part of the results subgroup
        #q_values = np.array(res.get('q')) # Make array out of q values
        
        
        return (data, q_values)


data, q = read_nxs(WAXS_path, 'ESRF')

data1, q1 = read_nxs(SAXS_path, 'DLS')



