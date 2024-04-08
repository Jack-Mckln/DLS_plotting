# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 16:08:41 2024

@author: jackm
"""
import h5py
import numpy as np 
import os
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib as mpl
import seaborn as sns
import glob
import math
from pathlib import Path
import re
import ciso8601 #for parsing time strings quickly
import datetime # for unix epoch conversion
import socket #
from matplotlib.colors import LogNorm, Normalize


### Determining which computer I am on to set the correct paths (not having to edit paths each time I change computer) ###

hostname = socket.gethostname()

if hostname == 'UDPC-1608-0104':
    user_path = Path(r'C:\Users\jpm93')

else:
    user_path = Path(r'C:\Users\jackm')


### Function to open and read the correct data ###
def read_nxs(filename, source, timing_path): 
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
            time = np.array(hdf.get('entry_0000/PyFAI/result_ave/t'))
            start_time = datetime.datetime.timestamp(
                ciso8601.parse_datetime(
                str(
                np.array(hdf.get('entry_0000/start_time')))
                .strip("b'")
                )
                )
            
        elif source == 'DLS':
            data = np.array(hdf.get('processed/result/data'))
            q_values = np.array(hdf.get('processed/result/q'))
            
            with h5py.File(timing_path) as time_hdf:
                
                #endtime from file
                end_time = np.array(time_hdf.get('entry1/end_time')) 
                #convert end time to string and strip superfluous characters
                end_time = str(end_time).strip("b'")
                #convert string into datetime object
                end_time = ciso8601.parse_datetime(end_time)
                #convert datetime object to seconds from unix epoch start to make usable
                end_time = datetime.datetime.timestamp(end_time)
                
                
                #same thing for starttime in one line
                start_time = datetime.datetime.timestamp(
                    ciso8601.parse_datetime(
                    str(
                    np.array(time_hdf.get('entry1/start_time')))
                    .strip("b'")
                    )
                    )
                
                
                    
                time_dif = end_time - start_time
                time_step = time_dif/len(data)
                time = np.array([time_step*x for x in range(len(data)+1)])


        else:
            raise Exception("Error: Source must be either 'DLS' or 'ESRF'")
 
        
        
        return (data, q_values, time, start_time)


### TREE PRINTING FUNCTION ###
def h5_tree(val, pre =''):
    items = len(val)
    for key, val in val.items():
        items -= 1
        if items == 0:
            #considering the last items
            if type(val) == h5py._hl.group.Group:
                print(pre + '¦_____' + key)
                h5_tree(val, pre + '     ')  #recursively through all groups
            try:
                print(pre + '¦----- ' + key + '(%d)' %len(val))
            #If the dataset is sacalar it doesn't have a length so this is to account for that error
            except:
                print(pre + '¦----- ' + key + '(scalar)')
        else:
            if type(val) == h5py._hl.group.Group:
                print(pre + '¦----- ' + key)
                h5_tree(val, pre + '¦      ')
            else:
                try:
                    print(pre + '¦----- ' + key + '(%d)' %len(val))
                #If the dataset is sacalar it doesn't have a length so this is to account for that error
                except:
                    print(pre + '¦----- ' + key + '(scalar)')



def plot_one_file_all_pos(filename, source, n_slices=1, step=1):
    data, q_values, t, _ = read_nxs(filename, source, timing)
    
    
    
    if source == 'DLS':
        if len(data.shape) == 3: #Check if there are multiple frames per scan
            data  = np.sum(data, axis=1)  #reducing the data by summing it along columns (i.e. along frames)      
        else:
            pass
        
        ### SETTING UP FIGURE ###
        fig_one_all = plt.figure()
        ax_one_all = fig_one_all.add_subplot(1,1,1)
        ax_one_all.set_yscale("log")
        ax_one_all.set_xlabel('q (1/nm)')
        ax_one_all.set_ylabel('Intensity (arb. units)')
        ax_one_all.set_yticklabels([])
        ax_one_all.set_prop_cycle(plt.cycler('color', 
                                           plt.cm.jet(np.linspace(0, 1, len(data))))) #setting the colour map to cycle through                                                
                                                                                             #can use different colour maps, here jet goes
        for x in range(0,len(data),step):
            plt.plot(q_values,(data[x]*10**x))
            
    elif source == 'ESRF':
        q_values = q_values/10
        
        slice_size=(len(data)/n_slices) # cut the data up into the vertical slices
        
        if slice_size.is_integer() == True:
            
        
            slice_start = 0
            while slice_start < (len(data)):
                
                ### SETTING UP FIGURE ###
                fig_one_all = plt.figure()
                ax_one_all = fig_one_all.add_subplot(1,1,1)
                ax_one_all.set_yscale("log")
                ax_one_all.set_xlabel('q (1/nm)')
                ax_one_all.set_ylabel('Intensity (arb. units)')
                ax_one_all.set_yticklabels([])
                ax_one_all.set_prop_cycle(plt.cycler('color', 
                                                   plt.cm.jet(np.linspace(0, 1, len(data))))) #setting the colour map to cycle through                                                
                                                                                                     #can use different colour maps, here jet goes
                for x in range(int(slice_start), int(slice_start+slice_size),step): 
                    ax_one_all.set_title(f'{x}')
    
                    plt.plot(q_values, data[x]*10**x,
                                            #s=0.1,
                                            label=f"{x}")
                
                plt.show(fig_one_all) #Unsure why this cant ber fig one all but then it doesnt work, and this means that I dont have the correct set up for scales etc
                plt.close(fig_one_all) # For some reason wihtout this it also prints the previous plot

                slice_start += slice_size
                
                                    
            
        
        else:
            raise Exception("Length of data / slice size is not an integer. Check your slice size")
        
    
    
    
    else:
        print("Error: Source must be either 'DLS' or 'ESRF'")
    
    return 
    

def plot_multi_file_one_pos_time(source, scatter_type, 
                                 position, n_pos,
                                 start_file_no = 0, end_file_no = 'end',
                                 range_min = 0, range_max = 2):
    
    #data, q_values, t, start_time = read_nxs(filename, source)
    t_vals = np.array([[]]) #array to be populated with the time values corresponding to each scan
    filter_data=[]
    
    
    
    
  
    if source == 'DLS' and scatter_type == 'SAXS':
       
        try:
            del(DLS_SAXS_file_list[0:start_file_no]) # remove the first files that we are not interested in (from 0 to start_file_no)
                                                     # Timing then start at the first file retained
            if end_file_no == 'end':
                pass
            else:
                del(DLS_SAXS_file_list[end_file_no:])
            
            for file in DLS_SAXS_file_list:
                timing_path = DLS_drop_dir / (file.split('_')[0]+'.nxs')
                data, q_values, t, start_time = read_nxs(DLS_drop_dir / file, source, timing_path)
                t_vals = np.append(t_vals, t+start_time)
                
                if source == 'DLS':
                    if len(data.shape) == 3: #Check if there are multiple frames per scan
                        data  = np.sum(data, axis=1) #reducing the data by summing it along columns (i.e. along frames)
                               
                    else:
                        pass
                
                filter_data.append(data[position])   # list of arratys each containing the data at the position requested summed over all frames  
       
        except TypeError:
            print(f"Incompatible file {file}, data missing")
            
            
        fig,ax = plt.subplots()
        for file in range(len(filter_data)):
            
            ax.plot(q_values, filter_data[file]*10**file) 
            ax.set_yscale('log')
            #need to change axes
                
        fig,ax = plt.subplots()
        c = ax.imshow(np.log(filter_data), vmin = range_min, vmax = range_max, aspect = 'auto')
        ax.set_xlabel('q (1/nm)')
        ax.set_ylabel('Time (s) ')
        plt.xticks(ticks = np.arange(175, len(q_values), (len(q_values)/5)), labels = ['{0:.1f}'.format(q_values[int(x)]) for x in list(np.arange(175, len(q_values), (len(q_values)/5)))])
        
        #### NEED TO FIX THE X TICKS SECTION HERE TO AMKE SURE IT CORRESPONDS TO THE Q_VALUES i WANT IT TO CORRESOND TO ###
        
        ax.invert_yaxis() #making sure that the longest time is plotted at the top
        #plt.colorbar(c)
        

        t_vals = t_vals-t_vals[0]
        t_vals = t_vals[position::n_pos] # take only the time correspoinding to a certain position
        
        
        ax_heat = sns.heatmap(np.log10(filter_data), 
                               vmin = range_min, vmax = range_max, 
                               xticklabels = 'auto', yticklabels = 10,
                               )
        ax_heat.set_xlabel('q (1/nm)')
        ax_heat.set_ylabel('Time (s) ')            
        ax_heat.invert_yaxis()
        plt.yticks(ticks = np.arange(0, len(t_vals), len(t_vals)/5), labels = ['{0:.0f}'.format(t_vals[int(x)]) for x in list(np.arange(0, len(t_vals), len(t_vals)/5))])
        plt.xticks(ticks = np.arange(175, len(q_values), (len(q_values)/5)), labels = ['{0:.1f}'.format(q_values[int(x)]) for x in list(np.arange(175, len(q_values), (len(q_values)/5)))])


        pass
    
    elif source == 'DLS' and scatter_type == 'WAXS':
        
        try:
            del(DLS_WAXS_file_list[0:start_file_no]) # remove the first files that we are not interested in (from 0 to start_file_no)
                                                     # Timing then start at the first file retained
                                                     
            if end_file_no == 'end': #remove final files we are not interested in
                pass
            else:
                del(DLS_WAXS_file_list[end_file_no:])
                
                
            for file in DLS_WAXS_file_list:
                timing_path = DLS_drop_dir / (file.split('_')[0]+'.nxs')
                data, q_values, t, start_time = read_nxs(DLS_drop_dir / file, source, timing_path)
                t_vals = np.append(t_vals, t+start_time)
                
                if source == 'DLS':
                    if len(data.shape) == 3: #Check if there are multiple frames per scan
                        data  = np.sum(data, axis=1) #reducing the data by summing it along columns (i.e. along frames)
                               
                    else:
                        pass
                
                filter_data.append(data[position])   # list of arratys each containing the data at the position requested summed over all frames  
        
        except TypeError:
            print(f"Incompatible file {file}, data missing")
        
        fig,ax = plt.subplots()
        for file in range(len(filter_data)):
            
            ax.plot(q_values, filter_data[file]*10**file) 
            ax.set_yscale('log')
            #need to change axes
                
        fig,ax = plt.subplots()
        c = ax.imshow(np.log(filter_data), vmin = range_min, vmax = range_max, aspect = 'auto')
        ax.set_xlabel('q (1/nm)')
        ax.set_ylabel('Time (s) ')
        plt.xticks(ticks = np.arange(175, len(q_values), (len(q_values)/5)), labels = ['{0:.1f}'.format(q_values[int(x)]) for x in list(np.arange(175, len(q_values), (len(q_values)/5)))])
        
        #### NEED TO FIX THE X TICKS SECTION HERE TO AMKE SURE IT CORRESPONDS TO THE Q_VALUES i WANT IT TO CORRESOND TO ###
        
        ax.invert_yaxis() #making sure that the longest time is plotted at the top
        #plt.colorbar(c)
        

        t_vals = t_vals-t_vals[0]
        t_vals = t_vals[position::n_pos] # take only the time correspoinding to a certain position
        
        
        ax_heat = sns.heatmap(np.log10(filter_data), 
                               vmin = range_min, vmax = range_max, 
                               xticklabels = 'auto', yticklabels = 10,
                               )
        ax_heat.set_xlabel('q (1/nm)')
        ax_heat.set_ylabel('Time (s) ')            
        ax_heat.invert_yaxis()
        plt.yticks(ticks = np.arange(0, len(t_vals), len(t_vals)/5), labels = ['{0:.0f}'.format(t_vals[int(x)]) for x in list(np.arange(0, len(t_vals), len(t_vals)/5))])
        plt.xticks(ticks = np.arange(175, len(q_values), (len(q_values)/5)), labels = ['{0:.1f}'.format(q_values[int(x)]) for x in list(np.arange(175, len(q_values), (len(q_values)/5)))])

        

        
        pass
    
    
    elif source == 'ESRF' and scatter_type == 'SAXS':
        
        try:
            del(ESRF_SAXS_file_list[0:start_file_no]) # remove the first files that we are not interested in (from 0 to start_file_no)
                                                     # Timing then start at the first file retained
            if end_file_no == 'end': # remove final files
                pass
            else:
                del(ESRF_SAXS_file_list[end_file_no:])
            
            for file in ESRF_SAXS_file_list:
                timing_path = DLS_drop_dir / (file.split('_')[0]+'.nxs')
                data, q_values, t, start_time = read_nxs(ESRF_drop_dir / file, source, timing_path)
                q_values = q_values/10
                t_vals = np.append(t_vals, t+start_time)
                
                
                filter_data.append(data[position])   # list of arratys each containing the data at the position requested summed over all frames  
           
        
        except TypeError:
           print(f"Incompatible file {file}, data missing")
           
           
        fig,ax = plt.subplots()
        
        for file in range(len(filter_data)):
            
            ax.plot(q_values, filter_data[file]*10**file) 
            ax.set_yscale('log')
            #data[x]*10**x
            #need to change axes
        
        
        
        
        fig,ax = plt.subplots()
        c = ax.imshow(np.log(filter_data), vmin = range_min, vmax = range_max, aspect = 'auto')
        ax.set_xlabel('q (1/nm)')
        ax.set_ylabel('Time (s) ')
        plt.xticks(ticks = np.arange(175, len(q_values), (len(q_values)/5)), labels = ['{0:.1f}'.format(q_values[int(x)]) for x in list(np.arange(175, len(q_values), (len(q_values)/5)))])
        
        #### NEED TO FIX THE X TICKS SECTION HERE TO AMKE SURE IT CORRESPONDS TO THE Q_VALUES i WANT IT TO CORRESOND TO ###
        
        ax.invert_yaxis() #making sure that the longest time is plotted at the top
        #plt.colorbar(c)
        
        t_vals = t_vals-t_vals[0]
        t_vals = t_vals[position::n_pos] # take only the time correspoinding to a certain position
        
        
        ax_heat = sns.heatmap(np.log10(filter_data), 
                               vmin = range_min, vmax = range_max, 
                               xticklabels = 'auto', yticklabels = 10,
                               )
        ax_heat.set_xlabel('q (1/nm)')
        ax_heat.set_ylabel('Time (s) ')            
        ax_heat.invert_yaxis()
        plt.yticks(ticks = np.arange(0, len(t_vals), len(t_vals)/5), labels = ['{0:.0f}'.format(t_vals[int(x)]) for x in list(np.arange(0, len(t_vals), len(t_vals)/5))])
        plt.xticks(ticks = np.arange(175, len(q_values), (len(q_values)/5)), labels = ['{0:.1f}'.format(q_values[int(x)]) for x in list(np.arange(175, len(q_values), (len(q_values)/5)))])
        
           
        
        pass
    
    elif source == 'ESRF' and scatter_type == 'WAXS':

        try:
            del(ESRF_WAXS_file_list[0:start_file_no]) # remove the first files that we are not interested in (from 0 to start_file_no)
                                                     # Timing then start at the first file retained
            if end_file_no == 'end': # remove final files
                pass
            else:
                del(ESRF_WAXS_file_list[end_file_no:])
            
            for file in ESRF_WAXS_file_list:
                timing_path = DLS_drop_dir / (file.split('_')[0]+'.nxs')
                data, q_values, t, start_time = read_nxs(ESRF_drop_dir / file, source, timing_path)
                q_values = q_values/10
                t_vals = np.append(t_vals, t+start_time)
    
    
                
                filter_data.append(data[position])   # list of arratys each containing the data at the position requested summed over all frames  
            
            
        except TypeError:
            print(f"Incompatible file {file}, data missing")
            
        fig,ax = plt.subplots()
        for file in range(len(filter_data)):
            
            ax.plot(q_values, filter_data[file]*10**file) 
            ax.set_yscale('log')
            #need to change axes
                
        fig,ax = plt.subplots()
        c = ax.imshow(np.log(filter_data), vmin = range_min, vmax = range_max, aspect = 'auto')
        ax.set_xlabel('q (1/nm)')
        ax.set_ylabel('Time (s) ')
        plt.xticks(ticks = np.arange(175, len(q_values), (len(q_values)/5)), labels = ['{0:.1f}'.format(q_values[int(x)]) for x in list(np.arange(175, len(q_values), (len(q_values)/5)))])
        
        #### NEED TO FIX THE X TICKS SECTION HERE TO AMKE SURE IT CORRESPONDS TO THE Q_VALUES i WANT IT TO CORRESOND TO ###
        
        ax.invert_yaxis() #making sure that the longest time is plotted at the top
        #plt.colorbar(c)
        
       # t_vals = t_vals-t_vals[0]
        t_vals = t_vals[position::n_pos] # take only the time correspoinding to a certain position
        
        
        ax_heat = sns.heatmap(np.log10(filter_data), 
                               vmin = range_min, vmax = range_max, 
                               xticklabels = 'auto', yticklabels = 10,
                               )
        ax_heat.set_xlabel('q (1/nm)')
        ax_heat.set_ylabel('Time (s) ')            
        ax_heat.invert_yaxis()
        plt.yticks(ticks = np.arange(0, len(t_vals), len(t_vals)/5), labels = ['{0:.0f}'.format(t_vals[int(x)]) for x in list(np.arange(0, len(t_vals), len(t_vals)/5))])
        plt.xticks(ticks = np.arange(175, len(q_values), (len(q_values)/5)), labels = ['{0:.1f}'.format(q_values[int(x)]) for x in list(np.arange(175, len(q_values), (len(q_values)/5)))])
            
    
    
    else: 
        raise Exception('Make sure source is one of "DLS" or "ESRF" and scatter type "WAXS" or "SAXS"')

        pass
    
    
    return (t_vals, filter_data, q_values)



    
def plot_one_file_one_pos(filename, source, position):
    
    data, q_values, t, _ = read_nxs(filename, source, timing)
    
    if source == 'DLS':
        if len(data.shape) == 3: #Check if there are multiple frames per scan
            data  = np.sum(data, axis=1)  #reducing the data by summing it along columns (i.e. along frames)      
        else:
            pass
        
        
                
    elif source == 'ESRF':
        q_values = q_values/10
        
    fig_one_all = plt.figure()
    ax_one_all = fig_one_all.add_subplot(1,1,1)
    ax_one_all.set_yscale("log")
    ax_one_all.set_xlabel('q (1/nm)')
    ax_one_all.set_ylabel('Intensity (arb. units)')
    ax_one_all.set_yticklabels([])
    ax_one_all.set_prop_cycle(plt.cycler('color', 
                                       plt.cm.jet(np.linspace(0, 1, len(data))))) #setting the colour map to cycle through                                                
       
                                                                                  #can use different colour maps, here jet goes
    plt.plot(q_values,data[position],
                #s=0.1
                )
    
    
    pass

### FUNCTION TO CREATE LISTS WITH THE REELVANT SAXS/WAXS FILES FOR THE MULTI TIME PLOT FUNCTIONS ###
def glob_re(pattern, strings):
    return list(filter(re.compile(pattern).match, strings)) #This isn't the best way to do this as it doesn't take advantage of lazy evaluation but it works and I dont really understand laziness yet


### Setting up the file paths ###

ESRF_mother_directory = user_path / 'OneDrive - University of Bath\PhD\Acoustic levtitation beamtimes\Data ESRF February 2024\X-Ray data' #specify the directory as a path object
DLS_mother_directory = user_path / 'OneDrive - University of Bath\PhD\Acoustic levtitation beamtimes\Diamond beamtime February 2024\Diamond X-ray data\processed'


ESRF_drop_dir = ESRF_mother_directory / 'Drop14_DPPC_5perc_2-pentanol'
DLS_drop_dir = DLS_mother_directory / 'Drop9_2-1_DPPC-POPC_7-5mgml_ves_80RH'
#background_name="i22-636979_saxs_Transmission_IvsQ_processed.nxs"


ESRF_scan_no = '00008' 
ESRF_SAXS_path = ESRF_drop_dir / f'drop14_DPPC_5perc_pentanol_eiger2_{ESRF_scan_no}_00_ave.h5'
ESRF_WAXS_path = ESRF_drop_dir / f'db02_waxs_{ESRF_scan_no}_00_ave.h5'



DLS_scan_no = '736088'
DLS_SAXS_path = DLS_drop_dir / f"i22-{DLS_scan_no}_saxs_Transmission_IvsQ_processed.nxs"
DLS_WAXS_path = DLS_drop_dir / f"i22-{DLS_scan_no}_waxs_Transmission_IvsQ_processed.nxs"


timing= DLS_drop_dir / f'i22-{DLS_scan_no}.nxs'


#background_path = droplet_directory / background_name # these two notations are equivalent

ESRF_SAXS_file_list=glob_re(r'.*eiger2_\d+.*ave\.h5', os.listdir(ESRF_drop_dir)) #this populates a list with al the names of files 
                                                                                #in the directory that end like this
                                                                                # and sorts them
ESRF_WAXS_file_list=glob_re(r'.*waxs_\d+.*ave\.h5', os.listdir(ESRF_drop_dir))


DLS_SAXS_file_list = glob_re(r'.*saxs_Transmission_IvsQ_processed.nxs', os.listdir(DLS_drop_dir))
DLS_WAXS_file_list = glob_re(r'.*waxs_Transmission_IvsQ_processed.nxs', os.listdir(DLS_drop_dir))



### FUNCTION CALLS ###

#t_vals = plot_multi_file_one_pos_time(DLS_SAXS_path, 'DLS', 5, 5, 5)

t_vals, filter_data, q_values = plot_multi_file_one_pos_time('DLS', 'SAXS', 
                                                             position = 9, n_pos = 21, 
                                                             start_file_no = 5,
                                                             range_min = 0 ,
                                                             range_max =1
                                                             )

#timing = DLS_drop_dir / DLS_SAXS_file_list[0]
#plot_one_file_all_pos(DLS_drop_dir / DLS_SAXS_file_list[0], 'DLS', step = 1)

#plot_one_file_one_pos(ESRF_SAXS_path, 'ESRF', 5)

#t_vals, t = plot_multi_file_one_pos_time(DLS_SAXS_path, 'DLS', 1, 1, 1)

#data, q_values, t, start_time= read_nxs(DLS_drop_dir / DLS_SAXS_file_list[0], 'ESRF', timing)
