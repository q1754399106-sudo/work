#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
This is the Python code for reading the catalog and the raw DAS spectra and plotting them.

Hilary Chang 2026.3.13

'''
# =============================================================================
iev=1466 # Select an target event (iev=1466 is Figure 3b)

# Please point the directory to this supplementary folder
opda_dir='/mnt/data0/hilarych/FORGE_2024/save/open_data_JGR' 

# =============================================================================

import pandas as pd
from obspy.core import UTCDateTime
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# Read catalog
infile=os.path.join(opda_dir,'event_catalog.csv')
cat=pd.read_csv(infile)             
ev=cat.iloc[iev]                    # Target event catalog
t0=UTCDateTime(ev['t0'])                        # Origin time
evid=t0.strftime('%Y%m%d%H%M%S_%f') # evid is formatted using t0
infile=os.path.join(opda_dir,'spectra','spec_D-1OB_P_arr_%s.npz'%evid)
r=np.load(infile,allow_pickle=True)                   # Read spectra file
ch_list_partial=r['ch_list_partial'] # Selected channels with good picks
spec_list_partial=r['spec_list_partial']   # P-arrival of raw DAS spectra of selected channels
spec_n_list_partial=r['spec_n_list_partial'] # Noise of raw DAS spectra (0.3 s before the signal window) of selected channels 
f_series=r['f_series'] # Frequency
num_all_ch=1600 # Number of channels in D-1B
fst=1000 # Sampling rate (Hz)
twin=0.1 # Time window for calculating the spectra

def pad_data_2_correct_channel_number(data_partial,ch_partial,num_all_ch):
    '''
    The spectra files only contain data of selected channels.
    This function pads nans to non-selected channels for the loaded data so that 
    the data sizes are the same for all events: (number of DAS channels in D-1OB 
                                                 (1600), number of time samples)
    
    '''
    # Initialize
    data_full=np.full([num_all_ch, data_partial.shape[1]], np.nan)
    data_full[ch_partial]=data_partial
    return data_full

spec_list=pad_data_2_correct_channel_number(spec_list_partial,ch_list_partial,num_all_ch)
spec_n_list=pad_data_2_correct_channel_number(spec_n_list_partial,ch_list_partial,num_all_ch)  

x=f_series
y1=spec_list
y2=spec_n_list

colormap = plt.cm.viridis
colorst = [colormap(i/len(y1)) for i in range(len(y1))]

fig,ax=plt.subplots(figsize=[6,4])
ax.set_facecolor('gray')
# Initialize
spec_all_stack=[]
for ist in range(len(y1)):
    # ch=mid_ch_list[ist]
    
    indsort=np.argsort(x)
    ax.semilogy(x[indsort],
                y1[ist][indsort],color=colorst[ist],lw=0.5,zorder=5)

    ax.semilogy(x[indsort],
                y2[ist][indsort],color='k',lw=0.2,zorder=4)        

ax.set_xlabel('Frequency (Hz)')
ax.set_ylabel('Strain rate amplitude spectra\n(radian per gauge length (10 m) per second)')
ax.set_xscale('log')
ax.set_yscale('log')
ax.grid(lw=0.5)
ax.set_xlim(1/twin,fst/2)
plt.title('D-1OB raw DAS spectra\nfor event t0:\n%s'%(t0))
norm = mpl.colors.Normalize(vmin=0, vmax=len(colorst))
cmap=plt.get_cmap("viridis", len(colorst))
sm = plt.cm.ScalarMappable(norm=norm, cmap=colormap)
cbar=fig.colorbar(sm,shrink=0.8,aspect=12,ax=ax)
cbar.ax.invert_yaxis()
cbar.ax.set_ylabel('Channel')
cbar.ax.tick_params(labelsize=9)

plt.show()
