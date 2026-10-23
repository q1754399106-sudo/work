#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
This is the Python code for reading and plotting the Q results and other 
measurements of the well logs.

Hilary Chang 2026.3.11

'''
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# Read the Q regults and well logs
# =============================================================================
# Please point the directory to this supplementary folder
opda_dir='/mnt/data0/hilarych/FORGE_2024/save/open_data_JGR'


# =============================================================================
infile=os.path.join(opda_dir,'Qp_profile_20-200Hz.npz')
r=np.load(infile,allow_pickle=True) 
med_Q_chs=r['med_Q_chs']
mad_Q_chs=r['mad_Q_chs']
dep_chs=r['dep_chs']
dep_interp=r['dep_interp']
Q_interp=r['Q_interp']
med_wlogRho_chs=r['med_wlogRho_chs']
mad_wlogRho_chs=r['mad_wlogRho_chs']
med_wlogVp_chs=r['med_wlogVp_chs']
mad_wlogVp_chs=r['mad_wlogVp_chs']
med_wlogG_chs=r['med_wlogG_chs']
mad_wlogG_chs=r['mad_wlogG_chs']


#% (Plot) Q vs density vs Vp vs shear modulus
fig=plt.figure(figsize=[9,4])
gs=GridSpec(1, 4,wspace=0.1,width_ratios=[1,0.5,0.5,0.5])
ylabel='Depth (m)'

# Initialize
axs=[]
axs.append(fig.add_subplot(gs[0]))

ax=axs[-1]
xlabel='$Q_p$'
x=med_Q_chs
xerr=mad_Q_chs
y=dep_chs
ax.scatter(x,y,zorder=5,alpha=0.5,edgecolor='none',lw=0.1,color='tab:blue',
           s=20,label='Measured')
ax.errorbar(x,y,xerr=xerr,capsize=1,lw=0,ms=0,capthick=0.3,color='gray',
            elinewidth=0.2,alpha=0.5)

x1=Q_interp
y1=dep_interp

ax.scatter(x1,y1,zorder=3,edgecolor='k',lw=0.2,color='none',
           label='Interpolated/\nExtrapolated',s=10)
ax.legend(fontsize=8)

ax.set_xscale('log')
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
ax.set_xlim([5,200])

axs.append(fig.add_subplot(gs[1]))
ax=axs[-1]
xlabel=r'$\rho$ ($g/cm^3$)'
x=med_wlogRho_chs
xerr=mad_wlogRho_chs
ax.scatter(x,y,zorder=5,alpha=0.5,edgecolor='none',lw=0.1,
           color='tab:blue',s=20)
ax.errorbar(x,y, xerr=xerr,capsize=1,lw=0,ms=0,capthick=0.3,color='gray',
            elinewidth=0.2)
ax.set_xlabel(xlabel)


axs.append(fig.add_subplot(gs[2]))
ax=axs[-1]
xlabel=r'$V_p$ ($km/s$)'
x=med_wlogVp_chs/1000
xerr=mad_wlogVp_chs/1000
ax.scatter(x,y,zorder=5,alpha=0.5,edgecolor='none',lw=0.1,
           color='tab:blue',s=20)
ax.errorbar(x,y, xerr=xerr,capsize=1,lw=0,ms=0,capthick=0.3,color='gray',
            elinewidth=0.2)
ax.set_xlabel(xlabel)

axs.append(fig.add_subplot(gs[3]))
ax=axs[-1]
xlabel='$\mu$ (GPa)'
x=med_wlogG_chs/1e9
xerr=mad_wlogG_chs/1e9
ax.scatter(x,y,zorder=5,alpha=0.5,edgecolor='none',lw=0.1,
           color='tab:blue',s=20)
ax.errorbar(x,y, xerr=xerr,capsize=1,lw=0,ms=0,capthick=0.3,color='gray',
            elinewidth=0.2)
ax.set_xlabel(xlabel)

ylim=[0,3000]
for i,ax in enumerate(axs):
    ax.grid(lw=0.5)
    ax.set_ylim(ylim)
    ax.invert_yaxis()
    if i!=0:
        ax.set_yticklabels([])

plt.show()