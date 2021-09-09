#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


###############################################################################
# STORM TO STOPMOTION (swan cases)

def plot_polar_coords(xr, yr, xp, yp, res, #xlen, ylen, 
                      sreg=0.5, spol=2.5, creg='grey', cpol='r', 
                      ttl_domain='15º', ttl_res='15km'):
    '''
    plot custom polar coordinates, stopmotion methodology
    
    xr, yr      - regular coordinates
    xp, yp      - polar coordinates
#    xlen, ylen  - domain length (x,y)
    res         - regular resolution 
    '''
    
    fac = (xr.max()-xr.min()) / (yr.max()-yr.min())

    plt.figure(figsize=(10*fac, 10))
    
    plt.scatter(xr, yr, c=creg, s=sreg, label='Regular ({0})'.format(xr.size))
    plt.scatter(xp, yp, c=cpol, s=spol, label='Polar ({0})'.format(xp.size))
    
    plt.title('Domain {0}x{0} ({1}),  Grid: regular ({2}), polar({3})'.format(
            ttl_domain, ttl_res, xr.size, xp.size), fontweight='bold')
    
    plt.xlabel('Longitude (m)', fontweight='bold'); 
    plt.ylabel('Latitude (m)', fontweight='bold')
    
    plt.xlim([xr.min()-res, xr.max()+res]); 
    plt.ylim([yr.min()-res, yr.max()+res]); # plt.legend();
#    plt.axis('equal')
    
    return plt

def plot_grid_segments(st_list, df_seg, xlen1, xlen2, ylen1, ylen2, st=None, 
                       n_rows=4, n_cols=6, width=12, height=6.4,
                       cwarm='crimson', cseg='limegreen', ctrack='silver',
                       cfacecolor1='blue', cfacecolor2='yellow'):
                        #n_rows=5, n_cols=4, width=8, height=8):
    '''
    st_list         - list of swan stopmotion cases 
    df_seg          - (pandas.Dataframe) stopmotion segments parameters
    st              - (pandas.Dataframe) real, filled, interpolated track
    xlen,ylen       - swan extent in cartesian domain [m]
    
    returns:  figure 1, stopmotion configuration for SWAN cases
              figure 2, real storm track analogue segments
    '''
    
    # figure 1: SWAN cases
    fig1 = plt.figure(figsize=(width, height))
    gs = gridspec.GridSpec(n_rows, n_cols, wspace=0, hspace=0)
    gr, gc = 0, 0

    for i,sti in enumerate(st_list):
        
        # number of steps for 1 day time
        dt = int(sti.attrs['override_dtcomp'].split()[0])
        n_steps = int(24 * 60 / dt)

        ax = plt.subplot(gs[gr, gc])
        
        # stopmotion segment event
        ax.plot(sti['x'][:n_steps], sti['y'][:n_steps], c=cwarm, linewidth=2)
        ax.plot(sti['x'][n_steps:], sti['y'][n_steps:], c=cseg, linewidth=3); 
        
        ax.set_xticklabels([]);             ax.set_yticklabels([]);
        ax.set_xlim([-xlen1, xlen2]);       ax.set_ylim([-ylen1, ylen2]);
        spc = (xlen1+xlen2)/38
        ax.text(-xlen1+spc, -ylen1+spc, i, fontsize=11, fontweight='bold')
#        ax.set_xlim([-xlen/2, xlen/2]);     ax.set_ylim([-ylen/2, ylen/2]);
#        spc = xlen/38
#        ax.text(-xlen/2+spc, -ylen/2+spc, i, fontsize=11, fontweight='bold')
        ax.patch.set_facecolor(cfacecolor1)
        ax.patch.set_alpha(0.08)        
        
        # counter
        gc += 1
        if gc >= n_cols:
            gc = 0
            gr += 1

    # figure 2: real storm segments
    if isinstance(st, pd.DataFrame):

        fig2 = plt.figure(figsize=(width, height))
        gs = gridspec.GridSpec(n_rows, n_cols, wspace=0, hspace=0)
        gr, gc = 0, 0
        
        lo0 = df_seg['lon0'].values[:]     # warmup origin coordinate
        la0 = df_seg['lat0'].values[:]
        lo = df_seg['lon'].values[:]       # target origin coordinate
        la = df_seg['lat'].values[:]
    
        xlim1 = np.floor(st['lon'].values.min()) -1
        xlim2 = np.ceil(st['lon'].values.max()) +1
        ylim1 = np.floor(st['lat'].values.min())
        ylim2 = np.ceil(st['lat'].values.max())
    
        for i,sti in enumerate(st_list):
    
            ax = plt.subplot(gs[gr, gc])
            
            # track
            ax.plot(st['lon'], st['lat'], '-', c=ctrack);
            # note: first segment skipped (has no preceding warmup)
            ax.plot([lo0[i+1], lo[i+1]], [la0[i+1], la[i+1]], c=cwarm, linewidth=2)
            ax.plot([lo[i+1], lo[i+2]], [la[i+1], la[i+2]], c=cseg, linewidth=3)
    
            ax.set_xticklabels([]);         ax.set_yticklabels([]);
            ax.set_xlim([xlim1, xlim2]);    ax.set_ylim([ylim1, ylim2])
            ax.text(xlim1 +0.35, ylim1 +0.35, i, fontsize=11, fontweight='bold')
            ax.patch.set_facecolor(cfacecolor2)
            ax.patch.set_alpha(0.08)        
            
            # counter
            gc += 1
            if gc >= n_cols:
                gc = 0
                gr += 1
                
    else:   fig2 = None

    return fig1, fig2

###############################################################################
# PLOT DATABASE (shytcwaves library)

def axplot_scatter(axs, x, y, s=1, color='grey'):
    
    axs.scatter(x, y, marker=".", s=s, color=color)
    
    return axs

def axplot_kde(axs, x, y, s=0.5, cmap='plasma'):
    
    xy = np.vstack([x, y])
    z = gaussian_kde(xy)(xy)
    
    axs.scatter(x, y, s=s, c=z, cmap=cmap)

    return axs
    
def axplot_multi_scatter(axs, df, labels, units, color='grey', 
                          xmin=None, xmax=None, ymin=None, ymax=None, 
                          s=1, fontsize=12, cmap='plasma', mode_kde=False):
    '''
    multiple variables scatter plot
    '''
    # get variables
    var = df.keys()
    
    for i, item in enumerate(var):
        
        var_copy = list(var.values).copy()
        labels_copy = labels.copy()
        units_copy = units.copy()
                
        if xmin:# & xmax & ymin & ymax:
            xmin_copy, xmax_copy = xmin.copy(), xmax.copy()
            
        for p in range(i+1):
            var_copy.pop(0)
            labels_copy.pop(0)
            units_copy.pop(0)
            if xmin:# & xmax & ymin & ymax:
                xmin_copy.pop(0)
                xmax_copy.pop(0)
                        
        for j, jtem in enumerate(var_copy):
            
            # subplot
            if mode_kde:  axplot_kde(axs[i,j+i], 
                                     df[jtem], df[item], s=s, cmap=cmap)
            else:         axplot_scatter(axs[i,j+i], 
                                         df[jtem], df[item], s=s, color=color)
            
            # attributes
            if xmin:# & xmax & ymin & ymax:
                axs[i,j+i].set_xlim(xmin_copy[j], xmax_copy[j])
                axs[i,j+i].set_ylim(ymin[i], ymax[i])
            if i==j+i:
                axs[i,j+i].set_xlabel('{0} ({1})'.format(
                        labels_copy[j], units_copy[j]), fontweight='bold')
                if i==0:
                    axs[i,j+i].set_ylabel('{0} ({1})'.format(
                            labels[i], units[i]), fontweight='bold')
            else:
                axs[i,j+i].set_xticks([])
                axs[i,j+i].set_yticks([])
                
            # axis off
            len_var = len(var)-len(var_copy)-1
            if len_var != 0:
                for k in range(len_var):    axs[i,k].axis('off')
                    
    return axs

def aux_modeSI(dfplot, var, units):
    'Converts velocities and distances into international system units'

    for i,v in enumerate(var):
        if v in ['dwseg','wseg','dvseg','vseg']:
            units[i] = 'km/h'
            dfplot[var[i]] *= 1.852  # [kt to km/h]
        if v in ['drseg','rseg']:
            units[i] = 'km'
            dfplot[var[i]] *= 1.852  # [nmile to km]
    
    return dfplot, units

def plot_params_scatter(df, var=['daseg','dpseg','pseg','dwseg','wseg', \
                                 'dvseg','vseg','drseg','rseg','laseg'], 
                        labels=['dA','dP','P','dW','W','dV','V','dR','R', \
                                'absLat'],
                        units=['º','mbar','mbar','kt','kt','kt','kt', \
                               'nmile','nmile','º'], 
                        color='grey', cmap='plasma', s=0.5, ttl='', 
                        modeSI=False, mode_kde=False):
    '''
    df      - (pandas.Dataframe)
    var     - list of variables to plot
    label   - list of variable labels
    units   - list of unit labels
    
    (optional)
    color   - color code
    ttl     - title
    modeSI  - activate for international system units [km/h, km]
    
    Plots parameterized stopmotion units
    '''
    
    # get variables for plot
    dfplot = df[var].dropna()
    
    if modeSI:  dfplot, units = aux_modeSI(dfplot, var, units)
        
    # figure
    N = dfplot.keys().size - 1
    fig, axs = plt.subplots(N, N, figsize=(15, 15))
    plt.subplots_adjust(wspace=0, hspace=0)

    # scatter
    if not mode_kde:
        axplot_multi_scatter(axs, dfplot, labels, units, color=color, s=1); 
    else:
        axplot_multi_scatter(axs, dfplot, labels, units, s=s, cmap=cmap, 
                             mode_kde=True); 
    # title
    plt.suptitle(ttl +' ({0})'.format(dfplot.shape[0]), 
                 y=0.90, fontweight='bold')
        
    return fig
    
#def plot_scatter_all(df, var=['daseg','dpseg','pseg','dwseg','wseg', \
#                              'dvseg','vseg','drseg','rseg','laseg'], 
#                     labels=['dA','dP','P','dW','W','dV','V','dR','R', \
#                             'absLat'],
#                     units=['º','mbar','mbar','kt','kt','kt','kt','nmile', \
#                            'nmile','º'], 
#                     color='grey', ttl='', modeSI=False):
#    '''
#    df      - (pandas.Dataframe)
#    var     - list of variables to plot
#    label   - list of variable labels
#    units   - list of unit labels
#    
#    (optional)
#    color   - color code
#    ttl     - title
#    modeSI  - activate for international system units [km/h, km]
#    
#    Plots parameterized stopmotion units
#    '''
#    
#    # get variables for plot
#    dfplot = df[var].dropna()
#    
#    if modeSI:  dfplot, units = aux_modeSI(dfplot, var, units)
#        
#    # figure
#    N = dfplot.keys().size - 1
#    fig, axs = plt.subplots(N, N, figsize=(15, 15))
#    plt.subplots_adjust(wspace=0, hspace=0)
#
#    # scatter
#    axplot_multi_scatter(axs, dfplot, labels, units, color=color, s=1); 
#
#    plt.suptitle(ttl +' ({0})'.format(dfplot.shape[0]), 
#                 y=0.90, fontweight='bold')
#        
#    return fig
#    
#def plot_kde_all(df, var=['daseg','dpseg','pseg','dwseg','wseg', \
#                              'dvseg','vseg','drseg','rseg','laseg'], 
#                     labels=['dA','dP','P','dW','W','dV','V','dR','R', \
#                             'absLat'],
#                     units=['º','mbar','mbar','kt','kt','kt','kt','nmile', \
#                            'nmile','º'], 
#                     cmap='plasma', s=0.5, ttl='', modeSI=False):
#        
#    # get variables for plot
#    dfplot = df[var].dropna()
#    
#    if modeSI:  dfplot, units = aux_modeSI(dfplot, var, units)
#        
#    # figure
#    N = dfplot.keys().size - 1
#    fig, axs = plt.subplots(N, N, figsize=(15, 15))
#    plt.subplots_adjust(wspace=0, hspace=0)
#
#    # scatter
#    axplot_multi_scatter(axs, dfplot, labels, units, s=s, cmap=cmap, 
#                         mode_kde=True); 
#
#    plt.suptitle(ttl +' ({0})'.format(dfplot.shape[0]), 
#                 y=0.90, fontweight='bold')
#        
#    return fig

def plot_params_mda(df_dataset, df_subset, 
                    var=['daseg','dpseg','pseg','dwseg','wseg', \
                         'dvseg','vseg','drseg','rseg','laseg'], 
                    labels=['dA','dP','P','dW','W','dV','V','dR','R','absLat'],
                    units=['º','mbar','mbar','kt','kt','kt','kt','nmile', \
                           'nmile','º'],
                    c_dataset='grey', c_subset='purple', 
                    s=0.5, ttl='', modeSI=False):
    
    # get variables for plot
    dfplot_1 = df_dataset[var].dropna()
    dfplot_2 = df_subset[var].dropna()
    
    if modeSI:  
        dfplot_1, units = aux_modeSI(dfplot_1, var, units)
        dfplot_2, units = aux_modeSI(dfplot_2, var, units)
        
    # figure
    N = dfplot_1.keys().size - 1
    fig, axs = plt.subplots(N, N, figsize=(15, 15))
    plt.subplots_adjust(wspace=0, hspace=0)

    # scatter
    axplot_multi_scatter(axs, dfplot_1, labels, units, color=c_dataset, s=1); 
    axplot_multi_scatter(axs, dfplot_2, labels, units, color=c_subset, s=s); 

    plt.suptitle(ttl + 'Dataset ({0}), MDA ({1})'.format(
            dfplot_1.shape[0], dfplot_2.shape[0]
            ), y=0.90, fontweight='bold')
        
    return fig
    
def plot_scatter_basins(df, var=['basin','daseg','dpseg','pseg','dwseg','wseg',\
                                 'dvseg','vseg','drseg','rseg','laseg'], 
                        labels=['dA','dP','P','dW','W','dV','V','dR','R', \
                                'absLat'],
                        units=['º','mbar','mbar','kt','kt','kt','kt', \
                               'nmile','nmile','º'], 
                        basins=['NA','WP','SP','SI','EP','NI','SA'], 
                        ttl='', modeSI=False):
    
    # get variables for plot
    dfplot = df[var].dropna()
    
    if modeSI:  dfplot, units = aux_modeSI(dfplot, var, units)
        
    # color dictionary
    df_color = pd.DataFrame([['b','m','r','c','k','y','g']], 
                            columns=['NA','SA','WP','EP','SP','NI','SI'])
    # figure
    N = dfplot.keys().size - 2
    fig, axs = plt.subplots(N, N, figsize=(15, 15))
    plt.subplots_adjust(wspace=0, hspace=0)

    # loop basins
    for bas in basins:
        df_basin = dfplot.loc[dfplot['basin'] == bas]
        df_plot = df_basin.drop(columns=['basin'])
        # scatter
        axplot_scatter(axs, df_plot, 
                       labels, units, color=df_color[bas]); 

    plt.suptitle(ttl +' ({0})'.format(dfplot.shape[0]), 
                 y=0.90, fontweight='bold')

    return fig


###############################################################################
# ENSEMBLE (shytcwaves postprocessing)
    

