#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os.path as op
import numpy as np
import pandas as pd
from scipy.ndimage.filters import gaussian_filter
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

def plot_grid_segments(st_list, df_seg, xlen1, xlen2, ylen1, ylen2, 
                       st=None, st_list_0=None, df_seg_0=None,
                       n_rows=4, n_cols=6, width=12, height=6.4,
                       cwarm='crimson', cseg='limegreen', ctrack='silver',
                       cfacecolor1='blue', cfacecolor2='yellow'):
                        #n_rows=5, n_cols=4, width=8, height=8):
    '''
    st_list         - list of swan stopmotion cases (analogue/real)
    df_seg          - (pandas.Dataframe) stopmotion segments parameters
    xlen,ylen       - swan extent in cartesian domain [m]
    st              - (pandas.Dataframe) real, filled, interpolated track
    
    optional:
    st_list_0       - list of swan stopmotion cases (real)
    df_seg_0        - (pandas.Dataframe) stopmotion segments parameters (real)
    
    returns:  figure 1, stopmotion configuration for SWAN cases
              figure 2, real storm track analogue segments
    '''
    
    # remove NaNs
    if isinstance(st, pd.DataFrame) and df_seg.isna().any(axis=1).any(): 
        st = st[~df_seg.isna().any(axis=1)]
    df_seg = df_seg[~df_seg.isna().any(axis=1)]
    if isinstance(df_seg_0, pd.DataFrame):
        df_seg_0 = df_seg_0[~df_seg_0.isna().any(axis=1)]
    
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
        ax.plot(sti['x'][:n_steps], sti['y'][:n_steps], c=cwarm, linewidth=3)
        ax.plot(sti['x'][n_steps:], sti['y'][n_steps:], c=cseg, linewidth=3.5); 

        # stopmotion segment event (REAL)
        if st_list_0:
            sti0 = st_list_0[i]
            ax.plot(sti0['x'][:n_steps], sti0['y'][:n_steps], c='k', 
                    linewidth=.75, linestyle='dashed')
            ax.plot(sti0['x'][n_steps:], sti0['y'][n_steps:], c='k', 
                    linewidth=.75, linestyle='dashed'); 

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
        
    # subtitle
    if isinstance(st, pd.DataFrame):
        sign = st.lat.values[0]    # hemisphere sign
        if sign > 0:    ttl = '(North Hemisphere)'
        elif sign < 0:  ttl = '(South Hemisphere)'
        plt.suptitle('segments, ' + ttl, fontweight='bold', y=0.95)  # CHECK SPACE VERTICAL

    # figure 2: real storm segments
    if isinstance(st, pd.DataFrame):

        fig2 = plt.figure(figsize=(width, height))
        gs = gridspec.GridSpec(n_rows, n_cols, wspace=0, hspace=0)
        gr, gc = 0, 0
        
        lo0 = df_seg['lon_w'].values[:]     # warmup origin coordinate
        la0 = df_seg['lat_w'].values[:]
        lo1 = df_seg['lon_i'].values[:]     # target origin coordinate
        la1 = df_seg['lat_i'].values[:]
        lo2 = df_seg['lon_t'].values[:]     # target endpoint coordinate
        la2 = df_seg['lat_t'].values[:]
#        lo = df_seg['lon'].values[:]       # target origin coordinate
#        la = df_seg['lat'].values[:]
        if isinstance(df_seg_0, pd.DataFrame):
            lo0_0 = df_seg_0['lon_w'].values[:]     # warmup origin coordinate
            la0_0 = df_seg_0['lat_w'].values[:]
            lo1_0 = df_seg_0['lon_i'].values[:]     # target origin coordinate
            la1_0 = df_seg_0['lat_i'].values[:]
            lo2_0 = df_seg_0['lon_t'].values[:]     # target endpoint coordinate
            la2_0 = df_seg_0['lat_t'].values[:]
    
        xlim1 = np.floor(st['lon'].values.min()) -1
        xlim2 = np.ceil(st['lon'].values.max()) +1
        ylim1 = np.floor(st['lat'].values.min()) -1
        ylim2 = np.ceil(st['lat'].values.max()) +1
    
        for i,sti in enumerate(st_list):
    
            ax = plt.subplot(gs[gr, gc])
            
            # track
            ax.plot(st['lon'], st['lat'], '-', c=ctrack);
            
            # note: first segment skipped (has no preceding warmup)
            ax.plot([lo0[i], lo1[i]], [la0[i], la1[i]], c=cwarm, linewidth=3)
            ax.plot([lo1[i], lo2[i]], [la1[i], la2[i]], c=cseg, linewidth=3.5)
#            ax.plot([lo0[i+1], lo[i+1]], [la0[i+1], la[i+1]], c=cwarm, linewidth=2)
#            ax.plot([lo[i+1], lo[i+2]], [la[i+1], la[i+2]], c=cseg, linewidth=3)
    
            if isinstance(df_seg_0, pd.DataFrame):
                ax.plot([lo0_0[i], lo1_0[i]], [la0_0[i], la1_0[i]], c='k', 
                        linewidth=1, linestyle='dashed')
                ax.plot([lo1_0[i], lo2_0[i]], [la1_0[i], la2_0[i]], c='k', 
                        linewidth=1, linestyle='dashed')
                
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

def axplot_analogue(axs, x, y, s=10, marker='.', color=''):
    
    axs.scatter(x, y, marker=marker, s=s, color=color)
    
    return axs

def axplot_multi_scatter(axs, df, labels, units, color='grey', 
                          xmin=None, xmax=None, ymin=None, ymax=None, 
                          s=1, fontsize=12, cmap='plasma', marker='.',
                          mode_kde=False, mode_analogue=False):
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
            elif mode_analogue:
                          axplot_analogue(axs[i,j+i], 
                                          df[jtem], df[item], s=s, marker=marker)
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
    
def plot_params_mda(df_dataset, df_subset, 
                    var=['daseg','dpseg','pseg','dwseg','wseg', \
                         'dvseg','vseg','drseg','rseg','laseg'], 
                    labels=['dA','dP','P','dW','W','dV','V','dR','R','absLat'],
                    units=['º','mbar','mbar','kt','kt','kt','kt','nmile', \
                           'nmile','º'],
                    c_dataset='grey', c_subset='purple', 
                    s=0.5, ttl1='', ttl2='', modeSI=False):
    
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

    plt.suptitle('{0} ({1}), {2} ({3})'.format(ttl1, dfplot_1.shape[0], 
                 ttl2, dfplot_2.shape[0]
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
# ANALOGUES (shytcwaves initialization)

def plot_params_analogue(df_dataset, df_seg, df_analogue,
                        var=['daseg','dpseg','pseg','dwseg','wseg', \
                             'dvseg','vseg','drseg','rseg','laseg'], 
                        labels=['dA','dP','P','dW','W','dV','V','dR', \
                                'R','absLat'],
                        units=['º','mbar','mbar','kt','kt','kt','kt', \
                               'nmile','nmile','º'],
                        c_dataset='thistle', c_seg='slateblue', 
                        c_analogue='crimson', s=0.5, ttl1='', ttl2='', 
                        modeSI=False, mode_2vars=False):
    
    # get variables for plot
    dfplot_1 = df_dataset[var].dropna()
    dfplot_2 = df_seg[var].dropna()
    dfplot_3 = df_analogue[var].dropna()
    
    if modeSI:  
        dfplot_1, units = aux_modeSI(dfplot_1, var, units)
        dfplot_2, units = aux_modeSI(dfplot_2, var, units)
        dfplot_3, units = aux_modeSI(dfplot_3, var, units)
        
    # figure
    N = dfplot_1.keys().size - 1
    fig, axs = plt.subplots(N, N, figsize=(15, 15))
    plt.subplots_adjust(wspace=0, hspace=0)

    # scatter
    if mode_2vars:
        axplot_multi_scatter_2vars(axs, dfplot_2, dfplot_3, labels, units, s=s)
    else:
        axplot_multi_scatter(axs, dfplot_1, labels, units, color=c_dataset, s=1); 
        axplot_multi_scatter(axs, dfplot_2, labels, units, s=s+4, color=c_seg)
        axplot_multi_scatter(axs, dfplot_3, labels, units, s=s, color=c_analogue)

    plt.suptitle('{0} ({1}), {2} ({3})'.format(ttl1, dfplot_1.shape[0], 
                 ttl2, dfplot_2.shape[0]
                 ), y=0.90, fontweight='bold')
        
    return fig

def axplot_multi_scatter_2vars(axs, df1, df2, labels, units, 
                              xmin=None, xmax=None, ymin=None, ymax=None, 
                              s=1, fontsize=12, c1='slateblue', c2='red'):
    '''
    multiple variables scatter plot
    '''
    # get variables
    var = df1.keys()    # both dataframes have the same vars

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
            axplot_analogue(axs[i,j+i], df1[jtem], df1[item], # square (1)
                            s=s, marker='s', color='silver')
            axplot_analogue(axs[i,j+i], df1[jtem], df1[item], # dot (1)
                            s=s, marker='.', color=c1)
            axplot_analogue(axs[i,j+i], df2[jtem], df2[item], # dot (2)
                            s=s, marker='.', color=c2)
            axs[i,j+i].plot([df1[jtem], df2[jtem]], [df1[item], df2[item]], 
                            '--', c='silver')                 # line (1-2)
            
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

def plot_params_anom_series(df_seg, df_analogue, 
                            var=['laseg','daseg','dpseg','pseg','dwseg', \
                                 'wseg','dvseg','vseg','drseg','rseg',],
                            labels=['absLat','dA','dP','P','dW','W','dV','V',\
                                    'dR','R'],
                            units=['º','º','mbar','mbar','kt','kt','kt','kt',\
                                   'nmile', 'nmile'],
                            n_rows=10, n_cols=2, width=25, height=10):
    
    # remove NaN
    df_seg = df_seg[~df_seg.isna().any(axis=1)]
    df_analogue = df_analogue[~df_analogue.isna().any(axis=1)]
    
    # get variables
    df_seg = df_seg[var]
    df_analogue = df_analogue[var]

    # figure
    fig = plt.figure(figsize=(width, height))
    gs = gridspec.GridSpec(n_rows, n_cols, wspace=0.1, hspace=0)
    
    for i,vari in enumerate(var):
        var_seg = df_seg[vari].values
        var_analogue = df_analogue[vari].values
        var_anom = df_analogue[vari].values - df_seg[vari].values
        var_x = np.arange(0,df_seg.shape[0])
        
        # left column
        ax = plt.subplot(gs[i, 0])
        
        ax.plot(var_seg, '-', c='slateblue', linewidth=0.5)
        ax.plot(var_seg, 's', c='slateblue', markersize=2)
        ax.plot(var_analogue, '-', c='red', linewidth=0.5)
        ax.plot(var_analogue, '^', c='red', markersize=2)
        ax.plot([0,df_seg.shape[0]-1], [0,0], '--', c='silver')
        
        ax.set_xticks(ticks=np.arange(0,df_seg.shape[0]))
        ax.set_ylabel('{0} ({1})'.format(
                        labels[i], units[i]), fontweight='bold')
        
        # right column
        ax = plt.subplot(gs[i, 1])
        
        # variables: real / analogue / anomaly
        ax.plot([0,df_seg.shape[0]-1], [0,0], '--', c='silver')
        ax.plot(var_anom, '-', c='grey', linewidth=1)
        ax.scatter(var_x, var_anom, c=var_anom, cmap='seismic', 
                   vmin=-20, vmax=20, edgecolors='k', s=np.abs(var_anom)*10)
        
        ax.set_ylim([-20,20])
        ax.set_xticks(ticks=np.arange(0,df_seg.shape[0]))
        
    # titles
    ax = plt.subplot(gs[0, 0]); 
    ax.set_title('Parameter values', fontweight='bold')
    ax = plt.subplot(gs[0, 1]); 
    ax.set_title('Parameter anomalies (analogue - real)', fontweight='bold')
    
    return fig

###############################################################################
# SHYTCWAVES LIBRARY RESULTS
from scipy.io import loadmat
import xarray as xr
from .common import custom_cmap
from datetime import datetime

def outmat2xr(path_proj, case_i):
    'read output .mat file and returns xarray.Dataset'

    # case path
    p_mat = op.join(path_proj, '{0:04d}/output_main.mat'.format(case_i))
    
    # matlab dictionary
    dmat = loadmat(p_mat)

    # find output variable keys inside .mat file
    ks = list(set([x.split('_')[0] for x in dmat.keys()]))
    ks = [x for x in ks if x]  # remove empty values
    if 'Windv' in ks: ks.remove('Windv'); ks.append('Windv_x'); ks.append('Windv_y')

    # get dates from one key
    hsfs = sorted([x for x in dmat.keys() if ks[0] == x.split('_')[0]])
    dates_str = ['_'.join(x.split('_')[1:]) for x in hsfs]
    dates = [datetime.strptime(s,'%Y%m%d_%H%M%S') for s in dates_str]

    # iterate times and variables 
    l_times = []
    for ix_t, ds in enumerate(dates_str):
        xds_t = xr.Dataset()
        for vn in ks:
            xds_t[vn] = (('Y','X',), dmat['{0}_{1}'.format(vn, ds)])
        l_times.append(xds_t)

    # join at times dim
    xds_out = xr.concat(l_times, dim='time')
    xds_out = xds_out.assign_coords(time=dates)

    return xds_out

def axplot_var_map(ax, XX, YY, vv,
                   vmin = None, vmax = None,
                   cmap = plt.get_cmap('seismic'),
                  ):
    'plot 2D map with variable data'

    # cplot v lims
    if vmin == None: vmin = np.nanmin(vv)
    if vmax == None: vmax = np.nanmax(vv)

    # plot variable 2D map
    pm = ax.pcolormesh(
        XX, YY, vv,
        cmap=cmap,
        vmin=vmin, vmax=vmax,
        shading='auto',
        #zorder=0,
    )

    # fix axes
    ax.set_xlim(XX[0], XX[-1])
    ax.set_ylim(YY[0], YY[-1])

    # return pcolormesh
    return pm

def axplot_grafiti(ax, xds_case, case_number, var_name,
                   storm_track_list=[], 
                   vmin=None, vmax=None):
    '''
    # TODO: documentar
    '''

    # plot grafiti
    xds_var = xds_case[var_name]

    # get mesh data from output dataset
    xa, ya = 'X', 'Y'
    X = xds_var[xa].values[:]
    Y = xds_var[ya].values[:]

    # grafiti
    xds_var_max = xds_var.max(dim='time')
    var_max = xds_var_max.values[:]

    # colormap
    if var_name=='W': 	 
        ccmap = custom_cmap(100, 'plasma_r', 0.05, 0.9, 'viridis', 0.2, 1)
    if var_name=='Hsig': 
        ccmap = custom_cmap(100, 'YlOrRd', 0.09, 0.9, 'YlGnBu_r', 0, 0.88)
    
    pc = axplot_var_map(
        ax, X, Y, var_max,  # TODO output.T, aqui no
        vmin = vmin, vmax = vmax,
        cmap = ccmap,
    )

    # plot storm track
    if storm_track_list:
        st = storm_track_list[0]  # select storm track for this case
        ax.plot(st.lon, st.lat, '-', linewidth=4, 
                color='black', label='Storm Track')
        
    # red dashed line that bounds storage domain
    xlen, ylen = X.size, Y.size
    if xlen > ylen:
        ix0 = xlen-ylen
        ax.plot([X[ix0], X[ix0]], [Y[0], Y[-1]], '--', c='grey')

    # number
    ax.text(0.02, 0.02, case_number, color='fuchsia', fontweight='bold', 
            fontsize=20, transform=ax.transAxes)

    plt.axis('scaled')
    plt.xlim(X[0], X[-1])
    plt.ylim(Y[0], Y[-1])
    plt.axis('off')

    # equal axis
    ax.set_aspect('equal', 'box')

    return pc

def get_grid_sizes(xds_mda, resolution):
    '''
    resolution  - computational grid extent [m]
    '''
    
    res = resolution
    xr_ls, yr_ls = [], []
    
    for num in [0,1,2]:

        xsize_1 = xds_mda.x_size_left.values[num]*10**3
        xsize_2 = xds_mda.x_size_right.values[num]*10**3
        ysize_1 = xds_mda.y_size_left.values[num]*10**3
        ysize_2 = xds_mda.y_size_right.values[num]*10**3
        
        # create domain centered at (x,y)=(0,0)
        lon = np.arange(-xsize_1, xsize_2, res)
        lat = np.arange(-ysize_1, ysize_2, res)
        
        # get XY
        lon_len = lon[-1]-lon[0]
        lon = np.linspace(lon[0], 
                          lon[0] + lon_len + lon_len / int(lon_len/res), 
                          int(lon_len/res))
        lat_len = lat[-1]-lat[0]
        lat = np.linspace(lat[0], 
                          lat[0] + lat_len + lat_len / int(lat_len/res), 
                          int(lat_len/res))
        
        xr_ls.append(lon)
        yr_ls.append(lat)
    
    return xr_ls, yr_ls    # list small/medium/large

def add_grid_coords(xds, icase, xds_mda, xr_ls, yr_ls):
    '''
    xds         - (xr.Dataset) mat output
    resolution  - computational grid extent [m]
    icase       - number of mda cases
    '''
    
    # get grid size
    pos_small = np.where(icase == xds_mda.indices_small)[0]
    pos_medium = np.where(icase == xds_mda.indices_medium)[0]
    pos_large = np.where(icase == xds_mda.indices_large)[0]
    
    if len(pos_small)==1:   num=0
    if len(pos_medium)==1:  num = 1
    if len(pos_large)==1:   num = 2
    
    if len(pos_small)==0 and len(pos_medium)==0 and len(pos_large)==0:   num=0
    
    # set X and Y values
    X, Y = xr_ls[num], yr_ls[num]
    xds = xds.assign_coords(X=X)
    xds = xds.assign_coords(Y=Y)

    return xds

def plot_multi_grid(path_proj, resolution, storm_track_list=[], 
                      case_list=np.arange(0,10), var_name='Hsig', 
                      var_min=0, var_max=50,
                      n_rows=5, n_cols=5, width=34, height=30, show=True):
    '''
    # TODO: documentar
    '''

    # load mda file
    p_mda = r'/media/administrador/HD8/BlueMath/data/stopmotion/'
    xds_mda = xr.open_dataset(op.join(p_mda, 'shytcwaves_mda_indices.nc'))
    
    # get size coords (small, medium, large)
    xr_ls, yr_ls = get_grid_sizes(xds_mda, resolution)

    # figure
    fig = plt.figure(figsize=(width, height))

    gs = gridspec.GridSpec(n_rows, n_cols, wspace=0, hspace=0)
    gr, gc = 0, 0

    for ic, case_i in enumerate(case_list):
        
        sys.stdout.write('\rCase %d' %case_i)
        sys.stdout.flush()

        # load case .mat to .nc
        xds_out = outmat2xr(path_proj, case_i)
        
        # add case coords
        xds_out = add_grid_coords(xds_out, case_i, xds_mda, xr_ls, yr_ls)

        # plot variable times
        ax = plt.subplot(gs[gr, gc])
        
        pc = axplot_grafiti(
            ax, xds_out, case_i, var_name,
            storm_track_list=[storm_track_list[ic]],
            vmin=var_min, vmax=var_max
        )

        # get lower positions
        if gr==n_rows-1 and gc==0:
            pax_l = ax.get_position()
        elif gr==n_rows-1 and gc==n_cols-1:
            pax_r = ax.get_position()

        # counter
        gc += 1
        if gc >= n_cols:
            gc = 0
            gr += 1

    cbar_ax = fig.add_axes([pax_l.x0, pax_l.y0-0.05, pax_r.x1 - pax_l.x0, 0.02])
    cb = fig.colorbar(pc, cax=cbar_ax, orientation='horizontal')
#    cb.set_label(label='{0} ({1})'.format(var_name, var_units), size=20, weight='bold')
    cb.set_label(label='{0}'.format(var_name), size=20, weight='bold')
    cb.ax.tick_params(labelsize=15)

    # show and return figure
    if show:    plt.show()
    else:       plt.close()
    
    return fig

import sys
def plot_multi_grid_W(path_wind, resolution, storm_track_list=[], 
                      case_list=np.arange(0,10), var_name='W', 
                      var_min=0, var_max=50,
                      n_rows=5, n_cols=5, width=34, height=30, show=True):
    '''
    # TODO: documentar
    '''

    # figure
    fig = plt.figure(figsize=(width, height))

    gs = gridspec.GridSpec(n_rows, n_cols, wspace=0, hspace=0)
    gr, gc = 0, 0

    for ic, case_i in enumerate(case_list):
        
        sys.stdout.write('\rCase %d' %case_i)
        sys.stdout.flush()

        
        # load case .nc
        path = op.join(path_wind, '{0:04d}'.format(case_i))
        xds_wind = xr.open_dataset(op.join(path, 'vortex_wind_main.nc'))
        
        # plot variable times
        ax = plt.subplot(gs[gr, gc])
        
        pc = axplot_grafiti(
            ax, xds_wind, case_i, var_name,
            storm_track_list=[storm_track_list[ic]],
            vmin=var_min, vmax=var_max
        )

        # get lower positions
        if gr==n_rows-1 and gc==0:
            pax_l = ax.get_position()
        elif gr==n_rows-1 and gc==n_cols-1:
            pax_r = ax.get_position()

        # counter
        gc += 1
        if gc >= n_cols:
            gc = 0
            gr += 1

    cbar_ax = fig.add_axes([pax_l.x0, pax_l.y0-0.05, pax_r.x1 - pax_l.x0, 0.02])
    cb = fig.colorbar(pc, cax=cbar_ax, orientation='horizontal')
#    cb.set_label(label='{0} ({1})'.format(var_name, var_units), size=20, weight='bold')
    cb.set_label(label='{0}'.format(var_name), size=20, weight='bold')
    cb.ax.tick_params(labelsize=15)

    # show and return figure
    if show:    plt.show()
    else:       plt.close()
    
    return fig


###############################################################################
# ENSEMBLE (shytcwaves postprocessing)
  
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib.patheffects as pe

# colormaps
cmap_wind = custom_cmap(100, 'plasma_r', 0.05, 0.9, 'viridis', 0.2, 1)
cmap_hs = custom_cmap(15, 'YlOrRd', 0.15, 0.9, 'YlGnBu_r', 0, 0.85)
cmap_tp = 'nipy_spectral_r'
cmap_sp = 'twilight_shifted'
cmap_probs = ListedColormap(['forestgreen','limegreen','chartreuse','yellow',
                             'gold','goldenrod','darkorange','crimson','brown','darkmagenta'])

def get_cmap_bmu(xs_shy_bulk):
    
    color_list = ['darkblue','steelblue','lightsteelblue','hotpink','violet', \
                  'pink','limegreen','lawngreen','greenyellow','firebrick', \
                  'crimson','darksalmon','deepskyblue','lightblue', \
                  'saddlebrown','rosybrown','y','gold','khaki','green', \
                  'mediumaquamarine','palegreen','grey','silver','darkorange',\
                  'burlywood','peachpuff','blueviolet','mediumpurple','thistle',]
    
    ncolor = xs_shy_bulk.case.size
    if ncolor > len(color_list):
        color_list *= int(np.ceil(ncolor / len(color_list)))
        
    return ListedColormap(color_list[:ncolor])

###############################################################################
# ENSEMBLE  (shytcwaves)

def aux_cbar_var(fig, im, ttl='Hs [m]', rect=[0.02, 0, 0.19, 0.15], 
                 extend=None, xlabel=False, orientation='vertical'):
    'plot colorbar'
    
    gs = gridspec.GridSpec(1,1)
    ax0=fig.add_subplot(gs[0])
    if xlabel:
        plt.colorbar(im, cax=ax0, orientation=orientation, extend=extend)
        ax0.set_xlabel(ttl, fontsize=14)
    else:
        plt.colorbar(im, cax=ax0, orientation=orientation, extend=extend, 
                     label=ttl)
    gs.tight_layout(fig, rect=rect)
    
    return gs

def aux_cbar_bmu(fig, im, ttl='bmu', rect=[0.02, 0, 0.19, 0.15], 
                 num_ticks=10, orientation='vertical'):
    'plot colorbar'

    gs = gridspec.GridSpec(1,1)
    ax0 = fig.add_subplot(gs[0])
    
    cb = plt.colorbar(im, cax=ax0, orientation=orientation, label=ttl, 
                      ticks = np.arange(num_ticks)+.5)
    cb.set_ticklabels(range(0, num_ticks + 1))
#     ax0.set_xlabel(ttl, fontsize=14, fontweight='bold')
    gs.tight_layout(fig, rect=rect)

    return gs

def aux_series_swath(fig, xs_shy_bulk, lonpts, latpts, numpt=0,
                     var=['hsbmu','tpbmu'], cmap='hot', sameylim=None,
                     xs_dyn=None, var_dyn=['Hsig','Tp'], 
                     rect=[0, 0.2, 0.2, .9]):
    '''
    Plot of time series.
    xs_shy_bulk     - (xr.Dataset) shytcwaves bulk params
    lonpts,latpts   - (array) control point coordinates
    numpt           - (float) index of control points for plotting
    var             - (list) strings of variables for series

    xs_dyn          - (xr.Dataset) optional: dynamic results
    var_dyn         - (list) optional: strings of variables for series
    '''

    xds = xs_shy_bulk

    # var max
    pargmin = np.argmin(np.abs(xds.lon.values - lonpts[numpt]) + \
                        np.abs(xds.lat.values - latpts[numpt]))
    vmax = [np.nanmax(xds[var[0]][pargmin]), np.nanmax(xds[var[1]][pargmin])]
    vmin = [0,0]
    cmapmax = xds.case.size
    
    if sameylim:    # common ylim for all point series
        vmax = sameylim

    if xs_dyn:      # dynamic simulation
        xs_sel = xs_dyn.interp(lon=lonpts[numpt], lat=latpts[numpt])
        vmax_dyn = [np.nanmax(xs_sel[var_dyn[0]]), np.nanmax(xs_sel[var_dyn[1]])]
        vmax = [np.ceil(np.nanmax([vmax[0], vmax_dyn[0]])), 
                np.ceil(np.nanmax([vmax[1], vmax_dyn[1]]))]
        
    # figure
    gs = gridspec.GridSpec(len(var)+1, 1, hspace=0, wspace=0)
    
    for ivar in range(len(var)):   # shytcwaves hs,tp

        ax = fig.add_subplot(gs[ivar])
        
        for i in np.arange(0, xds.time.size, 6):  # every 6h
            plt.plot([xds.time.values[i], xds.time.values[i]], 
                     [vmin[ivar], vmax[ivar]], '--', c='silver', linewidth=1)
        ax.scatter(xds.time, xds[var[ivar]][pargmin], c=xds['bmu'][pargmin], 
                   cmap=cmap, vmin=0, vmax=cmapmax, alpha=.75); 
        ax.set_xlim([xds.time[0], xds.time[-1]])
        ax.set_ylim([vmin[ivar], vmax[ivar]]); 
        ax.set_ylabel(var[ivar], fontweight='bold')
        ax.set_xticklabels([]); ax.set_xticks([])
        
        if xs_dyn:      # dynamic simulation
#            xs_sel = xs_dyn.interp(lon=lonpts[numpt], lat=latpts[numpt])
            ax.plot(xs_dyn.time, xs_sel[var_dyn[ivar]], c='b', linewidth=1.1)

    ax2 = fig.add_subplot(gs[ivar+1]) # bmu
    
    im2 = ax2.scatter(xds.time, xds['bmu'][pargmin], c=xds['bmu'][pargmin], 
                      cmap=cmap, vmin=0, vmax=cmapmax)
    ax2.text(xds.time[0], xs_shy_bulk.case.size-5, 'point {0}'.format(numpt), 
            c='b', fontweight='bold', fontsize=12);
    ax2.set_xlim([xds.time[0], xds.time[-1]])
    ax2.set_ylim([0, xds.case.size])
    
    gs.tight_layout(fig, rect=rect, h_pad=0)

    return gs, im2

def aux_var_shy(fig, st, xs_shy, mesh_lo, mesh_la, lonpts, latpts, 
                var='hsbmu', itime=0, swath=False, cmap='hot', vmax=1, 
                rect=[0, 0.2, 0.2, .9], text=False, spine=False):
    '''
    Plot 2d map.
    st              - (pd.Dataframe) track coordinates
    xs_shy          - (xr.Dataset) shytcwaves bulk params
    meshlo,meshla   - (array) grid control points
    lonpts,latpts   - (array) control point coordinates
    var             - (string) variable for plotting
    itime           - (float) index of time array

    swath           - (bool) activates swath of maximum variable
    text            - (bool) activates control point labels
    spine           - (bool) activates spine color for shytcwaves plot
    '''
       
    if swath:   # maximum swath
        xds = xs_shy
        z = xds[var].max(dim='time').values.reshape((mesh_lo.shape))
#        z = np.nanmax(xds[var].values, axis=1).reshape((mesh_lo.shape))
    else:       # time data
        xds = xs_shy.isel(time=itime)
        z = xds[var].values.reshape((mesh_lo.shape))#[0], mesh_lo.shape[1])
            
    # figure
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])
    
    im = ax.pcolor(mesh_lo, mesh_la, z, cmap=cmap, vmin=0, vmax=vmax)
    
    ax.plot(st.lon, st.lat, 'k')
    ax.plot(lonpts, latpts, 'o', c='fuchsia', mec='k', mew=1)
    if text:
        for i in range(np.array(lonpts).size):
            ax.text(lonpts[i]-1, latpts[i], i, size=12, c='fuchsia',
                    path_effects=[pe.withStroke(linewidth=1.5, foreground="k")])
    if spine:
        plt.setp(ax.spines.values(), color='violet', linewidth=5)
        plt.setp([ax.get_xticklines(), ax.get_yticklines()], color='violet')
        plt.setp(ax.spines.values(), color='violet')

    ax.set_title('{0}'.format(var), fontweight='bold')
    ax.set_ylim([xds.lat.values[0], xds.lat.values[-1]])
    
    gs.tight_layout(fig, rect=rect)

    return gs, im

def plot_ensemble_series(st, xs_shy_bulk, #xs_shy_spec, xs_pts, 
                             mesh_lo, mesh_la, lonpts, latpts, sameylim=None,
                             text=True, width=20, height=10):
    '''
    Panel of swath shytcwaves results.
    st              - (pd.Dataframe) track coordinates
    xs_shy_bulk     - (xr.Dataset) shytcwaves bulk params
    meshlo,meshla   - (array) grid control points
    lonpts,latpts   - (array) control point coordinates
    '''
    
    # colormap
    cmap_bmu = get_cmap_bmu(xs_shy_bulk)
    
    # maximum variables
#    wind_max = np.nanmax(xs_wind.W.values)
    hs_max = np.nanmax(xs_shy_bulk.hswath.values)
    tp_max = np.nanmax(xs_shy_bulk.tswath.values)
    cmapmax = xs_shy_bulk.case.size

    # figure
    fig = plt.figure(figsize=[width,height], dpi=200)

    # column nº1: vortex, hsbmu, tpbmu, bmu
    left, right = 0, .15
#     gs, im1 = aux_wind(fig, xs_wind, st, swath=True, vmax=wind_max, 
#                        cmap=cmap_wind, rect=[left, .75, right, 1])
    gs, im2 = aux_var_shy(fig, st, xs_shy_bulk, mesh_lo, mesh_la, lonpts, latpts, 
                          var='hsbmu', swath=True, cmap=cmap_hs, vmax=hs_max, 
                          text=text, rect=[left, .66, right, 1])
    gs, im3 = aux_var_shy(fig, st, xs_shy_bulk, mesh_lo, mesh_la, lonpts, latpts, 
                          var='tpbmu', swath=True, cmap=cmap_tp, vmax=tp_max, 
                          text=text, rect=[left, .33, right, .66])
    gs, im4 = aux_var_shy(fig, st, xs_shy_bulk, mesh_lo, mesh_la, lonpts, latpts, 
                          var='bmu', swath=True, cmap=cmap_bmu, vmax=cmapmax, 
                          text=text, rect=[left, 0, right, .33])

    # column nº2: colorbars
#     aux_cbar_var(fig, im1, ttl='W [m/s]', rect=[right, .75, .2, 1-.02],
#                  xlabel=True)
    aux_cbar_var(fig, im2, ttl='Hs [m]', rect=[right, .66, .2, 1-.02],
                 xlabel=True)
    aux_cbar_var(fig, im3, ttl='Tp [s]', rect=[right, .33, .2, .66-.02],
                 xlabel=True)
    aux_cbar_var(fig, im4, ttl='bmu []', rect=[right, 0, .2, .33-.02],
                 xlabel=True)

    # row nº3: series
    if sameylim:
        pargmin = [np.argmin(np.abs(xs_shy_bulk.lon.values - lonpts[i]) + \
                             np.abs(xs_shy_bulk.lat.values - latpts[i])) \
                   for i in range(len(lonpts))]
        sameylim = [np.nanmax(xs_shy_bulk['hsbmu'][pargmin]), 
                    np.nanmax(xs_shy_bulk['tpbmu'][pargmin])]  # vmax
        
    gs, im11 = aux_series_swath(fig, xs_shy_bulk, lonpts, latpts, numpt=0,
                                sameylim=sameylim, cmap=cmap_bmu, rect=[.2, .5, .6, 1])
    gs, im11 = aux_series_swath(fig, xs_shy_bulk, lonpts, latpts, numpt=1,
                                sameylim=sameylim, cmap=cmap_bmu, rect=[.2, 0, .6, .5])
    gs, im11 = aux_series_swath(fig, xs_shy_bulk, lonpts, latpts, numpt=2,
                                sameylim=sameylim, cmap=cmap_bmu, rect=[.6, .5, 1, 1])
    gs, im11 = aux_series_swath(fig, xs_shy_bulk, lonpts, latpts, numpt=3,
                                sameylim=sameylim, cmap=cmap_bmu, rect=[.6, 0, 1, .5])
    
#    aux_cbar_bmu(fig, im11, ttl='', rect=[.48, 0, .53, .66], 
#                 num_ticks=xs_shy_bulk.case.size)

    return fig

def aux_wind(fig, st, xs_wind, var='W', cmap='hot', vmax=1, 
             itime=0, swath=False, rect=[0, 0.2, 0.2, .9]):
    '''
    Plot 2d map.
    st              - (pd.Dataframe) track coordinates
    xs_wind         - (xr.Dataset) wind vortex
    var             - (string) variable for plotting
    itime           - (float) index of time array

    swath           - (bool) activates swath of maximum variable
    '''
    
    if swath:   # maximum swath
        xds = xs_wind
        z = xds[var].max(dim='time')
        ttl = 'Wind vortex  ({0} - {1})'.format(xds.time.values[0].astype('datetime64[h]'), 
                                                xds.time.values[1].astype('datetime64[h]'))
    else:       # time data
        xds = xs_wind.isel(time=itime)
        z = xds[var]
        ttl = 'Wind vortex  {0}'.format(xds.time.values)
            
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    im = ax.pcolor(xds.lon, xds.lat, z, cmap=cmap, vmin=0, vmax=vmax)
    
    ax.plot(st.lon, st.lat, 'k')
    ax.set_title(ttl, fontweight='bold')
    gs.tight_layout(fig, rect=rect)

    return gs, im

def aux_spectrum_vars(xds_ls, itime=0, swath=False):
    '''
    Gets spectra variables for time index (itime) or max (swath). If the list
    has two variables, calculates the difference 
    
    xds_ls          - (list) xr.Dataset of directional wave spectra
    itime           - (float) index of time array
    swath           - (bool) activates swath of maximum variable
    '''
    
    if len(xds_ls)==1:      # single variable
        xds = xds_ls[0]
        x, y = np.deg2rad(xds.dir), xds.freq, 
        if swath:   z = xds.efth.max(dim='time')
        else:       z = xds.isel(time=itime).efth
        return x, y, [z]
            
    elif len(xds_ls)==2:    # difference of two variables
        xds1, xds2 = xds_ls[0], xds_ls[1]
        x, y = np.deg2rad(xds2.dir), xds2.freq
        if swath:   z2, z1 = xds2.efth.max(dim='time'), xds1.efth.max(dim='time')
        else:       z2, z1 = xds2.isel(time=itime).efth, xds1.isel(time=itime).efth
        return x, y, [z2,z1]

#    if len(xds_ls)==1:      # single variable
#        xds = xds_ls[0]
#        x, y = np.deg2rad(xds.dir), xds.freq, 
#        if swath:   z = xds.efth.max(dim='time')
#        else:       z = xds.isel(time=itime).efth
#            
#    elif len(xds_ls)==2:    # difference of two variables
#        xds1, xds2 = xds_ls[0], xds_ls[1]
#        x, y = np.deg2rad(xds2.dir), xds2.freq
#        if swath:   z = xds2.efth.max(dim='time') - xds1.efth.max(dim='time')
#        else:       z = xds2.isel(time=itime).efth - xds1.isel(time=itime).efth
#        
#    return x, y, z

def aux_spectrum(fig, xds_ls=[], itime=0, swath=False, cmap='hot', 
                 vmin=0, vmax=0.2, ymin=0.03, ylim=0.49, 
                 rect=[0, 0, 1, 1], spine=False):
    '''
    Plot 2d map.
    xds_ls          - (list) xr.Dataset of directional wave spectra
    itime           - (float) index of time array

    swath           - (bool) activates swath of maximum variable
    spine           - (bool) activates spine color for shytcwaves plot
    '''
    
#    # calculate spectrum variables
#    if len(xds_ls)==1:      # single variable
#        xds = xds_ls[0]
#        x, y = np.deg2rad(xds.dir), xds.freq, 
#        if swath:   z = xds.efth.max(dim='time')
#        else:       z = xds.isel(time=itime).efth
#            
#    elif len(xds_ls)==2:    # difference of two variables
#        xds1, xds2 = xds_ls[0], xds_ls[1]
#        x, y = np.deg2rad(xds2.dir), xds2.freq
#        if swath:   z2, z1 = xds2.efth.max(dim='time'), xds1.efth.max(dim='time')
#        else:       z2, z1 = xds2.isel(time=itime).efth, xds1.isel(time=itime).efth
        
#    x,y,z = aux_spectrum_vars(xds_ls, itime=itime, swath=swath)
    x,y,z = aux_spectrum_vars(xds_ls, itime=itime, swath=swath)
    
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0], projection='polar')

    xx = np.append(x, x[0])
    yy = np.append(0, y)
    
    if len(z)==1:      # single variable
        z = z[0]
        zz = np.column_stack((z[:,:], z[:,-1]))
        im = ax.pcolormesh(xx, yy, np.sqrt(zz), vmin=vmin, vmax=vmax)  
    
    elif len(z)==2:    # difference of two variables
        z2, z1 = z[0], z[1]
        zz2 = np.column_stack((z2[:,:], z2[:,-1]))
        zz1 = np.column_stack((z1[:,:], z1[:,-1]))
        im = ax.pcolormesh(xx, yy, np.sqrt(zz2) - np.sqrt(zz1), vmin=vmin, vmax=vmax)  
        
    im.set_cmap(cmap)
    ax.set_theta_zero_location('N', offset=0)
    ax.set_theta_direction(-1)
    ax.set_ylim(ymin, ylim)
    ax.tick_params(axis='y', colors='plum', labelsize=10, 
                   grid_linestyle=':', grid_alpha=0.75, grid_color='plum')
    ax.tick_params(axis='x', colors='purple', labelsize=10, 
                   grid_linestyle=':', grid_alpha=0.75, grid_color='plum')
    ax.grid(color='plum', linestyle='--', linewidth=0.7, alpha=1)
    
    if spine:
        plt.setp(ax.spines.values(), color='violet', linewidth=5)
        plt.setp([ax.get_xticklines(), ax.get_yticklines()], color='violet')
        plt.setp(ax.spines.values(), color='violet')
    gs.tight_layout(fig, rect=rect)

    return gs, im

def plot_ensemble_sp(st, xs_wind, xs_shy_bulk, xs_shy_sp, #xs_shy_spec, 
                     mesh_lo, mesh_la, time_i=0, spec_max=0.2, #lonpts, latpts, 
                     swath=False, text=False, width=20, height=10):
    '''
    Plot panel with 2d map variables and superpoint directional wave spectra
    st              - (pd.Dataframe) track coordinates
    xs_wind         - (xr.Dataset) vortex winds
    xs_shy_bulk     - (xr.Dataset) shytcwaves bulk mesh
#    xs_shy_spec     - (xr.Dataset) shytcwaves spectra cps
    xs_shy_sp       - (xr.Dataset) shytcwaves superpoint
    meshlo,meshla   - (array) grid control points
    lonpts,latpts   - (array) control point coordinates
    itime           - (float) index of time array
    spec_max        - (float) maximum energy to plot spectrum

    swath           - (bool) activates swath of maximum variable
    text            - (bool) activates control point labels
    '''
    # TODO add landmask
    
    # colormaps
    cmap_bmu = get_cmap_bmu(xs_shy_bulk)
    
    # varmax
    wind_max = np.nanmax(xs_wind['W'].values)
    hs_max = np.nanmax(xs_shy_bulk['hswath'].values)
    tp_max = np.nanmax(xs_shy_bulk['tswath'].values)
    bmu_max = xs_shy_bulk['case'].size
    
    lonpts = xs_shy_sp.lon.values
    latpts = xs_shy_sp.lat.values

    # figure
    fig = plt.figure(figsize=[width,height], dpi=200)
    
    # row nº1: vortex, bmu
    bottom, up = .55, 1
    gs, im1 = aux_wind(fig, st, xs_wind, var='W', cmap=cmap_wind, vmax=wind_max, 
                       itime=time_i, swath=swath, rect=[0, bottom, 0.25, up])
    aux_cbar_var(fig, im1, ttl='W [m/s]', orientation='horizontal',
                 rect=[0.01, bottom-0.1, .25, bottom])

    gs, im2 = aux_var_shy(fig, st, xs_shy_bulk, mesh_lo, mesh_la, #  xs_shy_spec, 
                          lonpts, latpts, var='bmu', itime=time_i, 
                          swath=swath, cmap=cmap_bmu, vmax=bmu_max, 
                          text=text, spine=False, rect=[.25, bottom, .5, up])
    aux_cbar_bmu(fig, im2, ttl='segment case', rect=[.26, bottom-0.1, .5, bottom], 
                 num_ticks=xs_shy_bulk.case.size, orientation='horizontal')

    # row nº2: Hsbmu, Tpbmu
    bottom, up = 0, .45
    gs, im3 = aux_var_shy(fig, st, xs_shy_bulk, mesh_lo, mesh_la, #   xs_shy_spec, 
                          lonpts, latpts, var='hsbmu', itime=time_i, 
                          swath=swath, cmap=cmap_hs, vmax=hs_max, 
                          text=text, spine=False, rect=[0, bottom, .25, up])
    aux_cbar_var(fig, im3, ttl='Hs [m]', rect=[.01, bottom-.1, .25, bottom], 
                 orientation='horizontal')

    gs, im4 = aux_var_shy(fig, st, xs_shy_bulk, mesh_lo, mesh_la, #   xs_shy_spec, 
                          lonpts, latpts, var='tpbmu', itime=time_i, 
                          swath=swath, cmap=cmap_tp, vmax=tp_max, 
                          text=text, spine=False, rect=[.25, bottom, .5, up])
    aux_cbar_var(fig, im4, ttl='Tp [s]', orientation='horizontal', 
                 rect=[.26, bottom-.1, .5, bottom])

    # right: spectra
    bottom, up = 0, .95
    gs, im8 = aux_spectrum(fig, xds_ls=[xs_shy_sp], itime=time_i, vmax=spec_max, swath=swath, 
                           cmap=cmap_sp, rect=[.5, bottom, 1, up], spine=False)
    aux_cbar_var(fig, im8, ttl='sqrt(Variance density)', orientation='horizontal', 
                 rect=[.5, bottom-.1, 1, bottom], extend='max', xlabel=True)

    return fig     



###############################################################################
# ENSEMBLE  (validation with dynamic simulations)

def aux_var(fig, st, xs_main, lonpts, latpts,  #xs_shy_sp, 
            var='Hsig', itime=0, swath=False, cmap='hot', vmax=1, 
            rect=[0, 0.2, 0.2, .9], text=False):
    '''
    Plot 2d map.
    st              - (pd.Dataframe) track coordinates
    xs_main         - (xr.Dataset) dynamic simulation
    xs_shy          - (xr.Dataset) shytcwaves bulk params
    lonpts,latpts   - (array) control point coordinates
    
    var             - (string) variable for plotting
    itime           - (float) index of time array

    swath           - (bool) activates swath of maximum variable
    text            - (bool) activates control point labels
    '''
       
    if swath:   # maximum swath
        xds = xs_main
        z = xds[var].max(dim='time')
    else:       # time data
        xds = xs_main.isel(time=itime)
        z = xds[var].values
            
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])
    
    im = ax.pcolor(xds.lon, xds.lat, z, cmap=cmap, vmin=0, vmax=vmax)
    
    ax.plot(st.lon, st.lat, 'k')
    ax.plot(lonpts, latpts, 'o', c='fuchsia', mec='k', mew=1)  # control points
    if text:
        for i in range(np.array(lonpts).size):
            ax.text(lonpts[i]-1, latpts[i], i, size=12, c='fuchsia',
                    path_effects=[pe.withStroke(linewidth=1.5, foreground="k")])
#    ax.plot(xs_shy_sp.lon.values, xs_shy_sp.lat.values,   # control points
#            'o', c='fuchsia', mec='k', mew=1)
#    if text:
#        for i in range(xs_shy_sp.lon.size):
#            ax.text(xs_shy_sp.lon[i]-1, xs_shy_sp.lat[i], i,
#                    size=12, c='fuchsia',
#                    path_effects=[pe.withStroke(linewidth=1.5, foreground="k")])

    ax.set_title('{0}'.format(var), fontweight='bold')
    ax.set_ylim([xds.lat.values[0], xds.lat.values[-1]])
    
    gs.tight_layout(fig, rect=rect)

    return gs, im


def plot_ensemble_val_series(st, xs_shy_bulk, xs_dyn, xs_wind, #xs_pts, #xs_shy_spec, 
                             mesh_lo, mesh_la, lonpts, latpts, sameylim=None,
                             text=True, width=20, height=10):
    '''
    Panel of shytcwaves series and 2d maps (dynamic and shy).
    st              - (pd.Dataframe) track coordinates
    xs_shy_bulk     - (xr.Dataset) shytcwaves bulk params
    xs_dyn          - (xr.Dataset) dynamic simulation
    xs_wind         - (xr.Dataset) wind vortex

    meshlo,meshla   - (array) grid control points
    lonpts,latpts   - (array) control point coordinates
    '''
    
    # colormaps
    cmap_bmu = get_cmap_bmu(xs_shy_bulk)
    
    # varmax
    wind_max = np.nanmax(xs_wind.W.values)
    hs_max = np.max([np.nanmax(xs_dyn.Hsig.values), 
                       np.nanmax(xs_shy_bulk.hswath.values)])
    tp_max = np.max([np.nanmax(xs_dyn.Hsig.values), 
                     np.nanmax(xs_shy_bulk.tswath.values)])
    cmapmax = xs_shy_bulk.case.size

    # figure
    fig = plt.figure(figsize=[width,height], dpi=200)

    # row nº1: vortex, Hsmain, Tpmain, hsbmu, tpbmu, bmu
    bottom, up = 0.72, 1
    gs, im1 = aux_wind(fig, st, xs_wind, cmap=cmap_wind, vmax=wind_max, 
                       swath=True, rect=[0, bottom, 0.16, up])
    gs, im2 = aux_var_shy(fig, st, xs_shy_bulk, mesh_lo, mesh_la, lonpts, latpts, #xs_pts, 
                          var='bmu', swath=True, cmap=cmap_bmu, vmax=cmapmax, 
                          text=text, rect=[.17, bottom, .33, up], spine=True)
    gs, im3 = aux_var(fig, st, xs_dyn, lonpts, latpts, #xs_pts, 
                      var='Hsig', swath=True, cmap=cmap_hs, vmax=hs_max, 
                      text=text, rect=[.34, bottom, 0.5, up])
    gs, im4 = aux_var_shy(fig, st, xs_shy_bulk, mesh_lo, mesh_la, lonpts, latpts, #xs_pts
                          var='hsbmu', swath=True, cmap=cmap_hs, vmax=hs_max, 
                          text=text, rect=[.5, bottom, .66, up], spine=True)
    gs, im5 = aux_var(fig, st, xs_dyn, lonpts, latpts, #xs_pts, 
                      var='Tp', swath=True, cmap=cmap_tp, vmax=tp_max, 
                      text=text, rect=[.67, bottom, .83, up])
    gs, im6 = aux_var_shy(fig, st, xs_shy_bulk, mesh_lo, mesh_la, lonpts, latpts, #xs_pts
                          var='tpbmu', swath=True, cmap=cmap_tp, vmax=tp_max, 
                          text=text, rect=[.84, bottom, 1, up], spine=True)

    # row nº2: colorbars
    bottom, up = 0.65, 0.75
    aux_cbar_var(fig, im1, ttl='W [m/s]', rect=[0.02, bottom, .16, up],
                 xlabel=True, orientation='horizontal')
    aux_cbar_var(fig, im2, ttl='bmu []', rect=[.19, bottom, .33, up],
                 xlabel=True, orientation='horizontal')
    aux_cbar_var(fig, im3, ttl='Hs [m]', rect=[.36, bottom, .66, up],
                 xlabel=True, orientation='horizontal')
    aux_cbar_var(fig, im5, ttl='Tp [s]', rect=[.69, bottom, 1, up],
                 xlabel=True, orientation='horizontal')

    # row nº3: series
    if sameylim:
        pargmin = [np.argmin(np.abs(xs_shy_bulk.lon.values - lonpts[i]) + \
                             np.abs(xs_shy_bulk.lat.values - latpts[i])) \
                   for i in range(len(lonpts))]
        sameylim = [np.nanmax(xs_shy_bulk['hsbmu'][pargmin]), 
                    np.nanmax(xs_shy_bulk['tpbmu'][pargmin])]   # vmax
    
    gs, im11 = aux_series_swath(fig, xs_shy_bulk, lonpts, latpts, numpt=0, # xs_pts, pargmin, 
                                xs_dyn=xs_dyn, cmap=cmap_bmu, sameylim=sameylim, 
                                rect=[0, .33, .48, .66])
    gs, im11 = aux_series_swath(fig, xs_shy_bulk, lonpts, latpts, numpt=1, # xs_pts, pargmin, 
                                xs_dyn=xs_dyn, cmap=cmap_bmu, sameylim=sameylim, 
                                rect=[0, 0, .48, .33])
    gs, im11 = aux_series_swath(fig, xs_shy_bulk, lonpts, latpts, numpt=2, # xs_pts, pargmin, 
                                xs_dyn=xs_dyn, cmap=cmap_bmu, sameylim=sameylim, 
                                rect=[.53, .33, 1, .66])
    gs, im11 = aux_series_swath(fig, xs_shy_bulk, lonpts, latpts, numpt=3, # xs_pts, pargmin, 
                                xs_dyn=xs_dyn, cmap=cmap_bmu, sameylim=sameylim, 
                                rect=[.53, 0, 1, .33])
    
    aux_cbar_bmu(fig, im11, ttl='', rect=[.48, 0, .53, .66], 
                 num_ticks=xs_shy_bulk.case.size)

    return fig

def plot_ensemble_val_sp(st, xs_shy_bulk, xs_shy_sp, #xs_shy_spec, 
                         xs_wind, xs_dyn, xs_dyn_sp, 
                         mesh_lo, mesh_la, time_i=0, spec_max=0.2, #lonpts, latpts, 
                         swath=False, text=False, width=20, height=10):
    '''
    Panel of shytcwaves superpoint spectra and 2d maps (dynamic and shy).
    st              - (pd.Dataframe) track coordinates
    xs_shy_bulk     - (xr.Dataset) shytcwaves bulk params
#    xs_shy_spec     - (xr.Dataset) shytcwaves spectra cps
    xs_shy_sp       - (xr.Dataset) shytcwaves superpoint
    xs_wind         - (xr.Dataset) wind vortex
    xs_dyn          - (xr.Dataset) dynamic simulation
    xs_dyn_sp       - (xr.Dataset) dynamic superpoint

    meshlo,meshla   - (array) grid control points
    lonpts,latpts   - (array) control point coordinates

    itime           - (float) index of time array
    spec_max        - (float) maximum energy to plot spectrum
    swath           - (bool) activates swath of maximum variable
    text            - (bool) activates control point labels
    '''
    # TODO add landmask
    
    # colormaps
    cmap_bmu = get_cmap_bmu(xs_shy_bulk)
    
    # varmax
    wind_max = np.nanmax(xs_wind['W'].values)
    hs_max = np.max([np.nanmax(xs_dyn['Hsig'].values), 
                       np.nanmax(xs_shy_bulk['hswath'].values)])
    tp_max = np.max([np.nanmax(xs_dyn['Tp'].values), 
                     np.nanmax(xs_shy_bulk['tswath'].values)])
    bmu_max = xs_shy_bulk['case'].size
    
    lonpts = xs_shy_sp.lon.values
    latpts = xs_shy_sp.lat.values

    # figure
    fig = plt.figure(figsize=[width,height], dpi=200)
    
    # row nº1: vortex, Hsmain, Tpmain
    bottom, up = 0.72, 1
    gs, im1 = aux_wind(fig, st, xs_wind, cmap=cmap_wind, vmax=wind_max, 
                       itime=time_i, swath=swath, rect=[0, bottom, 0.28, up])
    aux_cbar_var(fig, im1, ttl='W [m/s]', rect=[.28, bottom, .33, up-0.01])

    gs, im2 = aux_var(fig, st, xs_dyn, lonpts, latpts, var='Hsig', itime=time_i, #xs_shy_spec, 
                      swath=swath, cmap=cmap_hs, vmax=hs_max, text=text, 
                      rect=[0.33, bottom, 0.61, up])

    gs, im3 = aux_var(fig, st, xs_dyn, lonpts, latpts, var='Tp', itime=time_i, #xs_shy_spec, st, 
                      swath=swath, cmap=cmap_tp, vmax=tp_max, text=text, 
                      rect=[0.66, bottom, 0.95, up])

    
    # row nº2: bmu, Hsshy, Tpshy
    bottom, up = 0.44, 0.72
    gs, im4 = aux_var_shy(fig, st, xs_shy_bulk, mesh_lo, mesh_la, #xs_shy_spec, 
                          lonpts, latpts, var='bmu', itime=time_i, swath=swath, 
                          cmap=cmap_bmu, vmax=bmu_max, text=text, spine=True,
                          rect=[0, bottom, .28, up])
    aux_cbar_bmu(fig, im4, ttl='segment case', rect=[.28, bottom, .33, up-.01], 
                 num_ticks=xs_shy_bulk.case.size, orientation='vertical')

    gs, im5 = aux_var_shy(fig, st, xs_shy_bulk, mesh_lo, mesh_la, #xs_shy_spec, 
                          lonpts, latpts, var='hsbmu', itime=time_i, swath=swath, 
                          cmap=cmap_hs, vmax=hs_max, text=text, spine=True,
                          rect=[.33, bottom, .61, up])
    aux_cbar_var(fig, im5, ttl='Hs [m]', rect=[.61, bottom, .66, .99])

    gs, im6 = aux_var_shy(fig, st, xs_shy_bulk, mesh_lo, mesh_la, #xs_shy_spec, 
                          lonpts, latpts, var='tpbmu', itime=time_i, swath=swath,
                          cmap=cmap_tp, vmax=tp_max, text=text, spine=True,
                          rect=[.66, bottom, .95, up])
    aux_cbar_var(fig, im6, ttl='Tp [s]', rect=[.95, bottom, 1, .99])

    # row nº3: spmain, spshy
    bottom, up = 0, 0.44
    gs, im7 = aux_spectrum(fig, xds_ls=[xs_dyn_sp], itime=time_i, swath=swath, 
                           vmax=spec_max, cmap=cmap_sp, rect=[0, bottom, .3, up])

    gs, im8 = aux_spectrum(fig, xds_ls=[xs_shy_sp], itime=time_i, swath=swath, 
                           vmax=spec_max, cmap=cmap_sp, rect=[.3, bottom, .6, up], spine=True)
    aux_cbar_var(fig, im8, ttl='sqrt(Variance density)', 
                 rect=[.6, bottom, .65, up], extend='max', xlabel=True)

    gs, im9 = aux_spectrum(fig, xds_ls=[xs_dyn_sp, xs_shy_sp], itime=time_i, 
                           swath=swath, cmap='seismic', vmin=-spec_max, vmax=spec_max,
                           rect=[.65, bottom, .95, up])
    aux_cbar_var(fig, im9, ttl='Diff sqrt(Variance density)', 
                 rect=[.95, bottom, 1, up], extend='both', xlabel=True)

    return fig               
                       

def plot_ensemble_val_panel(st, xs_shy_bulk, xs_dyn, xs_wind, #xs_pts, #xs_shy_spec, 
                            xs_dyn_sp, xs_shy_sp, spec_max,
                             mesh_lo, mesh_la, lonpts, latpts, 
                             text=False, width=25, height=5):
    '''
    Panel of shytcwaves series and 2d maps (dynamic and shy).
    st              - (pd.Dataframe) track coordinates
    xs_shy_bulk     - (xr.Dataset) shytcwaves bulk params
    xs_dyn          - (xr.Dataset) dynamic simulation
    xs_wind         - (xr.Dataset) wind vortex
    xs_shy_sp       - (xr.Dataset) shytcwaves superpoint
    xs_dyn_sp       - (xr.Dataset) dynamic superpoint
    spec_max        - (float) maximumvariance energy

    meshlo,meshla   - (array) grid control points
    lonpts,latpts   - (array) control point coordinates
    '''
    
    # varmax
    wind_max = np.nanmax(xs_wind.W.values)
    hs_max = np.max([np.nanmax(xs_dyn.Hsig.values), 
                       np.nanmax(xs_shy_bulk.hswath.values)])

    # figure
    fig = plt.figure(figsize=[width,height], dpi=200)

    # row nº1: vortex, Hsmain, Tpmain, hsbmu, tpbmu, bmu
    bottom, up = 0.2, 0.95
    gs, im1 = aux_wind(fig, st, xs_wind, cmap=cmap_wind, vmax=wind_max, 
                       swath=True, rect=[0, bottom, 0.18, up])
    
    gs, im2 = aux_var(fig, st, xs_dyn, lonpts, latpts, #xs_pts, 
                      var='Hsig', swath=True, cmap=cmap_hs, vmax=hs_max, 
                      text=text, rect=[.18, bottom, 0.36, up])
    gs, im3 = aux_var_shy(fig, st, xs_shy_bulk, mesh_lo, mesh_la, lonpts, latpts, #xs_pts
                          var='hsbmu', swath=True, cmap=cmap_hs, vmax=hs_max, 
                          text=text, rect=[.36, bottom, .54, up], spine=True)

    bottom, up = 0.15, 1
    gs, im4 = aux_spectrum(fig, xds_ls=[xs_dyn_sp], swath=True, 
                           vmax=spec_max, cmap=cmap_sp, rect=[.55, bottom, .7, up])
    gs, im5 = aux_spectrum(fig, xds_ls=[xs_shy_sp], swath=True, spine=True, 
                           vmax=spec_max, cmap=cmap_sp, rect=[.7, bottom, .85, up])
    gs, im6 = aux_spectrum(fig, xds_ls=[xs_dyn_sp, xs_shy_sp], swath=True, 
                           cmap='seismic', vmin=-spec_max, vmax=spec_max,
                           rect=[.85, bottom, 1, up])

    # row nº2: colorbars
    bottom, up = 0, 0.2
    aux_cbar_var(fig, im1, ttl='W [m/s]', rect=[0.02, bottom, .18, up],
                 xlabel=True, orientation='horizontal')
    aux_cbar_var(fig, im3, ttl='Hs [m]', rect=[.2, bottom, .54, up],
                 xlabel=True, orientation='horizontal')
    aux_cbar_var(fig, im5, ttl='sqrt(Variance density)', orientation='horizontal', 
                 rect=[.57, bottom, .85, up], extend='max', xlabel=True)
    aux_cbar_var(fig, im6, ttl='Diff sqrt(Variance density)', orientation='horizontal', 
                 rect=[.86, bottom, .99, up], extend='both', xlabel=True)

    return fig

###############################################################################
# FORECAST
    
import cartopy
from cartopy import crs as ccrs
from .nonstationary import get_storm_color
from ..storms import get_category

def axplot_storm_track(ax, st, name, cat_colors=True, mksize=10, 
                       title_on=False, legend=False, linewidth=1):
    'plots storm track, category and target point' 
    
    # storm track parameters 
    xt = st.lon.values
    yt = st.lat.values
#    tc_name=str(name)
    
    # plot track
    ax.plot(xt, yt, '-', linewidth=linewidth,
            color='black', label='Storm Track', transform=ccrs.PlateCarree()) 
    
    # plot categories
    if cat_colors: # get category
        categ = np.array(get_category(st.p0)) 
        for c in range(7):
            if len(np.where(categ==c)[0])>0:
                lonc = xt[np.where(categ==c)[0]]
                latc = yt[np.where(categ==c)[0]]
                label = 'cat {0}'.format(c)
                if c==6: label = 'unknown'
                ax.plot(
                    lonc, latc, '.', color=get_storm_color(c),
                    transform=ccrs.PlateCarree(), markersize=mksize, label=label)

    if legend:  ax.legend()
    if title_on:
        plt.title('{0} TC: Pmin: {1:.2f} hPa / Vmean: {2:.2f} km/h'.format(
                    name, np.min(st.p0), np.mean(st.vf)*1.852), 
                    fontsize=16, fontweight='bold')

#def plot_predicted_tracks(st, name, figsize=[20,9]):
#    
#    fig = plt.figure(figsize=figsize) 
#    ax = plt.axes(projection = ccrs.PlateCarree(central_longitude=190))
#    axplot_storm_track(ax, st, name)
#    ax.stock_img() # cartopy land feature
#    land_50m = cartopy.feature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='darkgrey', facecolor='gainsboro', zorder=1)
#    ax.add_feature(land_50m)
#    ax.gridlines()
#    return fig, ax

def plot_forecast_grid_shy(p_ensemble, tcname, lon, lat, 
                           lo01=165, lo02=210, la01=-40, la02=-14, 
                           var='hswath', unit='m', 
                           nstorms=81, nrows=9, ncols=9, width=25, height=15):
    'plots predicted tracks, shytcwaves swath'
    
    # plot prediction results
    
    # lon0, lon1, lat0, lat1 = 170, 200, -30, -16
    mesh_lo, mesh_la = np.meshgrid(lon, lat)
#    lo01, lo02, la01, la02 = 165, 210, -40, -14
    
    # maxvalues
    varmax = []
    for ipred in range(nstorms):
        p_bulk = op.join(p_ensemble, 
                         '{0}_forecast_{1}_xds_shy_bulk.nc'.format(
                                 tcname, ipred))
        if op.isfile(p_bulk):  
            xds_shy_bulk = xr.open_dataset(p_bulk)
            varmax.append(np.nanmax(xds_shy_bulk[var].values))
        else: varmax.append(0)
    vmax_pos = np.argmax(varmax)
    varmax = np.nanmax(varmax)
    vmin=0
    fig, ax = plt.subplots(nrows, ncols, figsize=(width,height), 
                           sharex=True, sharey=True, 
                           subplot_kw=dict(projection=ccrs.PlateCarree()))
    
    for ipred in range(nstorms):
        ax = plt.subplot(nrows, ncols, ipred+1, 
                         projection=ccrs.PlateCarree(central_longitude=190))
    
        path_pkl = op.join(p_ensemble, 
                            '{0}_forecast_{1}_track.pkl'.format(tcname, ipred))
                
        p_bulk = op.join(p_ensemble, '{0}_forecast_{1}_xds_shy_bulk.nc'.format(
                tcname, ipred))

        if op.isfile(path_pkl):
            
            # cluster python version is upgraded --> protocol 5
            import pickle5 as pickle
            with open(path_pkl, "rb") as fh:
                df_storm = pickle.load(fh)
          
        # df_storm = pd.read_pickle(op.join(p_ensemble, 
        #                                   '{0}_forecast_{1}_track.pkl'.format(
        #                                           tcname, ipred)))
    
            axplot_storm_track(ax, df_storm, tcname, mksize=5)
    
        # cartopy land feature
        land_50m = cartopy.feature.NaturalEarthFeature('physical', 'land', '50m', 
                                                       alpha=0.4, edgecolor='k', 
                                                       facecolor='gainsboro', 
                                                       zorder=1)
        ax.add_feature(land_50m)
        # ax.gridlines()
        ax.plot([lon.min(),lon.max(),lon.max(),lon.min(),lon.min()], 
                 [lat.min(),lat.min(),lat.max(),lat.max(),lat.min()], 
                '--', c='silver', transform=ccrs.PlateCarree(), 
                linewidth=2, zorder=0);
    
        if op.isfile(p_bulk):  
            xds_shy_bulk = xr.open_dataset(p_bulk)
    
            mask=np.where(~np.isnan(xds_shy_bulk[var].values))
    #         pc = ax.scatter(xds_shy_bulk.lon.values[mask], xds_shy_bulk.lat.values[mask], c=xds_shy_bulk['hswath'].values[mask], #30, 
    #                            cmap=cmap_hs, vmin=vmin, vmax=varmax, transform=ccrs.PlateCarree())
            pc = ax.tricontourf(xds_shy_bulk.lon.values[mask], 
                                xds_shy_bulk.lat.values[mask], 
                                xds_shy_bulk[var].values[mask], #30, 
    #                             levels=np.linspace(vmin, 6, 50), cmap=cmap_hs, transform=ccrs.PlateCarree())
                                30, cmap=cmap_hs, vmin=vmin, vmax=varmax, 
                                transform=ccrs.PlateCarree())
            if ipred==vmax_pos: im=pc
    
        ax.set_extent([lo01, lo02, la01, la02])
        plt.xticks([]); plt.yticks([])
        
    fig.subplots_adjust(right=0.85, wspace=0.005, hspace=0.005)
    cbar_ax = fig.add_axes([0.88, 0.15, 0.02, 0.7])
    fig.colorbar(im, cax=cbar_ax, label='{0} [{1}]'.format(var, unit));
    
    return fig

def plot_forecast_probs(p_ensemble, tcname, mesh_lo, dict_tc, 
#                        hspred, pred_varmax, var_ls=[], ttl_ls=[], dt_ls=[],
                        lon0=184, lon1=190, lat0=-23, lat1=-16, # shy domain
                        lo01=165, lo02=210, la01=-40, la02=-14, # extent
                        nstorms=81, nrows=2, ncols=5, width=15, height=7, 
                        label=['hswath','probability'], unit=['m','%'],
                        contour_smooth=True):
    
#    var_ls = [None, pred_mean, pred_pct95, pred_pct99, pred_max, pred_occur, pred_2m, pred_5m, pred_8m, pred_12m]
#    ttl_ls = ['JTWC predicted tracks', 'Mean swath [m]', 'Pctl95 swath [m]', 'Pctl99 swath [m]', 'Max swath [m]', 
#              'Prediction probability [%]', 'Prob >2m [%]', 'Prob >5m [%]', 'Prob >8m [%]', 'Prob >12m [%]']
#    nrows, ncols = 2, 5
    xds_shy_bulk = dict_tc['xds_shy_bulk']
    pred_varmax = dict_tc['pred_varmax']
    var_ls = dict_tc['var_ls']
    ttl_ls = dict_tc['ttl_ls']
    dt_ls = dict_tc['dt_ls']

    fig, ax = plt.subplots(nrows, ncols, figsize=(width,height), 
                           sharex=True, sharey=True, 
                           subplot_kw=dict(projection=ccrs.PlateCarree()))
    
    for ipred, (ivar, idt) in enumerate(zip(var_ls, dt_ls)):
        
        ax = plt.subplot(nrows, ncols, ipred+1, 
                         projection = ccrs.PlateCarree(central_longitude=190))
    
        # cartopy land feature
        land_50m = cartopy.feature.NaturalEarthFeature('physical', 'land', '50m', 
                                                       alpha=0.4, edgecolor='k', 
                                                       facecolor='gainsboro', 
                                                       zorder=5)
        ax.add_feature(land_50m)
        # ax.gridlines()
        ax.plot([lon0,lon1,lon1,lon0,lon0], [lat0,lat0,lat1,lat1,lat0], 
                '--', c='silver', transform=ccrs.PlateCarree(), 
                linewidth=2, zorder=0);
        
        if ipred==0:
            for sti in range(nstorms):

                path_pkl = op.join(p_ensemble, 
                                    '{0}_forecast_{1}_track.pkl'.format(tcname, sti))
                        
                # p_bulk = op.join(p_ensemble, '{0}_forecast_{1}_xds_shy_bulk.nc'.format(
                #         tcname, ipred))
        
                if op.isfile(path_pkl):
                    
                    # cluster python version is upgraded --> protocol 5
                    import pickle5 as pickle
                    with open(path_pkl, "rb") as fh:
                        df_storm = pickle.load(fh)
                  
                # df_storm = pd.read_pickle(op.join(p_ensemble, 
                #                          '{0}_forecast_{1}_track.pkl'.format(
                #                                  tcname, sti)))
                        axplot_storm_track(ax, df_storm, tcname, mksize=7)
    
        if ipred>0 and ipred <=4:
            mask = np.where(~np.isnan(ivar))
            pc1 = ax.tricontourf(xds_shy_bulk.lon.values[mask], 
                                 xds_shy_bulk.lat.values[mask], ivar[mask], 
                                 30, cmap=cmap_hs, vmin=0, vmax=pred_varmax, 
                                 transform=ccrs.PlateCarree(), zorder=2); 
    
        if ipred >4:
            cmap = plt.get_cmap(cmap_probs)
            cmap.set_bad('k')
#            var_ = ivar / hspred.shape[0] *100
            data = np.ma.masked_equal(ivar, 0)
            pc2 = ax.pcolor(xds_shy_bulk.lon.values.reshape(mesh_lo.shape), 
                            xds_shy_bulk.lat.values.reshape(mesh_lo.shape), 
                            data.reshape(mesh_lo.shape), cmap=cmap, vmax=100, 
                            transform=ccrs.PlateCarree(), zorder=2);
        
        if ipred==1 or ipred >5:
#        if (ipred>0 and ipred <=4) or ipred >4:
            if ipred==1:    ccontour = 'w'
            elif ipred >5:  ccontour = 'k'
            data = np.ma.masked_equal(idt, 0).reshape(mesh_lo.shape)
            if contour_smooth: # resample data using cubic spline interpolation
                data = gaussian_filter(data, sigma=.6)
            CS = ax.contour(xds_shy_bulk.lon.values.reshape(mesh_lo.shape), 
                            xds_shy_bulk.lat.values.reshape(mesh_lo.shape), 
                            data, levels=np.arange(1,8), colors=ccontour, 
                            transform=ccrs.PlateCarree(), zorder=3);
            fmt = {}
            strs = ['1', '2', '3', '4', '5', '6', '7']
            for l, s in zip(CS.levels, strs):   fmt[l] = s
            ax.clabel(CS, CS.levels, inline=True, fmt=fmt, fontsize=12)

        ax.set_title(ttl_ls[ipred], fontweight='bold')
        ax.set_extent([lo01, lo02, la01, la02])
        plt.xticks([]); plt.yticks([])
    
    
    fig.subplots_adjust(right=0.87, wspace=0.005, hspace=0.15)
    cbar_ax = fig.add_axes([0.88, 0.55, 0.01, 0.3])
    fig.colorbar(pc1, cax=cbar_ax, label='{0} [{1}]'.format(label[0], unit[0]));
    cbar_ax = fig.add_axes([0.88, 0.15, 0.01, 0.3])
    fig.colorbar(pc2, cax=cbar_ax, label='{0} [{1}]'.format(label[1], unit[1]));
    
    return fig


















#def aux_series_swath0(fig, xs_pts, pargmin, xs_shy_bulk, numpt=0,
#                     var_main=['Hsig','Tp'], var=['hsbmu','tpbmu'], cmap='hot', 
#                     rect=[0, 0.2, 0.2, .9]):
#    
#    xs_pts = xs_pts.isel(point=numpt)
#    pargmin = pargmin[numpt]
#    vmin = [0,0]
#    vmax = [np.nanmax(xs_shy_bulk[var[0]][pargmin]), 
#            np.nanmax(xs_shy_bulk[var[1]][pargmin])]
#    cmapmax = xs_shy_bulk.case.size
#    
#    gs = gridspec.GridSpec(3,1, hspace=0, wspace=0)
#    
#    # hs
#    ax = fig.add_subplot(gs[0])
#    for i in np.arange(0, xs_pts.time.size, 6): 
#        plt.plot([xs_pts.time.values[i], xs_pts.time.values[i]], [vmin, vmax], 
#                 '--', c='silver', linewidth=1)
#    ax.plot(xs_pts.time, xs_pts[var_main[0]], c='b', linewidth=1.1)
#    ax.scatter(xs_shy_bulk.time, xs_shy_bulk[var[0]][pargmin], 
#               c=xs_shy_bulk['bmu'][pargmin], cmap=cmap, 
#               vmin=0, vmax=cmapmax, alpha=.75); 
#    ax.set_xlim([xs_pts.time[0], xs_pts.time[-1]])
#    ax.set_ylim([vmin[0], vmax[0]]); ax.set_ylabel(var[0], fontweight='bold')
#    ax.set_xticklabels([]); ax.set_xticks([])
#    
#    # tp
#    ax1 = fig.add_subplot(gs[1])
#    for i in np.arange(0, xs_pts.time.size, 6): 
#        plt.plot([xs_pts.time.values[i], xs_pts.time.values[i]], [vmin, vmax], 
#                 '--', c='silver', linewidth=1)
#    ax1.plot(xs_pts.time, xs_pts[var_main[1]], c='b', linewidth=1.1)
#    ax1.scatter(xs_shy_bulk.time, xs_shy_bulk[var[1]][pargmin], 
#                c=xs_shy_bulk['bmu'][pargmin], cmap=cmap, 
#                vmin=0, vmax=cmapmax, alpha=.75); 
#    ax1.set_xlim([xs_pts.time[0], xs_pts.time[-1]])
#    ax1.set_ylim([vmin[1], vmax[1]]); ax1.set_ylabel(var[1], fontweight='bold')
#    ax1.set_xticklabels([]); ax.set_xticks([])
#
#    # bmu
#    ax2 = fig.add_subplot(gs[2])
#    im2 = ax2.scatter(xs_shy_bulk.time, xs_shy_bulk['bmu'][pargmin], 
#                      c=xs_shy_bulk['bmu'][pargmin], cmap=cmap, 
#                      vmin=0, vmax=cmapmax)
#    ax2.text(xs_pts.time[0], 5, 'point {0}'.format(numpt), 
#            c='b', fontweight='bold', fontsize=12);
#    ax2.set_ylim([0, xs_shy_bulk.case.size])
#    ax2.set_xlim([xs_pts.time[0], xs_pts.time[-1]])
#    
#    gs.tight_layout(fig, rect=rect, h_pad=0)
#
#    return gs, im2
#
#def aux_var_shy0(fig, mesh_lo, mesh_la, xs_shy, xs_shy_sp, st, var='hsbmu', 
#                itime=0, swath=False, cmap='hot', vmax=1, 
#                rect=[0, 0.2, 0.2, .9], text=False, spine=False):
#       
#    if swath:   # maximum swath
#        xds = xs_shy
#        z = xds[var].max(dim='time').values.reshape((mesh_lo.shape))
##        z = np.nanmax(xds[var].values, axis=1).reshape((mesh_lo.shape))
#    else:       # time data
#        xds = xs_shy.isel(time=itime)
#        z = xds[var].values.reshape((mesh_lo.shape))#[0], mesh_lo.shape[1])
#            
#    gs = gridspec.GridSpec(1,1)
#    ax = fig.add_subplot(gs[0])
#    
#    im = ax.pcolor(mesh_lo, mesh_la, z, cmap=cmap, vmin=0, vmax=vmax)
#    
#    ax.plot(st.lon, st.lat, 'k')
#    ax.plot(xs_shy_sp.lon.values, xs_shy_sp.lat.values, 
#            'o', c='fuchsia', mec='k', mew=1)
#    if text:
#        for i in range(xs_shy_sp.lon.size):
#            ax.text(xs_shy_sp.lon[i]-1, xs_shy_sp.lat[i], i,
#                    size=12, c='fuchsia',
#                    path_effects=[pe.withStroke(linewidth=1.5, foreground="k")])
#
#    ax.set_title('{0}'.format(var), fontweight='bold')
#    ax.set_ylim([xds.lat.values[0], xds.lat.values[-1]])
#    
#    if spine:
#        plt.setp(ax.spines.values(), color='violet', linewidth=5)
#        plt.setp([ax.get_xticklines(), ax.get_yticklines()], color='violet')
#        plt.setp(ax.spines.values(), color='violet')
#    gs.tight_layout(fig, rect=rect)
#
#    return gs, im
#
#def plot_ensemble_series0(xs_pts, st, xs_shy_bulk, #xs_shy_spec, 
#                             mesh_lo, mesh_la, lonpts, latpts, 
#                             text=True, width=20, height=10):
#
##     xs_pts = xs_pts.isel(point=pts_ind)
#    
#    # colormaps
##    cmap_wind = custom_cmap(100, 'plasma_r', 0.05, 0.9, 'viridis', 0.2, 1)
#    cmap_bmu = get_cmap_bmu(xs_shy_bulk)
#    
#    # varmax
##    wind_max = np.nanmax(xs_wind.W.values)
#    hs_max = np.nanmax(xs_shy_bulk.hswath.values)
#    tp_max = np.nanmax(xs_shy_bulk.tswath.values)
#    cmapmax = xs_shy_bulk.case.size
#
#    # figure
#    fig = plt.figure(figsize=[width,height], dpi=200)
#
#    # column nº1: vortex, hsbmu, tpbmu, bmu
#    left, right = 0, .15
##     gs, im1 = aux_wind(fig, xs_wind, st, swath=True, vmax=wind_max, 
##                        cmap=cmap_wind, rect=[left, .75, right, 1])
#    gs, im2 = aux_var_shy(fig, mesh_lo, mesh_la, xs_shy_bulk, xs_pts, st, 
#                          var='hsbmu', swath=True, cmap=cmap_hs, vmax=hs_max, 
#                          text=text, rect=[left, .66, right, 1])
#    gs, im3 = aux_var_shy(fig, mesh_lo, mesh_la, xs_shy_bulk, xs_pts, st, 
#                          var='tpbmu', swath=True, cmap=cmap_tp, vmax=tp_max, 
#                          text=text, rect=[left, .33, right, .66])
#    gs, im4 = aux_var_shy(fig, mesh_lo, mesh_la, xs_shy_bulk, xs_pts, st, 
#                          var='bmu', swath=True, cmap=cmap_bmu, vmax=cmapmax, 
#                          text=text, rect=[left, 0, right, .33])
#
#    # column nº2: colorbars
##     aux_cbar_var(fig, im1, ttl='W [m/s]', rect=[right, .75, .2, 1-.02],
##                  xlabel=True)
#    aux_cbar_var(fig, im2, ttl='Hs [m]', rect=[right, .66, .2, 1-.02],
#                 xlabel=True)
#    aux_cbar_var(fig, im3, ttl='Tp [s]', rect=[right, .33, .2, .66-.02],
#                 xlabel=True)
#    aux_cbar_var(fig, im4, ttl='bmu []', rect=[right, 0, .2, .33-.02],
#                 xlabel=True)
#
#    # row nº3: series
#    pargmin = [np.argmin(np.abs(xs_shy_bulk.lon.values - lonpts[i]) + \
#                         np.abs(xs_shy_bulk.lat.values - latpts[i])) for i in range(len(lonpts))]
#    
#    gs, im11 = aux_series_swath(fig, xs_pts, pargmin, xs_shy_bulk, numpt=0,
#                                cmap=cmap_bmu, rect=[.2, .5, .6, 1])
#    gs, im11 = aux_series_swath(fig, xs_pts, pargmin, xs_shy_bulk, numpt=1,
#                                cmap=cmap_bmu, rect=[.2, 0, .6, .5])
#    gs, im11 = aux_series_swath(fig, xs_pts, pargmin, xs_shy_bulk, numpt=2,
#                                cmap=cmap_bmu, rect=[.6, .5, 1, 1])
#    gs, im11 = aux_series_swath(fig, xs_pts, pargmin, xs_shy_bulk, numpt=3,
#                                cmap=cmap_bmu, rect=[.6, 0, 1, .5])
#    
##    aux_cbar_bmu(fig, im11, ttl='', rect=[.48, 0, .53, .66], 
##                 num_ticks=xs_shy_bulk.case.size)
#
#    return fig
