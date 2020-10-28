#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import pandas as pd
import numpy as np
from scipy import interpolate
from math import sqrt

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.dates as mdates
import matplotlib.colors as colors
from matplotlib.colors import ListedColormap
from matplotlib import cm

#from .storms import get_category

def get_storm_color(categ):

    dcs = {
        0 : 'lime',
        1 : 'yellow',
        2 : 'orange',
        3 : 'red',
        4 : 'purple',
        5 : 'black',
        6 : 'gray',
    }

    return dcs[categ]


def plot_trackcategory_basemap(lo0, lo1, la0, la1, lo_tc, la_tc, categ, x0, y0, 
                               x_area, y_area, name, storm_year, target, zoom=False):
    'Plot Historical TCs tracks map, requires basemap module'

    try:
        from mpl_toolkits.basemap import Basemap
    except:
        print('basemap module required.')
        return

    # general map figure
    fig = plt.figure(figsize = (20,9))
    
    # map bound limits
    if zoom:    res = 'h'
    else:       res = 'l'
    map = Basemap(llcrnrlon = lo0, llcrnrlat = la0, urcrnrlon = lo1, urcrnrlat = la1,
                 resolution = res, projection = 'cyl', lat_0 = la0, lon_0 = lo0, area_thresh=0.01)
    
    # map format
    map.drawcoastlines()
    map.fillcontinents(color='darkgray')
    map.drawmapboundary(fill_color='lightcyan')
    map.drawparallels(np.arange(-60,60,20),labels=[1,1,0,1])
    map.drawmeridians(np.arange(-180,180,20),labels=[1,1,0,1])
    
    # plot storm track
    if zoom:    ini, num = 0, np.shape(lo_tc)[0]-1
    else:       ini, num = 1, np.shape(lo_tc)[0]
    
#    for j in list(range(1,np.shape(lo_tc)[0])):
    for j in list(range(ini,num)):
        if zoom:
            map.plot([lo_tc[j], lo_tc[j+1]], [la_tc[j], la_tc[j+1]],'-', 
                     color=get_storm_color(categ[j+1]), linewidth=3)
            map.plot(lo_tc[j], la_tc[j],'.', color=get_storm_color(categ[j]), markersize=12)
        else:
            map.plot([lo_tc[j-1], lo_tc[j]], [la_tc[j-1], la_tc[j]],'-', 
                     color=get_storm_color(categ[j-1]), linewidth=1)
            map.plot(lo_tc[j], la_tc[j],'.', color=get_storm_color(categ[j-1]), linewidth=1)
    
    # plot target 
    if zoom:    map.plot(x0, y0,'+', color='dodgerblue', mew=4, ms=22, label=target)
    else:       map.plot(x0, y0,'+', color='dodgerblue', mew=3, ms=15, label=target)
    
    # plot numerical domain
    map.plot(x_area, y_area,'-', color='dodgerblue', linewidth=3, markersize=12, label='Target area')
    
    # category labels
    map.plot([],[],'.-', color='lime', linewidth=3, markersize=12, label='Cat 0')
    map.plot([],[],'.-', color='yellow', linewidth=3, markersize=12, label='Cat 1')
    map.plot([],[],'.-', color='orange', linewidth=3, markersize=12, label='Cat 2')
    map.plot([],[],'.-', color='red', linewidth=3, markersize=12, label='Cat 3')
    map.plot([],[],'.-', color='darkmagenta', linewidth=3, markersize=12, label='Cat 4')
    map.plot([],[],'.-', color='black', linewidth=3, markersize=12, label='Cat 5')
    map.plot([],[],'.-', color='gray', linewidth=3, markersize=12, label='Unknown')
    
    plt.legend(loc='lower right', framealpha=0.4)
    plt.title('Target location & {0} track data ({1})'.format(name, storm_year), 
              fontsize=14, fontweight='bold')
    plt.show()
    
    return fig


def plot_trackcategory_basemap_list(lo0, lo1, la0, la1, lo_tc_ls, la_tc_ls, categ_ls, x0, y0, 
                                    x_area, y_area, year1, year2, target, zoom=False):
    'Plot Historical TCs tracks map, requires basemap module'

    try:
        from mpl_toolkits.basemap import Basemap
    except:
        print('basemap module required.')
        return

    # general map figure
    fig = plt.figure(figsize = (20,9))
    
    # map bound limits
    if zoom:    res = 'h'
    else:       res = 'l'
    map = Basemap(llcrnrlon = lo0, llcrnrlat = la0, urcrnrlon = lo1, urcrnrlat = la1,
                 resolution = res, projection = 'cyl', lat_0 = la0, lon_0 = lo0, area_thresh=0.01)
    
    # map format
    map.drawcoastlines()
    map.fillcontinents(color='darkgray')
    map.drawmapboundary(fill_color='lightcyan')
    map.drawparallels(np.arange(-60,60,20),labels=[1,1,0,1])
    map.drawmeridians(np.arange(-180,180,20),labels=[1,1,0,1])
    
    for i in range(len(lo_tc_ls)):
        lo_tc, la_tc, categ = lo_tc_ls[i], la_tc_ls[i], categ_ls[i]
        
        # plot storm track
        if zoom:    ini, num = 0, np.shape(lo_tc)[0]-1
        else:       ini, num = 1, np.shape(lo_tc)[0]
        
    #    for j in list(range(1,np.shape(lo_tc)[0])):
        for j in list(range(ini,num)):
            if zoom:
                map.plot([lo_tc[j], lo_tc[j+1]], [la_tc[j], la_tc[j+1]],'-', 
                         color=get_storm_color(categ[j+1]), linewidth=3)
                map.plot(lo_tc[j], la_tc[j],'.', color=get_storm_color(categ[j]), markersize=12)
            else:
                map.plot([lo_tc[j-1], lo_tc[j]], [la_tc[j-1], la_tc[j]],'-', 
                         color=get_storm_color(categ[j-1]), linewidth=1)
                map.plot(lo_tc[j], la_tc[j],'.', color=get_storm_color(categ[j-1]), linewidth=1)
    
    # plot target 
    if zoom:    map.plot(x0, y0,'+', color='dodgerblue', mew=4, ms=22, label=target)
    else:       map.plot(x0, y0,'+', color='dodgerblue', mew=3, ms=15, label=target)
    
    # plot numerical domain
    map.plot(x_area, y_area,'-', color='dodgerblue', linewidth=3, markersize=12, label='Target area')
    
    # category labels
    map.plot([],[],'.-', color='lime', linewidth=3, markersize=12, label='Cat 0')
    map.plot([],[],'.-', color='yellow', linewidth=3, markersize=12, label='Cat 1')
    map.plot([],[],'.-', color='orange', linewidth=3, markersize=12, label='Cat 2')
    map.plot([],[],'.-', color='red', linewidth=3, markersize=12, label='Cat 3')
    map.plot([],[],'.-', color='darkmagenta', linewidth=3, markersize=12, label='Cat 4')
    map.plot([],[],'.-', color='black', linewidth=3, markersize=12, label='Cat 5')
    map.plot([],[],'.-', color='gray', linewidth=3, markersize=12, label='Unknown')
    
    plt.legend(loc='lower right', framealpha=0.4)
    plt.title('Target location & track data ({0}-{1})'.format(year1, year2), 
              fontsize=14, fontweight='bold')
    plt.show()
    
    return fig


def plot_bathymetry(xds_bathy, vmin, vmax, x0, y0, np_shore, target, zoom=False):
    
    fig = plt.figure(figsize=(12,8))
    
    # plot depth, target, shoreline
    xds_bathy.elevation.plot(vmin=vmin, vmax=vmax, cmap='gist_earth')
    plt.plot(np_shore[:,0], np_shore[:,1], color='crimson')
    plt.plot(x0, y0, '+', color='k', mew=4, ms=22)
    
    plt.gca().set_aspect('equal', adjustable='box')
    plt.gca().set_xlabel('Longitude [º]', fontweight='bold')
    plt.gca().set_ylabel('Latitude [º]', fontweight='bold')
    
    # title
    if zoom: plt.title('GEBCO bathymetry (450m) - zoom {0}'.format(
                        target), fontsize=16, fontweight='bold')
    else:    plt.title('GEBCO bathymetry (450m resolution) {0}'.format(
                        target), fontsize=14, fontweight='bold')
    
    return fig

def custom_cmap(numcolors, map1, m1ini, m1end, map2, m2ini, m2end):
    '''
    Generate custom colormap
    Example: Red-Orange-Yellow-Green-Blue -- map1='YlOrRd' map2='YlGnBu_r'
    mXini, mXend:   colormap range of colors 
    numcolors:      number of colors (100-continuous, 15-discretization)
    '''
    
    # color maps
    cmap1 = plt.get_cmap(map1, numcolors)
    cmap2 = plt.get_cmap(map2, numcolors)
    
    # custom color ranges
    cmap1v = colors.LinearSegmentedColormap.from_list(
            'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap1.name, a=m1ini, b=m1end),
            cmap1(np.linspace(m1ini,m1end,100)))
    cmap2v = colors.LinearSegmentedColormap.from_list(
            'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap2.name, a=m2ini, b=m2end),
            cmap2(np.linspace(m2ini,m2end,100)))
    
    top = cm.get_cmap(cmap1v, 128)
    bottom = cm.get_cmap(cmap2v, 128)
    
    newcolors = np.vstack((bottom(np.linspace(0,1,128)),top(np.linspace(0,1,128))))
    newcmp = ListedColormap(newcolors, name='OrangeBlue')
    
    return newcmp


def plot_vortex_times(name, xds_out_case, var_name, p_export_case, num_nests, var_max, 
                      var_min, lon, lat, quiver=False, np_shore=np.array([]), cmap='jet'):
    '''
    Plots non-stationary SWAN execution output for selected var and case

    name       - name of TC
    xds_out_case - swan output (xarray.Dataset) 
    var_name   - 'Hsig', 'Tm02', 'Tpsmoo'
    p_export   - path for exporting figures
    num_nests  - number of nested grids
    lon, lat   - nested grids boundary limits

    opt. args
    quiver     - True for adding directional quiver plot
    np_shore   - shoreline, np.array x = np_shore[:,0] y = np.shore[:,1]
    cmap       - matplotlib colormap
    '''

    if not op.isdir(p_export_case): os.makedirs(p_export_case)

    # iterate case output over time
    for t in xds_out_case.time.values[:]:
        xds_oct = xds_out_case.sel(time=t)

        # time string
        t_str = pd.to_datetime(str(t)).strftime('%d-%b-%Y_%H:%M%p')

        # get mesh data from output dataset
        X = xds_oct.X.values[:]
        Y = xds_oct.Y.values[:]
        delta_x = X[-1]-X[0]
        num_quiver = 60  # number of quiver in x-axis
        step = delta_x/num_quiver
        scale_quiver = 1/step * var_max
#        scale_quiver = 2 * var_max

        # get variable and units
        var = xds_oct[var_name].values[:]
        var_units = xds_oct[var_name].attrs['units']

        # new figure
        fig, ax0 = plt.subplots(nrows=1, figsize=(12, 12))
        var_title = '{0}'.format(var_name)  # title

        # pcolormesh
        newcmap = custom_cmap(100, 'plasma_r', 0.05, 0.9, 'viridis', 0.2, 1)
        im = plt.pcolormesh(X, Y, var, cmap=newcmap, vmin= var_min, vmax= var_max)

        # add quiver plot 
        if quiver:
            var_title += '-Dir'
            var_dir = xds_oct.Dir.values[:]
            aux_quiver(X, Y, var, var_dir, scale_quiver, step, ax0)
            
            # try using streamplot arrows
#            speed = np.sqrt(var*np.sin(np.deg2rad(var_dir))**2 + var*np.cos(np.deg2rad(var_dir))**2)
#            lw = 5 * speed / speed.max() 
#            ax0.streamplot(X, Y, -var*np.sin(np.deg2rad(var_dir)), -var*np.cos(np.deg2rad(var_dir)), density=1, color='b', linewidth=lw) 

        # shoreline
        if np_shore.any():
            x_shore = np_shore[:,0]
            y_shore = np_shore[:,1]
            plt.plot(x_shore, y_shore,'.', color='dimgray', markersize=3)

        # plot NESTED meshgrids
        for i in np.arange(0,num_nests):
            plt.plot([lon[i][0], lon[i][0]], [lat[i][0], lat[i][1]], '-k', linewidth=2, label='')
            plt.plot([lon[i][1], lon[i][1]], [lat[i][0], lat[i][1]], '-k', linewidth=2, label='')
            plt.plot([lon[i][0], lon[i][1]], [lat[i][0], lat[i][0]], '-k', linewidth=2, label='')
            plt.plot([lon[i][0], lon[i][1]], [lat[i][1], lat[i][1]], '-k', linewidth=2, label='')

        # customize pcolormesh
        plt.title('{0}     {1}'.format(name, t_str), fontsize = 14, fontweight='bold')
        plt.xlabel('Longitude (º)', fontsize = 12)
        plt.ylabel('Latitude (º)', fontsize = 12)

        plt.axis('scaled')
        plt.xlim(X[0], X[-1])
        plt.ylim(Y[0], Y[-1])

        # add custom colorbar
        divider = make_axes_locatable(ax0)
        cax = divider.append_axes("right", size="5%", pad=0.1)
        plt.colorbar(im, cax=cax)
        plt.ylabel('{0} ({1})'.format(var_name, var_units), fontsize = 12)

        # export fig
        p_ex = op.join(p_export_case, 'outmap_{0}_{1}.png'.format(var_name, t_str))
        fig.savefig(p_ex)

        # close fig 
        plt.close()

def plot_output_nonstat(name, xds_out, var_name, p_export, num_nests,  
                        lon, lat, quiver=False, np_shore=np.array([])):
    '''
    Plots non-stationary SWAN execution output for selected var, for every case

    name       - name of TC
    xds_out    - swan output (xarray.Dataset)
    var_name   - 'Hsig', 'Tm02', 'Tpsmoo'
    p_export   - path for exporting figures
    num_nests  - number of nested grids
    lon, lat   - nested grids boundary limits

    opt. args
    quiver     - True for adding directional quiver plot
    np_shore   - shoreline, np.array x = np_shore[:,0] y = np.shore[:,1]
    '''

    # make export dir
    if not op.isdir(p_export): os.makedirs(p_export)

    for case_ix in xds_out.case.values[:]:

        # select case
        xds_out_case = xds_out.sel(case=case_ix)

        # output case subfolder
        case_id = '{0:04d}'.format(case_ix)
        p_export_case = op.join(p_export, case_id)

        # plot variable times
        plot_var_times(
            name, xds_out_case, var_name, p_export_case,
            num_nests, lon, lat, quiver=quiver, np_shore=np_shore)

def plot_var_times(name, xds_out_case, var_name, p_export_case, num_nests, lon, lat, 
                   quiver=False, np_shore=np.array([]), cmap='jet'):
    '''
    Plots non-stationary SWAN execution output for selected var and case

    name       - name of TC
    xds_out_case - swan output (xarray.Dataset) 
    var_name   - 'Hsig', 'Tm02', 'Tpsmoo'
    p_export   - path for exporting figures
    num_nests  - number of nested grids
    lon, lat   - nested grids boundary limits

    opt. args
    quiver     - True for adding directional quiver plot
    np_shore   - shoreline, np.array x = np_shore[:,0] y = np.shore[:,1]
    cmap       - matplotlib colormap
    '''

    if not op.isdir(p_export_case): os.makedirs(p_export_case)

    # maximum value of output variable over time (GENERAL MESH) --->  axis = (case, time, Y, X)
    var_max = float(xds_out_case[var_name].max(axis=0).max())
    var_min = float(xds_out_case[var_name].min(axis=0).min())

    # iterate case output over time
    for t in xds_out_case.time.values[:]:
        xds_oct = xds_out_case.sel(time=t)

        # time string
        t_str = pd.to_datetime(str(t)).strftime('%d-%b-%Y %H:%M%p') 

        # get mesh data from output dataset
        X = xds_oct.X.values[:]
        Y = xds_oct.Y.values[:]
        delta_x = X[-1]-X[0]
        num_quiver = 60  # number of quiver in x-axis
        step = delta_x/num_quiver
        scale_quiver = 1/step * var_max
#        scale_quiver = 2 * var_max

        # get variable and units
        var = xds_oct[var_name].values[:]
        var_units = xds_oct[var_name].attrs['units']

        # new figure
        fig, ax0 = plt.subplots(nrows=1, figsize=(12, 12))
        var_title = '{0}'.format(var_name)  # title

        # pcolormesh
        newcmap = custom_cmap(100, 'YlOrRd', 0.09, 0.9, 'YlGnBu_r', 0, 0.88)
        im = plt.pcolormesh(X, Y, var, cmap=newcmap, vmin= var_min, vmax= var_max)

        # add quiver plot 
        if quiver:
            var_title += '-Dir'
            var_dir = xds_oct.Dir.values[:]
            aux_quiver(X, Y, var, var_dir, scale_quiver, step, ax0)

        # shoreline
        if np_shore.any():
            x_shore = np_shore[:,0]
            y_shore = np_shore[:,1]
            plt.plot(x_shore, y_shore,'.', color='dimgray', markersize=3)

        # plot NESTED meshgrids
        for i in np.arange(0,num_nests):
            plt.plot([lon[i][0], lon[i][0]], [lat[i][0], lat[i][1]], '-k', linewidth=2, label='')
            plt.plot([lon[i][1], lon[i][1]], [lat[i][0], lat[i][1]], '-k', linewidth=2, label='')
            plt.plot([lon[i][0], lon[i][1]], [lat[i][0], lat[i][0]], '-k', linewidth=2, label='')
            plt.plot([lon[i][0], lon[i][1]], [lat[i][1], lat[i][1]], '-k', linewidth=2, label='')

        # customize pcolormesh
#        plt.title('{0} (t={1})'.format(var_title, t_str),
        plt.title('{0}     {1}'.format(name, t_str),
                  fontsize = 14, fontweight='bold')
        plt.xlabel(xds_oct.attrs['xlabel'], fontsize = 12)
        plt.ylabel(xds_oct.attrs['ylabel'], fontsize = 12)

        plt.axis('scaled')
        plt.xlim(X[0], X[-1])
        plt.ylim(Y[0], Y[-1])

        # add custom colorbar
        divider = make_axes_locatable(ax0)
        cax = divider.append_axes("right", size="5%", pad=0.1)
        plt.colorbar(im, cax=cax)
        plt.ylabel('{0} ({1})'.format(var_name, var_units), fontsize = 12)

        # export fig
        p_ex = op.join(p_export_case, 'outmap_{0}_{1}.png'.format(var_name, t_str))
        fig.savefig(p_ex)

        # close fig 
        plt.close()
        
def aux_quiver(X, Y, var, vdir, scale_quiver, step, ax=None):
    '''
    interpolates var and plots quiver with var_dir. Requires open figure

    var  - variable module
    vdir - variable direction (º clockwise relative to North)
    scale_quiver - scale for proportional size
    '''
    ax=ax
#    size = 30  # quiver mesh size
    # quiver mesh definition (proportional in both axis)
#    step = 0.25  # quiver each 0.5º
    degree_x = X[-1]-X[0]    
    degree_y = Y[-1]-Y[0]  
    size_x = round(degree_x / step)
    size_y = round(degree_y / step)

    # var and dir interpolators 
    vdir_f = vdir.copy()
    vdir_f[np.isnan(vdir_f)] = 0
#    f_dir = interpolate.interp2d(X, Y, vdir_f, kind='linear')
    vdir_f_y = np.cos(np.deg2rad(vdir_f))
    vdir_f_x = np.sin(np.deg2rad(vdir_f))
    f_dir_x = interpolate.interp2d(X, Y, vdir_f_x, kind='linear')
    f_dir_y = interpolate.interp2d(X, Y, vdir_f_y, kind='linear')
    
    var_f = var.copy()
    var_f[np.isnan(var_f)] = 0
    f_var = interpolate.interp2d(X, Y, var_f, kind='linear')

    # generate quiver mesh
    x_q = np.linspace(X[0], X[-1], num = size_x)
    y_q = np.linspace(Y[0], Y[-1], num = size_y)

    # interpolate data to quiver mesh
    var_q = f_var(x_q, y_q)
#    vdir_q = f_dir(x_q, y_q)
    vdir_q_x = f_dir_x(x_q, y_q)
    vdir_q_y = f_dir_y(x_q, y_q)
    vdir_q = np.rad2deg(np.arctan(vdir_q_x/vdir_q_y))
    # sign correction
    vdir_q[(vdir_q_x>0) & (vdir_q_y<0)] = vdir_q[(vdir_q_x>0) & (vdir_q_y<0)] + 180
    vdir_q[(vdir_q_x<0) & (vdir_q_y<0)] = vdir_q[(vdir_q_x<0) & (vdir_q_y<0)] + 180  

    # u and v dir components
    u = np.sin(np.deg2rad(vdir_q))
    v = np.cos(np.deg2rad(vdir_q))

    # plot quiver
    ax.quiver(x_q, y_q, -u*var_q, -v*var_q, width=0.0015, scale=scale_quiver, scale_units='x')#width=0.003, scale=1, scale_units='inches'

    return ax

def plot_output_points(name, xds_out, p_export):
    '''
    Plots SWAN execution output table points time series

    name       - name of TC
    xds_out    - swan points output (xarray.Dataset)
    p_export   - path for exporting figures
    '''

    # make export dir
    if not op.isdir(p_export): os.makedirs(p_export)

    for case_ix in xds_out.case.values[:]:

        # select case
        xds_out_case = xds_out.sel(case=case_ix)

        # output case subfolder
        case_id = '{0:04d}'.format(case_ix)
        p_export_case = op.join(p_export, case_id)

        # plot variable times
        plot_points_times(name, xds_out_case, p_export_case)

def plot_points_times(name, xds_out_case, p_export_case):
    '''
    Plots non-stationary SWAN points output for selected case

    name           - name of TC
    xds_out_case   - swan case output (xarray.Dataset)
    p_export_case  - path for exporting figures
    '''

    if not op.isdir(p_export_case): os.makedirs(p_export_case)

    # get initial and final dates (in datetime.datetime format)
    dates = xds_out_case.time.values[:]    
    date_0 = xds_out_case.time.values[0]
    date_1 = xds_out_case.time.values[-1]
    t_str_ini = pd.to_datetime(str(date_0)).strftime('%d-%b-%Y %H:%M%p')
    t_str_fin = pd.to_datetime(str(date_1)).strftime('%d-%b-%Y %H:%M%p')  

    # iterate over points
    n_pts = len(xds_out_case.point)
    for i in range(n_pts):

        # get point variables
        xds_pt = xds_out_case.isel(point=i)

        hs = xds_pt.HS.values[:]
        tm = xds_pt.TM02.values[:]
        tp = xds_pt.RTP.values[:]
        dr = xds_pt.DIR.values[:]

        # plot and save figure series of each output point
        fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(12, 12), sharex=True)
        fig.suptitle('Storm {0}   POINT nº{1}   ({2})--({3})'.format(
                name, i, t_str_ini, t_str_fin), fontweight = 'bold', fontsize = 14)
        
        ax1.plot(dates, hs, '.', color = 'b', markersize=2, label="Hs [m]")
        ax1.set_xlim(dates[0], dates[-1])
        plt.setp(ax1.get_xticklabels(), visible=False)
        ax1.set_title('Significant Wave Height [m]', fontweight = 'bold')

        ax2.plot(dates, tm, '.', color = 'b', markersize=2)
        plt.setp(ax2.get_xticklabels(), visible=False)
        ax2.set_title('Mean Period [s]', fontweight = 'bold')

        ax3.plot(dates, tp, '.', color = 'b', markersize=2)
        plt.setp(ax3.get_xticklabels(), visible=False)
        ax3.set_title('Peak Period [s]', fontweight = 'bold')

        ax4.plot(dates, dr, '.', color = 'b', markersize=2)
        ax4.set_ylim([0, 360])
        ax4.xaxis.set_major_formatter(mdates.DateFormatter('%d/%b %Hh'))
        plt.setp(ax4.get_xticklabels(), visible=True)
        ax4.set_title('Wave direction [º]', fontweight = 'bold')
        
        fig.autofmt_xdate()

        # export fig
        p_ex = op.join(p_export_case, 'point_{0}.png'.format(i))
        fig.savefig(p_ex)

        # close fig 
        plt.close()

def plot_pts_grids(name, xds_out, p_export, show=False, vortex=False):
    '''
    Plots SWAN execution output table points time series

    name       - name of TC
    xds_out    - swan points output (xarray.Dataset)
    p_export   - path for exporting figures
    '''

    # make export dir
    if not op.isdir(p_export): os.makedirs(p_export)

    for case_ix in xds_out[0].case.values[:]:

#        # select case
#        ls_xds_out_case = [] 
#        for j in range(len(xds_out)):
#            ls_xds_out_case.append(xds_out[j].sel(case=case_ix))

        # output case subfolder
        case_id = '{0:04d}'.format(case_ix)
        p_export_case = op.join(p_export, case_id)

        # plot variable times
#        plot_grids_times(name, ls_xds_out_case, p_export_case)
        if vortex:  plot_grids_times_vortex(name, xds_out, p_export_case, show)
        else:       plot_grids_times(name, xds_out, p_export_case, show)

def plot_grids_times(name, xds_out_case, p_export_case, show=False):
    '''
    Plots non-stationary SWAN points output for selected case

    name           - name of TC
    xds_out_case   - swan case output (xarray.Dataset)
    p_export_case  - path for exporting figures
    '''

    if not op.isdir(p_export_case): os.makedirs(p_export_case)

    # get initial and final dates (in datetime.datetime format)
    dates = xds_out_case[0].time.values[:]    
    date_0 = xds_out_case[0].time.values[0]
    date_1 = xds_out_case[0].time.values[-1]
    t_str_ini = pd.to_datetime(str(date_0)).strftime('%d-%b-%Y %H:%M%p')
    t_str_fin = pd.to_datetime(str(date_1)).strftime('%d-%b-%Y %H:%M%p')  
    
    # get variables y-limit for plot
    n_pts = len(xds_out_case[0].point)
    hs_max, tm_max, tp_max = [], [], []
    for j in np.arange(0, len(xds_out_case)):
        hs_max.append(xds_out_case[j].HS.values.max())
        tm_max.append(xds_out_case[j].TM02.values.max())
        tp_max.append(xds_out_case[j].RTP.values.max())    
    hs_max = np.nanmax(hs_max)
    tm_max = np.nanmax(tm_max)
    tp_max = np.nanmax(tp_max)
    
    # iterate over points
    for i in range(n_pts):

        # get point variables
        hs, tm, tp, dr, dspr, wx = [], [], [], [], [], []
        hs = np.empty((len(xds_out_case),len(xds_out_case[0].time)))
        tm = np.empty((len(xds_out_case),len(xds_out_case[0].time)))
        tp = np.empty((len(xds_out_case),len(xds_out_case[0].time)))
        dr = np.empty((len(xds_out_case),len(xds_out_case[0].time)))
        dspr = np.empty((len(xds_out_case),len(xds_out_case[0].time)))
        #wx = np.empty((len(xds_out_case),len(xds_out_case[0].time)))
        for j in np.arange(0, len(xds_out_case)):
            hs[j,:] = xds_out_case[j].isel(point=i).HS.values[:][0]
            tm[j,:] = xds_out_case[j].isel(point=i).TM02.values[:][0]
            tp[j,:] = xds_out_case[j].isel(point=i).RTP.values[:][0]
            dr[j,:] = xds_out_case[j].isel(point=i).DIR.values[:][0]
            dspr[j,:] = xds_out_case[j].isel(point=i).DSPR.values[:][0]
            #wx[j,:] = xds_out_case[j].isel(point=i).WIND.values[:][0]

            
        # plot and save figure series of each output point
        fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, 1, figsize=(20, 10), sharex=True)
        fig.suptitle('{0}   POINT nº{1}   ({2})--({3})'.format(
                name, i, t_str_ini, t_str_fin), fontweight = 'bold', fontsize = 14)
        
        c = ['k', 'r', 'deepskyblue', 'gold', 'lime']
#        label = ["Grid_0", "Nest_1", "Nest_2", "Nest_3", "Nest_4"]
        label = ["15 km", "5 km", "1 km", "250 m", "75 m"]
        for j, xds in enumerate(zip(xds_out_case)):       
            ax1.plot(dates, hs[j,:], '-', color = c[j], markersize=2, label=label[j])    
            ax2.plot(dates, tm[j,:], '-', color = c[j], markersize=2, label=label[j])    
            ax3.plot(dates, tp[j,:], '-', color = c[j], markersize=2, label=label[j])    
            ax4.plot(dates, dr[j,:], '-', color = c[j], markersize=2, label=label[j])
            ax5.plot(dates, dspr[j,:], '-', color = c[j], markersize=2, label=label[j])
            #ax6.plot(dates, wx[j,:], '-', color = c[j], markersize=2, label=label[j])
         
        ax1.set_xlim(dates[0], dates[-1])
#        ax1.set_ylim(0, hs_max+0.5)
        plt.setp(ax1.get_xticklabels(), visible=False)
        ax1.set_title('Significant Wave Height [m]')
        ax1.legend()
        
#        ax2.set_ylim(0, tp_max+0.5)
        plt.setp(ax2.get_xticklabels(), visible=False)
        ax2.set_title('Mean Period [s]')

#        ax3.set_ylim(0, tp_max+0.5)
        plt.setp(ax3.get_xticklabels(), visible=False)
        ax3.set_title('Peak Period [s]')
        
        ax4.set_ylim([0, 360])
        plt.setp(ax4.get_xticklabels(), visible=False)
        ax4.set_title('Wave direction [º]')

#        ax5.set_ylim(0, dspr_max+0.5)
        plt.setp(ax5.get_xticklabels(), visible=True)
        ax5.set_title('Directional spread ')

#        ax6.set_ylim(0, wx_max+0.5)
        ax5.xaxis.set_major_formatter(mdates.DateFormatter('%d/%b %Hh'))
        #plt.setp(ax6.get_xticklabels(), visible=True)
        #ax6.set_title('Wind x-axis [m/s]')

        fig.autofmt_xdate()

        # export fig
        p_ex = op.join(p_export_case, 'point_{0}.png'.format(i))
        fig.savefig(p_ex)
        
        if show:    plt.show()

        # close fig 
        plt.close()

def plot_grids_times_vortex(name, xds_out_case, p_export_case, show=False):
    '''
    Plots non-stationary SWAN points output for selected case

    name           - name of TC
    xds_out_case   - swan case output (xarray.Dataset)
    p_export_case  - path for exporting figures
    '''

    if not op.isdir(p_export_case): os.makedirs(p_export_case)

    # get initial and final dates (in datetime.datetime format)
    dates = xds_out_case[0].time.values[:]    
    date_0 = xds_out_case[0].time.values[0]
    date_1 = xds_out_case[0].time.values[-1]
    t_str_ini = pd.to_datetime(str(date_0)).strftime('%d-%b-%Y %H:%M%p')
    t_str_fin = pd.to_datetime(str(date_1)).strftime('%d-%b-%Y %H:%M%p')  
    
    # get variables y-limit for plot
    n_pts = len(xds_out_case[0].point)
    
    # iterate over points
    for i in range(n_pts):

        # get point variables
        w, wd = [], []
        w = np.empty((len(xds_out_case),len(xds_out_case[0].time)))
        wd = np.empty((len(xds_out_case),len(xds_out_case[0].time)))
        
        for j in np.arange(0, len(xds_out_case)):
            w[j,:] = xds_out_case[j].isel(point=i).W.values[:][0]
            wd[j,:] = xds_out_case[j].isel(point=i).Dir.values[:][0]
            
        # plot and save figure series of each output point
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(20, 4), sharex=True)
        fig.suptitle('{0}   POINT nº{1}   ({2})--({3})'.format(
                name, i, t_str_ini, t_str_fin), fontweight = 'bold', fontsize = 14)
        
        c = ['k', 'r', 'deepskyblue', 'gold', 'lime']
        for j, xds in enumerate(zip(xds_out_case)):       
            ax1.plot(dates, w[j,:], '-', color = c[j], markersize=2)    
            ax2.plot(dates, wd[j,:], '-', color = c[j], markersize=2)    
         
        ax1.set_xlim(dates[0], dates[-1])
        plt.setp(ax1.get_xticklabels(), visible=False)
        ax1.set_title('Wind speed [m/s]')
        
        ax2.set_ylim([0, 360])
        ax2.xaxis.set_major_formatter(mdates.DateFormatter('%d/%b %Hh'))
        plt.setp(ax2.get_xticklabels(), visible=True)
        ax2.set_title('Wind direction [º]')

        fig.autofmt_xdate()

        # export fig
        p_ex = op.join(p_export_case, 'point_vortex_{0}.png'.format(i))
        fig.savefig(p_ex)
        
        if show:    plt.show()

        # close fig 
        plt.close()

def plot_pts_mda(name, xds_out, p_export):
    '''
    Plots SWAN execution output table points time series

    name       - name of TC
    xds_out    - swan points output (xarray.Dataset)
    p_export   - path for exporting figures
    '''

    # make export dir
    if not op.isdir(p_export): os.makedirs(p_export)

    for case_ix in xds_out[0].case.values[:]:

        # select case
        xds_out_case = [xds_out[0].sel(case=case_ix)]

        # output case subfolder
        case_id = '{0:04d}'.format(case_ix)
        p_export_case = op.join(p_export, case_id)

        # plot variable times
        plot_pts_mda_times(name, xds_out_case, p_export_case)

def plot_pts_mda_times(name, xds_out_case, p_export_case):
    '''
    Plots non-stationary SWAN points output for selected case

    name           - name of TC
    xds_out_case   - swan case output (xarray.Dataset)
    p_export_case  - path for exporting figures
    '''

    if not op.isdir(p_export_case): os.makedirs(p_export_case)

    # get initial and final dates (in datetime.datetime format)
    dates = xds_out_case[0].time.values[:]    
    date_0 = xds_out_case[0].time.values[0]
    date_1 = xds_out_case[0].time.values[-1]
    t_str_ini = pd.to_datetime(str(date_0)).strftime('%d-%b-%Y %H:%M%p')
    t_str_fin = pd.to_datetime(str(date_1)).strftime('%d-%b-%Y %H:%M%p')  
    
    # get variables y-limit for plot
    n_pts = len(xds_out_case[0].point)
    for i in range(n_pts):
        # get point variables
        hs, tm, tp, dr = [], [], [], []
        hs = np.empty((n_pts,len(xds_out_case[0].time)))
        tm = np.empty((n_pts,len(xds_out_case[0].time)))
        tp = np.empty((n_pts,len(xds_out_case[0].time)))
        dr = np.empty((n_pts,len(xds_out_case[0].time)))
        hs_max, tm_max, tp_max = [], [], []
        for j in np.arange(0, len(xds_out_case)):
            hs[j,:] = xds_out_case[j].sel(point=i).HS.values[:][0]
            tm[j,:] = xds_out_case[j].sel(point=i).TM02.values[:][0]
            tp[j,:] = xds_out_case[j].sel(point=i).RTP.values[:][0]
            dr[j,:] = xds_out_case[j].sel(point=i).DIR.values[:][0]
        
        # variables y-limit for plot
        hs_max.append(hs.max())
        tm_max.append(tm.max())
        tp_max.append(tp.max())
    
    hs_max = max(hs_max)
    tm_max = max(tm_max)
    tp_max = max(tp_max)
    
    # iterate over points
    for i in range(n_pts):

        # get point variables
        hs, tm, tp, dr = [], [], [], []
        hs = np.empty((n_pts,len(xds_out_case[0].time)))
        tm = np.empty((n_pts,len(xds_out_case[0].time)))
        tp = np.empty((n_pts,len(xds_out_case[0].time)))
        dr = np.empty((n_pts,len(xds_out_case[0].time)))
        for j in np.arange(0, len(xds_out_case)):
            hs[j,:] = xds_out_case[j].sel(point=i).HS.values[:][0]
            tm[j,:] = xds_out_case[j].isel(point=i).TM02.values[:][0]
            tp[j,:] = xds_out_case[j].isel(point=i).RTP.values[:][0]
            dr[j,:] = xds_out_case[j].isel(point=i).DIR.values[:][0]
            
#        # remove nan
#        hs = hs[~np.isnan(hs)]
#        tm = tm[~np.isnan(tm)]
#        tp = tp[~np.isnan(tp)]
#        dr = dr[~np.isnan(dr)]
            
        # plot and save figure series of each output point
        fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(12, 12), sharex=True)
        fig.suptitle('Storm {0}   POINT nº{1}   ({2})--({3})'.format(
                name, i, t_str_ini, t_str_fin), fontweight = 'bold', fontsize = 14)
        
        c = ['k', 'r', 'b', 'g', 'pink']
        label = ["Grid_0", "Nest_1", "Nest_2", "Nest_3", "Nest_4"]
        for j, xds in enumerate(zip(xds_out_case)):       
            ax1.plot(dates, hs[j,:], '-', color = c[j], markersize=2, label=label[j])
            ax1.set_xlim(dates[0], dates[-1])
            ax1.set_ylim(0, hs_max+0.5)
            plt.setp(ax1.get_xticklabels(), visible=False)
            ax1.set_title('Significant Wave Height [m]', fontweight = 'bold')
            ax1.legend()
    
            ax2.plot(dates, tm[j,:], '-', color = c[j], markersize=2, label=label[j])
            ax2.set_ylim(0, tm_max+0.5)
            plt.setp(ax2.get_xticklabels(), visible=False)
            ax2.set_title('Mean Period [s]', fontweight = 'bold')
    
            ax3.plot(dates, tp[j,:], '-', color = c[j], markersize=2, label=label[j])
            ax3.set_ylim(0, tp_max+0.5)
            plt.setp(ax3.get_xticklabels(), visible=False)
            ax3.set_title('Peak Period [s]', fontweight = 'bold')
    
            ax4.plot(dates, dr[j,:], '-', color = c[j], markersize=2, label=label[j])
            ax4.set_ylim([0, 360])
            ax4.xaxis.set_major_formatter(mdates.DateFormatter('%d/%b %Hh'))
            plt.setp(ax4.get_xticklabels(), visible=True)
            ax4.set_title('Wave direction [º]', fontweight = 'bold')
        
            fig.autofmt_xdate()

        # export fig
        p_ex = op.join(p_export_case, 'point_{0}.png'.format(i))
        fig.savefig(p_ex)

        # close fig 
        plt.close()

def plot_storm_track(name, pd_storm, case_id, lon0, lon1, lat0, lat1, num_nests, 
                     lon, lat, p_export, np_shore=np.array([])):
    '''
    Plots SWAN execution output table points time series

    name        - storm name
    pd_storm    - storm track pandas.DataFrame (x0, y0, R as metadata)
    case_id     - case number
    lon0, lon1  - longitude axes limits (lon0, lon1)
    lat0, lat1  - latitude axes limits (lat0, lat1)
    num_nests   - number of nested grids
    lon, lat    - nested grids boundary limits
    p_export    - path for exporting figure

    opt. args
    np_shore   - shoreline, np.array x = np_shore[:,0] y = np.shore[:,1]
    '''

    # make export dir
    if not op.isdir(p_export): os.makedirs(p_export)

    # get storm track data
    xt = pd_storm.lon
    yt = pd_storm.lat
    pmin = np.min(pd_storm.p0)#[0]
    vmean = np.mean(pd_storm.vf)#[0]
    gamma = pd_storm.move[0]

    # get storm metadata
    x0 = pd_storm.x0
    y0 = pd_storm.y0

    # plot and save figure
    fig = plt.figure(figsize=(12, 12))

    # plot shoreline
    if np_shore.any():
        x_shore = np_shore[:,0]
        y_shore = np_shore[:,1]
        plt.plot(x_shore, y_shore,'.', color='dimgray', markersize=3, label='')

    # plot track
    plt.plot(xt, yt, 'o-', linewidth=2, color='orangered', label='Great Circle')

    # plot target
    plt.plot(x0, y0, '+', mew=3, ms=15, color='dodgerblue', label='')
    
    # plot NESTED mesh
    for i in np.arange(0,num_nests):
        plt.plot([lon[i][0], lon[i][0]], [lat[i][0], lat[i][1]], '-k', linewidth=2, label='')
        plt.plot([lon[i][1], lon[i][1]], [lat[i][0], lat[i][1]], '-k', linewidth=2, label='')
        plt.plot([lon[i][0], lon[i][1]], [lat[i][0], lat[i][0]], '-k', linewidth=2, label='')
        plt.plot([lon[i][0], lon[i][1]], [lat[i][1], lat[i][1]], '-k', linewidth=2, label='')

    # plot parameters
    plt.axis('scaled')
    plt.xlim([lon0, lon1])
    plt.ylim([lat0, lat1])
    plt.xlabel('Longitude (º)')
    plt.ylabel('Latitude (º)')
    plt.title('{0}    Pmin: {1} hPa  /  Vmean: {2} km/h  /  Gamma: {3}º'.format(
            name, round(pmin,2),  round(vmean*1.872,2), round(gamma,2)), 
            fontsize = 14, fontweight='bold')
    plt.legend()

    # export fig
    p_save = op.join(p_export, 'track_coords_{0}.png'.format(case_id))
    fig.savefig(p_save)

    # close fig
    plt.close()
    
def plot_points_loc(p_export, xds_depth, np_shore, x_out, y_out, 
                    lon0, lon1, lat0, lat1, show=False):
    '''
    Plots SWAN execution output table points time series

    p_export    - path for exporting figure
    xds_depth   - elevation data 
    np_shore    - shoreline, np.array x = np_shore[:,0] y = np.shore[:,1]
    x_out, y_out- output points
    lon0, lon1  - longitude axes limits (lon0, lon1)
    lat0, lat1  - latitude axes limits (lat0, lat1)
    '''
    fig = plt.figure(figsize=(12, 12))
    
    xds_depth.elevation.plot(vmin=-500, vmax=0)
#    plt.plot(np_shore[:,0], np_shore[:,1],'.')
#    plt.plot(x_out,y_out,'.', markersize=12)
    
    for i, xi in enumerate(x_out):
        plt.plot(xi, y_out[i],'or', markersize = 22)
        plt.text(xi, y_out[i], u'{0}'.format(i), fontsize = 16, fontweight='bold', 
                 ha='center', va='center',color='w')

    # plot parameters
    plt.axis('scaled')
    plt.xlim(lon0, lon1)
    plt.ylim(lat0, lat1)
    plt.xlabel('Longitude (º)', fontsize = 12)
    plt.ylabel('Latitude (º)', fontsize = 12)
    plt.title('Output points location\n', fontsize = 16, fontweight = 'bold')
    
    # export fig
    p_save = op.join(p_export, 'points_location.png')
    fig.savefig(p_save)
    
    if show:    plt.show()

    # close fig
    plt.close()

def plot_grafiti_nonstat(name, xds_out, var_name, st_list, p_export,  
                        quiver=False, mda=False, np_shore=np.array([]), show=False):
    '''
    Plots non-stationary SWAN execution output for selected var, for every case

    name       - name of TC
    xds_out    - swan output (xarray.Dataset)
    var_name   - 'Hsig', 'Tm02', 'Tpsmoo'
    p_export   - path for exporting figures

    opt. args
    quiver     - True for adding directional quiver plot
    mda        - 
    np_shore   - shoreline, np.array x = np_shore[:,0] y = np.shore[:,1]
    '''

    # make export dir
    if not op.isdir(p_export): os.makedirs(p_export)

    var_max_all = xds_out[var_name].max(axis=0).max()
    var_min_all = xds_out[var_name].min(axis=0).min()

    for case_n, (case_ix) in enumerate(xds_out.case.values[:]):

        # select case
        xds_out_case = xds_out.sel(case=case_ix)
#        pos_lon = st_list[case_n]['lon'].values
#        pos_lat = st_list[case_n]['lat'].values
        pos_lon = st_list[case_ix]['lon'].values
        pos_lat = st_list[case_ix]['lat'].values

        # output case subfolder
        case_id = '{0:04d}'.format(int(case_ix))
        p_export_case = op.join(p_export, case_id)

        # plot variable times
        grafiti_nonstat(
            name, xds_out_case, var_name, pos_lon, pos_lat, p_export_case, show,
#            quiver=quiver, mda=mda, np_shore=np_shore, st_list=st_list[case_n], case_n=case_n,
            quiver=quiver, mda=mda, np_shore=np_shore, st_list=st_list[case_ix],
            case_n=case_ix, var_max_all=var_max_all, var_min_all=var_min_all)

def grafiti_nonstat(name, xds_out_case, var_name, pos_lon, pos_lat, p_export_case, show,
                    quiver=False, mda=False, np_shore=np.array([]), cmap='jet', 
                    st_list=None, case_n=None, var_max_all=None, var_min_all=None):
    '''
    Plots non-stationary SWAN execution output for Hsmax during the storm event

    name       - name of TC
    xds_out    - swan output (xarray.Dataset)
    var_name   - 'Hsig', 'Tm02', 'Tpsmoo'
    p_export   - path for exporting figures

    opt. args
    quiver     - True for adding directional quiver plot
    np_shore   - shoreline, np.array x = np_shore[:,0] y = np.shore[:,1]
    cmap       - matplotlib colormap
    '''

    if not op.isdir(p_export_case): os.makedirs(p_export_case)
    
    # maximum value of output variable over time --->  axis = (case, time, Y, X)
#    var_max = float(xds_out_case[var_name].max(axis=0).max())
#    var_min = float(xds_out_case[var_name].min(axis=0).min())
    var_max = var_max_all
    var_min = var_min_all

    # maximum value of output variable over time --->  axis = (case, time, Y, X)
    xds_out_max = xds_out_case.max(axis=0)

    # get mesh data from output dataset
    X = xds_out_case.X.values[:]
    Y = xds_out_case.Y.values[:]

    # get variable and units
    var = xds_out_max[var_name].values[:]
    var_units = xds_out_case[var_name].attrs['units']
    
    # get initial and final dates (in datetime.datetime format)
    date_0 = xds_out_case.time.values[0]
    date_1 = xds_out_case.time.values[-1]
    t_str_ini = pd.to_datetime(str(date_0)).strftime('%d-%b-%Y %H:%M%p')
    t_str_fin = pd.to_datetime(str(date_1)).strftime('%d-%b-%Y %H:%M%p')  
    
    # new figure
    fig, ax0 = plt.subplots(nrows=1, figsize=(12, 12))
    var_title = '{0}'.format(name)  # title

    # pcolormesh
    newcmap = custom_cmap(100, 'YlOrRd', 0.09, 0.9, 'YlGnBu_r', 0, 0.88)
    im = plt.pcolormesh(X, Y, var, cmap=newcmap, vmin= var_min, vmax= var_max)
    
    # add quiver plot 
    if quiver:
        var_title += '-Dir'
        var_dir = xds_out_max.Dir.values[:]
        aux_quiver(X, Y, var, var_dir)

    # shoreline
    if np_shore.any():
        x_shore = np_shore[:,0]
        y_shore = np_shore[:,1]
        plt.plot(x_shore, y_shore,'.', color='dimgray', markersize=3)

    # track
    plt.plot(pos_lon, pos_lat,'.-', color='black', markersize=12)

    # customize pcolormesh
    if mda:
        pmin = min(st_list['p0'].values)
        vmean = np.mean(st_list['vf'].values)
        gamma = np.mean(st_list['move'].values)
        plt.title('MDA_{0}     Pmin={1}mbar,  Vmean={2}km/h,  Gamma={3}º'.format(
                case_n, round(pmin,2), round(vmean*1.872,2), round(gamma,2)), 
                  fontsize = 14, fontweight='bold')
    else:
        plt.title('{0}   ({1})--({2})'.format(var_title, t_str_ini, t_str_fin),
                  fontsize = 14, fontweight='bold')
        
    plt.xlabel(xds_out_case.attrs['xlabel'], fontsize = 12)
    plt.ylabel(xds_out_case.attrs['ylabel'], fontsize = 12)
    plt.axis('scaled')
    plt.xlim(X[0], X[-1])
    plt.ylim(Y[0], Y[-1])

    # add custom colorbar
    divider = make_axes_locatable(ax0)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    plt.colorbar(im, cax=cax)
    plt.ylabel('{0}max ({1})'.format(var_name, var_units), fontsize = 12)

    # export fig
    p_ex = op.join(p_export_case, 'grafiti_{0}.png'.format(var_name))
    fig.savefig(p_ex)
    
    if show:    plt.show()

    # close fig 
    plt.close()
   
def plot_grafiti_discret_nonstat(name, xds_out, var_name, st_list, p_export, 
                                 step_grafiti, quiver=False, np_shore=np.array([]), show=False):
    '''
    Plots non-stationary SWAN execution output for selected var, for every case

    name         - name of TC
    xds_out      - swan output (xarray.Dataset)
    var_name     - 'Hsig', 'Tm02', 'Tpsmoo'
    p_export     - path for exporting figures
    step_grafiti - discretization step [º]

    opt. args
    quiver       - True for adding directional quiver plot
    np_shore     - shoreline, np.array x = np_shore[:,0] y = np.shore[:,1]
    '''

    # make export dir
    if not op.isdir(p_export): os.makedirs(p_export)

    for case_n, (case_ix) in enumerate(xds_out.case.values[:]):

        # select case
        xds_out_case = xds_out.sel(case=case_ix)
        pos_lon = st_list[case_n]['lon'].values
        pos_lat = st_list[case_n]['lat'].values

        # output case subfolder
        case_id = '{0:04d}'.format(case_ix)
        p_export_case = op.join(p_export, case_id)

        # plot variable times
        grafiti_nonstat_discret(
            name, xds_out_case, var_name, pos_lon, pos_lat, p_export_case, 
            step_grafiti, show, quiver=quiver, np_shore=np_shore)

def grafiti_nonstat_discret(name, xds_out_case, var_name, pos_lon, pos_lat, 
                            p_export_case, d_step, show,
                            quiver=False, np_shore=np.array([]), cmap='jet'):
    '''
    Plots non-stationary SWAN execution output for Hsmax during the storm event

    name       - name of TC
    xds_out    - swan output (xarray.Dataset)
    var_name   - 'Hsig', 'Tm02', 'Tpsmoo'
    p_export   - path for exporting figures
    d_step     - mesh discretization step (º)

    opt. args
    quiver     - True for adding directional quiver plot
    np_shore   - shoreline, np.array x = np_shore[:,0] y = np.shore[:,1]
    cmap       - matplotlib colormap
    '''

    if not op.isdir(p_export_case): os.makedirs(p_export_case)
    
    # maximum value of output variable over time --->  axis = (case, time, Y, X)
    var_max = float(xds_out_case[var_name].max(axis=0).max())
    var_min = float(xds_out_case[var_name].min(axis=0).min())

    # maximum value of output variable over time --->  axis = (case, time, Y, X)
    xds_out_max = xds_out_case.max(axis=0)

    # get mesh data from output dataset
    X = xds_out_case.X.values[:]
    Y = xds_out_case.Y.values[:]

    # get variable and units
    var = xds_out_max[var_name].values[:]
    var_units = xds_out_case[var_name].attrs['units']
    
    # discretization of meshgrid data (every "d_step"º)
    XX,YY = np.meshgrid(X, Y)
    dy,dx = np.mgrid[slice(Y[0],Y[-1],d_step), slice(X[0],X[-1],d_step)]
    var_discret = np.zeros(np.shape(dx))
    for j,xj in enumerate(dx[0][:-1]):
        for i,yi in enumerate(dy[:,0][:-1]):
            kx, ky = np.where((XX>=xj) & (XX<=dx[0][j+1]) & (YY>=yi) & (YY<=dy[:,0][i+1]))
            var_discret[i,j] = max(var[kx,ky])
    
    # get initial and final dates (in datetime.datetime format)
    date_0 = xds_out_case.time.values[0]
    date_1 = xds_out_case.time.values[-1]
    t_str_ini = pd.to_datetime(str(date_0)).strftime('%d-%b-%Y %H:%M%p')
    t_str_fin = pd.to_datetime(str(date_1)).strftime('%d-%b-%Y %H:%M%p')  
    
    # new figure
    fig, ax0 = plt.subplots(nrows=1, figsize=(12, 12))
    var_title = '{0}'.format(name)  

    # pcolormesh
    newcmap = custom_cmap(15, 'YlOrRd', 0.15, 0.9, 'YlGnBu_r', 0, 0.85)
    im = plt.pcolormesh(dx, dy, var_discret, cmap=newcmap, vmin= var_min, vmax= var_max)

    # add quiver plot 
    if quiver:
        var_title += '-Dir'
        var_dir = xds_out_max.Dir.values[:]
        aux_quiver(X, Y, var, var_dir)

    # shoreline
    if np_shore.any():
        x_shore = np_shore[:,0]
        y_shore = np_shore[:,1]
        plt.plot(x_shore, y_shore,'.', color='dimgray', markersize=3)
        
    # track
    plt.plot(pos_lon, pos_lat,'.-', color='black', markersize=12)

    # customize pcolormesh
    plt.title('{0}   ({1})--({2})'.format(var_title, t_str_ini, t_str_fin),
              fontsize = 14, fontweight='bold')
    plt.xlabel(xds_out_case.attrs['xlabel'], fontsize = 12)
    plt.ylabel(xds_out_case.attrs['ylabel'], fontsize = 12)

    plt.axis('scaled')
    plt.xlim(X[0], X[-1])
    plt.ylim(Y[0], Y[-1])

    # add custom colorbar
    divider = make_axes_locatable(ax0)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    plt.colorbar(im, cax=cax)
    plt.ylabel('{0}max ({1})'.format(var_name, var_units), fontsize = 12)

    # export fig
    p_ex = op.join(p_export_case, 'grafiti_{0}_discret.png'.format(var_name))
    fig.savefig(p_ex)
    
    if show:    plt.show()

    # close fig 
    plt.close()
   
def plot_output_panel(name, xds_panel, var_name, p_export, num_nests,  
                        lon, lat, quiver=False, np_shore=np.array([])):
    '''
    Plots non-stationary SWAN execution output for selected var, for every case

    name       - name of TC
    xds_panel  - swan output (xarray.Dataset) (list of all dataset outputs)
    var_name   - 'Hsig', 'Tm02', 'Tpsmoo'
    p_export   - path for exporting figures
    num_nests  - number of nested grids
    lon, lat   - nested grids boundary limits

    opt. args
    quiver     - True for adding directional quiver plot
    np_shore   - shoreline, np.array x = np_shore[:,0] y = np.shore[:,1]
    '''

    # make export dir
    if not op.isdir(p_export): os.makedirs(p_export)

    for case_ix in xds_panel[0].case.values[:]:

        # output case subfolder
        case_id = '{0:04d}'.format(case_ix)
        p_export_case = op.join(p_export, case_id)

        # plot variable times
        plot_panel_times(
            name, xds_panel, case_ix, var_name, p_export_case,
            num_nests, lon, lat, quiver=quiver, np_shore=np_shore)

def plot_panel_times(name, xds_panel, case_ix, var_name, p_export_case, num_nests, 
                     lon, lat, quiver=False, np_shore=np.array([]), cmap='jet'):
    '''
    Plots non-stationary SWAN output panel for selected case
    Designed for Majuro (4 meshgrids): 15km, 5km, 1km, 200m

    name           - name of TC
    case_ix        - case number of output dataset 
    xds_panel      - swan cases output (xarray.Dataset), list of cases/meshgrids
    p_export_case  - path for exporting figures
    '''

    if not op.isdir(p_export_case): os.makedirs(p_export_case)

    # select colorbar limits of xds_nonstat case_ix
    var_max = float(xds_panel[0][var_name].sel(case=case_ix).max(axis=0).max())
    var_min = float(xds_panel[0][var_name].sel(case=case_ix).min(axis=0).min())

    # iterate case_ix output over time
    for t in xds_panel[0].time.values[:]:
        xds_oct_0 = xds_panel[0].sel(case=case_ix)     # mesh 15km
        xds_oct_1 = xds_panel[1].sel(case=case_ix)     # mesh 5km
        xds_oct_2 = xds_panel[2].sel(case=case_ix)     # mesh 1km
        xds_oct_3 = xds_panel[3].sel(case=case_ix)     # mesh 200m
        
        # time string
        t_str = pd.to_datetime(str(t)).strftime('%d-%b-%Y %H:%M%p')

        # new figure
        xds_oct = [xds_oct_0, xds_oct_1, xds_oct_2, xds_oct_3]
        color = ['fuchsia','fuchsia','fuchsia','fuchsia']
        res = ['15 km', '5 km', '1 km', '250 m']

        fig, axes = plt.subplots(2, 2, figsize=(19, 12), constrained_layout=True)
        
        for c, (ax) in enumerate(axes.flat):
            im = plot_panel_var_t(name, var_max, var_min, xds_oct[c].sel(time=t), 
                                  var_name, num_nests, lon, lat, quiver=True, 
                                  np_shore=np_shore[c], cmap='jet', ax=ax)
            x_text = xds_oct[c].X.values[3]
            y_text = xds_oct[c].Y.values[-4]
            ax.text(x_text, y_text, res[c], color=color[c], va='top', ha='left', 
                    fontsize='26', fontweight='bold')

        cb = fig.colorbar(im, ax=axes.ravel().tolist(), shrink=1, aspect=35) #, label='Hs (m)') 
        cb.ax.set_title('Hs (m)', size=16)
        #cb.set_label(label='Hs (m)', size=14, rotation=0, labelpad=25, y=0.5)
#        cb.ax.set_ylabel('data',labelsize=14)
        cb.ax.tick_params(labelsize=14)
        cb.ax.yaxis.set_ticks_position('left')
        cb.ax.yaxis.set_label_position('left')
        
        fig.suptitle('MAJURO     {0}     {1}'.format(name, t_str), fontsize = 18, 
                     fontweight='bold') 

        # export fig
        p_ex = op.join(op.join(p_export_case, 'panel_{0}.png'.format(t_str)))
        fig.savefig(p_ex)

        # close fig 
        plt.close()
           
def plot_panel_var_t(name, var_max, var_min, xds_out_case, var_name, num_nests, lon, lat, 
                   quiver=False, np_shore=np.array([]), cmap='jet', ax=None):
    '''
    Plots non-stationary SWAN execution output for selected var and case

    name       - name of TC
    xds_out_0       - swan output of general meshgrid
    xds_out_case    - swan output (xarray.Dataset) of nested meshgrid
    var_name        - 'Hsig', 'Tm02', 'Tpsmoo'
    p_export        - path for exporting figures

    opt. args
    quiver     - True for adding directional quiver plot
    np_shore   - shoreline, np.array x = np_shore[:,0] y = np.shore[:,1]
    cmap       - matplotlib colormap
    '''
    ax = ax

    xds_oct = xds_out_case

    # get mesh data from output dataset
    X = xds_oct.X.values[:]
    Y = xds_oct.Y.values[:]
    flechas = 60
    delta_x = X[-1]-X[0]
    step = delta_x/flechas
    scale_quiver = 1/step * var_max

    # get variable and units
    var = xds_oct[var_name].values[:]
    
    # pcolormesh
    newcmap = custom_cmap(100, 'YlOrRd', 0.09, 0.9, 'YlGnBu_r', 0, 0.88)
    im = ax.pcolormesh(X, Y, var, cmap=newcmap, vmin= var_min, vmax= var_max)

    # add quiver plot 
    if quiver:
        var_dir = xds_oct.Dir.values[:]
        aux_quiver(X, Y, var, var_dir, scale_quiver, step, ax=ax)

    # shoreline
    if np_shore.any():
        x_shore = np_shore[:,0]
        y_shore = np_shore[:,1]
        ax.plot(x_shore, y_shore,'.', color='dimgray', markersize=3)

    # plot NESTED meshgrids
#    for i in np.arange(0,num_nests):
    for i in np.arange(0,num_nests-1):  # when not plotting nest nº 4 
        ax.plot([lon[i][0], lon[i][0]], [lat[i][0], lat[i][1]], '-k', linewidth=2, label='')
        ax.plot([lon[i][1], lon[i][1]], [lat[i][0], lat[i][1]], '-k', linewidth=2, label='')
        ax.plot([lon[i][0], lon[i][1]], [lat[i][0], lat[i][0]], '-k', linewidth=2, label='')
        ax.plot([lon[i][0], lon[i][1]], [lat[i][1], lat[i][1]], '-k', linewidth=2, label='')

    ax.axis('scaled')
    ax.set_xlim(X[0], X[-1])
    ax.set_ylim(Y[0], Y[-1])
    
    ax.tick_params(axis="x", labelsize=12)
    ax.tick_params(axis="y", labelsize=12)
    ax.set_xlabel('Longitude (º)', fontsize = 12)
    ax.set_ylabel('Latitude (º)', fontsize = 12)

    return im

def plot_output_panel_vortex(name, xds_panel, xds_vortex, var_name, p_export, num_nests,  
                        lon, lat, quiver=False, np_shore=np.array([])):
    '''
    Plots non-stationary SWAN execution output for selected var, for every case

    name       - name of TC
    xds_panel  - swan output (xarray.Dataset) (list of all dataset outputs)
    var_name   - 'Hsig', 'Tm02', 'Tpsmoo'
    p_export   - path for exporting figures
    num_nests  - number of nested grids
    lon, lat   - nested grids boundary limits

    opt. args
    quiver     - True for adding directional quiver plot
    np_shore   - shoreline, np.array x = np_shore[:,0] y = np.shore[:,1]
    '''

    # make export dir
    if not op.isdir(p_export): os.makedirs(p_export)

    for case_ix in xds_panel[0].case.values[:]:

        # output case subfolder
        case_id = '{0:04d}'.format(case_ix)
        p_export_case = op.join(p_export, case_id)

        # plot variable times
        plot_panel_vortex_times(name, xds_panel, xds_vortex, case_ix, var_name, 
                                p_export_case, num_nests, lon, lat, 
                                quiver=quiver, np_shore=np_shore)

def plot_panel_vortex_times(name, xds_panel, xds_vortex, case_ix, var_name, 
                            p_export_case, num_nests, lon, lat, 
                            quiver=False, np_shore=np.array([]), cmap='jet'):
    '''
    Plots non-stationary SWAN output panel for selected case
    Designed for general meshgrid: i.e. 15km

    name           - name of TC
    xds_panel      - swan cases output (xarray.Dataset), list of cases/meshgrids
    xds_vortex     - vortex cases output (xarray.Dataset)
    case_ix        - case number of output dataset 
    p_export_case  - path for exporting figures
    '''

    if not op.isdir(p_export_case): os.makedirs(p_export_case)

    # select colorbar limits of xds_nonstat case_ix
    var_max = float(xds_panel[0][var_name].sel(case=case_ix).max(axis=0).max())
    var_min = float(xds_panel[0][var_name].sel(case=case_ix).min(axis=0).min())
    
#    W_max = float(xds_vortex['W'].sel(case=xds_vortex.case.values[case_ix]).max(axis=0).max())
#    W_min = float(xds_vortex['W'].sel(case=xds_vortex.case.values[case_ix]).min(axis=0).max())
    W_max = float(xds_vortex['W'].max(axis=0).max())
    W_min = float(xds_vortex['W'].min(axis=0).max())
    
    cbar_max = [var_max, W_max]
    cbar_min = [var_min, W_min]
    var_name = [var_name, 'W']
    cbar_max = [W_max, var_max]
    cbar_min = [W_min, var_min]
    var_name = ['W', 'Hsig']

    # iterate case_ix output over time
    for t in xds_panel[0].time.values[:]:
        xds_oct_0 = xds_panel[0].sel(case=case_ix)
        xds_oct_v = xds_vortex
        
        # time string
        t_str = pd.to_datetime(str(t)).strftime('%d-%b-%Y_%H:%M%p')

        # new figure
        xds_oct = [xds_oct_v, xds_oct_0]
        color = ['fuchsia','fuchsia']
        res = ['vortex', '15 km']

        fig, axes = plt.subplots(1, 2, figsize=(20, 8), constrained_layout=True)
        for c, (ax) in enumerate(axes.flat):
            plot_panel_vortex_var_t(name, cbar_max[c], cbar_min[c], xds_oct[c].sel(time=t), 
                                  var_name[c], num_nests, lon, lat, quiver=True, 
                                  np_shore=np_shore, cmap='jet', ax=ax)
            x_text = xds_oct[c].X.values[3]
            y_text = xds_oct[c].Y.values[-4]
            ax.text(x_text, y_text, res[c], color=color[c], va='top', ha='left', 
                    fontsize='26', fontweight='bold')

        fig.suptitle('{0}     {1}'.format(name, t_str), fontsize=18, fontweight='bold') 

        # export fig
        p_ex = op.join(op.join(p_export_case, 'panel_{0}.png'.format(t_str)))
        fig.savefig(p_ex)

        # close fig 
        plt.close()
           
def plot_panel_vortex_var_t(name, var_max, var_min, xds_out_case, var_name, num_nests, lon, lat, 
                            quiver=False, np_shore=np.array([]), cmap='jet', ax=None):
    '''
    Plots non-stationary SWAN execution output for selected var and case

    name            - name of TC
    xds_out_0       - swan output of general meshgrid
    xds_out_case    - swan output (xarray.Dataset) of nested meshgrid
    var_name        - 'Hsig', 'Tm02', 'Tpsmoo'
    p_export        - path for exporting figures

    opt. args
    quiver     - True for adding directional quiver plot
    np_shore   - shoreline, np.array x = np_shore[:,0] y = np.shore[:,1]
    cmap       - matplotlib colormap
    '''
    ax = ax

    xds_oct = xds_out_case

    # get mesh data from output dataset
    X = xds_oct.X.values[:]
    Y = xds_oct.Y.values[:]
    flechas = 60
    delta_x = X[-1]-X[0]
    step = delta_x/flechas
    scale_quiver = 1/step * var_max

    # get variable and units
    var = xds_oct[var_name].values[:]
    var_units = xds_oct[var_name].attrs['units']

    # cmap
    if var_name=='Hsig':    
        newcmap = custom_cmap(100, 'YlOrRd', 0.09, 0.9, 'YlGnBu_r', 0, 0.88)
    elif var_name=='W':     
        newcmap = custom_cmap(100, 'plasma_r', 0.05, 0.9, 'viridis', 0.2, 1)
        
    # pcolormesh
    im = ax.pcolormesh(X, Y, var, cmap=newcmap, vmin= var_min, vmax= var_max)

    # add quiver plot 
    if quiver:
        var_dir = xds_oct.Dir.values[:]
        aux_quiver(X, Y, var, var_dir, scale_quiver, step, ax=ax)

    # shoreline
    if np_shore.any():
        x_shore = np_shore[:,0]
        y_shore = np_shore[:,1]
        ax.plot(x_shore, y_shore,'.', color='dimgray', markersize=3)

    # plot NESTED meshgrids
    for i in np.arange(0,num_nests):
#    for i in np.arange(0,num_nests-1):  # when not plotting nest nº 4 
        ax.plot([lon[i][0], lon[i][0]], [lat[i][0], lat[i][1]], '-k', linewidth=2, label='')
        ax.plot([lon[i][1], lon[i][1]], [lat[i][0], lat[i][1]], '-k', linewidth=2, label='')
        ax.plot([lon[i][0], lon[i][1]], [lat[i][0], lat[i][0]], '-k', linewidth=2, label='')
        ax.plot([lon[i][0], lon[i][1]], [lat[i][1], lat[i][1]], '-k', linewidth=2, label='')

    ax.axis('scaled')
    ax.set_xlim(X[0], X[-1])
    ax.set_ylim(Y[0], Y[-1])
    
    ax.tick_params(axis="x", labelsize=12)
    ax.tick_params(axis="y", labelsize=12)
    ax.set_xlabel('Longitude (º)', fontsize = 12)
    ax.set_ylabel('Latitude (º)', fontsize = 12)
    
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    plt.colorbar(im, cax=cax)
    plt.ylabel('{0} ({1})'.format(var_name, var_units), fontsize = 12)

    return im

def plot_panel_grafiti(name, xds_panel, var_name, pos_lon, pos_lat, p_export, num_nests,  
                        lon, lat, quiver=False, np_shore=np.array([])):
    '''
    Plots non-stationary SWAN execution output for selected var, for every case

    name       - name of TC
    xds_panel  - swan output (xarray.Dataset) (list of all dataset outputs)
    var_name   - 'Hsig', 'Tm02', 'Tpsmoo'
    p_export   - path for exporting figures
    num_nests  - number of nested grids
    lon, lat   - nested grids boundary limits

    opt. args
    quiver     - True for adding directional quiver plot
    np_shore   - shoreline, np.array x = np_shore[:,0] y = np.shore[:,1]
    '''

    # make export dir
    if not op.isdir(p_export): os.makedirs(p_export)

    for case_ix in xds_panel[0].case.values[:]:

        # output case subfolder
        case_id = '{0:04d}'.format(case_ix)
        p_export_case = op.join(p_export, case_id)

        # plot variable times
        grafiti_nonstat_panel(
            name, xds_panel, case_ix, var_name, pos_lon, pos_lat, p_export_case,
            num_nests, lon, lat, quiver=quiver, np_shore=np_shore)

def grafiti_nonstat_panel(name, xds_panel, case_ix, var_name, pos_lon, pos_lat, 
                          p_export_case, num_nests, lon, lat, 
                          quiver=False, np_shore=np.array([]), cmap='jet'):
    '''
    Plots non-stationary SWAN execution output for Hsmax during the storm event

    name       - name of TC
    xds_out    - swan output (xarray.Dataset)
    var_name   - 'Hsig', 'Tm02', 'Tpsmoo'
    p_export   - path for exporting figures

    opt. args
    quiver     - True for adding directional quiver plot
    np_shore   - shoreline, np.array x = np_shore[:,0] y = np.shore[:,1]
    cmap       - matplotlib colormap
    '''

    if not op.isdir(p_export_case): os.makedirs(p_export_case)
    
    # select colorbar limits of xds_nonstat case_ix
    var_max = float(xds_panel[0][var_name].sel(case=case_ix).max(axis=0).max())
    var_min = float(xds_panel[0][var_name].sel(case=case_ix).max(axis=0).min())

    # get initial and final dates (in datetime.datetime format)
    t = xds_panel[0].time.values[0]
    t1 = xds_panel[0].time.values[-1]
    t_str_ini = pd.to_datetime(str(t)).strftime('%d-%b-%Y %H:%M%p')
    t_str_fin = pd.to_datetime(str(t1)).strftime('%d-%b-%Y %H:%M%p')  
    t_str = pd.to_datetime(str(t)).strftime('%d-%b-%Y %H:%M%p')
    
    # maximum value of output variable over time --->  axis = (case, time, Y, X)
    xds_out_max_0 = xds_panel[0].sel(case=case_ix)    # mesh 15km
    xds_out_max_1 = xds_panel[1].sel(case=case_ix)    # mesh 5km
    xds_out_max_2 = xds_panel[2].sel(case=case_ix)    # mesh 1km
    xds_out_max_3 = xds_panel[3].sel(case=case_ix)    # mesh 200m

    # new figure
    xds_oct = [xds_out_max_0, xds_out_max_1, xds_out_max_2, xds_out_max_3]
    color = ['fuchsia','fuchsia','fuchsia','fuchsia']
    res = ['15km', '5km', '1km', '250m']

    fig, axes = plt.subplots(2, 2, figsize=(19, 12), constrained_layout=True)
    for c, (ax) in enumerate(axes.flat):
 #       im = plot_panel_var_t(name, var_max, var_min, xds_oct[c].max(axis=0), var_name, num_nests, lon, lat, 
        im = plot_var_t_panel(name, var_max, var_min, xds_oct[c], var_name, 
                              pos_lon, pos_lat, num_nests, lon, lat, 
                              quiver=False, np_shore=np_shore[c], cmap='jet', ax=ax)
        x_text = xds_oct[c].X.values[3]
        y_text = xds_oct[c].Y.values[-4]
        ax.text(x_text, y_text, res[c], color=color[c], va='top', ha='left', 
                fontsize='26', fontweight='bold')

    cb = fig.colorbar(im, ax=axes.ravel().tolist(), shrink=1, aspect=35) #, label='Hs (m)') 
    cb.ax.set_title('Hsmax (m)', size=16)
    #cb.set_label(label='Hs (m)', size=14, rotation=0, labelpad=25, y=0.5)
#        cb.ax.set_ylabel('data',labelsize=14)
    cb.ax.tick_params(labelsize=14)
    cb.ax.yaxis.set_ticks_position('left')
    cb.ax.yaxis.set_label_position('left')
    
    fig.suptitle('{0}     ({1})--({2})'.format(name, t_str_ini, t_str_fin), 
                 fontsize = 18, fontweight='bold') 

    # export fig
    p_ex = op.join(op.join(p_export_case, 'panel_{0}.png'.format(t_str)))
    fig.savefig(p_ex)

    # close fig 
    plt.close()


def plot_var_t_panel(name, var_max, var_min, xds_out_case, var_name, 
                     pos_lon, pos_lat, num_nests, lon, lat, 
                     quiver=False, np_shore=np.array([]), cmap='jet', ax=None):
    '''
    Plots non-stationary SWAN execution output for selected var and case

    name       - name of TC
    xds_out_0       - swan output of general meshgrid
    xds_out_case    - swan output (xarray.Dataset) of nested meshgrid
    var_name        - 'Hsig', 'Tm02', 'Tpsmoo'
    p_export        - path for exporting figures

    opt. args
    quiver     - True for adding directional quiver plot
    np_shore   - shoreline, np.array x = np_shore[:,0] y = np.shore[:,1]
    cmap       - matplotlib colormap
    '''
    ax = ax

#    xds_oct = xds_out_case.sel(time=t)
    xds_oct = xds_out_case.max(axis=0)

    # get mesh data from output dataset
    X = xds_out_case.X.values[:]
    Y = xds_out_case.Y.values[:]
    flechas = 60
    delta_x = X[-1]-X[0]
    step = delta_x/flechas
    scale_quiver = 1/step * var_max

    # get variable and units
    var = xds_oct[var_name].values[:]

    # cmap
    newcmap = custom_cmap(100, 'YlOrRd', 0.09, 0.9, 'YlGnBu_r', 0, 0.88)
    
    # pcolormesh
    im = ax.pcolor(X, Y, var, cmap=newcmap, vmin= var_min, vmax= var_max)

    # add quiver plot 
    if quiver:
#        var_title += '-Dir'
        var_dir = xds_oct.Dir.values[:]
        aux_quiver(X, Y, var, var_dir, scale_quiver, step, ax=ax)

    # shoreline
    if np_shore.any():
        x_shore = np_shore[:,0]
        y_shore = np_shore[:,1]
        ax.plot(x_shore, y_shore,'.', color='dimgray', markersize=3)

    # plot NESTED meshgrids
    for i in np.arange(0,num_nests-1):
        ax.plot([lon[i][0], lon[i][0]], [lat[i][0], lat[i][1]], '-k', linewidth=2, label='')
        ax.plot([lon[i][1], lon[i][1]], [lat[i][0], lat[i][1]], '-k', linewidth=2, label='')
        ax.plot([lon[i][0], lon[i][1]], [lat[i][0], lat[i][0]], '-k', linewidth=2, label='')
        ax.plot([lon[i][0], lon[i][1]], [lat[i][1], lat[i][1]], '-k', linewidth=2, label='')

    # track
    ax.plot(pos_lon, pos_lat,'.-', color='black', markersize=12)

    ax.axis('scaled')
    ax.set_xlim(X[0], X[-1])
    ax.set_ylim(Y[0], Y[-1])
    
    ax.tick_params(axis="x", labelsize=12)
    ax.tick_params(axis="y", labelsize=12)
    ax.set_xlabel('Longitude (º)', fontsize = 12)
    ax.set_ylabel('Latitude (º)', fontsize = 12)

    return im



def GetDivisors(x):
    l_div = []
    i = 1
    while i<x:
        if x%i == 0:
            l_div.append(i)
        i = i + 1
    return l_div

def GetBestRowsCols(n):
    'try to square number n, used at gridspec plots'

    sqrt_n = sqrt(n)
    if sqrt_n.is_integer():
        n_r = int(sqrt_n)
        n_c = int(sqrt_n)
    else:
        l_div = GetDivisors(n)
        n_c = l_div[int(len(l_div)/2)]
        n_r = int(n/n_c)

    return n_r, n_c


import matplotlib.gridspec as gridspec
def plot_grafiti_nonstat_matrix100(name, xds_out_ls, var_name, st_list, p_export, 
                                   fig_id, show=False, quiver=False, np_shore=np.array([])):
    '''
    Plots non-stationary SWAN execution output for selected var, for every case

    name       - name of TC
    xds_out    - swan output (xarray.Dataset)
    var_name   - 'Hsig', 'Tm02', 'Tpsmoo'
    p_export   - path for exporting figures

    opt. args
    quiver     - True for adding directional quiver plot
    np_shore   - shoreline, np.array x = np_shore[:,0] y = np.shore[:,1]
    '''

    # make export dir
    if not op.isdir(p_export): os.makedirs(p_export)

    var_max_all, var_min_all, n_clusters = [], [], []
    for i in range(len(xds_out_ls)):
        var_max_all.append(xds_out_ls[i][var_name].max(axis=0).max())
        var_min_all.append(xds_out_ls[i][var_name].min(axis=0).min())
        n_clusters.append(xds_out_ls[i].case.values.size)
    var_max_all = np.max(var_max_all)
    var_min_all = np.min(var_min_all)

    # get number of rows and cols for gridplot 
    n_clusters = np.sum(n_clusters)
    n_rows, n_cols = GetBestRowsCols(n_clusters)

    # plot figure
    fig = plt.figure(figsize=(23*1.5, 20*1.5))

    gs = gridspec.GridSpec(n_rows, n_cols, wspace=0, hspace=0)
    gr, gc = 0, 0


    for j in range(len(xds_out_ls)):
        xds_out = xds_out_ls[j]
        for case_n, (case_ix) in enumerate(xds_out.case.values[:]):
            # select case
            xds_out_case = xds_out.sel(case=case_ix)
            var_units = xds_out_case[var_name].attrs['units']
    
            pos_lon = st_list[case_ix]['lon'].values
            pos_lat = st_list[case_ix]['lat'].values
    
            # output case subfolder
            case_id = '{0:04d}'.format(int(case_ix))
            p_export_case = op.join(p_export, case_id)
    
            # plot variable times
            ax = plt.subplot(gs[gr, gc])
            pc = grafiti_matrix(
                ax, case_ix, name, xds_out_case, var_name, pos_lon, pos_lat, p_export_case, 
                quiver=quiver, np_shore=np_shore, st_list=st_list[case_ix], case_n=case_ix,
                var_max_all=var_max_all, var_min_all=var_min_all, 
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
    cb.set_label(label='{0} ({1})'.format(var_name, var_units), size=20, weight='bold')
    cb.ax.tick_params(labelsize=15)

    # show and return figure
    if show: plt.show()
    
    p_ex = op.join(p_export, 'grafiti_{0}_matrix_{1}.png'.format(var_name, fig_id))
    fig.savefig(p_ex)
    
    plt.close()

#def plot_grafiti_nonstat_matrix100(name, xds_out, var_name, st_list, p_export, fig_id, show=False, 
#                        quiver=False, np_shore=np.array([])):
#    '''
#    Plots non-stationary SWAN execution output for selected var, for every case
#
#    name       - name of TC
#    xds_out    - swan output (xarray.Dataset)
#    var_name   - 'Hsig', 'Tm02', 'Tpsmoo'
#    p_export   - path for exporting figures
#
#    opt. args
#    quiver     - True for adding directional quiver plot
#    np_shore   - shoreline, np.array x = np_shore[:,0] y = np.shore[:,1]
#    '''
#
#    # make export dir
#    if not op.isdir(p_export): os.makedirs(p_export)
#
#    var_max_all = xds_out[var_name].max(axis=0).max()
#    var_min_all = xds_out[var_name].min(axis=0).min()
#
#    # get number of rows and cols for gridplot 
#    n_clusters = xds_out.case.values.size
#    n_rows, n_cols = GetBestRowsCols(n_clusters)
#
#    # plot figure
#    fig = plt.figure(figsize=(23*1.5, 20*1.5))
#
#    gs = gridspec.GridSpec(n_rows, n_cols, wspace=0, hspace=0)
#    gr, gc = 0, 0
#
##    for ic in range(n_clusters):
#
#    for case_n, (case_ix) in enumerate(xds_out.case.values[:]):
#        # select case
#        xds_out_case = xds_out.sel(case=case_ix)
#        var_units = xds_out_case[var_name].attrs['units']
#
##        pos_lon = st_list[case_n]['lon'].values
##        pos_lat = st_list[case_n]['lat'].values
#        pos_lon = st_list[case_ix]['lon'].values
#        pos_lat = st_list[case_ix]['lat'].values
#
#        # output case subfolder
#        case_id = '{0:04d}'.format(int(case_ix))
#        p_export_case = op.join(p_export, case_id)
#
#        # plot variable times
#        ax = plt.subplot(gs[gr, gc])
#        pc = grafiti_matrix(
#            ax, case_ix, name, xds_out_case, var_name, pos_lon, pos_lat, p_export_case, 
##            quiver=quiver, np_shore=np_shore, st_list=st_list[case_n], case_n=case_n,
#            quiver=quiver, np_shore=np_shore, st_list=st_list[case_ix], case_n=case_ix,
#            var_max_all=var_max_all, var_min_all=var_min_all, 
#        )
#
#        # get lower positions
#        if gr==n_rows-1 and gc==0:
#            pax_l = ax.get_position()
#        elif gr==n_rows-1 and gc==n_cols-1:
#            pax_r = ax.get_position()
#
#        # counter
#        gc += 1
#        if gc >= n_cols:
#            gc = 0
#            gr += 1
#
##    # add a colorbar  
##    cb = fig.colorbar(gs[:,-1])
###    cb = fig.colorbar(im, ax=axes.ravel().tolist(), shrink=1, aspect=35) #, label='Hs (m)') 
##    cb.ax.set_title('Hs (m)', size=16)
##    #cb.set_label(label='Hs (m)', size=14, rotation=0, labelpad=25, y=0.5)
###        cb.ax.set_ylabel('data',labelsize=14)
##    cb.ax.tick_params(labelsize=14)
##    cb.ax.yaxis.set_ticks_position('left')
##    cb.ax.yaxis.set_label_position('left')
#            
##    cb = plt.colorbar(im, orientation='horizontal', pad=0.15)
##    cb.set_label(label='{0} ({1})'.format(var_name, var_units), size=15, weight='bold')
##    cb.ax.tick_params(labelsize='large')
#
#    cbar_ax = fig.add_axes([pax_l.x0, pax_l.y0-0.05, pax_r.x1 - pax_l.x0, 0.02])
#    cb = fig.colorbar(pc, cax=cbar_ax, orientation='horizontal')
#    cb.set_label(label='{0} ({1})'.format(var_name, var_units), size=20, weight='bold')
#    cb.ax.tick_params(labelsize=15)
#
##    ticklabs = cb.get_xticklabels()
##    cb.set_yticklabels(ticklabs, fontsize=10)
#    
#    # show and return figure
#    if show: plt.show()
#    
#    p_ex = op.join(p_export, 'grafiti_{0}_matrix_{1}.png'.format(var_name, fig_id))
#    fig.savefig(p_ex)
#    
#    plt.close()


def grafiti_matrix(ax, case_ix, name, xds_out_case, var_name, pos_lon, pos_lat, p_export_case,
                    quiver=False, np_shore=np.array([]), cmap='jet', 
                    st_list=None, case_n=None, var_max_all=None, var_min_all=None):
    '''
    Plots non-stationary SWAN execution output for Hsmax during the storm event

    name       - name of TC
    xds_out    - swan output (xarray.Dataset)
    var_name   - 'Hsig', 'Tm02', 'Tpsmoo'
    p_export   - path for exporting figures

    opt. args
    quiver     - True for adding directional quiver plot
    np_shore   - shoreline, np.array x = np_shore[:,0] y = np.shore[:,1]
    cmap       - matplotlib colormap
    '''

    var_max = var_max_all
    var_min = var_min_all

    # maximum value of output variable over time --->  axis = (case, time, Y, X)
    xds_out_max = xds_out_case.max(axis=0)

    # get mesh data from output dataset
    X = xds_out_case.X.values[:]
    Y = xds_out_case.Y.values[:]

    # get variable and units
    var = xds_out_max[var_name].values[:]
    
    # pcolormesh
    if var_name=='W': 	 newcmap = custom_cmap(100, 'plasma_r', 0.05, 0.9, 'viridis', 0.2, 1)
    if var_name=='Hsig': newcmap = custom_cmap(100, 'YlOrRd', 0.09, 0.9, 'YlGnBu_r', 0, 0.88)

    pc = ax.pcolormesh(X, Y, var, cmap=newcmap, vmin= var_min, vmax= var_max)
    
    # add quiver plot 
    if quiver:
        var_dir = xds_out_max.Dir.values[:]
        aux_quiver(X, Y, var, var_dir)

    # shoreline
    if np_shore.any():
        x_shore = np_shore[:,0]
        y_shore = np_shore[:,1]
        ax.plot(x_shore, y_shore,'.', color='dimgray', markersize=3)

    # track
    ax.plot(pos_lon, pos_lat,'.-', color='black', markersize=3)
    
    # number
    ax.text(X[5], Y[5], case_ix, color='fuchsia', fontweight='bold', fontsize=20)

    plt.axis('scaled')
    plt.xlim(X[0], X[-1])
    plt.ylim(Y[0], Y[-1])
    plt.axis('off')
    
    return pc

def plot_grafiti_matrix_hist(name, season, xds_out_ls, var_name, st_list, p_export, 
                                   fig_id, show=False, quiver=False, np_shore=np.array([])):
    '''
    Plots non-stationary SWAN execution output for selected var, for every case

    name       - name of TC
    xds_out    - swan output (xarray.Dataset)
    var_name   - 'Hsig', 'Tm02', 'Tpsmoo'
    p_export   - path for exporting figures

    opt. args
    quiver     - True for adding directional quiver plot
    np_shore   - shoreline, np.array x = np_shore[:,0] y = np.shore[:,1]
    '''

    # make export dir
    if not op.isdir(p_export): os.makedirs(p_export)

    var_max_all, var_min_all, n_clusters = [], [], []
    for i in range(len(xds_out_ls)):
        var_max_all.append(xds_out_ls[i][var_name].max(axis=0).max())
        var_min_all.append(xds_out_ls[i][var_name].min(axis=0).min())
        n_clusters.append(xds_out_ls[i].case.values.size)
    var_max_all = np.max(var_max_all)
    var_min_all = np.min(var_min_all)

    # get number of rows and cols for gridplot 
    n_clusters = np.sum(n_clusters)
    n_rows, n_cols = GetBestRowsCols(n_clusters)

    # plot figure
    fig = plt.figure(figsize=(23*1.5, 20*1.5))

    gs = gridspec.GridSpec(n_rows, n_cols, wspace=0, hspace=0)
    gr, gc = 0, 0


    for j in range(len(xds_out_ls)):
        xds_out = xds_out_ls[j]
        for case_n, (case_ix) in enumerate(xds_out.case.values[:]):
            # select case
            xds_out_case = xds_out.sel(case=case_ix)
            var_units = xds_out_case[var_name].attrs['units']
    
            pos_lon = st_list[case_ix]['lon'].values
            pos_lat = st_list[case_ix]['lat'].values
    
            # output case subfolder
            case_id = '{0:04d}'.format(int(case_ix))
            p_export_case = op.join(p_export, case_id)
    
            # plot variable times
            ax = plt.subplot(gs[gr, gc])
            pc = grafiti_matrix_hist(
                ax, case_ix, name[case_ix], season[case_ix], xds_out_case, var_name, 
                pos_lon, pos_lat, p_export_case, quiver=quiver, np_shore=np_shore, 
                var_max_all=var_max_all, var_min_all=var_min_all, 
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
    cb.set_label(label='{0} ({1})'.format(var_name, var_units), size=20, weight='bold')
    cb.ax.tick_params(labelsize=15)

    # show and return figure
    if show: plt.show()
    
    p_ex = op.join(p_export, 'grafiti_{0}_matrix_{1}.png'.format(var_name, fig_id))
    fig.savefig(p_ex)
    
    plt.close()

def grafiti_matrix_hist(ax, case_ix, name, season, xds_out_case, var_name, pos_lon, pos_lat, p_export_case,
                    quiver=False, np_shore=np.array([]), cmap='jet', 
                    var_max_all=None, var_min_all=None):
    '''
    Plots non-stationary SWAN execution output for Hsmax during the storm event

    name       - name of TC
    xds_out    - swan output (xarray.Dataset)
    var_name   - 'Hsig', 'Tm02', 'Tpsmoo'
    p_export   - path for exporting figures

    opt. args
    quiver     - True for adding directional quiver plot
    np_shore   - shoreline, np.array x = np_shore[:,0] y = np.shore[:,1]
    cmap       - matplotlib colormap
    '''

    var_max = var_max_all
    var_min = var_min_all

    # maximum value of output variable over time --->  axis = (case, time, Y, X)
    xds_out_max = xds_out_case.max(axis=0)

    # get mesh data from output dataset
    X = xds_out_case.X.values[:]
    Y = xds_out_case.Y.values[:]

    # get variable and units
    var = xds_out_max[var_name].values[:]
    
    # pcolormesh
    if var_name=='W': 	 newcmap = custom_cmap(100, 'plasma_r', 0.05, 0.9, 'viridis', 0.2, 1)
    if var_name=='Hsig': newcmap = custom_cmap(100, 'YlOrRd', 0.09, 0.9, 'YlGnBu_r', 0, 0.88)

    pc = ax.pcolormesh(X, Y, var, cmap=newcmap, vmin= var_min, vmax= var_max)
    
    # add quiver plot 
    if quiver:
        var_dir = xds_out_max.Dir.values[:]
        aux_quiver(X, Y, var, var_dir)

    # shoreline
    if np_shore.any():
        x_shore = np_shore[:,0]
        y_shore = np_shore[:,1]
        ax.plot(x_shore, y_shore,'.', color='dimgray', markersize=3)

    # track
    ax.plot(pos_lon, pos_lat,'.-', color='black', markersize=3)
    
    # number
    ax.text(X[5], Y[5], case_ix, color='fuchsia', fontweight='bold', fontsize=20)
    ax.text(X[20], Y[5], '{0} ({1})'.format(name, season), color='fuchsia', fontweight='bold', fontsize=16)

    plt.axis('scaled')
    plt.xlim(X[0], X[-1])
    plt.ylim(Y[0], Y[-1])
    plt.axis('off')
    
    return pc

