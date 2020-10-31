#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op

import numpy as np
import pandas as pd
import xarray as xr

import matplotlib.pyplot as plt

from .common import GetBestRowsCols, calc_quiver, custom_cmap

# import constants
from .config import _faspect, _fsize, _fdpi

from ..storms import get_category


# TODO:
#    PLOT TIME SERIES POINT OUTPUT
#    PLOT VIDEO PANEL WINDS INPUT - VAR OUTPUT
# 
#    REFACTOR plots/nonstationary (y mesh/io/wrap)
#    REPASO IMPORTs
#    DOC 

# aux.functions

def mesh2np(mesh):
    'generates np depth from swan mesh'

    # TODO MESH.GET_DG_XY()
    depth = mesh.depth
    xpc = mesh.dg['xpc']
    ypc = mesh.dg['ypc']
    xlenc = mesh.dg['xlenc']
    ylenc = mesh.dg['ylenc']
    mxc = depth.shape[1]
    myc = depth.shape[0]
    XX = np.linspace(xpc, xpc+xlenc, mxc)
    YY = np.linspace(ypc, ypc+ylenc, myc)

    return XX, YY, depth

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


# axes generation

def axplot_shore(ax, np_shore):
    'adds a shore (numpy) to axes'

    ax.plot(
        np_shore[:,0], np_shore[:,1], '.', color='dimgray',
        markersize=3, label=''
    )

def axplot_labels(ax, coords_mode):
    'sets xlabel and ylabel from swan coords_mode'

    xlab, ylab = 'X (m)', 'Y (m)'
    if coords_mode == 'SPHERICAL':
        xlab, ylab = 'Longitude (º)', 'Latitude (º)'
    ax.set_xlabel(xlab, fontweight='bold')
    ax.set_ylabel(ylab, fontweight='bold')

def axplot_quiver(ax, XX, YY, vv, vd):
    'Plot 2D quiver map'

    x_q, y_q, var_q, u, v = calc_quiver(XX, YY, vv, vd, size=12)
    ax.quiver(
        x_q, y_q, -u*var_q, -v*var_q,
        width=0.003,
        #scale = 0.5,
        scale_units='inches',
    )

def axplot_var_map(ax, XX, YY, vv,
                   vmin = None, vmax = None,
                   cmap = plt.get_cmap('seismic'),
                  ):
    'plot 2D map with variable data'

    # cplot v lims
    if vmin == None: vmin = np.nanmin(vv)
    if vmax == None: vmax = np.nanmax(vv)

    # plot variable 2D map
    pm = ax.pcolormesh(
        XX, YY, vv,
        cmap=cmap,
        vmin=vmin, vmax=vmax,
        shading='auto',
    )

    # fix axes
    ax.set_xlim(XX[0], XX[-1])
    ax.set_ylim(YY[0], YY[-1])

    # return pcolormesh
    return pm

def axplot_storm_track(ax, st, cat_colors=True):
    'plots storm track, category and target point'

    # storm trac parameters
    xt = st.lon
    yt = st.lat

    # plot track
    plt.plot(
        st.lon, st.lat, '-', linewidth=2,
        color='darkgrey', label='Storm Track'
    )

    # plot categories
    if cat_colors:

        # get category
        categ = np.array(get_category(st.p0))

        for c in range(7):
            lonc = st.lon[np.where(categ==c)[0]]
            latc = st.lat[np.where(categ==c)[0]]
            label = 'cat {0}'.format(c)
            if c==6: label = 'unknown'
            ax.plot(
                lonc, latc, '.', color=get_storm_color(c),
                markersize=10, label=label,
            )

    # target point
    ax.plot(
        st.x0, st.y0, '+', mew=3, ms=15,
        color='dodgerblue', label='target',
    )


# Plot functions list

def plot_project_site(swan_proj):
    '''
    Plots SwanProject site
        - bathymetry
        - control points
        - shoreline
        - nested meshes locations
    '''

    # figure
    fig, (axs) = plt.subplots(
        nrows=1, ncols=1,
        figsize=(_fsize*_faspect, _fsize),
    )

    # plot bathymetry
    mesh = swan_proj.mesh_main
    XX, YY, depth = mesh2np(mesh)

    pm = axplot_var_map(
        axs, XX, YY, depth,
        cmap = 'gist_earth_r',
    )
    cbar = fig.colorbar(pm, ax=axs)
    cbar.ax.set_ylabel('depth (m)', rotation=90, va="bottom", fontweight='bold')

    # mesh coordinates labels
    axplot_labels(axs, swan_proj.params['coords_mode'])

    # plot shoreline
    shore = swan_proj.shore
    if shore.any():
        axplot_shore(axs, np_shore=shore)

    # plot output points (control points)
    x_out = swan_proj.params['output_points_x']
    y_out = swan_proj.params['output_points_y']
    if x_out and len(x_out)==len(y_out):
        axs.plot(
            x_out, y_out, '.', color='red',
            markersize=16, label='control points'
        )

    # plot nested meshes
    for mn in swan_proj.mesh_nested_list:
        XX_m, YY_m, _ = mesh2np(mn)
        xr = [XX_m[0], XX_m[-1], XX_m[-1], XX_m[0], XX_m[0]]
        yr = [YY_m[0], YY_m[0], YY_m[-1], YY_m[-1], YY_m[0]]

        axs.plot(xr, yr, '-', color='w', linewidth='4', label=mn.ID)

    # title
    axs.set_title('SWAN Project Site: {0}'.format(swan_proj.name),
                  fontsize=16, fontweight='bold')

    # turn on legend
    plt.legend(loc='upper right', prop={'size':10})

    return fig

def plot_case_input(swan_proj, storm_track_list=[], case_number=0):
    '''
    Plots case input
        - boundary waves time series
        - TCs storm track
        - shoreline
        - control point
    '''


    # mesh values
    mesh = swan_proj.mesh_main
    XX, YY, _ = mesh2np(mesh)

    # figure
    fig, (axs) = plt.subplots(
        nrows=1, ncols=1,
        figsize=(_fsize*_faspect, _fsize),
    )

    # mesh values
    mesh = swan_proj.mesh_main
    XX, YY, _ = mesh2np(mesh)

    # mesh coordinates labels
    axplot_labels(axs, swan_proj.params['coords_mode'])

    # plot shoreline
    shore = swan_proj.shore
    if shore.any():
        axplot_shore(axs, np_shore=shore)

    # TODO: add plot wave boundaries series

    # plot storm track
    if storm_track_list:
        st = storm_track_list[case_number]  # select storm track for this case
        axplot_storm_track(axs, st)

        # add text to title
        ttl_st = '\nPmin: {0:.2f} hPa / Vmean: {1:.2f} km/h / Gamma: {2:.2f}º'.format(
            np.min(st.p0), np.mean(st.vf)*1.872, st.move[0])

    # title
    axs.set_title(
        'SWAN Project: {0}, Case: {1:04d} {2}'.format(
            swan_proj.name, case_number, ttl_st),
        fontsize=16, fontweight='bold')

    # fix axes
    axs.set_xlim(XX[0], XX[-1])
    axs.set_ylim(YY[0], YY[-1])

    # legend
    plt.legend(loc='upper right', prop={'size':10})

    axs.set_facecolor('lightcyan')

    return fig

def plot_case_vortex_grafiti(swan_wrap, storm_track_list=[], case_number=0,
                            mesh=None):
    '''
    # TODO: documentar
    '''
    # TODO: refactor con el siguiente grafiti

    # swan project
    swan_proj = swan_wrap.proj

    # default to main mesh
    if mesh == None: mesh = swan_proj.mesh_main

    # read vortex Wind module and direction
    case_id = '{0:04d}'.format(case_number)
    p_case = op.join(swan_proj.p_cases, case_id)
    code = 'wind_{0}'.format(mesh.ID)
    p_vortex = op.join(p_case, 'vortex_{0}.nc'.format(code))

    xds_vortex = xr.open_dataset(p_vortex)
    var_name = 'W'


    # figure
    fig, (axs) = plt.subplots(
        nrows=1, ncols=1,
        figsize=(_fsize*_faspect, _fsize),
    )

    # plot grafiti
    xds_var = xds_vortex[var_name]
    var_units = xds_var.units

    # get mesh data from output dataset
    coords_mode = swan_proj.params['coords_mode']
    if coords_mode == 'SPHERICAL':
        xa, ya = 'lon', 'lat'
    else:
        xa, ya = 'X', 'Y'
    X = xds_var[xa].values[:]
    Y = xds_var[ya].values[:]

    # maximum and minimum values 
    vmax = float(xds_var.max().values)
    vmin = float(xds_var.min().values)

    # grafiti
    xds_var_max = xds_var.max(dim='time')
    var_max = xds_var_max.values[:]

    ccmap = custom_cmap(100, 'plasma_r', 0.05, 0.9, 'viridis', 0.2, 1)
    pm = axplot_var_map(
        axs, X, Y, var_max,  # TODO output.T, aqui no
        vmin = vmin, vmax = vmax,
        cmap = ccmap,
    )
    cbar = fig.colorbar(pm, ax=axs)
    cbar.ax.set_ylabel(
        '{0} max. ({1})'.format(var_name, var_units),
        rotation=90, va="bottom", fontweight='bold',
        labelpad=15,
    )

    # plot shoreline
    shore = swan_proj.shore
    if shore.any():
        axplot_shore(axs, np_shore=shore)

    # plot storm track
    if storm_track_list:
        st = storm_track_list[case_number]  # select storm track for this case
        axplot_storm_track(axs, st, cat_colors=False)

    # mesh coordinates labels
    axplot_labels(axs, swan_proj.params['coords_mode'])

    # title
    axs.set_title(
        'SWAN Project: {0}, Case: {1:04d} - Vortex Model'.format(
            swan_proj.name, case_number ),
        fontsize=16, fontweight='bold')

    return fig


def plot_case_output_grafiti(
    swan_wrap, var_name='Hsig', case=0, mesh=None,
    storm_track_list = []):
    '''
    # TODO DOC
    Plots case output "grafiti": variable max
        - todo1
        - todo2
    '''

    # swan project
    swan_proj = swan_wrap.proj

    # default to main mesh
    if mesh == None: mesh = swan_proj.mesh_main

    # load output
    xds_out = swan_wrap.extract_output(
        case_ini=case, case_end=case+1,
        mesh=mesh,
    ).squeeze(drop=True)


    # figure
    fig, (axs) = plt.subplots(
        nrows=1, ncols=1,
        figsize=(_fsize*_faspect, _fsize),
    )

    # plot grafiti
    xds_var = xds_out[var_name]
    var_units = xds_var.units

    # get mesh data from output dataset
    coords_mode = swan_proj.params['coords_mode']
    if coords_mode == 'SPHERICAL':
        xa, ya = 'lon', 'lat'
    else:
        xa, ya = 'X', 'Y'
    X = xds_var[xa].values[:]
    Y = xds_var[ya].values[:]

    # maximum and minimum values 
    vmax = float(xds_var.max().values)
    vmin = float(xds_var.min().values)

    # grafiti
    xds_var_max = xds_var.max(dim='time')
    var_max = xds_var_max.values[:]

    ccmap = custom_cmap(15, 'YlOrRd', 0.15, 0.9, 'YlGnBu_r', 0, 0.85)
    pm = axplot_var_map(
        axs, X, Y, var_max.T,
        vmin = vmin, vmax = vmax,
        cmap = ccmap,
    )
    cbar = fig.colorbar(pm, ax=axs)
    cbar.ax.set_ylabel(
        '{0} max. ({1})'.format(var_name, var_units),
        rotation=90, va="bottom", fontweight='bold',
        labelpad=15,
    )

    # plot shoreline
    shore = swan_proj.shore
    if shore.any():
        axplot_shore(axs, np_shore=shore)

    # plot storm track
    if storm_track_list:
        st = storm_track_list[case]  # select storm track for this case
        axplot_storm_track(axs, st, cat_colors=False)

    # mesh coordinates labels
    axplot_labels(axs, swan_proj.params['coords_mode'])

    # title
    date_0 = xds_out.time.values[0]
    date_1 = xds_out.time.values[-1]
    fmt = '%d-%b-%Y %H:%M%p'
    t_str_ini = pd.to_datetime(str(date_0)).strftime(fmt)
    t_str_fin = pd.to_datetime(str(date_1)).strftime(fmt)
    ttl_t = '\n{0} : {1}'.format(t_str_ini, t_str_fin)
    axs.set_title(
        'SWAN Project: {0}, Case: {1:04d} {2}'.format(
            swan_proj.name, case, ttl_t),
        fontsize=16, fontweight='bold')

    return fig



