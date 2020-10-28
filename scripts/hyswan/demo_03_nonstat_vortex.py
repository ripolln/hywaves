#!/usr/bin/env python
# -*- coding: utf-8 -*-

# TODO: add nested mesh to demo
# TODO: add plots
# TODO: check input and output grid orientation (caution transpose)

# common 
import sys
import os
import os.path as op

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

# dev
sys.path.insert(0, op.join(op.dirname(__file__), '..', '..'))

# swan wrap module 
from hywaves.swan.storms import track_from_parameters
from hywaves.swan.wrap import SwanProject, SwanMesh, SwanWrap_NONSTAT


# TODO: actualizar a track_site_parameters

# --------------------------------------
# data
p_data = op.abspath(op.join(op.dirname(__file__), '..', '..', 'data'))
p_demo = op.join(p_data, 'demo')

# bathymetry (from .nc file)
p_bathy = op.join(p_demo, 'depth.nc')

xds_bathy = xr.open_dataset(p_bathy)
lon = xds_bathy.lon.values[:]
lat = xds_bathy.lat.values[:]

depth = xds_bathy.elevation.values[:] * -1  # elevation to depth 


# --------------------------------------
# SWAN case input: waves_event (empty) + storm_track (from MDA parameters)

# time array for SWAN case input
date_ini = '2000-01-02 00:00'
hours = 6
time = pd.date_range(date_ini, periods=hours, freq='H')


# generate empty wave event
we = pd.DataFrame(index=time, columns=['hs', 'per', 'dir', 'spr', 'U10', 'V10'])
we['level'] = 0
we['tide'] = 0


# generate storm track from MDA parameters 
pmin = 924.9709      # center pressure 
vmean = 69.0352      # translation velocity (km/h)
delta = 87.8432      # azimut
gamma = 92.7126      # translation angle (nautical convention)
x0, y0 = 167.5, 9.5  # Roi-Namur coordinates
x1 = 175             # enter point in the computational grid
R = 4                # smaller radius in degrees

st = track_from_parameters(
    pmin, vmean, delta, gamma, x0, y0, x1, R, date_ini, hours,
    great_circle=True
)

# TODO vf a vfx, vfy

print('\ninput storm track')
print(st)


# --------------------------------------
# SWAN project (config bathymetry, parameters, computational grid)

p_proj = op.join(p_data, 'projects')  # swan projects main directory
n_proj = '03_vortex'                  # project name

sp = SwanProject(p_proj, n_proj)


# --------------------------------------
# SWAN main mesh
main_mesh = SwanMesh()

# depth grid description (input bathymetry grid)
main_mesh.dg = {
    'xpc': lon[0],                             # x origin
    'ypc': lat[0],                             # y origin
    'alpc': 0,                                 # x-axis direction 
    'xlenc': lon[-1]-lon[0],                   # grid length in x
    'ylenc': lat[-1]-lat[0],                   # grid length in y
    'mxc': depth.shape[1]-1,                   # number mesh x (TODO -1?)
    'myc': depth.shape[0]-1,                   # number mesh y (TODO -1?)
    'dxinp': (lon[-1]-lon[0])/depth.shape[1],  # size mesh x
    'dyinp': (lat[-1]-lat[0])/depth.shape[0],  # size mesh y
}

# depth value (from file)
main_mesh.depth = depth

# computational grid description
main_mesh.cg = {
    'xpc': 160,
    'ypc': 2,
    'alpc': 0,
    'xlenc': 15,
    'ylenc': 13,
    'mxc': 15,
    'myc': 13,
    'dxinp': 1,
    'dyinp': 1,
}

sp.set_main_mesh(main_mesh)


# SWAN parameters (sea level, jonswap gamma)
input_params = {
    'set_level': 0,
    'set_convention': 'NAUTICAL',
    'set_cdcap': 2.5*10**-3,

    'coords_mode': 'SPHERICAL',
    'coords_projection': 'GCM',

    'boundw_jonswap': 3.3,
    'boundw_period': 'MEAN',

    'physics':[
        'WIND DRAG WU',
        'GEN3 ST6 5.7E-7 8.0E-6 4.0 4.0 UP HWANG VECTAU TRUE10',
        'QUAD iquad=8',
        'WCAP',
        #'SETUP',  # not compatible with spherical coords
        'TRIADS',
        'DIFFRAC',
    ],

    'numerics':[
        'PROP BSBT',
    ]
}
sp.set_params(input_params)


# SWAN output points
sp.x_out = [167.5, 167.5, 167]
sp.y_out = [9.5, 9.6, 9.45]


# --------------------------------------
# SWAN wrap NONSTAT (create case files, launch SWAN num. model, extract output)

sw = SwanWrap_NONSTAT(sp)

# build non-stationary cases from wave_events list and storm_tracks list
sw.build_cases([we], storm_track_list=[st], make_waves=False)

# run SWAN
sw.run_cases()

# extract grid output from non-stationary cases
xds_out_main = sw.extract_output()
print('\noutput main mesh')
print(xds_out_main)

