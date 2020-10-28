#!/usr/bin/env python
# -*- coding: utf-8 -*-

# TODO: add nested mesh to demo
# TODO: add plots
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
from hywaves.swan.storms import track_site_parameters
from hywaves.swan.wrap import SwanProject, SwanMesh, SwanWrap_NONSTAT


# --------------------------------------
# data
p_data = op.abspath(op.join(op.dirname(__file__), '..', '..', 'data'))
p_demo = op.join(p_data, 'demo')

# bathymetry (from .nc file)
p_bathy = op.join(p_demo, 'depth_samoa.nc')

xds_bathy = xr.open_dataset(p_bathy)

# sign convention [0º,360º]
xds_lon = xds_bathy.lon.values
xds_lon[xds_lon<0] = xds_lon[xds_lon<0] + 360
xds_bathy.lon.values[:] = xds_lon
xds_bathy = xds_bathy.sel(lon=slice(182,194), lat=slice(-20,-7.5))

lon = xds_bathy.lon.values[:]
lat = xds_bathy.lat.values[:]

depth = xds_bathy.elevation.values[:] * -1  # elevation to depth 

# shoreline (from .nc file)
p_shore = op.join(p_demo, 'shore_samoa.npy')

np_shore = np.load(p_shore)

# sign convention [0º,360º]
lon_shore = np_shore[:,0]
lon_shore[lon_shore<0] = lon_shore[lon_shore<0] + 360
np_shore[:,0] = lon_shore


# --------------------------------------
# target location
target = 'Samoa'
x0, y0 = -172.8, -13.5      # coordinates
if x0 < 0:  x0 = x0 + 360   # sign convention [0º,360º]


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
#x1 = 175             # enter point in the computational grid
R = 4                # smaller radius in degrees

tstep = 30           # computational time step (minutes) for track interpolation

# TODO: genera storm track de un unico instante de tiempo
st = track_site_parameters(
    tstep, pmin, vmean, delta, gamma, x0, y0, lon[0], lon[-1], lat[0], lat[-1],
    R, date_ini
)

print('\ninput storm track')
print(st)


# --------------------------------------
# SWAN project (config bathymetry, parameters, computational grid)

p_proj = op.join(p_data, 'projects')  # swan projects main directory
n_proj = '04_vortex_samoa'            # project name

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
    'mxc': depth.shape[1]-1,                   # number mesh x
    'myc': depth.shape[0]-1,                   # number mesh y
    'dxinp': (lon[-1]-lon[0])/depth.shape[1],  # size mesh x
    'dyinp': (lat[-1]-lat[0])/depth.shape[0],  # size mesh y
}

# depth value (from file)
main_mesh.depth = depth

# computational grid description
main_mesh.cg = {
    'xpc': 182,
    'ypc': -20,
    'alpc': 0,
    'xlenc': 12,
    'ylenc': 12.5,
    'mxc': int(round(12/0.136)),    # grid resolution of 15km (=0.136º)
    'myc': int(round(12.5/0.136)),
    'dxinp': 12/int(round(12/0.136)),
    'dyinp': 12.5/int(round(12.5/0.136)),
}

sp.set_main_mesh(main_mesh)

# SWAN parameters (sea level, jonswap gamma, ...)
input_params = {
    'set_level': 0,
    'set_convention': 'NAUTICAL',
    'set_cdcap': 2.5*10**-3,

    'boundw_jonswap': 3.3,
    'boundw_period': 'MEAN',

    'boundn_mode': 'CLOSED',

    'wind_deltinp': '30 MIN',
    'level_deltinp': '1 HR',

    'compute_deltc': '5 MIN',
    'output_deltt': '30 MIN',

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

