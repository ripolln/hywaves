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

# dev
sys.path.insert(0, op.join(op.dirname(__file__), '..', '..'))

# swan wrap module
from hywaves.swan.wrap import SwanProject, SwanMesh, SwanWrap_NONSTAT


# --------------------------------------
# data
p_data = op.abspath(op.join(op.dirname(__file__), '..', '..', 'data'))
p_demo = op.join(p_data, 'demo')

# test data: csiro point 
p_waves_demo = op.join(p_demo, 'waves_csiro_demo.nc')


# --------------------------------------
# SWAN case input: waves_event 

# non-stationary case requires a wave_event (dim: time). variables:
# time
# waves: hs, t02, dir, spr
# wind: U10, V10 (units_ m/s)
# water level: level, tide     (not in csiro?)

xds_waves = xr.open_dataset(p_waves_demo)
xds_waves = xds_waves.squeeze()   # remove lon,lat dim (len=1)
waves = xds_waves.to_dataframe()  # xarray --> pandas

# now we generate the wave event 
vs = ['hs', 't02', 'dir', 'spr', 'U10', 'V10']
we = waves['2000-01-02 00:00':'2000-01-02 03:00'][vs]
we['level'] = 0 # no water level data 
we['tide'] = 0 # no tide data 
we.rename(columns={'t02': 'per'}, inplace=True)  # rename for swan

print('\ninput wave event')
print(we)


# --------------------------------------
# SWAN project (config bathymetry, parameters, computational grid)

p_proj = op.join(p_data, 'projects')  # swan projects main directory
n_proj = '02_nonstat'                 # project name

sp = SwanProject(p_proj, n_proj)


# --------------------------------------
# SWAN main mesh
main_mesh = SwanMesh()

# depth grid description (input bathymetry grid)
main_mesh.dg = {
    'xpc': 0,       # x origin
    'ypc': 0,       # y origin
    'alpc': 0,      # x-axis direction 
    'xlenc': 400,   # grid length in x
    'ylenc': 400,  # grid length in y
    'mxc': 1,       # number mesh x
    'myc': 1,       # number mesh y
    'dxinp': 400,   # size mesh x
    'dyinp': 400,  # size mesh y
}

# depth value
main_mesh.depth = np.ones((2,2)) * 155

# computational grid description
main_mesh.cg = {
    'xpc': 0,
    'ypc': 0,
    'alpc': 0,
    'xlenc': 400,
    'ylenc': 400,
    'mxc': 40,
    'myc': 20,
    'dxinp': 10,
    'dyinp': 20,
}

sp.set_main_mesh(main_mesh)

# --------------------------------------
# SWAN nest1 mesh
mesh_nest1 = SwanMesh()

# depth grid description
mesh_nest1.dg = {
    'xpc': 50,
    'ypc': 100,
    'alpc': 0,
    'xlenc': 80,
    'ylenc': 100,
    'mxc': 8,
    'myc': 10,
    'dxinp': 10,
    'dyinp': 10,
}

# depth value
mesh_nest1.depth = np.ones((10,8)) * 158

# computational grid description
mesh_nest1.cg = {
    'xpc': 50,
    'ypc': 100,
    'alpc': 0,
    'xlenc': 80,
    'ylenc': 100,
    'mxc': 8,
    'myc': 10,
    'dxinp': 10,
    'dyinp': 10,
}

sp.set_nested_mesh_list([mesh_nest1])


# --------------------------------------
# SWAN parameters (sea level, jonswap gamma)
input_params = {
    'set_level': 4,
    'set_convention': 'NAUTICAL',

    'boundw_jonswap': 1.9,
    'boundw_period': 'MEAN',

    'boundn_mode': 'CLOSED',

    'wind_deltinp': '1 HR',
    'level_deltinp': '1 HR',

    'compute_deltc': '5 MIN',
    'output_deltt': '10 MIN',

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


# --------------------------------------
# SWAN wrap NONSTAT (create case files, launch SWAN num. model, extract output)

sw = SwanWrap_NONSTAT(sp)

# build non-stationary cases from wave_events list
sw.build_cases([we])  # test one event

# run SWAN
sw.run_cases()

# extract output from non-stationary cases
xds_out_main = sw.extract_output()
print('\noutput main mesh')
print(xds_out_main)

# extract output from nest1 mesh 
xds_out_nest1 = sw.extract_output(mesh=sp.mesh_nested_list[0])
print('\noutput nest1 mesh')
print(xds_out_nest1)
