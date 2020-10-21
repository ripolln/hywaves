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
we = waves['2000-01-02 00:00':'2000-01-02 02:00'][vs]
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
    'xlenc': 700,   # grid length in x
    'ylenc': 1000,  # grid length in y
    'mxc': 1,       # number mesh x
    'myc': 1,       # number mesh y
    'dxinp': 700,   # size mesh x
    'dyinp': 1000,  # size mesh y
}

# depth value
main_mesh.depth = np.ones((2,2)) * 155

# computational grid description
main_mesh.cg = {
    'xpc': 0,
    'ypc': 0,
    'alpc': 0,
    'xlenc': 700,
    'ylenc': 1000,
    'mxc': 100,
    'myc': 50,
    'dxinp': 7,
    'dyinp': 20,
}

sp.set_main_mesh(main_mesh)


# SWAN parameters (sea level, jonswap gamma)
sp.params = {
    'sea_level': 4,
    'jonswap_gamma': 1.9,
    'cdcap': None,
    'coords_spherical': None,
    'waves_period': 'MEAN',
    'maxerr': None,
}


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

