#!/usr/bin/env python
# -*- coding: utf-8 -*-

# common 
import sys
import os
import os.path as op

import numpy as np
import pandas as pd
import xarray as xr

# dev library 
sys.path.insert(0, op.join(op.dirname(__file__), '..', '..'))

# swan wrap module
from hywaves.swan.wrap import SwanProject, SwanWrap_STAT


# --------------------------------------
# data
p_data = op.abspath(op.join(op.dirname(__file__), '..', '..', 'data'))
p_demo = op.join(p_data, 'demo')

# test data: csiro waves point 
p_waves_demo = op.join(p_demo, 'waves_csiro_demo.nc')


# --------------------------------------
# SWAN case input: waves_event 

# waves: hs, per, dir, spr     (per can be PEAK or MEAN wave period)
# water level: level, tide     (not in csiro?)

xds_waves = xr.open_dataset(p_waves_demo)
xds_waves = xds_waves.squeeze()   # remove lon,lat dim (len=1)
waves = xds_waves.to_dataframe()  # xarray --> pandas

# now we generate the wave event 
vs = ['hs', 't02', 'dir', 'spr']
we = waves['2000-01-02 00:00':'2000-01-03 00:00'][vs]
we['level'] = 0 # no water level data 
we['tide'] = 0 # no tide data 
we.rename(columns={'t02': 'per'}, inplace=True)  # rename for swan

print('\ninput wave events')
print(we)


# --------------------------------------
# SWAN project 

p_proj = op.join(p_data, 'projects')  # swan projects main directory
n_proj = '01_stat'                    # project name

sp = SwanProject(p_proj, n_proj)


# --------------------------------------
# SWAN main mesh

# depth grid description (input bathymetry grid)
sp.mesh_main.dg = {
    'xpc': 0,      # x origin
    'ypc': 0,      # y origin
    'alpc': 0,     # x-axis direction 
    'xlenc': 400,  # grid length in x
    'ylenc': 400,  # grid length in y
    'mxc': 1,      # number mesh x
    'myc': 1,      # number mesh y
    'dxinp': 400,  # size mesh x
    'dyinp': 400,  # size mesh y
}

# depth value
sp.mesh_main.depth = np.ones((2,2)) * 155

# computational grid description
sp.mesh_main.cg = {
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


# --------------------------------------
# SWAN nest1 mesh (allows nest1, nest2 and nest3)

# depth grid description
sp.mesh_nest1.dg = {
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
sp.mesh_nest1.depth = np.ones((10,8)) * 158

# computational grid description
sp.mesh_nest1.cg = {
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

# activate nest1
sp.run_nest1 = True


# --------------------------------------
# SWAN parameters (sea level, jonswap gamma)

sp.params = {
    'sea_level': 4,
    'jonswap_gamma': 1.9,
    'coords_spherical': None,   # None, 'GCM', 'CCM' 
    'waves_period': 'MEAN',     # 'PEAK / MEAN'
    'nested_bounds': 'CLOSED',  # 'CLOSED' / 'OPEN'
}


# --------------------------------------
# SWAN wrap STAT (create case files, launch SWAN num. model, extract output)

sw = SwanWrap_STAT(sp)

# build stationary cases from waves data
sw.build_cases(we)

# run SWAN
sw.run_cases()

# extract output from main mesh 
xds_out_main = sw.extract_output()
print('\noutput main mesh')
print(xds_out_main)
print()

# extract output from nest1 mesh 
xds_out_nest1 = sw.extract_output(mesh=sp.mesh_nest1)
print('\noutput nest1 mesh')
print(xds_out_nest1)

