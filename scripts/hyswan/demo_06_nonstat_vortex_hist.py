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
from hywaves.swan.storms import historic_track_preprocessing, historic_track_interpolation
from hywaves.swan.wrap import SwanProject, SwanMesh, SwanWrap_NONSTAT


# --------------------------------------

# data
p_data = op.abspath(op.join(op.dirname(__file__), '..', '..', 'data'))
p_demo = op.join(p_data, 'demo')

# bathymetry (from .nc file)
p_bathy = op.join(p_demo, 'depth_majuro.nc')

xds_bathy = xr.open_dataset(p_bathy)

# sign convention [0º,360º]
xds_lon = xds_bathy.lon.values
xds_lon[xds_lon<0] = xds_lon[xds_lon<0] + 360
xds_bathy.lon.values = xds_lon

lon = xds_bathy.lon.values[:]
lat = xds_bathy.lat.values[:]

depth = xds_bathy.elevation.values[:] * -1  # elevation to depth 

# shoreline (from .nc file)
p_shore = op.join(p_demo, 'shore_majuro.npy')

np_shore = np.load(p_shore)

# sign convention [0º,360º]
lon_shore = np_shore[:,0]
lon_shore[lon_shore<0] = lon_shore[lon_shore<0] + 360
np_shore[:,0] = lon_shore


# historic storm
storm = xr.open_dataset(op.join(p_demo, 'storm_ibtracs_paka.nc'))

# target coordinates
target = 'Kwajalein'
x0, y0 = 167.73, 8.72
if x0 < 0:  x0 = x0 + 360   # sign convention [0º,360º]



# --------------------------------------
# SWAN project (config bathymetry, parameters, computational grid)

p_proj = op.join(p_data, 'projects')  # swan projects main directory
n_proj = '03_vortex_hist'                  # project name

sp = SwanProject(p_proj, n_proj)


# --------------------------------------
# SWAN main mesh
main_mesh = SwanMesh()

# depth grid description (input bathymetry grid)
# grid resolution of 15km (=0.136º)
res = 0.136
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
    'xpc': 163.5,
    'ypc': 0.5,
    'alpc': 0,
    'xlenc': 15,
    'ylenc': 13,
    'mxc': int(round(15/res)),    # grid resolution of 15km (=0.136º)
    'myc': int(round(13/res)),
    'dxinp': 15/int(round(15/res)),
    'dyinp': 13/int(round(13/res)),
}

sp.set_main_mesh(main_mesh)



# --------------------------------------
# SWAN nested meshes


# NEST 1
mesh_nest1 = SwanMesh()

# depth grid description (input bathymetry grid)
# grid resolution of 5km (=0.0453º)
res1 = 0.04533
mesh_nest1.dg = main_mesh.dg

# depth value (from file)
mesh_nest1.depth = main_mesh.depth

# computational grid description
mesh_nest1.cg = {
    'xpc': 168.5,
    'ypc': 5.5,
    'alpc': 0,
    'xlenc': 5.5,
    'ylenc': 3.5,
    'mxc': int(round(5.5/res1)),
    'myc': int(round(3.5/res1)),
    'dxinp': 5.5/int(round(5.5/res1)),
    'dyinp': 3.5/int(round(3.5/res1)),
}


# NEST 2
mesh_nest2 = SwanMesh()

# depth grid description (input bathymetry grid)
# grid resolution of 1km (=0.009º)
res2 = 0.009
mesh_nest2.dg = main_mesh.dg

# depth value (from file)
mesh_nest2.depth = main_mesh.depth

# computational grid description
mesh_nest2.cg = {
    'xpc': 170.9,
    'ypc': 6.8,
    'alpc': 0,
    'xlenc': 1.2,
    'ylenc': 0.7,
    'mxc': int(round(1.2/res2)),
    'myc': int(round(0.7/res2)),
    'dxinp': 1.2/int(round(1.2/res2)),
    'dyinp': 0.7/int(round(0.7/res2)),
}


# set project nested mesh list
sp.set_nested_mesh_list([mesh_nest1, mesh_nest2])


# --------------------------------------
# SWAN parameters (sea level, jonswap gamma)
sp.params = {
    'sea_level': 0,
    'jonswap_gamma': 3.3,
    'cdcap': 2.5*10**-3,
    'coords_spherical': 'GCM',
    'waves_period': 'MEAN',
    'maxerr': None,
}


# SWAN output points
sp.x_out = [172.5, 172.5, 171]
sp.y_out = [8.5, 9.6, 7.4]

# add shoreline for plotting and check-ups
sp.np_shore = np_shore
sp.np_shore_nested = [np_shore, np_shore]

# computation time step (equal to the storm dt interpolation)
sp.dt_comp = 30
sp.dt_comp_nested = 20
# TODO: si se define dt_comp_nested hay que generar st_nest con ese paso de tiempo

## 
#sp.vortex_max = None
#sp.vortex_min = None
#
## save plots
#sp.plot_vortex = None
#sp.plot_track = None


# --------------------------------------
# SWAN case input: waves_event (empty) + storm_track (from historic storm)

# variables names for v3.10 Ibtracs
d_vns = {
    'longitude':'lon_wmo',
    'latitude':'lat_wmo',
    'time': 'time_wmo',
    'pressure':'pres_wmo',
    'maxwinds':'wind_wmo',
}

# preprocess storm variables
st_time, ylat_tc, ylon_tc, ycpres, ywind, ts, categ = historic_track_preprocessing(storm, d_vns)

# generate interpolated storm track  
st, time_input = historic_track_interpolation(st_time, ylon_tc, ylat_tc, ycpres, ywind,
                                              y0, x0, lat[0], lon[0], lat[1], lon[1],  
                                              ts, sp.dt_comp, wind=ywind,
                                              great_circle=True)

# generate wave event (empty)
we = pd.DataFrame(index=time_input, columns=['hs', 't02', 'dir', 'spr', 'U10', 'V10'])
we['level'] = 0
we['tide'] = 0

print('\ninput storm track')
print(st)


# --------------------------------------
# SWAN wrap NONSTAT (create case files, launch SWAN num. model, extract output)

sw = SwanWrap_NONSTAT(sp)

# build non-stationary cases from wave_events list and storm_tracks list
sw.build_cases([we], storm_track_list=[st], make_waves=False)

# run SWAN
sw.run_cases()

sys.exit()
# extract grid output from non-stationary cases
xds_out_main = sw.extract_output()
print('\noutput main mesh')
print(xds_out_main)
xds_out_main.to_netcdf(op.join(sw.proj.p_cases, 'output_nonstat.nc'))

xds_out_nest1 = sw.extract_output(mesh = sp.mesh_nest1)
print('\noutput nest1 mesh')
print(xds_out_nest1)
xds_out_nest1.to_netcdf(op.join(sw.proj.p_cases, 'output_nonstat_nest1.nc'))

xds_out_nest2 = sw.extract_output(mesh = sp.mesh_nest2)
print('\noutput nest2 mesh')
print(xds_out_nest2)
xds_out_nest2.to_netcdf(op.join(sw.proj.p_cases, 'output_nonstat_nest2.nc'))

# extract point output from non-stationary cases
xds_out_pts = sw.extract_output_points()
# TODO: to integrate mesh


