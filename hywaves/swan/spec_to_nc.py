#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os.path as op
import netCDF4 as nc4
from datetime import datetime
import numpy as np
import xarray as xr


def specdat2nc(p_dat, p_nc):
    'Parse spectral output .dat file to .nc file'

    # get time dimension length 
    with open(p_dat) as f:
        len_time = sum(l.count("date and time") for l in f)

    with open(p_dat) as f:

        # skip header
        for i in range(3): f.readline()

        # skip TIME lines
        for i in range(2): f.readline()

        # skip LONLAT line
        f.readline()

        # get number of points and coordinates
        n_points = int(f.readline().split()[0])
        lon_points, lat_points = [], []
        for i in range(n_points):
            lonlat = f.readline()
            lon_points.append(float(lonlat.split()[0]))
            lat_points.append(float(lonlat.split()[1]))

        # skip AFREQ line
        f.readline()

        # get FREQs
        n_freq = int(f.readline().split()[0])
        freqs = []
        for i in range(n_freq):
            freqs.append(float(f.readline()))

        # skip NDIR line
        f.readline()

        # get DIRs
        n_dir = int(f.readline().split()[0])
        dirs = []
        for i in range(n_dir):
            dirs.append(float(f.readline()))

        # skip QUANT lines
        for i in range(2): f.readline()

        # skip vadens lines
        for i in range(2): f.readline()

        # get exception value
        ex = float(f.readline().split()[0])

        # start reading spectral output data
        times = []

        # generate output netcdf4 file
        ncf = nc4.Dataset(p_nc, mode='w')

        ncf.createDimension('time', len_time)
        ncf.createDimension('freq', len(freqs))
        ncf.createDimension('dir', len(dirs))
        ncf.createDimension('point', len(lon_points))

        ncv_time = ncf.createVariable('time', np.float32, ('time',))
        ncv_freq = ncf.createVariable('freq', np.float32, ('freq',))
        ncv_dir = ncf.createVariable('dir', np.float32, ('dir',))
        ncv_lat = ncf.createVariable('lat_pts', np.float32, ('point',))
        ncv_lon = ncf.createVariable('lon_pts', np.float32, ('point',))
        ncv_spec= ncf.createVariable('spec', np.float32, ('time', 'freq', 'dir', 'point'), zlib=True, complevel=9)

        ncv_time.units = 'hours since 1800-01-01'
        ncv_freq.units = 'Hz'
        ncv_dir.units = 'degr'
        #ncv_lat.units = 'degr'
        #ncv_lon.units = 'degr'
        ncv_spec.units = 'm2/Hz/degr'

        # read time by time
        ct = 0
        tl = f.readline()
        while tl:
            print(tl.split()[0])
            time = datetime.strptime(tl.split()[0], '%Y%m%d.%H%M%S')
            times.append(time)

            # spectral output numpy storage
            spec_pts_t = np.ones((n_freq, n_dir, n_points)) * np.nan

            # read all points
            for p in range(n_points):

                # check we have result
                fac_line = f.readline()

                # read spectra values 
                if 'FACTOR' in fac_line:

                    # get factor
                    fac = float(f.readline())

                    # read matrix
                    for i in range(n_freq):
                        spec_pts_t[i,:,p] = np.fromstring(f.readline(), sep=' ')

                    # multiply spec by factor
                    spec_pts_t[:,:,p] = spec_pts_t[:,:,p] * fac

            # append spectra to .nc file
            ncv_spec[ct, :, :, :] = spec_pts_t

            # read next time line (if any)
            tl = f.readline()
            ct += 1

        # fill .nc data
        ncv_time[:] = nc4.date2num(times, ncv_time.units)
        ncv_lat[:] = lat_points
        ncv_lon[:] = lon_points
        ncv_freq[:] = freqs
        ncv_dir[:] = dirs

        # close .nc file
        ncf.close()


# ---------------------------------

# test folder
p_case = '/Users/admin/Projects/BlueMath/HybridModels/HyWaves/temp_quarantine/spec_9g'

# .dat input file
p_dat = op.join(p_case, 'spec_compgrid_main.dat')

# .nc output file
p_nc = op.join(p_case, 'spec_compgrid_main.nc')

# parse file
specdat2nc(p_dat, p_nc)


# load file with xarray
specs = xr.open_dataset(p_nc)
print(specs)


