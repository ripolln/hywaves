#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import shutil as su
from datetime import datetime

import numpy as np
import pandas as pd
import xarray as xr
from scipy.io import loadmat
from scipy import interpolate

from .geo import geo_distance_azimuth


# SWAN INPUT/OUTPUT STAT LIBRARY

class SwanIO(object):
    'SWAN numerical model input/output'

    def __init__(self, swan_proj):

        # needs SwanProject 
        self.proj = swan_proj

    def make_project(self):
        'makes swan project folder and subfolders'

        if not op.isdir(self.proj.p_main): os.makedirs(self.proj.p_main)
        if not op.isdir(self.proj.p_cases): os.makedirs(self.proj.p_cases)


class SwanIO_STAT(SwanIO):
    'SWAN numerical model input/output - STATIONARY cases'

    def make_input(self, p_file, id_run, ws, bnd):
        '''
        Writes input.swn file from waves sea state for stationary execution

        p_file  - input.swn file path
        ws      - wave sea state (hs, per, dr, spr)
        bnd     - wave sea state active boundaries

        more info: http://swanmodel.sourceforge.net/online_doc/swanuse/node23.html
        '''
        # TODO: check readinp idla

        # .swn file parameters
        sea_level = self.proj.params['sea_level']
        jonswap_gamma = self.proj.params['jonswap_gamma']
        coords_spherical = self.proj.params['coords_spherical']
        waves_period = self.proj.params['waves_period']

        # main mesh
        mm = self.proj.mesh_main

        # .swn text file
        t = "PROJ '{0}' '{1}'\n$\n".format(self.proj.name, id_run)
        t += 'MODE STAT\n'

        # spherical coordinates (mercator) switch
        if coords_spherical != None:
            t += 'COORDINATES SPHER {0}\n'.format(coords_spherical)

        # sea level
        t += 'SET level={0}  NAUTICAL\n$\n'.format(sea_level)

        # computational grid
        t += 'CGRID REGULAR {0} {1} {2} {3} {4} {5} {6} CIRCLE 72 0.0345 1.00  34\n$\n'.format(
            mm.cg['xpc'], mm.cg['ypc'], mm.cg['alpc'], mm.cg['xlenc'],
            mm.cg['ylenc'], mm.cg['mxc']-1, mm.cg['myc']-1)

        # bathymetry
        t += 'INPGRID BOTTOM REGULAR {0} {1} {2} {3} {4} {5} {6}\n'.format(
            mm.dg['xpc'], mm.dg['ypc'], mm.dg['alpc'], mm.dg['mxc'],
            mm.dg['myc'], mm.dg['dxinp'], mm.dg['dyinp'])

        t += "READINP BOTTOM 1 '{0}' {1} 0 FREE\n$\n".format(
            mm.depth_fn, mm.dg_idla)

        # waves boundary conditions
        t += 'BOUND SHAPespec JONswap {0} {1} DSPR DEGR\n'.format(
            jonswap_gamma, waves_period)
        for ic in bnd:
            t += "BOUN SIDE {0} CONstant PAR {1:.3f} {2:.3f} {3:.3f} {4:.3f}\n".format(
                ic, ws.hs, ws.per, ws.dir, ws.spr)
        t += "$\n"

        # numerics
        t += 'OFF QUAD\n'
        # t += 'PROP BSBT\n'
        # t += 'WCAP\n'
        t += 'BREA\n'
        t += 'FRICTION JONSWAP\n$\n'

        # optional nested mesh
        r_ns = [self.proj.run_nest1, self.proj.run_nest2, self.proj.run_nest3]
        m_ns = [self.proj.mesh_nest1, self.proj.mesh_nest2, self.proj.mesh_nest3]
        nout_0 = ['nest1', 'nest2', 'nest3']
        nout_1 = ['bounds_nest1.dat', 'bounds_nest2.dat', 'bounds_nest3.dat']

        for r_n, m_n, n0, n1 in zip(r_ns, m_ns, nout_0, nout_1):
            if r_n:
                t += "NGRID '{0}' {1} {2} {3} {4} {5} {6} {7}\n".format(
                    n0, m_n.cg['xpc'], m_n.cg['ypc'], m_n.cg['alpc'],
                    m_n.cg['xlenc'], m_n.cg['ylenc'],
                    np.int32(m_n.cg['xlenc']/mm.cg['dxinp']),
                    np.int32(m_n.cg['ylenc']/mm.cg['dyinp'])
                )
                t += "NESTOUT '{0}' '{1}'\n".format(n0, n1)

        # output
        t += "BLOCK 'COMPGRID' NOHEAD '{0}' LAY 3 HSIGN TM02 DIR TPS DSPR\n$\n".format(
            mm.output_fn,
        )

        # compute
        t += 'TEST  1,0\n'
        t += 'COMPUTE \n'
        t += 'STOP\n$\n'

        # write file:
        with open(p_file, 'w') as f:
            f.write(t)

        # log    
        fmt2 = ' 7.2f'
        print(
            'SWAN CASE: {1} ---> hs {2:{0}}, per {3:{0}}, dir {4:{0}}, spr {5:{0}}'.format(
                fmt2, id_run, ws.hs, ws.per, ws.dir, ws.spr
            )
        )

    def make_input_nested(self, p_file, id_run):
        '''
        Writes input_nested.swn file from waves sea state for stationary execution

        p_file  - input_nestedN.swn file path
        '''

        # TODO check myc-1, mxc -1 

        # .swn file parameters
        sea_level = self.proj.params['sea_level']
        coords_spherical = self.proj.params['coords_spherical']
        nested_bounds = self.proj.params['nested_bounds']

        # SWAN nested Computacional grid
        mn1 = self.proj.mesh_nest1


        # .swn text file
        t = "PROJ '{0}' '{1}'\n$\n".format(self.proj.name, id_run)
        t += 'MODE STAT\n'

        # spherical coordinates (mercator) switch
        if coords_spherical != None:
            t += 'COORDINATES SPHER {0}\n'.format(coords_spherical)

        t += 'SET level={0}  NAUTICAL\n$\n'.format(sea_level)

        # computational grid
        t += 'CGRID REGULAR {0} {1} {2} {3} {4} {5} {6} CIRCLE 72 0.03558410 1.00  35\n$\n'.format(
            mn1.cg['xpc'], mn1.cg['ypc'], mn1.cg['alpc'], mn1.cg['xlenc'],
            mn1.cg['ylenc'], mn1.cg['mxc']-1, mn1.cg['myc']-1)

        # bathymetry
        t += 'INPGRID BOTTOM REGULAR {0} {1} {2} {3} {4} {5} {6}\n'.format(
            mn1.dg['xpc'], mn1.dg['ypc'], mn1.dg['alpc'], mn1.dg['mxc']-1,
            mn1.dg['myc']-1, mn1.dg['dxinp'], mn1.dg['dyinp'])

        t += "READINP BOTTOM 1 '{0}' {1} 0 FREE\n$\n".format(
            mn1.depth_fn, mn1.dg_idla)

        # Boundary Conditions
        t += "BOUN NEST '{0}' {1}\n".format('bounds_nest1.dat', nested_bounds)

        #  wind file
        t += "$\n"

        # numerics
        t += 'OFF QUAD\n'
        # t += 'GEN1\n'
        # t += 'PROP BSBT\n'
        # t += 'WCAP\n'
        t += 'BREA\n'
        t += 'FRICTION JONSWAP\n$\n'

        # output
        t += "BLOCK 'COMPGRID' NOHEAD '{0}' LAY 3 HSIGN TM02 DIR TPS DSPR\n$\n".format(
            mn1.output_fn,
        )

        # compute
        t += 'TEST  1,0\n'
        t += 'COMPUTE \n'
        t += 'STOP\n$\n'

        # write file:
        with open(p_file, 'w') as f:
            f.write(t)

    def build_case(self, case_id, waves_ss, bnd=['N', 'E', 'W', 'S']):
        '''
        Build SWAN STAT case input files for given wave sea state (hs, per, dir, spr)

        ix_case  - SWAN case index (int)
        waves_ss - wave sea state (hs, per, dir, spr)
        bnd      - wave sea state active boundaries
        '''

        # SWAN case path
        p_case = op.join(self.proj.p_cases, case_id)

        # make execution dir
        if not op.isdir(p_case): os.makedirs(p_case)

        # make depth file for main mesh
        self.proj.mesh_main.export_dat(p_case)

        # make input.swn file
        self.make_input(op.join(p_case, 'input.swn'), case_id, waves_ss, bnd)

        # optional nested mesh depth and input files
        r_ns = [self.proj.run_nest1, self.proj.run_nest2, self.proj.run_nest3]
        m_ns = [self.proj.mesh_nest1, self.proj.mesh_nest2, self.proj.mesh_nest3]
        i_ns = ['input_nest1.swn', 'input_nest2.swn', 'input_nest3.swn']

        for r_n, m_n, i_n in zip(r_ns, m_ns, i_ns):
            if r_n:
                m_n.export_dat(p_case)
                self.make_input_nested(op.join(p_case, i_n), case_id)

    def outmat2xr(self, p_mat):

        # matlab dictionary
        dmat = loadmat(p_mat)

        # return dataset
        xds_out = xr.Dataset(
            {
                'Hsig':   (('X','Y',), dmat['Hsig'].T,   {'units':'m'}),
                'Tm02':   (('X','Y',), dmat['Tm02'].T,   {'units':'s'}),
                'Dir':    (('X','Y',), dmat['Dir'].T,    {'units':'º'}),
                'Dspr':    (('X','Y',), dmat['Dspr'].T,  {'units':'º'}),
                'TPsmoo': (('X','Y',), dmat['TPsmoo'].T, {'units':'s'}),
            }
        )

        return xds_out

    def output_case(self, p_case, mesh):
        'read .mat output file from stationary and returns xarray.Dataset'

        # extract output from selected mesh
        p_mat = op.join(p_case, mesh.output_fn)
        xds_out = self.outmat2xr(p_mat)

        # set X and Y values
        X, Y = mesh.get_XY()
        xds_out = xds_out.assign_coords(X=X)
        xds_out = xds_out.assign_coords(Y=Y)

        # rename to longitude latitude in spherical coords cases
        coords_spherical = self.proj.params['coords_spherical']
        if coords_spherical != None:
            xds_out = xds_out.rename({'X':'lon', 'Y':'lat'})

        return xds_out


class SwanIO_NONSTAT(SwanIO):
    'SWAN numerical model input/output - NON STATIONARY cases'

    def make_out_points(self, p_file):
        'Generates desired output-points coordinates file'

        # define and save output points
        x_out = self.proj.x_out
        y_out = self.proj.y_out

        if not x_out or not y_out:
            return

        else:
            points = np.vstack((x_out,y_out)).T
            np.savetxt(p_file, points, fmt='%.2f')

    def make_wave_files(self, p_case, waves_event, time, bnd):
        'Generate event wave files (swan compatible)'

        # wave variables
        hs = waves_event.hs.values[:]
        per = waves_event.per.values[:]
        direc = waves_event.dir.values[:]
        spr = waves_event.spr.values[:]

        # csv file 
        num_data = len(time)
        data = np.zeros((num_data, 5))
        data[:, 0] = time
        data[:, 1] = hs
        data[:, 2] = per
        data[:, 3] = direc
        data[:, 4] = spr

        # Copy file for all boundaries
        save = op.join(p_case, 'series_waves.dat')
        np.savetxt(save, data, header='TPAR', comments='', fmt='%8.4f %2.3f %2.3f %3.2f %3.1f')
        for i in bnd:
            su.copyfile(save, op.join(p_case, 'series_waves_{0}.dat'.format(i)))

    def make_wind_files(self, p_case, waves_event):
        '''
        Generate event wind mesh files (swan compatible)

        uses wave_event U10 and V10 values at the entire SWAN comp. grid
        '''

        # wind variables
        u10 = waves_event.U10.values[:]
        v10 = waves_event.V10.values[:]

        # main mesh
        mm = self.proj.mesh_main

        # each time needs 2D (mesh) wind files (U,V) 
        mxc = mm.cg['mxc']  # number mesh x
        myc = mm.cg['myc']  # number mesh y

        txt = ''
        for c, (u, v) in enumerate(zip(u10,v10)):

            # single point wind -> entire SWAN comp.grid wind
            aux = np.ones((mxc, myc))

            # TODO: wind has to be rotated if alpc != 0

            # csv file 
            u_2d = aux * u
            v_2d = aux * v
            u_v_stack = np.vstack((u_2d, v_2d))
            save = op.join(p_case, 'wind_{0:06}.dat'.format(c))
            np.savetxt(save, u_v_stack, fmt='%.2f')

            # wind list file
            txt += 'wind_{0:06}.dat\n'.format(c)

        # winds file path
        save = op.join(p_case, 'series_wind.dat')
        with open(save, 'w') as f:
            f.write(txt)

    def make_vortex_files(self, p_case, case_id, storm_track, 
                          nested=False, num_nest=None):
#    def make_vortex_files(self, p_case, storm_track):
        '''
        Generate event wind mesh files (swan compatible)

        uses wave_event storm path data over SWAN computational grid
        needs SPHERICAL COORDINATES
        '''

        # parameters
        RE = 6378.135 * 1000            # Earth radius [m]
        beta = 0.9                      # conversion factor of wind speed
        rho_air = 1.15                  # air density
        w = 2 * np.pi / 86184.2         # Earth's rotation velocity (rad/s)
        pifac = np.arccos(-1) / 180     # pi/180
        one2ten = 0.8928                # conversion from 1-min to 10 min
            
        # wind variables
        storm_move = storm_track.move.values[:]
#        storm_vf =   storm_track.vf.values[:]
        storm_vfx =  storm_track.vfx.values[:]
        storm_vfy =  storm_track.vfy.values[:]
        storm_lon =  storm_track.lon.values[:]
        storm_lat =  storm_track.lat.values[:]
        storm_p0 =   storm_track.p0.values[:]
        storm_pn =   storm_track.pn.values[:]
        times =      storm_track.index[:]
        storm_vmax = storm_track.vmax.values[:]       
        storm_lat_orig = storm_lat

        # main mesh
        mm = self.proj.mesh_main

        # comp. grid for generating vortex wind files
        mxc = mm.cg['mxc']  # number mesh x
        myc = mm.cg['myc']  # number mesh y

        # comp. grid lat, lon limits 
        lon0 = mm.cg['xpc']
        lat0 = mm.cg['ypc']
        lon1 = mm.cg['xpc'] + mm.cg['xlenc']
        lat1 = mm.cg['ypc'] + mm.cg['ylenc']

        # general meshgrid limits
        lon00, lat00, lon01, lat01 = lon0, lat0, lon1, lat1       
        
#        if nested == True:
#            # comp. grid for generating vortex wind files
#            mxc = self.proj.cc_mxn[num_nest]        # number mesh x
#            myc = self.proj.cc_myn[num_nest]        # number mesh y
#    
#            # comp. grid lat, lon limits 
#            lon0 = self.proj.cc_xpn[num_nest]
#            lat0 = self.proj.cc_ypn[num_nest]
#            lon1 = self.proj.cc_xpn[num_nest] + self.proj.cc_xlenn[num_nest]
#            lat1 = self.proj.cc_ypn[num_nest] + self.proj.cc_ylenn[num_nest]            

        cg_lon = np.linspace(lon0, lon1, mxc)
        cg_lat = np.linspace(lat0, lat1, myc)
        mg_lon, mg_lat = np.meshgrid(cg_lon, cg_lat)

        # wind output holder
        hld_W = np.zeros((len(cg_lat), len(cg_lon), len(storm_move)))
        hld_D = np.zeros((len(cg_lat), len(cg_lon), len(storm_move)))

        # Correction when track is in south hemisphere for vortex generation 
        if any (i < 0 for i in storm_lat) == True:
                south_hemisphere = True
        else:   south_hemisphere = False

        # each time needs 2D (mesh) wind files (U,V)
        txt = ''
        for c, (lo, la, la_orig, p0, pn, ut, vt, vmax) in enumerate(zip(
            storm_lon, storm_lat, storm_lat_orig, storm_p0, storm_pn, 
            storm_vfx, storm_vfy, storm_vmax)):
            
            # Wind model code (from ADCIRC, transcribed by Antonio Espejo) and 
            # later slightly modified by Sara Ortega to include TCs at southern 
            # hemisphere

            # get distance and angle between points 
            arcl, theta = geo_distance_azimuth(mg_lat, mg_lon, la, lo)
            r = arcl * np.pi / 180.0 * RE
            
            # angle correction for southern hemisphere
            if south_hemisphere:        
                thet = np.arctan2((mg_lat-la)*pifac, -(mg_lon-lo)*pifac)
            if not south_hemisphere:    
                thet = np.arctan2((mg_lat-la)*pifac, (mg_lon-lo)*pifac)
            
            # ADCIRC MODEL 
            CPD = (pn - p0) * 100    # central pressure deficit [Pa]
            if CPD < 100: CPD = 100  # limit central pressure deficit
            
            # Wind model 
            f = 2 * w * np.sin(abs(la)*np.pi/180)  # Coriolis
            
            # Substract the translational storm speed from the observed maximum 
            # wind speed to avoid distortion in the Holland curve fit. 
            # The translational speed will be added back later
            vkt = vmax - np.sqrt(np.power(ut,2) + np.power(vt,2))   # [kt]
            
            # Convert wind speed from 10m altitude to wind speed at the top of 
            # the atmospheric boundary layer
            vgrad = vkt / beta    # [kt]
            v = vgrad 
            vm = vgrad * 0.52     # [m/s]
            
            # Knaff et al. (2016) - Radius of maximum wind (RMW)
            rm = 218.3784 - 1.2014*v + np.power(v/10.9844,2) - np.power(v/35.3052,3) - 145.509*np.cos(la*pifac)  # nautical mile
            rm = rm * 1.852 * 1000   # from [n mi] to [m]
            rn = rm / r              # []
            
            # Holland B parameter with upper and lower limits
            B = rho_air * np.exp(1) * np.power(vm,2) / CPD
            if B > 2.5: B = 2.5       
            elif B < 1: B = 1
            
            # Wind velocity at each node and time step   [m/s]
            vg = np.sqrt(np.power(rn,B) * np.exp(1-np.power(rn,B)) * np.power(vm,2) + np.power(r,2)*np.power(f,2)/4) - r*f/2
            
            # Determine translation speed that should be added to final storm  
            # wind speed. This is tapered to zero as the storm wind tapers to 
            # zero toward the eye of the storm and at long distances from the storm
            vtae = (abs(vg) / vgrad) * ut    # [m/s]
            vtan = (abs(vg) / vgrad) * vt
            
            # Find the velocity components and convert from wind at the top of the 
            # atmospheric boundary layer to wind at 10m elevation
            if south_hemisphere:        hemisphere_sign = 1
            if not south_hemisphere:    hemisphere_sign = -1
            ve = hemisphere_sign * vg * beta * np.sin(thet)     # [m/s]
            vn = vg * beta * np.cos(thet)
            
            # Convert from 1 minute averaged winds to 10 minute averaged winds
            ve = ve * one2ten    # [m/s]
            vn = vn * one2ten
            
            # Add the storm translation speed
            vfe = ve + vtae      # [m/s]
            vfn = vn + vtan
            
            # wind module
            W = np.sqrt(np.power(vfe,2) + np.power(vfn,2))    # [m/s]
            
            # Surface pressure field
            pr = p0 + (pn-p0) * np.exp(- np.power(rn,B))      # [mbar]
            py, px = np.gradient(pr)
            ang = np.arctan2(py, px) + np.sign(la_orig) * np.pi/2.0
            
            # Wind field
            u_2d = W * np.cos(ang)    # m/s
            v_2d = W * np.sin(ang)    # m/s
            u_v_stack = np.vstack((u_2d, v_2d))
            
#            if nested == True:
#                save = op.join(p_case, 'wind_nest_{0}_{1:06}.dat'.format(num_nest+1, c))
#                np.savetxt(save, u_v_stack, fmt='%.2f')
#                # wind list file
#                txt += 'wind_nest_{0}_{1:06}.dat\n'.format(num_nest+1, c)
#            else:
#                save = op.join(p_case, 'wind_{0:06}.dat'.format(c))
#                np.savetxt(save, u_v_stack, fmt='%.2f')
#                # wind list file
#                txt += 'wind_{0:06}.dat\n'.format(c)

            # csv file 
            save = op.join(p_case, 'wind_{0:06}.dat'.format(c))
            np.savetxt(save, u_v_stack, fmt='%.2f')

            # wind list file
            txt += 'wind_{0:06}.dat\n'.format(c)

            # hold wind data (m/s)
            hld_W[:,:,c] = W 
            hld_D[:,:,c] =  270 - np.rad2deg(ang)  # direction (º clock. rel. north)

#        # winds file path
#        if nested == True:
#            save = op.join(p_case, 'series_wind_nest_{0}.dat'.format(num_nest+1))
#            with open(save, 'w') as f:
#                f.write(txt)
#        else:
#            save = op.join(p_case, 'series_wind.dat')
#            with open(save, 'w') as f:
#                f.write(txt)

        # winds file path
        save = op.join(p_case, 'series_wind.dat')
        with open(save, 'w') as f:
            f.write(txt)

        # aux. save vortex wind fields
        if south_hemisphere == True:
            cg_lon = np.linspace(lon00, lon01, mxc)
            cg_lat = np.linspace(lat00, lat01, myc)

#        # save vortex dataset
#        xds_vortex = xr.Dataset(
#            {
#                'W': (('Y','X','time'), hld_W, {'units':'m/s'}),
#                'Dir': (('Y','X','time'), hld_D, {'units':'º'})
#            },
#            coords={
#                'time' : times,
##                'case' : int(case_id),
#                'Y' : cg_lat,
#                'X' : cg_lon,
#            }
#        )
#        xds_vortex.attrs['xlabel'] = 'Longitude (º)'
#        xds_vortex.attrs['ylabel'] = 'Latitude (º)'
#
#        if nested == True:
#            p_vortex = op.join(p_case, 'vortex_wind_nest_{0}.nc'.format(num_nest+1))
#            xds_vortex_nest = xds_vortex
#            xds_vortex_nest.to_netcdf(p_vortex)
#        else:
#            p_vortex = op.join(p_case, 'vortex_wind.nc')
#            xds_vortex.to_netcdf(p_vortex)
#            self.proj.vortex_max = float(xds_vortex['W'].max(axis=0).max())
#            self.proj.vortex_min = float(xds_vortex['W'].min(axis=0).min())

        # aux. save vortex wind fields
        p_vortex = op.join(p_case, 'vortex_wind.nc')
        xds_vortex = xr.Dataset(
            {
                'W': (('lat','lon','time'), hld_W, {'units':'m/s'}),
                'Dir': (('lat','lon','time'), hld_D, {'units':'º'})
            },
            coords={
                'Y' : cg_lat,
                'X' : cg_lon,
                'time' : times,
            }
        )
        xds_vortex.attrs['xlabel'] = 'Longitude (º)'
        xds_vortex.attrs['ylabel'] = 'Latitude (º)'
        xds_vortex.to_netcdf(p_vortex)
        
#        # save plot vortex wind files
#        if self.proj.plot_vortex == True:
#            if nested == True:
#                p_export = op.join(p_case, 'vortex_figs_nest_{0}'.format(num_nest+1))
#                p_vortex_figs = op.join(p_export, case_id)
#                plot_vortex_times(
#                    self.proj.storm, xds_vortex, 'W', p_vortex_figs, self.proj.num_nest, 
#                    self.proj.vortex_max, self.proj.vortex_min, 
#                    self.proj.lon_nested, self.proj.lat_nested, 
#                    cmap='hot_r', quiver=True, np_shore=np_shore)
#            else:
#                p_export = op.join(p_case, 'vortex_figs')
#                p_vortex_figs = op.join(p_export, case_id)
#                plot_vortex_times(
#                    self.proj.storm, xds_vortex, 'W', p_vortex_figs, self.proj.num_nest, 
#                    self.proj.vortex_max, self.proj.vortex_min, 
#                    self.proj.lon_nested, self.proj.lat_nested, 
#                    cmap='hot_r', quiver=True, np_shore=np_shore)
#
#        # plot storm track
#        if self.proj.plot_track == True:
#            plot_storm_track(self.proj.storm, storm_track, case_id, lon00, lon11, 
#                             lat00, lat11, self.proj.num_nest, self.proj.lon_nested, 
#                             self.proj.lat_nested, p_case, np_shore=np_shore)


    def make_level_files(self, p_case, wave_event):
        'Generate event level mesh files (swan compatible)'

        # parse pandas time index to swan iso format
        swan_iso_fmt = '%Y%m%d.%H%M'
        time = pd.to_datetime(wave_event.index).strftime(swan_iso_fmt).values[:]

        # level variables
        zeta = wave_event.level.values[:]
        tide = wave_event.tide.values[:]

        # main mesh
        mm = self.proj.mesh_main

        # each time needs 2D (mesh) level 
        mxc = mm.cg['mxc']  # number mesh x
        myc = mm.cg['myc']  # number mesh y

        txt = ''
        for c, (z, t) in enumerate(zip(zeta, tide)):

            # single point level -> entire SWAN comp.grid level
            aux = np.ones((mxc, myc)).T

            # csv file 
            l = z + t  # total level
            l_2d = aux * l
            save = op.join(p_case, 'level_{0:06}.dat'.format(c))
            np.savetxt(save, l_2d, fmt='%.2f')

            # level list file
            txt += 'level_{0:06}.dat\n'.format(c)

        # waves file path
        save = op.join(p_case, 'series_level.dat')
        with open(save, 'w') as f:
            f.write(txt)

    def make_input(self, p_file, id_run, time, make_waves=True,
                   make_winds=True, wvs_bnd=['N', 'E', 'W', 'S']):
        '''
        Writes input.swn file from waves event for non-stationary execution

        p_file  - input.swn file path
        time    - event time at swan iso format

        make_waves - activates waves input files generation (at waves_bnd)
        make_winds - activates wind input files generation

        more info: http://swanmodel.sourceforge.net/online_doc/swanuse/node23.html
        '''

        # event time (swan iso format)
        t0_iso = time[0]
        t1_iso = time[-1]

        # .swn file parameters
        sea_level = self.proj.params['sea_level']
        jonswap_gamma = self.proj.params['jonswap_gamma']
        cdcap = self.proj.params['cdcap']
        maxerr = self.proj.params['maxerr']
        coords_spherical = self.proj.params['coords_spherical']
        waves_period = self.proj.params['waves_period']

        # main mesh
        mm = self.proj.mesh_main

        # nested mesh
        if self.proj.nested:    mn = self.proj.meshes_nested[0]

        # output points
        x_out = self.proj.x_out
        y_out = self.proj.y_out

        # computational time step data (minutes)
        if type(self.proj.dt_comp) == list:   dt_comp = self.proj.dt_comp[int(float(id_run))]
        if type(self.proj.dt_comp) != list:   dt_comp = self.proj.dt_comp

        # .swn text file
        t = "PROJ '{0}' '{1}' 'GENERAL'\n$\n".format(self.proj.name, id_run)
        t += 'MODE NONSTAT\n'

        # spherical coordinates (mercator) swich
        if coords_spherical:
            t += 'COORDINATES SPHER CCM\n'

        # cdcap
        cdcap_str = ''
        if cdcap: cdcap_str = 'cdcap={0}'.format(cdcap)

        # max error (caution)
        maxerr_str = ''
        if maxerr: maxerr_str = 'maxerr={0}'.format(maxerr)

        # set level and cdcap (if available)
        t += 'SET level={0} {1} {2}  NAUTICAL\n$\n'.format(
            sea_level, cdcap_str, maxerr_str
        )

        # computational grid
        t += 'CGRID REGULAR {0} {1} {2} {3} {4} {5} {6} CIRCLE 72 0.03 1.00 \n$\n'.format(
            mm.cg['xpc'], mm.cg['ypc'], mm.cg['alpc'], mm.cg['xlenc'],
            mm.cg['ylenc'], mm.cg['mxc'], mm.cg['myc'])

        # bathymetry
        t += 'INPGRID BOTTOM REGULAR {0} {1} {2} {3} {4} {5} {6}\n'.format(
            mm.dg['xpc'], mm.dg['ypc'], mm.dg['alpc'], mm.dg['mxc'],
            mm.dg['myc'], mm.dg['dxinp'], mm.dg['dyinp'])

        t += "READINP BOTTOM 1 '{0}' {1} 0 FREE\n$\n".format(
            mm.depth_fn, mm.dg_idla)

        # wind
        t += 'INPGRID  WIND  REGULAR {0} {1} {2} {3} {4} {5} {6} NONSTAT {7} {8} MIN {9}\n'.format(
            mm.cg['xpc'], mm.cg['ypc'], mm.cg['alpc'], mm.cg['mxc'],
            mm.cg['myc'], mm.cg['dxinp'], mm.cg['dyinp'], t0_iso, dt_comp, t1_iso)
        t += "READINP  WIND 1. SERIES '{0}' 3 0 FREE\n$\n".format('series_wind.dat')

        # level
        t += 'INPGRID  WLEV  REGULAR {0} {1} {2} {3} {4} {5} {6} NONSTAT {7} {8} MIN {9}\n'.format(
            mm.cg['xpc'], mm.cg['ypc'], mm.cg['alpc'], mm.cg['mxc'],
            mm.cg['myc'], mm.cg['dxinp'], mm.cg['dyinp'], t0_iso, dt_comp, t1_iso)
        t += "READINP  WLEV 1. SERIES '{0}' 3 0 FREE\n$\n".format('series_level.dat')

        # waves boundary conditions
        if make_waves:
            t += 'BOUND SHAPespec JONswap {0} {1} DSPR DEGR\n'.format(
                jonswap_gamma, waves_period)
            for ic in wvs_bnd:
                t += "BOUN SIDE {0} CONstant FILE 'series_waves_{0}.dat'\n".format(ic)

        # numerics & physics
        t += 'WIND DRAG WU\n'
        t += 'GEN3 ST6 5.7E-7 8.0E-6 4.0 4.0 UP HWANG VECTAU TRUE10\n'
        t += 'QUAD iquad=8\n'
        t += 'WCAP\n'
        if not coords_spherical:
            t += 'SETUP\n'  # not compatible with spherical 
        t += 'BREA\n'
        t += 'FRICTION JONSWAP\n$\n'
        t += 'TRIADS\n'
#        t += 'DIFFRAC\n'

        # numerics 
        t += 'PROP BSBT\n$\n'

        # optional nested mesh
        if self.proj.nested:
            t += "NGRID '{0}' {1} {2} {3} {4} {5} {6} {7}\n".format(
                    'nest1', mn.cg['xpc'], mn.cg['ypc'], mn.cg['alpc'], 
                    mn.cg['xlenc'], mn.cg['ylenc'], mn.cg['mxc'], mn.cg['myc'])
            t += "NESTOUT '{0}' '{1}' OUT {2} {3} MIN\n$\n".format(
                    'nest1', 'bounds_nest1.dat', t0_iso, dt_comp)
        
        # output
        t += "BLOCK 'COMPGRID' NOHEAD '{0}' LAY 3 HSIGN TM02 DIR TPS DSPR OUT {1} {2} MIN\n$\n".format(
            mm.output_fn, t0_iso, dt_comp)

        # output points
        if not x_out or not y_out:
            pass
        else:
            t += "POINTS 'outpts' FILE 'points_out.dat'\n"
            t += "TABLE 'outpts' NOHEAD 'table_outpts.dat' DEP HS HSWELL DIR RTP TM02 DSPR WIND WATLEV  OUT {0} {1} MIN\n$\n".format(t0_iso, dt_comp)

        # compute
        t += 'TEST  1,0\n'
        t += 'COMPUTE NONSTAT {0} {1} MIN {2}\n'.format(t0_iso, dt_comp, t1_iso)
        t += 'STOP\n$\n'

        # write file:
        with open(p_file, 'w') as f:
            f.write(t)

    def make_input_nested(self, p_file, id_run, time, num_nest):
        '''
        Writes input_nested.swn file from waves sea state for non-stationary execution

        p_file  - input_nestedN.swn file path
        time    - event time at swan iso format
        '''

        # event time (swan iso format)
        t0_iso = time[0]
        t1_iso = time[-1]

        # .swn file parameters
        sea_level = self.proj.params['sea_level']
#        jonswap_gamma = self.proj.params['jonswap_gamma']
        cdcap = self.proj.params['cdcap']
        maxerr = self.proj.params['maxerr']
        coords_spherical = self.proj.params['coords_spherical']
#        waves_period = self.proj.params['waves_period']

        # main mesh
        mm = self.proj.meshes_nested[num_nest]

        # nested mesh
        if num_nest < self.proj.num_nest-1:
            mn = self.proj.meshes_nested[num_nest+1]

        # output points
        x_out = self.proj.x_out
        y_out = self.proj.y_out

        # computational time step data (minutes)
        if type(self.proj.dt_comp_nested) == list:   dt_comp = self.proj.dt_comp_nested[int(float(id_run))]
        if type(self.proj.dt_comp_nested) != list:   dt_comp = self.proj.dt_comp_nested

        # .swn text file
        t = "PROJ '{0}' '{1}' 'NESTED_{2}'\n$\n".format(self.proj.name, id_run, num_nest+1)
        t += 'MODE NONSTAT\n'

        # spherical coordinates (mercator) swich
        if coords_spherical:
            t += 'COORDINATES SPHER CCM\n'

        # cdcap
        cdcap_str = ''
        if cdcap: cdcap_str = 'cdcap={0}'.format(cdcap)

        # max error (caution)
        maxerr_str = ''
        if maxerr: maxerr_str = 'maxerr={0}'.format(maxerr)

        # set level and cdcap (if available)
        t += 'SET level={0} {1} {2}  NAUTICAL\n$\n'.format(
            sea_level, cdcap_str, maxerr_str
        )

        # computational grid
        t += 'CGRID REGULAR {0} {1} {2} {3} {4} {5} {6} CIRCLE 72 0.03 1.00 \n$\n'.format(
            mm.cg['xpc'], mm.cg['ypc'], mm.cg['alpc'], mm.cg['xlenc'],
            mm.cg['ylenc'], mm.cg['mxc'], mm.cg['myc'])

        # bathymetry
        t += 'INPGRID BOTTOM REGULAR {0} {1} {2} {3} {4} {5} {6}\n'.format(
            mm.dg['xpc'], mm.dg['ypc'], mm.dg['alpc'], mm.dg['mxc'],
            mm.dg['myc'], mm.dg['dxinp'], mm.dg['dyinp'])

        t += "READINP BOTTOM 1 '{0}' {1} 0 FREE\n$\n".format(
            mm.depth_fn, mm.dg_idla)

        # wind
        t += 'INPGRID  WIND  REGULAR {0} {1} {2} {3} {4} {5} {6} NONSTAT {7} {8} MIN {9}\n'.format(
            mm.cg['xpc'], mm.cg['ypc'], mm.cg['alpc'], mm.cg['mxc'],
            mm.cg['myc'], mm.cg['dxinp'], mm.cg['dyinp'], t0_iso, dt_comp, t1_iso)
        t += "READINP  WIND 1. SERIES '{0}' 3 0 FREE\n$\n".format(
                'series_wind_nest{0}.dat'.format(num_nest+1))

        # level
        t += 'INPGRID  WLEV  REGULAR {0} {1} {2} {3} {4} {5} {6} NONSTAT {7} {8} MIN {9}\n'.format(
            mm.cg['xpc'], mm.cg['ypc'], mm.cg['alpc'], mm.cg['mxc'],
            mm.cg['myc'], mm.cg['dxinp'], mm.cg['dyinp'], t0_iso, dt_comp, t1_iso)
        t += "READINP  WLEV 1. SERIES '{0}' 3 0 FREE\n$\n".format(
                'series_level_nest{0}.dat'.format(num_nest+1))

        # no waves input for nested SWAN runs

        # Boundary Conditions for previous meshgrids
        t += "BOUN NEST '{0}'\n".format('bounds_nest{0}.dat'.format(num_nest+1))

        # physics
        t += 'WIND DRAG WU\n'
        t += 'GEN3 ST6 5.7E-7 8.0E-6 4.0 4.0 UP HWANG VECTAU TRUE10\n'
        t += 'QUAD\n'
        t += 'WCAP\n'
        if not coords_spherical:
            t += 'SETUP\n'  # not compatible with spherical 
        t += 'BREA\n'
        t += 'FRICTION JONSWAP\n'
        t += 'TRIADS\n'
#        t += 'DIFFRAC\n'
        
        # numerics 
        t += 'PROP BSBT\n$\n'

        # optional nested mesh
        if num_nest < self.proj.num_nest-1:
            t += "NGRID '{0}' {1} {2} {3} {4} {5} {6} {7}\n".format(
                    'nest{0}'.format(num_nest+2), mn.cg['xpc'], mn.cg['ypc'], mn.cg['alpc'], 
                    mn.cg['xlenc'], mn.cg['ylenc'], mn.cg['mxc'], mn.cg['myc'])
            t += "NESTOUT '{0}' '{1}' OUT {2} {3} MIN\n$\n".format(
                    'nest{0}'.format(num_nest+2), 'bounds_nest{0}.dat'.format(num_nest+2), 
                    t0_iso, dt_comp)
        
        # output
        t += "BLOCK 'COMPGRID' NOHEAD '{0}' LAY 3 HSIGN TM02 DIR TPS DSPR OUT {1} {2} MIN\n$\n".format(
            mm.output_fn, t0_iso, dt_comp)

        # output points
        if not x_out or not y_out:
            pass
        else:
            t += "POINTS 'outpts' FILE 'points_out.dat'\n"
            t += "TABLE 'outpts' NOHEAD {0} DEP HS HSWELL DIR RTP TM02 DSPR WIND WATLEV  OUT {1} {2} MIN\n$\n".format("'table_outpts_nest_{0}.dat'".format(num_nest+1),t0_iso, dt_comp)


        # compute
        t += 'TEST  1,0\n'
        t += 'COMPUTE \n'
        t += 'STOP\n$\n'

        # write file:
        with open(p_file, 'w') as f:
            f.write(t)

    def build_case(self, case_id, waves_event, storm_track=None,
                   make_waves=True, make_winds=True, waves_bnd=['N', 'E', 'W', 'S']):
        '''
        Build SWAN NONSTAT case input files for given wave dataset

        case_id  - SWAN case index (int)

        waves_event - waves event time series (pandas.Dataframe)
        also contains level, tide and wind (not storm track) variables
        [n x 8] (hs, per, dir, spr, U10, V10, level, tide)

        storm_track - None / storm track time series (pandas.Dataframe)
        storm_track generated winds have priority over waves_event winds
        [n x 6] (move, vf, lon, lat, pn, p0)
        '''

        # SWAN case path
        p_case = op.join(self.proj.p_cases, case_id)

        # make execution dir
        if not op.isdir(p_case): os.makedirs(p_case)

        # make depth file for main mesh
        self.proj.mesh_main.export_dat(p_case)
        
        # make depth file for nested meshes (optional)
        if self.proj.nested:
            for i in range(self.proj.num_nest):
                self.proj.meshes_nested[i].export_dat(p_case)
        # TODO: depth files for nested meshes
        
#        self.make_depth(op.join(p_case, 'depth.dat'), self.proj.depth)
#        if self.proj.nested == True:
#            for i in np.arange(0,self.proj.num_nest):
#                self.make_depth(op.join(p_case, 'depth_nest_{0}.dat'.format(i+1)), 
#                                self.proj.depth_nested[i])

        # make output points file
        self.make_out_points(op.join(p_case, 'points_out.dat'))

        # parse pandas time index to swan iso format
        swan_iso_fmt = '%Y%m%d.%H%M'
        time_swan = pd.to_datetime(waves_event.index).strftime(swan_iso_fmt).values[:]

        # make wave files
        if make_waves:
            self.make_wave_files(p_case, waves_event, time_swan, waves_bnd)

        # make wind files
        if isinstance(storm_track, pd.DataFrame):
            # vortex model from storm tracks
            self.make_vortex_files(p_case, case_id, storm_track, 
                                   nested=False, num_nest=None)
            if self.proj.nested == True:
                for i in np.arange(0,self.proj.num_nest):
                    self.make_vortex_files(p_case, case_id, storm_track, 
                                           nested=True, num_nest=i)
        else:
            # meshgrid wind
            self.make_wind_files(p_case, waves_event)


        # make water level files
        self.make_level_files(p_case, waves_event)

        # make input.swn file
        self.make_input(
            op.join(p_case, 'input.swn'), case_id, time_swan,
            make_waves = make_waves, make_winds = make_winds,
        )
        
        # optional make input_nested.swn file(s)
        if self.proj.nested:
            for i in range(self.proj.num_nest):
                save = op.join(p_case, 'input_nest{0}.swn'.format(i+1))
                self.make_input_nested(save, case_id, time_swan, i)


    def outmat2xr(self, p_mat):

        # matlab dictionary
        dmat = loadmat(p_mat)

        # get dates from one key
        hsfs = sorted([x for x in dmat.keys() if 'Hsig' in x])
        dates_str = ['_'.join(x.split('_')[1:]) for x in hsfs]
        dates = [datetime.strptime(s,'%Y%m%d_%H%M%S') for s in dates_str]

        # read times
        l_times = []
        for ds in dates_str:
            xds_t = xr.Dataset(
               {
                   'Hsig':   (('X','Y',), dmat['Hsig_{0}'.format(ds)].T,   {'units':'m'}),
                   'Tm02':   (('X','Y',), dmat['Tm02_{0}'.format(ds)].T,   {'units':'s'}),
                   'Dir':    (('X','Y',), dmat['Dir_{0}'.format(ds)].T,    {'units':'º'}),
                   'Dspr':   (('X','Y',), dmat['Dspr_{0}'.format(ds)].T,   {'units':'º'}),
                   'TPsmoo': (('X','Y',), dmat['TPsmoo_{0}'.format(ds)].T, {'units':'s'}),
               }
            )
            l_times.append(xds_t)

        # join at times dim
        xds_out = xr.concat(l_times, dim='time')
        xds_out = xds_out.assign_coords(time=dates)

        return xds_out

    def output_case(self, p_case, mesh):
        'read .mat output file from non-stationary and returns xarray.Dataset'

        # extract output from selected mesh
        p_mat = op.join(p_case, mesh.output_fn)
        xds_out = self.outmat2xr(p_mat)

        # set X and Y values
        X, Y = mesh.get_XY()
        xds_out = xds_out.assign_coords(X=X)
        xds_out = xds_out.assign_coords(Y=Y)

        # rename to longitude latitude in spherical coords cases
        coords_spherical = self.proj.params['coords_spherical']
        if coords_spherical != None:
            xds_out = xds_out.rename({'X':'lon', 'Y':'lat'})

        return xds_out

    def get_t0_dt(self, p_input):
        'gets output points time_ini and delta_time (min) from SWAN input.swn file'

        # read input.swn and file data
        with open(p_input, 'r') as fR:
            ls = fR.readlines()

        lx = [x for x in ls if x.startswith('TABLE')][0].split(' ')
        t0_str = lx[-3]  # start date
        dt_min = lx[-2]  # dt (minutes)

        swan_iso_fmt = '%Y%m%d.%H%M'
        t0 = datetime.strptime(t0_str, swan_iso_fmt)

        return t0, dt_min

    def output_points(self, p_case):
        'read table_outpts.dat output file and returns xarray.Dataset'

        p_dat = op.join(p_case, 'table_outpts.dat')

        # variable names
        names = ['DEP', 'HS', 'HSWELL', 'DIR', 'RTP', 'TM02', 'DSPR', 'WIND',
                 'WATLEV', 'OUT' ]

        x_out = self.proj.x_out
        y_out = self.proj.y_out

        # points are mixed at output file
        np_pts = np.genfromtxt(p_dat)
        n_rows = np_pts.shape[0]

        # number of points
        n_pts = len(x_out)

        l_xds_pts = []
        for i in range(n_pts):
            ix_p = np.arange(i, n_rows, n_pts)

            np_pti = np_pts[ix_p, :]
            xds_pti = xr.Dataset({}) #, coords='time')
            for c, n in enumerate(names):
                xds_pti[n] = (('time'), np_pti[:,c])

            l_xds_pts.append(xds_pti)

        xds_out = xr.concat(l_xds_pts, dim='point')

        # add point x and y
        xds_out['x_point'] = (('point'), x_out)
        xds_out['y_point'] = (('point'), y_out)

        # add times dim values
        t0, dt_min = self.get_t0_dt(op.join(p_case, 'input.swn'))
        time_out = pd.date_range(t0, periods=len(xds_out.time), freq='{0}min'.format(dt_min))
        xds_out = xds_out.assign_coords(time=time_out)

        return xds_out

