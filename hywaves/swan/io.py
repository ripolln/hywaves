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

# output variables metadatai 'CODE':('name', 'units', 'matname')
# TODO: completar
meta_out = {
    'HSIGN':('Hsig', 'm', 'Hsig'),
    'DIR':('Dir', 'º', 'Dir'),
    'PDIR':('PkDir', 'º', 'PkDir'),
    'TM02':('Tm02', 's', 'Tm02'),
    'TPS':('Tp', 's', 'TPsmoo'),
    'RTP':('RTpeak', '?', 'RTpeak'),
    'FSPR':('FSpr', '?', 'FSpr'),
    'DSPR':('Dspr', 'º', 'Dspr'),
    'DEPTH':('Depth', 'm', 'Depth'),
    'WATLEV':('WaterLevel', 'm', 'Watlev'),
    'WIND':('Wind', 'm/s', 'Windv'), # Windv_u, Windv_v
    'OUT':('OUT', '-', 'OUT'),
}


# input.swn TEMPLATES - COMMON

def swn_coordinates(proj):
    'COORDINATES .swn block'

    coords_mode = proj.params['coords_mode']
    coords_projection = proj.params['coords_projection']

    proj_str = ''
    if coords_projection: proj_str = '{0}'.format(coords_projection)

    t = ''
    if coords_mode:
        t = 'COORDINATES {0} {1}\n'.format(coords_mode, proj_str)

    return t

def swn_set(proj):
    'SET .swn block'

    set_level = proj.params['set_level']
    set_maxerr = proj.params['set_maxerr']
    set_cdcap = proj.params['set_cdcap']
    set_convention = proj.params['set_convention']

    level_str = 'level={0}'.format(set_level)
    cdcap_str = ''
    if set_cdcap: cdcap_str = 'cdcap={0}'.format(set_cdcap)
    maxerr_str = ''
    if set_maxerr: maxerr_str = 'maxerr={0}'.format(set_maxerr)
    conv_str = ''
    if set_convention: conv_str = '{0}'.format(set_convention)

    return 'SET {0} {1} {2} {3}\n$\n'.format(
        level_str, cdcap_str, maxerr_str, conv_str,
    )

def swn_computational(proj, mesh):
    'COMPUTATIONAL GRID .swn block'
    # only regular grid and full circle spectral directions !!!

    cgrid_mdc = proj.params['cgrid_mdc']
    cgrid_flow = proj.params['cgrid_flow']
    cgrid_fhigh = proj.params['cgrid_fhigh']

    # TODO: parametro [msc] se deja ausente?
    return 'CGRID REGULAR {0} {1} {2} {3} {4} {5} {6} CIRCLE {7} {8} {9} \n$\n'.format(
        mesh.cg['xpc'], mesh.cg['ypc'], mesh.cg['alpc'],
        mesh.cg['xlenc'], mesh.cg['ylenc'], mesh.cg['mxc']-1, mesh.cg['myc']-1,
        cgrid_mdc, cgrid_flow, cgrid_fhigh,
    )

def swn_bathymetry(mesh):
    'BATHYMETRY GRID .swn block'

    mxc = mesh.dg['mxc']
    myc = mesh.dg['myc']

    # fix nested / main different behaviour
    if mesh.is_nested:
        mxc = mxc - 1
        myc = myc - 1

    t = ''
    t += 'INPGRID BOTTOM REGULAR {0} {1} {2} {3} {4} {5} {6}\n'.format(
        mesh.dg['xpc'], mesh.dg['ypc'], mesh.dg['alpc'],
        mxc, myc, mesh.dg['dxinp'], mesh.dg['dyinp'])

    t += "READINP BOTTOM 1 '{0}' {1} 0 FREE\n$\n".format(
        mesh.fn_depth, mesh.dg_idla)

    return t

def swn_physics(proj):
    'PHYSICS .swn block'

    list_physics = proj.params['physics']
    t = ''
    for l in list_physics: t += '{0}\n'.format(l)
    t += '$\n'
    return t

def swn_numerics(proj):
    'NUMERICS .swn block'

    list_numerics = proj.params['numerics']
    t = ''
    for l in list_numerics: t += '{0}\n'.format(l)
    t += '$\n'
    return t

def swn_bound_waves_nested(proj, boundn_file):
    'boundary waves (NESTED) .swn block'

    boundn_mode = proj.params['boundn_mode']

    t = ''
    t += "BOUN NEST '{0}' {1}\n".format(boundn_file, boundn_mode)
    t += '$\n'

    return t

def swn_nestout(proj, t0_iso=None, compute_deltc=None):
    'output for nested boundary waves .swn block'

    mesh_main = proj.mesh_main

    t = ''
    for mesh_n in proj.mesh_nested_list:

        t += "NGRID '{0}' {1} {2} {3} {4} {5} {6} {7}\n".format(
            mesh_n.ID, mesh_n.cg['xpc'], mesh_n.cg['ypc'], mesh_n.cg['alpc'],
            mesh_n.cg['xlenc'], mesh_n.cg['ylenc'],
            np.int32(mesh_n.cg['xlenc'] / mesh_main.cg['dxinp']),
            np.int32(mesh_n.cg['ylenc'] / mesh_main.cg['dyinp'])
        )

        # prepare nonstat times str
        nonstat_str = ''
        if t0_iso != None:
            nonstat_str = 'OUT {0} {1}'.format(t0_iso, compute_deltc)

        t += "NESTOUT '{0}' '{1}' {2} \n$\n".format(mesh_n.ID, mesh_n.fn_boundn, nonstat_str)

    return t

# input.swn TEMPLATES - STAT

def swn_bound_waves_stat(proj, ws, waves_bnd):
    'boundary waves (STAT) .swn block'

    boundw_jonswap = proj.params['boundw_jonswap']
    boundw_period = proj.params['boundw_period']

    t = ''
    t += 'BOUND SHAPespec JONswap {0} {1} DSPR DEGR\n'.format(
        boundw_jonswap, boundw_period)
    for ic in waves_bnd:
        t += "BOUN SIDE {0} CONstant PAR {1:.3f} {2:.3f} {3:.3f} {4:.3f}\n".format(
            ic, ws.hs, ws.per, ws.dir, ws.spr)
    t += '$\n'

    return t

# input.swn TEMPLATES - NONSTAT

def swn_bound_waves_nonstat(proj, waves_bnd):
    'boundary waves (NONSTAT) .swn block'

    boundw_jonswap = proj.params['boundw_jonswap']
    boundw_period = proj.params['boundw_period']

    t = ''
    t += 'BOUND SHAPespec JONswap {0} {1} DSPR DEGR\n'.format(
        boundw_jonswap, boundw_period)
    for ic in waves_bnd:
        t += "BOUN SIDE {0} CONstant FILE 'series_waves_{0}.dat'\n".format(ic)
    t += '$\n'

    return t

def swn_inp_levels_nonstat(proj, mesh, t0_iso, t1_iso):
    'input level files (NONSTAT) .swn block'

    level_deltinp = proj.params['level_deltinp']

    level_fn = 'series_level_{0}.dat'.format(mesh.ID)
    level_idla = 3 # TODO: comprobar archivos level se generan acorde

    t = ''
    t += 'INPGRID  WLEV  REGULAR {0} {1} {2} {3} {4} {5} {6} NONSTAT {7} {8} {9}\n'.format(
        mesh.cg['xpc'], mesh.cg['ypc'], mesh.cg['alpc'],
        mesh.cg['mxc']-1, mesh.cg['myc']-1, mesh.cg['dxinp'], mesh.cg['dyinp'],
        t0_iso, level_deltinp, t1_iso)
    t += "READINP  WLEV 1. SERIES '{0}' {1} 0 FREE\n$\n".format(
        level_fn, level_idla)

    return t

def swn_inp_winds_nonstat(proj, mesh, t0_iso, t1_iso, wind_deltinp):
    'input level files (NONSTAT) .swn block'

    wind_fn = 'series_wind_{0}.dat'.format(mesh.ID)
    wind_idla = 3

    t = ''
    t += 'INPGRID  WIND  REGULAR {0} {1} {2} {3} {4} {5} {6} NONSTAT {7} {8} {9}\n'.format(
        mesh.cg['xpc'], mesh.cg['ypc'], mesh.cg['alpc'],
        mesh.cg['mxc']-1, mesh.cg['myc']-1, mesh.cg['dxinp'], mesh.cg['dyinp'],
        t0_iso, wind_deltinp, t1_iso)
    t += "READINP  WIND 1. SERIES '{0}' {1} 0 FREE\n$\n".format(
        wind_fn, wind_idla)

    return t


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

    def output_case(self, p_case, mesh):
        'read .mat output file from non-stationary and returns xarray.Dataset'

        # extract output from selected mesh
        p_mat = op.join(p_case, mesh.fn_output)
        xds_out = self.outmat2xr(p_mat)

        # set X and Y values
        X, Y = mesh.get_XY()
        xds_out = xds_out.assign_coords(X=X)
        xds_out = xds_out.assign_coords(Y=Y)

        # rename to longitude latitude in spherical coords cases
        coords_mode = self.proj.params['coords_mode']
        if coords_mode == 'SPHERICAL':
            xds_out = xds_out.rename({'X':'lon', 'Y':'lat'})

        return xds_out


class SwanIO_STAT(SwanIO):
    'SWAN numerical model input/output - STATIONARY cases'

    def make_input(self, p_file, id_run,
                   mesh, ws,
                   waves_bnd=['N', 'E', 'W', 'S'],
                   ttl_run = ''):
        '''
        Writes input.swn file from waves sea state for stationary execution

        p_file      - input.swn file path
        mesh        - SwanMesh instance
        ws          - wave sea state (hs, per, dr, spr)
        waves_bnd   - wave sea state active boundaries

        ttl_run     - project title (optional)

        more info: http://swanmodel.sourceforge.net/online_doc/swanuse/node23.html
        '''

        # -- PROJECT --
        t = "PROJ '{0}' '{1}' '{2}'\n$\n".format(self.proj.name, id_run, ttl_run)

        # -- MODE STATIONARY --
        t += 'MODE STAT\n'

        # -- COORDINATES --
        t += swn_coordinates(self.proj)

        # -- SET -- 
        t += swn_set(self.proj)

        # -- COMPUTATIONAL GRID --
        t += swn_computational(self.proj, mesh)

        # -- BATHYMETRY --
        t += swn_bathymetry(mesh)

        # -- SWAN STATIONARY -- INPUT WAVES --

        # MAIN mesh - boundary waves
        if not mesh.is_nested:
            t += swn_bound_waves_stat(self.proj, ws, waves_bnd)

        # NESTED mesh - nested waves
        else:
            t += swn_bound_waves_nested(self.proj, mesh.fn_boundn)

        # -- PHYSICS --
        t += swn_physics(self.proj)

        # -- NUMERICS --
        t += swn_numerics(self.proj)

        # -- OUTPUT: NESTED MESHES  -- 
        if not mesh.is_nested:
            t += swn_nestout(self.proj)

        # -- OUTPUT: BLOCK  -- 
        t += "BLOCK 'COMPGRID' NOHEAD '{0}' LAY 3 HSIGN TM02 DIR TPS DSPR\n$\n".format(
            mesh.fn_output,
        )

        # -- COMPUTE --
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

    def build_case(self, case_id, waves_ss, waves_bnd=['N', 'E', 'W', 'S']):
        '''
        Build SWAN STAT case input files for given wave sea state (hs, per, dir, spr)

        ix_case  - SWAN case index (int)
        waves_ss - wave sea state (hs, per, dir, spr)
        bnd      - wave sea state active boundaries
        '''

        # SWAN case path
        p_case = op.join(self.proj.p_cases, case_id)
        if not op.isdir(p_case): os.makedirs(p_case)

        # MAIN mesh
        self.proj.mesh_main.export_depth(p_case)  # export main depth file
        p_swn = op.join(p_case, self.proj.mesh_main.fn_input)

        # make input.swn file
        self.make_input(
            p_swn, case_id,
            self.proj.mesh_main,
            waves_ss,
            waves_bnd = waves_bnd,
        )

        # NESTED meshes 
        for mesh_nested in self.proj.mesh_nested_list:

            mesh_nested.export_depth(p_case)  # export nested depth file
            p_swn = op.join(p_case, mesh_nested.fn_input)

            # make input_nestX.swn file
            self.make_input(
                p_swn, case_id,
                mesh_nested,
                waves_ss,
            )

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


class SwanIO_NONSTAT(SwanIO):
    'SWAN numerical model input/output - NON STATIONARY cases'

    def make_out_points(self, p_file):
        'Generates desired output-points coordinates file'

        # define and save output points
        x_out = self.proj.params['output_points_x']
        y_out = self.proj.params['output_points_y']

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

    def make_wind_files(self, p_case, waves_event, mesh):
        '''
        Generate event wind mesh files (swan compatible)

        uses wave_event U10 and V10 values at the entire SWAN comp. grid
        '''

        # wind variables
        u10 = waves_event.U10.values[:]
        v10 = waves_event.V10.values[:]

        # each time needs 2D (mesh) wind files (U,V) 
        mxc = mesh.cg['mxc']  # number mesh x
        myc = mesh.cg['myc']  # number mesh y
        code = 'wind_{0}'.format(mesh.ID)

        txt = ''
        for c, (u, v) in enumerate(zip(u10,v10)):

            # single point wind -> entire SWAN comp.grid wind
            aux = np.ones((mxc, myc))

            # TODO: wind has to be rotated if alpc != 0

            # csv file 
            u_2d = aux * u
            v_2d = aux * v
            u_v_stack = np.vstack((u_2d, v_2d))
            save = op.join(p_case, '{0}_{1:06}.dat'.format(code, c))
            np.savetxt(save, u_v_stack, fmt='%.2f')

            # wind list file
            txt += '{0}_{1:06}.dat\n'.format(code, c)

        # winds file path
        save = op.join(p_case, 'series_{0}.dat'.format(code))
        with open(save, 'w') as f:
            f.write(txt)

    def make_vortex_files(self, p_case, case_id, mesh,
                          storm_track):
        '''
        Generate event wind mesh files (swan compatible)

        uses wave_event storm path data over SWAN computational grid
        needs SPHERICAL COORDINATES

        mesh      - mesh (main or nested)
        '''

        # TODO: llevarlo a hywaves/swan/storms 

        code = 'wind_{0}'.format(mesh.ID)

        # parameters
        RE = 6378.135 * 1000            # Earth radius [m]
        beta = 0.9                      # conversion factor of wind speed
        rho_air = 1.15                  # air density
        w = 2 * np.pi / 86184.2         # Earth's rotation velocity (rad/s)
        pifac = np.arccos(-1) / 180     # pi/180
        one2ten = 0.8928                # conversion from 1-min to 10 min

        # wind variables
        # TODO: vf?
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

        # select main mesh or nested mesh
        mm = mesh

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

            # csv file 
            save = op.join(p_case, '{0}_{1:06}.dat'.format(code, c))
            np.savetxt(save, u_v_stack, fmt='%.2f')

            # wind list file
            txt += '{0}_{1:06}.dat\n'.format(code, c)

            # hold wind data (m/s)
            hld_W[:,:,c] = W
            hld_D[:,:,c] =  270 - np.rad2deg(ang)  # direction (º clock. rel. north)

        # winds file path
        save = op.join(p_case, 'series_{0}.dat'.format(code))
        with open(save, 'w') as f:
            f.write(txt)

        # aux. save vortex wind fields
        if south_hemisphere == True:
            cg_lon = np.linspace(lon00, lon01, mxc)
            cg_lat = np.linspace(lat00, lat01, myc)

        # aux. save vortex wind fields
        p_vortex = op.join(p_case, 'vortex_{0}.nc'.format(code))
        xds_vortex = xr.Dataset(
            {
                'W': (('lat','lon','time'), hld_W, {'units':'m/s'}),
                'Dir': (('lat','lon','time'), hld_D, {'units':'º'})
            },
            coords={
                'lat' : cg_lat,
                'lon' : cg_lon,
                'time' : times,
            }
        )
        xds_vortex.attrs['xlabel'] = 'Longitude (º)'
        xds_vortex.attrs['ylabel'] = 'Latitude (º)'

        xds_vortex.to_netcdf(p_vortex)

    def make_level_files(self, p_case, wave_event, mesh):
        'Generate event level mesh files (swan compatible)'

        # parse pandas time index to swan iso format
        swan_iso_fmt = '%Y%m%d.%H%M'
        time = pd.to_datetime(wave_event.index).strftime(swan_iso_fmt).values[:]

        # level variables
        zeta = wave_event.level.values[:]
        tide = wave_event.tide.values[:]

        # each time needs 2D (mesh) level 
        mxc = mesh.cg['mxc']  # number mesh x
        myc = mesh.cg['myc']  # number mesh y
        code = 'level_{0}'.format(mesh.ID)

        txt = ''
        for c, (z, t) in enumerate(zip(zeta, tide)):

            # single point level -> entire SWAN comp.grid level
            aux = np.ones((mxc, myc)).T

            # csv file 
            l = z + t  # total level
            l_2d = aux * l
            save = op.join(p_case, '{0}_{1:06}.dat'.format(code, c))
            np.savetxt(save, l_2d, fmt='%.2f')

            # level list file
            txt += '{0}_{1:06}.dat\n'.format(code, c)

        # waves file path
        save = op.join(p_case, 'series_{0}.dat'.format(code))
        with open(save, 'w') as f:
            f.write(txt)

    def make_input(self, p_file, id_run,
                   mesh,
                   time, compute_deltc, wind_deltinp=None,
                   ttl_run='',
                   make_waves=True, make_winds=True, make_levels=True,
                   waves_bnd=['N', 'E', 'W', 'S']):
        '''
        Writes input.swn file from waves event for non-stationary execution

        p_file     - input.swn file path
        time       - event time at swan iso format
        compute_deltc - computational delta time (swan project parameter)

        ttl_run    - execution title that will appear in the output

        make_waves - activates waves input files generation (at waves_bnd)
        make_winds - activates wind input files generation
        make_levels - activates level input files generation

        more info: http://swanmodel.sourceforge.net/online_doc/swanuse/node23.html
        '''

        # -- PROJECT --
        t = "PROJ '{0}' '{1}' '{2}'\n$\n".format(self.proj.name, id_run, ttl_run)

        # -- MODE NONSTATIONARY --
        t += 'MODE NONSTAT\n'

        # -- COORDINATES --
        t += swn_coordinates(self.proj)

        # -- SET -- 
        t += swn_set(self.proj)

        # -- COMPUTATIONAL GRID --
        t += swn_computational(self.proj, mesh)

        # -- BATHYMETRY --
        t += swn_bathymetry(mesh)

        # -- SWAN NON STATIONARY -- INPUT GRIDS --
        t0_iso = time[0]   # initial time (SWAN ISOFORMAT)
        t1_iso = time[-1]  # end time (SWAN ISOFORMAT)

        # level series files
        if make_levels:
            t += swn_inp_levels_nonstat(self.proj, mesh, t0_iso, t1_iso)

        # wind series files
        if make_winds:
            t += swn_inp_winds_nonstat(self.proj, mesh, t0_iso, t1_iso,
                                       wind_deltinp)

        # -- BOUNDARY WAVES CONDITIONS --

        # MAIN mesh - boundary waves
        if not mesh.is_nested:
            if make_waves:
                t += swn_bound_waves_nonstat(self.proj, waves_bnd)

        # NESTED mesh - nested waves
        else:
            t += swn_bound_waves_nested(self.proj, mesh.fn_boundn)

        # -- PHYSICS --
        t += swn_physics(self.proj)

        # -- NUMERICS --
        t += swn_numerics(self.proj)

        # -- OUTPUT: NESTED MESHES  -- 
        if not mesh.is_nested:
            t += swn_nestout(self.proj, t0_iso=t0_iso, compute_deltc=compute_deltc)

        # output variables
        out_vars = ' '.join(self.proj.params['output_variables'])

        # -- OUTPUT: BLOCK  -- 
        dt_out = self.proj.params['output_deltt']
        t += "BLOCK 'COMPGRID' NOHEAD '{0}' LAY 3 {1} {2} {3}\n$\n".format(
            mesh.fn_output, out_vars, t0_iso, dt_out)

        # -- OUTPUT: POINTS  -- 
        x_out = self.proj.params['output_points_x']
        y_out = self.proj.params['output_points_y']

        if not x_out or not y_out:
            pass
        else:
            t += "POINTS 'outpts' FILE 'points_out.dat'\n"
            t += "TABLE 'outpts' NOHEAD '{0}' {1} {2} {3}\n$\n".format(
                mesh.fn_output_points, out_vars, t0_iso, dt_out)

        # -- COMPUTE --
        t += 'TEST  1,0\n'
        t += 'COMPUTE NONSTAT {0} {1} {2}\n'.format(t0_iso, compute_deltc, t1_iso)
        t += 'STOP\n$\n'

        # write file:
        with open(p_file, 'w') as f:
            f.write(t)

    def build_case(self, case_id, waves_event, storm_track=None,
                   make_waves=True, make_winds=True, make_levels=True,
                   waves_bnd=['N', 'E', 'W', 'S']):
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

        # parse pandas time index to swan iso format
        swan_iso_fmt = '%Y%m%d.%H%M'
        time_swan = pd.to_datetime(waves_event.index).strftime(swan_iso_fmt).values[:]

        # project computational and winds_input delta time
        compute_deltc = self.proj.params['compute_deltc']
        wind_deltinp = proj.params['wind_deltinp']

        # SWAN case path
        p_case = op.join(self.proj.p_cases, case_id)
        if not op.isdir(p_case): os.makedirs(p_case)

        # MAIN mesh
        self.proj.mesh_main.export_depth(p_case)  # export depth file

        # make water level files
        if make_levels: self.make_level_files(p_case, waves_event, self.proj.mesh_main)

        # make wave files
        if make_waves: self.make_wave_files(p_case, waves_event, time_swan, waves_bnd)

        # make wind files
        if make_winds:

            # vortex model from storm tracks  //  meshgrind wind
            if isinstance(storm_track, pd.DataFrame):
                self.make_vortex_files(p_case, case_id, self.proj.mesh_main, storm_track)

                # optional: override computational/winds dt with storm track attribute 
                if 'override_dtcomp' in storm_track.attrs:
                    compute_deltc = storm_track.attrs['override_dtcomp']
                    wind_deltinp = storm_track.attrs['override_dtcomp']
                    print('CASE {0} - compute_deltc, wind_deltinp override with storm track: {1}'.format(
                        case_id, compute_deltc))
            else:
                self.make_wind_files(p_case, waves_event, self.proj.mesh_main)

        # make output points file
        self.make_out_points(op.join(p_case, 'points_out.dat'))

        # make input.swn file
        p_swn = op.join(p_case, self.proj.mesh_main.fn_input)

        self.make_input(
            p_swn, case_id,
            self.proj.mesh_main,
            time_swan,
            compute_deltc, wind_deltinp,
            make_waves = make_waves,
            make_winds = make_winds,
            make_levels = make_levels,
            waves_bnd = waves_bnd
        )

        # NESTED mesh depth and input (wind, level) files
        for mesh_n in self.proj.mesh_nested_list:

            mesh_n.export_depth(p_case)  # export nested depth file
            p_swn = op.join(p_case, mesh_n.fn_input)

            if make_levels:
                self.make_level_files(p_case, waves_event, mesh_n)

            if make_winds:

                # vortex model from storm tracks  //  meshgrind wind
                if isinstance(storm_track, pd.DataFrame):
                    self.make_vortex_files(p_case, case_id, mesh_n, storm_track)
                else:
                    self.make_wind_files(p_case, waves_event, mesh_n)

            # make input_nestX.swn file
            self.make_input(
                p_swn, case_id,
                mesh_n,
                time_swan,
                compute_deltc, wind_deltinp,
                make_waves = make_waves,
                make_winds = make_winds,
                make_levels = make_levels,
            )

    def outmat2xr(self, p_mat):

        # matlab dictionary
        dmat = loadmat(p_mat)

        # get dates from one key
        hsfs = sorted([x for x in dmat.keys() if 'Hsig' in x])
        dates_str = ['_'.join(x.split('_')[1:]) for x in hsfs]
        dates = [datetime.strptime(s,'%Y%m%d_%H%M%S') for s in dates_str]

        # variables to extract
        # TODO detectarlas del .mat automaticamente?
        names = self.proj.params['output_variables']
        not_proc = ['WIND', 'OUT']  # filter variables

        # read times
        l_times = []
        for ds in dates_str:

            xds_t = xr.Dataset()
            for vn in names:
                if vn in meta_out.keys() and vn not in not_proc:
                    vn_ni = meta_out[vn][0]
                    vn_un = meta_out[vn][1]
                    mat_code = meta_out[vn][2]
                    xds_t[vn_ni] = (('Y','X',), dmat['{0}_{1}'.format(mat_code, ds)], {'units':vn_un})

            l_times.append(xds_t)

        # join at times dim
        xds_out = xr.concat(l_times, dim='time')
        xds_out = xds_out.assign_coords(time=dates)

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

    def output_points(self, p_case, mesh):
        'read table_outpts_meshID.dat output file and returns xarray.Dataset'

        # extract output from selected mesh
        p_dat = op.join(p_case, mesh.fn_output_points)

        # variable names
        names = self.proj.params['output_variables']

        x_out = self.proj.params['output_points_x']
        y_out = self.proj.params['output_points_y']

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

                # search metadata for good description and override n
                if n in meta_out.keys():
                    n='{0} ({1})'.format(meta_out[n][0], meta_out[n][1])

                xds_pti[n] = (('time'), np_pti[:,c])

            l_xds_pts.append(xds_pti)

        xds_out = xr.concat(l_xds_pts, dim='point')

        # add point x and y
        xds_out['x_point'] = (('point'), x_out)
        xds_out['y_point'] = (('point'), y_out)

        # mesh ID
        xds_out.attrs['mesh_ID'] = mesh.ID

        # add times dim values
        p_swn = op.join(p_case, self.proj.mesh_main.fn_input)
        t0, dt_min = self.get_t0_dt(p_swn)
        time_out = pd.date_range(t0, periods=len(xds_out.time), freq='{0}min'.format(dt_min))
        xds_out = xds_out.assign_coords(time=time_out)

        return xds_out

