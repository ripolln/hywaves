#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import subprocess as sp
import sys

import numpy as np
import xarray as xr

# SWAN STAT LIBRARY
from .io import SwanIO_STAT, SwanIO_NONSTAT


# grid description template
d_grid_template = {
    'xpc': None,      # x origin
    'ypc': None,      # y origin
    'alpc': None,     # x-axis direction 
    'xlenc': None,    # grid length in x
    'ylenc': None,    # grid length in y
    'mxc': None,      # number mesh x
    'myc': None,      # number mesh y
    'dxinp': None,    # size mesh x
    'dyinp': None,    # size mesh y
}

# swan input parameters template
d_params_template = {
    'sea_level': None,
    'jonswap_gamma': None,
    'coords_spherical': None,  # None, 'GCM', 'CCM'
    'cdcap': None,
    'maxerr': None,            # None, 1, 2, 3
    'waves_period': None,      # 'PEAK', 'MEAN'
    'nested_bounds': None,     # 'CLOSED', 'OPEN'
}


class SwanMesh(object):
    'SWAN numerical model mesh'

    def __init__(self):

        self.depth = None   # bathymetry depth value (2D numpy.array)
        self.depth_fn = ''  # filename used in SWAN execution
        self.output_fn = ''  # output .mat file for mesh comp. grid

        # grid parameters
        self.cg = d_grid_template.copy()  # computational grid
        self.dg = d_grid_template.copy()  # depth grid
        self.dg_idla = 3  # 1/3 input swan parameter, handles read order at readinp

    def export_dat(self, p_case):
        'exports depth values to .dat file'
        
        # TODO: it needs to be able to write also "self.depth_nested" in 

        p_export = op.join(p_case, self.depth_fn)
        np.savetxt(p_export, self.depth, fmt='%.2f')

    def get_XY(self):
        'returns mesh X, Y arrays from computational grid'

        # computational grid
        cg = self.cg

        x0 = cg['xpc']
        x1 = cg['xlenc'] + cg['xpc'] - cg['dxinp']
        xN = cg['mxc']
        X = np.linspace(x0, x1, xN)

        y0 = cg['ypc']
        y1 = cg['ylenc'] + cg['ypc'] - cg['dyinp']
        yN = cg['myc']
        Y = np.linspace(y0, y1, yN)

        return X, Y


class SwanProject(object):
    'SWAN numerical model project parameters, grids and information'

    def __init__(self, p_proj, n_proj):
        '''
        SWAN project information will be stored here

        http://swanmodel.sourceforge.net/online_doc/swanuse/node25.html
        '''

        self.p_main = op.join(p_proj, n_proj)    # project path
        self.name = n_proj                       # project name
        self.storm = None                        # storm name

        # sub folders 
        self.p_cases = op.join(self.p_main, 'cases')  # project cases

        # nested parameters
        self.nested = None          # activation
        self.num_nest = None        # number of nested meshes to run
        self.lon_nested = None      # longitude of nested meshgrids
        self.lat_nested = None      # latitude of nested meshgrids

        # SWAN mesh: main 
        self.mesh_main = SwanMesh()
        self.mesh_main.depth_fn = 'depth_main.dat'
        self.mesh_main.output_fn = 'output_main.mat'


        # SWAN mesh: nest1
        self.mesh_nest1 = SwanMesh()
        self.mesh_nest1.depth_fn = 'depth_nest1.dat'
        self.mesh_nest1.output_fn = 'output_nest1.mat'
        self.run_nest1 = False

        # SWAN mesh: nest2
        self.mesh_nest2 = SwanMesh()
        self.mesh_nest2.depth_fn = 'depth_nest2.dat'
        self.mesh_nest2.output_fn = 'output_nest2.mat'
        self.run_nest2 = False

        # SWAN mesh: nest3
        self.mesh_nest3 = SwanMesh()
        self.mesh_nest3.depth_fn = 'depth_nest3.dat'
        self.mesh_nest3.output_fn = 'output_nest3.mat'
        self.run_nest3 = False
        
        self.meshes_nested = [self.mesh_nest1, self.mesh_nest2, self.mesh_nest3]
        # TODO: generate "self.mesh_nest_i" for range("num_nest")
        
#        # SWAN mesh: nested
#        if self.nested == True:
#            i = 0
#            while i < self.num_nest:
                

        # swan execution parameters
        self.params = d_params_template.copy()
        
        # set nested meshes to run
        if self.num_nest == 1:
            self.run_nest1 = True
        if self.num_nest == 2:
            self.run_nest1 = True
            self.run_nest2 = True
        if self.num_nest == 3:
            self.run_nest1 = True
            self.run_nest2 = True
            self.run_nest3 = True

        # output points
        self.x_out = []
        self.y_out = []

        # shoreline for plotting
        self.np_shore = np.array([])
        self.np_shore_nested = np.array([])
        
        # computation time step (equal to the storm dt interpolation)
        self.dt_comp = None
        self.dt_comp_nested = None
        
#        # bound limits of the general mesh vortex
#        self.vortex_max = None
#        self.vortex_min = None
#        
#        # save plots
#        self.plot_vortex = None
#        self.plot_track = None
        

class SwanWrap(object):
    'SWAN numerical model wrap for multi-case handling'

    def __init__(self, swan_proj, swan_io):
        '''
        swan_proj - SwanProject() instance, contains project parameters
        swan_io   - SwanIO_STAT / SwanIO_NONSTAT modules (auto from children)
        '''

        # set project and IO module
        self.proj = swan_proj           # swan project parameters
        self.io = swan_io(self.proj)     # swan input/output 

        # swan bin executable
        p_res = op.join(op.dirname(op.realpath(__file__)), 'resources')
        self.bin = op.abspath(op.join(p_res, 'swan_bin', 'swan_ser.exe'))

    def get_run_folders(self):
        'return sorted list of project cases folders'

        # TODO: will find previously generated cases... fix it 

        ldir = sorted(os.listdir(self.proj.p_cases))
        fp_ldir = [op.join(self.proj.p_cases, c) for c in ldir]

        return [p for p in fp_ldir if op.isdir(p)]

    def run_cases(self):
        'run all cases inside project "cases" folder'

        # TODO: improve log / check execution ending status

        # get sorted execution folders
        run_dirs = self.get_run_folders()

        for p_run in run_dirs:

            # run case main mesh
            self.run(p_run)

            # TODO: integrate "num_nest"

            # run case nested mesh (optional)
            r_ns = [
                self.proj.run_nest1,
                self.proj.run_nest2,
                self.proj.run_nest3,
            ]
            i_ns = [
                'input_nest1.swn',
                'input_nest2.swn',
                'input_nest3.swn'
            ]
            for r_n, i_n in zip(r_ns, i_ns):
                if r_n:
                    self.run(p_run, input_file=i_n)

            # log
            p = op.basename(p_run)
            print('SWAN CASE: {0} SOLVED'.format(p))

    def run(self, p_run, input_file='input.swn'):
        'Bash execution commands for launching SWAN'

        # aux. func. for launching bash command
        def bash_cmd(str_cmd, out_file=None, err_file=None):
            'Launch bash command using subprocess library'

            _stdout = None
            _stderr = None

            if out_file:
                _stdout = open(out_file, 'w')
            if err_file:
                _stderr = open(err_file, 'w')

            s = sp.Popen(str_cmd, shell=True, stdout=_stdout, stderr=_stderr)
            s.wait()

            if out_file:
                _stdout.flush()
                _stdout.close()
            if err_file:
                _stderr.flush()
                _stderr.close()

        # check if windows OS
        is_win = sys.platform.startswith('win')

        if is_win:
            # WINDOWS - use swashrun command
            cmd = 'cd {0} && swashrun input && {1} input'.format(
                p_run, self.bin)

        else:
            # LINUX/MAC - ln input file and run swan case
            cmd = 'cd {0} && ln -sf {1} INPUT && {2} INPUT'.format(
                p_run, input_file, self.bin)

        bash_cmd(cmd)

#        # run nested meshes if activated
#        if self.proj.nested == True:
#            for i in np.arange(0,self.proj.num_nest):
#                cmdn = 'cd {0} && ln -sf input_nest{1}.swn INPUT && {2} INPUT'.format(
#                    p_run, i+1, self.bin)
#                bash_cmd(cmdn)


    def extract_output(self, mesh=None):
        '''
        exctract output from all cases generated by "build_cases_stat"

        return xarray.Dataset (uses new dim "case" to join output)
        '''

        # select main or nested mesh
        if mesh == None: mesh = self.proj.mesh_main

        # get sorted execution folders
        run_dirs = self.get_run_folders()

        # exctract output case by case and concat in list
        l_out = []
        for p_run in run_dirs:

            # read output file
            xds_case_out = self.io.output_case(p_run, mesh)
            l_out.append(xds_case_out)

        # concatenate xarray datasets (new dim: case)
        xds_out = xr.concat(l_out, dim='case')

        return(xds_out)

    def extract_output_nonstat(self, output_run, case_ini=None, case_end=None, num=None):
        '''
        exctract output from all cases generated by "build_cases_nonstat"
        (it is possible to choose which cases to extract)

        return xarray.Dataset (uses new dim "case" to join output)
        '''
        
        # TODO: integrate mesh input as for "extract_output"
        # in the code below I can choose which cases to extract separately 
        
        # get sorted execution folders
        run_dirs = self.get_run_folders()
        if (case_ini != None) & (case_end != None):   
            run_dirs=run_dirs[case_ini:case_end]

        # exctract output case by case and concat in list
        l_out = []
        for p_run in run_dirs:

            # read output file
            xds_case_out = self.io.output_case_nonstat(p_run, output_run, num)
            l_out.append(xds_case_out)

        # concatenate xarray datasets (new dim: case)
        xds_out = xr.concat(l_out, dim='case')
        if (case_ini != None) & (case_end != None):   
            xds_out = xds_out.assign(case=np.arange(case_ini, case_end))

        return(xds_out)

    def extract_output_points(self, case_ini=None, case_end=None):
        '''
        extract output from points all cases table_outpts.dat
        (it is possible to choose which cases to extract)

        return xarray.Dataset (uses new dim "case" to join output)
        '''
        
        # TODO: integrate mesh, same as for "extract_output"

        # get sorted execution folders
        run_dirs = self.get_run_folders()
        if (case_ini != None) & (case_end != None):   
            run_dirs=run_dirs[case_ini:case_end]

        # exctract output case by case and concat in list
        l_out = []
        for p_run in run_dirs:

            # read output file
            xds_case_out = self.io.output_points(p_run)
            l_out.append(xds_case_out)

        # concatenate xarray datasets (new dim: case)
        xds_out = xr.concat(l_out, dim='case')
        if (case_ini != None) & (case_end != None):   
            xds_out = xds_out.assign(case=np.arange(case_ini, case_end))

        return(xds_out)


class SwanWrap_STAT(SwanWrap):
    'SWAN numerical model wrap for STATIONARY multi-case handling'

    def __init__(self, swan_proj):
        super().__init__(swan_proj, SwanIO_STAT)

    def build_cases(self, waves_dataset):
        '''
        generates all files needed for swan stationary multi-case execution

        waves_dataset - pandas.dataframe with "n" boundary conditions setup
        [n x 4] (hs, per, dir, spr)
        '''

        # make main project directory
        self.io.make_project()

        # one stat case for each wave sea state
        for ix, (_, ws) in enumerate(waves_dataset.iterrows()):

            # build stat case 
            case_id = '{0:04d}'.format(ix)
            self.io.build_case(case_id, ws)


class SwanWrap_NONSTAT(SwanWrap):
    'SWAN numerical model wrap for NON STATIONARY multi-case handling'

    def __init__(self, swan_proj):
        super().__init__(swan_proj, SwanIO_NONSTAT)

    def build_cases(self, waves_event_list, storm_track_list=None,
                    make_waves=True, make_winds=True):
        '''
        generates all files needed for swan non-stationary multi-case execution

        waves_event_list - list waves events time series (pandas.DataFrame)
        also contains level, tide and wind (not storm track) variables
        [n x 8] (hs, per, dir, spr, U10, V10, level, tide)

        storm_track_list - list of storm tracks time series (pandas.DataFrame)
        storm_track generated winds have priority over waves_event winds
        [n x 6] (move, vf, lon, lat, pn, p0)
        '''

        # check user input: no storm tracks
        if storm_track_list == None:
            storm_track_list = [None] * len(waves_event_list)

        # make main project directory
        self.io.make_project()

        # one non-stationary case for each wave time series
        for ix, (wds, sds) in enumerate(
            zip(waves_event_list, storm_track_list)):

            # build stat case 
            case_id = '{0:04d}'.format(ix)
            self.io.build_case(
                case_id, wds, storm_track=sds,
                make_waves=make_waves, make_winds=make_winds
            )

