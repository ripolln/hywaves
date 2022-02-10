#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import numpy as np
import pandas as pd
import xarray as xr
import gc
from datetime import timedelta
import datetime
import sys
from tqdm import tqdm
from wavespectra.specarray import SpecArray
from wavespectra.specdataset import SpecDataset
#from math import radians, degrees, sin, cos, atan2, pi  #asin, acos, sqrt, 

from .geo import shoot, gc_distance, GeoAzimuth
from .storms import Extract_basin_storms, get_vmean, stopmotion_trim_circle, get_category,\
d_vns_basinscenter, historic_track_preprocessing, historic_track_interpolation
from ..statistical.mda import nearest_indexes, nearest_indexes_weighted

# auxiliar function to generate custom polar coordinates to extract SWAN output

def generate_polar_coords(xlon, ylat, res, 
                          rings=28, radii_ini=5000, angle_inc=5):
    '''
    xlon,ylat  - domain dimensions [m]
    res        - spatial resolution [m]
    
    rings      - dimension of radii array
    radi_ini   - initial radii [m]
    angle_inc  - increment interval [º]
    
    returns:  (xxr,yyr) regular and (xxp,yyp) polar grid points 
    '''
    
    # A: regular cartesian grid
    xxr, yyr = np.meshgrid(xlon, ylat) 

    # B: custom polar grid calculation

    # radii array
    rc = np.zeros(rings)
    rc[0] = radii_ini
    for i in range(1, rings):       rc[i] = ((i+1)/i) ** 1.5 * rc[i-1]
    
    # radii (m), angle (º) meshgrid
    rc, thetac = np.meshgrid(rc, np.arange(0, 360, angle_inc))
    
    # polar grid
    xxp = rc * np.sin(thetac * np.pi/180)
    yyp = rc * np.cos(thetac * np.pi/180)
    
    return xxr, yyr, xxp, yyp, rc, thetac  # 2d-arrays


###############################################################################
# STOPMOTION library
# functions to preprocess and interpolate coordinates at swan computational 
# timesteps for any storm track according to the SHyTCWaves methodology units:
# 6-hour target segment preceded by 24-hour warmup segments 
###############################################################################

def pres2wind(xds_coef, st_pres, st_center, st_basin):
    '''
    xds_coef    - (xarray.Dataset) fitting coefficients for Pmin-Vmax
    st_pres     - storm pressure [mbar]
    st_center   - storm center (RSMC, ibtracs)
    st_basin    - storm basin (NA,SA,WP,EP,SP,NI,SI)
    
    returns:  empirical maximum wind speed [kt, 1-min avg]
    '''
    
    # select Pmin-Wmax coefficients
    coefs = xds_coef.sel(center=st_center, basin=st_basin
                         ).coef.values[:]
    
    # prediction
    wmax_pred = np.polyval(coefs, st_pres)
    
    return wmax_pred    # [kt]

def storm2stopmotion(df_storm):
    '''
    df_storm    - (pandas.DataFrame)  Variables from storm preprocessing: 
                  (lon,lat,p0,vmax,rmw, vmaxfill,rmwfill)
                  * _fill inform of historical gaps filled with estimates
                  * p0,vmax,rmw: mean values for 6h forward segment
                  
    Generation of stopmotion units methodology from storm track segments of 6h:
        
    A.  Warmup segment (24h):
            4 segments to define start/end coordinates, {Vmean, relative angle}
            last 4th segment, mean {Pmin, Wmax, Rmw}, endpoint {lat}
            
    B.  Target segment (6h): {dP,dV,dW,dR,dAng}
        
    * Absolute value of latitude is stored, start of target segment
    * Relative angle is referenced to the geographic north (southern 
      hemisphere is multiplied with a factor of -1)
    
    returns:  df_ (pandas.Dataframe)
    '''
    
    # no need to remove NaTs -> historic_track_preprocessing
    # no need to remove NaNs -> historic_track_preprocessing
    # no need to remove wind/rmw -> filled data at historic_track_interpolation
    
    # remove NaN
    df_ = df_storm.dropna()
    df_['time'] = df_.index.values

#    # round time to hours
#    df['time'] = pd.to_datetime(df['time'], format='%Y-%m-%d %H:%M:%S'
#                                ).dt.round('1h')
#
#    # only keep 0,6,12,18 hours    (REDUNDANT ???)
#    hr = df['time'].dt.hour.values
#    pos = np.where((hr==0) | (hr==6) | (hr==12) | (hr==18))[0]
#    df_ = df.iloc[pos]
#    df_.index = df_['time']
    
    # constant segments variables
    lon  = df_['lon'].values[:]
    lat  = df_['lat'].values[:]
    pres = df_['p0'].values[:]      # [mbar]
    wind = df_['vmax'].values[:]    # [kt]
    rmw  = df_['rmw'].values[:]     # [nmile]

    # generate stopmotion segments: 24h warmup + 6h target segments
    
    # timestep [hours]
    dt = np.diff(df_.index) / np.timedelta64(1, 'h')
    
    # warmup 4-segments (24h) variables
    lo0   = np.full(df_.shape[0], np.nan)   # warmup x-coordinate
    la0   = np.full(df_.shape[0], np.nan)   # warmup y-coordinate
    vseg  = np.full(df_.shape[0], np.nan)   # mean translational speed
    vxseg = np.full(df_.shape[0], np.nan)   # (dirx)
    vyseg = np.full(df_.shape[0], np.nan)   # (diry)
    aseg  = np.full(df_.shape[0], np.nan)   # azimuth, geographic North

    # warmup last-segment (6h) variables
    pseg  = np.full(df_.shape[0], np.nan)   # segment pressure
    wseg  = np.full(df_.shape[0], np.nan)   # segment maxwinds
    rseg  = np.full(df_.shape[0], np.nan)   # segment radii rmw
    lseg  = np.full(df_.shape[0], np.nan)   # latitude (north hemisphere)
    laseg = np.full(df_.shape[0], np.nan)   # (absolute)
    
    # target segment (6h) variables
    lo1 = np.full(df_.shape[0], np.nan)     # target origin x-coordinate
    la1 = np.full(df_.shape[0], np.nan)     # idem y-coordinate
    dv  = np.full(df_.shape[0], np.nan)     # translational speed variation
    dvx = np.full(df_.shape[0], np.nan)     # (dirx)
    dvy = np.full(df_.shape[0], np.nan)     # (diry)
    da  = np.full(df_.shape[0], np.nan)     # azimuth variation
    dp  = np.full(df_.shape[0], np.nan)     # pressure variation
    dw  = np.full(df_.shape[0], np.nan)     # maxwinds variation
    dr  = np.full(df_.shape[0], np.nan)     # radii rmw variation
    dl  = np.full(df_.shape[0], np.nan)     # latitude variation
    dla = np.full(df_.shape[0], np.nan)     # (absolute)
    lo2 = np.full(df_.shape[0], np.nan)     # target endpoint x-coordinate
    la2 = np.full(df_.shape[0], np.nan)     # idem y-coordinate
    
    # loop
    for i in np.arange(1, dt.size):
        
        # get stopmotion endpoints coordinates (24h+6h)
        if i < 4:   # < four preceding segments

            # number of "missing" preceding segments to last 24h
            n_missing = 4 - i
            
            # last available preceding segment
            lon1, lon2 = lon[1], lon[0]
            lat1, lat2 = lat[1], lat[0]
            
            # distance of last available preceding segment
            arcl_h, gamma_h = gc_distance(lat2, lon2, lat1, lon1)
            RE = 6378.135                       # earth radius [km]
            r = arcl_h * np.pi / 180.0 * RE     # distance [km]

            # shoot backwards to calculate (lo0,la0) of 24h preceding warmup 
            dist = r * n_missing
            glon, glat, baz = shoot(lon2, lat2, gamma_h + 180, dist)
            
            # endpoint coordinates (-24h, 0h, 6h)
            lon_0, lon_i, lon_i1 = glon, lon[i], lon[i+1]
            lat_0, lat_i, lat_i1 = glat, lat[i], lat[i+1]
#            lo0, lo1, lo2 = glon, lon[i], lon[i+1]
#            la0, la1, la2 = glat, lat[i], lat[i+1]
        
        if i >= 4:  # >= four preceding segments
            
            # endpoint coordinates (-24h, 0h, 6h)
            lon_0, lon_i, lon_i1 = lon[i-4], lon[i], lon[i+1]
            lat_0, lat_i, lat_i1 = lat[i-4], lat[i], lat[i+1]
#            lo0, lo1, lo2 = lon[i-4], lon[i], lon[i+1]
#            la0, la1, la2 = lat[i-4], lat[i], lat[i+1]

        # segment endpoints
        lo0[i], lo1[i], lo2[i] = lon_0, lon_i, lon_i1
        la0[i], la1[i], la2[i] = lat_0, lat_i, lat_i1
        
        # warmup 4-segments (24h) variables
#        lo0[i], la0[i] = lon_0, lat_0
        _, vseg[i], vxseg[i], vyseg[i] = get_vmean(lat_0, lon_0, 
                                                   lat_i, lon_i, 
                                                   24) #dt[i-4:i].sum())
        aseg[i] = GeoAzimuth(lat_0, lon_0, lat_i, lon_i)
#        aseg[i] = calculate_azimut(lon_0, lat_0, lon_i, lat_i)
        
        # warmup last-segment (6h) variables
        pseg[i]  = pres[i-1]      
        wseg[i]  = wind[i-1]
        rseg[i]  = rmw[i-1]
        lseg[i]  = lat_i 
        laseg[i] = np.abs(lat_i)  
    
        # target segment (6h) variables
        _, v, vx, vy = get_vmean(lat_i, lon_i, lat_i1, lon_i1, 
                                 dt[i:i+1].sum())
        dv[i]  = v - vseg[i]             # [km/h]
        dvx[i] = vx - vxseg[i]
        dvy[i] = vy - vyseg[i]
        dp[i]  = pres[i] - pres[i-1]     # [mbar]
        dw[i]  = wind[i] - wind[i-1]     # [kt]
        dr[i]  = rmw[i] - rmw[i-1]       # [nmile]
        dl[i]  = lat_i1 - lat_i          # [º]
        dla[i] = np.abs(dl[i])
        
        # angle variation
        ang1 = aseg[i]
        ang2 = GeoAzimuth(lat_i, lon_i, lat_i1, lon_i1)
#        ang2 = calculate_azimut(lon_i, lat_i, lon_i1, lat_i1)
        dt_ang = ang2 - ang1            # [º]
        sign = np.sign(lseg[i])         # hemisphere: north (+), south (-)
        
        if (ang2 > ang1) & (dt_ang < 180):      da[i] = sign * (dt_ang)
        elif (ang2 > ang1) & (dt_ang > 180):    da[i] = sign * (dt_ang - 360)
        elif (ang2 < ang1) & (dt_ang > -180):   da[i] = sign * (dt_ang)
        elif (ang2 < ang1) & (dt_ang < -180):   da[i] = sign * (dt_ang + 360) 

    # add to dataframe
    df_['vseg'] = vseg / 1.852      # [kt]
#    df_['vxseg'] = vxseg / 1.852
#    df_['vyseg'] = vyseg / 1.852
    df_['dvseg'] = dv / 1.852
#    df_['dvxseg'] = dvx / 1.852
#    df_['dvyseg'] = dvy / 1.852
    df_['pseg'] = pseg              # [mbar]
    df_['dpseg'] = dp
    df_['wseg'] = wseg              # [kt, 1-min avg]
    df_['dwseg'] = dw
    df_['rseg'] = rseg              # [nmile]
    df_['drseg'] = dr
    df_['aseg'] = aseg              # [º]
    df_['daseg'] = da 
    df_['lseg'] = lseg              # [º]
    df_['laseg'] = laseg
    df_['dlaseg'] = dla
    df_['lon_w'] = lo0             # warmup origin
    df_['lat_w'] = la0
    df_['lon_i'] = lo1             # target origin
    df_['lat_i'] = la1
    df_['lon_t'] = lo2             # target endpoint
    df_['lat_t'] = la2

    return df_ 

def stopmotion_interpolation(df_seg, st=None,
                             t_warm=24, t_seg=6, t_prop=42):
    '''
    df_seg      - (pandas.DataFrame)  Stopmotion parameterized units: 
                  (vseg,pseg,wseg,rseg,laseg,dvseg,dpseg,dwseg,drseg,daseg)
    st          - (pandas.DataFrame)  real storm
                - "None" for MDA segments (unrelated to historic tracks)
    t_warm      - warmup period [hour]
    t_seg       - target period [hour]
    t_prop      - propagation period [hour]
                  
    Generation of SWAN cases, cartesian coordinates (SHyTCWaves configuration)
    A.  Warmup period (24h): over the negative x-axis ending at (x,y)=(0,0)
    B.  Target period (6h): starting at (x,y)=(0,0)
    C.  Propagation period (42h): no track coordinates (no wind forcing)
        
    returns:  st_list (pandas.DataFrame)    
    '''
    
    # sign hemisphere (+north, -south)
    if isinstance(st, pd.DataFrame):    # historic track
        sign = np.sign(st['lat'][0])
        method = st.attrs['method']
        center = st.attrs['center']
    else:                               # mda cases 
        sign = 1
        method = 'mda'
        center = 'mda'

    # remove NaNs
    df = df_seg.dropna()
    
    # number of stopmotion events
    N = df.shape[0]

    # list of SWAN cases (paramterized events)
    st_list, we_list = [], []
    for i in range(N):
        
        seg_i = df.iloc[i]

        # stopmotion unit parameters
        vseg  = seg_i['vseg']*1.852     # [km/h]
        dvseg = seg_i['dvseg']*1.852
        pseg  = seg_i['pseg']           # [mbar]
        dpseg = seg_i['dpseg']
        wseg  = seg_i['wseg']           # [kt, 1-min avg]
        dwseg = seg_i['dwseg']
        rseg  = seg_i['rseg']           # [nmile]
        drseg = seg_i['drseg']
        daseg = seg_i['daseg']          # [º]
        laseg = seg_i['laseg']          # [º]
        
        # vmean criteria for SWAN computational timestep [minutes]
        seg_vmean = vseg + dvseg
        if (vseg > 20) or (seg_vmean > 20):     dt_comp = 10
        else:                                   dt_comp = 20
#        if seg_vmean < 20:      dt_comp = 20
#        elif seg_vmean >= 20:   dt_comp = 10

        # time array for SWAN input
        ts = t_warm + t_seg + t_prop        # [h] simulation duration 
        ts = np.asarray(ts) * 60 / dt_comp  # [] intervals of computation
        
        ts_warmup = int(t_warm * 60 / dt_comp)
        ts_segment = int(t_seg * 60 / dt_comp)
    
        # random initial date
        date_ini = pd.Timestamp(1999, 12, 31, 0)
        time_input = pd.date_range(date_ini, periods=int(ts), 
                                    freq='{0}MIN'.format(dt_comp))
        time_input = np.array(time_input)

        # vortex input variables
        x =     np.full(int(ts), np.nan)    # [m]
        y =     np.full(int(ts), np.nan)    # [m]
        vmean = np.full(int(ts), np.nan)    # [km/h]
        ut =    np.full(int(ts), np.nan)
        vt =    np.full(int(ts), np.nan)
        p0 =    np.full(int(ts), np.nan)    # [mbar]
        vmax =  np.full(int(ts), np.nan)    # [kt, 1-min avg]
        rmw =   np.full(int(ts), np.nan)    # [nmile]
        lat =   np.full(int(ts), np.nan)    # [º]

        # (A) preceding 24h segment: over negative x-axis ending at (x,y)=(0,0)
        
        for j in np.arange(0, ts_warmup):
            
            if j==0:    x[j] = - vseg * 24 * 10**3
            else:       x[j] = x[j-1] + vseg * (dt_comp/60) *10**3
            y[j]    = 0
            vmean[j]= vseg
            ut[j]   = vseg
            vt[j]   = 0
            p0[j]   = pseg
            vmax[j] = wseg
            rmw[j]  = rseg
            lat[j]  = laseg
            
        # (B) target 6h segment: starting at (x,y)=(0,0)
        
        for j in np.arange(ts_warmup, ts_warmup+ts_segment):
        
            vel  = vseg + dvseg     # [km/h]
            velx = vel * np.sin((daseg * sign + 90) * np.pi/180)
            vely = vel * np.cos((daseg * sign + 90) * np.pi/180)

            x[j]    = x[j-1] + velx * (dt_comp/60) *10**3
            y[j]    = y[j-1] + vely * (dt_comp/60) *10**3
            vmean[j]= vel
            ut[j]   = velx
            vt[j]   = vely
            p0[j]   = pseg + dpseg
            vmax[j] = wseg + dwseg
            rmw[j]  = rseg + drseg
            lat[j]  = laseg
        
        # (C) propagation 42h segment: remaining values of data arrays
        
        # store dataframe
        st_seg = pd.DataFrame(index=time_input,
                              columns=['x','y','vf','vfx','vfy',\
                                       'pn','p0','vmax','rmw','lat'])
        
        st_seg['x']        = x              # [m]
        st_seg['y']        = y               
        st_seg['lon']      = x              # (idem for plots)
        st_seg['lat']      = y             
        st_seg['vf']       = vmean/1.852    # [kt]
        st_seg['vfx']      = ut/1.852      
        st_seg['vfy']      = vt/1.852      
        st_seg['pn']       = 1013           # [mbar]
        st_seg['p0']       = p0             
        st_seg['vmax']     = vmax           # [kt]
        st_seg['rmw']      = rmw            # [nmile]
        st_seg['latitude'] = lat*sign       # [º]
        
        # add metadata
        st_seg.attrs = {
                'method':   method,
                'center':   center,
                'override_dtcomp': '{0} MIN'.format(dt_comp),
                'x0':       0,
                'y0':       0,
                'p0':       'mbar',
                'vf':       'kt',
                'vmax':     'kt, 1-min avg',
                'rmw':      'nmile',
                }

        # append to stopmotion event list
        st_list.append(st_seg)
        
        # generate wave event (empty)
        we = pd.DataFrame(index=time_input, 
                          columns=['hs', 't02', 'dir', 'spr', 'U10', 'V10'])
        we['level'] = 0
        we['tide'] = 0
        we_list.append(we)
        
    return st_list, we_list


###############################################################################
# STOPMOTION database
# generation of a historical database of stopmotion units (SHyTCWaves 
# methodology) using all available IBTrACS official centers and basins.
# Each storm is preprocessed, gaps filled with estimates (maxwinds, radii),
# interpolated at 6h timesteps, and parameterized as 24h preceding segment
# {v,p,w,r,lat} and 6h target segment {dv,dp,dw,dr,da}
###############################################################################

def stopmotion_database(allcenters=['USA','TOKYO','CMA','HKO','NEWDELHI', \
                                    'REUNION','BOM','WELLINGTON','NADI','WMO'],
                        p_save='', filename=''):
    '''
    Historical IBTrACS storm database is loaded, to loop over official RMSCs
    and their official basins, to generate for each storm the proposed 
    stopmotion methodology parameterized segment units (24h warmup + 6h target)
    
    (A) STORM PREPROCESSING:
        remove NaT, remove NaN, longitude convention [0,360º], timestep, 
        category, basin, convert (X)-min to 1-min avg winds (depends on center)
        
    (B) STORM INTERPOLATION:
        gaps filled with maxwinds/radii estimates, interpolate every 6h,
        mean constant values per segment
        
    (C) STOPMOTION SEGMENTS:
        parameterize 24h preceding warmup and 6h target segments
        Warmup segment (24h):
            4 segments to define start/end coordinates, {Vmean, relative angle}
            last 4th segment, mean {Pmin, Wmax, Rmw}, endpoint {lat}
        Target segment (6h): {dP,dV,dW,dR,dAng}
            
        * Absolute value of latitude is stored, start of target segment
        * Relative angle is referenced to the geographic north (southern 
          hemisphere is multiplied with a factor of -1)
    
    returns:  df_database (pandas.Dataframe)
    '''
    
    # get ibtracs database
    filename = 'ibtracs.nc'
    p_ibtracs = op.join(op.dirname(op.realpath(__file__)),'resources',filename)
    xds_ibtracs = xr.open_dataset(p_ibtracs)
    ibtracs_version = xds_ibtracs.attrs['product_version']
    
    # set storms id as coordinate
    xds_ibtracs['stormid'] = (('storm'), xds_ibtracs.storm.values)
#    xds_ibtracs = xds_ibtracs.set_coords('stormid')

    # loop for all centers and official basins
#    allcenters = ['USA','TOKYO','CMA','HKO','NEWDELHI','REUNION','BOM',\
#                  'WELLINGTON','NADI','WMO']
    import datetime

    df_list = []
    for center in allcenters:
                
        # get center's official basins ['NA','SA','WP','EP','SP','NI','SI']
        basin_ids = d_vns_basinscenter[center]

        # loop for all basins
        for basin in basin_ids:
            
            print(center, basin, '  start time:', datetime.datetime.now())

            # extract storms at basin X
            xds_basin = Extract_basin_storms(xds_ibtracs, basin)
            
            for i in range(xds_basin.storm.size):
                storm = xds_basin.sel(storm=i)
                stid = storm.stormid.values
                                
                # (A) -- STORM PREPROCESSING
                # returns: {center,basin,dist,lo,la,move,vf,p,w,r,cat,ts}
                df = historic_track_preprocessing(storm, center=center, 
                                                  database_on=True)
                
                # separate storm blocks (independent if gaps > 24h)
                pos = np.where(df['timestep'].values > 24)[0]
                df_blocks = []
                
                if pos.size == 0:
                    df_blocks.append(df)
                    
                if pos.size == 1:
                    df_blocks.append(df.iloc[:pos[0]+1])
                    df_blocks.append(df.iloc[pos[0]+1:])
                
                elif pos.size > 1:
                    df_blocks.append(df.iloc[0:pos[0]+1])
                    for i in range(pos.size-1):
                        loc_ini, loc_end = pos[i]+1, pos[i+1]+1
                        df_blocks.append(df.iloc[loc_ini:loc_end])
                    df_blocks.append(df.iloc[loc_end:])


                for df in df_blocks:
                    # when storm is not empty
                    if df.shape[0] > 1:
                        
#                        print(stid, df['timestep'].values)
                        
                        # (B) -- STORM INTERPOLATION
                        # generate filled, interpolated track (constant segments)
                        # {center,basin,dist,lo,la,move,vf,pn,p0,w,r,_fill}
                        
                        dt_interp = 6*60    # stopmotion segment interval (6h)
                        
                        st, _ = historic_track_interpolation(df, dt_interp, 
                                                             interpolation=False, 
                                                             mode='mean',
                                                             )
                        # (C) -- STOPMOTION SEGMENTS
                        # generate stopmotion segments (24h warmup + 6h target)
                        # + {lo0,la0,vseg,pseg,wseg,rseg,laseg,dv,dp,dw,dr,da}
                        df_seg = storm2stopmotion(st)
                        
                        # (option) remove on-land data ("dist2land"==0)
                        # (option) remove vmaxfill, rmwfill to ignore estimates
        
                        # add storm id 
                        df_seg['stid'] = stid 
                        
                        # keep stopmotion parameters
                        df_ = df_seg[['center','basin','stid','dist2land', \
                                      'vseg', 'dvseg', 'pseg', 'dpseg', \
                                      'wseg', 'dwseg', 'rseg', 'drseg', \
                                      'aseg', 'daseg', 'lseg', 'laseg', \
                                      'vmaxfill','rmwfill','lon_i','lat_i']]
                        
                        # append to list
                        df_list.append(df_)
    
    # concatenate dataframes
    df_database = pd.concat(df_list)
            
    # add medatada
    df_database.attrs = {
            'version':      ibtracs_version,
            'center':       'ibtracs source',
            'stid':         'storm identifier',
            'dist2land':    'km, distance to nearest land',
            'vseg':         'kt, translation velocity',
            'pseg':         'mbar, minimum central pressure',
            'wseg':         'kt, 1-min avg, maximum wind speed',
            'rseg':         'nmile, radii of maximum wind speed',
            'daseg':        'º, angle change of forward direction',
            'laseg':        'º, absolute latitude',
            }
    
    # save to pickle
    if p_save:      df_database.to_pickle(op.join(p_save, filename))
    
    return df_database


###############################################################################
# STOPMOTION ensemble
# functions to collect the 6h-segments from library cases to obtain the hybrid
# storm track (closest to the real track), those segments are rotated and 
# assigned time-geographical coordinates, and point-based ensemble is calculated
###############################################################################

def find_analogue(df_library, df_case, ix_weights):
    '''
    Finds the minimum distance in a n-dimensional normalized space, 
    corresponding to the shytcwaves 10 parameters: 
        {pseg, vseg, wseg, rseg, dp, dv, dw, dr, dA, lat}
    
    df_library  - (pandas.Dataframe) library parameters
    df_case     - (pandas.Dataframe) target case parameters
    ix_weights  - (array) columns indices weigth factors
    '''
    
    # remove NaNs from storm segments
    df_case = df_case.dropna()
    data_case = df_case[['daseg','dpseg','pseg','dwseg','wseg', \
                         'dvseg','vseg','drseg','rseg','laseg']].values
    
    # library segments parameter
    data_lib = df_library[['daseg','dpseg','pseg','dwseg','wseg', \
                           'dvseg','vseg','drseg','rseg','laseg']].values
    ix_scalar = [1,2,3,4,5,6,7,8,9]
    ix_directional = [0]

    # get indices of nearest n-dimensional point
    ix_near = nearest_indexes_weighted(data_case, data_lib, 
                                       ix_scalar, ix_directional, ix_weights)

    return ix_near
    
def analogue_endpoints(df_seg, df_analogue):
    '''
    adds segments endpoints from original target origin coordinates
    
    df_seg          - (pandas.DataFrame) historic parameterized storm
    df_analogue     - (pandas.DataFrame) analogue segments
    '''
    
    # remove NaNs
    df_seg = df_seg[~df_seg.isna().any(axis=1)]
    
    # get hemishpere
    # relative angles are multiplied by "sign" to account for hemisphere
    sign = np.sign(df_seg.lat.values[0])
    
    # get historic variables
    lon0 = df_seg.lon.values                    # target origin coords
    lat0 = df_seg.lat.values
    aseg = df_seg.aseg.values                   # warmup azimuth
    daseg = df_seg.daseg.values * sign          # target azimuth variation
    
    # get analogue variables
    vseg1_an = df_analogue.vseg.values          # warmup velocity [kt]
    dvseg_an = df_analogue.dvseg.values
    vseg2_an = vseg1_an + dvseg_an              # target velocity
    daseg_an = df_analogue.daseg.values * sign  # target azimuth variation
    
    # new variables
    az1_an = np.full(lon0.shape, np.nan)
    glon1_an = np.full(lon0.shape, np.nan)
    glon2_an = np.full(lon0.shape, np.nan)
    glat1_an = np.full(lon0.shape, np.nan)
    glat2_an = np.full(lon0.shape, np.nan)
    
    for i in range(df_seg.shape[0]):
        
        # azimuth angles for ensemble
        az2 = aseg[i] + daseg[i]        # target azimuth (historic-fixed)
        az1 = az2 - daseg_an[i]         # warmup azimuth (stopmotion analogue)
        
        # shoot backwards to warmup origin
        dist1 = vseg1_an[i] * 1.852 * 24  # [km]
        glon1, glat1, baz = shoot(lon0[i], lat0[i], az1 + 180, dist1)
        
        # shoot forwards to target endpoint
        dist2 = vseg2_an[i] * 1.852 * 6   # [km]
        glon2, glat2, baz = shoot(lon0[i], lat0[i], az2 + 180 - 180, dist2)
        
        # store 
        az1_an[i] = az1
        glon1_an[i] = glon1
        glon2_an[i] = glon2
        glat1_an[i] = glat1
        glat2_an[i] = glat2
            
    # longitude convention
    glon1_an[glon1_an<0] += 360
    glon2_an[glon2_an<0] += 360
    
    # add to dataframe
    df_analogue['aseg'] = az1_an
    df_analogue['lon_w'] = glon1_an    # warmup origin
    df_analogue['lat_w'] = glat1_an
    df_analogue['lon_i'] = df_seg.lon.values
    df_analogue['lat_i'] = df_seg.lat.values
    df_analogue['lon_t'] = glon2_an
    df_analogue['lat_t'] = glat2_an
    
    return df_analogue


def find_analogue_grid_coords(path_mda, icase):
    '''
    Finds the analogue segment case grid size (small/medium/large) and
    calculates the polar coordinates grid
    
    rings
    res
    '''
    
    # load MDA indices (grid sizes)
    # TODO: path
#    p_mda = r'/media/administrador/Disco/GitLab/BlueMath/data/stopmotion/'
#    xds_ind_mda = xr.open_dataset(op.join(p_mda, 'shytcwaves_mda_indices.nc'))
    xds_ind_mda = xr.open_dataset(op.join(path_mda, 'shytcwaves_mda_indices_clean.nc'))
    
    # get grid code
    pos_small = np.where(icase == xds_ind_mda.indices_small)[0]
    pos_medium = np.where(icase == xds_ind_mda.indices_medium)[0]
    pos_large = np.where(icase == xds_ind_mda.indices_large)[0]
    
    # TODO: number of rings variable
    if len(pos_small)==1:       num, rings = 0, 29
    elif len(pos_medium)==1:    num, rings = 1, 33
    elif len(pos_large)==1:     num, rings = 2, 33
    
    # get sizes
    xsize_1 = xds_ind_mda.x_size_left.values[num]*10**3
    xsize_2 = xds_ind_mda.x_size_right.values[num]*10**3
    ysize_1 = xds_ind_mda.y_size_left.values[num]*10**3
    ysize_2 = xds_ind_mda.y_size_right.values[num]*10**3

    #  (cartesian convention)
    # computational grid extent
    res = 15 *10**3   # km to m
    # TODO: resolution variable

    # create domain centered at (x,y)=(0,0)
    lon = np.arange(-xsize_1, xsize_2, res)
    lat = np.arange(-ysize_1, ysize_2, res)

    # custom output coordinates
    _, _, _, _, rc, thetac = generate_polar_coords(lon, lat, res, 
                              rings=rings, radii_ini=5000, angle_inc=5)
        
    return rc, thetac

#def analogue2location(path_library, df_analogue, path, df_seg, 
#                      st, #rc, thetac, 
#                      storm_name='OFA'):
#    '''
#    Function to access the library analogue cases for a given storm track,
#    to allocate them with absolute coordinates (geographical+rotation and time)
#    and to store them for postprocessing
#    
#    path_library- (path) library of prerun segments
#    path        - (path) project directory
#    
#    df_analogue - (pandas.DataFrame) analogue prerun segments from library
#    df_seg      - (pandas.DataFrame) storm 6h-segments parameters 
#    st          - (pandas.DataFrame) storm track interpolated every 6h
#    
##    rc          - (xarray) radii of polar control points (shape)
##    thetac      - (xarray) angle of polar control points 
#    
#    returns:    (xarray.Dataset) saved files
#    '''
#    
#    # folder path for project ensemble
#    p_ensemble = op.join(path, 'ensemble')
#    if not op.isdir(p_ensemble):  os.mkdir(p_ensemble)
#    
#    # remove NaN
#    df_seg = df_seg[~df_seg.isna().any(axis=1)]
#    df_analogue = df_analogue[~df_analogue.isna().any(axis=1)]
#    
#    # assign time
#    df_seg['time'] = df_seg.index.values
#    
#    # get hemisphere
#    sign = np.sign(df_seg.lat.values[0])
#    
#    # print time
#    print('Initial time: ', datetime.datetime.now())
#    
#    # -------------------------------------------------------------------------
#    # evaluate each REAL time coordinate of the storm (st)
#    # (index segments inform about the time coordinates)
#
#    for iseg in range(df_seg.shape[0]):
#        
#        # get storm segment 'i'
#        df_case = df_seg.iloc[iseg]
#        seg_time = np.datetime64(df_case.time)  # timestamp to datetime64
##        seg_ang = df_case.aseg                  # relative angle   REMOVE!!!!
#    
#        # ---------------------------------------------------------------------
#        # get storm coordinates at "seg_time"
#        ind_time_st = np.where(st.index.values == seg_time)[0][0]
#        # storm coordinates
#        st_time = st.index.values[ind_time_st]
#        # target origin coordinates (real=fixed)
#        st_lon = st.iloc[ind_time_st].lon
#        st_lat = st.iloc[ind_time_st].lat
#    
#        # ---------------------------------------------------------------------
#        # get analogue segment from library
#        iseg_analogue = df_analogue.index.values[iseg]
#        seg_ang = df_analogue.iloc[iseg].aseg                  # RELATIVE ANGLE
#
##         df_analogue, pos_analogue = find_analogue(df_library, df_case)
#        filename = 'spec_outpts_main.nc'
#        p_analogue = op.join(path_library, 
##                              'cases/{0:04d}'.format(iseg_analogue), filename)
#                             '{0:04d}'.format(iseg_analogue), filename)
#        # load file (*time:48, *freq:38, *dir:72, point)
#        seg_sim = xr.open_dataset(p_analogue)
#        rc, thetac = find_analogue_grid_coords(iseg_analogue)
#
#
#        # ---------------------------------------------------------------------
##        # PLOT 1
##        seg_sim1 = seg_sim.copy()
##        seg_sim1 = seg_sim1.rename({'frequency': 'freq'})
##        seg_sim1 = seg_sim1.rename({'direction': 'dir'})
##        
##        hs0 = seg_sim1.spec.hs()
##        hs1 = seg_sim1.spec.hs().max(axis=1)
##        
##        plt.subplot(131)
##        plt.scatter(seg_sim1.lon_pts, seg_sim1.lat_pts, c=hs0.values[:,20], cmap='hot');
##        plt.plot(st_mda[numseg].lon, st_mda[numseg].lat); plt.axis('equal')
##        plt.title('real', fontweight='bold'); plt.colorbar();
#
#        # if southern hemisphere, mirror axis
#        if sign == -1:
#            
#            # reorder points (swan stores radii outwards, angle clockwise)
#            n_point = seg_sim.point.shape[0]
#            n_freq = rc.shape[1]
#            n_dir_half = int(rc.shape[0] / 2)
#            
#            def reorder_polar_points(n_freq, n_dir_half, xo, p=[]):
#                
#                for i in range(n_dir_half):
#                    p.append(np.arange(xo, xo + n_freq))
#                    xo -= n_freq
#                
#                return p
#
#            xo = n_point/2 - n_freq    # first half
#            p = reorder_polar_points(n_freq, n_dir_half, xo, p=[])
#            xo = n_point - n_freq      # second half
#            p = reorder_polar_points(n_freq, n_dir_half, xo, p=p)
#            
#            ind_order = np.reshape(np.array(p, dtype=int), -1)  # 1-dim points
#        
#            seg_sim = xr.Dataset(
#                    {'efth': (('freq','dir','point','time'), seg_sim.efth.values),
#                     'lon_pts': ('point', seg_sim.lon_pts.values[ind_order]),
#                     'lat_pts': ('point', seg_sim.lat_pts.values[ind_order]),
#                            },
#                    {'freq': seg_sim.frequency.values,
#                     'dir': seg_sim.direction.values,
#                     'time': seg_sim.time.values,
#                            })
#            
#        seg_sim['hs'] = seg_sim.spec.hs()#.max(axis=1)
#        
##        # plot
##        plt.subplot(132)
##        plt.scatter(seg_sim.lon_pts, seg_sim.lat_pts, 
##                    c=seg_sim.hs.values[:,20], cmap='hot'); 
##        plt.plot(st_list_analogue[iseg].x, st_list_analogue[iseg].y); 
##        plt.axis('equal')
##        plt.title('flipped', fontweight='bold'); plt.colorbar();
#
#        # ---------------------------------------------------------------------
#        # time array
#        hour_intervals = seg_sim.time.size
#        time = [st_time + np.timedelta64(1,'h')*i for i in range(hour_intervals)]
#        time_array = np.array(time)
#    
#        # spec array
#        # spec array
#        spec_array = np.reshape(seg_sim.efth.values, 
#                                (seg_sim.freq.size, seg_sim.dir.size, 
#                                 rc.T.shape[0], rc.T.shape[1],
#                                 seg_sim.time.size)
#                                )
##        spec_array = np.reshape(seg_sim.efth.values, 
##                                (seg_sim.time.size, seg_sim.frequency.size, 
##                                 seg_sim.direction.size, 
##                                 rc.shape[0], rc.shape[1])
##                                )
#
#        # get lon,lat arrays
#        xds_x = np.reshape(seg_sim.lon_pts.values, rc.T.shape)  # to (r,ang)
#        xds_y = np.reshape(seg_sim.lat_pts.values, rc.T.shape)
#        xds_rad = (np.sqrt(xds_x**2 + xds_y**2))      # radii array
#        xds_ang = (np.full(xds_rad.shape, np.nan))    # angle array
#        
#        deg = np.arctan2(xds_y, xds_x) *180/np.pi
#        deg[deg<0] += 360
#        xds_ang = np.round(np.mod(-deg + 90, 360), 0)
#        
#        lon_arr = np.full(xds_rad.shape, np.nan)
#        lat_arr = np.full(xds_rad.shape, np.nan)
#        for i in range(xds_rad.shape[0]):
#            for j in range(xds_rad.shape[1]):
#                lon_arr[i,j], lat_arr[i,j],_ = shoot(st_lon, st_lat, 
#                                                     xds_ang[i,j] +seg_ang -90, 
#                                                     xds_rad[i,j]/1000)  # [km]
#        lon_arr[lon_arr<0] += 360    # [0º,360º] longitude convention
#        
#    
#        # ---------------------------------------------------------------------
#        # dataset with absolute coordinates
#        spec_arr = np.expand_dims(spec_array, axis=0)
#
#        xd = xr.Dataset({
#            'efth': (('case','freq','dir','radii','angle','time'), spec_arr),
#            'lon': (('radii','angle',), lon_arr),
#            'lat': (('radii','angle',), lat_arr),
#        }, coords={
#            'case': np.atleast_1d(iseg),
#            'freq': seg_sim.freq.values,
#            'dir': seg_sim.dir.values,
#            'angle': thetac[:,0],
#            'radii': rc[0,:],
#            'time': time_array,
#        })
##        xd = xr.Dataset({
##            'efth': (('case','time','freq','dir','angle','radii'), spec_arr),
##            'lon': (('angle','radii'), lon_arr),
##            'lat': (('angle','radii'), lat_arr),
##        }, coords={
##            'case': np.atleast_1d(iseg),
##            'time': time_array,
##            'freq': seg_sim.frequency.values,
##            'dir': seg_sim.direction.values,
##            'radii': rc[0,:],
##            'angle': thetac[:,0],
##        })
#    
##        # ---------------------------------------------------------------------
##        # PLOT 3
##        hs3 = xd.spec.hs()#.max(axis=-1)
##        
##        plt.subplot(133)
##        # plt.pcolor(xd.lon, xd.lat, hs3.values[0,20,:,:], cmap='hot')
##        plt.scatter(xd.lon, xd.lat, c=hs3.values[0,:,:,20], cmap='hot'); plt.axis('equal')
##        plt.plot([df_analogue.iloc[iseg].lon_w, df_analogue.iloc[iseg].lon_i, df_analogue.iloc[iseg].lon_t],
##                 [df_analogue.iloc[iseg].lat_w, df_analogue.iloc[iseg].lat_i, df_analogue.iloc[iseg].lat_t])
##        # plt.plot(st_mda[numseg].lon, st_mda[numseg].lat); plt.axis('equal')
##        plt.title('rotated', fontweight='bold'); plt.colorbar();
##        
##        plt.gcf().set_size_inches(15,5)
##        # ---------------------------------------------------------------------
#
#        # store
#        filename = 'storm_{0}_segment{1}.nc'.format(storm_name, iseg)
#        xd.to_netcdf(op.join(p_ensemble, filename))
#        
#        print(iseg, datetime.datetime.now())
#    print('done')


#def analogue_point_spectra(p_ensemble, df_seg, p_lon, p_lat,
#                           max_dist = 60, storm_name='OFA'):
#    '''
#    Function to extract the wave spectra from each parameterized segment of the
#    storm track, at the closest grid control point to the target point.
#    
#    p_ensemble  - (path) ensemble project directory
#    df_seg      - (pandas.DataFrame) storm 6h-segments parameters 
#    p_lon,p_lat - (float) target point for ensemble    
#    
#    max_dist    - (float) maximum distance to assign interpolated spectra [km]
#    
#    returns:    (xarray.Dataset) closest point-based raw wave spectra segments
#    '''
#        
#    # earth radius [km]
#    RE = 6378.135   
#    dd = 1.5    # interval [degree]
#    
#    # remove NaN
#    df_seg = df_seg[~df_seg.isna().any(axis=1)]
##     # storm parameterized segments indices
##     indices_df_seg = df_seg.dropna().index.values
#
#    # extract spectra at point (case,time,freq,dir)
#    for iseg in range(df_seg.shape[0]):
##     for pos in range(indices_df_seg.size - 1):
#    
#        sys.stdout.write('\rCase %d' %iseg)
#        sys.stdout.flush()
#        
#        # load storm segment in absolute coordinates 
#        filename = 'storm_{0}_segment{1}.nc'.format(storm_name, iseg)
#        xds = xr.open_dataset(op.join(p_ensemble, filename))
#    
#        pos = np.where((xds.lon.values < p_lon+dd) & (xds.lon.values > p_lon-dd) & \
#                       (xds.lat.values < p_lat+dd) & (xds.lat.values > p_lat-dd))
#        
#        xds_lon = xds.lon.values[pos[0],pos[1]]
#        xds_lat = xds.lat.values[pos[0],pos[1]]
#        
#        # ---------------------------------------------------------------------
#        # interpolate target point in output grid control points (i,j)
#        # method: nearest.... weighted average also???
#    
#        # calculate distance between consecutive coordinates [km]
#        dist = np.full(xds_lat.shape, np.nan)
#        for ii in range(xds_lat.shape[0]):
#            arcl_h, _ = gc_distance(p_lat, p_lon, 
#                                    xds_lat[ii], xds_lon[ii])
#            dist[ii] = arcl_h * np.pi / 180.0 * RE
##        dist = np.full(xds.lat.values.shape, np.nan)
##        for ii in range(xds_lat.shape[0]):
##            for jj in range(xds_lat.shape[1]):
##                arcl_h, _ = gc_distance(p_lat, p_lon, 
##                                        xds_lat[ii,jj], xds_lon[ii,jj])
##                dist[ii,jj] = arcl_h * np.pi / 180.0 * RE
#    
#        # TODO: remove duplicate angle in the first np.arange(0,360,5)
#        
#        if np.amin(dist) < max_dist:  # assign closest node
#        
#            # indices for closest node
#            ind = np.where(dist==np.amin(dist))[0][0]   # (radii,angle)
#            
#            xds_sel = xds.isel(radii=pos[0][ind], angle=pos[1][ind], 
##                               ).squeeze(['angle','radii']
#                               ).drop(['angle','radii','lon','lat'])
#        
#            
##            ind_dist_angle, ind_dist_radii = pos[np.where(dist==np.amin(dist))]
##    
##            # dataset selection
##            if ind_dist_angle.size > 1:   
##                ind_dist_angle = ind_dist_angle[0]
##                ind_dist_radii = ind_dist_radii[0]
##    
##            xds_sel = xds.isel(angle=ind_dist_angle, 
##                               radii=ind_dist_radii
##                               ).squeeze(['angle','radii']
##                               ).drop(['angle','radii','lon','lat'])
#        
#        # TODO check it works --> remove squeeze
#        else:  # assign NaNs when there is no closest point for interpolation
#            xds_sel = xds.isel(angle=0, radii=0
##                               ).squeeze(['angle','radii']
#                               ).drop(['angle','radii','lon','lat'])
#            xds_sel['efth'] = (('case','freq','dir','time'),
#                                np.full(xds_sel['efth'].shape, np.nan))
#            
#        # store
#        if iseg == 0:  xds_segs = xds_sel
#        if iseg > 0:   xds_segs = xr.merge([xds_segs, xds_sel])
#        
#        print(iseg, datetime.datetime.now(), xds_segs.dims, xds_sel.dims)
#
#    return xds_segs
   
#def analogue_point_reconstruction(xds_segs, mask_mode='none'):
#    '''
#    Generate the point-based ensemble using the closest nodes wave spectra over 
#    time, apply a dynamic obstacle mask, and calculate the envelope of 
#    reconstructed directional wave spectra over time.
#    
#    xds_segs    - (xarray.Dataset) closest point-based raw wave spectra segments
#    mask_mode   - 'none' = no mask
#                  'fixed' = obstacle sector mask (constant in time)
#                  'dynamic' = obstacle mask plus distance criteria over time
#    
#    returns:    (xarray.Dataset) closest point-based masked spectra envelope
#    '''
#    
#    # WARNING: "wavespectra" library expect names to be 'efth','freq','dir'
#    # rename dataset
##    xds_segs = xds_segs.rename({'spectra': 'efth'})
#    
#    # TODO: reorder variables to make sure order
#    
#    # reconstruct at every case, over time 
#    xds_segs['hs'] = (('case','point','time'), xds_segs.spec.hs())
#    
#    xds_pt_ls = []
#    for ip in range(xds_segs.point.size):
#        
#        xds_pt = xds_segs.sel(point=ip)
#        xds_pt = xds_pt.dropna(dim='case', how='all')
#        xds_pt = xds_pt.dropna(dim='time', how='all')
#    
#        time_indices = np.arange(0, xds_pt.time.size)
#        
## TODO
##    elif mask_mode=='fixed':
##        # apply obstacle mask
##    
##    elif mask_mode=='dynamic':
##        # apply obstacle mask
##        
##        # distance criteria: distance to storm center < distance to obstacle
##        
#
#        # select segment case indices (maximum along dim='case')
#        seg_bmu_new = xds_pt.hs.argmax(dim='case').values
#    
#        efth_bmu = xds_pt.efth.values[seg_bmu_new,:,:,time_indices]
#        hs_bmu = xds_pt.hs.values[seg_bmu_new,time_indices]
#        seg_bmu = xds_pt.case.values[xds_pt.hs.argmax(dim='case').values]
#        
#        # store
#        xds_pti = xr.Dataset({
#            'efth': (('point','time','freq','dir'), np.expand_dims(efth_bmu, axis=0)),
#            'hs': (('point','time'), np.expand_dims(hs_bmu, axis=0)),
#            'lon': (('point'), [xds_segs.lon.values[ip]]),
#            'lat': (('point'), [xds_segs.lat.values[ip]]),
#            'bmu': (('point','time'), np.expand_dims(seg_bmu, axis=0)),
#        }, coords={
#            'point': [ip],
#            'time': xds_pt.time.values,
#            'freq': xds_pt.freq.values,
#            'dir': xds_pt.dir.values,
#        })
#        
#        xds_pt_ls.append(xds_pti)
#
#    # merge datasets    
#    xds_pts = xr.merge(xds_pt_ls)
#
#    return xds_pts


def stopmotion_st_bmu(path_library, df_analogue, df_seg, path_mda,
                      st, cp_lon_ls, cp_lat_ls, max_dist=60, 
                      list_out=False, tqdm_out=True, text_out=True):
    '''
    Function to access the library analogue cases for a given storm track,
    calculate distance and angle from the target segment origin to the control
    point (relative coordinate system), and extract the directional wave 
    spectra at the closest node (for every analogue segment)
    
    path_library    - (path) library of prerun segments
    path_mda        - (path) mda sample indices (grid sizes)
    
    df_analogue     - (pandas.DataFrame) analogue prerun segments from library
    df_seg          - (pandas.DataFrame) storm 6h-segments parameters 
    st              - (pandas.DataFrame) storm track interpolated every 6h
    cp_lon/lat_ls   - (list) control point geographical coordinates
    max_dist        - (float) maximum distance [km] to extract closest node
    
    returns:    (xarray.Dataset) wave directional spectra (dim 'case')
    '''

    # remove NaN
    df_seg = df_seg[~df_seg.isna().any(axis=1)]
    df_analogue = df_analogue[~df_analogue.isna().any(axis=1)]
    
    # assign time
    df_seg['time'] = df_seg.index.values
    
    # get hemisphere
    sign = np.sign(df_seg.lat.values[0])

    xds_list = []
#    print('loop start...', datetime.datetime.now())

    if tqdm_out:    array = tqdm(range(df_seg.shape[0]))
    else:           array = range(df_seg.shape[0])

    for iseg in array:     # each segment
        
        # get storm segment 'i'
        df_icase = df_seg.iloc[iseg]
        df_ianalogue = df_analogue.iloc[iseg]
        iseg_analogue = df_ianalogue.name       # analogue id
        aseg = df_ianalogue.aseg                # analogue warmseg azimuth (N)
    
        # ---------------------------------------------------------------------
        # get storm coordinates at "seg_time"
        seg_time = np.datetime64(df_icase.time)  # timestamp to datetime64
        ind_time_st = np.where(st.index.values == seg_time)[0][0]
        
        # storm coordinates (real storm eye)
        st_lon = st.iloc[ind_time_st].lon       # target origin longitude
        st_lat = st.iloc[ind_time_st].lat       # target origin latitude
        st_time = st.index.values[ind_time_st]  # time coordinates
    
        # ---------------------------------------------------------------------
        # get cp coordinates in relative system (radii, angle)
        cp_dist_ls, cp_ang_ls = get_cp_radii_angle(st_lat, st_lon, 
                                                    cp_lat_ls, cp_lon_ls, 
                                                    sign, aseg)
        
        # get SWAN output mask indices (rad, ang)
        mask_rad, mask_ang = get_mask_radii_angle(path_mda, iseg_analogue)
        
        # find closest index
        ix_near = find_nearest(cp_dist_ls, cp_ang_ls, mask_rad, mask_ang)
#        pos_nonan = np.where(np.abs(mask_rad[ix_near] - cp_dist_ls) \
#                             <= max_dist)[0]
        pos_nonan = np.abs(mask_rad[ix_near] - cp_dist_ls) <= max_dist
        
        # ---------------------------------------------------------------------
        # load library Hs reconstructed
        xds_rec = xr.open_dataset(op.join(path_library, 
                                          'library_shytcwaves_bulk_params.nc'))
        
        # extract HS,TP at closest points (case,point,time)
        hs_arr = np.full((ix_near.size, xds_rec.time.values.size), np.nan)
        hs_arr[pos_nonan,:] = xds_rec.hs.values[iseg_analogue,ix_near[pos_nonan],:]
        
        tp_arr = np.full((ix_near.size, xds_rec.time.values.size), np.nan)
        tp_arr[pos_nonan,:] = xds_rec.tp.values[iseg_analogue,ix_near[pos_nonan],:]

        # time array
        hour_intervals = xds_rec.time.size
        time = [st_time + np.timedelta64(1,'h')*i for i in range(hour_intervals)]
        time_array = np.array(time)

        # store dataset
        xds = xr.Dataset({
                'hs': (('case','point','time'), np.expand_dims(hs_arr, axis=0)),
                'tp': (('case','point','time'), np.expand_dims(tp_arr, axis=0)),
                'lon': (('point'), cp_lon_ls),
                'lat': (('point'), cp_lat_ls),
                'ix_near': (('case','point'), np.expand_dims(ix_near, axis=0)),
                'pos_nonan': (('case','point'), np.expand_dims(pos_nonan, axis=0)),
                  },
                  coords={
                          'case': [iseg],
                          'time': time_array,
                          })
        xds_list.append(xds)

#        sys.stdout.flush()
#        sys.stdout.write('\rSegment: {0} / {1}\r'.format(iseg+1, df_seg.shape[0]))
##        print('segment', iseg, datetime.datetime.now())
#
#    print('Segment: ', iseg+1, '/', df_seg.shape[0])
#    print('loop end.....', datetime.datetime.now())
    if list_out:    xds_out = xds_list
    else:       
        # merge
        xds_out = xr.merge(xds_list)
#        print('loop merge...', datetime.datetime.now())
        if text_out: print('Merging bulk envelope...', datetime.datetime.now())
        
        
        # add envelope variables
        xds_bmu = xds_out.copy()
        hsval = xds_bmu.hs.values[:]
        hsval[np.isnan(hsval)] = 0                      # remove nan
        bmu = np.argmax(hsval, axis=0).astype(float)    # bmu indices
        hsmax = np.sort(hsval, axis=0)[-1,:,:]          # max over 'case'
    
        # bmu, hs
        bmu[hsmax==0] = np.nan                          # restitute nans
        hsmax[hsmax==0] = np.nan
        
        xds_out['bmu'] = (('point','time'), bmu)
        xds_out['hsbmu'] = (('point','time'), hsmax)
    
        # tp
        tpmax = np.full(xds_out.hsbmu.shape, np.nan)
        nanmask = ~np.isnan(bmu)
        mesht, meshp = np.meshgrid(np.arange(0,xds_out.time.size), 
                                   np.arange(0,xds_out.point.size))
        
        tpmax[nanmask] = xds_out.tp.values[bmu.ravel()[nanmask.ravel()].astype('int64'), 
                                          meshp.ravel()[nanmask.ravel()], 
                                          mesht.ravel()[nanmask.ravel()]]
        xds_out['tpbmu'] = (('point','time'), tpmax)

        # add swath variables
        hh = hsmax.copy()
        hh[np.isnan(hh)] = 0
        posw = np.argmax(hh, axis=1)
        
        xds_out['hswath'] = (('point'), hsmax[np.arange(0,xds_out.point.size), posw])
        xds_out['tswath'] = (('point'), tpmax[np.arange(0,xds_out.point.size), posw])

    return xds_out

def stopmotion_st_spectra(path_library, df_analogue, df_seg, path_mda,
                         st, cp_lon_ls, cp_lat_ls, cp_names=[], 
                         max_dist=60, list_out=False, tqdm_out=True, text_out=True):
    '''
    Function to access the library analogue cases for a given storm track,
    calculate distance and angle from the target segment origin to the control
    point (relative coordinate system), and extract the directional wave 
    spectra at the closest node (for every analogue segment)
    
    path_library    - (path) library of prerun segments
    path_mda        - (path) mda sample indices (grid sizes)
    
    df_analogue     - (pandas.DataFrame) analogue prerun segments from library
    df_seg          - (pandas.DataFrame) storm 6h-segments parameters 
    st              - (pandas.DataFrame) storm track interpolated every 6h
    cp_lon/lat_ls   - (list) control point geographical coordinates
    max_dist        - (float) maximum distance [km] to extract closest node
    
    returns:    (xarray.Dataset) wave directional spectra (dim 'case')
    '''

#    # folder path for project ensemble
#    p_ensemble = op.join(path, 'ensemble')
#    if not op.isdir(p_ensemble):  os.mkdir(p_ensemble)
    
    # remove NaN
    df_seg = df_seg[~df_seg.isna().any(axis=1)]
    df_analogue = df_analogue[~df_analogue.isna().any(axis=1)]
    
    # assign time
    df_seg['time'] = df_seg.index.values
    
    # get hemisphere
    sign = np.sign(df_seg.lat.values[0])

    # get bmu (wavespectra reconstructed)
    # it provides 'time,bmu,ix_near,pos_nonan' (point,time)
    xds_bmu = stopmotion_st_bmu(path_library, df_analogue, df_seg, path_mda,
                                 st, cp_lon_ls, cp_lat_ls, 
                                 max_dist=max_dist, list_out=list_out, 
                                 tqdm_out=tqdm_out, text_out=text_out)
#    print('bmu identifier obtained for all analogue segments.')
        
#    shy_time = xds_bmu.time.values      # (dim time)
#    shy_bmu = xds_bmu.bmu.values        # (point,time)
    
    # spectral energy
    seg_sim = xr.open_dataset(op.join(path_library, '0000/spec_outpts_main.nc'))
    efth_arr = np.full((seg_sim.frequency.size, seg_sim.direction.size, 
                        xds_bmu.point.size, xds_bmu.time.size), np.nan) # 38,72,cp,t
    
#    print('loop start...', datetime.datetime.now())

    if tqdm_out:    array = tqdm(range(xds_bmu.case.size))
    else:           array = range(xds_bmu.case.size)
    for iseg in array:
    
        # get storm segment 'i'
        df_icase = df_seg.iloc[iseg]
        df_ianalogue = df_analogue.iloc[iseg]
        iseg_analogue = df_ianalogue.name       # analogue id
    
        # ---------------------------------------------------------------------
        # get storm coordinates at "seg_time"
        seg_time = np.datetime64(df_icase.time)  # timestamp to datetime64
        st_time = st.index.values[st.index.values == seg_time][0]
    
        # ---------------------------------------------------------------------
        # get analogue segment from library
        filename = 'spec_outpts_main.nc'
        p_analogue = op.join(path_library, 
                         '{0:04d}'.format(iseg_analogue), filename)
        # load file
        seg_sim = xr.open_dataset(p_analogue)   # freq,dir,2088,48
    
        # time array
        hour_intervals = seg_sim.time.size
        time = [st_time + np.timedelta64(1,'h')*i for i in range(hour_intervals)]
        time_array = np.array(time)
    
        # get intersect time iseg vs xds_bmu
        _, ix_time_st, ix_time_shy = np.intersect1d(time_array, xds_bmu.time.values, 
                                                    return_indices=True)
    
        # find all closest grid points
        shy_inear = xds_bmu.ix_near.values[iseg,:].astype('int64') # case,point
    
        # find bmu indices for iseg
        in_pt, in_t = np.where(xds_bmu.bmu.values==iseg)
    
        # get indices of casei
        in_t_ = in_t - ix_time_shy[0]
        
        # reorder spectral directions
        base=5
        efth_case = seg_sim.isel(point=shy_inear)#.isel(point=in_pt, time=in_t_)
        if sign < 0: 
            efth_case['direction'] = 360 - seg_sim.direction.values
            new_dirs = np.round(efth_case.direction.values + \
                                base * round(df_icase.aseg/base) + 90, 1)
        else:
            new_dirs = np.round(efth_case.direction.values + \
                                base * round(df_icase.aseg/base) - 90, 1)  # TODO CHECK WITH KATRINA
        new_dirs = np.mod(new_dirs, 360)
        new_dirs[new_dirs > 270] -= 360
        efth_case['direction'] = new_dirs
        efth_case = efth_case.sel(direction=seg_sim.direction.values)
    
        # insert spectral values for bmu=iseg
        efth_arr[:,:,in_pt,in_t] = efth_case.efth.values[:,:,in_pt, in_t_]
    
#        sys.stdout.flush()
#        sys.stdout.write('\rInserting spectra segment: {0} / {1}\r'.format(
#                 iseg+1, xds_bmu.case.size))
##        print('segment', iseg, 'bmu spectra inserted.', datetime.datetime.now())
#
#    print('Inserting spectra segment: ', iseg+1, '/', xds_bmu.case.size)
#    print('loop end.....', datetime.datetime.now())
#    print('Inserted envelope spectra.')
    if text_out:  print('Inserting envelope spectra...', datetime.datetime.now())
        
    # store dataset
    xds_spec = xr.Dataset({'efth': (('freq','dir','point','time'), efth_arr),
                          'lon': (('point'), np.array(cp_lon_ls)), 
                          'lat': (('point'), np.array(cp_lat_ls)), 
                          'station': (('point'), np.array(cp_names)),
                         },
                        coords={
                            'freq': seg_sim.frequency.values,
                            'dir': seg_sim.direction.values,
                            'point': xds_bmu.point.values,
                            'time': xds_bmu.time.values,
                        })
    return xds_spec, xds_bmu

def get_cp_spectra(path_library, df_analogue, path, df_seg, path_mda,
                   st, cp_lon_ls, cp_lat_ls, max_dist=60, list_out=False):
    '''
    Function to access the library analogue cases for a given storm track,
    calculate distance and angle from the target segment origin to the control
    point (relative coordinate system), and extract the directional wave 
    spectra at the closest node (for every analogue segment)
    
    path_library    - (path) library of prerun segments
    path            - (path) project directory
    path_mda        - (path) mda sample indices (grid sizes)
    
    df_analogue     - (pandas.DataFrame) analogue prerun segments from library
    df_seg          - (pandas.DataFrame) storm 6h-segments parameters 
    st              - (pandas.DataFrame) storm track interpolated every 6h
    cp_lon/lat_ls   - (list) control point geographical coordinates
    max_dist        - (float) maximum distance [km] to extract closest node
    
    returns:    (xarray.Dataset) wave directional spectra (dim 'case')
    '''

    # folder path for project ensemble
    p_ensemble = op.join(path, 'ensemble')
    if not op.isdir(p_ensemble):  os.mkdir(p_ensemble)
    
    # remove NaN
    df_seg = df_seg[~df_seg.isna().any(axis=1)]
    df_analogue = df_analogue[~df_analogue.isna().any(axis=1)]
    
    # assign time
    df_seg['time'] = df_seg.index.values
    
    # get hemisphere
    sign = np.sign(df_seg.lat.values[0])

    xds_list = []
    for iseg in range(df_seg.shape[0]):
        
        # get storm segment 'i'
        df_icase = df_seg.iloc[iseg]
        df_ianalogue = df_analogue.iloc[iseg]
        iseg_analogue = df_ianalogue.name       # analogue id
        aseg = df_ianalogue.aseg                # analogue warmseg azimuth (N)
    
        # ---------------------------------------------------------------------
        # get storm coordinates at "seg_time"
        seg_time = np.datetime64(df_icase.time)  # timestamp to datetime64
        ind_time_st = np.where(st.index.values == seg_time)[0][0]
        
        # storm coordinates (real storm eye)
        st_lon = st.iloc[ind_time_st].lon       # target origin longitude
        st_lat = st.iloc[ind_time_st].lat       # target origin latitude
        st_time = st.index.values[ind_time_st]  # time coordinates
    
        # ---------------------------------------------------------------------
        # get cp coordinates in relative system (radii, angle)
        cp_dist_ls, cp_ang_ls = get_cp_radii_angle(st_lat, st_lon, 
                                                    cp_lat_ls, cp_lon_ls, 
                                                    sign, aseg)

        # get SWAN output mask indices (rad, ang)
        mask_rad, mask_ang = get_mask_radii_angle(path_mda, iseg_analogue)

        # find closest index
        ix_near = find_nearest(cp_dist_ls, cp_ang_ls, mask_rad, mask_ang)
        
        # ---------------------------------------------------------------------
        # get analogue segment from library
        filename = 'spec_outpts_main.nc'
        p_analogue = op.join(path_library, 
                             '{0:04d}'.format(iseg_analogue), filename)
        # load file
        seg_sim = xr.open_dataset(p_analogue)
        
        # time array
        hour_intervals = seg_sim.time.size
        time = [st_time + np.timedelta64(1,'h')*i for i in range(hour_intervals)]
        time_array = np.array(time)

        # extract spectra for closest index
        spec_arr = seg_sim.efth.values[:,:,ix_near,:]   # dims:freq,dir,point,time
        spec_arr = np.expand_dims(spec_arr, axis=0)
        
        # assign NaN when closest node is more distant than "max_dist"
        pos_nan = np.where(np.abs(mask_rad[ix_near] - cp_dist_ls) > max_dist)[0]
        spec_arr[:,:,:,pos_nan,:] = np.nan
                
        # store dataset
        xds = xr.Dataset({'efth': (('case','freq','dir','point','time'), spec_arr),
                          'lon': (('point'), cp_lon_ls),
                          'lat': (('point'), cp_lat_ls)},
                          coords={
                                  'case': [iseg],
                                  'freq': seg_sim.frequency.values[:],
                                  'dir': seg_sim.direction.values[:],
                                  'time': time_array,
                                  })
        xds_list.append(xds)
 
        print(iseg, datetime.datetime.now(), pos_nan, 'nan',
              [int(ele) for ele in np.abs(mask_rad[ix_near] - cp_dist_ls)],
              [int(ele) for ele in cp_dist_ls])
        
        
#        xds = xr.Dataset({'efth': (('case','freq','dir','point','time'), spec_arr),
#                          'lon': (('point'), cp_lon_ls),
#                          'lat': (('point'), cp_lat_ls)},
#                          coords={
#                                  'case': [iseg],
#                                  'freq': seg_sim.frequency.values[:],
#                                  'dir': seg_sim.direction.values[:],
#                                  'time': time_array,
#                                  })
#        xds_list.append(xds)
# 
#        print(iseg, datetime.datetime.now(), np.abs(mask_rad[ix_near] - cp_dist_ls),
#              cp_dist_ls)#, cp_ang_ls, mask_rad[ix_near], mask_ang[ix_near])

#        # extract spectra for closest index
#        spec_arr = seg_sim.efth.values[:,:,ix_near,:]   # dims:freq,dir,point,time
#        spec_arr = np.expand_dims(spec_arr, axis=0)
#        for ipt in range(ix_near.size):
#            if mask_rad[ix_near[ipt]] > max_dist: 
#                spec_arr[:,:,ipt,:] *= np.nan
#        
#        xds = xr.Dataset({'efth': (('case','freq','dir','point','time'), spec_arr),
#                          'lon': (('point'), cp_lon_ls),
#                          'lat': (('point'), cp_lat_ls)},
#                          coords={
#                                  'case': [iseg],
#                                  'freq': seg_sim.frequency.values[:],
#                                  'dir': seg_sim.direction.values[:],
#                                  'time': time_array,
#                                  })
#        xds_list.append(xds)
# 
#        print(iseg, datetime.datetime.now(), 
#              cp_dist_ls, cp_ang_ls, mask_rad[ix_near], mask_ang[ix_near])

    if list_out:    xds_out = xds_list
    else:       
        # merge
        xds_out = xr.merge(xds_list)
        print('merged', datetime.datetime.now())
    
    return xds_out

def find_nearest(cp_rad_ls, cp_ang_ls, mask_rad, mask_ang):
    '''
    Finds the minimum distance in a n-dimensional normalized space.
    
    cp_rad_ls, cp_ang_ls    - (list) control point
    mask_rad, mask_and      - (array) SWAN output point dimension
    '''
    
    # create dataframes
    df_cp = pd.DataFrame({'radii': cp_rad_ls,
                          'angle': cp_ang_ls,
                          })
                
    df_mask = pd.DataFrame({'radii': mask_rad,
                            'angle': mask_ang,
                            })
    # indices
    ix_scalar = [0]
    ix_directional = [1]
    
    # get indices of nearest n-dimensional point
    ix_near = nearest_indexes(df_cp.values, df_mask.values, 
                              ix_scalar, ix_directional)
    
    return ix_near


def get_cp_radii_angle(st_lat, st_lon, cp_lat_ls, cp_lon_ls, sign, aseg):
    '''
    Extracts the distance and angle of the control point in the relative 
    coordinate system (analogue segment)
    
    st_lat, st_lon  - (float) storm coordinates
    cp_lat, cp_lon  - (list) control point coordinates 
    sign            - (1) north hemisphere / (-1) south hemisphere
    aseg            - azimuth of the analogue warm segment (geographic north)
    '''

    # earth radius [km]
    RE = 6378.135
    
    cp_dist_ls, cp_ang_ls = [], []
    for i in range(len(cp_lat_ls)):
        
        cp_lat, cp_lon = cp_lat_ls[i], cp_lon_ls[i]

        # get point polar reference
        # azimut is refered to geographical north (absolute system)
        arcl_h, ang_abs = gc_distance(st_lat, st_lon, cp_lat, cp_lon)
        cp_dist_ls.append(arcl_h * np.pi / 180.0 * RE)   # [km]
    
        # change of coordinate system (absolute to relative)
        ang_rel = ang_abs - (aseg - 90)
        if ang_rel < 0:    ang_rel = np.mod(ang_rel, 360)
        
        # south hemisphere effect
        if sign == -1:
            if (ang_rel>=0) and (ang_rel <=180):    ang_rel = 180 - ang_rel
            elif (ang_rel>=180) and (ang_rel<=360): ang_rel = 360 - (ang_rel-180)
        
        cp_ang_ls.append(ang_rel)

    return cp_dist_ls, cp_ang_ls   # [km], [º]
        
def get_mask_radii_angle(path_mda, icase):
    '''
    Extracts the indices for radii/angle at the output points dimension
    
    path_mda    - directory
    icase       - analogue case id
    '''
               
    # load output indices associated with distances/angles to target origin
    p_mask = op.join(path_mda, 'mda_mask_indices.nc')
    if not op.isfile(p_mask):  
        xds_mask_ind = get_mask_indices(path_mda)
    else:
        xds_mask_ind = xr.open_dataset(p_mask)

    # load MDA indices (grid sizes)
#    xds_ind_mda = xr.open_dataset(op.join(path_mda, 'shytcwaves_mda_indices.nc'))
    xds_ind_mda = xr.open_dataset(op.join(path_mda, 'shytcwaves_mda_indices_clean.nc'))
    
    # get grid code
    pos_small = np.where(icase == xds_ind_mda.indices_small)[0]
    pos_medium = np.where(icase == xds_ind_mda.indices_medium)[0]
    pos_large = np.where(icase == xds_ind_mda.indices_large)[0]
    
    if len(pos_small)==1:   
        rad = xds_mask_ind.radii_sma.values[:] / 1000
        ang = xds_mask_ind.angle_sma.values[:]
        
    elif len(pos_medium)==1:
        rad = xds_mask_ind.radii_med.values[:] / 1000
        ang = xds_mask_ind.angle_med.values[:]
        
    elif len(pos_large)==1:
        rad = xds_mask_ind.radii_lar.values[:] / 1000
        ang = xds_mask_ind.angle_lar.values[:]

    return rad, ang   # [km], [º]

def get_mask_indices(path_mda):
    
    # load
#    xds_mda = xr.open_dataset(op.join(path_mda, 'shytcwaves_mda_indices.nc'))
    xds_mda = xr.open_dataset(op.join(path_mda, 'shytcwaves_mda_indices_clean.nc'))
    
    # get one case for each grid size
    case_sma = xds_mda.indices_small.values[0]
    case_med = xds_mda.indices_medium.values[0]
    case_lar = xds_mda.indices_large.values[0]

    rr, aa = [], []
    # get analogue case
    for ic,case_id in enumerate([case_sma, case_med, case_lar]):

        # get polar coords (radii/distance, angle)
        rc, thetac = find_analogue_grid_coords(path_mda, case_id)

        # to point dimension
        rr.append(np.reshape(rc, -1))       # [m]
        aa.append(np.reshape(thetac, -1))   # [º]

    # store all grid sizes mask indices
    xd = xr.Dataset({
        'radii_sma': (('point_sma',), rr[0]),
        'angle_sma': (('point_sma',), aa[0]),
        'radii_med': (('point_med',), rr[1]),
        'angle_med': (('point_med',), aa[1]),
        'radii_lar': (('point_lar',), rr[2]),
        'angle_lar': (('point_lar',), aa[2]),
    })
    
    # store dataset
    xd.to_netcdf(op.join(path_mda, 'mda_mask_indices.nc'))

    return xd


def aux_rec_wavespectra(p_analogue, path_library, iseg, spec_hs, spec_tp):

    # load file
    seg_sim = xr.open_dataset(p_analogue)
    seg_sim = seg_sim.rename(frequency='freq', direction='dir')
    spec_time = seg_sim.time.values
    
    # bulk hs
    sp_hs = seg_sim.efth.spec.hs()
    size_pts = sp_hs.shape[0]
    spec_hs[iseg,:size_pts,:] = sp_hs[:,:]
    
    # bulk tp
    sp_tp = seg_sim.efth.spec.tp()
    spec_tp[iseg,:size_pts,:] = sp_tp[:,:]
    
    print('case', iseg, 'wavespectra reconstruction done. ', 
          datetime.datetime.now(), spec_hs.shape, 'gridsize', sp_hs.shape[0])
    
    gc.collect()   # free wavespectra memory

    return spec_hs, spec_tp, spec_time
    
def save_library_bulk_hs_tp(path_library, path_mda, analogue_ls, 
                            num_cases=5000, num_time=48):
    '''
    Access the library analogue cases (or list), to calculate and store 
    the reconstructed bulk hs
    
    path_library    - (path) library of prerun segments
    path_mda        - (path) mda sample indices (grid sizes)
    analogue_ls     - (list) cases selected for reconstruction
    
    num_cases       - (float) total number of prerun cases
    num_time        - (float) total number of time steps [hour]
    
    returns:    (xarray.Dataset) reconstructed wave directional spectra
    '''
    
    # get grid largest size
    xds_mda_ind = xr.open_dataset(op.join(path_mda, 'mda_mask_indices.nc'))
    num_large = xds_mda_ind.point_lar.size   # largest size
        
    # variables for reconstruction (case,point,time)
    spec_hs = np.full((num_cases, num_large, num_time), np.nan)
    spec_tp = np.full((num_cases, num_large, num_time), np.nan)
    
    print('starting loop...', datetime.datetime.now())
    
    for iseg in np.unique(analogue_ls):     # each analogue case
        
        # get analogue segment from library
        filename = 'spec_outpts_main.nc'
        p_analogue = op.join(path_library, 
                             '{0:04d}'.format(iseg), filename)
        
        # wavespectra reconstruction
        spec_hs, spec_tp, time = aux_rec_wavespectra(p_analogue, path_library,
                                                     iseg, spec_hs, spec_tp)
        
    # store dataset
    xd = xr.Dataset({
        'hs': (('case','point','time'), spec_hs),
        'tp': (('case','point','time'), spec_tp),
    }, coords={
        'case': np.arange(0, num_cases),
        'point': np.arange(0, num_large),
        'time': time, #seg_sim.time.values,
    })
    
    # save at library path
    xd.to_netcdf(op.join(path_library, 'library_shytcwaves_bulk_params.nc'))
    
    print('Bulk variables file stored.')
 
    
###############################################################################
# STOPMOTION forecast predicted tracks
###############################################################################

from .storms import track_triming
from .superpoint import SuperPoint_Superposition

def TCforecast_process(xds, center='WMO', ipred=None):
    '''
    center    - (string) 'WMO' uses 1-min avg wind data
    ipred     - (float/None) prediction index
    '''
    
    xds = xds.rename({'longitude':'lon', 
                      'latitude':'lat', 
                      'pressure':'wmo_pres', 
                      'maxwinds':'wmo_wind', 
#                       'Basin':'basin', 
                     })
    
    # select predicted track / real track is given
    if type(ipred)==int:  xds = xds.isel(track=ipred)
    
    # processing
    xds = xds.rename({'time': 'time_real'})
    xds = xds.rename_dims({'position':'time'})
#     if ipred:  xds = xds.rename_dims({'position':'time'})    # Laura: unificar label
#     else:      xds = xds.rename_dims({'index':'time'})
        
    xds = xds.assign_coords({"time": xds.time_real.values})
    if '' not in xds.basin.values:  basin_id = np.unique(xds.basin.values)[0]
    else:                           basin_id = np.unique(xds.basin.values)[1]
        

    if type(ipred)==int:  
        xds = xds.isel(time=np.where(np.isnat(xds['time']) == False)[0]
                                ).resample(time='6H').interpolate('linear')
    else: 
        _, index = np.unique(xds['time'], return_index=True)
        xds = xds.isel(time=index).resample(time='6H').interpolate('linear')

    xds['basin'] = (('time'), [basin_id]*xds.time.size)
    xds['center'] = (('time'), [center]*xds.time.size)
    xds['dist2land'] = (('time'), [100]*xds.time.size)

    return xds

def TCforecast2shytcwaves(p_save, tc_name, iprediction, storm, center, lon, lat, 
                          dict_site={}):
    '''
    Function that for a storm track and a target domain, provides the shytcwaves 
    estimated induced spectral wave energy over time
    
    p_save     - (path) store for results
    tc_name    - (str) name
    iprediction- (float) predicted number
    storm      - (xarray.Dataset) with standard variables extracted from IBTrACS
                 compulsory: 'longitude,latitude', pressure[mbar], maxwinds[kt]
                 optional: rmw[nmile], dist2land[], basin[str]
    center     - (str) IBTrACS center track data (many are available)
    lon,lat    - (array) longitude/latitude nodes for calculating swath
    dict_site  - (dict) site data for superpoint building
    '''

    # A: stopmotion segmentation, 6h interval 
    df = historic_track_preprocessing(storm, center=center, forecast_on=True)

    dt_int_seg = 6*60       # [minutes] constant segments
    st, _ = historic_track_interpolation(df, dt_int_seg,
                                         interpolation=False, mode='mean')

    # skip when only NaN or 0
    lons, lats = st.lon.values, st.lat.values
    if (np.unique(lons[~np.isnan(lons)]).all() == 0) & \
        (np.unique(lats[~np.isnan(lats)]).all() == 0): 
        print('Forecast: {0}, No predicted track coordinates'.format(iprediction))

    else: 
        st_trim = track_triming(st, lat[0], lon[0], lat[-1], lon[-1])
        
        # store tracks for shytcwaves
        st.to_pickle(op.join(p_save, '{0}_forecast_{1}_track.pkl'.format(
                tc_name, iprediction)))
        st_trim.to_pickle(op.join(p_save, '{0}_forecast_{1}_track_trim.pkl'.format(
                tc_name, iprediction)))
    
        # parameterized segts (24h warmup + 6htarget)
        df_seg = storm2stopmotion(st_trim)   
        
        if df_seg.shape[0] > 2:
    
            print('st:', st.shape[0], 'df_seg:', df_seg.shape[0])
            st_list, we_list = stopmotion_interpolation(df_seg, st=st_trim)
    
            #######################################################################
            # B: analogue segments from library
    #        df_mda = pd.read_pickle(r'/media/administrador/HD8/BlueMath/data/stopmotion/shytcwaves_mda.pkl')
            df_mda = pd.read_pickle(r'/media/administrador/HD8/BlueMath/data/stopmotion/shytcwaves_mda_clean.pkl')
            ix_weights = [1]*10   # equal weights
            ix = find_analogue(df_mda, df_seg, ix_weights)
    
            df_analogue = df_mda.iloc[ix]
            df_analogue = analogue_endpoints(df_seg, df_analogue)
    
            st_list_analogue, we_list_analogue = stopmotion_interpolation(
                    df_analogue, st=st_trim)
    
            #######################################################################
            # C: extract bulk envelope (to plot swaths)
            mesh_lo, mesh_la = np.meshgrid(lon, lat)
            print('Forecast: {0}, number of segments: {1}, number of swath nodes: {2}'.format(
                iprediction, len(st_list), mesh_lo.size))
    #        print('Extracting bulk envelope...')
    
            p_library = r'/media/administrador/HD8/conf_drm4g/output'
            p_stopmotion = r'/media/administrador/HD8/BlueMath/data/stopmotion'
            xds_shy_bulk = stopmotion_st_bmu(p_library, df_analogue, 
                                             df_seg, p_stopmotion, st_trim,
                                             list(np.ravel(mesh_lo)), 
                                             list(np.ravel(mesh_la)), max_dist=60)
            # store
            xds_shy_bulk.to_netcdf(op.join(p_save, 
                '{0}_forecast_{1}_xds_shy_bulk.nc'.format(tc_name, iprediction)))
    
            #######################################################################
            # D: extract spectra envelope
            if type(dict_site)==dict:
    #            print('Extracting spectra envelope...')
                xds_shy_spec, _ = stopmotion_st_spectra(p_library, df_analogue, 
                                                        df_seg, p_stopmotion, 
                                                        st_trim, dict_site['lonpts'], 
                                                        dict_site['latpts'], 
                                                        cp_names=dict_site['namepts'], 
                                                        max_dist=60, list_out=False)
                # store
                xds_shy_spec.to_netcdf(op.join(p_save, 
                    '{0}_forecast_{1}_xds_shy_spec_{2}.nc'.format(
                            tc_name, iprediction, dict_site['site'])))
    
                ###################################################################
                # E: build superpoint
                stations = list(np.arange(0, xds_shy_spec.point.size))
                xds_shy_sp = SuperPoint_Superposition(xds_shy_spec, stations, 
                                                      dict_site['sectors'], 
                                                      dict_site['deg_superposition'])
                # store
                xds_shy_sp.to_netcdf(op.join(p_save, 
                    '{0}_forecast_{1}_xds_shy_sp_{2}.nc'.format(
                            tc_name, iprediction, dict_site['site'])))
                        
        print('Predicted forecast files stored.\n')

def SHy_forecast_probs(p_ensemble, tc_target, 
                       ntracks=10, threshold=[2,5,8,12]):
    
    # get a single bulk file
    p_bulk = op.join(p_ensemble, '{0}_forecast_{1}_xds_shy_bulk.nc'.format(
            tc_target, 2))
    # TODO: get the first file (not empty)
    xds_shy_bulk = xr.open_dataset(p_bulk)
    
    # get predicted swaths and times
    hspred = np.full((ntracks, xds_shy_bulk.hswath.size), np.nan)
    tpred = np.full((ntracks, xds_shy_bulk.point.size), np.nan, 
                    dtype='datetime64[h]')
    i_, varmax = [], []

    for ipred in range(ntracks):
        p_bulk = op.join(p_ensemble, '{0}_forecast_{1}_xds_shy_bulk.nc'.format(
                tc_target, ipred))
        
        if op.isfile(p_bulk):  
            xds_shy_bulk = xr.open_dataset(p_bulk)
            hspred[ipred,:] = xds_shy_bulk.hswath.values
            varmax.append(np.nanmax(xds_shy_bulk.hswath.values))
            i_.append(ipred)
            
            hsbmu = xds_shy_bulk.hsbmu.values
            hsbmu[np.isnan(hsbmu)] = -1           # replace nan
            posbmu = np.nanargmax(hsbmu, axis=1)  # max indices
            valbmu = np.nanmax(hsbmu, axis=1)     # max values
            posmask = np.where(valbmu>=0)[0]      # mask to add nans
            # valbmu[valbmu==-1] = np.nan
            tpred[ipred,posmask] = xds_shy_bulk.time.values.astype(
                    'datetime64[h]')[posbmu[posmask]]
            
    # remove empty rows
    hspred = hspred[i_,:]
    tpred = tpred[i_,:]
    varmax = np.nanmax(varmax)
    dtpred = (tpred - np.nanmin(tpred)) / np.timedelta64(1, 'D')   # days
        
    # probabilities for plotting
    pred_mean = np.nanmean(hspred, axis=0)
    pred_pct95 = np.nanpercentile(hspred, 95, axis=0)
    pred_pct99 = np.nanpercentile(hspred, 99, axis=0)
    pred_max = np.nanmax(hspred, axis=0)
    pred_varmax = np.nanmax(pred_max)
    pred_occur = np.nansum(~np.isnan(hspred), axis=0) / ntracks *100
    pred_prob = []
    for i in range(len(threshold)):
        pred_prob.append(np.nansum(hspred > threshold[i], axis=0) \
                         / np.nansum(hspred >= 0, axis=0) * 100)
        
#    pred_2 = np.nansum(hspred >5, axis=0) / np.nansum(hspred>=0, axis=0)*100
#    pred_3 = np.nansum(hspred >8, axis=0) / np.nansum(hspred>=0, axis=0)*100
#    pred_4 = np.nansum(hspred >12, axis=0) / np.nansum(hspred>=0, axis=0)*100
#    pred_2m = np.nansum(hspred > 2, axis=0) / np.nansum(hspred >= 0, axis=0)*100
#    pred_5m = np.nansum(hspred > 5, axis=0) / np.nansum(hspred >= 0, axis=0)*100
#    pred_8m = np.nansum(hspred > 8, axis=0) / np.nansum(hspred >= 0, axis=0)*100
#    pred_12m = np.nansum(hspred > 12, axis=0) / np.nansum(hspred >= 0, axis=0)*100
    
    # times for contour plotting
    dtpred_mean = np.nanmean(dtpred, axis=0)
    dt_prob = []
    for i in range(len(threshold)):
        dtpred_ = dtpred.copy()
        dtpred_[hspred <= threshold[i]] = np.nan
        dt_prob.append(np.nanmean(dtpred_, axis=0))
        
#    dtpred_2m, dtpred_5m, dtpred_8m, dtpred_12m = dtpred.copy(), dtpred.copy(), dtpred.copy(), dtpred.copy()
#    dtpred_2m[hspred <= 2] = np.nan
#    dtpred_5m[hspred <= 5] = np.nan
#    dtpred_8m[hspred <= 8] = np.nan
#    dtpred_12m[hspred <= 12] = np.nan
#    dtpred_2m = np.nanmean(dtpred_2m, axis=0)
#    dtpred_5m = np.nanmean(dtpred_5m, axis=0)
#    dtpred_8m = np.nanmean(dtpred_8m, axis=0)
#    dtpred_12m = np.nanmean(dtpred_12m, axis=0)
    
    var_ls = [None, pred_mean, pred_pct95, pred_pct99, pred_max, 
              pred_occur, pred_prob[0], pred_prob[1], pred_prob[2], pred_prob[3]]
    ttl_ls = ['JTWC predicted tracks', 'Mean swath', 'Pctl95 swath', 
              'Pctl99 swath', 'Max swath', 'Occurence probability', 
              'Prob Hs>{0}m'.format(threshold[0]), 
              'Prob Hs>{0}m'.format(threshold[1]), 
              'Prob Hs>{0}m'.format(threshold[2]), 
              'Prob Hs>{0}m'.format(threshold[3])]
#    ttl_ls = ['JTWC predicted tracks', 'Mean swath [m]', 'Pctl95 swath [m]', 
#              'Pctl99 swath [m]', 'Max swath [m]', 'Prediction probability [%]', 
#              'Prob >2m [%]', 'Prob >5m [%]', 'Prob >8m [%]', 'Prob >12m [%]']
    dt_ls = [None, dtpred_mean, None, None, None, 
             None, dt_prob[0], dt_prob[1], dt_prob[2], dt_prob[3]]
#             None, dtpred_2m, dtpred_5m, dtpred_8m, dtpred_12m]
    
    dict_tc = {'var_ls': var_ls,
               'ttl_ls': ttl_ls,
               'dt_ls': dt_ls, 
               'xds_shy_bulk': xds_shy_bulk,
#               'hspred': hspred,
               'pred_varmax': pred_varmax,
               }
    
    return dict_tc

###############################################################################
# STOPMOTION tesla emulator tracks
###############################################################################

def TCemulator2shytcwaves(p_library, path_mda, p_save, storm, 
                          loc, lac, radii, num=0, name='historic',
                          center='WMO', dict_site={}):
    '''
    Function that for a storm track and a target domain, provides the shytcwaves 
    estimated induced spectral wave energy over time
    
    p_save     - (path) store for results, site folder
    storm      - (xarray.Dataset) with standard variables extracted from IBTrACS
                 compulsory: 'longitude,latitude', pressure[mbar], maxwinds[kt]
                 optional: rmw[nmile], dist2land[], basin[str]
    lon,lat    - (array) longitude/latitude nodes for calculating swath
    center     - (str) IBTrACS center track data (many are available)
    dict_site  - (dict) site data for superpoint building
    '''

    if not type(dict_site)==dict:  print('dict_site must be a dictionary')
        
    # A: stopmotion segmentation, 6h interval 
    df = historic_track_preprocessing(storm, center=center)

    dt_int_seg = 6*60       # [minutes] constant segments
    st, _ = historic_track_interpolation(df, dt_int_seg,
                                         interpolation=False, mode='mean')
    
    # parameterized segments (24h warmup + 6htarget)
    df_seg = storm2stopmotion(st)   
    
    if not df_seg.dropna().empty:   # not empty

        # segments/storm trim within domain (optional)
        df_trim = stopmotion_trim_circle(df_seg, loc, lac, radii)
        
        if not df_trim.empty:
            if df_trim.shape[0] > 3:
                
                st_trim = st.iloc[(st.index>=np.min(df_trim.index)) & \
                                  (st.index<=np.max(df_trim.index))]
#                 st_trim = track_triming_circle(st, loc, lac, radii)
                st_cat = np.max(get_category(st_trim.p0))
                st_id = int(storm.storm.values)

                # generate stopmotion segments
                print(num, datetime.datetime.now(), '/ id:', st_id, '/ cat:', 
                      st_cat, '/ df_trim:', df_trim.shape[0])

                st_list, we_list = stopmotion_interpolation(df_trim, st=st_trim)

                # store
                st_trim.to_pickle(op.join(p_save, 
                                          'cat{0}/track_{1}_cat{0}_id{2}_{3}.pkl'.format(
                                              st_cat, dict_site['site'], st_id, name)))
#                 print(df_trim)
#                 print(st_trim)

                #######################################################################
                # B: analogue segments from library
#                df_mda = pd.read_pickle(r'/media/administrador/HD8/BlueMath/data/stopmotion/shytcwaves_mda.pkl')
                df_mda = pd.read_pickle(path_mda)
                ix_weights = [1]*10   # equal weights
                ix = find_analogue(df_mda, df_trim, ix_weights)

                df_analogue = df_mda.iloc[ix]
                df_analogue = analogue_endpoints(df_trim, df_analogue)

                st_list_analogue, we_list_analogue = stopmotion_interpolation(
                        df_analogue, st=st_trim)

                #######################################################################
                # D: extract spectra envelope
#                p_library = r'/media/administrador/HD8/conf_drm4g/output'
#                p_stopmotion = r'/media/administrador/HD8/BlueMath/data/stopmotion'

        #            print('Extracting spectra envelope...')
                xds_shy_spec, _ = stopmotion_st_spectra(p_library, df_analogue, 
                                                        df_trim, path_mda, 
                                                        st_trim, dict_site['lonpts'], 
                                                        dict_site['latpts'], 
                                                        cp_names=dict_site['namepts'], 
                                                        max_dist=60, list_out=False, 
                                                        tqdm_out=False, text_out=False)
                # store
                xds_shy_spec.to_netcdf(op.join(p_save, 
                                               'cat{0}/spectra_{1}_cat{0}_id{2}_{3}.nc'.format(
                                                   st_cat, dict_site['site'], st_id, name)))

                ###################################################################
                # E: build superpoint
                stations = list(np.arange(0, xds_shy_spec.point.size))
                xds_shy_sp = SuperPoint_Superposition(xds_shy_spec, stations, 
                                                      dict_site['sectors'], 
                                                      dict_site['deg_superposition'])

                # store
                xds_shy_sp.to_netcdf(op.join(p_save, 
                                             'cat{0}/superpoint_{1}_cat{0}_id{2}_{3}.nc'.format(
                                                 st_cat, dict_site['site'], st_id, name)))
        
    return st_id

def TCforecast2shytcwaves_cluster(p_save, path_mda, p_library, p_stopmotion,
                                  tc_name, iprediction, 
                                  storm, center, lon, lat, dict_site={}):
    '''
    Function that for a storm track and a target domain, provides the shytcwaves 
    estimated induced spectral wave energy over time
    
    p_save     - (path) store for results
    tc_name    - (str) name
    iprediction- (float) predicted number
    storm      - (xarray.Dataset) with standard variables extracted from IBTrACS
                 compulsory: 'longitude,latitude', pressure[mbar], maxwinds[kt]
                 optional: rmw[nmile], dist2land[], basin[str]
    center     - (str) IBTrACS center track data (many are available)
    lon,lat    - (array) longitude/latitude nodes for calculating swath
    dict_site  - (dict) site data for superpoint building
    '''

    # A: stopmotion segmentation, 6h interval 
    df = historic_track_preprocessing(storm, center=center, forecast_on=True)

    dt_int_seg = 6*60       # [minutes] constant segments
    st, _ = historic_track_interpolation(df, dt_int_seg,
                                         interpolation=False, mode='mean')

    # skip when only NaN or 0
    lons, lats = st.lon.values, st.lat.values
    if (np.unique(lons[~np.isnan(lons)]).all() == 0) & \
        (np.unique(lats[~np.isnan(lats)]).all() == 0): 
    
        print('Forecast: {0}, No predicted track coordinates'.format(iprediction))

    else: 
        st_trim = track_triming(st, lat[0], lon[0], lat[-1], lon[-1])
        
        # store tracks for shytcwaves
        st.to_pickle(op.join(p_save, '{0}_forecast_{1}_track.pkl'.format(
                tc_name, iprediction)))
        st_trim.to_pickle(op.join(p_save, '{0}_forecast_{1}_track_trim.pkl'.format(
                tc_name, iprediction)))
    
        # parameterized segts (24h warmup + 6htarget)
        df_seg = storm2stopmotion(st_trim)   
        
        if df_seg.shape[0] > 2:
    
            print('st:', st.shape[0], 'df_seg:', df_seg.shape[0])
            st_list, we_list = stopmotion_interpolation(df_seg, st=st_trim)
    
            #######################################################################
            # B: analogue segments from library
    #        df_mda = pd.read_pickle(r'/media/administrador/HD8/BlueMath/data/stopmotion/shytcwaves_mda.pkl')
            # df_mda = pd.read_pickle(r'/media/administrador/HD8/BlueMath/data/stopmotion/shytcwaves_mda_clean.pkl')
            df_mda = pd.read_pickle(op.join(path_mda, 'shytcwaves_mda_clean.pkl'))
            ix_weights = [1]*10   # equal weights
            ix = find_analogue(df_mda, df_seg, ix_weights)
    
            df_analogue = df_mda.iloc[ix]
            df_analogue = analogue_endpoints(df_seg, df_analogue)
    
            st_list_analogue, we_list_analogue = stopmotion_interpolation(
                    df_analogue, st=st_trim)
    
            #######################################################################
            # C: extract bulk envelope (to plot swaths)
            mesh_lo, mesh_la = np.meshgrid(lon, lat)
            print('Forecast: {0}, number of segments: {1}, number of swath nodes: {2}'.format(
                iprediction, len(st_list), mesh_lo.size))
    #        print('Extracting bulk envelope...')
    
            # p_library = r'/media/administrador/HD8/conf_drm4g/output'
            # p_stopmotion = r'/media/administrador/HD8/BlueMath/data/stopmotion'
            xds_shy_bulk = stopmotion_st_bmu(p_library, df_analogue, 
                                             df_seg, p_stopmotion, st_trim,
                                             list(np.ravel(mesh_lo)), 
                                             list(np.ravel(mesh_la)), max_dist=60)
            # store
            xds_shy_bulk.to_netcdf(op.join(p_save, 
                '{0}_forecast_{1}_xds_shy_bulk.nc'.format(tc_name, iprediction)))
    
            #######################################################################
            # D: extract spectra envelope
            if type(dict_site)==dict:
    #            print('Extracting spectra envelope...')
                xds_shy_spec, _ = stopmotion_st_spectra(p_library, df_analogue, 
                                                        df_seg, p_stopmotion, 
                                                        st_trim, dict_site['lonpts'], 
                                                        dict_site['latpts'], 
                                                        cp_names=dict_site['namepts'], 
                                                        max_dist=60, list_out=False)
                # store
                xds_shy_spec.to_netcdf(op.join(p_save, 
                    '{0}_forecast_{1}_xds_shy_spec_{2}.nc'.format(
                            tc_name, iprediction, dict_site['site'])))
    
                ###################################################################
                # E: build superpoint
                stations = list(np.arange(0, xds_shy_spec.point.size))
                xds_shy_sp = SuperPoint_Superposition(xds_shy_spec, stations, 
                                                      dict_site['sectors'], 
                                                      dict_site['deg_superposition'])
                # store
                xds_shy_sp.to_netcdf(op.join(p_save, 
                    '{0}_forecast_{1}_xds_shy_sp_{2}.nc'.format(
                            tc_name, iprediction, dict_site['site'])))
            
    print('Predicted forecast files stored.\n')

