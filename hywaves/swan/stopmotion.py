#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os.path as op
import numpy as np
import pandas as pd
import xarray as xr
from datetime import timedelta
#from math import radians, degrees, sin, cos, atan2, pi  #asin, acos, sqrt, 

from .geo import shoot, gc_distance, GeoAzimuth
from .storms import Extract_basin_storms, get_vmean, \
d_vns_basinscenter, historic_track_preprocessing, historic_track_interpolation


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
#    xxr, yyr = np.meshgrid(np.arange(-xlength/2, xlength/2, res), 
#                           np.arange(-xlength/2, xlength/2, res))

    # B: custom polar grid calculation

    # radii array
    rc = np.zeros(rings)
    rc[0] = radii_ini
    for i in range(1, rings):       rc[i] = ((i+1)/i) ** 1.5 * rc[i-1]
    
    # radii (m), angle (º) meshgrid
    rc, thetac = np.meshgrid(rc, 
                             np.arange(0, 360, angle_inc))
    
    # polar grid
    xxp = rc * np.sin(thetac * np.pi/180)
    yyp = rc * np.cos(thetac * np.pi/180)
    
    return xxr, yyr, xxp, yyp  # 2d-arrays


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
    dv  = np.full(df_.shape[0], np.nan)     # translational speed variation
    dvx = np.full(df_.shape[0], np.nan)     # (dirx)
    dvy = np.full(df_.shape[0], np.nan)     # (diry)
    da  = np.full(df_.shape[0], np.nan)     # azimuth variation
    dp  = np.full(df_.shape[0], np.nan)     # pressure variation
    dw  = np.full(df_.shape[0], np.nan)     # maxwinds variation
    dr  = np.full(df_.shape[0], np.nan)     # radii rmw variation
    dl  = np.full(df_.shape[0], np.nan)     # latitude variation
    dla = np.full(df_.shape[0], np.nan)     # (absolute)
    
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
        
        if i >= 4:  # >= four preceding segments
            
            # endpoint coordinates (-24h, 0h, 6h)
            lon_0, lon_i, lon_i1 = lon[i-4], lon[i], lon[i+1]
            lat_0, lat_i, lat_i1 = lat[i-4], lat[i], lat[i+1]

        # warmup 4-segments (24h) variables
        lo0[i], la0[i] = lon_0, lat_0
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
    df_['lon0'] = lo0
    df_['lat0'] = la0
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
                                    'REUNION','BOM','WELLINGTON','NADI','WMO']):
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
    
    returns:  df_database (pandas.Datafrate)
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
                                      'vmaxfill', 'rmwfill']]
                        
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
#    path = op.join(op.dirname(op.realpath(__file__)), '..', '..', '..', \
#                   'stopmotion', 'stopmotion_database_{0}.pkl'.format(
#                           ibtracs_version))
#    df_database.to_pickle(path)
    
    return df_database

