#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from datetime import timedelta

from .geo import shoot, gc_distance


# STORM TRACK LIBRARY

# TODO descatalogar? 
def track_from_parameters(
    pmin, vmean, delta, gamma,
    x0, y0, x1, R,
    date_ini, hours,
    great_circle=False):
    '''
    Calculates storm track variables from storm track parameters

    pmin, vmean, delta, gamma  - storm track parameters
    x0, y0                     - site coordinates (longitude, latitude)
    x1                         - enter point in computational grid
    R                          - radius (º)
    date_ini                   - initial date 'yyyy-mm-dd HH:SS'
    hours                      - number of hours to generate
    great_circle               - True for using great circle lon,lat calculation
    '''

    RE = 6378.135  # earth radius

    # generation of storm track 
    xc = x0 + R * np.sin(delta * np.pi/180)  # enter point in the smaller radius
    yc = y0 + R * np.cos(delta * np.pi/180)

    d = (x1 - xc) / np.sin(gamma * np.pi/180)
    y1 = yc + d * np.cos(gamma * np.pi/180)

    # time array for SWAN input
    time_input = pd.date_range(date_ini, periods=hours, freq='H')

    # storm track (pd.DataFrame)
    st = pd.DataFrame(index=time_input, columns=['move', 'vf', 'pn', 'p0', 'lon', 'lat'])
    st['move'] = gamma
    st['vf'] = vmean
    st['pn'] = 1013
    st['p0'] = pmin

    # calculate lon and lat
    if not great_circle:
        st['lon'] = x1 - (st['vf']*180/(RE*np.pi)) * np.sin(gamma*np.pi/180) * list(range(len(st)))
        st['lat'] = y1 - (st['vf']*180/(RE*np.pi)) * np.cos(gamma*np.pi/180) * list(range(len(st)))

    else:
        x2 = x1 - (st['vf']*180/(RE*np.pi)) * np.sin(gamma*np.pi/180) * list(range(len(st)))
        y2 = y1 - (st['vf']*180/(RE*np.pi)) * np.cos(gamma*np.pi/180) * list(range(len(st)))
        xt, yt = [], []
        for i in list(range(0,hours)):
            glon, glat, baz = shoot(x1, y1, gamma+180, vmean * i)
            xt.append(glon)
            yt.append(glat)
        st['lon'] = xt
        st['lat'] = yt

    # add some metadata
    st.x0 = x0
    st.y0 = y0
    st.R = R

    return st

def get_category(ycpres):
    'Defines storm category according to minimum pressure centers'

    categ = []
    for i in range(len(ycpres)):
        if (ycpres[i] == 0) or (np.isnan(ycpres[i])):
            categ.append(6)
        elif ycpres[i] < 920:  categ.append(5)
        elif ycpres[i] < 944:  categ.append(4)
        elif ycpres[i] < 964:  categ.append(3)
        elif ycpres[i] < 979:  categ.append(2)
        elif ycpres[i] < 1000: categ.append(1)
        elif ycpres[i] >= 1000: categ.append(0)

    return categ

def historic_track_preprocessing(xds, d_vns):
    '''
    Historic track is preprocessed, by removing NaN data, apply longitude
    sign convention, change time format and define storm category

    xds: historic track dataset (storm dimension)
    d_vns: dictionary to set longitude, latitude, time, pressure and wind varnames
    '''

    # get names of vars
    nm_lon = d_vns['longitude']
    nm_lat = d_vns['latitude']
    nm_prs = d_vns['pressure']
    nm_tim = d_vns['time']
    nm_win = d_vns['maxwinds']

    # get var time
    ytime = xds[nm_tim].values         # dates format: datetime64

    # remove time nans
    ycpres = xds[nm_prs].values[~np.isnat(ytime)]     # minimum pressure
    ylat_tc = xds[nm_lat].values[~np.isnat(ytime)]    # latitude
    ylon_tc = xds[nm_lon].values[~np.isnat(ytime)]    # longitude
    ywind = xds[nm_win].values[~np.isnat(ytime)]      # wind speed [kt]
    ytime = ytime[~np.isnat(ytime)]

    # remove pmin nans
    ylat_tc = ylat_tc[~np.isnan(ycpres)]
    ylon_tc = ylon_tc[~np.isnan(ycpres)]
    ywind = ywind[~np.isnan(ycpres)]
    ytime = ytime[~np.isnan(ycpres)]
    ycpres = ycpres[~np.isnan(ycpres)]

    # sign convention: [0º,360º]
    ylon_tc[ylon_tc<0] = ylon_tc[ylon_tc<0] + 360

    # round dates to hour
    st_time = []
    for i in range(len(ytime)):
        round_to = 3600
        # each gauge has different time format
        dt = ytime[i].astype('datetime64[s]').tolist()
        seconds = (dt - dt.min).seconds
        rounding = (seconds+round_to/2) // round_to * round_to
        out = dt + timedelta(0,rounding-seconds,-dt.microsecond)
        st_time.append(out)
    st_time = np.asarray(st_time)

    # time step data [hours]
    ts = (st_time[1:] - st_time[:-1])
    ts = [ts[i].total_seconds() / 3600 for i in range(ts.size)]

    # storm category centers
    categ = get_category(ycpres)

    return st_time, ylat_tc, ylon_tc, ycpres, ywind, ts, categ

def ibtrac_basin_fitting(x0, y0):
    '''
    Assigns cubic polynomial fitting curve coefficient for each basin of
    historical TCs data (IBTrACS)
    '''

    # determination of the location basin 
    if y0 < 0:                  basin = 5
    elif (y0 > 0) & (x0 > 0):   basin = 3
    else:                       print('Basin not defined')

    # cubic polynomial fitting curve for Ibtracs and each basin
    # TODO: obtain all basin fitting coefficients

    if basin == 3:      # West Pacific
        p1 = -7.77328602747578e-06
        p2 = 0.0190830514629838
        p3 = -15.9630945598490
        p4 = 4687.76462404360

    elif basin == 5:    # South Pacific
        p1 = -4.70481986864773e-05
        p2 = 0.131052968357409
        p3 = -122.487981649828
        p4 = 38509.7575283218

    return p1, p2, p3, p4

def historic_track_interpolation(st_time, ylon_tc, ylat_tc, ycpres, ywind, y0, x0,
                                 lat00, lon00, lat01, lon01, ts, dt_comp,
                                 wind=None, great_circle=False, fit=False):
    '''
    Calculates storm track variables from storm track parameters and interpolates
    track points in between historical data (for "dt_comp" time step)

    st_time                    - time track
    lat, lon                   - track coordinates (longitude, latitude)
    pmin, ywind                - storm track parameters
    x0, y0                     - target coordinates (longitude, latitude)
    lat0, lon0, lat1, lon1     - numerical domain bound limits
    ts                         - track time step data [hours]
    dt_comp                    - simulation computation time step [minutes]
    wind                       - False for using vmax approximation eq.
    great_circle               - True for using great circle lon,lat calculation
    fit                        - True for fitting vmax when ywind=0 (start of storm)
    '''

    RE = 6378.135   # earth radius [km]

    # cubic polynomial fitting curve for each IBTrACS basin
    p1, p2, p3, p4 = ibtrac_basin_fitting(x0, y0)

    # generate lists
    time_storm = list(st_time)  # datetime format
    pmin = list(ycpres)
    lat = list(ylat_tc)
    lon = list(ylon_tc)
    if wind.any() != None:
        mwind = wind
        if fit:
            wind_fitting = p1 * np.power(pmin,3) + p2 * np.power(pmin,2) + p3 * np.power(pmin,1) + p4
            pos = np.where(mwind==0)
            mwind[pos] = wind_fitting[pos]

    # number of time steps between consecutive interpolated track points in order
    # to  match SWAN computational time step
    ts = np.asarray(ts) * 60 / dt_comp

    # initialize
    move, vmean, pn, p0, lon_t, lat_t, vmax = [], [], [], [], [], [], []
    vu, vy = [], []
    time_input = np.empty((0,),dtype='datetime64[ns]')

    for i in range(0, len(time_storm)-1):
        # time array for SWAN input
        date_ini = time_storm[i]
        time_input0 = pd.date_range(
                date_ini, periods=int(ts[i]), freq='{0}MIN'.format(dt_comp))
        time_input = np.append(np.array(time_input), np.array(time_input0))

        # track pair of successive coordinates
        lon1 = lon[i]
        lat1 = lat[i]
        lon2 = lon[i+1]
        lat2 = lat[i+1]

        # translation speed 
        arcl_h, gamma_h = gc_distance(lat2, lon2, lat1, lon1)
        r = arcl_h * np.pi / 180.0 * RE     # distance between consecutive track points (km)
        dx = r / ts[i]                      # interpolation distance 
        vx = float(dx) /3.6                 # translation speed (m/s)
        vx = vx /0.52                       # translation speed (kt)

        for j in range(int(ts[i])):
            # append track parameters
            move.append(gamma_h)
            vmean.append(vx)
            vu.append(vx * np.sin((gamma_h+180)*np.pi/180))
            vy.append(vx * np.cos((gamma_h+180)*np.pi/180))
            pn.append(1013)
            p0.append(pmin[i] + j* (pmin[i+1]-pmin[i])/ts[i])
            if wind.any() != None:
                vmax.append(mwind[i] + j* (mwind[i+1]-mwind[i])/ts[i])    #[kt]

            # calculate lon, lat
            if not great_circle:
                lon_h = lon1 - (dx*180/(RE*np.pi)) * np.sin(gamma_h*np.pi/180) * j
                lat_h = lat1 - (dx*180/(RE*np.pi)) * np.cos(gamma_h*np.pi/180) * j
            else:
                xt, yt = [], []
                glon, glat, baz = shoot(lon1, lat1, gamma_h + 180, float(dx) * j)
                xt = np.append(xt,glon)
                yt = np.append(yt,glat)
                lon_h = xt
                lat_h = yt
            lon_t.append(lon_h)
            lat_t.append(lat_h)

    # to array
    move = np.array(move)
    vmean = np.array(vmean)
    vu = np.array(vu)
    vy = np.array(vy)
    p0 = np.array(p0)
    vmax = np.array(vmax)
    lon_t = np.array(lon_t)
    lat_t = np.array(lat_t)

    # longitude sign convention --> (0º,360º)
    lon_t[lon_t<0]= lon_t[lon_t<0] + 360

    # select interpolation data within the target domain area
    loc = []
    for i, (lo,la) in enumerate(zip(lon_t, lat_t)):
        if (lo<=lon01) & (lo>=lon00) & (la<=lat01) & (la>=lat00):
            loc.append(i)

    # storm track (pd.DataFrame)
    st = pd.DataFrame(index=time_input[loc],
                      columns=['move','vf','vfx','vfy','pn','p0','lon','lat','vmax'])
    st['move'] = move[loc]
    st['vf'] = vmean[loc]
    st['vfx'] = vu[loc]
    st['vfy'] = vy[loc]
    st['pn'] = 1013
    st['p0'] = p0[loc]
    st['lon'] = lon_t[loc]
    st['lat'] = lat_t[loc]
    # vmax is calculated from Pmin-Vmax basin-fitting when the value is not given 
    if wind.any() != None:
        st['vmax'] = vmax[loc]
    else:
        st['vmax'] = p1 * np.power(p0[loc],3) + p2 * np.power(p0[loc],2) + p3 * np.power(p0[loc],1) + p4   # [kt]

    # add some metadata
    st.x0 = x0
    st.y0 = y0
    st.R = 4

    return st, time_input[loc]

def entrance_coords(delta, gamma, x0, y0, R, lon0, lon1, lat0, lat1):
    '''
    Calculates storm track first coordinates

    delta, gamma               - storm track parameters
    x0, y0                     - site coordinates (longitude, latitude)
    R                          - radius (º)
    lon0, lon1, lat0, lat1     - computational coordinates (outer grid)
    '''

    # enter point in the radius
    xc = x0 + R * np.sin(delta * np.pi/180)
    yc = y0 + R * np.cos(delta * np.pi/180)

    # calculate angles that determine the storm boundary entrance  [degrees]
    ang_1 = np.arctan((lon1-xc)/(lat1-yc)) *180/np.pi       # upper right corner
    ang_2 = np.arctan((lon1-xc)/(lat0-yc)) *180/np.pi +180  # lower right
    ang_3 = np.arctan((lon0-xc)/(lat0-yc)) *180/np.pi +180  # lower left
    ang_4 = np.arctan((lon0-xc)/(lat1-yc)) *180/np.pi +360  # upper left

    if (gamma > ang_1) & (gamma < ang_2):
        x1 = lon1
        d = (x1 - xc) / np.sin(gamma * np.pi/180)
        y1 = yc + d * np.cos(gamma * np.pi/180)

    elif (gamma > ang_2) & (gamma < ang_3):
        y1 = lat0
        d = (y1 - yc) / np.cos(gamma * np.pi/180)
        x1 = xc + d * np.sin(gamma * np.pi/180)

    elif (gamma > ang_3) & (gamma < ang_4):
        x1 = lon0
        d = (x1 - xc) / np.sin(gamma * np.pi/180)
        y1 = yc + d * np.cos(gamma * np.pi/180)

    elif (gamma > ang_4) | (gamma < ang_1):
        y1 = lat1
        d = (y1 - yc) / np.cos(gamma * np.pi/180)
        x1 = xc + d * np.sin(gamma * np.pi/180)

    return x1, y1

def track_site_parameters(step, pmin, vmean, delta, gamma,
                          x0, y0, lon0, lon1, lat0, lat1, R, date_ini):
    '''
    Calculates storm track variables from storm track parameters within the study area
    (uses great circle)

    step                       - computational time step (in minutes)
    pmin, vmean, delta, gamma  - storm track parameters   (NOTE: vmean in [kt])
    x0, y0                     - site coordinates (longitude, latitude)
    lon0, lon1, lat0, lat1     - enter point in computational grid
    R                          - radius (º)
    date_ini                   - initial date 'yyyy-mm-dd HH:SS'
    great_circle               - default option
    '''

    # cubic polynomial fitting curve for each IBTrACS basin
    # to calculate vmax from relation Pmin-Vmax of the specific basin
    p1, p2, p3, p4 = ibtrac_basin_fitting(x0, y0)

    # storm boundary entrance coordinates
    x1, y1 = entrance_coords(delta, gamma, x0, y0, R, lon0, lon1, lat0, lat1)

    # calculate lon,lat storm coordinates
    xt, yt = [x1], [y1]
    i = 1
    glon, glat, baz = shoot(x1, y1, gamma+180, vmean*1.872 * i*step/60)  # velocity in [km/h]
    if glon < 0: glon += 360

    while (glon < lon1) & (glon > lon0) & (glat < lat1) & (glat > lat0):
        xt.append(glon)
        yt.append(glat)
        i += 1
        glon, glat, baz = shoot(x1, y1, gamma+180, vmean*1.872 * i*step/60)  # velocity in [km/h]
        if glon < 0: glon += 360

    frec = len(xt)

    # time array for SWAN input
    time_input = pd.date_range(date_ini, periods=frec, freq='{0}min'.format(step))

    # storm track (pd.DataFrame)
    st = pd.DataFrame(index=time_input,
                      columns=['move','vf','vfx','vfy','pn','p0','lon','lat','vmax'])
    st['move'] = gamma
    st['vf'] = vmean   # [kt]    
    st['pn'] = 1013
    st['p0'] = pmin
    st['vfx'] = vmean * np.sin((gamma+180) * np.pi/180)   # [kt]
    st['vfy'] = vmean * np.cos((gamma+180) * np.pi/180)   # [kt]
    st['vmax'] = p1 * np.power(pmin,3) + p2 * np.power(pmin,2) + p3 * np.power(pmin,1) + p4   # [kt]
    st['lon'] = xt
    st['lat'] = yt

    # add some metadata
    st.x0 = x0
    st.y0 = y0
    st.R = R

    return st


# VORTEX LIBRARY
# TODO 

