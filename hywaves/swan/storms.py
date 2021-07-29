#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os.path as op
import numpy as np
import pandas as pd
import xarray as xr
from datetime import timedelta
from numpy import polyfit

from .geo import shoot, gc_distance


###############################################################################
# IBTrACS database is used to extract historical storm track data. However, 
# specialized centers (RSMCs) employ different procedures to compute the maximum
# sustained winds (Vmax). Therefore a conversion factor is applied to obtain the 
# estimated 1-min average maximum sustained winds (input variable for "make_vortex").
# Harper et al. (2010) recommends 0.93 to convert from 1-min to 10-min (at sea conditions)
# instead of the traditional 0.88 factor.
###############################################################################

# dictionary of conversion factors (1-min to X-min)
d_vns_windfac = {
        'USA':      1,        # 1-min avg, (RSMC Miami, RSMC Honolulu, NHC, JTWC, CPHC)
        'TOKYO':    0.93,     # 10-min avg, (RSMC Tokyo (JMA))
        'CMA':      0.99,     # 2-min avg, (Chinese Meteorological Administration Shanghai Typhoon Institute)
        'HKO':      0.93,     # 10-min avg, (Hong Kong Observatory (operated by JMA))
        'NEWDELHI': 0.99,     # 3-min avg, (RSMC New Delhi)
        'REUNION':  0.93,     # 10-min avg, (RSMC La Reunion (MeteoFrance))
        'BOM':      0.93,     # 10-min avg, (Australian Bureau of Meteorology)
        'WELLINGTON': 0.93,   # 10-min avg, (New Zealand MetService)
        'NADI':     0.93,     # 10-min avg, (RSMC Nadi, Fiji Meteorological Service)
        'WMO':      1,        # mixed, it is assumed 1-min!!!
        }

# dictionary of centers names for calling variables
d_vns_idcenter = {
        'WMO':          'wmo',
        'USA':          'usa',
        'TOKYO':        'tokyo',
        'CMA':          'cma',
        'HKO':          'hko',
        'NEWDELHI':     'newdelhi',
        'REUNION':      'reunion',
        'BOM':          'bom',
        'WELLINGTON':   'wellington',
        'NADI':         'nadi',
        }

# dictionary of centers names for calling variables
d_vns_basinscenter = {
        'WMO':          ['NA','SA','EP','WP','SP','NI','SI'],
        'USA':          ['NA','SA','EP','WP','SP','NI','SI'],
        'TOKYO':        ['WP'],
        'CMA':          ['WP'],
        'HKO':          ['WP'],
        'NEWDELHI':     ['NI'],
        'REUNION':      ['SI'],
        'BOM':          ['SP','SI'],
        'WELLINGTON':   ['SP'],
        'NADI':         ['SP'],
        }

def get_dictionary_center(d_vns_idcenter, center='WMO'):
    'Returns the dictionary of a center (IBTrACS)'
    
    source = d_vns_idcenter[center]
    
    # lon, lat
    if source=='wmo':   dict_lon, dict_lat = 'lon','lat'
    else:               dict_lon, dict_lat = source+'_lon',source+'_lat'
    
    # radii RMW is provided by centers: USA, REUNION, BOM
    if source=='usa' or source=='reunion' or source=='bom':
        dict_radii = source+'_rmw'
    else: dict_radii = None
        
    # variables dictionary
    d_vns = {
        'source':       center,
        'time':         'time',
        'longitude':    dict_lon,
        'latitude':     dict_lat,
        'pressure':     source+'_pres',
        'maxwinds':     source+'_wind',
        'rmw':          dict_radii,
    }
    
    return d_vns


# generate pmin-vmax polynomial fit coefficients

def Extract_basin_storms(xds, id_basin):
    '''
    Selects storms originated in a given a basin
    '''
    
    # select genesis basin (np.bytes_)
    origin_basin = xds.sel(date_time=xds.date_time[0]).basin.values
    origin_basin = np.asarray([c.decode('UTF-8') for c in origin_basin])
    
    # select storms within given basin
    storm_pos = np.where(origin_basin == id_basin)[0]
    
    # extract storms for a given genesis basin
    xds_basin = xds.sel(storm=storm_pos)
    
    return xds_basin

#d_vns_nadi = {
#    'source':    'NADI',
#    'basins':    ['SP'],
#    'longitude': 'nadi_lon',
#    'latitude':  'nadi_lat',
#    'time':      'time',
#    'pressure':  'nadi_pres',
#    'wind':      'nadi_wind',
#    'rmw':       None,
#}
def ibtracs_fit_pmin_vmax(d_vns_idcenter):
    '''
    Generate third order polynomial fit coefficients for statistical relationship 
    between minimum central pressure and maximum wind speed (at 1-min avg).
    
    Maximum wind speed are converted to 1-min average for each RSMC center.
    '''
    
    p_ibtracs = op.join(op.dirname(op.realpath(__file__)), 'resources', 'IBTrACS.nc')
    xds_ibtracs = xr.open_dataset(p_ibtracs)
    
    # all basins, centers
    basin_all = ['NA','SA','WP','EP','SP','NI','SI']
    center_all = ['USA','TOKYO','CMA','HKO','NEWDELHI','REUNION','BOM','WELLINGTON','NADI','WMO']
    
    # set storms id as coordinate
    xds_ibtracs['stormid'] = (('storm'), xds_ibtracs.storm.values)
    xds_ibtracs.set_coords('stormid')
    
    # polynomial order
    N = 3   
    coef_fit = np.nan * np.zeros((len(center_all), len(basin_all), N+1))
#    pres_data = np.nan * np.zeros((len(center_all), len(basin_all), 50000))
#    wind_data = np.nan * np.zeros((len(center_all), len(basin_all), 50000))

    # loop for all centers
    for ic,center in enumerate(center_all):

        d_vns_center = get_dictionary_center(d_vns_idcenter, center=center)
        dict_pres = d_vns_center['pressure']    # mbar
        dict_wind = d_vns_center['maxwinds']    # kt
        basin_ids = d_vns_basinscenter[center]  # basins

        # loop for all basins
        for basin in basin_ids:
            
            ibasin = np.where(basin == np.array(basin_all))[0][0]
            
            # select storms at basin X
            xds_basin = Extract_basin_storms(xds_ibtracs, basin)
    
            # extract basin data: pressure, wind, landbasin
            Pbasin = xds_basin[dict_pres].values.reshape(xds_basin.storm.size * xds_basin.date_time.size)
            Wbasin = xds_basin[dict_wind].values.reshape(xds_basin.storm.size * xds_basin.date_time.size)
            Wbasin /= d_vns_windfac[center]     # winds are converted to 1-min
            landbasin = xds_basin['dist2land'].values.reshape(xds_basin.storm.size * xds_basin.date_time.size)
    
            PWbasin_s = np.column_stack((Pbasin, Wbasin, landbasin))
    
            # index for removing nans (including landmask)
            ix_nonan = ~np.isnan(PWbasin_s).any(axis=1)
            PWbasin_s = PWbasin_s[ix_nonan]
    
            # Fitting Polynomial Regression to the dataset 
            X, y = PWbasin_s[:,0], PWbasin_s[:,1]
            u = polyfit(X, y, deg=N)
    
            # store
            coef_fit[ic,ibasin,:] = u
#            pres_data[ic,ib,:PWbasin_s.shape[0]] = PWbasin_s[:,0].T
#            wind_data[ic,ib,:PWbasin_s.shape[0]] = PWbasin_s[:,1].T

    # store
    xds = xr.Dataset(
        {
            'coef': (('center','basin','polynomial'), coef_fit),
#            'pres': (('center','basin','data'), pres_data),
#            'wind': (('center','basin','data'), wind_data),
        },
    {
        'center': np.asarray(center_all),
        'basin': np.asarray(basin_all),
    })
    # save
    p_coef = op.join(op.dirname(op.realpath(__file__)), 'resources', 'ibtracs_coef_pmin_vmax.nc')
    xds.to_netcdf(p_coef)

#    xds.coef.attrs = 'Pressure (mbar), Wind speed (kt, 1-min avg)'
#    xds.pres.attrs['name'] = 'Pressure'
#    xds.pres.attrs['units'] = 'mbar'
#    xds.wind.attrs['name'] = 'MaxWinds'
#    xds.wind.attrs['units'] = 'kt (1-min)'
    
    return xds

###############################################################################



# STORM TRACK LIBRARY

def get_category(ycpres):
    'Defines storm category according to minimum pressure center'

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

def historic_track_preprocessing(xds, center='USA'): #d_vns):
    '''
    Historic track is preprocessed, by removing NaN data, apply longitude
    convention [0º-360º], change time format and define storm category

    xds:   historic track dataset (storm dimension)
    d_vns: dictionary for longitude, latitude, time, pressure and wind varnames
    
    returns variables:  
           time, lat, lon, pressure, wind, timestep, category, mean translational speed
    '''
    
    ###########################################################################
    # dictionary for IBTrACS center
    d_vns = get_dictionary_center(d_vns_idcenter, center=center)
    ###########################################################################

    # get names of vars
    nm_lon = d_vns['longitude']
    nm_lat = d_vns['latitude']
    nm_prs = d_vns['pressure']
    nm_tim = d_vns['time']
    nm_win = d_vns['maxwinds']
    ###########################################################################
    nm_rmw = d_vns['rmw']
    ###########################################################################

    # get var time
    ytime = xds[nm_tim].values         # datetime64

    # remove NaTs (time)
    ycpres = xds[nm_prs].values[~np.isnat(ytime)]     # minimum pressure
    ylat_tc = xds[nm_lat].values[~np.isnat(ytime)]    # latitude
    ylon_tc = xds[nm_lon].values[~np.isnat(ytime)]    # longitude
    ywind = xds[nm_win].values[~np.isnat(ytime)]      # wind speed [kt]
    ###########################################################################
    if nm_rmw:  yradii = xds[nm_rmw].values[~np.isnat(ytime)]  # RMW [nmile]
    ###########################################################################
    ytime = ytime[~np.isnat(ytime)]

    # remove NaNs (pressure)
    ylat_tc = ylat_tc[~np.isnan(ycpres)]
    ylon_tc = ylon_tc[~np.isnan(ycpres)]
    ywind = ywind[~np.isnan(ycpres)]
    ###########################################################################
    if nm_rmw:  yradii = yradii[~np.isnat(ycpres)]
    if not nm_rmw:  yradii = None
    ###########################################################################
    ytime = ytime[~np.isnan(ycpres)]
    ycpres = ycpres[~np.isnan(ycpres)]
    
    ###########################################################################
    # convert (X)min to 1-min avg winds (depends on center)
    ywind = ywind / d_vns_windfac[center]
    ###########################################################################

    # longitude convention: [0º,360º]
    ylon_tc[ylon_tc<0] = ylon_tc[ylon_tc<0] + 360

    # round dates to hour
    round_to = 3600
    st_time = []
    for i in range(len(ytime)):
        dt = ytime[i].astype('datetime64[s]').tolist()
        seconds = (dt - dt.min).seconds
        rounding = (seconds+round_to/2) // round_to * round_to
        out = dt + timedelta(0, rounding-seconds, -dt.microsecond)
        st_time.append(out)
    st_time = np.asarray(st_time)

    # storm coordinates timestep [hours]
    ts = (st_time[1:] - st_time[:-1])
    ts = [ts[i].total_seconds() / 3600 for i in range(ts.size)]

    # storm coordinates category
    categ = get_category(ycpres)
    
    # calculate Vmean
    RE = 6378.135   # earth radius [km]
    vmean = []
    for i in range(0, len(st_time)-1):

        # consecutive storm coordinates
        lon1, lon2 = ylon_tc[i], ylon_tc[i+1]
        lat1, lat2 = ylat_tc[i], ylat_tc[i+1]

        # translation speed 
        arcl_h, gamma_h = gc_distance(lat2, lon2, lat1, lon1)
        r = arcl_h * np.pi / 180.0 * RE     # distance between consecutive coordinates (km)
        vmean.append(r / ts[i] / 1.852)     # translation speed (km/h to kt) 
        
    # mean value
    vmean = np.mean(vmean)  # [kt]
    
    return st_time, ylat_tc, ylon_tc, ycpres, ywind, ts, categ, vmean, yradii

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
                                 wind=None, great_circle=False, fit=False, 
                                 interpolation=True, mode='first'):
    '''
    Calculates storm track variables from storm track parameters and interpolates
    track points in between historical data (for "dt_comp" time step)

    st_time                    - storm dates
    ylon_tc, ylat_tc           - storm coordinates (longitude, latitude)
    ycpres, ywind              - storm minimum pressure, maximum winds
    x0, y0                     - target coordinates (longitude, latitude)
    lat0, lon0, lat1, lon1     - bound limits for numerical domain (longitude, latitude)
    ts                         - storm coordinates time step [hours]
    dt_comp                    - computation time step [minutes]
    
    wind                       - historical ywind can be assigned or 'None' for vmax empirical estimate
    great_circle               - True for distances over great circle 
    fit                        - True for fitting empirical vmax when ywind=0 (start of storm)
    imterpolation              - True for storm variables interpolation (historic storms)
                                 False for mean storm variables (storm segments constant)
    mode ('first','mean')      - when interpolation is activated, chooses value for constant segments
    
    returns:  dataframe with interpolated storm coordinate variables 
    '''

    RE = 6378.135   # earth radius [km]
    
    # check file of IBTrACS fitting coefficients (Pmin-Vmax)
    p_coef = op.join(op.dirname(op.realpath(__file__)), 'resources', 'ibtracs_coef_pmin_vmax.nc')
    if not p_coef:  xds_coef = ibtracs_fit_pmin_vmax(d_vns_idcenter)
    else:           xds_coef = xr.open_dataset(p_coef)

    # select polynomial fitting coefficients (center, basin)
    p1, p2, p3, p4 = ibtrac_basin_fitting(x0, y0)

    # storm variables
    time_storm = list(st_time)  # datetime format
    pmin = list(ycpres)
    lat = list(ylat_tc)
    lon = list(ylon_tc)
    if wind.any() != None:  # wind variable, no data is filled with IBTrACS fitting coefficients
        mwind = wind
        if fit:
            wind_fitting = p1 * np.power(pmin,3) + p2 * np.power(pmin,2) + p3 * np.power(pmin,1) + p4
            pos = np.where(mwind==0)
            mwind[pos] = wind_fitting[pos]

    # number of time steps between consecutive interpolated storm coordinates in order
    # to  match SWAN computational time step
    ts_h = ts                               # hours
    ts = np.asarray(ts) * 60 / dt_comp      # number of intervals

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

        # consecutive storm coordinates
        lon1, lon2 = lon[i], lon[i+1]
        lat1, lat2 = lat[i], lat[i+1]

        # translation speed 
        arcl_h, gamma_h = gc_distance(lat2, lon2, lat1, lon1)
        r = arcl_h * np.pi / 180.0 * RE     # distance between consecutive storm coordinates [km]
        dx = r / ts[i]                      # distance during time step
        tx = ts_h[i] / ts[i]                # time period during time step
        vx = float(dx) / tx / 3.6           # translation speed [km to m/s]
        vx = vx /0.52                       # translation speed [m/s to kt]

        for j in range(int(ts[i])):
            # append storm track parameters
            move.append(gamma_h)
            vmean.append(vx)
            vu.append(vx * np.sin((gamma_h+180)*np.pi/180))
            vy.append(vx * np.cos((gamma_h+180)*np.pi/180))
            pn.append(1013)
            
            # append pmin, wind with/without interpolation along the storm track
            if interpolation:       p0.append(pmin[i] + j* (pmin[i+1]-pmin[i])/ts[i])
            if not interpolation:   
                if mode=='mean':    p0.append(np.mean((pmin[i], pmin[i+1])))
                elif mode=='first': p0.append(pmin[i])
            
            if wind.any() != None:
                if interpolation:   vmax.append(mwind[i] + j* (mwind[i+1]-mwind[i])/ts[i])    #[kt]
                if not interpolation:   
                    if mode=='mean': vmax.append(np.mean((mwind[i], mwind[i+1])))   
                    if mode=='first':vmax.append(mwind[i])     #[kt]

            # calculate timestep lon, lat
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

    # longitude convention [0º,360º]
    lon_t[lon_t<0]= lon_t[lon_t<0] + 360

    # select interpolation coordinates within the target domain area
    loc = []
    for i, (lo,la) in enumerate(zip(lon_t, lat_t)):
        if (lo<=lon01) & (lo>=lon00) & (la<=lat01) & (la>=lat00):
            loc.append(i)

    # storm track (pd.DataFrame)
    st = pd.DataFrame(index=time_input[loc],
                      columns=['move','vf','vfx','vfy','pn','p0','lon','lat','vmax'])
    
    st['move'] = move[loc]      # gamma, forward direction
    st['vf'] = vmean[loc]       # translational speed [kt]
    st['vfx'] = vu[loc]         # x-component
    st['vfy'] = vy[loc]         # y-component
    st['pn'] = 1013             # average pressure at the surface [mbar]
    st['p0'] = p0[loc]          # minimum central pressure [mbar]
    st['lon'] = lon_t[loc]      # longitude coordinate
    st['lat'] = lat_t[loc]      # latitude coordinate
    # maximum wind speed (if no value is given it is assigned the empirical Pmin-Vmax basin-fitting)
    if wind.any() != None:  st['vmax'] = vmax[loc]  # [kt]
    else:                   st['vmax'] = p1 * np.power(p0[loc],3) + p2 * np.power(p0[loc],2) + \
                                        p3 * np.power(p0[loc],1) + p4   # [kt]

    # add some metadata
    # TODO: move to st.attrs (this metada gets lost with any operation with st)
    st.x0 = x0
    st.y0 = y0
    st.R = 4

    return st, time_input[loc]

def entrance_coords(delta, gamma, x0, y0, R, lon0, lon1, lat0, lat1):
    '''
    Calculates storm track initial coordinates

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

    step                       - computational time step (minutes)
    pmin, vmean, delta, gamma  - storm track parameters   (NOTE: vmean in [kt])
    x0, y0                     - site coordinates (longitude, latitude)
    lon0, lon1, lat0, lat1     - enter point in computational grid
    R                          - radius (º)
    date_ini                   - initial date 'yyyy-mm-dd HH:SS'
    great_circle               - default option
    '''

    # cubic polynomial fitting coefficients for IBTrACS basins Pmin-Vmax relationship
    p1, p2, p3, p4 = ibtrac_basin_fitting(x0, y0)

    # storm entrance coordinates at the domain boundary
    x1, y1 = entrance_coords(delta, gamma, x0, y0, R, lon0, lon1, lat0, lat1)

    # calculate computation storm coordinates
    xt, yt = [x1], [y1]
    i = 1
    glon, glat, baz = shoot(x1, y1, gamma+180, vmean*1.852 * i*step/60)  # velocity in [km/h]
    if glon < 0: glon += 360
    while (glon < lon1) & (glon > lon0) & (glat < lat1) & (glat > lat0):
        xt.append(glon)
        yt.append(glat)
        i += 1
        glon, glat, baz = shoot(x1, y1, gamma+180, vmean*1.852 * i*step/60)  # velocity in [km/h]
        if glon < 0: glon += 360
    frec = len(xt)

    # time array for SWAN input
    time_input = pd.date_range(date_ini, periods=frec, freq='{0}min'.format(step))

    # storm track (pd.DataFrame)
    st = pd.DataFrame(index=time_input,
                      columns=['move','vf','vfx','vfy','pn','p0','lon','lat','vmax'])
    
    st['move'] = gamma      # gamma, forward direction
    st['pn'] = 1013         # average pressure at the surface [mbar]
    st['p0'] = pmin         # minimum central pressure [mbar]
    st['vf'] = vmean        # translation speed [kt]    
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
