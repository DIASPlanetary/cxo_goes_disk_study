# Script to find GOES solar X-ray flux data for the same time interval as the CXO observation, after accounting for light travel time.

# Authors: D. M. Weigt, Seán McEntee, Brad Snios

# Relevant packages
import pandas as pd
import numpy as np
import sunpy
import sunpy.timeseries as ts
from sunpy.net import Fido
from sunpy.net import attrs as a

import astropy.units as u
from astropy.time import Time                   #convert between different time coordinates
from astropy.time import TimeDelta              #add/subtract time intervals 

# Reading in the config.ini file containing any hard wired inputs
import configparser
config = configparser.ConfigParser()
config.read('config.ini')

obsID = str(config['inputs']['obsID'])

# Reading in excel catalogue file with all relevant info
catalogue = pd.read_excel('catalogue_all_data.xlsx')
index = np.where(catalogue['ObsID'] == int(obsID))[0][0]

start_date = catalogue['Start Date'][index]
end_date = catalogue['End Date'][index]
exp_evt = catalogue['Exposure Time (ks)'][index]
tstart_evt = catalogue['Tstart'][index]
tstop_evt = catalogue['Tstop'][index]

# Brad's Horizons code to extract the ephemeris file
from astroquery.jplhorizons import Horizons     #automatically download ephemeris 

# The start and end times are taken from the horizons file.
tstart_eph=Time(tstart_evt, format='cxcsec')
tstop_eph=Time(tstop_evt, format='cxcsec')
dt = TimeDelta(0.125, format='jd')
# Below sets the parameters of what observer the ephemeris file is generated form. For example, '500' = centre of the Earth, '500@-151' = CXO
obj_jup_cxo = Horizons(id=599,location='500@-151',epochs={'start':tstart_eph.iso, 'stop':tstop_eph.iso, 'step':'1m'}, id_type='majorbody')

eph_jup_cxo = obj_jup_cxo.ephemerides()

obj_sun_jup = Horizons(id=599,location='500@10',epochs={'start':tstart_eph.iso, 'stop':tstop_eph.iso, 'step':'1m'}, id_type='majorbody')

eph_sun_jup= obj_sun_jup.ephemerides()

"""Setting up event times"""
jup_cxo_lt = np.mean(eph_jup_cxo['lighttime'])
sun_jup_lt = np.mean(eph_sun_jup['lighttime'])
sun_earth_lt = np.mean(eph_sun_jup['earth_lighttime'])
# earth_jup_lt = np.mean(eph_earth_jup['lighttime'])
lt = (TimeDelta((sun_jup_lt + jup_cxo_lt - sun_earth_lt) * u.min)).datetime

cxo_tstart = Time(start_date, format='iso') - lt # taking away lt from Earth to Jupiter and then to CXO
cxo_exp = TimeDelta(exp_evt * u.s).datetime # replace <exp> with exposure time of event in seconds
cxo_tend = cxo_tstart + cxo_exp

# defining start and end times 10 minutes either side of the CXO observation after accounting for lt 
tstart = (cxo_tstart - 600 * u.s).iso
tend = (cxo_tend + 600 * u.s).iso

if (int(obsID) == 1862) or (int(obsID) == 2519):
    sat_num = 10 # changes based on date of observation
else:
    sat_num = 15
result = Fido.search(a.Time(tstart, tend), a.Instrument("XRS"))
result_goes_cxo = Fido.search(a.Time(tstart, tend), a.Instrument("XRS"), a.goes.SatelliteNumber(sat_num))
#Note: GOES satellite number will change depending on year! 'result' will show you what satellites are available
goes_cwd = str(config['inputs']['goes_cwd']) # will need to change this in config.ini file
files_cxo = Fido.fetch(result_goes_cxo, path=goes_cwd) # fetches the GOES data via sunpy from NOAA dtatabase
goes_cxo = ts.TimeSeries(files_cxo, concatenate=True) # creates time series of selected GOES data

# Converting goes lc to something similar to Jupiter light curve
goes_time_arr = Time(goes_cxo.index).cxcsec
goes_lc_cxo = goes_time_arr[np.where((goes_time_arr > cxo_tstart.cxcsec) & (goes_time_arr < cxo_tend.cxcsec))[0]]
goes_lc_cxo_flux = np.array(goes_cxo.quantity('xrsb')[np.where((goes_time_arr > cxo_tstart.cxcsec) & (goes_time_arr < cxo_tend.cxcsec))[0]])
goes_lc_cxo = [(x - cxo_tstart.cxcsec)/60 for x in goes_lc_cxo]

# Writing goes data to text file so we don't have to deal with sunpy each time.
np.savetxt(f'{obsID}/{obsID}_GOES_XRS_bhardwaj.txt', np.c_[goes_lc_cxo, goes_lc_cxo_flux], delimiter=',', header='time (min),flux (W/m^2)')

