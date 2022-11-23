# Script that has as its output a script comparing the net jovian disk light curve with GOES solar X-ray flux over the same time interval, after accounting for light travel time. 

# Relevant packages
import astropy.units as u
from astropy.time import Time                   #convert between different time coordinates
from astropy.time import TimeDelta              #add/subtract time intervals 
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from math import ceil

# Reading in the config.ini file containing any hard wired inputs
import configparser
config = configparser.ConfigParser()
config.read('config.ini')

obsID = str(config['inputs']['obsID'])

# Reading in excel catalogue file with all relevant info
catalogue = pd.read_excel('catalogue_all_data.xlsx')
index = np.where(catalogue['ObsID'] == int(obsID))[0][0]

disk_area = catalogue['Disk Area (px^2)'][index]
bg_area = catalogue['Background Area (px^2)'][index]

J_E_dist = catalogue['J-E Distance (AU)'][index]
ang_diam = catalogue['Ang Diam (arcsec)'][index]
esj_ang = catalogue['E-S-J Angle (deg)'][index]

start_date = catalogue['Start Date'][index]
tstart_evt = catalogue['Tstart'][index]
tstop_evt = catalogue['Tstop'][index]

exp_time = catalogue['Exposure Time (ks)'][index]

# Panel plot with Jup disk lc in first panel and goes lc in second panel with times overlapping
fig = plt.figure(figsize=(8,5))

# Reading in background file with PI filter
bg_file = pd.read_csv(f'{obsID}/{obsID}_bg_times.txt')
bg_t = bg_file['Time'].tolist()
t_bg_min = [(bg_t[i] - tstart_evt)/60 for i in range(len(bg_t))]

# Reading in jovian disk data PI filter 
PI_file = pd.read_csv(f'{obsID}/{obsID}_photonlist_PI_10_250.txt')
time = PI_file['# t(s)'].tolist()
time_min = [(time[i] - tstart_evt)/60 for i in range(len(time))]
tstop_min = (tstop_evt - tstart_evt)/60

bin_size = 5 # select 5-min bins
bin_arr = np.arange(ceil(tstop_min + 5), step=bin_size)

# generate histogram for disk photons
counts, bins = np.histogram(time_min, bins=bin_arr)
norm_cts = counts/bin_size # normalise by dividing by bin size
counts_per_area = norm_cts/disk_area # account for area of disk

# and now the same for the background
cts_bg, bins = np.histogram(t_bg_min, bins=bin_arr)
norm_cts_bg = cts_bg/bin_size
bg_cts_per_area = norm_cts_bg/bg_area # account for area of background region

# calculating net count rate
net_cr = counts_per_area - bg_cts_per_area

centre = (bins[:-1] + bins[1:])/2 # define bin centre to plot histogram

# making plot
ax = fig.add_subplot(211)
loc, label =  plt.xticks()
ax.plot(centre, net_cr, 'g', zorder=10, label='Net count rate')

major_xticks_time = np.arange(0, ceil(tstop_min + 6), step=100)
minor_xticks_time = np.arange(0, ceil(tstop_min + 6), step=10)

ax.set_xticks(major_xticks_time)
ax.set_xticks(minor_xticks_time, minor=True)
ax.tick_params(which='both', direction='in', top=True, right=True, labelbottom=False)
ax.set_xticks([])
plt.ylabel(r'$\mathrm{counts\;min^{-1}\;px^{-2}}$')
plt.xlim(-5, ceil(tstop_min) + 8)
ax.text(0.02, 0.95, '(a)', transform=ax.transAxes, fontsize=14, va='top')
plt.title(f'Observation start date: {start_date} \nExp time: {np.round(exp_time, 2)} ks - E-S-J angle: {np.round(esj_ang, 2)} deg - J-E dist: {np.round(J_E_dist, 2)} AU', pad=20, fontsize=12)
plt.ylim(-2.6e-5, 5.3e-5)
ax.legend(loc='upper right', fontsize=9, handlelength=1.0)

# goes data that aligns with cxo data
ax2 = fig.add_subplot(212)

# Reading in GOES data
GOES_file = pd.read_csv(f'{obsID}/{obsID}_GOES_XRS_bhardwaj.txt')
t_goes = np.array(GOES_file['# time (min)'].tolist())
flux_goes = np.array(GOES_file['flux (W/m^2)'].tolist())

# plotting GOES data
ax2.plot(t_goes, flux_goes, 'k', label = r'1.0 - 8.0 $\mathrm{\AA}$')
ax2.set_xticks(major_xticks_time)
ax2.set_xticks(minor_xticks_time, minor=True)

goes_D = {'X':1e-4, # sets up flare classes for plotting
          'M':1e-5,
          'C':1e-6,
          'B':1e-7,
          'A':1e-8,}

ax2.set_yscale('log')
ax2.set_ylim(1e-8, 1e-4)
plt.ylabel(r'GOES X-ray Flux ($\mathrm{W m^{-2}}$)')
ax2.tick_params(which='both', direction='in', top=True, right=True)
plt.xlabel('Time (min)')
ax2a = ax2.twinx()
ax2a.set_yscale('log')
ax2a.set_ylim(1e-8, 1e-4)

# Indicate on plot the different flare classifications
for ii in list(goes_D.values()):
    ax2a.axhline(ii, linestyle = '--', alpha=0.1)
ax2a.set_yticks(list(goes_D.values()))
ax2a.set_yticklabels(goes_D.keys())
ax2a.tick_params(axis='both', which='both')

plt.xlim(-5, ceil(tstop_min) + 8)
plt.subplots_adjust(hspace=0.05)
ax2.text(0.02, 0.95, '(b)', transform=ax2.transAxes, fontsize=14, va='top')
ax2.legend(loc='upper right', fontsize=9, handlelength=1.0)
plt.savefig(f'{obsID}/{obsID}_cxo_goes_comp.png', dpi=500, bbox_inches='tight')
plt.show()
