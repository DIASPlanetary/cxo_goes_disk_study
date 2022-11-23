# Script plotting the output of the go_chandra algorithm, showing the photons that are selected as jovian photons, and also the photons included in the background region.

#relevant packages 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.patches import Rectangle

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
exp_time = catalogue['Exposure Time (ks)'][index]

# JPL Horizons data
ang_diam = catalogue['Ang Diam (arcsec)'][index]
tilt_ang_deg = catalogue['North Pole Tilt Angle (deg)'][index]
tilt_ang_rad = np.deg2rad(tilt_ang_deg) # converting north pole tilt angle from degrees to radians

# Reading in ellipse data
ellipse_data = pd.read_csv(f'{obsID}/{obsID}_photonlist_full_obs_ellipse.txt')
x_ell = np.array(ellipse_data.iloc[:,1])
y_ell = np.array(ellipse_data.iloc[:,2])

# Reading in photons in full image
all_photons = pd.read_csv(f'{obsID}/{obsID}_all_photons.txt', header=None, delimiter=' ')
bg_reg = all_photons[(all_photons[0] > -100) & (all_photons[0] < 100) & (all_photons[1] > -100) & (all_photons[1] < 100)]
x = np.array(bg_reg.iloc[:,0])
y = np.array(bg_reg.iloc[:,1])

ellip = 0.3543163789650 # ellipticity of Jupiter
R_eq = (ang_diam/2.)/np.cos(tilt_ang_rad) # Equatorial radius of Jupiter
R_p = R_eq * np.sqrt(1 - ellip**2) # Polar radius of Jupiter

# Defining ellipse 1.5 times the size of Jupiter - photons outside of this region were included in the background
R_eq_bg = 1.5 * (ang_diam/2.)/np.cos(tilt_ang_rad)
R_p_bg = R_eq_bg * np.sqrt(1 - ellip**2)

ell = Ellipse((0,0), R_eq * 2., R_p * 2., tilt_ang_deg, facecolor='None', edgecolor='r', lw=1, zorder=50) # ellipse outline
ell_bg = Ellipse((0,0), R_eq_bg * 2., R_p_bg * 2., tilt_ang_deg, facecolor='None', edgecolor='k', lw=1, zorder=50) # ellipse outline

box = Rectangle((-100, -100), 200, 200, facecolor='None', edgecolor='k', lw=1, zorder=50) # boundary of background region

# get position of photons that lie in background region
x_bg = []
y_bg = []
for i in range(len(x)):
    if (x[i] * np.cos(tilt_ang_rad) + y[i] * np.sin(tilt_ang_rad)) ** 2./(R_eq_bg ** 2.) + (x[i] * np.sin(tilt_ang_rad) - y[i] * np.cos(tilt_ang_rad)) ** 2./(R_p_bg ** 2.) > 1.0:
        x_bg.append(x[i])
        y_bg.append(y[i])

# Schematic diagram showing location Jup and bg regions
fig, ax = plt.subplots(figsize=(7,7))
ax.add_artist(ell)
ax.add_artist(ell_bg)
ax.add_artist(box)
ax.plot(x, y, 'w.', markersize=1, zorder=0) # excluding photons that do not lie in either the Jupiter region or the bg region 
ax.plot(x_bg, y_bg, 'k.', markersize=1, zorder=2,)
ax.plot(x_ell, y_ell, 'r.', markersize=1, zorder=1)
ax.plot([-105, -105], [-105, 105], 'k', label = 'Background')
ax.plot([-105, -105], [-105, 105], 'r', label = 'Jupiter')
ax.set_xlim(-100, 100)
ax.set_ylim(-100, 100)
ax.set_xlabel('x (arcsec)'); ax.set_ylabel('y (arcsec)')
ax.set_title(f'ObsID: {obsID} \n{start_date} - {end_date} ({np.round(exp_time, 2)} ks)')
plt.legend(loc='upper right', facecolor='w', framealpha=0.8, markerscale=5.)
plt.tight_layout()
plt.savefig(f'{obsID}/{obsID}_bg_plus_ell.png', dpi=500)
plt.show()
