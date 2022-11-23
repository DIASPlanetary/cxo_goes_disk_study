# Script plotting JRM09 model over tessellation plots

# importing packages
import pandas as pd
import numpy as np
from scipy.io import readsav
import matplotlib.pyplot as plt
from astropy.time import Time
import matplotlib as mpl
import matplotlib.cm as cm

# Reading in the config.ini file containing any hard wired inputs
import configparser

config = configparser.ConfigParser()
config.read('config.ini')

# reading in JRM09 internal magnetic field model
f = readsav('map_jrm09.sav') # reading in .sav file containing JRM09 iso-contours
b = np.array(f['b']) # magnetic field magnitude
lon = np.arange(0, 360, 1) # longitude grid
theta = np.arange(-90, 91, 1) # latitude grid

eq_rad = 71492. # equatorial radius of Jupiter
pol_rad = 66854. # polar radiues of Jupiter
e = (71492 - 66854)/71492. # flattening of Jupiter
theta_graphic = np.arctan(np.tan(theta/180. * np.pi)/((1 - e) ** 2.))/np.pi * 180.
b_graphic = np.zeros((181,360)) # co-ordinate transformation from planetocentric to planetographic co-ords
for ilongitude in range(360):
    b_graphic[:, ilongitude] = np.interp(theta, theta_graphic, b[:, ilongitude]) # calculating new b values in planetographic co-ords

lon_grid, theta_grid = np.meshgrid(lon, theta) # creating grids to plot contours

obsID = str(config['inputs']['obsID']) # observation ID
folder_path = str(config['inputs']['folder_path'])

# Reading in catalogue containing information about each observation - used to get details like start date, exposure time, etc. to make figure title.
chandra_props = pd.read_excel('catalogue_all_data.xlsx')
index = np.where(chandra_props['ObsID'] == int(obsID))[0][0]
date_start = chandra_props['Start Date'][index]
exp_time = chandra_props['Exposure Time (ks)'][index]
J_E_dist = chandra_props['J-E Distance (AU)'][index]
esj_angle = chandra_props['E-S-J Angle (deg)'][index]

# Reading in data within PI filter for all of Jupiter 
PI_file = pd.read_csv(folder_path +  f'/{obsID}/{obsID}_photonlist_PI_filter_Jup_full_10_250.txt')
lat = np.array(PI_file['lat (deg)'].tolist())
lon = np.array(PI_file['SIII_lon (deg)'].tolist())
time = PI_file['# t(s)'].tolist()

# Making data periodic
lon2 = lon + 360.0
lon_neg = lon  - 360.0

lat2 = -lat + 180.
lat3 = -lat - 180.

lat_test = np.append(lat2, (lat, lat3, lat2, lat, lat3, lat2, lat, lat3))
lon_test = np.append(lon_neg, (lon_neg, lon_neg, lon, lon, lon, lon2, lon2, lon2))

# Plotting point maps for each obsID
# plotting contours
fig = plt.figure(figsize=(8,5))
ax = fig.add_subplot(1, 1, 1)
contours = plt.contour(lon_grid, theta_grid, b_graphic, colors='yellow', levels=[2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24], zorder=50)
levels = contours.levels
clabel = plt.clabel(contours, inline=True, fontsize=8, fmt='%1.0f', manual=False)
major_xticks = np.arange(350, -1, step=-50)
minor_xticks = np.arange(360, -1, step=-10)

plt.title(f'Observation start date: {date_start} \nExp time: {np.round(exp_time, 2)} ks - E-S-J angle: {np.round(esj_angle, 2)} deg - J-E dist: {np.round(J_E_dist, 2)} AU')
plt.scatter(lon_test, lat_test, s=2, color='k', zorder=20) #, alpha=0.3)
ax.set_xticks(major_xticks)
ax.set_xticks(minor_xticks, minor=True)
ax.tick_params(which='both', direction='in', top=True, right=True)
plt.xlabel('SIII Longitude (deg)'); plt.ylabel('SIII Latitude (deg)')

# Adding Voronoi regions to plot
points = np.column_stack((lon_test, lat_test))
from scipy.spatial import Voronoi, voronoi_plot_2d
vor = Voronoi(points)
fig_new = voronoi_plot_2d(vor, ax, show_vertices=False, show_points=False, point_size=2, color='w')
plt.xlim(362, -2)
plt.ylim(-90, 90)

# Function to calculate area of polygons.
def PolyArea(x,y):
    return 0.5 * abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))

area_list = []
lon_list = []
lat_list = []
norm = mpl.colors.Normalize(vmin=0, vmax=300, clip=True) # setting limits on colourbar
mapper = cm.ScalarMappable(norm=norm, cmap=cm.Blues_r)

# fixing regions so only calculate area within bounds
for i, point in enumerate(vor.points):
    if (point[0] >= 0) and (point[0] < 360) and (point[1] > -90) and (point[1] < 90):
        region = vor.regions[vor.point_region[i]]
        x_poly = np.array([vor.vertices[i][0] for i in region])
        y_poly = np.array([vor.vertices[i][1] for i in region])

        polygon = [vor.vertices[i] for i in region]

        # Find area of regions
        area_poly = PolyArea(x_poly, y_poly)
        area_list.append(area_poly)
        lon_list.append(point[0])
        lat_list.append(point[1])
        plt.fill(*zip(*polygon), color=mapper.to_rgba(area_poly), alpha=1.0, zorder=0) # colour polygons by area
cbar = fig.colorbar(mapper) # making colourbar
cbar.set_label(r'area ($\mathrm{deg}^2$)', rotation=270, labelpad=30)

# adding latitude constraints of planetary disk region
ax.hlines(-55, -2, 362, colors='r', zorder=60)
ax.hlines(45, -2, 362, colors='r', zorder=60)

# printing out results
print('# Polygons: ' + str(len(area_list)))
print('Max polygon area: ' + str(round(max(area_list), 2)) + ' deg\u00b2')
plt.tight_layout()
plt.savefig(str(folder_path) + f'/{obsID}/{obsID}_VRT_w_b_field.png', dpi=500)
plt.show()
