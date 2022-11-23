# Script to make total Chandra HRC-I Jupiter X-ray heat map. Results are displayed in 5 deg x 5 deg (lon, lat) bins. 

#relevant packages 
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

# Reading in excel catalogue file with all relevant info
catalogue = pd.read_excel('catalogue_all_data.xlsx')

# creating lon, lat grid to bin photons into for each observation
lon_arr = np.arange(361, step=5)
lat_arr = np.arange(start=-90, stop=91, step=5)
result = np.zeros((len(lon_arr) - 1, len(lat_arr) -1))

obsIDs = np.array(catalogue['ObsID'])
for obsID in obsIDs:
    # Reading in all Jupiter data with PI filter 
    raw_file  = pd.read_csv(f'{obsID}/{obsID}_photonlist_PI_filter_Jup_full_10_250.txt')
    lat = raw_file['lat (deg)'].tolist()
    lon = raw_file['SIII_lon (deg)'].tolist()
    result += np.histogram2d(lon, lat, bins=(lon_arr, lat_arr))[0] # adding photons

# making plot
fig2 = plt.figure(figsize=(9,4))
ax2 = fig2.add_subplot(111)

major_xticks = np.arange(350, -1, step=-50)
minor_xticks = np.arange(360, -1, step=-10)
major_yticks = np.arange(-80, 81, step=20)
minor_yticks = np.arange(-90, 91, step=10)

ax2.set_xticks(major_xticks)
ax2.set_xticks(minor_xticks, minor=True)
ax2.set_yticks(major_yticks)
ax2.set_yticks(minor_yticks, minor=True)

X, Y = np.meshgrid(lon_arr, lat_arr)
plt.pcolormesh(X, Y, result.T, cmap=plt.cm.magma, vmax=40) # plotting lon, lat grid with colourbar to express counts
cbar = plt.colorbar()
cbar.set_label(r'Number of photons / $5^{\circ}$ S3 lon x $5^{\circ}$ lat', rotation=270, labelpad=30)
plt.xlim(360, 0)
ax2.hlines(-55, 0, 360, colors='w')
ax2.hlines(45, 0, 360, colors='w')
plt.title('Total Chandra HRC-I Jupiter X-ray map')
plt.xlabel('SIII Longitude (degrees)'); plt.ylabel('SIII Latitude (degrees)')
plt.tight_layout()
plt.savefig('Jup_Xray_map_saturated_10_250.png', dpi=500)
plt.show()
