# Script that creates a 3 panel plot comparing sunspot number with jovian disk count rate (calculated from Chandra High Resolution Camera (HRC) observations of Jupiter).

#relevant packages 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.time import Time

# Reading in the config.ini file containing any hard wired inputs
import configparser
config = configparser.ConfigParser()
config.read('config.ini')

sun_data_path = str(config['inputs']['sun_data_path']) # make sure this entry in config file matches your description

# Reading in sunspot data
sun_data = pd.read_csv(sun_data_path, header=None, sep=";")
time_ssp = np.array(sun_data.iloc[:, 2])
ssp = np.array(sun_data.iloc[:, 3])


# Reading in excel catalogue file with all relevant info
catalogue = pd.read_excel('catalogue_all_data.xlsx')

start_date = np.array(catalogue['Start Date'])
dec_year = [Time(x).decimalyear for x in start_date]
exp_time = np.array(catalogue['Exposure Time (ks)'])

# disk data
disk_counts = np.array(catalogue['Disk Counts'])
disk_area = np.array(catalogue['Disk Area (px^2)'])
disk_cr = disk_counts/(disk_area * exp_time) # calculating disk count rate (counts/(ks px^2))

# background data
bg_counts = np.array(catalogue['Background Counts'])
bg_area = np.array(catalogue['Background Area (px^2)'])
bg_cr = bg_counts/(bg_area * exp_time) # calculating backgroung count rate (counts/(ks px^2))


net_cr = disk_cr - bg_cr # net count rate

# Sunspot plot comparison with filtered disk count rates (3 panel plot)
fig4 = plt.figure(figsize=(8,7))
ax4 = fig4.add_subplot(311)
loc, label =  plt.xticks()
ax4.plot(time_ssp, ssp, 'k', zorder=10) # plotting sunspot data
major_xticks = np.arange(2000, 2021, 5)
minor_xticks = np.arange(1999, 2022, 1)
plt.ylabel('Sunspot Number')
plt.xlim(1999, max(time_ssp)); plt.ylim(-10, 360)
plt.title('Chandra X-ray Observations')
ax4.set_xticks(major_xticks)
ax4.set_xticks(minor_xticks, minor=True)

major_yticks = np.arange(0, 351, 50)
minor_yticks = np.arange(0, 351, 10)
ax4.set_yticks(major_yticks)
ax4.set_yticks(minor_yticks, minor=True)

ax4.tick_params(which='both', direction='in', top=True, right=True, labelbottom=False)
ax4.text(0.02, 0.95, s='(a)', transform=ax4.transAxes, fontsize=14, va='top')
ax4.vlines(dec_year, ymin=-10, ymax=300, color='b', zorder=0, label='HRC-I') # adding times when HRC-I observations of Jupiter occurred

# 2nd panel - plotting disk and bg count rates
ax5 = fig4.add_subplot(312)

hrc_cr = ax5.scatter(dec_year, disk_cr, color='b', marker='^', zorder=0, label='Jovian Disk') # disk count rate
bg_cr = ax5.scatter(dec_year, bg_cr, color='r', marker='^', zorder=0, label='Background') # background count rate
plt.xlim(1999, max(time_ssp))
ax5.set_xticks(major_xticks)
ax5.set_xticks(minor_xticks, minor=True)
ax5.tick_params(which='both', direction='in', top=True, right=False, labelbottom=False)
ax5.text(0.02, 0.95, s='(b)', transform=ax5.transAxes, fontsize=14, va='top')
plt.ylabel(r'$\mathrm{counts\;ks^{-1}\;px^{-2}}$')
ax5.legend(loc='lower right', fontsize=9, handlelength=1.0, ncol=1)

# 3rd panel - plotting net count rate i.e. disk - bg
ax6 = fig4.add_subplot(313)
net_cr= ax6.scatter(dec_year, net_cr, color='g', marker='^', zorder=0, label='Net count rate')
plt.xlim(1999, max(time_ssp))
ax6.set_xticks(major_xticks)
ax6.set_xticks(minor_xticks, minor=True)
ax6.tick_params(which='both', direction='in', top=True, right=False)
ax6.text(0.02, 0.95, s='(c)', transform=ax6.transAxes, fontsize=14, va='top')
plt.xlabel('Year')
plt.ylabel(r'$\mathrm{counts\;ks^{-1}\;px^{-2}}$')
ax6.legend(loc='upper right', fontsize=9, handlelength=1.0, ncol=1)

plt.subplots_adjust(hspace=0)
plt.savefig('sunspot_plot_3_panel.png', dpi=500, bbox_inches='tight')
plt.show()
