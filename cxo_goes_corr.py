# Script that calulates the Pearson's Correlation Coefficient (PCC) between the net jovian disk light curves and the median GOES solar X-ray flux. Also produces a plot of the correlation as an output.

# relevant packages
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

# Reading in the config.ini file containing any hard wired inputs
import configparser
config = configparser.ConfigParser()
config.read('config.ini')

# Reading in excel catalogue file with all relevant info
catalogue = pd.read_excel('catalogue_all_data.xlsx')

obsIDs = np.array(catalogue['ObsID'])

start_date = np.array(catalogue['Start Date'])
exp_time = np.array(catalogue['Exposure Time (ks)'])

# disk data
disk_counts = np.array(catalogue['Disk Counts'])
disk_area = np.array(catalogue['Disk Area (px^2)'])
disk_cr = disk_counts/(disk_area * exp_time) # calculating disk count rate (counts/(ks px^2))

# background data
bg_counts = np.array(catalogue['Background Counts'])
bg_area = np.array(catalogue['Background Area (px^2)'])
bg_cr = bg_counts/(bg_area * exp_time) # calculating backgroung count rate (counts/(ks px^2))


net_cr = disk_cr - bg_cr # net jovian disk count rate

# read in median goes flux data
goes_data = pd.read_csv('median_flux_goes.txt')
med_flux_goes = np.array(goes_data['# flux (W/m^2)'])

# make correlation figure
plt.figure()
plt.scatter(med_flux_goes, net_cr, color='k', s=5, zorder=10)

# calculating PCC
from scipy.stats import pearsonr
corr, _ = pearsonr(med_flux_goes, net_cr)
m_corr, b_corr = np.polyfit(med_flux_goes, net_cr, 1)

# plotting line of best fit
X_corr = np.linspace(min(med_flux_goes), max(med_flux_goes), 1000)
corr_line, = plt.plot(X_corr, m_corr * X_corr + b_corr, 'b', zorder=0, label=f'y = {int(np.round(m_corr))}x + {np.round(b_corr, 2)}')
# plt.legend(handles=[corr_line], loc='upper left')
plt.title(f"""Net Disk Count Rate vs Median Solar X-ray Flux \nPearson's Correlation: {np.round(corr, 2)}""")
plt.ylabel(r'Net disk count rate ($\mathrm{counts\;min^{-1}\;px^{-2}}$)')
plt.xlabel(r'Solar X-ray flux ($\mathrm{W m^{-2}}$)')
plt.xscale('log')
plt.tight_layout()
plt.savefig('cxo_goes_corr.png', dpi=500)
plt.show()

