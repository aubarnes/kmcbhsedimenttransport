"""
Hovmoller Diagram of 95 transects from Kāneʻohe Marine Corps Base, Oʻahu, Hawaiʻi
Transects and data from University of Hawaiʻi at Mānoa Coastal Research Collaborative (CRC)
CoastSat and Planet-Scope derived shorelines

Austin Barnes
February 2025

Import necessary libraries
Data directory: "../data/"
Output directory: "../output/"
Read in csv files into separate pandas dataframes:
These are located in the data directory/CRCnorth_beach/
oahu0045_CoastSat_outliersremoved.csv
oahu0045_PlanetScope_outliersremoved_biascorrected.csv
oahu0046_CoastSat_outliersremoved.csv
oahu0046_PlanetScope_outliersremoved_biascorrected.csv

Data columns are: index, datetime, transect name (multiple), satellite name
The datetime column is in format "YYYY-MM-DD HH:MM:SS+00:00" (UTC)
The transect names in the columns are of format "oahu0045_0001" etc.
The data for each transect name is a shoreline position in meters

Create a Hovmoller diagram of the 95 transects
The x-axis is the datetime
The y-axis is the transect name (lowest at the top to highest at the bottom)
The color is the shoreline position in meters (use a colormap with red as negative, blue as positive, and white as zero)
"""
#%% Imports and directories
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
from scipy.stats import linregress

from sklearn.decomposition import PCA

data_dir = "../data/"
output_dir = "../output/"

#%% Read in csv files into separate pandas dataframes
oahu0045_coastsat = pd.read_csv(data_dir + "CRCnorth_beach/oahu0045_CoastSat_outliersremoved.csv")
oahu0045_planetscope = pd.read_csv(data_dir + "CRCnorth_beach/oahu0045_PlanetScope_outliersremoved_biascorrected.csv")
oahu0046_coastsat = pd.read_csv(data_dir + "CRCnorth_beach/oahu0046_CoastSat_outliersremoved.csv")
oahu0046_planetscope = pd.read_csv(data_dir + "CRCnorth_beach/oahu0046_PlanetScope_outliersremoved_biascorrected.csv")

## Drop the first (index) and last (satellite name) columns of the Coastsat dataframes
oahu0045_coastsat = oahu0045_coastsat.drop(columns="Unnamed: 0")
oahu0045_coastsat = oahu0045_coastsat.drop(columns="satname")
oahu0046_coastsat = oahu0046_coastsat.drop(columns="Unnamed: 0")
oahu0046_coastsat = oahu0046_coastsat.drop(columns="satname")

## Make the "dates" column the index
oahu0045_coastsat["dates"] = pd.to_datetime(oahu0045_coastsat["dates"])
oahu0045_coastsat = oahu0045_coastsat.set_index("dates")
oahu0045_planetscope["dates"] = pd.to_datetime(oahu0045_planetscope["dates"], utc=True)
oahu0045_planetscope = oahu0045_planetscope.set_index("dates")
oahu0046_coastsat["dates"] = pd.to_datetime(oahu0046_coastsat["dates"], utc=True)
oahu0046_coastsat = oahu0046_coastsat.set_index("dates")
oahu0046_planetscope["dates"] = pd.to_datetime(oahu0046_planetscope["dates"], utc=True)
oahu0046_planetscope = oahu0046_planetscope.set_index("dates")

## Transpose the dataframes
oahu0045_coastsat = oahu0045_coastsat.T
oahu0045_planetscope = oahu0045_planetscope.T
oahu0046_coastsat = oahu0046_coastsat.T
oahu0046_planetscope = oahu0046_planetscope.T

## Rename the planetscope rows to match the coastsat rows (insert an underscore) in the 5th to last character
oahu0045_planetscope.index = oahu0045_planetscope.index.str.slice_replace(-4, -4, "_")
oahu0046_planetscope.index = oahu0046_planetscope.index.str.slice_replace(-4, -4, "_")

## Compute mean shoreline position as the mean of each transect, combining coastsat and planetscope data
oahu0045_mean = (oahu0045_coastsat.mean(axis=1)*np.shape(oahu0045_coastsat)[0] + oahu0045_planetscope.mean(axis=1)*np.shape(oahu0045_planetscope)[0])/(np.shape(oahu0045_coastsat)[0] + np.shape(oahu0045_planetscope)[0])
oahu0046_mean = (oahu0046_coastsat.mean(axis=1)*np.shape(oahu0046_coastsat)[0] + oahu0046_planetscope.mean(axis=1)*np.shape(oahu0046_planetscope)[0])/(np.shape(oahu0046_coastsat)[0] + np.shape(oahu0046_planetscope)[0])

## Create new dataframes with the shoreline position anomalies
oahu0045_coastsat_demean = oahu0045_coastsat.apply(lambda x: x - oahu0045_mean)
oahu0045_planetscope_demean = oahu0045_planetscope.apply(lambda x: x - oahu0045_mean)
oahu0046_coastsat_demean = oahu0046_coastsat.apply(lambda x: x - oahu0046_mean)
oahu0046_planetscope_demean = oahu0046_planetscope.apply(lambda x: x - oahu0046_mean)

## Transpose the dataframes back for sorting
oahu0045_coastsat_demean = oahu0045_coastsat_demean.T
oahu0045_planetscope_demean = oahu0045_planetscope_demean.T
oahu0046_coastsat_demean = oahu0046_coastsat_demean.T
oahu0046_planetscope_demean = oahu0046_planetscope_demean.T

## Create a combined dataframe of the 95 transects, combining transects from coastsat and planetscope by name, sort by date
oahu0045_combined_demean = pd.concat([oahu0045_coastsat_demean, oahu0045_planetscope_demean], axis=0)
oahu0045_combined_demean = oahu0045_combined_demean.sort_index()
oahu0046_combined_demean = pd.concat([oahu0046_coastsat_demean, oahu0046_planetscope_demean], axis=0)
oahu0046_combined_demean = oahu0046_combined_demean.sort_index()

## Combine the two demeaned dataframes, fill with nans where data is missing. All indices from both dataframes are included
all_combined = oahu0045_combined_demean.combine_first(oahu0046_combined_demean)

## Monthly averages
all_combined_monthly = all_combined.resample("M").mean()

#%% Scatter plot of selected transect
## * for coastsat, . for planetscope
transect1 = "oahu0045_0010"
transect2 = "oahu0046_0056"
%matplotlib qt

fig, ax = plt.subplots(figsize=(12, 8))
# ax.scatter(oahu0045_coastsat_demean.index, oahu0045_coastsat_demean[transect1], marker="*", label="CoastSat")
# ax.scatter(oahu0045_planetscope_demean.index, oahu0045_planetscope_demean[transect1], marker=".", label="PlanetScope")
ax.scatter(oahu0046_coastsat_demean.index, oahu0046_coastsat_demean[transect2], marker="*", label="CoastSat")
ax.scatter(oahu0046_planetscope_demean.index, oahu0046_planetscope_demean[transect2], marker=".", label="PlanetScope")
ax.set_title("Shoreline Posistion Anomalies for Transect " + transect2)
ax.set_xlabel("Date")
ax.set_ylabel("Shoreline Position Anomaly (m)")
ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
ax.grid()
plt.legend(loc = "upper left")
plt.show()


#%% Hovmoller diagram
transects = all_combined.columns
transects = np.insert(transects, np.size(transects), 'S end')

dates = all_combined.index
## Add datetime64('2023-01-01 00:00:00+00:00') to the end of the dates array (for pcolormesh)
dates_extended = np.insert(dates, np.size(dates), pd.to_datetime('2023-01-01 00:00:00+00:00'))

# Create the figure and axis
fig, ax = plt.subplots(figsize=(12, 8))

# Define colormap and normalization
cmap = plt.get_cmap("coolwarm_r")  # Reverse the colormap
# norm = plt.Normalize(-25, 25)
norm = plt.Normalize(-50, 50)

mesh = ax.pcolormesh(dates_extended, transects, all_combined.T, cmap=cmap, norm=norm, shading='auto')

ax.text(pd.to_datetime('2023-01-01 00:00:00+00:00'), 'oahu0045_0001', 'Pyramid Rock', ha='left', va='bottom', fontsize=12)
ax.text(pd.to_datetime('2023-01-01 00:00:00+00:00'), 'S end', 'North Beach', ha='left', va='top', fontsize=12)
ax.text(pd.to_datetime('2023-01-01 00:00:00+00:00'), 'oahu0046_0001', 'Runway', ha='left', va='bottom', fontsize=10)

# Set axis labels and title
ax.set_title("Hovmöller Diagram of Pyramid Rock -> North Beach")
ax.set_xlabel("Date")
ax.set_ylabel("Transect Name")

# Reverse the y-axis direction
ax.invert_yaxis()

# Format x-axis to show dates properly
ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
plt.xticks(rotation=45)

# Add colorbar
cbar = plt.colorbar(mesh, ax=ax)
cbar.set_label("Shoreline Position Anomaly (m)")
# cbar.set_ticks(np.arange(-25, 30, 5))
cbar.set_ticks(np.arange(-50, 55, 10))

# ax.set_xlim([pd.to_datetime('1999-12-31 00:00:00+00:00'), pd.to_datetime('2023-01-01 00:00:00+00:00')])
ax.set_xlim([pd.to_datetime('2014-12-31 00:00:00+00:00'), pd.to_datetime('2023-01-01 00:00:00+00:00')])

plt.tight_layout()
plt.show()
#%% Hovmoller diagram of monthly averages
transects = all_combined_monthly.columns
transects = np.insert(transects, np.size(transects), 'S end')

dates = all_combined_monthly.index
## Add datetime64('2023-01-01 00:00:00+00:00') to the end of the dates array (for pcolormesh)
dates_extended = np.insert(dates, np.size(dates), pd.to_datetime('2023-01-01 00:00:00+00:00'))

# Create the figure and axis
fig, ax = plt.subplots(figsize=(12, 8))

# Define colormap and normalization
cmap = plt.get_cmap("coolwarm_r")  # Reverse the colormap
# norm = plt.Normalize(-25, 25)
norm = plt.Normalize(-50, 50)

mesh = ax.pcolormesh(dates_extended, transects, all_combined_monthly.T, cmap=cmap, norm=norm, shading='auto')

ax.text(pd.to_datetime('2023-01-01 00:00:00+00:00'), 'oahu0045_0001', 'Pyramid Rock', ha='left', va='bottom', fontsize=12)
ax.text(pd.to_datetime('2023-01-01 00:00:00+00:00'), 'S end', 'North Beach', ha='left', va='top', fontsize=12)
ax.text(pd.to_datetime('2023-01-01 00:00:00+00:00'), 'oahu0046_0001', 'Runway', ha='left', va='bottom', fontsize=10)

# Set axis labels and title
ax.set_title("Hovmöller Diagram of Pyramid Rock -> North Beach\nMonthly Averages")
ax.set_xlabel("Date")
ax.set_ylabel("Transect Name")

# Reverse the y-axis direction
ax.invert_yaxis()

# Format x-axis to show dates properly
ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
plt.xticks(rotation=45)

# Add colorbar
cbar = plt.colorbar(mesh, ax=ax)
cbar.set_label("Shoreline Position Anomaly (m)")
# cbar.set_ticks(np.arange(-25, 30, 5))
cbar.set_ticks(np.arange(-50, 55, 10))

# ax.set_xlim([pd.to_datetime('1999-12-31 00:00:00+00:00'), pd.to_datetime('2023-01-01 00:00:00+00:00')])
# ax.set_xlim([pd.to_datetime('2014-12-31 00:00:00+00:00'), pd.to_datetime('2023-01-01 00:00:00+00:00')])
ax.set_xlim([pd.to_datetime('2013-12-31 00:00:00+00:00'), pd.to_datetime('2023-01-01 00:00:00+00:00')])

plt.tight_layout()
plt.show()
#%% Compute EOFs of monthly averages

## Define number of EOFs to compute
n_eof = 4

# Define analysis period
subset_start = '2016-01-01'
subset_end = '2022-12-31'

# Subset the data
subset_data = all_combined_monthly.loc[subset_start:subset_end]

# Interpolate missing values using time-aware linear interpolation
subset_interp = subset_data.interpolate(method='time').ffill().bfill()

# Prepare data matrix
X = subset_interp.values

# Perform PCA/EOF analysis
pca = PCA(n_components=n_eof)
scores = pca.fit_transform(X)

# Calculate explained variance
explained_var = pca.explained_variance_ratio_ * 100

## Calculate scaling factors for physical units
max_abs_scores = np.max(np.abs(scores), axis=0)  # Used for amplitude normalization

# Create scaled variables
eof_amp = scores / max_abs_scores  # Range [-1, 1]
pca_spatialamp = pca.components_ * max_abs_scores[:, np.newaxis]

# Plot configuration
fig, ax = plt.subplots(n_eof, 2, figsize=(14, 12), gridspec_kw={'width_ratios': [3, 2]})

for i in range(n_eof):
    # Spatial pattern plot as bar plot
    bars = ax[i,0].bar(range(len(pca_spatialamp[i])), pca_spatialamp[i], color=['magenta' if val < 0 else 'black' for val in pca_spatialamp[i]])
    ax[i,0].axhline(0, color='black', linestyle='--')
    ax[i,0].axvline(23.5, color='green', linestyle='--') ## Runway Gap
    ax[i,0].axvspan(33, 40, color='lightgray', alpha=0.5) ## Rocky/Reefy Shoreline
    ax[i,0].axvspan(82, 94, color='lightgray', alpha=0.5) ## Rocky/Reefy Shoreline
    if i == 0:
        ax[i,0].text(24, ax[i,0].get_ylim()[1], 'Runway\nGap', ha='left', va='top', color='green')
        ax[i,0].text(33, ax[i,0].get_ylim()[1], 'Rocky/Reefy\nShoreline', ha='left', va='top', color='black')
    ax[i,0].set_title(f'EOF Mode {i+1} ({pca.explained_variance_ratio_[i]*100:.1f}% Variance)')
    ax[i,0].set_xlabel('Transect Number (West to East)')
    ax[i,0].set_ylabel('Spatial Amplitude (m)')
    
    # Time series plot
    ax[i,1].plot(subset_interp.index, eof_amp[:,i], color='darkred')
    ax[i,1].axhline(0, color='black', linestyle='--')
    ax[i,1].set_title(f'PC {i+1} Amplitude')
    ax[i,1].set_xlabel('Date')
    ax[i,1].set_ylabel('Amplitude')
    ax[i,1].set_ylim([-1, 1])
    
    ax[i,1].xaxis.set_major_locator(mdates.YearLocator())
    ax[i,1].xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    ax[i,1].grid(True, which='both', axis='x', linestyle='--', color='gray')

plt.tight_layout()
plt.show()
#%% Monthly Time Series of Spatially Averaged Width

## Transpose the dataframes
oahu0045_coastsat = oahu0045_coastsat.T
oahu0045_planetscope = oahu0045_planetscope.T
oahu0046_coastsat = oahu0046_coastsat.T
oahu0046_planetscope = oahu0046_planetscope.T

oahu0045_combined = pd.concat([oahu0045_coastsat, oahu0045_planetscope], axis=0)
oahu0045_combined = oahu0045_combined.sort_index()
oahu0046_combined = pd.concat([oahu0046_coastsat, oahu0046_planetscope], axis=0)
oahu0046_combined = oahu0046_combined.sort_index()

all_combined_shoreline = oahu0045_combined.combine_first(oahu0046_combined)

## Monthly average, then spatial average
shorelines_monthly = all_combined_shoreline.resample("M").mean()
shoreline_monthly = shorelines_monthly.mean(axis=1)

shoreline_monthly_subset = shoreline_monthly.loc['2016-01-01':'2022-12-31']

#%% Linear fit of shoreline_monthly from 2016-2023
## Perform linear regression
slope, intercept, r_value, p_value, std_err = linregress(np.arange(len(shoreline_monthly_subset)), shoreline_monthly_subset)
#%% Plot the ave_monthly and shoreline_monthly
fig, ax = plt.subplots(figsize=(12, 8))
ax.plot(shoreline_monthly_subset, label="Spatially Averaged Shoreline Position")
ax.plot(shoreline_monthly_subset.index, slope*np.arange(len(shoreline_monthly_subset)) + intercept, label="Linear Fit", linestyle="--")

ax.set_title("Monthly Time Series of Spatially Averaged Width")
ax.set_xlabel("Date")

ax.xaxis.set_major_locator(mdates.YearLocator())
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
ax.grid(True, which='both', linestyle='--', color='gray')
plt.legend()

plt.tight_layout()
plt.show()
#%% Monthly Time Series of West vs. East side average width
# east_west_div = 'oahu0046_0001' ## Runway
east_west_div = 'oahu0046_0028' ## Near Bill's Rocking Point
east_west_div_idx = np.where(all_combined_shoreline.columns == east_west_div)[0][0]
# pyramid_all = all_combined_shoreline.loc[:, :east_west_div]
# nb_all = all_combined_shoreline.loc[:, east_west_div:]

pyramid_all = all_combined_shoreline.iloc[:, :east_west_div_idx]
nb_all = all_combined_shoreline.iloc[:, east_west_div_idx:]

pyramid_all_monthly = pyramid_all.resample("M").mean()
nb_all_monthly = nb_all.resample("M").mean()

pyramid_ave_monthly = pyramid_all_monthly.mean(axis=1)
nb_ave_monthly = nb_all_monthly.mean(axis=1)

pyramid_ave_monthly_subset = pyramid_ave_monthly.loc['2016-01-01':'2022-12-31']
## Interpolate missing values
pyramid_ave_monthly_subset = pyramid_ave_monthly_subset.interpolate(method='time').ffill().bfill()
nb_ave_monthly_subset = nb_ave_monthly.loc['2016-01-01':'2022-12-31']

## Perform linear regression
slope_pyramid, intercept_pyramid, r_value_pyramid, p_value_pyramid, std_err_pyramid = linregress(np.arange(len(pyramid_ave_monthly_subset)), pyramid_ave_monthly_subset)
slope_nb, intercept_nb, r_value_nb, p_value_nb, std_err_nb = linregress(np.arange(len(nb_ave_monthly_subset)), nb_ave_monthly_subset)

#%% Plot Pyramid and NB ave_monthly with linear fits
fig, ax = plt.subplots(figsize=(12, 8))
ax.plot(pyramid_ave_monthly_subset, label="Pyramid Rock")
ax.plot(nb_ave_monthly_subset, label="North Beach")
ax.plot(pyramid_ave_monthly_subset.index, slope_pyramid*np.arange(len(pyramid_ave_monthly_subset)) + intercept_pyramid, label="Pyramid Rock Linear Fit", linestyle="--")
ax.plot(nb_ave_monthly_subset.index, slope_nb*np.arange(len(nb_ave_monthly_subset)) + intercept_nb, label="North Beach Linear Fit", linestyle="--")

ax.set_title("Monthly Time Series of West vs. East Side Average Width")
ax.set_xlabel("Date")

ax.xaxis.set_major_locator(mdates.YearLocator())
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
ax.grid(True, which='both', linestyle='--', color='gray')
plt.legend()

plt.tight_layout()
plt.show()
#%% Take subset of shorelines_monthly, find linear trend for each transect
shorelines_monthly_subset = shorelines_monthly.loc['2016-01-01':'2022-12-31']

## Interpolate missing values
shorelines_monthly_subset = shorelines_monthly_subset.interpolate(method='time').ffill().bfill()

## Perform linear regression for each transect
slope_all = np.zeros(np.shape(shorelines_monthly_subset)[1])
intercept_all = np.zeros(np.shape(shorelines_monthly_subset)[1])
r_value_all = np.zeros(np.shape(shorelines_monthly_subset)[1])
p_value_all = np.zeros(np.shape(shorelines_monthly_subset)[1])
std_err_all = np.zeros(np.shape(shorelines_monthly_subset)[1])

for i in range(np.shape(shorelines_monthly_subset)[1]):
    slope_all[i], intercept_all[i], r_value_all[i], p_value_all[i], std_err_all[i] = linregress(np.arange(len(shorelines_monthly_subset)), shorelines_monthly_subset.iloc[:, i])

#%% Plot the linear trends for each transect
fig, ax = plt.subplots(figsize=(12, 8))
ax.bar(range(len(slope_all)), slope_all*12, color=['magenta' if val < 0 else 'black' for val in slope_all])
ax.axhline(0, color='black', linestyle='--')
ax.axvline(23.5, color='green', linestyle='--') ## Runway Gap
ax.axvspan(33, 40, color='lightgray', alpha=0.5) ## Rocky/Reefy Shoreline
ax.axvspan(82, 94, color='lightgray', alpha=0.5) ## Rocky/Reefy Shoreline
ax.text(24, ax.get_ylim()[1], 'Runway\nGap', ha='left', va='top', color='green')
ax.text(33, ax.get_ylim()[1], 'Rocky/Reefy\nShoreline', ha='left', va='top', color='black')
ax.set_title("Shoreline Position Trends, Jan 2016 - Dec 2022")
ax.set_xlabel("Transect Number (West to East)")
ax.set_ylabel("Shoreline Position\nLinear Trend (m/year)")

plt.tight_layout()
plt.show()

# %%
