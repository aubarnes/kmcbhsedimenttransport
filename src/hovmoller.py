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
# %%
