'''
Function to calculate beach year (Oct 1 - Sep 30)
seasonally weighted annual means of shoreline position anomaly

Inputs:
    xdatetime: numpy array of datetime indices of satellite shoreline position observations
    xpos: numpy array of data values

Outputs:
    beach_years: numpy array of beach years
    annual_means: numpy array of seasonally weighted annual means
    sparse: numpy array of sparse data flags

Notes:
    - The function calculates the seasonally weighted annual means for each beach year in the input data.
    - The function uses the following algorithm:
        1. Calculate the mean of the last 3 months of the previous year and the first 3 months of the current year.
        2. Calculate the mean of the second quarter of the current year.
        3. Calculate the mean of the third quarter of the current year.
        4. Calculate the mean of the fourth quarter of the current year.
        5. Calculate the weighted mean of the second and fourth quarters.
        6. Calculate the weighted mean of the third and fourth quarters.
        7. Calculate the mean of the weighted means.
    - The function handles missing data by:
        - Filling missing data with the next available data point.
        - Adjusting the data points to ensure consistency.
        - Making sparse data adjustments.
    - The function returns the beach years, seasonally weighted annual means, and a sparse flag indicating if sparse data adjustments were made.

Reference:
Adapted from Matlab code written by William C. O'Reilly, Scripps Institution of Oceanography, UCSD
Date received from author: 2025-03-12
Adapted by Austin Barnes
'''

import numpy as np
import pandas as pd
from datetime import datetime, timedelta

def get_beach_year_seasonally_weighted_annual_means(xdatetime, xpos):
    vflag = 0  # verbose flag; 1 = print details to screen; 0 = quiet

    # Calculate beach years
    beach_years = np.unique(pd.to_datetime(xdatetime) + pd.DateOffset(months=3)).year
    beach_years = np.arange(beach_years[0], beach_years[-1] + 1)

    annual_means = np.full(len(beach_years), np.nan)
    avg_month_count = np.zeros(len(beach_years))
    sparse = np.zeros(len(beach_years))

    # Convert xdatetime to pandas datetime
    xdatetime = pd.to_datetime(xdatetime)

    # Calculate year-month means
    df = pd.DataFrame({'datetime': xdatetime, 'pos': xpos})
    df['year'] = df['datetime'].dt.year
    df['month'] = df['datetime'].dt.month
    ymmean = df.groupby(['year', 'month'])['pos'].mean().unstack()
    ymcount = df.groupby(['year', 'month'])['pos'].count().unstack()

    yrs = np.arange(ymmean.index[0], ymmean.index[-1] + 1)

    for y in yrs[1:]:
        iy = np.where(yrs == y)[0][0]
        ipy = np.where(yrs == y - 1)[0][0] if y - 1 in yrs else None
        iny = np.where(yrs == y + 1)[0][0] if y + 1 in yrs else None

        q1 = q2 = q3 = q4 = wmean = smean = np.nan

        if ipy is not None:
            if vflag == 1:
                print(f"{y} {' '.join([f'{x:4.1f}' for x in ymmean.loc[y-1, 10:12].tolist() + ymmean.loc[y, 1:10].tolist()])}")

            q1 = np.nanmean(ymmean.loc[y-1, 10:12])
            q2 = np.nanmean(ymmean.loc[y, 1:3])
            q3 = np.nanmean(ymmean.loc[y, 4:6])
            q4 = np.nanmean(ymmean.loc[y, 7:9])

            if np.isnan([q1, q2, q3, q4]).sum() > 0:
                idx = np.where(beach_years == y)[0][0]
                sparse[idx] = 1
                if vflag == 1:
                    print(f"{y} Making Sparse Data Adjustment of these 4 Qs: {q1:4.1f} {q2:4.1f} {q3:4.1f} {q4:4.1f}")

                if np.isnan(q4) and not np.isnan(ymmean.loc[y, 10]):
                    q4 = ymmean.loc[y, 10]
                    q1 = np.nanmean(ymmean.loc[y-1, 11:12])
                    if vflag == 1:
                        print(f"{y} 1. filling q4 using next Oct data {q1:3.1f} {q2:3.1f} {q3:3.1f} {q4:3.1f}")

                if np.isnan(q2) and not np.isnan(ymmean.loc[y, 4]):
                    q2 = ymmean.loc[y, 4]
                    q3 = np.nanmean(ymmean.loc[y, 5:6])
                    if vflag == 1:
                        print(f"{y} 2. adjusting q2 to use Apr data {q1:3.1f} {q2:3.1f} {q3:3.1f} {q4:3.1f}")

                if not np.isnan(q3) and q4 >= q3:
                    q3 = np.nan
                elif q3 > q4:
                    q4 = np.nan

                if not np.isnan(q1) and q2 <= q1:
                    q1 = np.nan
                elif q1 < q2:
                    q2 = np.nan

                wmean = np.nanmin([q1, q2])
                smean = np.nanmax([q3, q4])
            else:
                wmean = np.nanmean([q1, q2])
                smean = np.nanmean([q3, q4])

        if vflag == 1:
            print(f"{y} Q {q1:4.1f} {q2:4.1f} {q3:4.1f} {q4:4.1f}")
            print(f"{y} S {wmean:4.1f} {smean:4.1f}")

        ymean = np.nanmean([wmean, smean])
        idx = np.where(beach_years == y)[0][0]
        if idx is not None:
            annual_means[idx] = ymean

        if vflag == 1:
            print(f"{y} Y {ymean:4.1f}")

    return beach_years, annual_means, sparse

# Example usage:
# xdatetime = pd.date_range(start='2000-01-01', end='2022-12-31', freq='D')
# xpos = np.random.rand(len(xdatetime)) * 100
# beach_years, annual_means, sparse = get_beach_year_seasonally_weighted_annual_means(xdatetime, xpos)
