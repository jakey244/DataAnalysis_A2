WRITEUP
=======

PART 1:
-------
In this part of the assignment, the geostrophic ocean transport at 26°N is computed from temperature
and salinity profiles at the boundaries.

For this, the profiles are first converted to TEOS-10, then the dynamic height is calculated for each level.
Using the thermal wind relation, the vertical shear in meridional transport is calculated.
The thermal wind profile is then integrated from the surface down towards a maximum depth of 1100 meters.
We do not do a full depth integral because the level of 1100 meters is considered the "upper limb" of the AMOC which we want to isolate here.
This value is chosen because it corresponds to the depth of the AMOC maximum.
As a reference level for the dynamic height, 4820 dBar is chosen as the level of no motion.

Figure UMO_geostrophy.png shows the calculated time series from geostrophy compared to the published Upper Mid-Ocean transport from RAPID.
To make the time series comparable, the vertical profiles for the geostrophy calculation were filtered using a 10-day lowpass Tukey filter.

The mean of the two timeseries is different by about 6.25 Sv (-12.11 for geostrophy, -18.36 for UMO).
The standard deviation of the two timeseries is different by about 2.75 (6.16 for geostrophy, 3.41 for UMO).

This difference arises from the fact that the UMO transport includes more data, not just the boundary profiles.
Its calculation is more complex e.g. it splits the interior at the Mid-Atlantic Ridge and applies a mass-balance adjustment.
The geostrophy calculation is just the interior geostrophic transport, the UMO includes a barotropic compensation term
such that the net transport is mass balanced (when also considering Ekman and Florida Current).
The extra variance comes from the compensation term, as strong deviations from the mean in the interior geostrophic transport
are balanced by the compensation and therefore do not show up in the published UMO as strongly.
(doi.org/10.1029/2021GL093045)


For the unfiltered geostrophy time series, the correlation coefficient is r = 0.751.
For the 10-day lowpassed series, the correlation coefficient is r = 0.773.

Figure scatter_geostrophy_UMO.png shows a scatter plot for each data point in the timeseries.
A linear fit is calculated to show the correlated relationship and the confidence interval at 95% is shown.
The relationship is statistically significant at 95%.


PART 2A:
-------
Figures MOC_seasonal_cycle.png and UMO_seasonal_cycle.png show the mean monthly seasonal cycle for the UMO and the MOC transports.

Figures MOC_seasonal_cycle_comparison.png and UMO_seasonal_cycle_comparison.png show the original timeseries and the deseasonalized timeseries
overlayed for both transport datasets.

Comparing the two figures, it can be seen that for the UMO, the deseasonalized timeseries deviates more strongly from the original data than for MOC.

Both the original and deseasonalized timeseries are regridded on a monthly grid (using monthly means), as the seasonality is calculated on such a grid as well.

For both timeseries, we calculate the autocorrelation which is shown in Figure MOC_UMO_autocorrelation.png.
We furthermore calculate the effective degrees of freedom, using the integral timescale definition of the decorrelation timescale.
The data that is put into the autocorrelation is detrended but not deseasonalized.
For MOC, the decorrelation timescale is 18.5 days, for UMO it is 22.5 days.
Using the decorrelation timescale, the number of effective degrees of freedom can be calculated. For MOC, EDOF = 196. For UMO, EDOF = 161.
The effective degrees of freedom are rounded down to the nearest integer.



PART 2B:
-------
Hallo

PART 3:
-------
Hallo
