# PyOAE
Python tools for analysis of Ocean Alkalinity Enhancement and Ocean Acidification

by Greg Pelletier

Tools for analysis of Ocean Alkalinity Enhancement (OAE) and Ocean Acidification (OA), including the following functions and examples:

- **f_dTA** - Root-finding method to solve for the OAE treatment needed to restore any carbonate system variable (e.g. TA-DIC, CO3--, pH, Ωara, Ωcal) to pre-industrial conditions in the global oceans or any selected region.
- **etamax** - Calculation of the maximum hypothetical OAE efficiency ηmax (etamax) for any assumed addition of alkalinity. The ηmax is a dimensionless quantity that is the hypothetical maximum potential CDR (umol/kg) divided by the amount of added alkalinity (umol/kg).
- **dic_bio** - Calcuation of the biological component of DIC. The observed surface ocean DIC concentration (DIC_obs) is influenced by air–sea gas exchange, the biological production/consumption of organic matter, and calcium carbonate (CaCO3) formation/dissolution (Burt et al, 2016). To isolate the biological component of DIC (DIC_bio), a surface DIC concentration at atmospheric equilibrium (DIC_atm) is computed and subsequently removed from the observed DIC (DIC_obs), such that DIC_bio = DIC_obs - DIC_atm
- **sine_fit** - A sine-regression function to model time-series of variables with periodic seasonal cycles (e.g. DIC_bio)
- **pco2_tnorm** and **pco2_fass** - Temperature-normalization of pCO2 in seawater using the equations of Takahashi (2002), or Fassbender's method as described in Rodgers et al (2022). The thermal and non-thermal compoments of observed pCO2 are estimated, including an estimate of the seasonal cycle anomaly of the thermal component.

The PyOAE package requires that you have already installed numpy, scipy, and PyCO2SYS packages. We also recommend that you have installed xarray, cartopy, and matplotlib to analyze and plot maps using data from netcdf files. We also recommend the use of the multiprocessing package to apply the root-finding method to large datasets.

# Installation for Python or Jupyter Notebook

If you have not already installed PyOAE, enter the following with pip or !pip in your notebook or terminal:<br>
```
pip install git+https://github.com/gjpelletier/PyOAE.git
```

if you are upgrading from a previous installation of PyOAE, enter the following with pip pr !pip in your notebook or terminal:<br>
```
pip install git+https://github.com/gjpelletier/PyOAE.git --upgrade
```

Next import the PyOAE functions as follows in your notebook or python code:<br>
```
from PyOAE import f_dTA, etamax, dic_bio, sine_fit, pco2_tnorm
```

As an alternative to the commands above, you can download PyOAE.py from this github repository (https://github.com/gjpelletier/PyOAE) and copy and paste the functions into your own project.<br>

# Geospatial Data Analysis Project

#### *PyOAE - Geospatial analysis of the potential mitigation of ocean acidification (OA) using ocean alkalinity enhancement (OAE)*

**Calculation of the amount of OAE addition of TA needed to restore current ocean conditions to pre-industrial for OA indicator variables in the global oceans**

by Greg Pelletier, 16-Mar-2025

The "PyOAE_example_root_finding_multiprocessing.ipynb" notebook presents a new method to calculate the amount of Ocean Alkalinity Enhancement (OAE) addition of total alkalinity (TA) that is needed to restore the current untreated surface ocean condtions (e.g. in year 2010) to the same conditions that were present in pre-industrial times (e.g. year 1750), for a selected Ocean Acidification (OA) indicator variable (e.g. TA-DIC, CO3--, Ωara, Ωcal, or pH).   

**Background**

Anthropogenic emissions have significantly increased the atmospheric carbon dioxide (CO2) concentration, from a pre-industrial concentration of around 278 ppm, to over 420 ppm today (Lan et al 2023). Uptake of CO2 by the ocean leads to ocean acidification (OA). OA has already resulted in negative effects on marine ecosystems and organisms, especially marine calcifiers (Bednarsek et al 2023, Cornwall et al 2024). To mitigate the negative effects of global climate change, carbon dioxide removal (CDR) strategies, including marine carbon dioxide removal (mCDR), are considered to be necessary (Lee et al 2023). One important mCDR solution is ocean alkalinity enhancement (OAE), where alkalinity is added to the ocean to increase absorbtion of CO2 from the atmosphere into the ocean. OAE involve adding alkaline substances, such as sodium hydroxide (NaOH) or sodium carbonate (Na2CO3), to the ocean.

While the primary objective of OAE is to increase the uptake of CO2 from atmosphere into the ocean, a secondary benefit could be some amount of mitigation of OA. However, there is an inverse relationship between the OAE efficiency of increasing the uptake of CO2 from the atmosphere and the magntiude of possible OA mitigation.

Analysis of the response of OA indicator variables may be used to evaluate the potential for mitigation of OA by OAE. The difference between TA and DIC, also known as Alk* (Sarmiento & Gruber, 2006), can be used as a surrogate variable to interpret the response of other carbonate system variables such as carbonate ion concentrations (CO3--), pH, aragonite saturation state (Ωara), and calcite saturation state (Ωcal). Each of these variables can be used to interpret the response of calcification of sensitive marine organisms to to ocean acidification (OA) and ocean alkalinity enhancement (OAE) (Xue & Cai, 2020). Ocean acidification has caused a decreases in TA-DIC, CO3--, pH, Ωara, Ωcal since pre-industrial conditions (Sarmiento & Gruber, 2006), with corresponding decrease in calcification rates and increases in dissolution of sensitive organisms (e.g. Bednarsek et al 2025). 

**Problem statement**

The questions we attempt to address in ths notebook are stated as follows: 

- What is the quantity of OAE required to restore the recent historical OA indicators (e.g. TA-DIC, CO3--, pH, Ωara, Ωcal) to pre-industrial conditions?
- How does the quantity of OAE needed to mitigate OA vary spatially?
- How does the quantity of OAE needed to mitigate OA vary depending on which OA indicator variable is evaluated? 

**Methods and data set**

In this notebook we present a root-finding method to solve for the OAE treatment to restore any carbonate system variable (e.g. TA-DIC, CO3--, pH, Ωar, etc.) to pre-industrial conditions. The advantage of the root-finding method is that it can be used to solve for the OAE treatment required for restoration of any carbonate system variable with excellent precision. The disadvantage of the root-finding method is that it takes several hours of CPU time to solve for a water layer over the entire global ocean grid at a spatial resolution of 1$^\circ$ x 1$^\circ$.

Jiang et al (2023) published a global data set of historical and predicted future values of ocean acidification indicator variables for the surface layer of the ocean between the years 1750 to 2100 (https://doi.org/10.1029/2022MS003563). We will use a subset of those data from the years 1750 and 2010 to estimate the amount of total alkalinity (TA) that would need to be added in order to restore the conditions in 2010 to pre-industrial condtions (1750). 

A "root-finding method" is a mathematical technique used to determine the values of a variable (e.g. the OAE treatment ΔTA) that make a given function equal to zero, essentially finding the "root" or "zero" of that function. This is done through an iterative numerical algorithm when analytical solutions are not readily available. In this project we use Brent’s method for the numerical algorithm as implemented in the scipy optimize brentq function.

The general formula for the root-finding method with brentq is as follows, as brentq tries to find a point x where fun(x) = 0 using Brent’s method:

x = scipy.optimize.brentq(fun, x_lwr, x_upr)		

where

-  x = the root of the OAE treatment ΔTA (µmol kg-1) such that the the result of the function ‘fun’ is equal to zero
- fun = the function that calculates the difference between the treated condition at time t (e.g. year 2010 afer hypothetical OAE treatment) compared with the pre-industrial condition (e.g. year 1750) for whichever objective variable is being evaluated (e.g. TA-DIC, CO3--, pH, Ωar).
-  x_lwr = lower bound range of the OAE treatment ΔTA (µmol kg-1).
-  x_upr = upper bound range of the OAE treatment ΔTA (µmol kg-1).

The function ‘fun’ determines the difference between the treated condition at time t (e.g. year 2010 after hypothetical OAE treatment) compared with the pre-industrial condition (e.g. year 1750) for whichever objective variable is being evaluated (e.g. TA-DIC, CO3--, pH, Ωara, Ωcal) as follows:

fun(x) = y_t,trt – y_PI	

where

-  y_t,trt = the value objective variable at time t (e.g. 2010) after hypothetical OAE treatment
-  y_PI = the value of the objective variable in pre-industrial conditions (e.g. 1750)
-  y = the carbonate system variable that is used as the objective variable (e.g. TA-DIC, CO3--, pH, Ωara, Ωcal).

The function ‘fun’ performs the following calculations using a carbonate system calculator (e,g, PyCO2SYS) and the following equations:

-	Calculate ΔDICtrt and ΔDICcdr from trial values of x = dTA = ΔTAtrt using the following equations:

    - ΔTAtrt = dTA = x = the trial value of the hypothetical OAE treatment addition of total alkalinity (umol/kg)
    - ΔDICtrt = direct chemical addition of DIC due to the trial value of the OAE treatment addition (umol/kg)
    - ΔDICcdr = indirect CDR addition of DIC in response to the OAE treatment addition (umol/kg)
    - TAt,ctl = control TA at time t before OAE treatment (umol/kg)
    - TAt,trt = TAt,ctl + dTA = treatment TA at time t after OAE treatment (umol/kg)
    - DICt,ctl = control DIC at time t before OAE treatment (umol/kg)
    - DICt,trt = treatment DIC at time t after OAE treatment (umol/kg)
    - pCO2t,ctl = control pCO2 at time t before OAE treatment (uatm)

    - DICeq = CO2SYS calculation of DIC at equilibrium when input TA = TAt,trt and input pCO2 = pCO2t,ctl
    - CDRpot = DICeq - DICt,ctl = hypothetical maximum potential CDR (umol/kg)
    - etamax = CDRpot / dTA = hypothetical maximum OAE efficiency, typically in the range of 0.7-0.9 (Yankovsky et al 2024) (dimensionless)
    - eta(t) = CDR / dTA = realized OAE efficiency at  time t due to the realized CDR (dimensionless)
    - CDReff = eta(t) / etamax = realized CDR efficiency, typically in the range of 0.65-0.95 (Zhou et al 2024, Wang et al 2023) (dimensionless)  

    - If NaOH is used for OAE, then the following equations apply:
        - ΔDICtrt = 0
        - ΔDICcdr = CDReff * etamax * dTA
    - If Na2CO3 is used for OAE, then the following equations apply:
        - ΔDICtrt = 0.5 * dTA
        - ΔDICcdr = cdreff * (etamax - 0.5) * dTA

-	Calculate the DICt,trt and TAt,trt at time t (2010) after OAE treatment as follows:
    - TAt,trt = TAt,ctl + dTA
    - DICt,trt = DICt,ctl + ΔDICtrt + ΔDICcdr
-	Calculate yt,trt corresponding to DICt,trt and TAt,trt
-	Calculate yPI corresponding to DICPI and TAPI
-	Calculate fun(x) = yt,trt – yPI


**References**

- Bednaršek, N., B. R. Carter, R. M. McCabe, R. A. Feely, E. Howard, F. P. Chavez, M. Elliott, J. L. Fisher, J. Jahncke, Z. Siegrist, Pelagic calcifiers face increased mortality and habitat loss with warming and ocean acidification. Ecol. Appl. 32, e2674 (2022).
- Bednaršek, N., H. van de Mortel, G. Pelletier, M. García-Reyes, R. A. Feely, A. G. Dickson, Assessment framework to predict sensitivity of marine calcifiers to ocean alkalinity enhancement – identification of biological thresholds and importance of precautionary principle. Biogeosciences 22, 473–498 (2025).
- Cornwall, C.E., S. Comeau, B. P. Harvey, Are physiological and ecosystem-level tipping points caused by ocean acidification? A critical evaluation. Earth Syst. Dyn. 15, 671–687 (2024).
- Jiang, L., J. Dunne, B. R. Carter, J. F. Tjiputra, J. Terhaar, J. D. Sharp, A. Olsen, S. Alin, D. C. E. Bakker, R. A. Feely, J. Gattuso, P. Hogan, T. Ilyina, N. Lange, S. K. Lauvset, E. R. Lewis, T. Lovato, J. Palmieri, Y. Santana‐Falcón, J. Schwinger, R. Séférian, G. Strand, N. Swart, T. Tanhua, H. Tsujino, R. Wanninkhof, M. Watanabe, A. Yamamoto, T. Ziehn, Global Surface Ocean Acidification Indicators From 1750 to 2100. J. Adv. Model. Earth Syst. 15, e2022MS003563 (2023).
- Lan, X., P. Tans, K. Thoning, NOAA Global Monitoring Laboratory, Trends in globally-averaged CO2 determined from NOAA Global Monitoring Laboratory measurements., NOAA GML (2023); https://doi.org/10.15138/9N0H-ZH07.
- Lee et al 2023, “IPCC, 2023: Climate Change 2023: Synthesis Report. Contribution of Working Groups I, II and III to the Sixth Assessment Report of the Intergovernmental Panel on Climate Change [Core Writing Team, H. Lee and J. Romero (eds.)]. IPCC, Geneva, Switzerland.” (Intergovernmental Panel on Climate Change (IPCC), 2023); https://doi.org/10.59327/IPCC/AR6-9789291691647.
- Sarmiento. J.L., N. Gruber, Ocean Biogeochemical Dynamics (Princeton University Press, 2006; https://www.jstor.org/stable/j.ctt3fgxqx).
- Wang, H., D. J. Pilcher, K. A. Kearney, J. N. Cross, O. M. Shugart, M. D. Eisaman, B. R. Carter, Simulated Impact of Ocean Alkalinity Enhancement on Atmospheric CO2 Removal in the Bering Sea. Earths Future 11, e2022EF002816 (2023).
- Xue, L., W.-J. Cai, Total alkalinity minus dissolved inorganic carbon as a proxy for deciphering ocean acidification mechanisms. Mar. Chem. 222, 103791 (2020).
- Yankovsky, E., M. Zhou, M. Tyka, S. Bachman, D. Ho, A. Karspeck, M. Long, Impulse response functions as a framework for quantifying ocean-based carbon dioxide removal. EGUsphere, 1–26 (2024).
- Zhou, M., M. D. Tyka, D. T. Ho, E. Yankovsky, S. Bachman, T. Nicholas, A. R. Karspeck, M. C. Long, Mapping the global variation in the efficiency of ocean alkalinity enhancement for carbon dioxide removal. Nat. Clim. Change, 1–7 (2024)


# Example use of the root-finding method to analyze the OAE needed to restore OA indicators to pre-industrial conditions

The difference between TA and DIC, also known as Alk* (Sarmiento & Gruber, 2006), can be used as a surrogate variable to interpret the response of other carbonate system variables (e.g. CO3--, pH, Ωara, Ωcal). Ocean acidification (OA) has caused a decrease in TA-DIC since pre-industrial conditions (Sarmiento & Gruber, 2006).

In this example we will use the root-finding method to solve for the amount of OAE needed to restore the TA-DIC in 2010 in the global oceans to pre-industrial conditions. The example netcdf file "Jiang_data_for_PyOAE.nc" is available to download from this github repository. This data set provides the decadal average global carbonate system variable values in the year 2010 and in the year 1750. The original source of these data is from Jiang et al (2023) (https://doi.org/10.1029/2022MS003563)

The current version of f_dTA analyzes the ocean chemistry data from one grid cell at a time. Processing each grid cell takes about 1-3 seconds. To analyze all of the grid cells in a model domain, or a subset for a region of selected grid cells, we use the multiprocessing package in Python to loop through all of the grid cells that need to be evaluated, and solve for the root in each grid cell one at a time in the loop. 

The following code solves for and maps the amount of OAE needed to restore the TA-DIC to pre-industrial conditions, assuming that NaOH is used, and the CDR efficiency is 80%. This code takes about 2.5 hours to run on a laptop using 8 CPU cores:<br>
```
# import the packages that are needed
import scipy.optimize as opt
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from PyOAE import f_dTA
import time
import multiprocessing
import scipy.io
# read the global arrays of surface ocean data and assign to a dictionary
ds = xr.open_dataset("Jiang_data_for_PyOAE.nc", chunks={"x":60})
# copy the dask 2d arrays in ds values to numpy 2d arrays in ds_dict
ds_dict = {var: ds[var].values for var in ds.data_vars}
# Reshape all 2D arrays to 1D
ds_dict_1d = {key: arr.flatten() for key, arr in ds_dict.items()}
# initialize 1d array of the dTA_root values
ds_dict_1d["dTA_root"] = np.full((180*360), np.nan) # init out array 
# specify the options to use for this example
obj_var = 'alkstar'   # objective variable: 'alkstar', 'co3', phtot', 'omara', 'omcal'
oae_type = 'NaOH'     # chemical used for OAE: 'NaOH' or 'Na2CO3'
cdreff = 0.8          # CDR efficiency between 0-1 (e.g. 0.8 = 80%)
x_lwr = 0             # lower bound of possible dTA values (umol/kg)
x_upr = 500           # upper bound of possible dTA values (umol/kg)
# define the function to find the root in each grid cell i
def find_root(i):
    chem_pi = np.full(7, np.nan)
    chem_pi[0] = ds_dict_1d["talk_1750"][i]    # TA in 1750 (umol/kg)
    chem_pi[1] = ds_dict_1d["dic_1750"][i]     # DIC in 1750 (umol/kg)
    chem_pi[2] = ds_dict_1d["sio3"][i]         # SiO3 in 1750 (umol/kg)
    chem_pi[3] = ds_dict_1d["po4"][i]          # PO4 in 1750 (umol/kg)
    chem_pi[4] = ds_dict_1d["temp_1750"][i]    # Temperature in 1750 (degC)
    chem_pi[5] = ds_dict_1d["sal_1750"][i]     # Salinity in 1750 (psu)        
    chem_pi[6] = 0
    chem_ctl = np.full(7, np.nan)
    chem_ctl[0] = ds_dict_1d["talk_2010"][i]   # TA in 2010 (umol/kg)
    chem_ctl[1] = ds_dict_1d["dic_2010"][i]    # DIC in 2010 (umol/kg)
    chem_ctl[2] = ds_dict_1d["sio3"][i]        # SiO3 in 2010 (umol/kg)
    chem_ctl[3] = ds_dict_1d["po4"][i]         # PO4 in 2010 (umol/kg)
    chem_ctl[4] = ds_dict_1d["temp_2010"][i]   # Temperature in 2010 (degC)
    chem_ctl[5] = ds_dict_1d["sal_2010"][i]    # Salinity in 2010 (psu) 
    chem_ctl[6] = 0
    kwargs = {
    'chem_pi': chem_pi,
    'chem_ctl': chem_ctl,
    'oae_type': oae_type,
    'obj_var': obj_var,
    'cdreff': cdreff
    }
    nnn_pi = np.count_nonzero(~np.isnan(chem_pi))  # number of non-nan
    nnn_ctl = np.count_nonzero(~np.isnan(chem_ctl))  # number of non-nan
    if nnn_pi==7 and nnn_ctl==7:
        f_x = lambda x: f_dTA(x, **kwargs)
        root = opt.brentq(f_x, x_lwr, x_upr)
        return np.array([i,root])
# parallel processing loop through all grid cells (takes about 2.4 hours using 6 CPUs)
ncpu = 6   # number of CPU cores to use for parallel processing
with multiprocessing.Pool(processes=ncpu) as pool:    
    # Use imap_unordered to apply the function to a range of numbers
    results = pool.imap_unordered(find_root, range(ds_dict_1d["dTA_root"].shape[0]))    
    # Iterate over the results as they become available
    for result in results:
        if result is not None:
            i = int(result[0])
            root = result[1]
            ds_dict_1d["dTA_root"][i] = root
            # print(i, root) 
# reshape 1d dTA_root to 2d (180, 360)
ds_dict["dTA_root"] = ds_dict_1d["dTA_root"].reshape(ds_dict["talk_1750"].shape)
# Robinson map of the OAE needed to restore TA-DIC to pre-industrial
import cartopy.crs as ccrs
plt.figure(figsize=(8, 5))
X = np.linspace(0.5, 359.5, 360)
Y = np.linspace(-89.5, 89.5, 180)
Z = np.abs(ds_dict["dTA_root"]).copy()
zmin = np.nanpercentile(Z,0.1)
zmax = np.nanpercentile(Z,99.9)
Z[Z<zmin]=zmin
Z[Z>zmax]=zmax
ax = plt.axes(projection=ccrs.Robinson(central_longitude=180))
ax.set_global()
ax.coastlines()
# ax.gridlines()
plt.title('OAE needed to restore TA-DIC to pre-industrial\nusing NaOH with CDR efficiency = 80%')
plt.contourf(X,Y,Z,cmap='turbo',levels=256,transform=ccrs.PlateCarree());
plt.colorbar(orientation="horizontal", pad=0.03,label=r'$\Delta$TA needed to restore TA-DIC to PI $\mu$mol/kg',
            ticks=[0,25,50,75,100,125,150,175,200,225,250]);
plt.savefig('OAE_needed_for_OA_global_NaOH_CDReff80.png', format='png', dpi=300)
plt.show()
```
![OAE_needed_for_OA_global_NaOH_CDReff80](https://github.com/user-attachments/assets/52e554fa-e134-4e00-b99e-74eb79af9953)

Figure 1. The amount of OAE needed in the to restore the TA-DIC in 2010 to pre-industrial conditions, assuming that NaOH is used for OAE, and the CDR efficiency is 80%.

# Example use of the etamax function to analyze the maximum hypothetical OAE efficiency in the global oceans

In this example we will read global arrays of surface ocean chemistry data in 2010 from the Jiang_data_for_PyOAE.nc file, available in this repository, calculate the global array of ηmax for the year 2010, and plot a map of the results. The Jupyter Notebooks available at this repository provide examples of additional analysis of ηmax, and citations for the sources of the ocean chemistry data provided in the example netcdf file.<br>
```
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from PyOAE import etamax
# read the global arrays of surface ocean data and assign to a dictionary
ds = xr.open_dataset("Jiang_data_for_PyOAE.nc", chunks={"x":60})
ds_dict = {var: ds[var].values for var in ds.data_vars}
# extract the global arrays of chemistry data from the year 2010
kwargs = dict(
    TA_ctl = ds_dict["talk_2010"],    # TA (umol/kg)
    DIC_ctl = ds_dict["dic_2010"],    # DIC (umol/kg)
    SiO3_ctl = ds_dict["sio3"],       # SiO3 (umol/kg)
    PO4_ctl = ds_dict["po4"],         # PO4 (umol/kg)
    Temp_ctl = ds_dict["temp_2010"],  # temperatre (degC)
    Sal_ctl = ds_dict["sal_2010"],    # salinity (psu)
    Pres_ctl = np.zeros((180, 360))   # pressure (dbar)
)
dTA = 1    # assumed amount of added alkalinity by OAE (umol/kg)
# calculate the etamax for the dTA in all grid cells of the global array
result = etamax(dTA,**kwargs)
etamax = result["etamax"]
# plot a map of the results
plt.figure(figsize=(8, 5))  # Set the figure size (width, height)
plt.imshow(np.flipud(etamax), cmap='plasma', interpolation='none')
plt.colorbar(orientation="horizontal", pad=0.03)  # Add a colorbar for reference
plt.title(r'Figure 2. $\eta$max in 2010 with $\Delta$TA=1 $\mu$mol $kg^{-1}$')
plt.xticks([])
plt.yticks([])
plt.savefig('etamax_2010.png', format='png', dpi=300)
plt.show()
```
![etamax_2010](https://github.com/user-attachments/assets/ae40e53c-3c01-46c8-87a8-ffa209a2ad6f)

Figure 2. Values of the maximum theoretical OAE efficiency ηmax for the global oceans based on data from Jiang et al (2023) for a ∆TA perturbation of 1 umol/kg

# Example sensitivity of ηmax to the assumed amount of OAE addition of TA

In this example we will look at the difference in ηmax to the assumed amount of added alkalinity ∆TA (umol/kg). We will show the difference between ηmax assuming the ∆TA = 1 umol/kg, compared with ηmax assuming ∆TA = 100 umol/kg. 

```
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from PyOAE import etamax
# read the global arrays of surface ocean data and assign to a dictionary
ds = xr.open_dataset("jiang_data_for_jupyter_v12.nc", chunks={"lon":0})
ds_dict = {var: ds[var].values for var in ds.data_vars}
# extract the global arrays of chemistry data from the year 2010
kwargs = dict(
    TA_ctl = ds_dict["talk_2010"],    # TA (umol/kg)
    DIC_ctl = ds_dict["dic_2010"],    # DIC (umol/kg)
    SiO3_ctl = ds_dict["sio3"],       # SiO3 (umol/kg)
    PO4_ctl = ds_dict["po4"],         # PO4 (umol/kg)
    Temp_ctl = ds_dict["temp_2010"],  # temperatre (degC)
    Sal_ctl = ds_dict["sal_2010"],    # salinity (psu)
    Pres_ctl = np.zeros((180, 360))   # pressure (dbar)
)
# calculate the etamax for the dTA=1 and dTA=100 in all grid cells of the global array
result_dTA1 = etamax(1,**kwargs)
result_dTA100 = etamax(100,**kwargs)
etamax_dTA1 = result_dTA1["etamax"]
etamax_dTA100 = result_dTA100["etamax"]
# calculate difference between etamax at dTA=1 vs dTA=100
etamax_difference = etamax_dTA100 - etamax_dTA1
# plot a map of the results
plt.figure(figsize=(8, 5))  # Set the figure size (width, height)
plt.imshow(np.flipud(etamax_difference), cmap='viridis', interpolation='none')
plt.colorbar(orientation="horizontal", pad=0.03)  # Add a colorbar for reference
plt.title(r'Figure 3. Difference $\Delta\eta$max between $\Delta$TA 100 vs 1 umol/kg')
plt.xticks([])
plt.yticks([])
plt.savefig('etamax_sensitivity_to_dTA.png', format='png', dpi=300)
plt.show()
```
![etamax_sensitivity_to_dTA](https://github.com/user-attachments/assets/3e54a89d-d621-4be3-9cb0-80ecef264645)

Figure 3. Difference in ηmax comparing ∆TA = 100 umol/kg vs ∆TA = 1 umol/kg. 

# Example sensitivity of ηmax to different dissociation constants for carbonic acid with PyCO2SYS

The PyOAE functions use the PyCO2SYS package for the calculation of carbonate system variables. PyOAE uses default constants for the pH scale (total), carbonic acid of Lueker et al. (2000), bisulfate (HSO4−) of Dickson (1990), hydrofluoric acid (HF) of Perez and Fraga (1987), and the total borate content of Lee et al. (2010).

PyOAE allows the user to specify different constants using optional keyword arguments in "kwargs", using the same as keywords that are described in the PyCO2SYS documentation (https://github.com/mvdh7/PyCO2SYS) as follows:

- opt_pH_scale (Choice of pH scale, default=1)
- opt_k_carbonic (Choice of H2CO3 and HCO3- dissociation constants, default =10)
- opt_k_bisulfate (Choice of HSO4- dissociation constant KSO4, default=1)
- opt_total_borate (Choice of boron:sal, default=2)
- opt_k_fluoride (Choice of hydrogen fluoride dissociation constant, default=2)

Below is an example showing the sensitivity of ηmax, with ∆TA=1 umol/kg, for the following different dissociation constants for carbonic acid:

- Option 1: opt_k_carbonic=10 (default using Lueker et al 2000)
- Option 2: opt_k_carbonic=4 (refit of Mehrbach 1973 by Dickson and Millero 1987)

```
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from PyOAE import etamax
# read the global arrays of surface ocean data and assign to a dictionary
ds = xr.open_dataset("jiang_data_for_jupyter_v12.nc", chunks={"lon":0})
ds_dict = {var: ds[var].values for var in ds.data_vars}
# extract the global arrays of chemistry data from the year 2010
kwargs_opt1 = dict(
    TA_ctl = ds_dict["talk_2010"],    # TA (umol/kg)
    DIC_ctl = ds_dict["dic_2010"],    # DIC (umol/kg)
    SiO3_ctl = ds_dict["sio3"],       # SiO3 (umol/kg)
    PO4_ctl = ds_dict["po4"],         # PO4 (umol/kg)
    Temp_ctl = ds_dict["temp_2010"],  # temperatre (degC)
    Sal_ctl = ds_dict["sal_2010"],    # salinity (psu)
    Pres_ctl = np.zeros((180, 360))   # pressure (dbar)
)
kwargs_opt2 = dict(
    TA_ctl = ds_dict["talk_2010"],    # TA (umol/kg)
    DIC_ctl = ds_dict["dic_2010"],    # DIC (umol/kg)
    SiO3_ctl = ds_dict["sio3"],       # SiO3 (umol/kg)
    PO4_ctl = ds_dict["po4"],         # PO4 (umol/kg)
    Temp_ctl = ds_dict["temp_2010"],  # temperatre (degC)
    Sal_ctl = ds_dict["sal_2010"],    # salinity (psu)
    Pres_ctl = np.zeros((180, 360)),   # pressure (dbar)
    opt_k_carbonic = 4                # refit of Mehrbach 1973 by Dicskon and Millero 1987
)
# calculate the etamax for the dTA=1 and dTA=100 in all grid cells of the global array
result_opt1 = etamax(1,**kwargs_opt1)
result_opt2 = etamax(1,**kwargs_opt2)
etamax_opt1 = result_opt1["etamax"]
etamax_opt2 = result_opt2["etamax"]
# calculate difference between etamax between Option 1 and 2
etamax_difference = etamax_opt2 - etamax_opt1
# plot a map of the results
plt.figure(figsize=(8, 5))  # Set the figure size (width, height)
plt.imshow(np.flipud(etamax_difference), cmap='viridis', interpolation='none')
plt.colorbar(orientation="horizontal", pad=0.03)  # Add a colorbar for reference
plt.title(r'Figure 4. Difference $\Delta\eta$max comparing opt_k_carbonic = 4 vs 10')
plt.xticks([])
plt.yticks([])
plt.savefig('etamax_sensitivity_to_K1K2.png', format='png', dpi=300)
plt.show()
```
![etamax_sensitivity_to_K1K2](https://github.com/user-attachments/assets/d69deae9-bc96-415e-b7da-4eba695a6c72)

Figure 4. Difference in ηmax, with a ∆TA perturbation of 1 umol/kg, comparing two options for the dissociation constants for carbonic acid (option 4= refit of Mehrbach 1973 by Dickson and Millero 1987, option 10= Lueker et 2000)

# Example calculation of the biological component of DIC in the global oceans

The surface water DIC concentration of a given water mass is altered by air–sea gas exchange, the biological production/consumption of organic matter, and calcium carbonate (CaCO3) formation/dissolution (Burt et al, 2016). To isolate the biological component of DIC (DIC_bio), a surface DIC concentration at atmospheric equilibrium (DIC_atm) is computed, and subsequently removed from the observed DIC (DIC_obs). 

We use the following equation to represent DIC_bio:

DIC_bio = DIC_obs – DIC_atm				(eqn 1)


where 

- DIC_obs = observed DIC (OceanSODA-ETHZ)
- DIC_atm = DIC at equilibrium with atmospheric pCO2 and observed TA

We use the DIC reported in OceanSODA-ETHZ as the value of "DIC_obs" in each month at each grid
cell from 1982-2022. Next, we calculate "DIC_atm" with PyCO2SYS using
the atmospheric fCO2 from the SeaFlux data set, combined the observed
TA from OceanSODA-ETHZ in each grid cell for each month from
1982-2022. Therefore, "DIC_atm" represents the hypothetical DIC that
would be in equilibrium with the atmospheric pCO2. Finally, we
calculated the 1982-2022 monthly values of "DIC_bio" using eqn 1 as
DIC_bio = DIC_obs - DIC_atm

The repeating annual cycle of DIC_bio can be represented as a sine function of the following form:

y = mean + amplitude * sin(2π * (x - phase) / period)	(eqn 2)

where 

- y = DIC_bio = DIC_obs – DIC_atm
- x = time as decimal year fraction (1982-2022) 
- mean = mean from sine-regression
- amplitude = amplitude from sine-regression
- phase = phase shift from sine-regression
- period = assumed equal to 1 cycle per year

PyOAE incudes functions (dic_bio and sine_fit) to calculate DIC_bio using eqn 1, and to find the optimum parameters of the sine-regression model in eqn 2. We solve for the best-fit values of the parameters in eqn 2 at each grid cell including the mean,
amplitude, and phase of the sine function at each location. 

In this example we use two netcdf files that we need to do the analysis, OceanSODA_ETHZ_for_PyOAE.nc and SeaFlux_for_PyOAE.nc, available to download at the following link:

https://drive.google.com/drive/folders/1BGgVRk2Gf6mxNnX1Fxg0Q4GtZSAYMzef?usp=sharing

We use the 1982-2022 monthly time series of DIC_bio as the observed "y" values in
eqn 2 to estimate the best-fit parameters of eqn 2 at each location.

References:
- Clargo et al 2015 (https://doi.org/10.1016/j.marchem.2015.08.010)
- Burt et al 2016 (https://doi.org/10.1002/lno.10243):

In the following code block we loop through all of the grid cells in the global oceans to do the analysis. This takes about 45 minutes to calculate the monthly DIC_bio in each grid cell from 1982-2022, and solve for the best fit sine-regression model in each grid cell 

```
print('Computing DIC_bio, this takes about 45 minutes, please wait ...')
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from PyOAE import dic_bio, sine_fit
# read the data into an xarray dataset
ds1 = xr.open_dataset("OceanSODA_ETHZ_for_PyOAE.nc", chunks={"lon":0})
ds2 = xr.open_dataset("SeaFlux_for_PyOAE.nc", chunks={"lon":0})
# Convert ds1 to dictionary of numpy arrays for computations
ds_dict = {var: ds1[var].values for var in ds1.data_vars}
# append yearfrac,lon,lat,time,pco2atm,fco2atm to ds_dict
ds_dict["yearfrac"] = ds1.yearfrac.values
ds_dict["lon"] = ds1.lon.values
ds_dict["lat"] = ds1.lat.values
ds_dict["time"] = ds1.time.values
ds_dict["pco2atm"] = ds2.pco2atm.values
ds_dict["fco2atm"] = ds2.fco2atm.values
# initialize new output arrays for ds_dict
ds_dict["dic_atm"] = np.full_like(ds_dict["talk"], np.nan)
ds_dict["dic_bio"] = np.full_like(ds_dict["talk"], np.nan)
ds_dict["dic_bio_fit"] = np.full_like(ds_dict["talk"], np.nan)
ds_dict["dic_bio_mean"] = np.full_like(ds_dict["lon"], np.nan)
ds_dict["dic_bio_amplitude"] = np.full_like(ds_dict["lon"], np.nan)
ds_dict["dic_bio_phase"] = np.full_like(ds_dict["lon"], np.nan)
ds_dict["dic_bio_rmse"] = np.full_like(ds_dict["lon"], np.nan)
ds_dict["dic_bio_adj_rsquared"] = np.full_like(ds_dict["lon"], np.nan)
ds_dict["dic_bio_pvalue"] = np.full_like(ds_dict["lon"], np.nan)
for i in range(ds_dict["talk"].shape[2]):
    print("dic_bio computed at lon %.1f degE" % (i+0.5))
    for j in range(ds_dict["talk"].shape[1]):
        kwargs = {
            'alkalinity': ds_dict["talk"][:,j,i],
            'dic': ds_dict["dic"][:,j,i],
            'pco2atm': ds_dict["fco2atm"][:,j,i],
            'pco2atm_type': 5,  # The second parameter 4=pCO2, 5=fCO2
            'total_silicate': ds_dict["sio3"][:,j,i],
            'total_phosphate': ds_dict["po4"][:,j,i],
            'temperature': ds_dict["temperature"][:,j,i],
            'salinity': ds_dict["salinity"][:,j,i],
            'total_pressure': 0,
            'opt_pH_scale': 1,  # pH scale (1= total scale)
            'opt_k_carbonic': 10,  # Choice of H2CO3 and HCO3- K1 and K2 (10= Lueker et al 2000)
            'opt_k_bisulfate': 1,  # Choice of HSO4- dissociation constant KSO4 (1= Dickson)
            'opt_total_borate': 1,  # Choice for boron:sal 
            'opt_k_fluoride': 1   # Choice for fluoride
            }
        # nnn = number of non-nan values for each input variable to dic_bio         
        nnn_talk = np.count_nonzero(~np.isnan(kwargs["alkalinity"]))  # number of non-nan
        nnn_dic = np.count_nonzero(~np.isnan(kwargs["dic"]))  # number of non-nan
        nnn_pco2atm = np.count_nonzero(~np.isnan(kwargs["pco2atm"]))  # number of non-nan
        nnn_sio3 = np.count_nonzero(~np.isnan(kwargs["total_silicate"]))  # number of non-nan
        nnn_po4 = np.count_nonzero(~np.isnan(kwargs["total_phosphate"]))  # number of non-nan
        nnn_temp = np.count_nonzero(~np.isnan(kwargs["temperature"]))  # number of non-nan
        nnn_sal = np.count_nonzero(~np.isnan(kwargs["salinity"]))  # number of non-nan
        if (nnn_talk==nnn_dic and nnn_talk==nnn_pco2atm and nnn_talk==nnn_sio3 and 
            nnn_talk==nnn_po4 and nnn_talk==nnn_temp and nnn_talk==nnn_sal and nnn_talk > 0):
            # solve for dic_bio = dic_obs - dic_atm
            result_1 = dic_bio(**kwargs)
            ds_dict["dic_atm"][:,j,i]= result_1["dic_atm"]
            ds_dict["dic_bio"][:,j,i] = result_1["dic_bio"]
            # solve for dic_bio vs time sine-regression mean, amplitude, phase, rmse, y_fit
            result_2 = sine_fit(ds_dict["yearfrac"],result_1["dic_bio"])  
            ds_dict["dic_bio_fit"][:,j,i]= result_2["y_fit"]
            ds_dict["dic_bio_mean"][j,i]= result_2["mean"]
            ds_dict["dic_bio_amplitude"][j,i]= result_2["amplitude"]
            ds_dict["dic_bio_phase"][j,i]= result_2["phase"]
            ds_dict["dic_bio_rmse"][j,i]= result_2["rmse"]
            ds_dict["dic_bio_adj_rsquared"][j,i]= result_2["adj_rsquared"]
            ds_dict["dic_bio_pvalue"][j,i]= result_2["pvalue"]
print('Done')
```

Next we will make a map showing the results for the regression estimate of the mean DIC_bio. Positive values, shown in red, indicate that long-term mean DIC_obs > DIC_atm. Negative values, shown in blue, indicate that long-term mean DIC_obs < DIC_atm.

```
# Robinson map of the sine-regression mean values of DIC_bio
import cartopy.crs as ccrs
from matplotlib.colors import TwoSlopeNorm
plt.figure(figsize=(8, 5),dpi=150)
X = ds_dict['lon']
Y = ds_dict['lat']
Z = ds_dict['dic_bio_mean']
# Define the zero point
vmin = np.nanpercentile(Z,1)
vmax = np.nanpercentile(Z,99)
vcenter = 0.0
# Create a normalization instance
norm = TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
ax = plt.axes(projection=ccrs.Robinson(central_longitude=180))
ax.set_global()
ax.coastlines()
# ax.gridlines()
plt.title(r'Mean of 1982-2022 DIC_bio')
plt.contourf(X,Y,Z,cmap='seismic',levels=256,transform=ccrs.PlateCarree(), norm=norm);
plt.colorbar(orientation="horizontal", pad=0.03,label='Mean DIC_bio (umol/kg)',ticks=[-60,-50,-40,-30,-20,-10,0,10,20,30,40,50]);
plt.savefig('Fig5_map_of_DIC_bio_mean_using_fco2atm_as_fco2_v2.png', format='png', dpi=300)
plt.show()
```
![Fig5_map_of_DIC_bio_mean_using_fco2atm_as_fco2_v2](https://github.com/user-attachments/assets/b25339b4-be17-4744-967b-b9bab12d79e6)

Figure 5. Mean DIC_bio during 1982-2022 from sine-regression

Next we will make a map showing the results of the regression estimate of the amplitude of DIC_bio. Note that the amplitude represents half of the distance from the peak to the trough of each annual cycle. Therefore the amplitude is half of the annual range of DIC_bio

```
# Robinson map of the sine-regression amplitude values of DIC_bio
import cartopy.crs as ccrs
plt.figure(figsize=(8, 5),dpi=150)
X = ds_dict['lon']
Y = ds_dict['lat']
Z = np.abs(ds_dict["dic_bio_amplitude"])
ax = plt.axes(projection=ccrs.Robinson(central_longitude=180))
ax.set_global()
ax.coastlines()
# ax.gridlines()
plt.title(r'Amplitude of 1982-2022 DIC_bio')
plt.contourf(X,Y,Z,cmap='turbo',levels=256,transform=ccrs.PlateCarree());
plt.colorbar(orientation="horizontal", pad=0.03,label='Amplitude DIC_bio (umol/kg)',ticks=[0,5,10,15,20,25,30,35,40,45,50]);
plt.savefig('Fig6_map_of_DIC_bio_amplitude_using_fco2atm_as_fco2_v2.png', format='png', dpi=300)
plt.show()
```
![Fig6_map_of_DIC_bio_amplitude_using_fco2atm_as_fco2_v2](https://github.com/user-attachments/assets/b3e271ce-bb6c-43a6-89cb-664edd5cb7ae)

Figure 6. Amplitude of DIC_bio during 1982-2022 from sine-regression

Next we will make a map showing the p-values of the sine-regressions. Most of the grid cells have statistically signficant regressions (p<0.05) shown in blue

```
# Robinson map of the sine-regression rsquared values of DIC_bio
from matplotlib.colors import TwoSlopeNorm
from matplotlib.ticker import LinearLocator
plt.figure(figsize=(8, 5),dpi=150)
X = ds_dict['lon']
Y = ds_dict['lat']
Z = np.abs(ds_dict['dic_bio_pvalue'])
Z[Z>0.1]=0.1
# Define the zero point
vmin = 0
vmax = 0.1
vcenter = 0.05
# Create a normalization instance
norm = TwoSlopeNorm(vmin=vmin,vcenter=vcenter,vmax=vmax)
ax = plt.axes(projection=ccrs.Robinson(central_longitude=180))
ax.set_global()
ax.coastlines()
plt.title(r'p-value of 1982-2022 DIC_bio')
plt.contourf(X,Y,Z,cmap='seismic',transform=ccrs.PlateCarree(),norm=norm,levels=256);
cbar = plt.colorbar(orientation="horizontal", pad=0.03,label='p-value',
                    ticks=[0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]);
plt.savefig('Fig7_map_of_DIC_bio_pvalue_using_fco2atm_as_fco2_v2.png', format='png', dpi=300)
plt.show()
```
![Fig7_map_of_DIC_bio_pvalue_using_fco2atm_as_fco2_v2](https://github.com/user-attachments/assets/29674081-f0cb-4f42-8dea-c82e668a0d97)

Figure 7. p-values of DIC_bio during 1982-2022 from sine-regression

Next, we will make a map showing the adjusted r^2 values for the regressions in each grid cell. Note that the r^2 values tend to be greater in areas that have the largest amplitudes of DIC_bio. In other words, the regressions are best in areas where there is the greatest biogeochemical effect in DIC_bio

```
# Robinson map of the sine-regression rsquared values of DIC_bio
from matplotlib.colors import TwoSlopeNorm
plt.figure(figsize=(8, 5),dpi=150)
X = ds_dict['lon']
Y = ds_dict['lat']
Z = np.abs(ds_dict['dic_bio_adj_rsquared'])
# Define the zero point
vmin = 0
vmax = 1.0
vcenter = 0.5
# Create a normalization instance
norm = TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
ax = plt.axes(projection=ccrs.Robinson(central_longitude=180))
ax.set_global()
ax.coastlines()
# ax.gridlines()
plt.title(r'r-squared of 1982-2022 DIC_bio')
plt.contourf(X,Y,Z,cmap='seismic_r',transform=ccrs.PlateCarree(), norm=norm,levels=256);
# plt.contourf(X,Y,Z,cmap='plasma',levels=256,transform=ccrs.PlateCarree());
plt.colorbar(orientation="horizontal", pad=0.03,label='Adjusted r-squared',ticks=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]);
plt.savefig('Fig8_map_of_DIC_bio_rsquared_using_fco2atm_as_fco2_v2.png', format='png', dpi=300)
plt.show()
```
![Fig8_map_of_DIC_bio_rsquared_using_fco2atm_as_fco2_v2](https://github.com/user-attachments/assets/de91ce0a-7e79-4750-9d13-22bf58c8d121)

Figure 8. Adjusted r-squared of DIC_bio during 1982-2022 from sine-regression

Finally, we will show the time series of DICobs, DICatm, and DICbio at a selected location in the coastal California Current Ecosystem near the Columbia River

```
# Columbia River location i,j coordinates
i=234
j=136
fig, ax = plt.subplot_mosaic(
    '''
    A
    B
    ''',
    figsize = (9, 9), dpi=150
    # constrained_layout = True
)

# Plot at ['A']
ax['A'].plot(ds_dict["yearfrac"], ds_dict["dic"][:,j,i], label='DIC_obs (OceanSODA-ETHZ)', linestyle='-', marker='')
ax['A'].plot(ds_dict["yearfrac"], ds_dict["dic_atm"][:,j,i], label='DIC_atm (equilibrium with atmospheric pCO2)', linestyle='-')
# ax['A'].set_xlabel('year')
ax['A'].set_ylabel('DIC (umol/kg)')
ax['A'].legend(loc='upper left')
# ax['A'].grid(True)
ax['A'].set_title('a. DIC observed (DIC_obs), and at equilibrium with atmospheric pCO2 (DIC_atm)')
# ax['A'].set_xlim(2010, 2021)
ax['A'].set_xlim(1982, 2021)
ax['A'].set_ylim(1850, 2085)
# ax['A'].text(1983, 36.5, 'Mean: '+f"{A_fit:.1f}"+', Amplitude: '+f"{B_fit:.1f}"+', RMSE: '+f"{rmse:.1f}"+' umol/kg',
#         fontsize=10, color='black', ha='left', va='center')
# ax['A'].axhline(y=0, color='k', linestyle=':')

# Plot at ['B']
ax['B'].plot(ds_dict["yearfrac"], ds_dict["dic_bio"][:,j,i], label='DIC_bio = DIC_obs - DIC_atm', color='black', linestyle='-', marker='')
ax['B'].plot(ds_dict["yearfrac"], ds_dict["dic_bio_fit"][:,j,i], label='Regression', color='red', linestyle='--')
# ax['A'].set_xlabel('year')
ax['B'].set_ylabel('DIC_bio (umol/kg)')
ax['B'].legend(loc='upper left')
# ax['B'].grid(True)
ax['B'].set_title('b. DIC_bio = DIC_obs - DIC_atm')
# ax['B'].set_xlim(2010, 2021)
ax['B'].set_xlim(1982, 2021)
ax['B'].set_ylim(-35, 15)
# ax['B'].text(1983, -10, 'Mean: '+f"{A_fit:.1f}"+', Amplitude: '+f"{B_fit:.1f}"+', RMSE: '+f"{rmse:.1f}"+' umol/kg',
#         fontsize=10, color='black', ha='left', va='center')
ax['B'].text(1983, -32.5, 'Mean: '+f"{ds_dict["dic_bio_mean"][j,i]:.1f}"+', Amplitude: '+f"{ds_dict["dic_bio_amplitude"][j,i]:.1f}"+', RMSE: '+f"{ds_dict["dic_bio_rmse"][j,i]:.1f}"+' umol/kg',
        fontsize=10, color='black', ha='left', va='center')
ax['B'].axhline(y=0, color='k', linestyle=':')

fig.savefig('Fig9_DIC_obs_atm_bio_at_ColumbiaRiver_1982_2020_fco2atm_as_fco2.png', format='png');
```
![Fig9_DIC_obs_atm_bio_at_ColumbiaRiver_1982_2020_fco2atm_as_fco2](https://github.com/user-attachments/assets/f8aeaecf-998c-42b7-9244-ca427132b0d5)

Figure 9. DIC_obs, DIC_atm, and DIC_bio at a coastal location in the California Current Ecosystem near the Columbia River from 1982-2020. **a)** Observed DIC (DIC_obs), and the hypothetical DIC (DIC_atm) that would be in equilibrium with atmospheric pCO2 and observed TA. **b)** DIC_bio calculated from DIC_obs-DIC_atm, and the regression estimate of DIC_bio from the sine regression.












