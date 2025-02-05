# PyOAE
Python tools for analysis of Ocean Alkalinity Enhancement

by Greg Pelletier

Tools for analysis of Ocean Alkalinity Enhancement (OAE), including Jupyter Notebooks for the following examples:

- Root-finding method to solve for the OAE treatment needed to restore any carbonate system variable (e.g. TA-DIC, CO3--, pH, Ωara, Ωcal) to pre-industrial conditions in the coastal California Current Ecosystem. This notebook can also be used to solve for any other location in the global oceans
- Calculation of the maximum hypothetical OAE efficiency ηmax (etamax) for any assumed addition of alkalinity. The ηmax is a dimensionless quantity that is the hypothetical maximum potential CDR (umol/kg) divided by the amount of added alkalinity (umol/kg).

The PyOAE package requires that you have already installed numpy, scipy, and PyCO2SYS packages. We also recommend that you have installed xarray, cartopy, and matplotlib to analyze and plot maps using data from netcdf files.

# Installation for Python or Jupyter Notebook

First, if you have not already done so, install PyOAE as follows with pip or !pip in your notebook or terminal (append --upgrade to the following line if you are upgrading from a previous installation):<br>
```
pip install git+https://github.com/gjpelletier/PyOAE.git
```

Next import the f_dTA and etamax functions as follows in your notebook or python code:<br>
```
from PyOAE import f_dTA, etamax
```

As an alternative to the commands above, you can download PyOAE.py from this github repository (https://github.com/gjpelletier/PyOAE) and copy and paste the functions into your own project.<br>

# Example use of the root-finding method to analyze the OAE needed to restore OA indicators to pre-industrial conditions

The difference between TA and DIC, also known as Alk* (Sarmiento & Gruber, 2006), can be used as a surrogate variable to interpret the response of other carbonate system variables (e.g. CO3--, pH, Ωara, Ωcal). Ocean acidification (OA) has caused a decrease in TA-DIC since pre-industrial conditions (Sarmiento & Gruber, 2006).

In this example we will use the root-finding method to solve for the amount of OAE needed to restore the TA-DIC in 2010 in the coastal California Current Ecosystem (CCE) to pre-industrial conditions. A detailed explanation of the source of the ocean chemistry data used in this example is provided in the Jupyter Notebooks available in this repository

The current version of f_dTA analyzes the ocean chemistry data from one grid cell at a time. Processing each grid cell takes about 1-3 seconds. To analyze all of the grid cells in a model domain, or a subset for a region of selected grid cells, we use Python to loop through all of the grid cells that need to be evaluated, and solve for the root in each grid cell one at a time in the loop. 

First we will solve for the amount of OAE needed to restore the TA-DIC in the coastal CCE to pre-industrial as follows:<br>
```
# import the packages that are needed
import scipy.optimize as opt
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from PyOAE import f_dTA
# read the global arrays of surface ocean data and assign to a dictionary
ds = xr.open_dataset("jiang_data_for_jupyter_v12.nc", chunks={"lon":0})
ds_dict = {var: ds[var].values for var in ds.data_vars}
# initialize output array and specify options
ds_dict["dTA_root"] = np.full((180, 360), np.nan) # init out array 
obj_var = 'alkstar'   # objective variable: 'alkstar', 'co3', phtot', 'omara', 'omcal'
oae_type = 'NaOH'     # chemical used for OAE: 'NaOH' or 'Na2CO3'
cdreff = 0.8          # CDR efficiency between 0-1 (e.g. 0.8 = 80%)
x_upr = 0             # lower bound of possible dTA values (umol/kg)
x_lwr = 500           # upper bound of possible dTA values (umol/kg)
# main loop through all grid cells
for i, j in np.ndindex((180,360)):
    if ds_dict["dist2coast"][i,j]<=100 and ds_dict["LME"][i,j]==11:  # coastal CCE
        chem_pi = np.full(7, np.nan)
        chem_pi[0] = ds_dict["talk_1750"][i,j]    # TA in 1750 (umol/kg)
        chem_pi[1] = ds_dict["dic_1750"][i,j]     # DIC in 1750 (umol/kg)
        chem_pi[2] = ds_dict["sio3"][i,j]         # SiO3 in 1750 (umol/kg)
        chem_pi[3] = ds_dict["po4"][i,j]          # PO4 in 1750 (umol/kg)
        chem_pi[4] = ds_dict["temp_1750"][i,j]    # Temperature in 1750 (degC)
        chem_pi[5] = ds_dict["sal_1750"][i,j]     # Salinity in 1750 (psu)        
        chem_pi[6] = 0
        chem_ctl = np.full(7, np.nan)
        chem_ctl[0] = ds_dict["talk_2010"][i,j]   # TA in 2010 (umol/kg)
        chem_ctl[1] = ds_dict["dic_2010"][i,j]    # DIC in 2010 (umol/kg)
        chem_ctl[2] = ds_dict["sio3"][i,j]        # SiO3 in 2010 (umol/kg)
        chem_ctl[3] = ds_dict["po4"][i,j]         # PO4 in 2010 (umol/kg)
        chem_ctl[4] = ds_dict["temp_2010"][i,j]   # Temperature in 2010 (degC)
        chem_ctl[5] = ds_dict["sal_2010"][i,j]    # Salinity in 2010 (psu) 
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
            root = opt.brentq(f_x, x_upr, x_lwr)
            ds_dict["dTA_root"][i,j] = root
            print("i: %.0f, j: %.0f, root: %.4f" % (i,j,root))
```

Next we will make a map showing the results:
```
fig = plt.figure(figsize=(5.1, 5.5))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
# define lon and lat
lon = np.linspace(0.5, 359.5, 360)
lat = np.linspace(-89.5, 89.5, 180)
lon2d, lat2d = np.meshgrid(lon, lat)
# define ticks and labels
lon_ticks = range(-120, -100, 10)  # Longitude ticks every 30 degrees
lat_ticks = range(20, 55, 5)    # Latitude ticks every 30 degrees
xtick_labels = ['120°W', '110°W']
ytick_labels = ['20°N', '25°N', '30°N', '35°N', '40°N', '45°N', '50°N']
# summary stats
# plotdata = ds_dict["dTA_root"]
mean = np.nanmean(ds_dict["dTA_root"])
stdev = np.nanstd(ds_dict["dTA_root"])
# ax.set_title('Figure 1. OAE needed to\nrestore OA in coastal CCE\nto pre-industrial')
ax.set_title('Figure 1. OAE needed to restore\nTA-DIC to pre-industrial')
ax.projection = ccrs.PlateCarree()
ax.set_extent([-127.5+360, -107.5+360, 20, 50], crs=ccrs.PlateCarree())
ax.set_xticks(lon_ticks, crs=ccrs.PlateCarree())
ax.set_yticks(lat_ticks, crs=ccrs.PlateCarree())
ax.set_xticklabels(xtick_labels)
ax.set_yticklabels(ytick_labels)
ax.add_feature(cfeature.LAND, facecolor='lightgray')
ax.coastlines('10m', edgecolor='none', linewidth=0.5)
ax.text(-126.5, 24, 'Mean: '+f"{mean:.0f} umol/kg", transform=ccrs.PlateCarree(),
        fontsize=11, color='black', ha='left', va='center')
ax.text(-126.5, 22, 'Stdev: '+f"{stdev:.0f} umol/kg", transform=ccrs.PlateCarree(),
        fontsize=11, color='black', ha='left', va='center')
cmap = plt.get_cmap('plasma').reversed()
contour = ax.pcolor(lon, lat, ds_dict["dTA_root"], transform=ccrs.PlateCarree(), cmap=cmap)
cbar = plt.colorbar(contour, orientation='vertical', pad=0.05)
cbar.set_label(r'$\Delta$TA needed to restore TA-DIC to PI $\mu$mol/kg')
plt.savefig('OAE_needed_for_OA_in_CCE.png', format='png')
```
![OAE_needed_for_OA_in_CCE](https://github.com/user-attachments/assets/ebe2b86e-8184-4887-9378-ecd657b69f83)

Figure 1. The amount of OAE needed in the coastal California Current Ecosystem to restore the TA-DIC in 2010 to pre-industrial conditions, assuming that NaOH is used for OAE, and the CDR efficiency is 80%.

# Example use of the etamax function to analyze the maximum hypothetical OAE efficiency in the global oceans

In this example we will read global arrays of surface ocean chemistry data in 2010 from the jiang_data_for_jupyter_v12.nc file, available in this repository, calculate the global array of ηmax for the year 2010, and plot a map of the results. The Jupyter Notebooks available at this repository provide examples of additional analysis of ηmax, and citations for the sources of the ocean chemistry data provided in the example netcdf file.<br>
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

# Example sensitivity of etamax to the assumed amount of OAE addition of TA

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

Figure 3. Difference in ηmax comparing ∆TA perturbation of 100 umol/kg vs ∆TA perturbation of 1 umol/kg. 

# Example using different constants with PyCO2SYS

The PyOAE functions use the PyCO2SYS package for the calculation of carbonate system variables. PyOAE uses default constants for the pH scale (total), carbonic acid of Lueker et al. (2000), bisulfate (HSO4−) of Dickson (1990), hydrofluoric acid (HF) of Perez and Fraga (1987), and the total borate content of Lee et al. (2010).

PyOAE allows the user to specify different constants using optional keyword arguments in "kwargs", using the same as keywords that are described in the PyCO2SYS documentation (https://github.com/mvdh7/PyCO2SYS) as follows:

- opt_pH_scale (Choice of pH scale, default=1)
- opt_k_carbonic (Choice of H2CO3 and HCO3- dissociation constants, default =10)
- opt_k_bisulfate (Choice of HSO4- dissociation constant KSO4, default=1)
- opt_total_borate (Choice of boron:sal, default=2)
- opt_k_fluoride (Choice of hydrogen fluoride dissociation constant, default=2)

Below is an example showing the sensitivity of etamax, with a perturbation ∆TA=1 umol/kg, for the following different dissociation constants for carbonic acid:

- Scenario 1: opt_k_carbonic=10 (default using Lueker et al 2000)
- Scenario 2: opt_k_carbonic=4 (refit of Mehrbach 1973 by Dickson and Millero 1987)

```
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from PyOAE import etamax
# read the global arrays of surface ocean data and assign to a dictionary
ds = xr.open_dataset("jiang_data_for_jupyter_v12.nc", chunks={"lon":0})
ds_dict = {var: ds[var].values for var in ds.data_vars}
# extract the global arrays of chemistry data from the year 2010
kwargs_scen1 = dict(
    TA_ctl = ds_dict["talk_2010"],    # TA (umol/kg)
    DIC_ctl = ds_dict["dic_2010"],    # DIC (umol/kg)
    SiO3_ctl = ds_dict["sio3"],       # SiO3 (umol/kg)
    PO4_ctl = ds_dict["po4"],         # PO4 (umol/kg)
    Temp_ctl = ds_dict["temp_2010"],  # temperatre (degC)
    Sal_ctl = ds_dict["sal_2010"],    # salinity (psu)
    Pres_ctl = np.zeros((180, 360))   # pressure (dbar)
)
kwargs_scen2 = dict(
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
result_scen1 = etamax(1,**kwargs_scen1)
result_scen2 = etamax(1,**kwargs_scen2)
etamax_scen1 = result_scen1["etamax"]
etamax_scen2 = result_scen2["etamax"]
# calculate difference between etamax between Scenario 1 and 2
etamax_difference = etamax_scen2 - etamax_scen1
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





