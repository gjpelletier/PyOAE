# PyOAE
Python tools for analysis of Ocean Alkalinity Enhancement

by Greg Pelletier

Tools for analysis of Ocean Alkalinity Enhancement, including Jupyter Notebooks for the following examples:

- Root-finding method to solve for the OAE treatment needed to restore any carbonate system variable (e.g. TA-DIC, CO3--, pH, Ωar, etc.) to pre-industrial conditions in the coastal California Current Ecosystem. This notebook can also be used to solve for any other location in the global oceans
- Calculation of the maximum hypothetical OAE efficiency etamax in the global oceans

The PyOAE package requires that you have already installed numpy, scipy, and PyCO2SYS packages. We also recommend that you have installed xarray, cartopy, and matplotlib to analyze and plot maps using from input data from netcdf files.

# Installation for Python, Jupyter Notebook, and Google Colab

First install the new functions as follows with pip or !pip in your notebook or terminal:<br>
```
pip install git+https://github.com/gjpelletier/PyOAE.git
```

Next import the delta_method and parametric_bootstrap functions as follows in your notebook or python code:<br>
```
from PyOAE import f_dTA, etamax
```

As an alternative, you can also download PyOAE.py from this github repository (https://github.com/gjpelletier/PyOAE) and add both functions to your own project.<br>

# Example use of root-finding method for Jupyter Notebook or Google Colab

The first step is to install PyOAE from github as follows:<br>
```
!pip install git+https://github.com/gjpelletier/PyOAE.git
```

Next we need to import the f_dTA function, and also scipy.optimize and numpy as follows:<br>
```
from PyOAE import f_dTA
import scipy.optimize as opt
import numpy as np
```

Next we will analyze a simple example of data from a single model grid cell
```
# import the packages that are needed
import scipy.optimize as opt
import numpy as np
from PyOAE import f_dTA
# assign the inputs that are needed
x_upr = 0   # lower bound of the range of dTA values (umol/kg) to search for the root
x_lwr = 500 # upper bound of the range of dTA values (umol/kg) to search for the root
chem_pi = np.array([2232,1861,1.346,0.201,26.683,34.004,0])    # TA, DIC, SiO3, PO4, temp, sal, pres for PI
chem_ctl = np.array([2230,1915,1.346,0.201,27.391,33.914,0])   # TA, DIC, SiO3, PO4, temp, sal, pres for control
oae_type = 'NaOH'     # 'NaOH' or 'Na2CO3'
obj_var = 'alkstar'   # 'alkstar', 'co3', 'omara', 'omcal', or 'phtot'
cdreff = 0.8          # e.g. use 0.8 for 80% CDR efficiency
# make the kwargs for f_dTA:
kwargs = {
  'chem_pi': chem_pi,
  'chem_ctl': chem_ctl,
  'oae_type': oae_type,
  'obj_var': obj_var,
  'cdreff': cdreff
  }
# make the lambda function that will be used to allow brentq to use the kwargs for f_dTA
f_x = lambda x: f_dTA(x, **kwargs)
# use brentq to find the root of dTA that results in OAE treated condtions equal to pre-industrial
root = opt.brentq(f_x, x_upr, x_lwr)
print("The dTA needed for restoration to pre-industrial conditions is %.2f umol/kg" % (root))
```


