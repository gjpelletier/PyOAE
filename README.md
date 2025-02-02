# PyOAE
Python tools for analysis of Ocean Alkalinity Enhancement

by Greg Pelletier

Tools and examples for using Python and Jupyter for analysis of Ocean Alkalinity Enhancement, including the following:

- root-finding method to solve for the OAE treatment needed to restore any carbonate system variable (e.g. TA-DIC, CO3--, pH, Ωar, etc.) to pre-industrial conditions
- calculation of the maximum hypothetical OAE efficiency etamax
- calculation of the maximum hypothetical potential CDR

EXAMPLE 1 - use of f_dTA function for root-finding to process one value at a time

import scipy.optimize as opt
import numpy as np
from PyOAE import f_dTA
x_upr = 0   # lower bound of the range of dTA values (umol/kg) to search for the root
x_lwr = 500 # upper bound of the range of dTA values (umol/kg) to search for the root
chem_pi = np.array([2232,1861,1.346,0.201,26.683,34.004,0])    # TA, DIC, SiO3, PO4, temp, sal, pres for PI
chem_ctl = np.array([2230,1915,1.346,0.201,27.391,33.914,0])   # TA, DIC, SiO3, PO4, temp, sal, pres for control
oae_type = 'NaOH'     # 'NaOH' or 'Na2CO3'
obj_var = 'alkstar'   # 'alkstar', 'co3', 'omara', 'omcal', or 'phtot'
cdreff = 0.8          # e.g. use 0.8 for 80% CDR efficiency
kwargs = {
  'chem_pi': chem_pi,
  'chem_ctl': chem_ctl,
  'oae_type': oae_type,
  'obj_var': obj_var,
  'cdreff': cdreff
  }
f_x = lambda x: f_dTA(x, **kwargs)
root = opt.brentq(f_x, x_upr, x_lwr)
print("The dTA needed for restoration to pre-industrial conditions is %.2f umol/kg" % (root))

EXAMPLE 2 - use of etamax function to process arrays of any shape

import numpy as np
from PyOAE import etamax
dTA = 1   # lower bound of the range of dTA values (umol/kg) to search for the root
TA = 2232
DIC = 1861
SiO3 = 1.346
PO4 = 0.201
Temp = 26.683
Sal = 34.004
Pres = 0
result = etamax(dTA, TA, DIC, SiO3, PO4, Temp, Sal, Pres)
result



