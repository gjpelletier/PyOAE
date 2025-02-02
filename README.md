# PyOAE
Python tools for analysis of Ocean Alkalinity Enhancement

by Greg Pelletier

Tools for analysis of Ocean Alkalinity Enhancement, including Jupyter Notebooks for the following examples:

- Root-finding method to solve for the OAE treatment needed to restore any carbonate system variable (e.g. TA-DIC, CO3--, pH, Ωar, etc.) to pre-industrial conditions in the coastal California Current Ecosystem. This notebook can also be used to solve for any other location in the global oceans
- Calculation of the maximum hypothetical OAE efficiency etamax in the global oceans

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




