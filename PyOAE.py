# -*- coding: utf-8 -*-

__version__ = "1.0.20"

def f_dTA(dTA, **kwargs):

    """
    PURPOSE

    f_dTA calculates the difference between pre-industrial and treated conditions for the objective variable
  
    EXAMPLE USAGE

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

    INPUTS

    dTA is the trial value of the added concentration of TA due to OAE treatment (umol/kg)
    kwargs are the keyword arguments of other inputs to f_dTA as described below

    kwargs['chem_pi'] = vector of 7 values with chemistry variables for pre-industrial (e.g. year 1750)
    chem_pi[0] = TA umol/kg
    chem_pi[1] = DIC umol/kg
    chem_pi[2] = SiO3 umol/kg
    chem_pi[3] = PO4 umol/kg
    chem_pi[4] = temperature degC
    chem_pi[5] = salinity psu
    chem_pi[6] = pressure dbar

    kwargs['chem_ctl'] = vector of y values with chemistry variables for pre-treatment control (e.g. year 2010)
    chem_ctl[0] = TA umol/kg
    chem_ctl[1] = DIC umol/kg
    chem_ctl[2] = SiO3 umol/kg
    chem_ctl[3] = PO4 umol/kg
    chem_ctl[4] = temperature degC
    chem_ctl[5] = salinity (psu)
    chem_ctl[6] = pressure (dbar)

    kwargs['oae_type'] = chemical used for OAE, either 'NaOH' or 'Na2CO3'

    kwargs['obj_var'] = the objective variable, either 'alkstar', 'talk2dic', 'co3', 'hco3', 'omara', 'omcal', 'phtot', or 'revelle' 
    (which variable are we using to find the root)

    kwargs['cdreff'] = CDR efficiency = eta / etamax (e.g. cdreff=0.8 corresponds to 80% CDR efficiency)
 
    OUTPUTS 

    dVar = difference between pre-industrial and treated conditons for the objective variable (obj_var) 
    e.g. if obj_var is alkstar, then dVar = alkstar_trt - alkstar_pi 
    e.g. if obj_var is co3, then f_dTA = co3_trt - co3_pi 
    etc.
    (use Brent's method to find root where f_dTA equals zero)

    """

    import numpy as np
    import PyCO2SYS as pyco2

    # Define default values of input data arguments
    defaults = {
        'chem_pi': np.array([]),
        'chem_ctl': np.array([]),
        'oae_type': 'NaOH',
        'obj_var': 'alkstar',
        'cdreff': 0,
        'opt_pH_scale': 1,  # pH scale at which the input pH is reported ("1" means "Total Scale")
        'opt_k_carbonic': 10,  # Choice of H2CO3 and HCO3- dissociation constants K1 and K2 ("10" means "Lueker et al 2000")
        'opt_k_bisulfate': 1,  # Choice of HSO4- dissociation constant KSO4 ("1" means "Dickson")
        'opt_total_borate': 2,  # Choice of boron:sal ("2" means "Lee et al 2010")
        'opt_k_fluoride': 2   # "2" means Perez and Fraga 1987        
        }

    # Update input data argumements with any provided keyword arguments in kwargs
    data = {**defaults, **kwargs}

    # unpack updated input data keyword arguments
    chem_pi = data["chem_pi"]
    chem_ctl = data["chem_ctl"]
    oae_type = data["oae_type"]
    obj_var = data["obj_var"]
    cdreff = data["cdreff"]

    # error trapping for zero or negative input values of dTA
    if dTA <= 0:
        dTA = 1e-3
    
    # - - -
    # Step 1: find pCO2_ctl = pCO2 of control before treatment

    # Define input conditions for PyCO2SYS using same constants options as Jiang et al 2023
    TA_ctl = chem_ctl[0]
    DIC_ctl = chem_ctl[1]
    kwargs = dict(
        par1_type = 1,  # The first parameter supplied is of type "1", which means "alkalinity"
        par1 = TA_ctl,  # value of the first parameter is TA = TAt,ctl
        par2_type = 2,  # The second parameter supplied is of type "2", which means "DIC"
        par2 = DIC_ctl,  # value of the second parameter DIC = DICt,ctl
        total_silicate = chem_ctl[2],  # Concentration of silicate  in the sample (in umol/kg)
        total_phosphate = chem_ctl[3],  # Concentration of phosphate in the sample (in umol/kg)
        temperature = chem_ctl[4],  # Temperature at input conditions
        salinity = chem_ctl[5],  # Salinity of the sample
        pressure = chem_ctl[6],  # Pressure    at input conditions
        opt_pH_scale = data["opt_pH_scale"],  # Choice of pH scale
        opt_k_carbonic = data["opt_k_carbonic"],  # Choice of H2CO3 and HCO3- dissociation constants
        opt_k_bisulfate = data["opt_k_bisulfate"],  # Choice of HSO4- dissociation constant KSO4
        opt_total_borate = data["opt_total_borate"],  # Choice of boron:sal
        opt_k_fluoride = data["opt_k_fluoride"]   # Choice of hydrogen fluoride dissociation constant
        )
    # Run PyCO2SYS
    results = pyco2.sys(**kwargs)
    pCO2_ctl = results['pCO2']

    # - - -
    # Step 2: find etamax and DIC_trt = DIC of treated conditions 

    # Define input conditions for PyCO2SYS using same constants options as Jiang et al 2023
    TA_trt = TA_ctl + dTA
    kwargs = dict(
        par1_type = 1,  # The first parameter supplied is of type "1", which means "alkalinity"
        par1 = TA_trt,  # value of the first parameter is TA = TAt,trt
        par2_type = 4,  # The second parameter supplied is of type "4", which means "pCO2"
        par2 = pCO2_ctl,  # value of the second parameter is pCO2 = pCO2t,ctl
        total_silicate = chem_ctl[2],  # Concentration of silicate  in the sample (in umol/kg)
        total_phosphate = chem_ctl[3],  # Concentration of phosphate in the sample (in umol/kg)
        temperature = chem_ctl[4],  # Temperature at input conditions
        salinity = chem_ctl[5],  # Salinity of the sample
        pressure = chem_ctl[6],  # Pressure    at input conditions
        opt_pH_scale = data["opt_pH_scale"],  # Choice of pH scale
        opt_k_carbonic = data["opt_k_carbonic"],  # Choice of H2CO3 and HCO3- dissociation constants
        opt_k_bisulfate = data["opt_k_bisulfate"],  # Choice of HSO4- dissociation constant KSO4
        opt_total_borate = data["opt_total_borate"],  # Choice of boron:sal
        opt_k_fluoride = data["opt_k_fluoride"]   # Choice of hydrogen fluoride dissociation constant
        )
    # Run PyCO2SYS
    results = pyco2.sys(**kwargs)
    DICeq = results['dic']   # equilibrium DIC at TAt,trt and pCO2t,ctl (umol/kg)
    CDRpot = DICeq - DIC_ctl   # CDR potential (umol/kg)
    etamax = CDRpot / dTA
    if oae_type == "NaOH":
        dDICtrt = 0
        dDICcdr = cdreff * etamax * dTA
    elif oae_type == "Na2CO3":
        dDICtrt = 0.5 * dTA
        dDICcdr =  cdreff * (etamax - 0.5) * dTA
    DIC_trt = DIC_ctl + dDICtrt + dDICcdr         

    # - - -
    # Step 3: calculate the objective function f_dTA that should be zero

    if obj_var == "alkstar":

    	TA_pi = chem_pi[0]
    	DIC_pi = chem_pi[1]
    	alkstar_pi = TA_pi - DIC_pi;
    	alkstar_trt = TA_trt - DIC_trt;
    	dVar = alkstar_trt - alkstar_pi;

    else:

        # pre-industrial conditions
        kwargs = dict(
            par1_type = 1,  # The first parameter supplied is of type "1", which means "alkalinity"
            par1 = chem_pi[0],  # TA_pi
            par2_type = 2,  # The second parameter supplied is of type "4", which means "DIC"
            par2 = chem_pi[1],  # DIC_pi
            total_silicate = chem_pi[2],  # Concentration of silicate  in the sample (in umol/kg)
            total_phosphate = chem_pi[3],  # Concentration of phosphate in the sample (in umol/kg)
            temperature = chem_pi[4],  # Temperature at input conditions
            salinity = chem_pi[5],  # Salinity of the sample
            pressure = chem_pi[6],  # Pressure    at input conditions
            opt_pH_scale = data["opt_pH_scale"],  # Choice of pH scale
            opt_k_carbonic = data["opt_k_carbonic"],  # Choice of H2CO3 and HCO3- dissociation constants
            opt_k_bisulfate = data["opt_k_bisulfate"],  # Choice of HSO4- dissociation constant KSO4
            opt_total_borate = data["opt_total_borate"],  # Choice of boron:sal
            opt_k_fluoride = data["opt_k_fluoride"]   # Choice of hydrogen fluoride dissociation constant
            )
        # Run PyCO2SYS
        results = pyco2.sys(**kwargs)
        phtot_pi = results['pH_total'] 	# pH total scale for CO2SYS "input" condition
        co3_pi = results['carbonate'] 	# CO3^2- umol/kg
        omcal_pi = results['saturation_calcite']  	# Calcite saturation state Omega for CO2SYS "input" condition
        omara_pi = results['saturation_aragonite']  	# Aragonite saturation state Omega for CO2SYS "input" condition

        # treated condtion at time t
        TA_trt = TA_ctl + dTA
        kwargs = dict(
            par1_type = 1,  # The first parameter supplied is of type "1", which means "alkalinity"
            par1 = TA_trt, 
            par2_type = 2,  # The second parameter supplied is of type "2", which means "DIC"
            par2 = DIC_trt,  
            total_silicate = chem_ctl[2],  # Concentration of silicate  in the sample (in umol/kg)
            total_phosphate = chem_ctl[3],  # Concentration of phosphate in the sample (in umol/kg)
            temperature = chem_ctl[4],  # Temperature at input conditions
            salinity = chem_ctl[5],  # Salinity of the sample
            pressure = chem_ctl[6],  # Pressure    at input conditions
            opt_pH_scale = data["opt_pH_scale"],  # Choice of pH scale
            opt_k_carbonic = data["opt_k_carbonic"],  # Choice of H2CO3 and HCO3- dissociation constants
            opt_k_bisulfate = data["opt_k_bisulfate"],  # Choice of HSO4- dissociation constant KSO4
            opt_total_borate = data["opt_total_borate"],  # Choice of boron:sal
            opt_k_fluoride = data["opt_k_fluoride"]   # Choice of hydrogen fluoride dissociation constant
            )
        # Run PyCO2SYS
        results = pyco2.sys(**kwargs)
        phtot_trt = results['pH_total'] 	# pH total scale for CO2SYS "input" condition
        co3_trt = results['carbonate'] 	# CO3^2- umol/kg
        omcal_trt = results['saturation_calcite']  	# Calcite saturation state Omega for CO2SYS "input" condition
        omara_trt = results['saturation_aragonite']  	# Aragonite saturation state Omega for CO2SYS "input" condition

        # calculate dVar
        if obj_var == 'phtot':
            dVar = phtot_trt - phtot_pi
        elif obj_var == 'co3':
            dVar = co3_trt - co3_pi
        elif obj_var == 'omara':
            dVar = omara_trt - omara_pi
        elif obj_var == 'omcal':
            dVar = omcal_trt - omcal_pi

    return dVar

def etamax(dTA, **kwargs):

    """
    PURPOSE

    Calculate the hypothetical maximum OAE efficiency "etamax" for an assumed alkalinity treatment addition "dTA"
    and control conditions of TA, DIC, SiO3, PO4, temperature, and salinity
  
    EXAMPLE USAGE

    # import the packages that are needed
    import numpy as np
    from PyOAE import etamax
    # assign the inputs that are needed
    dTA = 1   # lower bound of the range of dTA values (umol/kg) to search for the root
    kwargs = dict(
        TA_ctl = 2232,
        DIC_ctl = 1861,
        SiO3_ctl = 1.346,
        PO4_ctl = 0.201,
        Temp_ctl = 26.683,
        Sal_ctl = 34.004,
        Pres_ctl = 0
        )
    # call the etamax function
    result = etamax(dTA, **kwargs)
    # display the result dictionary
    result

    INPUTS

    dTA = assumed value of the added concentration of TA due to OAE treatment (umol/kg)
    TA = control TA before OAE addition (umol/kg)
    DIC = control DIC before OAE addition (umol/kg)
    SiO3 = control SiO3 (umol/kg)
    PO4 = control PO4 (umol/kg)
    Temp = control temperature (degC)
    Sal = control salinity (psu)
    Pres = control pressure (dbar)

    OUTPUTS 

    - result = dictionary of output varlables with the following keys:
        - 'etamax': hypothetical maximum OAE efficiency (dimensionless)
        - 'CDR_pot': hypothetical maximum potential CDR (umol/kg)
        - 'DIC_eq': DIC at equilibrium when input TA= TA_ctl + dTA and input pCO2= pCO2t_ctl 
        - 'pCO2_ctl': pCO2 at control conditions (uatm)

    """

    import numpy as np
    import PyCO2SYS as pyco2

    # Define default values of input data arguments
    defaults = {
        'TA_ctl': np.array([]),
        'DIC_ctl': np.array([]),
        'SiO3_ctl': np.array([]),
        'PO4_ctl': np.array([]),
        'Temp_ctl': np.array([]),
        'Sal_ctl': np.array([]),
        'Pres_ctl': np.array([]),
        'opt_pH_scale': 1,  # pH scale at which the input pH is reported ("1" means "Total Scale")
        'opt_k_carbonic': 10,  # Choice of H2CO3 and HCO3- dissociation constants K1 and K2 ("10" means "Lueker et al 2000")
        'opt_k_bisulfate': 1,  # Choice of HSO4- dissociation constant KSO4 ("1" means "Dickson")
        'opt_total_borate': 2,  # Choice of boron:sal ("1" means "Lee et al 2010")
        'opt_k_fluoride': 2   # "2" means Perez and Fraga 1987        
        }

    # Update input data argumements with any provided keyword arguments in kwargs
    data = {**defaults, **kwargs}

    # unpack updated input data keyword arguments
    TA_ctl = data["TA_ctl"]
    DIC_ctl = data["DIC_ctl"]
    SiO3_ctl = data["SiO3_ctl"]
    PO4_ctl = data["PO4_ctl"]
    Temp_ctl = data["Temp_ctl"]
    Sal_ctl = data["Sal_ctl"]
    Pres_ctl = data["Pres_ctl"]

    # error trapping for zero or negative input values of dTA
    if dTA <= 0:
        dTA = 1e-3
    
    # - - -
    # Step 1: find pCO2_ctl = pCO2 of control before treatment

    # Define input conditions for PyCO2SYS using same constants options as Jiang et al 2023
    kwargs = dict(
        par1_type = 1,  # The first parameter supplied is of type "1", which means "alkalinity"
        par1 = TA_ctl,  # value of the first parameter is TA = TAt,ctl
        par2_type = 2,  # The second parameter supplied is of type "2", which means "DIC"
        par2 = DIC_ctl,  # value of the second parameter DIC = DICt,ctl
        total_silicate = SiO3_ctl,  # Concentration of silicate  in the sample (in umol/kg)
        total_phosphate = PO4_ctl,  # Concentration of phosphate in the sample (in umol/kg)
        temperature = Temp_ctl,  # Temperature at input conditions
        salinity = Sal_ctl,  # Salinity of the sample
        pressure = Pres_ctl,  # Pressure    at input conditions
        opt_pH_scale = data["opt_pH_scale"],  # Choice of pH scale
        opt_k_carbonic = data["opt_k_carbonic"],  # Choice of H2CO3 and HCO3- dissociation constants
        opt_k_bisulfate = data["opt_k_bisulfate"],  # Choice of HSO4- dissociation constant KSO4
        opt_total_borate = data["opt_total_borate"],  # Choice of boron:sal
        opt_k_fluoride = data["opt_k_fluoride"]   # Choice of hydrogen fluoride dissociation constant
        )
    # Run PyCO2SYS
    results = pyco2.sys(**kwargs)
    pCO2_ctl = results['pCO2']

    # - - -
    # Step 2: find etamax and DIC_trt = DIC of treated conditions 

    # Define input conditions for PyCO2SYS using same constants options as Jiang et al 2023
    TA_trt = TA_ctl + dTA
    kwargs = dict(
        par1_type = 1,  # The first parameter supplied is of type "1", which means "alkalinity"
        par1 = TA_trt,  # value of the first parameter is TA = TAt,trt
        par2_type = 4,  # The second parameter supplied is of type "4", which means "pCO2"
        par2 = pCO2_ctl,  # value of the second parameter is pCO2 = pCO2t,ctl
        total_silicate = SiO3_ctl,  # Concentration of silicate  in the sample (in umol/kg)
        total_phosphate = PO4_ctl,  # Concentration of phosphate in the sample (in umol/kg)
        temperature = Temp_ctl,  # Temperature at input conditions
        salinity = Sal_ctl,  # Salinity of the sample
        pressure = Pres_ctl,  # Pressure    at input conditions
        opt_pH_scale = data["opt_pH_scale"],  # Choice of pH scale
        opt_k_carbonic = data["opt_k_carbonic"],  # Choice of H2CO3 and HCO3- dissociation constants
        opt_k_bisulfate = data["opt_k_bisulfate"],  # Choice of HSO4- dissociation constant KSO4
        opt_total_borate = data["opt_total_borate"],  # Choice of boron:sal
        opt_k_fluoride = data["opt_k_fluoride"]   # Choice of hydrogen fluoride dissociation constant
        )
    # Run PyCO2SYS
    results = pyco2.sys(**kwargs)
    DIC_eq = results['dic']   # equilibrium DIC at TAt,trt and pCO2t,ctl (umol/kg)
    CDR_pot = DIC_eq - DIC_ctl   # CDR potential (umol/kg)
    etamax = CDR_pot / dTA

    # make the result dictionary for output
    result = {
            'etamax': etamax,
            'CDR_pot': CDR_pot,
            'DIC_eq': DIC_eq,
            'pCO2_ctl': pCO2_ctl
            }
    
    return result

def dic_bio(**kwargs):

    """
    PURPOSE

    Calculate the biological component of DIC (DIC_bio) 
    as the difference between the observed DIC (DIC_obs) 
    compared with the DIC at equilibrium with atmospheric pCO2 (DIC_atm)

    where

    DIC_bio = DIC_obs - DIC_atm
  
    INPUTS

    kwargs = dict(
        alkalinity = ,  # Concentration of alkainity in the sample (in umol/kg)
        dic = ,  # Concentration of DIC in the sample (in umol/kg)
        pco2atm = ,  # atmospheric pCO2 or fCO2 (in umol/kg)
        pco2atm_type = ,  # Atmospheric pco2atm is as 4=pCO2 or 5=fCO2
        total_silicate = ,  # Concentration of silicate  in the sample (in umol/kg)
        total_phosphate = ,  # Concentration of phosphate in the sample (in umol/kg)
        temperature = , # Temperature at input conditions
        salinity = , # Salinity of the sample
        pressure = , # Pressure    at input conditions
        opt_pH_scale = ,  # Choice of pH scale
        opt_k_carbonic = ,  # Choice of H2CO3 and HCO3- dissociation constants
        opt_k_bisulfate = ,  # Choice of HSO4- dissociation constant KSO4
        opt_total_borate = ,  # Choice of boron:sal
        opt_k_fluoride =    # Choice of hydrogen fluoride dissociation constant
        )

    OUTPUTS 

    - result = dictionary of output varlables with the following keys:
        - 'dic_atm': DIC at equilibrium with atmospheric pCO2 (umol/kg)
        - 'dic_bio': Biological component of DIC (umol/kg)

    """

    import numpy as np
    import PyCO2SYS as pyco2

    # Define default values of input data arguments
    defaults = {
        'pco2atm_type': 4,  # The second parameter 4=pCO2, 5=fCO2
        'alkalinity': np.array([]),
        'dic': np.array([]),
        'pco2atm': np.array([]),
        'total_silicate': np.array([]),
        'total_phosphate': np.array([]),
        'temperature': np.array([]),
        'salinity': np.array([]),
        'pressure': np.array([]),
        'opt_pH_scale': 1,  # pH scale at which the input pH is reported ("1" means "Total Scale")
        'opt_k_carbonic': 10,  # Choice of H2CO3 and HCO3- dissociation constants K1 and K2 ("10" means "Lueker et al 2000")
        'opt_k_bisulfate': 1,  # Choice of HSO4- dissociation constant KSO4 ("1" means "Dickson")
        'opt_total_borate': 2,  # Choice of boron:sal ("1" means "Lee et al 2010")
        'opt_k_fluoride': 2   # "2" means Perez and Fraga 1987        
        }

    # Update input data argumements with any provided keyword arguments in kwargs
    data = {**defaults, **kwargs}

    # - - -
    # Step 1: find pCO2_ctl = pCO2 of control before treatment

    # Define input conditions for PyCO2SYS using same constants options as Jiang et al 2023
    kwargs = dict(
        par1_type = 1,  # The first parameter supplied is of type "1", which means "alkalinity"
        par1 = data["alkalinity"],  # value of the first parameter
        par2_type = data["pco2atm_type"],  # The second parameter 4=pCO2, 5=fCO2
        par2 = data["pco2atm"],  # value of the second parameter
        total_silicate = data["total_silicate"],  # Concentration of silicate  in the sample (in umol/kg)
        total_phosphate = data["total_phosphate"],  # Concentration of phosphate in the sample (in umol/kg)
        temperature = data["temperature"],  # Temperature at input conditions
        salinity = data["salinity"],  # Salinity of the sample
        pressure = data["total_pressure"],  # Pressure    at input conditions
        opt_pH_scale = data["opt_pH_scale"],  # Choice of pH scale
        opt_k_carbonic = data["opt_k_carbonic"],  # Choice of H2CO3 and HCO3- dissociation constants
        opt_k_bisulfate = data["opt_k_bisulfate"],  # Choice of HSO4- dissociation constant KSO4
        opt_total_borate = data["opt_total_borate"],  # Choice of boron:sal
        opt_k_fluoride = data["opt_k_fluoride"]   # Choice of hydrogen fluoride dissociation constant
        )
    # Run PyCO2SYS
    results = pyco2.sys(**kwargs)
    dic_atm = results['dic']
    dic_bio = data["dic"] - dic_atm
    
    # make the result dictionary for output
    result = {
            'dic_atm': dic_atm,
            'dic_bio': dic_bio
            }
    
    return result
    
def sine_fit(x,y):

    """
    PURPOSE

    Sine-regression to estimate the mean and amplitude of the annual cycles of a periodic variable, 
    for example DIC_bio, or any other periodic variable with repeating cycles per unit of time
    assuming period = 1

    y = mean + amplitude * sin(2 * pi * (x - phase) / period)

    For example, in the case of annual cycles of DIC_bio:
    
    y = DICbio = DICobs – DICatm
    x = time as decimal year fraction (1982-2022) 
    mean = mean DICbio from regression
    amplitude = amplitude of DICbio from regression
    phase = phase shift of DICbio from regression
    period = assumed equal to 1 cycle per year for DICbio

    INPUT

    x = time (e.g. decimal years for annual periodic data)
    y = observed periodic data (e.g. y = DIC_bio = DIC_obs - DIC_atm)

    OUTPUT

    result = dictionary of output results of mean, amplitude, phase,
    x (echo of input), y (echo of input), y_fit, rmse, popt, pcov,
    and various other regression statistics

    """
    
    from scipy.optimize import curve_fit
    import numpy as np
    from scipy import stats
    import sys

    ctrl = np.isreal(x).all() and (not np.isnan(x).any()) and (not np.isinf(x).any()) and x.ndim==1
    if not ctrl:
      print('Check x: it needs be a vector of real numbers with no infinite or nan values!','\n')
      sys.exit()
    ctrl = np.isreal(y).all() and (not np.isnan(y).any()) and (not np.isinf(y).any()) and y.ndim==1
    if not ctrl:
      print('Check y: it needs be a vector of real numbers with no infinite or nan values!','\n')
      sys.exit()
    ctrl =  np.size(x)==np.size(y)
    if not ctrl:
      print('Check x and y: x and y need to be the same size!','\n')
      sys.exit()
    
    f_x = lambda x,a,b,c: a + b * np.sin(2 * np.pi * (x - c))

    # Initial guess for the parameters [A, B, C]
    # A= vertical offset (mean)
    # B= amplitude
    # C= phase shift
    initial_guess = [np.mean(y), np.std(y), 0.5]

    # Perform the curve fitting to find 
    # popt= optimum parameter set, and pcov = covariance matrix
    popt, pcov = curve_fit(f_x, x, y, p0=initial_guess)

    # Extract the fitted parameters
    A_fit, B_fit, C_fit = popt

    # Generate y values using the fitted parameters
    y_fit = f_x(x, A_fit, B_fit, C_fit)

    # - - -
    # additional outputs of regression statistics

    nobs = np.size(x)
    nparam = np.size(popt)
    df = nobs - nparam
      
    SSE = np.sum((y-y_fit) ** 2)                 # sum of squares (residual error)
    MSE = SSE / df                              # mean square (residual error)
    syx = np.sqrt(MSE)                          # std error of the estimate
    rmse = np.sqrt(SSE / nobs)                  # root mean squared error
        
    SST = np.sum(y **2) - np.sum(y) **2 / nobs  # sum of squares (total)
    SSR = SST - SSE                             # sum of squares (regression model)
    MSR = SSR / (np.size(popt)-1)              # mean square (regression model)
    Fstat = MSR / MSE           # F statistic
    dfn = np.size(popt) - 1    # df numerator = degrees of freedom for model = number of model parameters - 1
    dfd = df                    # df denomenator = degrees of freedom of the residual = df = nobs - nparam
    pvalue = 1-stats.f.cdf(Fstat, dfn, dfd)      # p-value of F test statistic
    rsquared = SSR / SST                                                        # ordinary rsquared
    adj_rsquared = 1-(1-rsquared)*(np.size(x)-1)/(np.size(x)-np.size(popt)-1)  # adjusted rsquared
    
    # make the result dictionary for output
    result = {
            'mean': A_fit,
            'amplitude': B_fit,
            'phase': C_fit,
            'x': x,
            'y': y,
            'y_fit': y_fit,
            'rmse': rmse,
            'popt': popt,
            'pcov': pcov,
            'SST': SST,
            'SSR': SSR,
            'SSE': SSE,
            'MSR': MSR,
            'MSE': MSE,
            'syx': syx,
            'nobs': nobs,
            'nparam': nparam,
            'df': df,
            'Fstat': Fstat,
            'dfn': dfn,
            'dfd': dfd,
            'pvalue': pvalue,
            'rsquared': rsquared,
            'adj_rsquared': adj_rsquared
            }

    return result

def nnn(x):

    """
    PURPOSE
    Count the number of non-nan values in the numpy array x
    USAGE
    result = nnn(x)
    INPUT
    x = any numpy array of any dimension
    OUTPUT
    result = number of non-nan values in the array x
    """
    
    import numpy as np

    result = np.count_nonzero(~np.isnan(x))
    
    return result

def pco2_tnorm(pco2_obs,temp_obs):

    """
    PURPOSE

    Temperature normalization of pCO2 using the equation of 
    Takahashi (2002) (http://dx.doi.org/10.1016/S0967-0645(02)00003-6)

    In order to remove the temperature effect from the observed pCO2,
    the observed pCO2 values are normalized to the long-term mean
    temperature of seawater at the location, using the following equation:

    pCO2_Tmean = pCO2_obs * exp(0.0423 * (Tmean - T_obs))     (eqn 1)

    In a similar manner, the effect of temperature changes on pCO2 is computed
    by perturbing the mean annual pCO2. The adjusted mean pCO2 value at the 
    set of observed temperatures, Tobs, is computed using the following equation:

    pCO2_Tobs = pCO2_mean * exp(0.0423 * (Tobs - Tmean))     (eqn 2)

    INPUT

    pco2_obs = observed pCO2 (e.g. monthly average pCO2 from OceanSODA-ETHZ) (uatm)
    temp_obs = observed temperature co-inciding with observed pCO2 (degC)

    OUTPUT

    result = temperature-normalized pCO2 using eqn 1 (uatm)

    """

    import numpy as np
    import sys
    
    x = pco2_obs
    y = temp_obs
    ctrl = np.isreal(x).all() and (not np.isnan(x).any()) and (not np.isinf(x).any()) and x.ndim==1
    if not ctrl:
      print('Check pco2 input: it needs be a vector of real numbers with no infinite or nan values!','\n')
      sys.exit()
    ctrl = np.isreal(y).all() and (not np.isnan(y).any()) and (not np.isinf(y).any()) and y.ndim==1
    if not ctrl:
      print('Check temperature input: it needs be a vector of real numbers with no infinite or nan values!','\n')
      sys.exit()
    ctrl =  np.size(x)==np.size(y)
    if not ctrl:
      print('Check pco2 and temperature: both need to be the same size!','\n')
      sys.exit()

    temp_mean = np.nanmean(temp_obs)
    pco2_tmean = pco2_obs * np.exp(0.0423 * (temp_mean - temp_obs))

    pco2_mean = np.nanmean(pco2_obs)
    pco2_tobs = pco2_mean * np.exp(0.0423 * (temp_obs - temp_mean))

    result = {
        'pco2_tmean': pco2_tmean,
        'pco2_tobs': pco2_tobs
        }

    return result

def find_nearest_index_1d(arr, target):

    """
    Find the index of the nearest value to the target in the 1-d array

    Parameters:
    arr (list of float): The list of numbers.
    target (float): The target number.

    Returns:
    int: The index of the nearest value.

    USAGE
    numbers = [10, 22, 14, 3, 7, 9]
    target_value = 8
    index = find_nearest_index_1d(numbers, target_value)
    print(f"The index of the nearest value to {target_value} is {index}.")
    """
    nearest_index = min(range(len(arr)), key=lambda i: abs(arr[i] - target))
    return nearest_index

def find_nearest_index_2d(array, target):

    """
    PURPOSE
    Find the row,col index of value in the input 2-d "array" that is nearest to the input "target"
    
    USAGE
    array = [
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9]
        ]
    target = 7.6
    row,col = find_nearest_index_2d(array, target)
    print("The (row,col) index of the nearest value to the target %.1f is (%.0f,%.0f)" % (value,row,col))    
    """
    import numpy as np
    
    # Convert the array to a numpy array for easier manipulation
    array = np.array(array)
    
    # Flatten the array and find the index of the nearest value
    flat_index = np.abs(array - target).argmin()
    
    # Convert the flat index back to 2D coordinates
    row, col = divmod(flat_index, array.shape[1])
    
    return row, col

def pco2_fass(**kwargs):

    """
    PURPOSE

    Fassbender's method of separating thermal and non-thermal components of pCO2
    Rodgers et al 2022 (https://doi.org/10.1029/2023GB007798)

    INPUT 

    kwargs dictionary with the following keyword arguments:
    'alkalinity' = observed monthly alkalinity (umol/kg)
    'dic' = observed monthly DIC (umol/kg)
    'total_silicate' = observed monthly total silicate (umol/kg)
    'total_phosphate' = observed monthly alkalinity (umol/kg)
    'temperature' = observed monthly temperature (degC)
    'salinity' = observed monthly salinity (psu)
    'pressure' = observed monthly pressure (dbar)
    'opt_pH_scale' = PyCO2SYS option
    'opt_k_carbonic' = PyCO2SYS option
    'opt_k_bisulfate' = PyCO2SYS option
    'pt_total_borate' = PyCO2SYS option
    'opt_k_fluoride' = PyCO2SYS option

    OUTPUT

    result dictionary with the following keywords:
    'pCO2' = observed monthly pCO2 (uatm)
    'pCO2_AM' = mean pCO2 estimated using 
        annual mean TA, DIC, SIO4, PO4, T, S (uatm)
    'pCO2_TFASS' = monthly thermally-driven pCO2 using 
        monthly T with mean TA, DIC, SIO4, PO4, S (uatm)
    'pCO2_Tanom' = thermal pCO2 component seasonal cycle anomaly (uatm)
    'pCO2_NTFASS' = monthly non-thermal component of pCO2 (uatm)
    """

    import numpy as np
    import PyCO2SYS as pyco2

    # Define default values of input data arguments
    defaults = {
        'alkalinity': np.array([]),
        'dic': np.array([]),
        'total_silicate': np.array([]),
        'total_phosphate': np.array([]),
        'temperature': np.array([]),
        'salinity': np.array([]),
        'pressure': np.array([]),
        'opt_pH_scale': 1,  # pH scale at which the input pH is reported ("1" means "Total Scale")
        'opt_k_carbonic': 10,  # Choice of H2CO3 and HCO3- dissociation constants K1 and K2 ("10" means "Lueker et al 2000")
        'opt_k_bisulfate': 1,  # Choice of HSO4- dissociation constant KSO4 ("1" means "Dickson")
        'opt_total_borate': 2,  # Choice of boron:sal ("1" means "Lee et al 2010")
        'opt_k_fluoride': 2   # "2" means Perez and Fraga 1987        
        }

    # Update input data argumements with any provided keyword arguments in kwargs
    data = {**defaults, **kwargs}

    # - - -
    # Step 1: solve for observed monthly pCO2 from observed monthly TA,DIc,SIO4,PO4,T,S
    kwargs = dict(
        par1_type = 1,  # The first parameter supplied is of type "1", which means "alkalinity"
        par1 = data["alkalinity"],  # value of the first parameter
        par2_type = 2,  # The second parameter 2=DIC, 4=pCO2, 5=fCO2
        par2 = data["dic"],  # value of the second parameter
        total_silicate = data["total_silicate"],  # Concentration of silicate  in the sample (in umol/kg)
        total_phosphate = data["total_phosphate"],  # Concentration of phosphate in the sample (in umol/kg)
        temperature = data["temperature"],  # Temperature at input conditions
        salinity = data["salinity"],  # Salinity of the sample
        pressure = data["total_pressure"],  # Pressure    at input conditions
        opt_pH_scale = data["opt_pH_scale"],  # Choice of pH scale
        opt_k_carbonic = data["opt_k_carbonic"],  # Choice of H2CO3 and HCO3- dissociation constants
        opt_k_bisulfate = data["opt_k_bisulfate"],  # Choice of HSO4- dissociation constant KSO4
        opt_total_borate = data["opt_total_borate"],  # Choice of boron:sal
        opt_k_fluoride = data["opt_k_fluoride"]   # Choice of hydrogen fluoride dissociation constant
        )
    # Run PyCO2SYS
    results = pyco2.sys(**kwargs)
    pCO2 = results['pCO2']

    # - - -
    # Step 2: solve for pCO2_AM (eqn 4 in Rodgers et al)
    Tbar = np.nanmean(data["temperature"])
    Sbar = np.nanmean(data["salinity"])
    DICbar = np.nanmean(data["dic"])
    TAbar = np.nanmean(data["alkalinity"])
    PO4bar = np.nanmean(data["total_phosphate"])
    SIO4bar = np.nanmean(data["total_silicate"])
    kwargs = dict(
        par1_type = 1,  # The first parameter supplied is of type "1", which means "alkalinity"
        par1 = TAbar,  # value of the first parameter
        par2_type = 2,  # The second parameter 2=DIC, 4=pCO2, 5=fCO2
        par2 = DICbar,  # value of the second parameter
        total_silicate = SIO4bar,  # Concentration of silicate  in the sample (in umol/kg)
        total_phosphate = PO4bar,  # Concentration of phosphate in the sample (in umol/kg)
        temperature = Tbar,  # Temperature at input conditions
        salinity = Sbar,  # Salinity of the sample
        pressure = data["total_pressure"],  # Pressure    at input conditions
        opt_pH_scale = data["opt_pH_scale"],  # Choice of pH scale
        opt_k_carbonic = data["opt_k_carbonic"],  # Choice of H2CO3 and HCO3- dissociation constants
        opt_k_bisulfate = data["opt_k_bisulfate"],  # Choice of HSO4- dissociation constant KSO4
        opt_total_borate = data["opt_total_borate"],  # Choice of boron:sal
        opt_k_fluoride = data["opt_k_fluoride"]   # Choice of hydrogen fluoride dissociation constant
        )
    # Run PyCO2SYS
    results = pyco2.sys(**kwargs)
    pCO2_AM = results['pCO2']
 
    # - - -
    # Step 3: solve for pCO2_TFASS (eqn 5 in Rodgers et al)
    kwargs = dict(
        par1_type = 1,  # The first parameter supplied is of type "1", which means "alkalinity"
        par1 = TAbar,  # value of the first parameter
        par2_type = 2,  # The second parameter 2=DIC, 4=pCO2, 5=fCO2
        par2 = DICbar,  # value of the second parameter
        total_silicate = SIO4bar,  # Concentration of silicate  in the sample (in umol/kg)
        total_phosphate = PO4bar,  # Concentration of phosphate in the sample (in umol/kg)
        temperature = data["temperature"],  # Temperature at input conditions
        salinity = Sbar,  # Salinity of the sample
        pressure = data["total_pressure"],  # Pressure    at input conditions
        opt_pH_scale = data["opt_pH_scale"],  # Choice of pH scale
        opt_k_carbonic = data["opt_k_carbonic"],  # Choice of H2CO3 and HCO3- dissociation constants
        opt_k_bisulfate = data["opt_k_bisulfate"],  # Choice of HSO4- dissociation constant KSO4
        opt_total_borate = data["opt_total_borate"],  # Choice of boron:sal
        opt_k_fluoride = data["opt_k_fluoride"]   # Choice of hydrogen fluoride dissociation constant
        )
    # Run PyCO2SYS
    results = pyco2.sys(**kwargs)
    pCO2_TFASS = pCO2

    # - - -
    # Step 4: solve for pCO2_Tanom (eqn 6 in Rodgers et al)
    pCO2_Tanom = pCO2_TFASS - pCO2_AM

    # - - -
    # Step 5: solve for pCO2_NTFASS (eqn 7 in Rodgers et al)
    pCO2_NTFASS = data["pCO2"] - pCO2_Tanom
 
    # make the result dictionary for output
    result = {
            'pCO2': pCO2,
            'pCO2_AM': pCO2_AM,
            'pCO2_TFASS': pCO2_TFASS,
            'pCO2_TFASS': pCO2_Tanom,
            'pCO2_TFASS': pCO2_NTFASS
            }

    return result





