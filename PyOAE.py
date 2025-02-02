# -*- coding: utf-8 -*-

__version__ = "1.0.12"

def f_dTA(dTA, kwargs):

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

    # unpack the kwargs
    chem_pi = kwargs["chem_pi"]
    chem_ctl = kwargs["chem_ctl"]
    oae_type = kwargs["oae_type"]
    obj_var = kwargs["obj_var"]
    cdreff = kwargs["cdreff"]

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
        opt_pH_scale = 1,  # pH scale at which the input pH is reported ("1" means "Total Scale")
        opt_k_carbonic = 10,  # Choice of H2CO3 and HCO3- dissociation constants K1 and K2 ("10" means "Lueker et al 2000")
        opt_k_bisulfate = 1,  # Choice of HSO4- dissociation constant KSO4 ("1" means "Dickson")
        opt_total_borate = 2,  # Choice of boron:sal ("1" means "Lee et al 2010")
        opt_k_fluoride = 2   # "2" means Perez and Fraga 1987
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
        opt_pH_scale = 1,  # pH scale at which the input pH is reported ("1" means "Total Scale")
        opt_k_carbonic = 10,  # Choice of H2CO3 and HCO3- dissociation constants K1 and K2 ("10" means "Lueker et al 2000")
        opt_k_bisulfate = 1,  # Choice of HSO4- dissociation constant KSO4 ("1" means "Dickson")
        opt_total_borate = 2,  # Choice of boron:sal ("1" means "Lee et al 2010")
        opt_k_fluoride = 2   # "2" means Perez and Fraga 1987
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
            opt_pH_scale = 1,  # pH scale at which the input pH is reported ("1" means "Total Scale")
            opt_k_carbonic = 10,  # Choice of H2CO3 and HCO3- dissociation constants K1 and K2 ("10" means "Lueker et al 2000")
            opt_k_bisulfate = 1,  # Choice of HSO4- dissociation constant KSO4 ("1" means "Dickson")
            opt_total_borate = 2,  # Choice of boron:sal ("1" means "Lee et al 2010")
            opt_k_fluoride = 2   # "2" means Perez and Fraga 1987
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
            opt_pH_scale = 1,  # pH scale at which the input pH is reported ("1" means "Total Scale")
            opt_k_carbonic = 10,  # Choice of H2CO3 and HCO3- dissociation constants K1 and K2 ("10" means "Lueker et al 2000")
            opt_k_bisulfate = 1,  # Choice of HSO4- dissociation constant KSO4 ("1" means "Dickson")
            opt_total_borate = 2,  # Choice of boron:sal ("1" means "Lee et al 2010")
            opt_k_fluoride = 2   # "2" means Perez and Fraga 1987
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

def etamax(dTA, TA_ctl, DIC_ctl, SiO3_ctl, PO4_ctl, Temp_ctl, Sal_ctl, Pres_ctl):

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
    TA = 2232
    DIC = 1861
    SiO3 = 1.346
    PO4 = 0.201
    Temp = 26.683
    Sal = 34.004
    Pres = 0
    # call the etamax function
    result = etamax(dTA, TA, DIC, SiO3, PO4, Temp, Sal, Pres)
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
        opt_pH_scale = 1,  # pH scale at which the input pH is reported ("1" means "Total Scale")
        opt_k_carbonic = 10,  # Choice of H2CO3 and HCO3- dissociation constants K1 and K2 ("10" means "Lueker et al 2000")
        opt_k_bisulfate = 1,  # Choice of HSO4- dissociation constant KSO4 ("1" means "Dickson")
        opt_total_borate = 2,  # Choice of boron:sal ("1" means "Lee et al 2010")
        opt_k_fluoride = 2   # "2" means Perez and Fraga 1987
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
        opt_pH_scale = 1,  # pH scale at which the input pH is reported ("1" means "Total Scale")
        opt_k_carbonic = 10,  # Choice of H2CO3 and HCO3- dissociation constants K1 and K2 ("10" means "Lueker et al 2000")
        opt_k_bisulfate = 1,  # Choice of HSO4- dissociation constant KSO4 ("1" means "Dickson")
        opt_total_borate = 2,  # Choice of boron:sal ("1" means "Lee et al 2010")
        opt_k_fluoride = 2   # "2" means Perez and Fraga 1987
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






