# -*- coding: utf-8 -*-

__version__ = "1.0.32"

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
        'opt_k_bisulfate': 1,  # KSO4 1= Dickson 1990, 2= Khoo et al 1977, 3= Waters and Millero 2013/Waters et al 2014
        'opt_total_borate': 1,  # boron:salinity 1= Uppstrom 1974, 2= Lee et al 2010, 3= Kulinski et al 2018
        'opt_k_fluoride': 1   # HF dissociation 1= Dickson and Riley 1979, 2= Perez and Fraga 1987        
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
        - 'DIC_ctl': DIC at control conditions (umol/kg)
        - 'TA_ctl': TA at control conditions (umol/kg)
        - 'CO3_ctl': Carbonate ion concentration at control conditions (umol/kg)
        - 'pHtotal_ctl': pH (total) at control conditions
        - 'OmegaAra_ctl': Aragonite saturation state at control conditions
        - 'OmegaCal_ctl': Calcite saturation state at control conditions
        - 'RevelleFactor_ctl': Revelle Factor at control conditions

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
        'opt_k_bisulfate': 1,  # KSO4 1= Dickson 1990, 2= Khoo et al 1977, 3= Waters and Millero 2013/Waters et al 2014
        'opt_total_borate': 1,  # boron:salinity 1= Uppstrom 1974, 2= Lee et al 2010, 3= Kulinski et al 2018
        'opt_k_fluoride': 1   # HF dissociation 1= Dickson and Riley 1979, 2= Perez and Fraga 1987        
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
    TA_ctl = results['alkalinity']
    DIC_ctl = results['dic']
    CO3_ctl = results['CO3']
    pHtotal_ctl = results['pH_total']
    OmegaAra_ctl = results['saturation_aragonite']
    OmegaCal_ctl = results['saturation_calcite']
    RevelleFactor_ctl = results['revelle_factor']

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
            'pCO2_ctl': pCO2_ctl,
            'DIC_ctl': DIC_ctl,
            'TA_ctl': TA_ctl,
            'CO3_ctl': CO3_ctl,
            'pHtotal_ctl': pHtotal_ctl,
            'OmegaAra_ctl': OmegaAra_ctl,
            'OmegaCal_ctl': OmegaCal_ctl,
            'RevelleFactor_ctl': RevelleFactor_ctl
            }
    
    return result

def dic_bio(**kwargs):

    """
    PURPOSE

    Calculate the biological component of DIC (DIC_bio) 
    as the difference between the observed DIC (DIC_obs) 
    compared with the DIC at equilibrium with atmospheric pCO2 (DIC_atm)
    using the method of Clargo et al (2015) and Burt et al (2016):

    - Clargo et al 2015 (https://doi.org/10.1016/j.marchem.2015.08.010)
    - Burt et al 2016 (https://doi.org/10.1002/lno.10243):

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
        'opt_k_bisulfate': 1,  # KSO4 1= Dickson 1990, 2= Khoo et al 1977, 3= Waters and Millero 2013/Waters et al 2014
        'opt_total_borate': 1,  # boron:salinity 1= Uppstrom 1974, 2= Lee et al 2010, 3= Kulinski et al 2018
        'opt_k_fluoride': 1   # HF dissociation 1= Dickson and Riley 1979, 2= Perez and Fraga 1987        
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

def pco2_tnorm(pCO2_obs,temp_obs):

    """
    PURPOSE

    Temperature normalization of pCO2 using the equation of 
    Takahashi (2002) (http://dx.doi.org/10.1016/S0967-0645(02)00003-6)

    In order to remove the temperature effect from the observed pCO2,
    the observed pCO2 values are normalized to the long-term mean
    temperature of seawater at the location. The following equation
    is used to estimate the non-thermal component of pCO2 (pCO2_NT):

    pCO2_NT = pCO2_obs * exp(0.0423 * (T_mean - T_obs))     (eqn 1)

    In a similar manner, the effect of temperature changes on pCO2 is computed
    by perturbing the mean annual pCO2. The thermal component of pCO2 (pCO2_T)
    is computed using the following equation:

    pCO2_T = pCO2_mean * exp(0.0423 * (T_obs - T_mean))     (eqn 2)

    INPUT

    pCO2_obs = observed pCO2 (e.g. monthly average pCO2 from OceanSODA-ETHZ) (uatm)
    temp_obs = observed temperature co-inciding with observed pCO2 (degC)

    OUTPUT

    result = dictionary of outputs including the following:
    'pCO2_NT': non-thermal component of pCO2 (uatm)
    'pCO2_T': thermal component of pCO2 (uatm)
    'pCO2_mean': long-term mean observed pCO2 (uatm)
    'pCO2_Tanom': thermal component seasonal cycle anomaly pCO2 (uatm)
        
    """

    import numpy as np
    import sys
    
    x = pCO2_obs
    y = temp_obs
    ctrl = np.isreal(x).all() and (not np.isnan(x).any()) and (not np.isinf(x).any()) and x.ndim==1
    if not ctrl:
      print('Check pCO2 input: it needs be a vector of real numbers with no infinite or nan values!','\n')
      sys.exit()
    ctrl = np.isreal(y).all() and (not np.isnan(y).any()) and (not np.isinf(y).any()) and y.ndim==1
    if not ctrl:
      print('Check temperature input: it needs be a vector of real numbers with no infinite or nan values!','\n')
      sys.exit()
    ctrl =  np.size(x)==np.size(y)
    if not ctrl:
      print('Check pCO2 and temperature: both need to be the same size!','\n')
      sys.exit()

    # - - -
    # pCO2_NT = non-thermal component pCO2 (uatm)
    temp_mean = np.nanmean(temp_obs)
    pCO2_NT = pCO2_obs * np.exp(0.0423 * (temp_mean - temp_obs))

    # - - -
    # pCO2_T = thermal component pCO2 (uatm)
    pCO2_mean = np.nanmean(pCO2_obs)
    pCO2_T = pCO2_mean * np.exp(0.0423 * (temp_obs - temp_mean))

    # - - -
    # pCO2_Tanom = thermal pCO2 component  seasonal cycle anomaly (uatm)
    pCO2_Tanom = pCO2_T - pCO2_mean

    result = {
        'pCO2_NT': pCO2_NT,
        'pCO2_T': pCO2_T,
        'pCO2_Tanom': pCO2_Tanom,
        'pCO2_mean': pCO2_mean
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
    described by eqn 4-7 in Rodgers et al 2022 (https://doi.org/10.1029/2023GB007798)

    INPUT 

    kwargs dictionary with the following keyword arguments:
    'alkalinity' = observed monthly alkalinity (umol/kg)
    'dic' = observed monthly DIC (umol/kg)
    'pCO2' = observed monthly pCO2 (uatm) (optional)
    'fCO2' = observed monthly fCO2 (uatm) (optional)
    'par2_type' = PyCO2SYS second parameter type, 2=DIC, 4=pCO2, 5=fCO2 (optional)
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
    'pCO2_obs' = observed monthly pCO2 (uatm)
    'pCO2_mean' = mean pCO2 estimated using 
        mean TA, DIC, SIO4, PO4, T, S (uatm)
    'pCO2_T' = monthly thermally-driven pCO2 using 
        monthly T with mean TA, DIC, SIO4, PO4, S (uatm)
    'pCO2_Tanom' = thermal pCO2 component seasonal cycle anomaly (uatm)
    'pCO2_NT' = monthly non-thermal component of pCO2 (uatm)
    """

    import numpy as np
    import PyCO2SYS as pyco2

    # Define default values of input data arguments
    defaults = {
        'alkalinity': np.array([]),
        'dic': np.array([]),
        'pCO2': np.array([]),
        'fCO2': np.array([]),
        'par2_type': 2,   # The second parameter 2=DIC, 4=pCO2, 5=fCO2
        'total_silicate': np.array([]),
        'total_phosphate': np.array([]),
        'temperature': np.array([]),
        'salinity': np.array([]),
        'pressure': np.array([]),
        'opt_pH_scale': 1,  # pH scale at which the input pH is reported ("1" means "Total Scale")
        'opt_k_carbonic': 10,  # Choice of H2CO3 and HCO3- dissociation constants K1 and K2 ("10" means "Lueker et al 2000")
        'opt_k_bisulfate': 1,  # Choice of HSO4- dissociation constant KSO4 ("1" means "Dickson")
        'opt_total_borate': 1,  # Choice of boron:sal ("1" means "Lee et al 2010")
        'opt_k_fluoride': 1   # "2" means Perez and Fraga 1987        
        }

    # Update input data argumements with any provided keyword arguments in kwargs
    data = {**defaults, **kwargs}

    if data["par2_type"]==2:
        data["par2"] = data["dic"]
    if data["par2_type"]==4:
        data["par2"] = data["pCO2"]
        pCO2_obs = data['pCO2']
    if data["par2_type"]==5:
        data["par2"] = data["fCO2"]

    # - - -
    # Step 1: solve for observed monthly pCO2 from observed monthly TA, DIC (or fCO2), SIO4, PO4, T, S
    if data["par2_type"]!=4:
        kwargs = dict(
            par1_type = 1,  # The first parameter supplied is of type "1", which means "alkalinity"
            par1 = data["alkalinity"],  # value of the first parameter
            par2_type = data["par2_type"],  # The second parameter 2=DIC, 4=pCO2, 5=fCO2
            par2 = data["par2"],  # value of the second parameter
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
        pCO2_obs = results['pCO2']

    # - - -
    # Step 2: solve for long-term average pCO2_mean (eqn 4 in Rodgers et al)
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
    pCO2_mean = results['pCO2']
 
    # - - -
    # Step 3: solve for thermal component pCO2_T (eqn 5 in Rodgers et al)
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
    pCO2_T = results['pCO2']

    # - - -
    # Step 4: thermal pCO2 component  seasonal cycle anomaly (pCO2_Tanom) (eqn 6 in Rodgers et al)
    pCO2_Tanom = pCO2_T - pCO2_mean

    # - - -
    # Step 5: solve for non-thermal component pCO2_NT (eqn 7 in Rodgers et al)
    pCO2_NT = pCO2_obs - pCO2_Tanom
 
    # make the result dictionary for output
    result = {
            'pCO2_obs': pCO2_obs,
            'pCO2_mean': pCO2_mean,
            'pCO2_T': pCO2_T,
            'pCO2_Tanom': pCO2_Tanom,
            'pCO2_NT': pCO2_NT
            }

    return result

def m_per_deg_lat(latitude):

    """
    Calculate the length per degree of latitude in meters
    Example usage:
    latitude = 47.0  # Example latitude
    lat_length = m_per_deg_lat(latitude)
    Reference:
    See https://gis.stackexchange.com/questions/75528/understanding-terms-in-length-of-degree-formula
    for using Taylor series for degree-to-meter on ellipsoid with m1-p3.
    """
    import math
    # Convert latitude to radians
    lat_rad = math.radians(latitude)
    # Calculate the length of a degree of latitude
    m1 = 111132.92  # latitude calculation term 1
    m2 = -559.82    # latitude calculation term 2
    m3 = 1.175      # latitude calculation term 3
    m4 = -0.0023    # latitude calculation term 4
    lat_length = m1 + (m2 * math.cos(2 * lat_rad)) + (m3 * math.cos(4 * lat_rad)) + (m4 * math.cos(6 * lat_rad))
    return lat_length

def m_per_deg_lon(latitude):

    """
    Calculate the length per degree of longitude in meters
    Example usage:
    latitude = 47.0  # Example latitude
    lon_length = m_per_deg_lon(latitude)
    Reference:
    See https://gis.stackexchange.com/questions/75528/understanding-terms-in-length-of-degree-formula
    for using Taylor series for degree-to-meter on ellipsoid with m1-p3.
    """
    import math
    # Convert latitude to radians
    lat_rad = math.radians(latitude)
    # Calculate the length of a degree of longitude
    p1 = 111412.84  # longitude calculation term 1
    p2 = -93.5      # longitude calculation term 2
    p3 = 0.118      # longitude calculation term 3
    lon_length = (p1 * math.cos(lat_rad)) + (p2 * math.cos(3 * lat_rad)) + (p3 * math.cos(5 * lat_rad))
    return lon_length

def haversine_distance(lat1, lon1, lat2, lon2):

    """
    Haversine distance in kilometers between any two points 
    defined by longitude and latitude coordinates
    used by centered_grid_cell_area
    """
        
    import math
    import numpy as np

    # Radius of the Earth in kilometers
    # R = 6371.0
    R = 6378.137
    
    # Convert degrees to radians
    # lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])
    lat1 = np.radians(lat1)
    lon1 = np.radians(lon1)
    lat2 = np.radians(lat2)
    lon2 = np.radians(lon2)
    
    # Differences in coordinates
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    
    # Haversine formula
    a = math.sin(dlat / 2)**2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2)**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    
    # Distance in kilometers
    distance = R * c
    return distance

def centered_grid_cell_area(lat, lon, dlat, dlon):

    """
    Grid cell area in km^2 given the central latitude and longitude degrees (lat, lon)
    and the cell height and width in degrees (dlat, dlon)
    using haversine_distance
    EXAMPLE USAGE
    lat = 47.5  # Example latitude
    lon = -122.5  # Example longitude
    dlat = 1.0  # 1 degree latitude
    dlon = 1.0  # 1 degree longitude
    area_km2 = centered_grid_cell_area(lat, lon, dlat, dlon)
    """
        
    import math

    # Calculate the four corners of the grid cell around central lat,lon
    lat1, lon1 = lat - dlat/2, lon - dlon/2
    lat2, lon2 = lat + dlat/2, lon - dlon/2
    lat3, lon3 = lat - dlat/2, lon + dlon/2
    lat4, lon4 = lat + dlat/2, lon + dlon/2
    
    # Calculate the lengths of the sides using the Haversine formula
    side1 = haversine_distance(lat1, lon1, lat2, lon2)
    side2 = haversine_distance(lat1, lon1, lat3, lon3)
    side3 = haversine_distance(lat2, lon2, lat4, lon4)
    side4 = haversine_distance(lat3, lon3, lat4, lon4)
    
    # Approximate the area of the grid cell (assuming it's roughly rectangular)
    area = ((side1 + side3) / 2) * ((side2 + side4) / 2)
    return area

def fCO2_to_pCO2(fCO2, TempC, pressure_atmosphere=1.0, RGas=83.1451):
    """
    Convert CO2 fugacity to partial pressure.
    adapted from PyCO2SYS
    """

    import numpy as np

    # RGas = 83.1451;  # % ml bar-1 K-1 mol-1, DOEv2
    TempK    = TempC + 273.15
    RT = RGas * TempK
    Delta = 57.7 - 0.118 * TempK
    b = (
        -1636.75
        + 12.0408 * TempK
        - 0.0327957 * TempK**2
        + 3.16528 * 0.00001 * TempK**3
    )
    # # For a mixture of CO2 and air at 1 atm (at low CO2 concentrations):
    # P1atm = 1.01325  # in bar
    p_bar = pressure_atmosphere * 1.01325  # units conversion
    FugFac = np.exp((b + 2 * Delta) * p_bar / RT)

    #% Generate the associated pCO2 from fCO2:
    pCO2 = fCO2 / FugFac

    return pCO2

def pCO2_to_fCO2(pCO2, TempC, pressure_atmosphere=1.0, RGas=83.1451):
    """
    Convert CO2 partial pressure to fugacity.
    adapted from PyCO2SYS
    """

    import numpy as np

    # RGas = 83.1451;  # % ml bar-1 K-1 mol-1, DOEv2
    TempK    = TempC + 273.15
    RT = RGas * TempK
    Delta = 57.7 - 0.118 * TempK
    b = (
        -1636.75
        + 12.0408 * TempK
        - 0.0327957 * TempK**2
        + 3.16528 * 0.00001 * TempK**3
    )
    # # For a mixture of CO2 and air at 1 atm (at low CO2 concentrations):
    # P1atm = 1.01325  # in bar
    p_bar = pressure_atmosphere * 1.01325  # units conversion
    FugFac = np.exp((b + 2 * Delta) * p_bar / RT)

    #% Generate the associated fCO2 from pCO2:
    fCO2 = pCO2 * FugFac

    return pCO2

# def fugacityfactor(TempC, WhichKs, RGas, pressure_atmosphere=1.0):
def fugacityfactor(TempC, pressure_atmosphere=1.0, RGas=83.1451):
    """
    Calculate the fugacity factor.
    adapted from PyCO2SYS
    """

    import numpy as np

    # This assumes that the pressure is at one atmosphere, or close to it.
    # Otherwise, the Pres term in the exponent affects the results.
    # Following Weiss, R. F., Marine Chemistry 2:203-215, 1974.
    # Delta and B are in cm**3/mol.
    # RGas = 83.1451;  # % ml bar-1 K-1 mol-1, DOEv2
    TempK    = TempC + 273.15
    RT = RGas * TempK
    Delta = 57.7 - 0.118 * TempK
    b = (
        -1636.75
        + 12.0408 * TempK
        - 0.0327957 * TempK**2
        + 3.16528 * 0.00001 * TempK**3
    )
    # # For a mixture of CO2 and air at 1 atm (at low CO2 concentrations):
    # P1atm = 1.01325  # in bar
    p_bar = pressure_atmosphere * 1.01325  # units conversion
    FugFac = np.exp((b + 2 * Delta) * p_bar / RT)
    # GEOSECS and Peng assume pCO2 = fCO2, or FugFac = 1
    # FugFac = np.where((WhichKs == 6) | (WhichKs == 7), 1.0, FugFac)
    return FugFac

def moving_average(x, w=12):
    """
    moving mean of vector x with window size w
    with default window size w=12 
    example:
    a = np.arange(20)
    a_n4 = moving_average(a, 4)
    print('a_n4= ',a_n4, a_n4.shape)
    # a_n4=  [ 1.5  2.5  3.5  4.5  5.5  6.5  7.5  8.5  9.5 10.5 11.5 12.5 13.5 14.5
    # 15.5 16.5 17.5] (17,)    
    """
    import numpy as np
    import sys
    ctrl = np.isreal(x).all() and (not np.isnan(x).any()) and (not np.isinf(x).any()) and x.ndim==1
    if not ctrl:
      print('Check x: it needs be a vector of real numbers with no infinite or nan values!','\n')
      sys.exit()
    return np.convolve(x, np.ones(w), 'valid') / w





