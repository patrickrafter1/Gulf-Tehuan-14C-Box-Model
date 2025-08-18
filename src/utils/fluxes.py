import numpy as np
import PyCO2SYS as pyco2
from . import constants

def alpha_to_epsilon(alpha):
    """Converts alpha value to epsilon."""
    return (alpha - 1) * 1000

def mol_to_umolkg_conversion(surface_mass):
    """Converts from units of moles to units of concentration."""
    return 1e6 / surface_mass

def calculate_equilibrium_constants(temp, salt):
    """
    Calculate equilibrium constants based on temperature and salinity.

    Parameters:
        temp (float): Temperature in degrees Celsius.
        salt (float): Salinity in PSU.

    Returns:
        dict: A dictionary of equilibrium constants.
    """
    # local variable definitions
    Tk = temp + 273.15
    sqrtS = np.sqrt(salt)
    #Tk = temp + 273.15
    centiTk = 0.01 * Tk
    invTk = 1.0 / Tk
    logTk = np.log(Tk)
    sqrtS = np.sqrt(salt)
    SO4 = 19.924 * salt / (1000 - 1.005 * salt)
    sqrtSO4 = np.sqrt(SO4)
    scl = salt / 1.80655

    K1 = 10.0 ** (62.008 - 3670.7 / Tk - 9.7944 * np.log(Tk) + salt * (0.0118 - salt * 0.000116))
    K2 = 10.0 ** (-4.777 - 1394.7 / Tk + salt * (0.0184 - salt * 0.000118))
    
    # warnings.filterwarnings("ignore", category=RuntimeWarning, message="divide by zero encountered in scalar divide")

    Kb = np.exp(
        -invTk
        * (
            8966.90
            + sqrtS * (2890.53 + sqrtS * (77.942 - sqrtS * (1.728 - sqrtS * 0.0996)))
        )
        - logTk * (24.4344 + sqrtS * (25.085 + sqrtS * 0.2474))
        + Tk * (sqrtS * 0.053105)
        + 148.0248
        + sqrtS * (137.1942 + sqrtS * 1.62142)
    )

    K1p = np.exp(
        115.525 - 4576.752 / Tk - 18.453 * np.log(Tk) 
        + sqrtS * (0.69171 - 106.736 / Tk) 
        - salt * (0.01844 + 0.65643 / Tk)
    )
    K2p = np.exp(
        172.0883 - 8814.715 / Tk - 27.927 * np.log(Tk) 
        + sqrtS * (1.3566 - 160.340 / Tk) 
        - salt * (0.05778 - 0.37335 / Tk)
    )
    K3p = np.exp(
        -18.141 - 3070.75 / Tk + sqrtS * (2.81197 + 17.27039 / Tk) 
        - salt * (0.09984 + 44.99486 / Tk)
    )

    Ksi = np.exp(
        117.385
        - invTk * 8904.2
        - logTk * 19.334
        + sqrtSO4 * (3.5913 - invTk * 458.79)
        - SO4 * (1.5998 - invTk * 188.74 - SO4 * (0.07871 - invTk * 12.1652))
        + np.log(1.0 - 0.001005 * salt)
    )

    Kw = np.exp(
        148.9652 - 13847.26 / Tk - 23.6521 * np.log(Tk) 
        - sqrtS * (5.977 - 118.67 / Tk - 1.0495 * np.log(Tk)) 
        - salt * 0.01615
    )

    Ks = np.exp(
        141.328
        - invTk * 4276.1
        - logTk * 23.093
        +  sqrtSO4* (324.57 - invTk * 13856.0 - logTk * 47.986 - SO4 * invTk * 2698.0)
        - SO4 * (771.54 - invTk * 35474.0 - logTk * 114.723 - SO4 * invTk * 1776.0)
        + np.log(1.0 - 0.001005 * salt)
    )

    Kf = np.exp(
        -12.641
        + invTk * 1590.2
        + sqrtSO4 * 1.525
        + np.log(1.0 - 0.001005 * salt)
        + np.log(1.0 + 0.1400 * scl / (96.062 * Ks))
    )

    constants = {
        "K1": K1,
        "K2": K2,
        "Kb": Kb,
        "K1p": K1p,
        "K2p": K2p,
        "K3p": K3p,
        "Ksi": Ksi,
        "Kw": Kw,
        "Ks": Ks,
        "Kf": Kf,
    }

    return constants

def h_total(
    H, K1, K1p, K12, K12p, K123p, Kf, Ks, Kw, invKb, invKs, invKsi, alk, borate, dic, fluoride, phos, sili, sulfate
):
    """
    Evaluate the root, f([H+])
    Eqn (18), the expression for AT, from Dickson 2007, in terms of total concentrations and [H+].
    """
    H3 = H * H * H
    invH = 1.0 / H
    A = H * (K12p + H * (K1p + H)) + K123p
    B = H * (K1 + H) + K12
    C = 1.0 / (1.0 + sulfate * invKs)

    res = (
        dic * (K1 * H + 2.0 * K12) / B
        + borate / (1.0 + H * invKb)
        + Kw * invH
        + phos * (K12p * H + 2.0 * K123p - H3) / A
        + sili / (1.0 + H * invKsi)
        - H * C
        - sulfate / (1.0 + Ks * invH * C)
        - fluoride / (1.0 + Kf * invH)
        - alk
    )
    return res

def calculate_pco2_iteratively(tic, talk, temp, salt, po4=0.0, sio3=0.0):
    """
    This routine computes equilibrium partial pressure of CO2 (pCO2) in the surface seawater.
    """
    # local variable definitions
    Tk = temp + 273.15
    centiTk = 0.01 * Tk
    SO4 = 19.924 * salt / (1000 - 1.005 * salt)
    scl = salt / 1.80655

    # converting from umol/kg or mmol/m3 to mol/kg
    alk = talk * 0.000001
    dic = tic * 0.000001
    phos = po4 * 0.000001
    sili = sio3 * 0.000001

    # Equilibrium constants
    eq_constants = calculate_equilibrium_constants(temp, salt)
    K1 = eq_constants["K1"]
    K2 = eq_constants["K2"]
    Kb = eq_constants["Kb"]
    K1p = eq_constants["K1p"]
    K2p = eq_constants["K2p"]
    K3p = eq_constants["K3p"]
    Ksi = eq_constants["Ksi"]
    Kw = eq_constants["Kw"]
    Ks = eq_constants["Ks"]
    Kf = eq_constants["Kf"]

    ff = np.exp(
        -162.8301
        + 218.2968 / centiTk
        + np.log(centiTk) * 90.9241
        - centiTk**2 * 1.47696
        + salt * (0.025695 - centiTk * (0.025225 - centiTk * 0.0049867))
    )

    borate = 0.000232 * scl / 10.811
    sulfate = 0.14 * scl / 96.062
    fluoride = 0.000067 * scl / 18.9984

    K12 = K1 * K2
    K12p = K1p * K2p
    K123p = K12p * K3p
    invKb = 1.0 / Kb

    invKs = 1.0 / Ks
    invKsi = 1.0 / Ksi

    IbrackMax = 30
    H_lo = 10.0 ** (-10.0)
    H_hi = 10.0 ** (-5.0)

    f_hi = h_total(
        H_hi, K1, K1p, K12, K12p, K123p, Kf, Ks, Kw, invKb, invKs, invKsi, alk, borate, dic, fluoride, phos, sili, sulfate
    )
    H_mid = 0.5 * (H_lo + H_hi)

    for Ibrack in range(1, IbrackMax + 1):
        H_mid = 0.5 * (H_lo + H_hi)
        f_mid = h_total(
            H_mid, K1, K1p, K12, K12p, K123p, Kf, Ks, Kw, invKb, invKs, invKsi, alk, borate, dic, fluoride, phos, sili, sulfate
        )

        if f_mid == 0:
            break
        else:
            ftest = f_hi / f_mid
            if ftest > 0:
                H_hi = H_mid
                f_hi = f_mid
            else:
                H_lo = H_mid
            H_mid = 0.5 * (H_lo + H_hi)

    H = H_mid

    H2 = H * H
    co2starik = dic * H2 / (H2 + K1 * H + K1 * K2)
    co3_insitu = dic * K1 * K2 / (H2 + K1 * H + K1 * K2)
    co3_insitu = co3_insitu

    ca_sat = 10.28 / 1e3
    ksp = 10.0 ** (
        -171.945
        - 0.077993 * Tk
        + 2903.293 / Tk
        + 71.595 * np.log10(Tk)
        + (-0.068393 + 0.0017276 * Tk + 88.135 / Tk) * salt**0.5
        - 0.10018 * salt
        + 0.0059415 * salt**1.5
    )
    co3_sat = ksp / ca_sat
    omega_aragik = co3_insitu / co3_sat
    phik = -np.log10(H)
    pco2 = co2starik * 1000000.0 / ff

    return pco2, phik, omega_aragik

def calculate_pco2_func(current_state, temp_celsius, salinity, use_pyco2sys=False):
    """
    Calculate pCO2 of the surface ocean.

    Parameters:
        current_state (array-like): Current state with DIC and TA concentrations.
        temp_celsius (float): Temperature in Celsius.
        salinity (float): Salinity.
        use_pyco2sys (bool, optional): If True, use PyCO2SYS for calculation (slower but accurate).
                                       If False, use manual calculation (faster but approximate).
                                       Defaults to True.

    Returns:
        float: pCO2 value in ppm.
    """
    dic_umolkg = current_state[0]
    alk_umolkg = current_state[1]

    if use_pyco2sys:
        carbon_chemistry = pyco2.sys(
            par1=alk_umolkg,
            par2=dic_umolkg,
            par1_type=1,
            par2_type=2,
            temperature=temp_celsius,
            salinity=salinity,
        )
        pco2_ocean = carbon_chemistry["pCO2"]
        pH = carbon_chemistry["pH"]
        omega = carbon_chemistry["saturation_aragonite"]
    else:
        pco2_ocean, pH, omega = calculate_pco2_iteratively(dic_umolkg, alk_umolkg, temp_celsius, salinity)

    return pco2_ocean

def gas_exchange(
    current_state, num_tracers, CO2_atm, d13C_atm, D14C_atm, surface_area, surface_mass, temp_celsius, salinity, day_of_year
):
    """
    Calculate gas exchange fluxes.

    Parameters:
        current_state (array-like): Current state with DIC, TA, d13C, and D14C concentrations.
        num_tracer (int): Number of tracers.
        CO2_atm (float): Atmospheric CO2 concentration.
        d13C_atm (float): Atmospheric d13C concentration.
        D14C_atm (float): Atmospheric D14C concentration.
        surface_area (float): Surface area of the mixed layer.
        surface_mass (float): Surface mass of the mixed layer.
        temp_celsius (float): Temperature in Celsius.
        salinity (float): Salinity.
        day_of_year (int): Day of the year (1-365).


    Returns:
        tuple: Rate of change of tracers.
    """

    # Determine wind speed based on day of the year
    if day_of_year < 79 or day_of_year > 273:  # Winter period (October through February)
        wind_speed = constants.WIND_SPEED_WINTER
    else:  # Summer period (March through September)
        wind_speed = constants.WIND_SPEED_SUMMER
    
    # Calculate pCO2 of the ocean using the previously defined function
    pco2_ocean = calculate_pco2_func(current_state, temp_celsius, salinity,use_pyco2sys=False)

    # Other variables
    d13C_ocean = current_state[2] / current_state[0]  # ppmil
    D14C_ocean = current_state[3] / current_state[0]
    pco2_atm = CO2_atm
    d13C_atm = d13C_atm
    D14C_atm = D14C_atm

    # Kinetic Fractionation factor for CO2 gas transfer across the air-sea interface
    kinetic_frac_c13 = 0.9995  #  0.99919 in Zhang 1995
    kinetic_frac_c13_permil = alpha_to_epsilon(kinetic_frac_c13)  # converted to per mil

    # Temperature dependent equilibrium fractionation factors
    eq_frac_sa_c13 = -9.866 / (temp_celsius + 273.15) + 1.02412
    eq_frac_sa_c13_permil = alpha_to_epsilon(eq_frac_sa_c13)

    eq_frac_as_c13 = -0.373 / (temp_celsius + 273.15) + 1.00019
    eq_frac_as_c13_permil = alpha_to_epsilon(eq_frac_as_c13)

    # Calculate piston velocity based on wind speed and temperature
    piston_velocity = constants.piston_velocity(wind_speed, temp_celsius, salinity)
    piston_velocity_year = piston_velocity  # mol/(m²·yr·atm)
    
    SeatoAir = piston_velocity_year * (pco2_ocean/1e6) * surface_area  # mol yr-1
    AirtoSea = piston_velocity_year * (pco2_atm/1e6) * surface_area  # mol yr-1

    fract_sa_c13 = eq_frac_sa_c13_permil + kinetic_frac_c13_permil
    fract_as_c13 = eq_frac_as_c13_permil + kinetic_frac_c13_permil

    dC13_SeatoAir = SeatoAir * (d13C_ocean + fract_sa_c13)
    dC13_AirtoSea = AirtoSea * (d13C_atm + fract_as_c13)
    DC14_SeatoAir = SeatoAir * (D14C_ocean + (2 * fract_sa_c13))
    DC14_AirtoSea = AirtoSea * (D14C_atm + (2 * fract_as_c13))

    d_dt = np.zeros((num_tracers))
    d_dt[0] += (AirtoSea - SeatoAir) * mol_to_umolkg_conversion(surface_mass)
    d_dt[2] += (dC13_AirtoSea - dC13_SeatoAir) * mol_to_umolkg_conversion(surface_mass)
    d_dt[3] += (DC14_AirtoSea - DC14_SeatoAir) * mol_to_umolkg_conversion(surface_mass)

    return d_dt

def vertical_mixing(current_state, num_tracers, day_of_year):
    """
    Calculate mixing fluxes. 
    Values based on Cai, WJ., Xu, YY., Feely, R.A. et al. Controls on surface water carbonate chemistry along North American ocean margins. Nat Commun 11, 2691 (2020). https://doi.org/10.1038/s41467-020-16530-z

    Parameters:
        current_state (array-like): Current state with DIC, TA, d13C, and D14C, NO3 concentrations.
        num_tracer (int): Number of tracers.
        surface_area (float): Surface area of the mixed layer.
        surface_mass (float): Surface mass of the mixed layer.
        day_of_year (int): Day of the year (1-365).

    Returns:
        tuple: Rate of change of tracers.
    """
    # Current DIC concentration
    current_dic = current_state[0]

    # Determine the mixing rate based on day of the year
    if day_of_year < 79 or day_of_year > 273:  # Winter period (October through February)
        fractional_mixing = constants.VERTICAL_MIXING_WINTER / current_dic # per day
    else:  # Summer period (March through September)
        fractional_mixing = constants.VERTICAL_MIXING_SUMMER / current_dic # per day

    DIC_vertical_mixing_in = fractional_mixing * constants.SUBSURFACE_DIC
    DIC_vertical_mixing_out = fractional_mixing * current_state[0]
    
    ALK_vertical_mixing_in = fractional_mixing * constants.SUBSURFACE_ALK
    ALK_vertical_mixing_out = fractional_mixing * current_state[1]

    d13C_DIC_vertical_mixing_in = DIC_vertical_mixing_in * constants.SUBSURFACE_d13C
    d13C_DIC_vertical_mixing_out = fractional_mixing * current_state[2] # current state already in delta * concentration units

    D14C_DIC_vertical_mixing_in = DIC_vertical_mixing_in * constants.SUBSURFACE_D14C
    D14C_DIC_vertical_mixing_out = fractional_mixing * current_state[3] # current state already in delta * concentration units

    # Calculate the mixing fluxes
    d_dt = np.zeros((num_tracers))
    d_dt[0] += DIC_vertical_mixing_in - DIC_vertical_mixing_out
    d_dt[1] += ALK_vertical_mixing_in - ALK_vertical_mixing_out
    d_dt[2] += d13C_DIC_vertical_mixing_in - d13C_DIC_vertical_mixing_out
    d_dt[3] += D14C_DIC_vertical_mixing_in - D14C_DIC_vertical_mixing_out

    return d_dt

def biology(current_state, num_tracers, ncp):
    """
    Calculate biological fluxes, including export production and CaCO₃ dynamics.

    Parameters:
        current_state (array-like): Current state with DIC, TA, d13C, D14C, and NO3 concentrations.
        num_tracers (int): Number of tracers.
        day_of_year (int): Day of the year (1-365).
        ncp (float): Net community production based on Cai et al., 2020.

    Returns:
        array: Rate of change of tracers.
    """

    # Export production fraction (portion of NCP exported as POC)
    poc_export = ncp
    pic_export = constants.RAIN_RATIO * poc_export  # PIC export

    # Calculate d13C and D14C values for POC and PIC
    current_dic = current_state[0]
    current_d13C_ocean = current_state[2] / current_dic
    d13C_poc = current_d13C_ocean + constants.OFFSET_ORG
    d13C_pic = current_d13C_ocean + constants.OFFSET_CC
    current_D14C_ocean = current_state[3] / current_dic
    D14C_poc = current_D14C_ocean + 2 * constants.OFFSET_ORG
    D14C_pic = current_D14C_ocean + 2 * constants.OFFSET_CC

    d_dt = np.zeros((num_tracers))

    # Update DIC and tracers due to POC remineralization and PIC dissolution
    d_dt[0] -= (poc_export + pic_export)
    d_dt[1] -= (pic_export * 2)
    d_dt[2] -= (poc_export * d13C_poc) + (pic_export * d13C_pic)
    d_dt[3] -= (poc_export * D14C_poc) + (pic_export * D14C_pic)

    return d_dt

def dilution(current_state, num_tracers, salinity_forcing):
    """
    We use salinity to calculate changes in tracer concentrations due to total volume changes.

    For simplicity, this function assumes that salinity changes (e.g., due to evaporation and 
    precipitation) do not alter the isotopic composition (δ13C, δ14C), and therefore no changes 
    in isotopic ratios are simulated during this process.

    Parameters:
    - current_state (array): Current state variables [DIC, ALK, d13C, D14C, Salinity]
    - num_tracers (int): Number of tracers
    - salinity_forcing (float): New salinity value (PSU)

    Returns:
    - d_dt (array): Rate of change of tracers due to salinity changes
    """
    # Current and new salinity
    current_salinity = current_state[4]
    new_salinity = salinity_forcing
    
    # Dilution factor
    dilution_factor = new_salinity / current_salinity
    
    # Changes in concentration due to salinity changes
    d_dt = np.zeros((num_tracers))
    for i in range(num_tracers):
        d_dt[i] = (dilution_factor * current_state[i]) - current_state[i]
    
    return d_dt