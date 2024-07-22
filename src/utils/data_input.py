import PyCO2SYS as pyco2
import numpy as np

def load_seasonal_forcings(simulation_length_years):
    """
    Load seasonal forcings for temperature and salinity for a given number of years with daily data.
    This function generates seasonal profiles for the US East Coast.

    Parameters:
    - simulation_length_years (int): Total number of years for the simulation.

    Returns:
    - temperature (array): Temperature profile (°C) for each day of the simulation period.
    - salinity (array): Salinity profile (PSU) for each day of the simulation period.
    """
    days_in_year = 365
    days = np.arange(days_in_year)

    # Generate temperature profile for one year
    temp_mean = 15  # Mean temperature
    temp_amplitude = 8  # Amplitude of temperature variation
    temp_phase_shift = -52  # Phase shift to align with peak temperature in late August (day 240)
    temperature_one_year = temp_mean + temp_amplitude * np.sin(2 * np.pi * (days - temp_phase_shift) / days_in_year)

    # Generate salinity profile for one year
    sal_mean = 33.5  # Mean salinity
    sal_amplitude = 0.5  # Amplitude of salinity variation
    sal_phase_shift = 180 - 31  # Phase shift to align with minimum salinity in July (day 180)
    salinity_one_year = sal_mean + sal_amplitude * np.sin(2 * np.pi * (days - sal_phase_shift) / days_in_year)
    
    # Repeat profiles for the total simulation length
    total_days = simulation_length_years * days_in_year
    temperature = np.tile(temperature_one_year, simulation_length_years)
    salinity = np.tile(salinity_one_year, simulation_length_years)

    return temperature, salinity


def calc_talk(salt, pref):
    """
    Calculate Total Alkalinity (TA) based on salinity and region.

    Parameters:
    - salt (float): Salinity used to calculate TA.
    - pref (str): Region of the US East Coast. Options are:
                  'se' (South Atlantic Bight),
                  'me' (Mid Atlantic Bight),
                  'ne' (North Atlantic outside Gulf of Maine).

    Returns:
    - talk (float): Calculated Total Alkalinity.
    """
    if pref == "se":
        talk = 48.7 * salt + 608.8
    elif pref == "me":
        talk = 46.6 * salt + 670.6
    elif pref == "ne":
        talk = 39.1 * salt + 932.7
    else:
        raise ValueError("Invalid region code. Use 'se', 'me', or 'ne'.")
    return talk

def compute_initial_conditions(pCO2, alk, temp, sal):
    """
    Compute chemical properties using PyCO2SYS.

    Parameters:
    - pCO2 (float): Surface ocean pCO2 concentration.
    - alk (float): Alkalinity concentration.

    Returns:
    - dict: A dictionary containing various chemical properties.
    """
    kwargs = {
        "par1": alk,
        "par2": pCO2,
        "par1_type": 1,  # Alkalinity
        "par2_type": 4,  # pCO2
        "temperature": temp,
        "salinity": sal,
    }
    return pyco2.sys(**kwargs)