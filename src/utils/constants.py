import numpy as np

# Experiment length
DEFAULT_SIMULATION_LENGTH_YEARS = 45  # years
DEFAULT_SPIN_UP_TIME = 5  # years
TOTAL_YEARS = DEFAULT_SIMULATION_LENGTH_YEARS + DEFAULT_SPIN_UP_TIME
TOTAL_DAYS = TOTAL_YEARS * 365

# Constants
DEFAULT_TEMP_CELSIUS = 27.3 #model defalt is 15 #this is WOA annual mean temperature for Gulf of Tehuantepec
DEFAULT_TEMP_KELVIN = DEFAULT_TEMP_CELSIUS + 273.15
DEFAULT_SALINITY = 34.2 #WOA annual mean salinity at Gulf of Tehuantepec is 33.8
# Total borate in mol/kg-sw following Uppstrom (1974)
BORON_COEFFICIENT = 0.0004157 / 35

# Atmospheric Properties
ATMOSPHERIC_CO2 = 400  # ppm
ATMOSPHERIC_d13C = -8  # per mil
ATMOSPHERIC_D14C = 100  # per mil

# Mixed Layer Properties
MIXED_LAYER_DEPTH = 20  # meters
SURFACE_AREA = 1  # m^2 
SWD = 1029  # Seawater density kg/m^3
MIXED_LAYER_INITIAL_D14C = ATMOSPHERIC_D14C # this is NEW from me
MIXED_LAYER_INITIAL_d13C = ATMOSPHERIC_d13C+8.5 # this is NEW from me and assumes an equilibrium isotope fractionation from the atmosphere of +8.5 per mil

# Seasonal Wind Speeds
WIND_SPEED= 5
WIND_SPEED_SUMMER = WIND_SPEED  # m/s
WIND_SPEED_WINTER = WIND_SPEED  # m/s

# Seasonal Vertical Mixing
#REORGANIZING THESE TERMS TO BE JUST ONE "VERTICAL_MIXING"
UPWELLING_RATE_WIND_CONVERSION = 0.1022 #Based on Kris Karnauska's spatial mean estimate of Gulf of Tehuantepec upwelling rate
#UPWELLING_RATE_WIND_CONVERSION = 1.86 #Based on Kris Karnauska's estimate of *peak* Gulf of Tehuantepec upwelling rate
VERTICAL_MIXING = (WIND_SPEED * UPWELLING_RATE_WIND_CONVERSION)/ MIXED_LAYER_DEPTH # in units of per day or (1/day)

VERTICAL_MIXING_WINTER = VERTICAL_MIXING  # 0.1 * the current DIC concentration per day was used in Cai et al., 2020 #supposed to be upwelling
#VERTICAL_MIXING_WINTER = WIND_SPEED_WINTER / 100 # this is an experiment with Dervla to see if it looks more like her code output
VERTICAL_MIXING_SUMMER = VERTICAL_MIXING  # no summer mixing used in Cai et al., 2020
#VERTICAL_MIXING_SUMMER = WIND_SPEED_SUMMER /100 # this is an experiment with Dervla to see if it looks more like her code output
# From the text:    Winter mixing in NWA is simulated by adding 0.1 mmol m−3 DIC
#                   daily from October to February (over 155 days).

# Subsurface boundary conditions
SUBSURFACE_DIC = 1935  # µmol/kg #calculated by Rafter for GLODAP data just upstream of Gulf of Tehuantepec
SUBSURFACE_ALK = 2160  # µmol/kg #calculate by Rafter from existing GLODAP DIC at 20m and Chapa-Balacorta 2015 pCO2 (via CO2SYS)
SUBSURFACE_d13C = 0.0  # per mil - based on correspondance with Wei-Jun Cai
SUBSURFACE_D14C = 0  # per mil - based on GLODAP transect and work to identify the upwelling water mass at 100m

# Seasonal Biology Flux
RAIN_RATIO = 0.07 # PIC : POC. in sinking particles-Sarmiento et al. (2002)

# Biology Parameters
OFFSET_ORG = -22  # isotopic fractionation into organic matter
OFFSET_CC = 2 # isotopic fractionation into CaCO3


# Gas Exchange Parameters
# Calculate K0 (Henry's constant for CO2 solubility in mol/kg-sw/atm) as a function of temperature and salinity
# Weiss, R. F., Marine Chemistry 2:203-215, 1974.
# This is in mol/kg-SW/atm.
def calculate_k0(temp_kelvin, salinity):
    temp_k100 = temp_kelvin / 100
    ln_k0 = (
        -60.2409
        + 93.4517 / temp_k100
        + 23.3585 * np.log(temp_k100)
        + salinity * (0.023517 - 0.023656 * temp_k100 + 0.0047036 * temp_k100**2)
    )
    return np.exp(ln_k0)

# Calculate the Schmidt number (Sc) for CO2
def schmidt(temp_celsius):
    return 2073.1 - 125.62 * temp_celsius + 3.6276 * temp_celsius**2 - 0.043219 * temp_celsius**3

# Calculate kappa based on wind speed, temperature, and salinity
def calculate_kappa(wind_speed, temp_celsius):
    sc = schmidt(temp_celsius)
    return 0.251 * wind_speed**2 * (sc / 660)**-0.5

# Function to calculate piston velocity
def piston_velocity(wind_speed, temp_celsius, salinity):
    return 0.24 * calculate_kappa(wind_speed, temp_celsius) * calculate_k0(temp_celsius + 273.15, salinity)

# Functions to calculate equilibrium constants based on temperature and salinity
def calc_kw(temp_kelvin, salinity):
    return np.exp(
        148.9802
        - 13847.26 / temp_kelvin
        - 23.6521 * np.log(temp_kelvin)
        + (-79.2447 + 3298.72 / temp_kelvin + 12.0408 * np.log(temp_kelvin)) * np.sqrt(salinity)
        - 0.019813 * salinity
    )

def calc_k0(temp_kelvin, salinity):
    """
    Henry's constant for CO2 solubility in mol/kg-sw/atm following W74.
    Weiss, R. F., Marine Chemistry 2:203-215, 1974.
    This is in mol/kg-SW/atm.
    """
    temp_k100 = temp_kelvin / 100
    ln_k0 = (
        -60.2409
        + 93.4517 / temp_k100
        + 23.3585 * np.log(temp_k100)
        + salinity * (0.023517 - 0.023656 * temp_k100 + 0.0047036 * temp_k100**2)
    )
    return np.exp(ln_k0)

def calc_k1(temp_kelvin, salinity):
    """
    Carbonic acid dissociation constants following SLH20
    Coefficients and their 95% confidence intervals from SLH20 Table 1.
    doi:10.5194/os-2020-19
    this is on the Total pH scale in mol/kg-SW
    """
    pk1 = (
        8510.63 / temp_kelvin
        - 172.4493
        + 26.32996 * np.log(temp_kelvin)
        - 0.011555 * salinity
        + 0.0001152 * salinity**2
    )
    return 10.0**-pk1

def calc_k2(temp_kelvin, salinity):
    """
    Carbonic acid dissociation constants following SLH20
    Coefficients and their 95% confidence intervals from SLH20 Table 1.
    doi:10.5194/os-2020-19
    this is on the Total pH scale in mol/kg-SW
    """
    pk2 = (
        4226.23 / temp_kelvin
        - 59.4636
        + 9.60817 * np.log(temp_kelvin)
        - 0.01781 * salinity
        + 0.0001122 * salinity**2
    )
    return 10.0**-pk2

def calc_kb(temp_kelvin, salinity):
    """
    Boric acid dissociation constant following D90b.
    Dickson, A. G., Deep-Sea Research 37:755-766, 1990.
    lnKB is on Total pH scale
    """
    sqr_sal = np.sqrt(salinity)
    ln_kb_top = (
        -8966.9
        - 2890.53 * sqr_sal
        - 77.942 * salinity
        + 1.728 * sqr_sal * salinity
        - 0.0996 * salinity**2
    )
    ln_kb = (
        ln_kb_top / temp_kelvin
        + 148.0248
        + 137.1942 * sqr_sal
        + 1.62142 * salinity
        + (-24.4344 - 25.085 * sqr_sal - 0.2474 * salinity) * np.log(temp_kelvin)
        + 0.053105 * sqr_sal * temp_kelvin
    )
    return np.exp(ln_kb)

def calc_boron(salinity):
    return BORON_COEFFICIENT * salinity
