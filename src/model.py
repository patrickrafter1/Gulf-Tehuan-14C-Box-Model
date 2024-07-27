"""
Carbon isotope box model for an idealized surface mixed layer.
Inspired by Lynch-Stieglitz et al. (1995) and Cei et al. (2020).
"""

import numpy as np
from scipy.integrate import solve_ivp
import time

from utils import constants
from utils import fluxes
from utils import data_input
from utils import data_output

class MixedLayerModel:
    """
    Single box model for a mixed layer of the ocean.
    The model includes DIC, ALK, d13C, and D14C tracers.

    Atmospheric conditions are set to:
    - CO2: 395 ppm
    - d13C: -8 per mil
    - D14C: 0 per mil

    Oceanic conditions are set to:
    - pCO2: 395 ppm (in equilibrium with atmospheric CO2)
    - Alkalinity: calculated based on salinty using calc_talk function
    - DIC: calculated with PyCO2SYS from pCO2 and alkalinity
    - d13C: 1 per mil
    - D14C: 0 per mil

    Parameters:
    - temp_celsius (array-like): Time-varying temperature in Celsius
    - salinity (array-like): Time-varying salinity
    """

    def __init__(self, temp_celsius, salinity, CO2_atm=constants.ATMOSPHERIC_CO2, d13C_atm=constants.ATMOSPHERIC_d13C, D14C_atm=constants.ATMOSPHERIC_D14C, surface_area=constants.SURFACE_AREA):
        """
        Initialize all variables needed for the box model.

        Parameters:
        - temp_celsius (array-like): Time-varying temperature in Celsius
        - salinity (array-like): Time-varying salinity
        """
        self.num_tracers = 5
        self.surface_area = surface_area  # m^2
        self.surface_volume = self.surface_area * constants.MIXED_LAYER_DEPTH  # m^3
        self.surface_mass = self.surface_volume * constants.SWD  # kg

        self.temp_celsius = np.array(temp_celsius) # Celsius
        self.temp_kelvin = self.temp_celsius + 273.15
        self.salinity_forcing = np.array(salinity)
        self.salinity = self.salinity_forcing[0] # PSU

        self.CO2_atm = CO2_atm  # ppm
        self.d13C_atm = d13C_atm  # per mil
        self.D14C_atm = D14C_atm  # per mil

        self.alkalinity = data_input.calc_talk(self.salinity_forcing[0], pref="ne")
        
        initial_carbonate_system = data_input.compute_initial_conditions(
            pCO2=self.CO2_atm, alk=self.alkalinity, temp=self.temp_celsius[0], sal=self.salinity_forcing[0]
        )

        self.DIC = initial_carbonate_system["dic"] # µmol/kg
        self.del_13C = 1 * self.DIC # d13C*DIC units
        self.del_14C = 0 * self.DIC # D14C*DIC units

        self.pco2_ocean = self.CO2_atm  # µatm
        self.initial_state = np.hstack((self.DIC, self.alkalinity, self.del_13C, self.del_14C,self.salinity))

        self.result = None
        self.time = None
        self.output = None
        self.last_printed_year = None

    def calculate_day_of_year(self, time):
        """
        Calculate the day of the year from the given time value.

        Parameters:
        - time (float): Time value in years.

        Returns:
        - int: Day of the year (1-365).
        """
        total_days = time * 365  # Convert time in years to days
        day_of_year = int(total_days % 365) + 1  # Ensure value cycles between 1 and 365
        return day_of_year
    
    def _calculate_gas_exchange(self, state_vector, day_of_year, temp_celsius, salinity):
        """
        Calculate gas exchange fluxes.

        Parameters:
        - state_vector (array): Current state variables
        - day_of_year (int): Day of the year
        - temp_celsius (float): Current temperature in Celsius
        - salinity (float): Current salinity

        Returns:
        - array: Changes in tracers due to gas exchange
        """
        return fluxes.gas_exchange(
            current_state=state_vector,
            num_tracers=self.num_tracers,
            CO2_atm=self.CO2_atm,
            d13C_atm=self.d13C_atm,
            D14C_atm=self.D14C_atm,
            surface_area=self.surface_area,
            surface_mass=self.surface_mass,
            temp_celsius=temp_celsius,
            salinity=salinity,
            day_of_year=day_of_year,
        )
    
    def _calculate_mixing(self, state_vector, day_of_year):
        """
        Calculate vertical mixing fluxes.

        Parameters:
        - state_vector (array): Current state variables
        - day_of_year (int): Day of the year

        Returns:
        - array: Changes in tracers due to mixing
        """
        return fluxes.vertical_mixing(
            current_state=state_vector,
            num_tracers=self.num_tracers,
            day_of_year=day_of_year,
        )
    
    def _calculate_biology(self, state_vector, day_of_year):
        """
        Calculate photosynthesis and respiration fluxes.

        Parameters:
        - state_vector (array): Current state variables
        - day_of_year (int): Day of the year

        Returns:
        - array: Changes in tracers due to biology
        """
        return fluxes.biology(
            current_state=state_vector,
            num_tracers=self.num_tracers,
            day_of_year=day_of_year,
        )  
    
    def _calculate_salinity_effects(self, state_vector, day_of_year, salinity):
        """
        Calculate changes in tracer concentrations due to changes in salinity.

        Parameters:
        - state_vector (array): Current state variables
        - day_of_year (int): Day of the year
        - salinity (float): Current salinity

        Returns:
        - array: Changes in tracers due to dilution
        """
        return fluxes.salinity_effects(
            current_state=state_vector,
            num_tracers=self.num_tracers,
            day_of_year=day_of_year,
            salinity_forcing=salinity,
        )

    def model(self, time, state_vector):
        """
        Calculate change in the model for one time step.

        Parameters:
        - time (float): Current time step
        - state_vector (array): Current state variables

        Returns:
        - d_dt (array): Rate of change of tracers
        """
        day_of_year = self.calculate_day_of_year(time)
        current_temp_celsius = np.interp(day_of_year, np.arange(len(self.temp_celsius)), self.temp_celsius)
        current_salinity = np.interp(day_of_year, np.arange(len(self.salinity_forcing)), self.salinity_forcing)

        current_year = int(time)
        experiment_year = current_year - constants.DEFAULT_SPIN_UP_TIME

        if experiment_year >= 1 and (self.last_printed_year is None or experiment_year > self.last_printed_year):
            print(f"Year: {experiment_year}")
            self.last_printed_year = experiment_year

    
        d_dt = np.zeros((self.num_tracers))
        
        d_dt_gasexchange = self._calculate_gas_exchange(state_vector, day_of_year, current_temp_celsius, current_salinity)
        d_dt_mixing = self._calculate_mixing(state_vector, day_of_year)
        d_dt_biology = self._calculate_biology(state_vector, day_of_year)
        d_dt_dilution = self._calculate_salinity_effects(state_vector, day_of_year, current_salinity)

        d_dt += d_dt_gasexchange
        d_dt += d_dt_mixing
        d_dt += d_dt_biology
        d_dt += d_dt_dilution

        return d_dt

    def run_model(self, simulation_length_years=constants.DEFAULT_SIMULATION_LENGTH_YEARS, 
                  num_steps=None, spin_up_time=constants.DEFAULT_SPIN_UP_TIME):
        """
        Run the model with ODE solver giving initial_state as initial condition.

        Parameters:
        - simulation_length_years (int): Total years for experiment
        - num_steps (int): Number of steps that output is calculated (daily)
        - spin_up_time (int): Number of years to exclude from the analysis (for model spin-up)
        """

        if num_steps is None:
            num_steps = 365 * (spin_up_time + simulation_length_years)

        start_time = time.time()

        # Set the times to store the computed solution
        self.time = np.linspace(0, simulation_length_years, num_steps) # in days
        
        print(f"Model spinning up for {spin_up_time} years...")

        # Solve the ODE
        self.result = solve_ivp(
            self.model,
            [0, simulation_length_years],
            self.initial_state,
            method="RK45",
            t_eval=self.time,
            vectorized=True,
            rtol=1e-6,
            atol=1e-8,
        )

        end_time = time.time()
        self.time = self.result.t
        self.output = self.result.y
        print("Simulation complete.")
        print(f"This {simulation_length_years} year run took {end_time - start_time:.2f} seconds.")

        # Drop spin-up years
        if spin_up_time > 0:
            one_year_steps = int(num_steps / simulation_length_years)
            spin_up_steps = spin_up_time * one_year_steps
            self.time = self.time[spin_up_steps:] - spin_up_time
            self.output = self.output[:, spin_up_steps:]
            self.temp_celsius = self.temp_celsius[spin_up_steps:]
            self.salinity_forcing = self.salinity_forcing[spin_up_steps:]

        data_output.save_output_to_file(
            time=self.time,
            output=self.output,
            salinity=constants.DEFAULT_SALINITY,
            temperature=constants.DEFAULT_TEMP_CELSIUS,
        )

if __name__ == "__main__":
    from utils.data_input import load_model_forcings
    from utils.plotting import plot_variables
    from utils import constants
    
    # Configuration
    spin_up_time = constants.DEFAULT_SPIN_UP_TIME
    simulation_length_years = constants.DEFAULT_SIMULATION_LENGTH_YEARS
    total_years = spin_up_time + simulation_length_years
    num_steps = 365 * total_years

    # Load forcings from North Atlantic outside Gulf of Maine based on Fig 6 in Cei et al. 2020
    temperature, salinity = load_model_forcings('ne', total_years)

    # Initialize the model
    model_instance = MixedLayerModel(temp_celsius=temperature, salinity=salinity)
    
    # Run the model
    model_instance.run_model(simulation_length_years=total_years, num_steps=num_steps, spin_up_time=spin_up_time)

    # Plot results
    plot_variables(
        time=model_instance.time,   
        output=model_instance.output,
        temp_celsius=model_instance.temp_celsius,
        salinity=model_instance.salinity_forcing, 
    )