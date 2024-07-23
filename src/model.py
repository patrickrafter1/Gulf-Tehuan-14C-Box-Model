"""
Carbon isotope box model for an idealized mixed layer.
Inspired by Lynch-Stieglitz et al. (1995)
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
    The model includes DIC, ALK, d13C, and D14C tracers which evolve through time due to air-sea gas exchange.

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

    def __init__(self, temp_celsius, salinity):
        """
        Initialize all variables needed for the box model.

        Parameters:
        - temp_celsius (array-like): Time-varying temperature in Celsius
        - salinity (array-like): Time-varying salinity
        """
        self.num_tracers = 4
        self.surface_area = constants.SURFACE_AREA  # m^2
        self.surface_volume = self.surface_area * constants.MIXED_LAYER_DEPTH  # m^3
        self.surface_mass = self.surface_volume * constants.SWD  # kg

        self.temp_celsius = np.array(temp_celsius)
        self.temp_kelvin = self.temp_celsius + 273.15
        self.salinity = np.array(salinity)

        self.CO2_atm = 395  # ppm
        self.d13C_atm = -8  # per mil
        self.D14C_atm = 0  # per mil

        self.alkalinity = data_input.calc_talk(self.salinity[0], pref="se")
        
        initial_carbonate_system = data_input.compute_initial_conditions(
            pCO2=self.CO2_atm, alk=self.alkalinity, temp=self.temp_celsius[0], sal=self.salinity[0]
        )

        self.DIC = initial_carbonate_system["dic"] # µmol/kg

        self.del_13C = 1 * self.DIC # delta * concenration tracer units
        self.del_14C = 0 * self.DIC # delta * concenration tracer units
        self.initial_state = np.hstack((self.DIC, self.alkalinity, self.del_13C, self.del_14C))

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
        current_salinity = np.interp(day_of_year, np.arange(len(self.salinity)), self.salinity)

        # Print time only once per year
        current_year = int(time)+1
        if self.last_printed_year is None or current_year > self.last_printed_year:
            print(f"Year: {time:.2f}")
            self.last_printed_year = current_year  

        d_dt = np.zeros((self.num_tracers))
        
        # Air-Sea Gas Exchange
        d_dt_gasexchange = fluxes.gas_exchange(
            current_state=state_vector,
            num_tracers=self.num_tracers,
            CO2_atm=self.CO2_atm,
            d13C_atm=self.d13C_atm,
            D14C_atm=self.D14C_atm,
            surface_area=self.surface_area,
            surface_mass=self.surface_mass,
            temp_celsius=current_temp_celsius,
            salinity=current_salinity,
            day_of_year=day_of_year,
        )
        
        # Mixing
        d_dt_mixing = fluxes.vertical_mixing(num_tracers=self.num_tracers, day_of_year=day_of_year)
        
        # Biology
        d_dt_biology = fluxes.biology(current_state=state_vector, num_tracers=self.num_tracers, day_of_year=day_of_year)
        
        d_dt += d_dt_mixing
        d_dt += d_dt_biology
        d_dt += d_dt_gasexchange

        return d_dt

    def run_model(self, simulation_length_years, num_steps, spin_up_time=0):
        """
        Run the model with ODE solver giving initial_state as initial condition.

        Parameters:
        - simulation_length_years (int): Total years for experiment
        - num_steps (int): Number of steps that output is calculated
        """
        start_time = time.time()

        # Set the times to store the computed solution
        self.time = np.linspace(0, simulation_length_years, num_steps)

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

        print(f"This {simulation_length_years} year run took {end_time - start_time:.2f} seconds.")

        # Drop spin-up years
        if spin_up_time > 0:
            one_year_steps = int(num_steps / simulation_length_years)
            spin_up_steps = spin_up_time * one_year_steps
            self.time = self.time[spin_up_steps:] - spin_up_time
            self.output = self.output[:, spin_up_steps:]
            self.temp_celsius = self.temp_celsius[spin_up_steps:]
            self.salinity = self.salinity[spin_up_steps:]

        # Save the output to file
        data_output.save_output_to_file(
            time=self.time,
            output=self.output,
            salinity=self.salinity,
            temperature=self.temp_celsius,
        )


if __name__ == "__main__":
    from utils.data_input import load_seasonal_forcings
    from utils.plotting import plot_variables
    
    spin_up_years = 5  # Number of years to spin up the model
    simulation_length_years = 5  # Set the length of the simulation

    num_steps = 365 * (spin_up_years+simulation_length_years)  # Number of steps for output calculation

    # Load seasonal forcings
    temperature, salinity = load_seasonal_forcings(spin_up_years+simulation_length_years)  # Load for one year

    # Initialize the model with realistic forcings
    model_instance = MixedLayerModel(temp_celsius=temperature, salinity=salinity)
    model_instance.run_model(spin_up_years+simulation_length_years, num_steps, spin_up_time=spin_up_years)

    # After running the model and saving output to file
    plot_variables(
        time=model_instance.time,   
        DIC=model_instance.output[0],
        ALK=model_instance.output[1],
        d13C=model_instance.output[2]/model_instance.output[0],
        temp_celsius=model_instance.temp_celsius,
        salinity=model_instance.salinity, 
        spin_up_time=spin_up_years,
    )
