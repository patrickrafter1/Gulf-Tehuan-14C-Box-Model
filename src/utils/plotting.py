import matplotlib.pyplot as plt
import numpy as np
from . import data_output
import os

def plot_variables(time, DIC, ALK, d13C, temp_celsius, salinity,spin_up_time, title="Model results"):
    """
    Plot the variables temperature, salinity, pCO2, pH, DIC, and d13C over time.

    Parameters:
    - time: Array of time values
    - DIC: Array of DIC values (µmol / kg)
    - ALK: Array of ALK values (µmol / kg)
    - d13C: Array of d13C values (per mil)
    - temp_celsius: Array of temperature values (°C)
    - salinity: Array of salinity values (PSU)
    - title: Title of the plot
    """
    # Extend temperature and salinity to match the time length
    num_years = int(np.ceil(len(time) / 365))
    extended_temperature = np.tile(temp_celsius, num_years)
    extended_salinity = np.tile(salinity, num_years)

    # Ensure temperature and salinity arrays match the length of time
    extended_temperature = extended_temperature[:len(time)]
    extended_salinity = extended_salinity[:len(time)]

    carb_chem = data_output.compute_carbonate_system(DIC, ALK, extended_temperature, extended_salinity)

    pCO2 = carb_chem["pCO2"]
    pH = carb_chem["pH"]
    omega = carb_chem["omega"]

    fig, axs = plt.subplots(4, 1, figsize=(14, 12), sharex=True)

    # Temperature and Salinity
    ax1 = axs[0]
    ax1.plot(time, extended_salinity, label="Salinity", color='blue')
    ax1.set_ylabel("SSS", color='blue')
    ax1.tick_params(axis='y', labelcolor='blue')
    for tl in ax1.get_yticklabels():
        tl.set_color('blue')

    ax2 = ax1.twinx()
    ax2.plot(time, extended_temperature, label="Temperature °C", color='red')
    ax2.set_ylabel("SST (°C)", color='red')
    ax2.tick_params(axis='y', labelcolor='red')
    for tl in ax2.get_yticklabels():
        tl.set_color('red')

    # pCO2 and pH
    ax3 = axs[1]
    ax3.plot(time, pCO2, label="pCO2 µatm", color='blue')
    ax3.set_ylabel("${p}$CO$_2$ (µatm)", color='blue')
    ax3.tick_params(axis='y', labelcolor='blue')
    for tl in ax3.get_yticklabels():
        tl.set_color('blue')

    ax4 = ax3.twinx()
    ax4.plot(time, pH, label="pH", color='red')
    ax4.set_ylabel("pH", color='red')
    ax4.tick_params(axis='y', labelcolor='red')
    for tl in ax4.get_yticklabels():
        tl.set_color('red')

    # DIC
    ax5 = axs[2]
    ax5.plot(time, DIC, label="DIC µmol / kg", color='blue')
    ax5.set_ylabel("DIC µmol kg$^{-1}$", color='blue')
    ax5.tick_params(axis='y', labelcolor='blue')
    for tl in ax5.get_yticklabels():
        tl.set_color('blue')

    ax6 = ax5.twinx()
    ax6.plot(time, omega, label="Ω", color='red')
    ax6.set_ylabel("Ω$_{arag}$", color='red')
    ax6.tick_params(axis='y', labelcolor='red')
    for tl in ax6.get_yticklabels():
        tl.set_color('red')

    # d13C
    ax7 = axs[3]
    ax7.plot(time, d13C, label="d13C", color='blue')
    ax7.set_ylabel("$\delta^{13}$C$_{DIC}$", color='blue')
    ax7.tick_params(axis='y', labelcolor='blue')
    for tl in ax7.get_yticklabels():
        tl.set_color('blue')


    # Custom x-axis labels
    num_ticks = (num_years) * 2 + 1  # For Jan, Jun each year
    # tick_positions = np.linspace(spin_up_time, num_years+spin_up_time, num_ticks)
    tick_positions = np.linspace(0, num_years, num_ticks)
    tick_labels = []
    for year in range(num_years):
        tick_labels.append(f'January of year {year+1}')
        tick_labels.append('')
    tick_labels.append(f'January of year {num_years+1}')
    
    plt.xticks(tick_positions, tick_labels, rotation=45, ha='right')

    # Title and labels
    ax1.set_title(title)
    plt.xlabel("Time")
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    # Save the figure
    plot_dir = os.path.join(os.getcwd(), "data/plots")
    os.makedirs(plot_dir, exist_ok=True)
    plot_path = os.path.join(plot_dir, "model_results.png")
    plt.savefig(plot_path)
    print(f"Plot saved to '{plot_path}'")
    plt.close()