import matplotlib.pyplot as plt
import numpy as np
from . import data_output
from . import constants
import os

def plot_variables(time, output, temp_celsius, salinity, title="Model results"):
    """
    Plot the variables temperature, salinity, pCO2, pH, DIC, and d13C over time.

    Parameters:
    - time: Array of time values
    - ouput: matrix of model output (µmol / kg)
    - temp_celsius: Array of temperature forcing (°C)
    - salinity: Array of salinity forcing (PSU) 
    - title: Title of the plot
    """
    # # Extend temperature and salinity to match the time length
    # num_years = int(np.ceil(len(time) / 365))
    num_years = constants.TOTAL_YEARS
    # extended_temperature = np.tile(temp_celsius, num_years)
    # extended_salinity = np.tile(salinity, num_years)

    # # Ensure temperature and salinity arrays match the length of time
    # extended_temperature = extended_temperature[:len(time)]
    # extended_salinity = extended_salinity[:len(time)]

    if len(temp_celsius) < len(time):
        temp_celsius = np.repeat(temp_celsius, len(time))
        salinity = np.repeat(salinity, len(time))

    # DIC, ALK, and d13C from model output
    DIC = output[0]
    ALK = output[1]
    d13C = output[2] / output[0]
    D14C = output[3] / output[0]
    

    extended_salinity = salinity
    extended_temperature = temp_celsius

    # Calculate carbonate system variables from PyCO2SYS
    carb_chem = data_output.compute_carbonate_system(DIC, ALK, extended_temperature, extended_salinity)
    pCO2 = carb_chem["pCO2"]
    pH = carb_chem["pH"]
    omega = carb_chem["omega"]

    fig, axs = plt.subplots(4, 1, figsize=(6, 11), sharex=True)

    # Temperature and Salinity
    ax1 = axs[0]
    ax1.plot(time, extended_salinity, label="Salinity", color='blue')
    ax1.set_ylabel("SSS", color='blue')
    ax1.tick_params(axis='y', labelcolor='blue')
    ax1.set_ylim(33.7,34)
    for tl in ax1.get_yticklabels():
        tl.set_color('blue')

    ax2 = ax1.twinx()
    ax2.plot(time, extended_temperature, label="Temperature °C", color='red')
    ax2.set_ylabel("SST (°C)", color='red')
    ax2.tick_params(axis='y', labelcolor='red')
    ax2.set_ylim(27, 29)
    for tl in ax2.get_yticklabels():
        tl.set_color('red')

    # pCO2 and pH
    ax3 = axs[1]
    ax3.plot(time, pCO2, label="pCO2 µatm", color='blue')
    ax3.set_ylabel("${p}$CO$_2$ (µatm)", color='blue')
    ax3.set_ylim(100, 1500)
    ax3.tick_params(axis='y', labelcolor='blue')
    for tl in ax3.get_yticklabels():
        tl.set_color('blue')

    ax4 = ax3.twinx()
    ax4.plot(time, pH, label="pH", color='red')
    ax4.set_ylabel("pH", color='red')
    ax4.tick_params(axis='y', labelcolor='red')
    ax4.set_ylim(7.5, 8.5)
    for tl in ax4.get_yticklabels():
        tl.set_color('red')

    # DIC & ALK
    ax5 = axs[2]
    ax5.plot(time, DIC, label="DIC µmol / kg", color='blue')
    ax5.set_ylabel("DIC µmol kg$^{-1}$", color='blue')
    ax5.tick_params(axis='y', labelcolor='blue')
    ax5.set_ylim(1600,2300)
    for tl in ax5.get_yticklabels():
        tl.set_color('blue')

    ax6 = ax5.twinx()
    ax6.plot(time, ALK, label="alkalinity", color='red')
    ax6.set_ylabel("ALK µmol kg$^{-1}$", color='red')
    ax6.tick_params(axis='y', labelcolor='red')
    ax6.set_ylim(2000, 2500)
    for tl in ax6.get_yticklabels():
        tl.set_color('red')

    
    # ax6 = ax5.twinx()
    # ax6.plot(time, omega, label="Ω", color='red')
    # ax6.set_ylabel("Ω$_{arag}$", color='red')
    # ax6.tick_params(axis='y', labelcolor='red')
    # ax6.set_ylim(1.3, 2.7)
    # for tl in ax6.get_yticklabels():
    #     tl.set_color('red')

    # d13C
    ax7 = axs[3]
    ax7.plot(time, d13C, label="d13C", color='blue')
    ax7.set_ylabel("$\delta^{13}$C$_{DIC}$", color='blue')
    ax7.tick_params(axis='y', labelcolor='blue')
    ax7.set_ylim(-1,2)
    for tl in ax7.get_yticklabels():
        tl.set_color('blue')

    # D14C
    ax8 = ax7.twinx()
    #ax8 = axs[4]
    ax8.plot(time, D14C, label="D14C", color='red')
    ax8.set_ylabel("$Δ^{14}$C$_{DIC}$", color='red')
    ax8.tick_params(axis='y', labelcolor='red')
    ax8.set_ylim(-300,300)
    for tl in ax8.get_yticklabels():
        tl.set_color('red')
    
    # Custom x-axis labels
    tick_positions = []
    tick_labels = []
    for year in range(num_years):
        for month in ["Jan", "Jul",]:
            tick_positions.append(year + (["Jan", "Jul"].index(month) / 2))
            tick_labels.append(month)
    tick_positions.append(num_years)
    tick_labels.append("Jan")

    ax7.set_xticks(tick_positions)
    ax7.set_xticklabels(tick_labels, rotation=45, ha='right')

    # Set x-axis limit to end at the second April
  # ax7.set_xlim([0, 1 + 1/3])

    # Add vertical grid lines for each tick
    for pos in tick_positions:
        for ax in axs:
            ax.axvline(pos, color='gray', linestyle='-', linewidth=0.5)

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