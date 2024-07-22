import pandas as pd
import os
import PyCO2SYS as pyco2

def save_output_to_file(time, output, label, salinity, temperature):
    """
    Save computed chemical properties to a file.

    Parameters:
    - time (array-like): Array of time values.
    - output (matrix): Matrix containing DIC, Alkalinity, d13C, and D14C tracers.
    - label (str): Label for the output file.
    """
    df = pd.DataFrame(
        {
            "time": time,
            "DIC": output[0],
            "ALK": output[1],
            "d13C": output[2] / output[0],
            "D14C": output[3] / output[0],
            "salinity": salinity,
            "temperature": temperature,
        }
    )

    current_directory = os.getcwd()
    output_file_path = os.path.join(current_directory, f"data/model_output/initial_conditions_{label}.txt")

    df.to_csv(output_file_path, sep="\t", float_format="%.2f", index=False)
    print(f"Output printed to '{output_file_path}'")

def compute_carbonate_system(DIC, ALK, temperature, salinity):
    """
    Compute the carbonate system properties using PyCO2SYS.

    Parameters:
    - DIC: Array of DIC values (µmol / kg)
    - ALK: Array of ALK values (µmol / kg)
    - temperature: Array of temperature values (°C)
    - salinity: Array of salinity values (PSU)

    Returns:
    - Dictionary containing pCO2, pH, and saturation state of aragonite (omega)
    """

    kwargs = {
        "par1": ALK,
        "par2": DIC,
        "par1_type": 1,  # Alkalinity
        "par2_type": 2,  # DIC
        "temperature": temperature,
        "salinity": salinity,
    }

    try:
        result = pyco2.sys(**kwargs)
    except Exception as e:
        print(f"Error in PyCO2SYS: {e}")
        return None

    return {
        "pCO2": result["pCO2"],
        "pH": result["pH"],
        "omega": result["saturation_aragonite"],
    }
