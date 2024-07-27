# Carbon Isotope Enabled Box Model for an Idealized Surface Mixed Layer 

### Model Description

Our model simulates a generic water column with a mixed layer depth of 30 m and a surface area of 1 m² (these values can be adjusted in the 'constants.py' file). The model assumes a fixed atmospheric pCO₂ of 395 μatm with a δ<sup>13</sup>C of -8‰ (these values can be adjusted in the constants.py file). It includes state variables for dissolved inorganic carbon (DIC; units µmol kg⁻¹), alkalinity (ALK; units µmol kg⁻¹), δ<sup>13</sup>C-DIC (δ notation, units per mil), and Δ<sup>14</sup>C-DIC (Δ notation, units per mil).

The model is originally based on the work by Lynch-Stieglitz et al. (1995):

> Lynch-Stieglitz, J., Stocker, T. F., Broecker, W. S., & Fairbanks, R. G. (1995). The influence of air-sea exchange on the isotopic composition of oceanic carbon: Observations and modeling. Global Biogeochemical Cycles, 9(4), 653–665. doi:10.1029/95GB02574

It has been extended to include realistic forcings based on Cai et al. (2020):

> Cai, W. J., Xu, Y. Y., Feely, R. A., et al. (2020). Controls on surface water carbonate chemistry along North American ocean margins. Nature Communications, 11, 2691. https://doi.org/10.1038/s41467-020-16530-z

Temperature and salinity forcings are derived from the Mercator 1/12° data-assimilated General Circulation Model, with a daily resolution. Wind speed is prescribed seasonally: 7 m/s from March through September and 10.5 m/s from October through February. A fractional change in tracers due to vertical mixing is based on a 0.1 µmol kg⁻¹ day⁻¹ DIC flux used in Cai et al. (2020) and the current DIC concentration at each time step. This fractional flux is applied to all tracers during winter months (October through February) with no vertical mixing in summer. A similar approach is used for biological production, with a fractional flux of 0.1 µmol kg⁻¹ day⁻¹ DIC from March through September and 0.01 µmol kg⁻¹ day⁻¹ DIC for the rest of the year, based on maximum ΔDIC<sub>bio</sub> (0.1 mmol C day⁻¹ m⁻³) from Cai et al. (2020).

Biological processes influence the carbon isotopic composition of the ocean (e.g., Fontugne and Duplessy, 1978). During photosynthesis, organisms preferentially utilize ¹²C (the lighter carbon isotope), enriching the surface ocean in ¹³C and relatively enriching the underlying waters in ¹²C upon remineralization. We use a constant isotopic fractionation of -22‰ for photosynthesis (Toggweiler and Sarmiento, 1985; Vogel et al., 1970) and a small fractionation of +2‰ for calcium carbonate formation, following Jahn et al. (2015).

We account for both kinetic and equilibrium isotopic fractionation during air-sea gas exchange. The kinetic fractionation factor for CO₂ transfer across the air-sea interface is set to -0.5‰, following Siegenthaler and Munnich (1981) and Zhang (1995). Equilibrium fractionation between atmospheric CO₂ and dissolved CO₂ in seawater is temperature-dependent, calculated using the fractionation factor α<sub>CO2</sub> = -0.373 / T<sub>K</sub> + 1.00019 (Vogel et al., 1970). The fractionation between dissolved CO₂ and HCO₃⁻ in seawater, α<sub>HCO3</sub>, is similarly temperature-dependent, calculated as α<sub>HCO3</sub> = -9.866 / T<sub>K</sub> + 1.02412. These fractionation factors (α) are converted to isotopic enrichment factors (ε) to correspond to the delta notation (units per mil) for δ<sup>13</sup>C tracers. Since the isotopic fractionation is approximately twice as large for ¹⁴C as for ¹³C (Zeebe and Wolf-Gladrow, 2001), all isotopic enrichment factors for ¹⁴C are doubled.

Given that chemical and isotopic equilibrium occurs on a timescale of seconds (Zeebe et al., 1999), much faster than the years to decades required for air-sea equilibrium of CO₂ and carbon isotopes (Broecker and Peng, 1982; Schmittner et al., 2013), we assume our carbon species in the mixed layer are in chemical and isotopic equilibrium for all timescales relevant to our experiments.

### Model Initiation and Iteration

The model is spun up for 5 years to ensure a steady state of isotopic equilibrium with the atmosphere. It is initialized with a δ<sup>13</sup>C-DIC of 1‰ and a Δ<sup>14</sup>C of 0‰. Following Cai et al. (2020), the carbonate system is initialized with values calculated from PyCO2SYS (Humphreys et al., 2022) based on a pCO₂(aq) in equilibrium with atmospheric values (395 μatm) and alkalinity calculated with salinity. ΔDIC<sub>bio</sub> and ΔDIC<sub>vertical</sub> are prescribed as fractional changes (values mentioned above), while ΔpCO₂<sub>air-sea</sub> is calculated at each time step using the relationship ΔpCO₂<sub>air-sea</sub> = 0.24 * k * K₀ * (pCO₂<sub>t(aq)</sub> - pCO₂<sub>atm</sub>). K₀ is the CO₂ gas solubility (Weiss, 1974) and k is defined as 0.251 * W² * (Sc/660)⁻⁰·⁵, where W is wind speed in m/s and Sc is the Schmidt number (Wanninkhof, 1992; 2014). For ΔpCO₂<sub>air-sea</sub>, pCO₂<sub>t(aq)</sub> is calculated using an iterative CO₂ solver. At each time step, DIC<sub>t+1</sub> is updated based on ΔDIC<sub>bio</sub>, ΔDIC<sub>vertical</sub>, and ΔpCO₂<sub>air-sea</sub>. At the end of the simulation, PyCO2SYS is used to back-calculate pCO₂, pH, and Ω<sub>aragonite</sub> based on the model's DIC, ALK, and temperature and salinity forcing.


## Running the Model

To run the model, follow these steps:

1. **Clone the repository:**

    ```bash
    git clone https://github.com/RyanAGreen/mld-box-model-c-isotopes.git
    cd mld-box-model-c-isotopes
    ```

2. **Create and activate a virtual environment:**

    - **Using `venv` (Python 3 standard library):**

      ```bash
      python3 -m venv mld-box-model
      ```

    - **Activate the virtual environment:**

      - On Windows:

        ```bash
        .\mld-box-model\Scripts\activate
        ```

      - On macOS and Linux:

        ```bash
        source mld-box-model/bin/activate
        ```

3. **Install the required dependencies:**

    ```bash
    pip install -r requirements.txt
    ```

4. **Adjust model parameters:**
   
    You can also modify the experiment length, mixed layer properties, atmospheric conditions, the parameterized seasonal forcings, and all other parameters within `src/utils/constants.py`.

5. **Run the model:**

    ```bash
    python src/model.py
    ```
    An overview figure will be generated and saved as `model_results.png` in the `data/plots` directory. Additionally, a text file containing the model results will be saved as `model_results.txt` in the `data` directory for further analysis.

6. **Deactivate the virtual environment:**

    When you're done working with the model, deactivate the virtual environment:

    ```bash
    deactivate
    ```

## Code Structure

The table below explains the files and directories included in this project.

| Directory/File | Description |
|---|---|
| `README.md` | Project documentation (this file) |
| `LICENSE` | License information |
| `requirements.txt` | Dependencies for the project |
| `data` | Directory for data used or generated by the model |
| `src` | Source code for the model |
| `src/model.py` | Main model code |
| `src/utils` | Utility functions or modules |
| `src/utils/constants.py` | Constants used in the model |
| `src/utils/fluxes.py` | Flux calculations for the model |
| `src/utils/data_input.py` | Data input functions |
| `src/utils/data_output.py` | Data output functions |
| `src/utils/plotting.py` | Plotting functions |
| `src/utils/indata/mercator_tseries.h5` | Temperature and salinity data |
| `src/__init__.py` | Python package initialization file |


## Output

The model generates:

- **Overview Figure**: Saved in the `data/plots` directory. This figure shows the seasonal salinity and SST forcing, DIC and δ13C model tracers, and pCO2, pH, and Ωarag generated with PyCO2SYS from the simulated DIC and ALK with the seasonal salinity and SST forcings. This figure is similar to Fig 6a,b,c from [Cai et al., 2020](https://doi.org/10.1038/s41467-020-16530-z).
- **Text File of Results**: Saved in the `data` directory. This only contains the model tracers (DIC, alkalinity, d13C, and D14C), salinity and SST forcing, and time in units of years. Other carbonate parameters can be calculated from these DIC, alkalinty, salinity and SST using PyCO2SYS/CO2SYS or another CO2 solver.  

Both outputs are designed for further data analysis and visualization.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Contact

For any questions or issues, please open an issue in the GitHub repository.
