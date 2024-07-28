# Carbon Isotope Enabled Box Model for an Idealized Surface Mixed Layer 

### Model Description

Our model simulates a generic water column with a mixed layer depth of 30 m and a surface area of 1 m² (these values can be adjusted in the `constants.py` file). The model assumes a fixed atmospheric pCO₂ of 395 μatm with a δ¹³C of -8‰ (adjustable in the `constants.py` file). It includes state variables for dissolved inorganic carbon (DIC; units µmol kg⁻¹), alkalinity (ALK; units µmol kg⁻¹), δ¹³C-DIC (δ notation, units per mil), and Δ¹⁴C-DIC (Δ notation, units per mil).

The model is originally based on the work by Lynch-Stieglitz et al. (1995):

> Lynch-Stieglitz, J., Stocker, T. F., Broecker, W. S., & Fairbanks, R. G. (1995). The influence of air-sea exchange on the isotopic composition of oceanic carbon: Observations and modeling. Global Biogeochemical Cycles, 9(4), 653–665. doi:10.1029/95GB02574

It has been extended to include realistic forcings based on Cai et al. (2020):

> Cai, W. J., Xu, Y. Y., Feely, R. A., et al. (2020). Controls on surface water carbonate chemistry along North American ocean margins. Nature Communications, 11, 2691. https://doi.org/10.1038/s41467-020-16530-z

**Forcings and Processes:**

- **Temperature and Salinity**: Derived from the Mercator 1/12° data-assimilated General Circulation Model, provided at daily resolution.
- **Wind Speed**: Prescribed seasonally: 7 m/s from March through September and 10.5 m/s from October through February.

**Vertical Mixing and Biological Production:**

1. **Vertical Mixing**: The fractional change in tracers due to vertical mixing is based on a 0.1 µmol kg⁻¹ day⁻¹ DIC flux, as used in Cai et al. (2020), normalized to the current DIC concentration. This fractional rate is then applied to the subsurface boundary conditions (DIC = 2200 µmol kg⁻¹; ALK = 2400 µmol kg⁻¹; δ¹³C-DIC = 1‰; Δ¹⁴C-DIC = -100‰), effectively scaling the flux according to the subsurface conditions. Vertical mixing is applied only during winter months (October through February), with no mixing in summer.

2. **Biological Production**: A similar normalization approach is used for biological production, with a maximum DIC flux of 0.1 µmol kg⁻¹ day⁻¹ from March through September and a reduced flux of 0.01 µmol kg⁻¹ day⁻¹ for the rest of the year. This flux is normalized to the average DIC concentration of 2050 µmol kg⁻¹, based on the maximum ΔDIC_bio (0.1 mmol C day⁻¹ m⁻³) reported by Cai et al. (2020).

3. **Export Production, CaCO₃ Production, and Dissolution**: The model includes the export of particulate organic carbon (POC) and the production and dissolution of CaCO₃. An organic export percentage of 20% is used, based on Henson et al. (2011). The PIC:POC rain ratio is set at 0.07, based on Sarmiento et al. (2002), and a daily CaCO₃ dissolution rate of 38% is applied, based on Hales and Emerson (1997).

**Isotopic Fractionation During Biology:**

Biological processes influence the carbon isotopic composition of the ocean (e.g., Fontugne and Duplessy, 1978). During photosynthesis, organisms preferentially utilize ¹²C (the lighter carbon isotope), enriching the surface ocean in ¹³C and relatively enriching the underlying waters in ¹²C upon remineralization. A constant isotopic fractionation of -22‰ for photosynthesis (Toggweiler and Sarmiento, 1985; Vogel et al., 1970) is used, and a small fractionation of +2‰ for calcium carbonate formation is applied, following Jahn et al. (2015).

**Isotopic Fractionation During Air-Sea Gas Exchange:**

The kinetic fractionation factor for CO₂ transfer across the air-sea interface is set to -0.5‰, following Siegenthaler and Munnich (1981) and Zhang (1995). Equilibrium fractionation between atmospheric CO₂ and dissolved CO₂ in seawater is temperature-dependent, calculated using the fractionation factor α_CO₂ = -0.373 / Tₖ + 1.00019 (Vogel et al., 1970). The fractionation between dissolved CO₂ and HCO₃⁻ in seawater, α_HCO₃, is similarly temperature-dependent, calculated as α_HCO₃ = -9.866 / Tₖ + 1.02412. These fractionation factors (α) are converted to isotopic enrichment factors (ε) to correspond to the delta notation (units per mil) for δ¹³C tracers. Since the isotopic fractionation is approximately twice as large for ¹⁴C as for ¹³C (Zeebe and Wolf-Gladrow, 2001), all isotopic enrichment factors for ¹⁴C are doubled.

**Equilibrium Assumptions:**

Given that chemical and isotopic equilibrium occurs on a timescale of seconds (Zeebe et al., 1999), much faster than the years to decades required for air-sea equilibrium of CO₂ and carbon isotopes (Broecker and Peng, 1982; Schmittner et al., 2013), we assume our carbon species in the mixed layer are in chemical and isotopic equilibrium for all timescales relevant to our experiments.

### Model Initiation and Iteration

The model is spun up for 5 years to ensure a steady state of isotopic equilibrium with the atmosphere. It is initialized with a δ¹³C-DIC of 1‰ and a Δ¹⁴C of 0‰. Following Cai et al. (2020), the carbonate system is initialized with values calculated from PyCO2SYS (Humphreys et al., 2022) based on a pCO₂(aq) in equilibrium with atmospheric values (395 μatm) and alkalinity calculated with salinity. ΔDIC_bio and ΔDIC_vertical are prescribed as fractional changes (values mentioned above), while ΔpCO₂_air-sea is calculated at each time step using the relationship ΔpCO₂_air-sea = 0.24 * k * K₀ * (pCO₂_t(aq) - pCO₂_atm). K₀ is the CO₂ gas solubility (Weiss, 1974) and k is defined as 0.251 * W² * (Sc/660)⁻⁰·⁵, where W is wind speed in m/s and Sc is the Schmidt number (Wanninkhof, 1992; 2014). For ΔpCO₂_air-sea, pCO₂_t(aq) is calculated using an iterative CO₂ solver. At each time step, DIC_t+1 is updated based on ΔDIC_bio, ΔDIC_vertical, and ΔpCO₂_air-sea. At the end of the simulation, PyCO2SYS is used to back-calculate pCO₂, pH, and Ω_aragonite based on the model's DIC, ALK, and temperature and salinity forcing.

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

## References

1. **Cai et al. (2020)**: Cai, W. J., Xu, Y. Y., Feely, R. A., et al. (2020). Controls on surface water carbonate chemistry along North American ocean margins. Nature Communications, 11, 2691. https://doi.org/10.1038/s41467-020-16530-z
2. **Henson et al. (2011)**: Henson, S. A., Sanders, R., & Madsen, E. (2012). Global patterns in efficiency of particulate organic carbon export and transfer to the deep ocean. Global Biogeochemical Cycles, 26(1). doi:10.1029/2011GB004099
3. **Sarmiento and Gruber (2002)**: Sarmiento, J. L., & Gruber, N. (2002). Sinks for anthropogenic carbon. Physics Today, 55(8), 30-36. doi:10.1063/1.1510279
4. **Hales and Emerson (1997)**: Hales, B., & Emerson, S. (1997). CaCO₃ dissolution in sediments of the Ceara Rise, western equatorial Atlantic. Geochimica et Cosmochimica Acta, 61(3), 501-514. doi:10.1016/S0016-7037(96)00363-3
5. **Fontugne and Duplessy (1978)**: Fontugne, M. R., & Duplessy, J. C. (1978). Carbon isotope ratio of marine plankton related to surface water masses. Earth and Planetary Science Letters, 41(3), 365-371.
6. **Toggweiler and Sarmiento (1985)**: Toggweiler, J. R., & Sarmiento, J. L. (1985). Glacial to interglacial changes in atmospheric carbon dioxide: The critical role of ocean surface water in high latitudes. In E. T. Sundquist & W. S. Broecker (Eds.), The carbon cycle and atmospheric CO2: Natural variations Archean to present (Vol. 32, pp. 163-184). American Geophysical Union.
7. **Vogel et al. (1970)**: Vogel, J. C., Grootes, P. M., & Mook, W. G. (1970). Isotopic fractionation between gaseous and dissolved carbon dioxide. Zeitschrift für Physik A Hadrons and Nuclei, 230(3), 225-238.
8. **Jahn et al. (2015)**: Jahn, A., Claussen, M., Brovkin, V., & Ganopolski, A. (2015). Quantifying the effect of CO2 concentration on the spread of arid zones during the late Miocene. Earth and Planetary Science Letters, 452, 223-231.
9. **Siegenthaler and Munnich (1981)**: Siegenthaler, U., & Munnich, K. O. (1981). 13C/12C fractionation during CO2 transfer from air to sea. In E. T. Sundquist & W. S. Broecker (Eds.), The Carbon Cycle and Atmospheric CO2: Natural Variations Archean to Present (Vol. 32, pp

## Contact

For any questions or issues, please open an issue in the GitHub repository.
