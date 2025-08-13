# Carbon Isotope Enabled Box Model for an Idealized Surface Mixed Layer 
#Hi dummy

### Model Description

Our model simulates a generic water column with a mixed layer depth of 30 m and a surface area of 1 m² (these values can be adjusted in the `constants.py` file). The model assumes a fixed atmospheric pCO₂ of 395 μatm with a δ¹³C of -8‰ (adjustable in the `constants.py` file). It includes state variables for dissolved inorganic carbon (DIC; units µmol kg⁻¹), alkalinity (ALK; units µmol kg⁻¹), δ¹³C-DIC (δ notation, units per mil), and Δ¹⁴C-DIC (Δ notation, units per mil).

The model is originally based on the work by Lynch-Stieglitz et al. (1995):

> Lynch-Stieglitz, J., Stocker, T. F., Broecker, W. S., & Fairbanks, R. G. (1995). The influence of air-sea exchange on the isotopic composition of oceanic carbon: Observations and modeling. Global Biogeochemical Cycles, 9(4), 653–665. doi:10.1029/95GB02574

It has been extended to include realistic forcings from Cai et al. (2020):

> Cai, W. J., Xu, Y. Y., Feely, R. A., et al. (2020). Controls on surface water carbonate chemistry along North American ocean margins. Nature Communications, 11, 2691. https://doi.org/10.1038/s41467-020-16530-z

**Forcings:** 

The following forcings are applied to the model to simulate realistic oceanic conditions:

- **Temperature and Salinity**: Derived from the Mercator 1/12° data-assimilated General Circulation Model, provided at daily resolution.
- **Wind Speed**: Prescribed seasonally: 7 m/s from March through September and 10.5 m/s from October through February.
- **Net Community Production**: Derived from weekly averages of daily changes in satellite-derived chlorophyll (Chl) for the year 2015 using the daily MODIS aqua 4km level 3 product. The NCP values are interpolated to daily time series, with a maximum NCP of ~0.1 mmol C day⁻¹ m⁻³. More details can be found in Cai et al. (2020).
- **Vertical Mixing**: Vertical mixing during winter months (October through February) is simulated based on the fixed vertical mixing flux from Cai et al. (2020) of 0.1 mmol m⁻³ day⁻¹. We convert the mixing flux to a unitless daily mixing rate by dividing it by the current surface DIC concentration. This rate is then used to calculate the flux of tracers (DIC, ALK, δ¹³C, and Δ¹⁴C) between the subsurface (based on subsurface boundary conditions of DIC = 2200 µmol kg⁻¹; ALK = 2400 µmol kg⁻¹; δ¹³C-DIC = 1‰; Δ¹⁴C-DIC = -100‰) and surface layers. Under the assumption of steady-state conditions, the same amount of tracers mixed into the surface layer from the subsurface is also removed. No mixing occurs during summer (March through September).

**Isotopic Fractionation During Biology:**

Biological processes influence the carbon isotopic composition of the ocean (e.g., Fontugne and Duplessy, 1978). During photosynthesis, organisms preferentially utilize ¹²C (the lighter carbon isotope), enriching the surface ocean in ¹³C and relatively enriching the underlying waters in ¹²C upon remineralization. A constant isotopic fractionation of -22‰ for photosynthesis (Toggweiler and Sarmiento, 1985; Vogel et al., 1970) is used, and a small fractionation of +2‰ for calcium carbonate formation is applied, following Jahn et al. (2015).

**Isotopic Fractionation During Air-Sea Gas Exchange:**

The kinetic fractionation factor (α<sub>k</sub>) for CO₂ transfer across the air-sea interface is set to 0.9995, following Siegenthaler and Munnich (1981) and Zhang (1995). The equilibrium fractionation between atmospheric CO₂ and dissolved CO₂ in seawater is temperature-dependent, calculated using the fractionation factor α<sub>CO₂</sub> = -0.373 / Tₖ + 1.00019 (Vogel et al., 1970). Therefore, the overall fractionation factor for air-sea carbon transfer is α<sub>k</sub> * α<sub>CO₂</sub>.

Following Lynch-Stieglitz et al. (1995), the temperature-dependent equilibrium fractionation between dissolved CO₂ and DIC is approximated using the fractionation factor between dissolved CO₂ and HCO₃⁻, given by α<sub>HCO₃</sub> = -9.866 / Tₖ + 1.02412. Thus, the fractionation factor for sea-air carbon transfer is calculated as α<sub>k</sub> * α<sub>HCO₃</sub>.

These fractionation factors (α) are converted to isotopic enrichment factors (ε) to align with the delta notation (units per mil) used for δ¹³C tracers. Considering that isotopic fractionation for ¹⁴C is approximately twice as large as for ¹³C (Zeebe and Wolf-Gladrow, 2001), all isotopic enrichment factors for ¹⁴C are doubled.

**Equilibrium Assumptions:**

Given that chemical and isotopic equilibrium occurs on a timescale of seconds (Zeebe et al., 1999), much faster than the years to decades required for air-sea equilibrium of CO₂ and carbon isotopes (Broecker and Peng, 1982; Schmittner et al., 2013), we assume our carbon species in the mixed layer are in chemical and isotopic equilibrium for all timescales relevant to our experiments.

### Model Initiation and Iteration

The model is spun up for 5 years to ensure a steady state of isotopic equilibrium with the atmosphere. It is initialized with a δ<sup>13</sup>C-DIC of 1‰ and a Δ<sup>14</sup>C of 0‰. Following Cai et al. (2020), the carbonate system is initialized with values calculated from PyCO2SYS (Humphreys et al., 2022) based on a pCO₂(aq) in equilibrium with atmospheric values (395 μatm) and alkalinity calculated based on salinity. ΔDIC<sub>bio</sub> and ΔDIC<sub>vertical</sub> are prescribed as fractional changes (values mentioned above), while ΔpCO₂<sub>air-sea</sub> is calculated at each time step using the relationship ΔpCO₂<sub>air-sea</sub> = 0.24 * k * K₀ * (pCO₂<sub>t(aq)</sub> - pCO₂<sub>atm</sub>). K₀ is the CO₂ gas solubility (Weiss, 1974) and k is defined as 0.251 * W² * (Sc/660)⁻⁰·⁵, where W is wind speed in m/s and Sc is the Schmidt number (Wanninkhof, 1992; 2014). For ΔpCO₂<sub>air-sea</sub>, pCO₂<sub>t(aq)</sub> is calculated using an iterative CO₂ solver. At each time step, DIC<sub>t+1</sub> is updated based on ΔDIC<sub>bio</sub>, ΔDIC<sub>vertical</sub>, and ΔpCO₂<sub>air-sea</sub>. At the end of the simulation, PyCO2SYS is used to back-calculate pCO₂, pH, and Ω<sub>aragonite</sub> based on the model's DIC, ALK, and temperature and salinity forcing.

## Running the Model

To run the model using `conda`, follow these steps:

1. **Ensure `conda` or `Anaconda` is installed:**

   If you don't have `conda` installed, you can download and install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/products/distribution).

2. **Clone the repository:**

    ```bash
    git clone https://github.com/RyanAGreen/mld-box-model-c-isotopes.git
    cd mld-box-model-c-isotopes
    ```

3. **Create and activate a `conda` environment:**

    - Create a new environment using the provided `environment.yml` file inside the cloned repository:

      ```bash
      conda env create -f environment.yaml
      ```

    - Activate the environment:

      ```bash
      conda activate mld-box-model
      ```

4. **Adjust model parameters:**
   
    You can modify the experiment length, mixed layer properties, atmospheric conditions, prescribed forcings, and other parameters within `src/utils/constants.py` as needed.

5. **Run the model:**

    ```bash
    python src/model.py
    ```

    Note that this will take a minute or two the first time the code is run. After that, each simulation should take only a few seconds (depending on the length of the experiment). An overview figure will be generated and saved as `model_results.png` in the `data/plots` directory. Additionally, a text file containing the model results will be saved as `model_results.txt` in the `data` directory for further analysis.

6. **Deactivate the `conda` environment:**

    When you're done working with the model, deactivate the environment:

    ```bash
    conda deactivate
    ```

7. **Note**: If you want to work on the model later, remember to reactivate the `conda` environment before running any code:

    ```bash
    conda activate mld-box-model
    ```


## Code Structure

The table below explains the files and directories included in this project.

| Directory/File | Description |
|---|---|
| `README.md` | Project documentation (this file) |
| `LICENSE` | License information |
| `environment.yml` | Conda environment file for dependencies |
| `data` | Directory for data used or generated by the model |
| `src` | Source code for the model |
| `src/model.py` | Main model code |
| `src/utils` | Utility functions or modules |
| `src/utils/constants.py` | Constants used in the model |
| `src/utils/fluxes.py` | Flux calculations for the model |
| `src/utils/data_input.py` | Data input functions |
| `src/utils/data_output.py` | Data output functions |
| `src/utils/plotting.py` | Plotting functions |
| `src/utils/indata/mercator_tseries.h5` | Temperature, salinity, and net community production forcing data |
| `src/__init__.py` | Python package initialization file |


## Output

The model generates:

- **Overview Figure**: Saved in the `data/plots` directory. This figure shows the seasonal salinity and SST forcing, DIC and δ13C model tracers, and pCO2, pH, and Ωarag generated with PyCO2SYS from the simulated DIC and ALK with the seasonal salinity and SST forcings. This figure is similar to Fig 6a,b,c from [Cai et al., 2020](https://doi.org/10.1038/s41467-020-16530-z).
- **Text File of Results**: Saved in the `data` directory. This only contains the model tracers (DIC, alkalinity, d13C, and D14C), salinity and SST forcing, and time in units of years. Other carbonate parameters can be calculated from these DIC, alkalinty, salinity and SST using PyCO2SYS/CO2SYS or another CO2 solver.  

Both outputs are designed for further data analysis and visualization.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## References

### References

1. **Cai et al. (2020)**: Cai, W. J., Xu, Y. Y., Feely, R. A., et al. (2020). Controls on surface water carbonate chemistry along North American ocean margins. Nature Communications, 11, 2691. https://doi.org/10.1038/s41467-020-16530-z
2. **Henson et al. (2011)**: Henson, S. A., Sanders, R., & Madsen, E. (2012). Global patterns in efficiency of particulate organic carbon export and transfer to the deep ocean. Global Biogeochemical Cycles, 26(1). doi:10.1029/2011GB004099
3. **Sarmiento and Gruber (2002)**: Sarmiento, J. L., & Gruber, N. (2002). Sinks for anthropogenic carbon. Physics Today, 55(8), 30-36. doi:10.1063/1.1510279
4. **Hales and Emerson (1997)**: Hales, B., & Emerson, S. (1997). CaCO₃ dissolution in sediments of the Ceara Rise, western equatorial Atlantic. Geochimica et Cosmochimica Acta, 61(3), 501-514. doi:10.1016/S0016-7037(96)00363-3
5. **Fontugne and Duplessy (1978)**: Fontugne, M. R., & Duplessy, J. C. (1978). Carbon isotope ratio of marine plankton related to surface water masses. Earth and Planetary Science Letters, 41(3), 365-371.
6. **Toggweiler and Sarmiento (1985)**: Toggweiler, J. R., & Sarmiento, J. L. (1985). Glacial to interglacial changes in atmospheric carbon dioxide: The critical role of ocean surface water in high latitudes. In E. T. Sundquist & W. S. Broecker (Eds.), The carbon cycle and atmospheric CO2: Natural variations Archean to present (Vol. 32, pp. 163-184). American Geophysical Union.
7. **Vogel et al. (1970)**: Vogel, J. C., Grootes, P. M., & Mook, W. G. (1970). Isotopic fractionation between gaseous and dissolved carbon dioxide. Zeitschrift für Physik A Hadrons and Nuclei, 230(3), 225-238.
8. **Jahn et al. (2015)**: Jahn, A., Claussen, M., Brovkin, V., & Ganopolski, A. (2015). Quantifying the effect of CO2 concentration on the spread of arid zones during the late Miocene. Earth and Planetary Science Letters, 452, 223-231.
9. **Siegenthaler and Munnich (1981)**: Siegenthaler, U., & Munnich, K. O. (1981). 13C/12C fractionation during CO2 transfer from air to sea. In E. T. Sundquist & W. S. Broecker (Eds.), The Carbon Cycle and Atmospheric CO2: Natural Variations Archean to Present (Vol. 32, pp. 259-269). American Geophysical Union. doi:10.1029/GM032p0259
10. **Zhang (1995)**: Zhang, J. (1995). The 13C/12C ratios of atmospheric CO2 in Beijing, China. Geochimica et Cosmochimica Acta, 59(4), 831-837. doi:10.1016/0016-7037(95)90522-Z
11. **Zeebe and Wolf-Gladrow (2001)**: Zeebe, R. E., & Wolf-Gladrow, D. A. (2001). CO2 in seawater: equilibrium, kinetics, isotopes. Elsevier Science.
12. **Broecker and Peng (1982)**: Broecker, W. S., & Peng, T. H. (1982). Tracers in the sea. Lamont-Doherty Geological Observatory.
13. **Schmittner et al. (2013)**: Schmittner, A., Mix, A. C., & Pisias, N. G. (2013). Regional patterns of deglacial warming. Quaternary Science Reviews, 68, 22-34.
14. **Zeebe et al. (1999)**: Zeebe, R. E., Sanyal, A., Ortiz, J. D., & Wolf-Gladrow, D. A. (1999). A theoretical study of the kinetics of the CO2 system in seawater. Marine Chemistry, 65(3-4), 135-152.
15. **Humphreys et al. (2022)**: Humphreys, M. P., Achterberg, E. P., & Briggs, E. M. (2022). PyCO2SYS v1.7: marine carbonate system calculations in Python. Geoscientific Model Development, 15(4), 1377-1392. doi:10.5194/gmd-15-1377-2022
16. **Weiss (1974)**: Weiss, R. F. (1974). Carbon dioxide in water and seawater: the solubility of a non-ideal gas. Marine Chemistry, 2, 203-215.
17. **Wanninkhof (1992)**: Wanninkhof, R. (1992). Relationship between wind speed and gas exchange over the ocean. Journal of Geophysical Research: Oceans, 97(C5), 7373-7382.
18. **Wanninkhof (2014)**: Wanninkhof, R. (2014). Relationship between wind speed and gas exchange over the ocean revisited. Limnology and Oceanography: Methods, 12(6), 351-362.

## Contact

For any questions or issues, please open an issue in the GitHub repository.
