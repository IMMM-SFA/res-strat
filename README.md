# res-strat

A Multi-layer Reservoir Thermal Stratification Module for Earth System Models

## Getting started with `res-strat`

- Clone this repository using `git clone https://github.com/IMMM-SFA/res-strat.git`
- Navigate into your `res-strat` repository directory
- Use the examples provided in the `data/inputs` directory or replace the contents of the input directories and files with your own.
- Compile the model by running:  `gfortran -o res-strat src/ReservoirStratification.F90`
- Run the model by simply executing `res-strat`. It may take an hour to simulate a reservoir for ten years on daily time step.
- Outputs will be created in the `outputs/depth` and `outputs/stratification` directories.

## Contents

Input data format is ASCII tab delimited and rows represent daily (or desired time step) data.  No header should be present.  File names correspond to the names in the `data/inputs/reservoirs.txt` file.  Currently, reservoir file names represent the coordinates represent the centroid of the grid where reservoirs are present in a `MOSART-heat-wm` simulation over the CONUS. This is only for ease of recognition; any name can be used.

| File or Directory | IO | Description | Fields (Units) |
| -- | -- | -- | -- |
| `data/inputs/reservoirs.txt` | Input | List of reservoirs to be simulated | file name of each reservoir file |
| `data/inputs/geometry` | Input | Directory containing geometry datasets | geometry code (-), mean depth (m), height (m), length (km),width (km), vol. error (mcm), area error (sqkm), vol. coefficient (-), area coefficient (-), vol. difference (mcm), area difference (sqkm), no. of layers (-) |
| `data/inputs/flow` | Input | Directory containing reservoir inflow, inflow temperature, and reservoir outflow | inflow (m3/s), inflow temperature (Kelvin), outflow (m3/s) |
| `data/inputs/forcing` | Input | Directory containing NLDAS2 forcing | cosine of zenith (-), long wave absorbed (w/m2), shortwave (w/m2), relative humidity (%), air temperature (Kelvin), wind speed (m/s) |
| `data/outputs/depth` | Output | Output directory for reservoir depth | Depth representation |
| `data/outputs/stratification` | Output | Output directory for stratified reservoir temperature daily simulation period 2001-2010 | one row representing one day (or desired time step) for all layers used. The last column represents the surface |

## Acknowledgement

This research was supported by the US Department of Energy, Office of Science, as part of research in Multi-Sector Dynamics, Earth and Environmental System Modeling Program.  
PNNL is operated for DOE by Battelle Memorial Institute under contract DE-AC05-76RL01830.
A complete reservoir geometry dataset is available at https://doi.org/10.5281/zenodo.1322884.

## Citation

Yigzaw, W., Hong-Yi Li, Xing Fang, L. Ruby Leung, Nathalie Voisin, Mohamad I. Hejazi, Yonas Demissie. 2019. 'A Multi-layer Reservoir Thermal Stratification Module for Earth System Models'. Journal of Advances in Modeling Earth Systems.

## Contact

For any questions about the database please contact: Dr. HongYi Li (hli57@uh.edu) or Dr. Wondmagegn Yigzaw (wyigzaw@uh.edu)
Wonder of Waters Lab at University of Houston.
webpage:  http://wowuoh.wixsite.com/home
