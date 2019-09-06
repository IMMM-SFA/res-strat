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

`src/ReservoirStratification.F90` - source code
`data/inputs/reservoirs.txt` - list of reservoirs to be simulated
`data/inputs/geometry` - directory containing geometry datasets
`data/inputs/flow` - directory containing reservoir inflow, inflow temperature, and reservoir outflow
`data/inputs/forcing` - directory containing NLDAS2 forcing
`data/outputs/depth` - output directory for reservoir depth
`data/outputs/stratification` - output directory for stratified reservoir temperature simulation period 2001-2010, daily

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
