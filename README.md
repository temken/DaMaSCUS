# DaMaSCUS

<img src="https://cloud.githubusercontent.com/assets/29034913/26834962/4092f75c-4ad7-11e7-86db-a359734ea2ef.png" width="425">

Dark Matter Simulation Code for Underground Scatterings

DaMaSCUS Version 1.0 06/06/2016

## GENERAL NOTES

- DaMaSCUS is a MC simulator of dark matter particles as they move through the Earth and scatter on terrestrial nuclei. 
- It allows to compute the local distortions of the DM phase space caused by collisions with nuclei. 
- The thusly distorted distribution functions and densities are used to give precise estimates of time-dependent signal rates for direct detection experiments and diurnal modulations.
- A full, realistic model of the Earth is implemented as well as the Earth's time-dependent velocity and orientation in the galactic frame.
- DaMaSCUS is written in C++ and fully parallelized (openMPI).

## CONTENT

The included folders are:

- `analysis`: This folder contains the relevant header and cpp files of the analysis code, which reads in the raw data from the data-folder and computes all relevant results, which then get stored in the results-folder.
- `data`: Once a simulation run is performed, the generated data will be stored here.
- `plots`: To visualize the results, created by the analysis module we include the small Mathematica package "DaMaSCUStoolbox" and an example notebook creating and saving plots.
- `results`: The analysis module saves its results and histograms here.
- `simulation`: This folder contains the simulation code and is the core of DaMaSCUS. You will find the relevant header and cpp files, together with the config file of input parameters.


## INSTALLATION AND USAGE

For installation and usage we refer to the manual.

## CITING DaMaSCUS

If you decide to use this code, please cite:

xxx

as well as the original publication:

xxx

## AUTHORS & CONTACT

The authors of DaMaSCUS are Timon Emken and Chris Kouvaris.

For questions, bug reports or other suggestions please contact Timon Emken (emken@cp3.sdu.dk).


## LICENCE

This project is licensed under the MIT License - see the LICENSE file.
