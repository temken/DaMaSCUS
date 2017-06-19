# DaMaSCUS

Dark Matter Simulation Code for Underground Scatterings

<a href="http://ascl.net/1706.003"><img src="https://img.shields.io/badge/ascl-1706.003-blue.svg?colorB=262255" alt="ascl:1706.003" /></a>

DaMaSCUS Version 1.0 06/06/2016

<img src="https://cloud.githubusercontent.com/assets/29034913/26834962/4092f75c-4ad7-11e7-86db-a359734ea2ef.png" width="425">

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

## Updates
The code will be updated continuously. Here I list only important updates and major bug fixes 
- `15.06.2017`: Release of v1.0.1: major bug fix concerning mostly very high cross-sections.
- `06.06.2017`: Release of v1.0


## CITING DaMaSCUS

If you decide to use this code, please cite

>Emken, T. & Kouvaris, C., 2017, DaMaSCUS, Astrophysics Source Code Library, record [ascl:1611.012](http://ascl.net/code/v/1702)

as well as the original publication,

>Emken, T. & Kouvaris, C., DaMaSCUS: The Impact of Underground Scatterings on Direct Detection of Light Dark Matter, (2017), [arXiv:1706.02249](https://arxiv.org/abs/1706.02249).

## AUTHORS & CONTACT

The authors of DaMaSCUS are Timon Emken and Chris Kouvaris.

For questions, bug reports or other suggestions please contact Timon Emken (emken@cp3.sdu.dk).


## LICENCE

This project is licensed under the MIT License - see the LICENSE file.
