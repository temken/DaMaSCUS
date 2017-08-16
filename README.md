# DaMaSCUS

Dark Matter Simulation Code for Underground Scatterings

<a href="http://ascl.net/1706.003"><img src="https://img.shields.io/badge/ascl-1706.003-blue.svg?colorB=262255" alt="ascl:1706.003" /></a>

DaMaSCUS Version 1.0 06/06/2016

<img src="https://cloud.githubusercontent.com/assets/29034913/26834962/4092f75c-4ad7-11e7-86db-a359734ea2ef.png" width="425">

## GENERAL NOTES

For the underlying physics we refer to the paper [arXiv:1706.02249](https://arxiv.org/abs/1706.02249).

- DaMaSCUS is a MC simulator of dark matter particles as they move through the Earth and scatter on terrestrial nuclei. 
- It allows to compute the local distortions of the DM phase space caused by collisions with nuclei. 
- The thusly distorted distribution functions and densities are used to give precise estimates of time-dependent signal rates for direct detection experiments and diurnal modulations.
- A full, realistic model of the Earth is implemented as well as the Earth's time-dependent velocity and orientation in the galactic frame.
- DaMaSCUS is written in C++ and fully parallelized (openMPI).

## INSTALLATION AND USAGE

For installation and usage we refer to the documentation.

## Updates
The code will be updated continuously. Here I list only important updates and major bug fixes 
- `16.08.2017`: Release of v1.0.2: Re-structuring of the folder structure and new documentation.
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
