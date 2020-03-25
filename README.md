[![Documentation Status](https://readthedocs.org/projects/damascus/badge/?version=latest)](http://damascus.readthedocs.io/en/latest/?badge=latest)
[![Build Status](https://travis-ci.com/temken/DaMaSCUS.svg?branch=master)](https://travis-ci.com/temken/DaMaSCUS)
[![codecov](https://codecov.io/gh/temken/DaMaSCUS/branch/master/graph/badge.svg)](https://codecov.io/gh/temken/DaMaSCUS)
  

# DaMaSCUS

Dark Matter Simulation Code for Underground Scatterings

<a href="http://ascl.net/1706.003"><img src="https://img.shields.io/badge/ascl-1706.003-blue.svg?colorB=262255" alt="ascl:1706.003" /></a>
[![arXiv](https://img.shields.io/badge/arXiv-1706.02249-B31B1B.svg)](https://arxiv.org/abs/1706.02249)

DaMaSCUS Version 1.1 25/03/2020

<img src="https://cloud.githubusercontent.com/assets/29034913/26834962/4092f75c-4ad7-11e7-86db-a359734ea2ef.png" width="425">

## GENERAL NOTES

For the underlying physics we refer to the paper [![arXiv](https://img.shields.io/badge/arXiv-1706.02249-B31B1B.svg)](https://arxiv.org/abs/1706.02249).

- DaMaSCUS is a MC simulator of dark matter particles as they move through the Earth and scatter on terrestrial nuclei. 
- It allows to compute the local distortions of the DM phase space caused by collisions with nuclei. 
- The thusly distorted distribution functions and densities are used to give precise estimates of time-dependent signal rates for direct detection experiments and diurnal modulations.
- A full, realistic model of the Earth is implemented as well as the Earth's time-dependent velocity and orientation in the galactic frame.
- DaMaSCUS is written in C++ and fully parallelized (openMPI).

## INSTALLATION AND USAGE

For installation and usage we refer to the [documentation](http://damascus.readthedocs.io/en/latest/).

## Updates
The code will be updated continuously. Here I list only important updates and major bug fixes 
- `25.03.2020`: Release of v1.1
- `16.08.2017`: Release of v1.0.2: Re-structuring of the folder structure and new documentation.
- `15.06.2017`: Release of v1.0.1: major bug fix concerning mostly very high cross-sections.
- `06.06.2017`: Release of v1.0

## CITING DaMaSCUS

If you decide to use this code, please cite

>Emken, T. & Kouvaris, C., 2017, DaMaSCUS, Astrophysics Source Code Library, record [[ascl:1706.003]](https://ascl.net/1706.003)

as well as the original publication,

>Emken, T. & Kouvaris, C., DaMaSCUS: The Impact of Underground Scatterings on Direct Detection of Light Dark Matter,
[![JCAP](https://img.shields.io/badge/JCAP-1710(2017)no.10,031-255773.svg)](http://iopscience.iop.org/article/10.1088/1475-7516/2017/10/031/meta), 
[![arXiv](https://img.shields.io/badge/arXiv-1706.02249-B31B1B.svg)](https://arxiv.org/abs/1706.02249)

In addition, version 1.1 of the code is archived under

> [[DOI:10.5281/zenodo.XXXXX]](https://doi.org/10.5281/zenodo.XXXXX)

## AUTHORS & CONTACT

The authors of DaMaSCUS are Timon Emken and Chris Kouvaris.

For questions, bug reports or other suggestions please contact Timon Emken (emken@chalmers.se).


## LICENSE

This project is licensed under the MIT License - see the LICENSE file.
