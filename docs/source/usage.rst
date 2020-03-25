==============
Using DaMaSCUS
==============

---------
Work Flow
---------

To simulate dark matter trajectories and analyze the generated data DaMaSCUS has a clear work-flow.

1. Adjust your input parameter (such as the dark matter mass) inside the configuration file and asign a simulation ID to identify this simulation run.
2. Run the simulation module to generate the raw data.
3. Run the analysis module to process the data.
4. The results can be plotted e.g. with the included Mathematica notebook.

^^^^^^^^^^^^^^^^^^^^^
1. Configuration File
^^^^^^^^^^^^^^^^^^^^^

You can find an example configuration file at

.. code-block:: bash

	/bin/DaMaSCUS.cfg 

Here you adjust all the input parameter for the next DaMaSCUS run. We go through it block by block.

.. code-block:: guess

	//DaMaSCUS Configuration File

	//Simulation input parameter
			simID		=	"exampleID";	//MC Simulation ID
			initialruns	=	100000000L;	//Number of particles in the initial MC run
			samplesize	=	10000;		//velocity sample size per isodetection ring
			vcutoff		=	1.0;		//velocity cutoff in cm/sec
			rings		=	36;				//number of isodetection rings

First you assign the simulation run a unique identifying ID. You also decide the number of particles you want to simulate in the initial run without scatterings (initialruns) and how many data points you need in each isodetection ring (samplesize). The velocity cut-off, below which a trajectory simulation is aborted.

New in version 1.1: The number of isodetection rings is now flexible and can be set in the configuration file (rings).

.. warning::

	The "L" behind the value for initialruns is necessary to denote that it is a **long int**.

.. code-block:: guess

	//Simulation Time:
			date		=	[15,02,2016];	//Date [dd,mm,yyyy]
			time		=	[0,0,0];	//Universal time [h,m,s]

Next you fix the simulation date and time, which is mostly used to determine the Earth's velocity in the galactic frame.

.. code-block:: guess

	//Dark Matter Data
		//Particle data
			mass		=	500.0;		//in MeV
			sigma 		=	1.0;		//in pb 
			formfactor	=	"None";		//Options: "None", "HelmApproximation"
		//DM Halo 
			halomodel	=	"SHM";		//Options: Standard Halo Model "SHM",...
			rho			=	0.3;	//DM halo energy density in GeV/cm^3

These are the essential input parameter. Here you fix the DM mass, its cross-section with nucleons and whether or not to use a Helm form factor. For light DM it is recommended to deactivate it, since this speeds up the simulation considerably. 

So far only the standard halo model is implemented, for which we can set the DM energy density. 

.. code-block:: guess

	//Detector depth:
			depth		=	1400.0;		//in meter

The depth of the isodetection rings, i.e. the underground depth of the detector of interest is determined here.

.. code-block:: guess

	//Analysis parameter
			experiment	=	"CRESST-II";	//Options: "LUX" for heavy DM,"CRESST-II" for light DM and "None"

And finally we decide the type of data analysis:
1.	Set "LUX" for a LUX-type detector. Use this option for heavy dark matter.
2. 	Use "CRESST-II" for a CRESST-type detector, which is sensitive to DM masses down to 500 MeV.
3. Set "None", if you are e.g. just interested in the resulting speed distributions across the globe.

.. warning::

   Note that the configuration file can be sensitive to the input parameter type. For example it might complain if an input parameter of type **double** is given as "1" instead of "1.0".

^^^^^^^^^^^^^^^^^^^^^^^^^
2. Running the simulation
^^^^^^^^^^^^^^^^^^^^^^^^^

After setting the input parameter and assigning a unique simulation ID we can start the MC simulation from the **/bin/** directory. To start run

.. code-block:: bash

	$ mpirun -n N DaMaSCUS-Simulator input.cfg

where *N* is the number of MPI processes and *input.cfg* is your configuration file.

After a successful run your terminal should show something like

.. code-block:: guess

	$ mpirun -n 4 DaMaSCUS-Simulator DaMaSCUS.cfg

	##############################
	DaMaSCUSv1.0 - Simulation

	Starting Time: Wed Aug 16 11:34:08 2017
	Simulation ID: exampleID

	Creating logfile.
	Creating folder for velocity data.
	Start initial MC simulation run without DM scatterings.
		Initial run finished	(1 s).

	Start main MC simulation run with scatterings.
		Main MC run finished	(4 s).

	Processing Time:	5.89347s (0:0:5:893).
	##############################

A copy of the used configuration file is stored in the **/data/** directory together with the raw data. In addition a logfile, which documents important input and output parameter is created in the **/results/** folder.

^^^^^^^^^^^^^^^^^^^
3. Analyze the data
^^^^^^^^^^^^^^^^^^^

Next we can analyze the generated data by running

.. code-block:: bash

	$ mpirun -n N DaMaSCUS-Analyzer SimulationID

in your terminal from the **/bin/** directory. Again *N* is the number of MPI processes. The analysis type is set inside the config file **/data/SimulationID.cfg** and can be adjusted after the simulation has finished. The terminal output of a successful analysis looks like

.. code-block:: guess

	$ mpirun -n 4 DaMaSCUS-Analyzer exampleID

	##############################
	DaMaSCUSv1.0 - Data Analysis

	Starting Time:	Wed Aug 16 11:49:05 2017
	Simulation ID:	exampleID
	Experiment:	None

	Creating folder for histograms.
	Done.

	Creating temporary files.
	Reading in local DM densities.
	Broadcast local DM densities to all MPI processes.
	Start data analysis...

	MPI rank	Isodetection ring	Local Progress	Computing time[s]	Residual time estimate[s]
	(*some infos about the progress*)

	Data analysis complete.
	Creating ASCII output.
	Delete temporary files and finish.

	Processing Time:	1.02264s (0:0:1:22).
	##############################


^^^^^^^^^^^^^^^^^^^
4. Plot the results
^^^^^^^^^^^^^^^^^^^

After the both modules have finished their computations you can use the included Mathematica notebook **/plots/plots.nb** to create and save plots of e.g. the speed histograms or the event rate modulation.