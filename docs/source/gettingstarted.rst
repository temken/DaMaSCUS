Getting started
===============

------------
Requirements
------------

^^^^^^^^^^^^
Dependencies
^^^^^^^^^^^^

These are the dependencies of DaMaSCUS:

* **libconfig**: To handle the input configuration files we use the \textit{libconfig} library.

	`http://www.hyperrealm.com/libconfig/ <http://www.hyperrealm.com/libconfig/>`_

* **Eigen**: DaMaSCUS relies heavily on this linear algebra C++ library.

	`http://eigen.tuxfamily.org/ <http://eigen.tuxfamily.org/>`_

* **openMPI**: For the parallelization we implemented our code using the open *Message Passing Interface*.
	
	`https://www.open-mpi.org <https://www.open-mpi.org>`_

------------
Download
------------

The DaMaSCUS code is available at

	`https://github.com/temken/DaMaSCUS/ <https://github.com/temken/DaMaSCUS/>`_ .

To download it via git simply run

.. code-block:: bash

	$ git clone https://github.com/temken/DaMaSCUS/

in your console or terminal. 

^^^^^^^^^^^^^^^^
Folder Structure
^^^^^^^^^^^^^^^^

You will now find the following folders in your destination directory:

* **/bin/**: After successful compilation this folder contains two executables as well as the configuration file.
* **/build/**: This folder contains all object files. Both the object files and the executables in **/bin/** are deleted via 

.. code-block:: bash

	$ make clean

* **/data/**: Once a simulation run is performed, the generated data will be stored here
* **/include/**: The DaMaSCUS header files are stored here. Necessary 3rd party libraries can also be placed here.
* **/plots/**: To visualize the results, created by the analysis module we include the small Mathematica package *DaMaSCUStoolbox* and an example notebook creating and saving plots.
* **/results/**: The analysis module saves its results and histograms here.
* **/src/**: All the source code files of the two DaMaSCUS modules can be found here.

------------
Installation
------------

DaMaSCUS consists of two more ore less independent modules:

1. **DaMaSCUS-Simulator**: Simulates the dark matter trajectories and genererates the raw data.
2. **DaMaSCUS-Analyzer**: Analyzes the raw data and calculates e.g. velocity histograms or detection rates.

The code is compiled using the Makefile. You might have to adjust the first lines 

.. code-block:: guess

	#Compiler and compiler flags
	CXX := mpic++
	CXXFLAGS := -Wall -std=c++11 
	LIB := -lconfig++
	INC := -I include
	(...)

to your local settings. Next to install DaMaSCUS and compile the code simply run

.. code-block:: bash

	$ make

from the root directory in your terminal. Alternatively you can also run

.. code-block:: bash

	$ make simulator

or

.. code-block:: bash

	$ make analyzer

to just compile one of the modules. 

Finally, running

.. code-block:: bash

	$ make clean

deletes all object files and executables.