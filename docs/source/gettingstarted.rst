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

in your console or terminal. You will now find the following folders in your destination directory:

* bin
* build
* data
* include
* plots
* results
* src

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
	CXXFLAGS := -Wall -std=c++11  -O2 
	LIB := -lconfig++
	INC := -I include
	(...)

to your local settings. Next to install DaMaSCUS and compile the code simply run

.. code-block:: bash

	$ make

in your terminal. Alternatively you can also run

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