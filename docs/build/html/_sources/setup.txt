**Setup**
=========

Which prerequisites are necessary?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The compilation of the program utilizes open source libraries as well as standard Linux tools.

* Check your Linux installation for a working installation of a **C/C++ build system** including **cmake** and **make**.
* Check your **C and C++ compiler**. We utilized successfully the Intel (v14.0 and newer) and the GNU (v4.8 and newer) compiler.
* You need a working installation of an MPI_ (Message Passing Interface) API library to be able to compile the program. 
* We utilized successfully both the OpenMPI (v1.10.4 and newer) and the IntelMPI (2017.0.12) implementation.
* The minimum threading support level of the MPI implementation required is **MPI_THREAD_FUNNELED**. 
* You need a working installation of the **Boost C++ libraries**. Details on where to get it and how to install are found here_.
* The shipped-with default version of Boost of at least Ubuntu 16.04 provided to me all functionality required.
* You need the **Eigen libraries** that come along with the TopologyTracer package. Further information to this library is provided_.


.. _MPI: http://www.mcs.anl.gov/research/projects/mpi/
.. _here: http://www.boost.org 
.. _provided: http://eigen.tuxfamily.org/
 
Which prerequisites are optional?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* It is well-known that the general purpose standard malloc(3) memory allocator class implementation does not assure that in particular many small allocations become locality-aware placed. Hence, the performance of the TopologyTracer can be improved by linking against an alternative memory allocator class, such as **Jason Evans jemalloc**. A documentation of how to obtain the library and how to compile it can be found online_.
* In the future, an OpenMP-parallelization and HDF5 extension of the TopologyTracer is planned. So far however HDF5 is not required to run the TopologyTracer.

.. _online: http://www.canonware.com/jemalloc/

Where to get source code from?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
TopologyTracer is available for Linux only. It is free software.

* Download the source from its git repository **https://github.com/mkuehbach/TopologyTracer**
* Eventually unpack the repository such that finally the following ends up in a single folder.
* This folder from now on will be called the **root** directory. You can give it any Linux-conformant name.
* Make sure in this root there is a **src** subdirectory with the cpp and the h source code files, 
* a **build** directory for storing the executable and an exemplary **XML control file**.
* Additionally, check whether there is a **CMakeLists.txt** file in the root folder.
* Next, utilize the top section of this CMakeList.txt file to switch on and off the compiler (GNU or Intel) 
* Next, open a console and dive into the **build** directory.
* If now its the first time you compile the TopologyTracer type **cmake ..**. 
* This inspects your system and generates a customized makefile for you.
* Next use this makefile by typing **make** to compile the program.
* Upon success you should now have a **binary** within build called **topotracer2d3d**

Placement of files
^^^^^^^^^^^^^^^^^^
The resulting executable expects the XML control file always in its current location folder!. Relative indexing is utilized. Other than that restriction, the executable can be renamed and relocated. The latter enables the scripting of TopologyTracer job queues and executing batch jobs.

Optimization
^^^^^^^^^^^^
If desired, adjust the level of compiler optimization via the OPTLEVEL variable in the CMakeLists.txt upper section. OPTLEVEL "-O0" means no optimization and should be utilized for debugging purposes only, while "-O3" is the maximum and recommended level for production tasks. Improvements between the two extremes vary between a factor of 2 - 5 faster with maximum optimization compared to without.

Troubleshooting?!
^^^^^^^^^^^^^^^^^
If in between the compilation process unrecoverable errors occur, attempt first a **make clean** command. 
If this does not help: Delete everything in the build folder except for the **TopologyTracer2D3D_Parameter.xml** control file and start over with **cmake ..**.


Utilizing the GNU compiler? If the __int64 definition is unknown to the preprocessor and the code therefore does not complains, does not compile: define it as typedef long __int64 in the **TopologyTracer_Topology.cpp**.