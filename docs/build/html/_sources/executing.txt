**Program execution**
=====================

MPI only
^^^^^^^^	
All set? Excellent! Then, the tracing is executed via one command line call::

   mpiexec -n <nprocesses> <topotracer> <simid> <TopologyTracer2D3D_Parameter.xml> <optional: arguments>
   
**Please note that the angle brackets must not be typed into the command line as they only mark the input parameter!**
   
In its current version the following input arguments are required:

* <nprocesses> How many MPI processes to utilize?
* <topotracer> The name of the executable.
* <simid> JobID, a positive integer to distinguish the datamining results from runs with other settings but the same raw data. 
* <TopologyTracer2D3D_Parameter.xml> a properly formatted XML control file. The name can be changed as long as the file remains a properly formatted XML file.

**Be careful: if the <simid> value is set to the same value during subsequent runs in the same folder, data will be overwritten without prompting!**

Below is a typical example call which executes the program with 10 MPI processes, reads from MySpecialSettings what should be done and tags all output with a consistent run ID, 1000 in this case::

   mpiexec -n 10 topotracer2d3d 1000 MySpecialSettings.xml

MPI/OpenMP
^^^^^^^^^^
Core functionalities within the analyses modes 4, 5, and 6 are additionally OpenMP thread parallelized.
For this, the OpenMP environment variable OMP_NUM_THREADS must be set to the desired number of threads spawned at most per MPI process.
Otherwise the program call is as above::

   export OMP_NUM_THREADS=10

When utilizing exemplarily the Intel compiler it is usually worthwhile to reconsider the thread placement on the computing cores.
One can guide this operating system decision by for instance the KMP AFFINITY_ environment.

 .. _AFFINITY: https://software.intel.com/en-us/node/522691