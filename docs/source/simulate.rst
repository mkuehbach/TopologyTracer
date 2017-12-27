What is TopologyTracer input?
=============================
A temporal sequence of microstructure dynamics snapshots. Currently, the TopologyTracer requests each of these snapshots as to be pairs of two binary files: a **Faces_*.bin** and a **Textures_*.bin** file. The asterisk placeholder (*) stands for a set of integer numbers separated by a constant offset, i.e. (1,2,3,4) or (10, 20, 30, 40) for instance, on an interval [SnapshotFirst,SnapshotOffset,SnapshotLast].
As the format of these files is generic, it is possible to feed either in-situ microstructural data from experiments or check pointing data from simulations.
Currently, the TopologyTracer is implemented to work with input from the GraGLeS level-set code (https://github.com/GraGLeS).
Feel free contacting me when there is interest in adding input routines for other tools. I am aware that conventions within our community differ so I am highly interested in offering such alternative entry to use my code.


Binary file format layout
=========================
The TopologyTracer utilizes binary data in order to provide accuracy, to minimize data storage costs, and to reduce I/O time via a one-file-per-process concept. During I/O operations all binary files are interpreted assuming

* LittleEndian byte order
* The default MPI view of the file is as of a linear byte stream
* Namely, displacement = 0, etype = MPI_BYTE, filetype = MPI_BYTE

In what follows a 64-bit C/C++ double (precision) floating point value reads as **d**, a 32-bit integer as **i**, a 32-bit unsigned integer as **ui**. 

Grain boundary face data via Faces files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The **Faces_*.bin** file is headerless and contains all disjoint grain boundary faces. It reads as a contiguous block of as many MPI structs as faces exist. Each struct specifies the segment length and the two disjoint integer IDs A, B (with A > B >= 0). These IDs reference two grains that share the boundary face. As such, each boundary is stored only once. Currently, its description occupies 16B on a 64-bit machine. A single entry reads as follows:: 

	d ui ui
	
Grain meta data via Texture files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The **Texture_*.bin** file is headerless as well. It contains the properties of all grains at the timestep, i.e. size, barycenter position to name but a few. The binary file encodes a contiguous block of as many MPI structs as grains remain. Each MPI struct encodes the properties of a disjoint grain. **The struct layout differs for 2D and 3D data in an effort to reduce the amount of data!** 

In **2D** the order is as follows:
Area, Perimeter, GBEnergy, Stored elastic energy, Bunge-Euler phi1, Bunge-Euler Psi, Bunge-Euler phi2, Barycenter X, Barycenter Y, 
GrainID, First-order neighbor count, Contact with boundary yes/or no, padding.

In **3D** the order is as follows:
Area, Perimeter, GBEnergy, Stored elastic energy, Bunge-Euler phi1, Bunge-Euler Psi, Bunge-Euler phi2, Barycenter X, Barycenter Y, Barycenter Z,
GrainID, First-order neighbor count, Contact with boundary yes/or no, padding.

The MPI struct encodes these pieces of information as follows::

	d d d d d d d d d ui ui ui ui
	
	d d d d d d d d d d ui ui ui ui

Consequently, a grain occupies 88B in 2D and 96B in 3D, respectively on a 64-bit architecture when utilizing the default byte order.

Illustrative example of expected datavolume
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
For a 3D polycrystal with 500,000 grains spanning a boundary network of 3,500,000 individual faces a single snapshot requires ((5.0e5 * 96B) + (3.5e6 * 16B)) / 1024  / 1024 KB/B MB/KB = 99.2MB disk space. For comparison, a practical peak bandwidth of a Lustre parallel file system for a one-file-per-process I/O scenario with a few MPI processes ranges typically between 1000-15000 MB/s.
	