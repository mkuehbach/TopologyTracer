**Output data format specification**
====================================
The TopologyTracer generates a number of binary results file which can be visualized and processed with Matlab. These files do not include the simulated or scaled real physical time. **That one has to be supported manually as a part of the Matlab analysis!**

In what follows, a 64-bit C/C++ double precision floating point value reads as **d**, a 32-bit integer as **i**, and a 32-bit unsigned integer as **ui**. The byte order is read and written assuming LittleEndian. MPI I/O views by default all files with displacement = 0, etype = MPI_TYPE and filetype = MPI_BYTE. The TopologyTracer follows this default. All files encode as 2D matrices of an elemental datatype (either d,i, or ui, or an MPI struct) as a number of nc columns of a contiguous block (row) of size nr each. All means are arithmetic mean values. 

**The results filenames are speaking with NR = nr and NC = nc!**

**Mind that in all what follows the statement all grains means explicitly all grains filtered for the analysis. This number is not necessarily equal to the number of grains in the population, because if the simulation is performed with open boundary conditions the grains touching the domain boundary have been expelled from the analysis!**
In general the TopologyTracer assumes that the grain with ID 0 is a special grain, namely a proxy ID for the RVE domain itself. Equivalently, we can state that the TopologyTracer utilizes Fortran ID labeling for grains and faces.

analyze_vol_quantiles()
^^^^^^^^^^^^^^^^^^^^^^^
Writes two files **VolQuantiles** which specifies the n = STATISTICS_VOL_NQUANTILES (currently 100) (1/n,2/n,...,n/n) quantile of the distribution of grain size (area (2D), volume (3D)). Specifically, nr = n and nc = (SnapshotLast - SnapshotFirst)/SnapshotOffset + 1.
The elemental datatype is d.

The file **VolMeta** provides descriptive statistics with nr = 1 and nc as above. The elemental type is a contiguous MPI struct of five d. They encode

* The total area covered by all grains considered
* The mean size
* The variance
* The total number of grains in the snapshot
* The total number of grains considered (for mean and var)


analyze_vol_forward()
^^^^^^^^^^^^^^^^^^^^^
Writes four 2D matrices, all with the same layout. Except for **FW.NF**, whose elemental datatype is ui, all others have d.

* **FW.NF** is the number of neighbors
* **FW.VOL** is the size of the grains (micron^Dimensionality)
* **FW.HAGBFRAC** gives fraction of boundary (length (2D), area (3D) with disorientation of target and neighbor >= HAGBDetectionThreshold
* **FW.MOBDSEE** gives area-weighted product of mobility * (stored elastic energy nbor - stored elastic energy target) (m/s)

Values for no longer existing grains get filled up with zeros!


analyze_vol_backward()
^^^^^^^^^^^^^^^^^^^^^^
Monitor for all targets (grains which survived the coarsening step up to SnapshotLast) backwards in time for the evolution of their faces, size, HAGB fraction, and instantaneous biased velocity. **Even though tracking is performed backwards, the organization of the matrix is positive in time, i.e. higher column vector indices mean later in time.**

The layout and content is the same as for analyze_vol_forward(). The only differences are the file name suffixes.

* **BK.NF**
* **BK.VOL**
* **BK.HAGBFRAC**
* **BK.MOBDSEE**

analyze_sizegain_vs_bk()
^^^^^^^^^^^^^^^^^^^^^^^^
Idea: for all survivors identify in each snapshot how their size is relative to the mean of the matrix.
The matrix constitutes of all grains excluding the survivors. Writes one 2D matrix of elemental type d. 
Each first entry in each row describes the number of matrix grains remaining, while the second entry the mean size of the matrix grains.

* **SIZEGAINBK** 


analyze_approx_rxfraction()
^^^^^^^^^^^^^^^^^^^^^^^^^^^
Calculate approximate recrystallized area/volume fraction assuming the following.
Matrix grains constitute all grains. The recrystallized grains are all those which survive the coarsening, i.e. all remaining up to SnapshotLast.
A 2D matrix with nr = 1 is written. The elemental datatype is a contiguous MPI_struct of d. These detail

* The total size of all grains considered in the analysis
* The total size covered by the targets, i.e. the considered as recrystallized
* The recrystallized coverage as TotalSizeTargets/TotalSizeAllGrains
* The total number of all grains considered in the analysis but were still detectable in the snapshot
* The total number of rxgrains


analyze_modf()
^^^^^^^^^^^^^^
Writes a 2D matrix of doubles (d) which encodes the disorientation distribution function (MDF) with equiangular adjustable binning, via nr = largest positive integer ( MaxDisoriAngle/DisoriAngleBinWidth) and nc = 1 + ((SnapshotLast-SnapshotFirst)/SnapshotOffset+1). The first row of the matrix gives the bin ends in degrees. Thereafter, the MDF for each snapshot follows.


analyze_see()
^^^^^^^^^^^^^
Same layout as analyze_modf. First row gives bin ends. Thereafter, each snapshot.

analyze_gsd()
^^^^^^^^^^^^^
Equal in mentality and data output structure to SEE and MODF.

analyze_abnormal_graingrowth()
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Writes a 2D matrix of structs containing 13 doubles (d) each which encode the following

The total size of all grains in the list of

* survivors
* matrix
* matrix exclusive survivors

The equivalent radius of

* survivors
* matrix
* matrix exclusive survivors

The remaining number of grains in 

* survivors
* matrix
* matrix exclusive survivors

The number of abnormal large grains, i.e. whose equivalent radius is larger than 3.0 the average of the individual populations

* matrix exclusive survivors abnormal when testing against mean of matrix
* matrix exclusive survivors abnormal when testing against mean of matrix exclusive survivors
* survivors abnormal when testing against mean of matrix
* survivors abnormal when testing against mean of matrix exclusive survivors

Mind that when approaching SnapshotLast the population of matrix grains ceases because this understood as the matrix into which the survivors grow.



**Currently undocumented**
^^^^^^^^^^^^^^^^^^^^^^^^^^
* analyze_grainsize_quantiles()
* analyze_drivingforce_see()
* analyze_topologydifference_forward()
* analyze_classical_nucmodels()
* analyze_knn_naive()
