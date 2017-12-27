**XML Control File Settings**
=============================
The entire datamining is controlled by only one control parameter settings file --- the **TopologyTracer2D3D_Parameter.xml** file. All angular values are expected in degrees!

Analysis mode
^^^^^^^^^^^^^
| **AnalysisMode**
|	Determines which execution model and fundamental tasks should be executed.
|	Option (1) (the default) is for tracing data in parallel. It distributes an ensemble of Face/Texture data on the MPI processes.
|	Option (2) tracking sequentially by identifying first all grains close to boundaries, thereafter loading successively the files for processing. 
|   Make sure to populate GrainBoundaryContact with a file to be generate/or existent that stores the deadtimes of the grains.
|   You should switch on **ProbeBoundaryContact** to generate such file and switch it off once done and the file exists.
|	Option (3) reads sequentially a TargetGrainIDs file and a KNNFile to correlate the evolution of a grain with a prognosis 
	of their evolution based on the k-th neighborhood in the initial structure (SnapshotFirst). The functionality is UNDOCUMENTED.
|	Option (4) reads binary GBContourPoints_<fid>.bin files to compute unbiased quantitative metrics for each nucleus to its long-range environment in 2d. The functionality is UNDOCUMENTED.
|	Option (5) is the currently only MonteCarlo discretely implemented extension of this in 3D. The functionality is UNDOCUMENTED.
|   Option (6) reads binary GBContourPoints_<fid>.bin files to approximate for each grain the capillary driving force by evaluating the turning angles at the piecewise segmented grain boundary contour. The functionality works at the moment only for 2d contours.
|   Option (7) is UNDOCUMENTED.

Fundamental settings
^^^^^^^^^^^^^^^^^^^^
| **SnapshotFirst**
|	Specifies the integer index of the first dataset to analyze.
| **SnapshotOffset**
|	Specifies the constant integer increment of datasets desired, 
|	i.e. tracing is performed on the interval [First, Last] in increments of Offset. 
| **SnapshotLast**
| 	Specifies the integer index of the last dataset to analyze.
|	All indices are integer, >= 0. First and Last can be set to the same value to analyze only one dataset namely the SnapshotFirst.
| **Dimensionality**
|	Is the dataset describing 2d (2) or 3d (3) simulations?
| **LargestGrainID**
|	The user has to specify the maximum grain ID in the entire dataset.
| **ProbeBoundaryContact**
|   You should switch this on only when in mode 2 to analyze once the dead times of the grains into a file named **GrainBoundaryContact**.

Physical properties
^^^^^^^^^^^^^^^^^^^
| **HAGBMobility**
|   Specifies the **maximum** mobility of all grain boundaries (m^4/Js).
| **HAGBEnergy**
|	Specifies the **maximum** specific energy of all grain boundaries (J/m^2).
| **DislocEnPerM**
|	Specifies the product 0.5Gbb,i.e the dislocation line energy per unit length dislocation (J/m).
| **PhysicalDomainEdgeLength**
|	Specifies the real distance which the simulation domain with InitialDomainEdgeLength represents (m). 
| **InitialDomainEdgeLength**
|	Specifies the integer edge length of the simulation domain to the beginning of the simulation (square in 2D, cube in 3D).

| **Mind that all these settings have to meet those in the coarsening simulation!**
| **Especially, LargestGrainID is required as the largest unsigned integer grainID**
| This is not necessarily the total number of grains in the simulation!

Datamining operations
^^^^^^^^^^^^^^^^^^^^^
The following options define which analyses should be executed (1) or not (0).

| **AnalyzeGrainSizeQuantiles**
|	Computes descriptive means of the grain size (area, vol) and the volume-average quantiles for each Snapshot.
| **AnalyzeSEE**
|	Approximates the CDF of the stored elastic energy in the grains in an area/volume averaged manner.
| **AnalyzeMODF**
|	Approximates the CDF of the MODF (disorientation angles in boundary network (face segment length/area averaged)
|	**This is not the MacKenzie plot (probability density of occurrence over disorientation angle) but its integrated form!**
|	Only considers boundaries between grains which both (instantaneously) do not make simulation domain boundary contact.
|	*GraGLeS outputs also an MODF, via an file MODF_* ASCII file. Numerical differences exist.*
| **AnalyzeGSD**
|	Approximated a histogram of the population arithmetic average normalized distribution of grain sizes (area, volume).
| **AnalyzeMaxSizeGainForward**
|	Computes for all grains the maximum size (area/volume) they obtain in all Snapshots.	
| **AnalyzeApproxRXFraction**
|	Determines two populations of grains: one with all grains which never touched the simulation domain boundaries 
|	and another from those remaining in SnapshotLast (targets). 
|	Therefrom it is computed the instantaneous coverage ratio (area, volume) of these targets with the population.
| **AnalyzeTrajectoriesBackward**
|	Extracts a list of those grains in the last timestep (SnapshotLast) which not made contact with the simulation domain boundaries.
|	These grains are then followed backwards in time. This is the suggested mode of operation to follow the evolution of persistent grains.
| **AnalyzeTrajectoriesForward**
|	Extracts first all grains which never during their existence, i.e. (size > 0) made contact with simulation domain boundaries if any. 
|	Then, their properties are tracked forward in time. If a grain disappeared earlier than SnapshotLast, 
|	all succeeding result values for within the column vector are set to zero.
|	**Mind that the resulting output can be very large as a total of 10 million grains initially tracked for 1000 snapshots will result in a matrix of 10 billion doubles!**
| **AnalyzeAbnormalGrains**
|	Determines three populations of grains, all of which never touched simulation domain boundaries. 
|	One from SnapshotLast (the possible AGG candidate grains).
|	One from all grains.
|	One of all the latter grains but excluding the candidates (the matrix). 
|	It is evaluated for each time step the (spherical equivalent radius ratio) of the candidate grains to the matrix.
| **AnalyzeKNN**
|	Determines the k-th order (nearest) neighbors to a target of the grains. Each grain of the dataset is inspected independently
| 	and resulting in a successful detection only if none of the neighbors made contact with domain boundaries to prevent bias. 
|	During this neighborhood analysis also the region boundary length is computed.
|	Between the 0-th and the 1-th order neighbor this is the grain boundary perimeter of the target itself. 
|	Also the fraction of high-angle grain boundaries is computed. WORKFLOW IS UNDOCUMENTED!


The following functions are currently UNDOCUMENTED

| **AnalyzeDrivingForceSEE**
| **AnalyzeMeanDrivingForceSEEForward**
| **AnalyzeSizeGainVsMatrixBackward**
| **AnalyzeSizeGainVsMatrixForward**
| **AnalyzeClassicalNucModels**
| **AnalyzeTopologyDifferenceForward**


Datamining operation specific settings
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Datamining is performed on the timestep interval [SnapshotFirst, SnapshotLast] in increments of SnapshotOffset. Grains with contact to the simulation domain boundaries (special grain ID 0) are expelled from the analyses in every case by default. The code allows for the implementation of utilizing results from simulations that were conducted under periodic boundary conditions to enable higher statistical significance. The latter is of particular interest for 3D simulations and **in particular when quantifying the higher-order neighbors**. The results are either stored in ASCII files (***.csv**) or as binary matrices 2D (***.bin**). The latter have no header but a **speaking filename** (<whatever>.F.<SnapshotFirst>.O.<SnapshotOffset>.L.<SnapshotLast>.NC.<ncols>.NR<nrows>.bin) which
specifies the number of columns (NC) and rows (NR), respectively. In general, columns are column vectors, each of which encoding all grains at one time step, while rows are row vectors which specify properties of one grain along all time steps.
With these definitions the binary results file read as a 2D matrix (grains-time) in implicit form by aligning entire column vectors according to their row ID contiguously in memory.

| **AnalyzeTrajectoriesForwardMode**
|   Modifies AnalyzeTrajectoriesForward, when unset (0) all grains are tracked, when set (1) the TopologyTracer requires TargetGrainIDs to be set in order to perform tracking only on these included IDs. Then an additional file needs to be supplied (**TargetGrainIDS**).
| **MaximumNumberOfKShells**
|	Specifies the maximum order of (long-range) neighbors during the analysis AnalyzeKNN. 
| **HAGBDetectionThreshold**
|	When is a grain boundary considered as to be of high-angle character? 
|	Sensible default in accordance with Read-Shockley theory is 15 degrees.
| **DisoriAngleBinningMin**
| **DisoriAngleBinningMax**
| **DisoriAngleBinningWidth**
|	The values specify the binning of the MODF.
| **StoredElasticEnergyBinningMin**
| **StoredElasticEnergyBinningMax**
| **StoredElasticEnergyBinningWidth**
|	The equivalent binning specification for the stored elastic energy density.
| **GSDBinningMin**
| **GSDBinningMax**
| **GSDBinningWidth**
| 	The equivalent binning specification for the mean-normalized size histogram.

| **PersistenceRadiusMin**
| **PersistenceRadiusIncr**
| **PersistenceRadiusMax**
|	Radius (normalized to RVE domain edge length) of the fixed size reference window about each target for which the intruding neighboring is computed. The incrementer allows to setup a basic loop to sample differently sized reference windows within [Min,Min+Incr,Min+Incr+Incr,...,Max].
| **UDSFile**
|	A uds file specifying the properties of all the grains in the synthetic structure. Is written by the microstructure generator.

| **ComputeCurvatureAlsoAtTJP**
|	Flag to overwrite the default behavior that for supporting points close to junctions the curvature is not computed into a scalar average of the grain. See supplementary material to the long-range environment reference paper for more details.
| **TranslateBinary2GNU**
|	Does what it reads at the cost of potentially significant disk space.

| **CapillaryActivityTargets**
|	Allows to reduce list of target grains to specific user-defined set of IDs rather than to operate on grains of all IDs in the input data.
| **CapillaryActivityMode**
|	Without digging into source code, leave with 0.

Auxiliary files
^^^^^^^^^^^^^^^
| **GrainBoundaryContact**
|	See explanation for analysis mode 2 and ProbeBoundaryContact.
| **TargetGrainIDs**
|   Enables to supply a headerless ASCII file of grainIDs (only positive 32-bit unsigned int IDs on interval 
|	[1,LargestGrainID] are interpreted) in order to restrict analyses to specific grains. 
|	Specifically, when AnalyzeTrajectoriesForward is chosen (1)
|	and AnalyzeTrajectoriesForwardMode set to tracking forward (1) in time the forward tracking is performed only for 
|	these targets and not the entire population to reduce the overall size of the output.
| **KNNFile**
|	enables to supply a single line header-equipped ASCII file with the layout grainID, x, y 
|	coordinate for a currently UNDOCUMENTED option.
| **GNUFile**
|	A gnu 2D grain boundary contour file. The functionality is UNDOCUMENTED and DEPRECATED.

Additional settings
^^^^^^^^^^^^^^^^^^^
| **OnlyProbeTheWorkPartitioning**
|	In mode 1 this option allows to execute only the dataset partitioning planning without reading the heavy data. 
|	This allows to plan the allocated load partitioning and to verify that all 
|	processes get snapshot data volume as equally as possible distributed.
| **LocalDatabaseMaximumSize**
|	Controls the maximum size of binary data (Faces, Texture) which a single MPI process accepts to store and process in AnalysisModes 1 and 2. 
|	At the moment the workpartitioning works as follows: i) it computes the total size of all binary file pairs, ii) then it partitions on the MPI processes.
|	These do usually **but not necessarily handle an equal number of Faces,Texture binary file pairs!** 
|	**In every case the processes neither handle incomplete snapshots nor do they share or duplicate snapshot data!**
|	Namely, either they handle as many datasets (pairs Faces/Texture) as their capacity allows. If their capacity is exceeded,
|	the next process accumulates snapshots, until either there are no more snapshots to distribute and the processing can start or the program has to exit in a controlled manner because all processes have already been exhausted and should no longer load snapshot data. Hence, the value must be chosen with care. The idea behind this concept is to enable the user to set an upper bound on how much memory the MPI processes have to supply at least to load snapshot data and thus prevent the program from flooding memory or running in above hard memory limits when running in batch queues.
|	The design implies that the more data per process are set the more processes with the higher rank IDs will idle as the partitioning
|	of snapshots proceeds in the order of the MPI process ID (rank). On the contrary, if the value is too low such as that they do not 
|	allow to partition the entire ensemble, the post-processing terminates.**
| **MPIReadBlockLength**
|	Sensible default of the MPI I/O buffer size.
| **MaxIDRange**
|	Sensible default for how many contiguous grain IDs are stored per IDBucket. The smaller this value is chosen,
|	the more buckets are required but also the proportionally lower it is the time to find a specific grain.
|	Clearly, the value is performance relevant, as the larger the value is set the proportionally more grains have to be 
|	scanned on average before finding a specific. On the contrary, a very small value will reduce the total
|	number of checks but eventually may not end up in pointing to aligned pieces of memory, unless additional measures of 
|	improving contiguity (with for instance the jemalloc allocator class) are enforced. 
|	Hence, the value should be in the order of the size of single cache block times size(unsigned int). 100 is a sensible default. 
|	The choice to go for a concept of IDBuckets rather than a hashfield of pointer is to avoid the 
|	storing of many NULL pointer for non-existent IDs for those snapshots in which only few IDs remain.
| **MemRegionsX**
| **MemRegionsY**
| **MemRegionsZ**
|	**Required to be 1**
|	An integer spatial partitioning of the domain into regions. All grains whose barycenter are in the region are
|	handled by an independent memory handler that is a local memory region class object. This improves memory locality for a future
|	OpenMP extension of the code, as then shared memory accessess within the region are more likely found in local caches.
| **DeveloperMode**
|	With this option enabled (1) each binary file is translated into an ASCII file.
|	**Make sure to disable this option (0) whenever running production jobs!**