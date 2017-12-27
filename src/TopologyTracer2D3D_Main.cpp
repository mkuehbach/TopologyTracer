/*
	Copyright Markus Kühbach, 2015-2017

	TopologyTracer is an MPI-parallel datamining tool for the conducting of spatiotemporal 
	analyses of microstructural element dynamics. Its purpose is the quantification of correlations
	--- spatial and temporal --- between individual microstructural elements and their 
	higher-order neighboring elements. Specifically, its current functionalities allow to 
	study the volume evolution of individual grains over time and set their growth history into 
	relation to the evolution of the higher-order neighbors. The tool is unique insofar as it 
	allows processing these individual surveys in a parallelized manner. Thus, enabling the 
	post-processing of so far intractable large datasets.

	The source code was developed by Markus Kühbach during his PhD time with Luis A. Barrales-Mora 
	and Günter Gottstein at the Institute of Physical Metallurgy and Metal Physics with RWTH Aachen University. 
	Being now with the Max-Planck-Institut fur Eisenforschung GmbH in Dusseldorf, I maintain the code, 
	though at disregular intervals. Nonetheless, feel free to utilize the tool, do not hesitate contacting 
	me for sharing thoughts, suggesting improvements, or reporting your experiences.
	markus.kuehbach at rwth-aachen.de and m.kuehbach at mpie.de


	The authors gratefully acknowledge the financial support from the Deutsche Forschungsgemeinschaft
	(DFG) within the Reinhart Koselleck-Project (GO 335/44-1) and computing time grants kindly provided
	by RWTH Aachen University and the FZ Jülich within the scope of the JARAHPC project JARA0076.


	This file is part of TopologyTracer.

	TopologyTracer is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	TopologyTracer is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with SCORE.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "TopologyTracer2D3D_Topohdl.h"


using namespace std;

//GENERAL STARTUP PARAMETER
#define	SIMID						1		//unique unsigned int identifier to distinguish analysis output
#define PARAMETERFILE				2		//name of XML control file

//additional or optional startup parameter for specific analysis modes
//optional for E_NUCSITES_LONGRANGE3D
#define FID							3
#define EXTENT						4

//additional for E_APPROX_CURVATURE_CONTOUR
#define SUFFIX						3		//additional naming suffix


#define VERSION_MAJOR				1
#define VERSION_MINOR				2
#define VERSION_REVISION			0
#define VERSION_BUILD				0
//#define IS_BETAVERSION_YES				//MK::no this source code version is not

void versioninfo( int myrank, int nranks ) {
		cout << "TopologyTracer2D3D starting rank " << myrank << " within " << nranks << endl;
		cout << "v" << VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_REVISION;
#ifdef IS_BETAVERSION_YES
		cout << "beta" << endl;
#else
		cout << endl;
#endif
//		cout << "SCORE the coordinate system is right-handed with x pointing right and ||RD, y pointing inwards and ||TD, and z pointing upwards and ||ND" << endl;
//#ifdef USE_QUATLIB_DEGRAEF
//		cout << "SCORE the mathlibrary is the Rowenhorst/Rollett/Konijnenberg/deGraef et. al. library" << endl;
//#endif
//#ifdef USE_QUATLIB_DAMASK
//		cout << "SCORE the mathlibrary is the DAMASK library" << endl;
//#endif
//#ifdef USE_QUATLIB_IMM
//		cout << "SCORE the mathlibrary is the IMM library" << endl;
//#endif
}


void TrackingInParallel(int pargc, char** pargv, unsigned int sid)
{
//INITIALIZATION
		//MK::Init MPI capable of hybrid MPI/OpenMP even though the TopologyTracer currently at no point calls MPI commands within an OpenMP parallel region
	int supportlevel_desired = MPI_THREAD_FUNNELED;
	int supportlevel_provided;
	MPI_Init_thread(&pargc, &pargv, supportlevel_desired, &supportlevel_provided);

	int nranks, rank;
	if (supportlevel_provided < supportlevel_desired) {
		versioninfo(MASTER, SINGLE_PROCESS);
		cerr << "ERROR::Worker::TrackingParallel Insufficient threading capabilities of the MPI library!"; return;
	}
	else {
		MPI_Comm_size(MPI_COMM_WORLD, &nranks);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	}

	double timer = MPI_Wtime();
	versioninfo(rank, nranks);
	if (rank == MASTER)
		cout << "Starting parallelized tracking of grains in " << Settings::Dimensionality << "d..." << endl;

	topoHdlP worker = new topoHdl;
	worker->set_MPICommunication(rank, nranks);
	worker->init_MPIDatatypes();

	//mode specific checks and reading of data before starting heavy I/O machinery
	if (Settings::AnalyzeTrajectoriesForward == true && Settings::AnalyzeTrajectoriesForwardMode == E_FORWARD_SELECTED) {
		if (worker->read_targets_fromlist() != true) {
			cerr << "ERROR::Worker::TrackingParallel " << worker->get_Rank() << " list of TargetGrainIDs was nonexistent or errornous!" << endl;
			delete worker; worker = NULL; MPI_Finalize();
			return;
		}
	}
	worker->myprofiler.logev("MPIStartup" + std::to_string(Settings::SimID), (double)(MPI_Wtime() - timer));

//DATASET MANAGEMENT
	timer = MPI_Wtime();
	unsigned int DatasetGood[1] = { 0 };
	if (worker->get_Rank() == MASTER) { //MASTER scans through all snapshots and determines file size to organize the workload distribution, then synchronizes with waiting partners
		cout << "...MASTER probing rawdata set for analysis run SimID " << sid << " with " << Settings::LocalDatabaseMaximumSize << " byte per process --- ";
		if (Settings::ProbeWorkPartMode == true)
			cout << "Probing mode enabled...\n";
		else
			cout << " Developer mode disabled...\n";

		if (worker->determineWorkload() == true)
			DatasetGood[0] = 1;
	}
	MPI_Barrier(MPI_COMM_WORLD); //tell all others whether data set is good, MK::probing done with only MASTER to reduce file sytem load 
	MPI_Bcast(&DatasetGood, 1, MPI_UNSIGNED, MASTER, MPI_COMM_WORLD);
	if (DatasetGood[0] == 0) {
		cerr << "ERROR::Worker::TrackingParallel " << worker->get_Rank() << " dataset is incomplete or errornous!" << endl;
		delete worker; worker = NULL; MPI_Finalize();
		return;
	}

	//ARRANGE WORKLOAD DISTRIBUTION	
	worker->commitWorkload();

	if (Settings::ProbeWorkPartMode == true) {
		cout << "...Worker " << worker->get_Rank() << " workpartitioning was probed, now terminating!" << endl;
		delete worker; worker = NULL; MPI_Finalize(); return;
	}

	worker->myprofiler.logev("DatasetManagement" + std::to_string(Settings::SimID), (double)(MPI_Wtime() - timer));
	timer = MPI_Wtime();

	//DATASET I/O COLLECTING PHASE
	unsigned int workerReadProcessGood = 0; //assume reading provokes an error false
	if (worker->readDatasets() == true) {
		workerReadProcessGood = 1;
	}

	worker->myprofiler.logev("DatasetIOHandling" + std::to_string(Settings::SimID), (double)(MPI_Wtime() - timer));
	timer = MPI_Wtime();

	unsigned int allReadProcessGood = 0; //communicate potential errors over process pool
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Allreduce(&workerReadProcessGood, &allReadProcessGood, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
	if (allReadProcessGood != worker->get_nRanks()) { //somebody unsuccessful?
		if (worker->get_Rank() == MASTER) {
			cout << "ERROR::Worker " << worker->get_Rank() << " reading rawdata was unsuccessful over the MPI process group!" << endl;
		}
		MPI_Finalize(); delete worker; worker = NULL;
		return;
	}

	//DATA ANALYSIS PHASE ALL PROCESS READ RAWDATA SUCCESSFULLY
	if (worker->get_Rank() == MASTER) {
		cout << "Worker " << worker->get_Rank() << " all datasets were read in successfully, now analyzing ..." << endl;
	}


	//PERFORM INDIVIDUAL DATA ANALYSES, START WITH DESCRIPTIVE COUNTING STUFF
	MPI_Barrier(MPI_COMM_WORLD); //wait for all workers becoming ready
	worker->myprofiler.logev("WaitingBeforeAnalyses" + std::to_string(Settings::SimID), (double)(MPI_Wtime() - timer)); timer = MPI_Wtime();

	if (Settings::AnalyzeGrainSizeQuantiles == true)
		worker->analyze_grainsize_quantiles_db();

	MPI_Barrier(MPI_COMM_WORLD);
	worker->myprofiler.logev("AnalyzeGSDQuantiles" + std::to_string(Settings::SimID), (double)(MPI_Wtime() - timer));  timer = MPI_Wtime();

	if (Settings::AnalyzeSEE == true)
		worker->analyze_see_db();

	MPI_Barrier(MPI_COMM_WORLD);
	worker->myprofiler.logev("AnalyzeSEEDF" + std::to_string(Settings::SimID), (double)(MPI_Wtime() - timer));  timer = MPI_Wtime();

	if (Settings::AnalyzeMODF == true)
		worker->analyze_modf_db();

	MPI_Barrier( MPI_COMM_WORLD );
	worker->myprofiler.logev("AnalyzeMODF" + std::to_string(Settings::SimID), (double)(MPI_Wtime() - timer));  timer = MPI_Wtime();

	if (Settings::AnalyzeGSD == true)
		worker->analyze_gsd_db();

	MPI_Barrier( MPI_COMM_WORLD );
	worker->myprofiler.logev("AnalyzeGSD" + std::to_string(Settings::SimID), (double)(MPI_Wtime() - timer));  timer = MPI_Wtime();

	if (Settings::AnalyzeDrvForceSEE == true)
		worker->analyze_drvforce_see_db();

	MPI_Barrier(MPI_COMM_WORLD);
	worker->myprofiler.logev("AnalyzeDrvForceSEE" + std::to_string(Settings::SimID), (double)(MPI_Wtime() - timer));  timer = MPI_Wtime();

	if ( Settings::AnalyzeMaxSizeGainForward == true )
		worker->analyze_maxsizegain_fw_db();
	
	MPI_Barrier( MPI_COMM_WORLD ); 
	worker->myprofiler.logev("AnalyzeMaxSzGainFW" + std::to_string(Settings::SimID), (double)(MPI_Wtime() - timer));  timer = MPI_Wtime();

	if ( Settings::AnalyzeApproxRXFraction == true ) 
		worker->analyze_rxfraction_db( Settings::SnapshotLast );

	MPI_Barrier( MPI_COMM_WORLD );
	worker->myprofiler.logev("AnalyzeRXFraction" + std::to_string(Settings::SimID), (double)(MPI_Wtime() - timer));  timer = MPI_Wtime();

	if ( Settings::AnalyzeTrajectoriesForward == true )
		worker->analyze_trajectories_fw_db();

	MPI_Barrier( MPI_COMM_WORLD );
	worker->myprofiler.logev("AnalyzeTrajectoriesFW" + std::to_string(Settings::SimID), (double)(MPI_Wtime() - timer));  timer = MPI_Wtime();

	if ( Settings::AnalyzeTrajectoriesBackward == true )
		worker->analyze_trajectories_bk_db();

	MPI_Barrier( MPI_COMM_WORLD ); 
	worker->myprofiler.logev("AnalyzeTrajectoriesBK" + std::to_string(Settings::SimID), (double)(MPI_Wtime() - timer));  timer = MPI_Wtime();

	if ( Settings::AnalyzeAbnormalGrains == true )
		worker->analyze_abnormalgraingrowth_db();

	MPI_Barrier( MPI_COMM_WORLD );
	worker->myprofiler.logev("AnalyzeAGG" + std::to_string(Settings::SimID), (double)(MPI_Wtime() - timer));  timer = MPI_Wtime();

	if ( Settings::AnalyzeSizeGainVsMatrixBackward == true )
		worker->analyze_sizegain_vs_matrix_bk_db();

	MPI_Barrier(MPI_COMM_WORLD);
	worker->myprofiler.logev("AnalyzeSzVsMatrixBK" + std::to_string(Settings::SimID), (double)(MPI_Wtime() - timer));  timer = MPI_Wtime();

	if (Settings::AnalyzeKNN == true)
		worker->analyze_knn_naive_db( Settings::MaxNumberOfKShells );

	MPI_Barrier(MPI_COMM_WORLD);
	worker->myprofiler.logev("AnalyzeKNN" + std::to_string(Settings::SimID), (double)(MPI_Wtime() - timer));  timer = MPI_Wtime();

//##MK::begin of undocumented functionality
	if ( Settings::AnalyzeClassicalNucModels == true )
		worker->analyze_classical_nucmodels_db();

	MPI_Barrier(MPI_COMM_WORLD);
	worker->myprofiler.logev("AnalyzeClassicalNucModels" + std::to_string(Settings::SimID), (double)(MPI_Wtime() - timer));  timer = MPI_Wtime();

//##MK::end of undocumented functionality
	cout << "...Worker " << worker->get_Rank() << " DONE WITH ALL ANALYSES, now clearing memory from rawdata..." << endl;
	worker->spit_profiling( "TrackParallel" );

	delete worker; worker = NULL;	
	MPI_Finalize();
}


void TrackingSequential (int pargc, char** pargv, unsigned int sid )
{
//INITIALIZATION
	//MK::Init MPI capable of hybrid MPI/OpenMP even though the TopologyTracer currently at no point calls MPI commands within an OpenMP parallel region
	int supportlevel_desired = MPI_THREAD_FUNNELED;
	int supportlevel_provided;
	MPI_Init_thread(&pargc, &pargv, supportlevel_desired, &supportlevel_provided);

	int nranks, rank;
	if (supportlevel_provided < supportlevel_desired) {
		versioninfo(MASTER, SINGLE_PROCESS);
		cerr << "ERROR::TrackingSeq Insufficient threading capabilities of the MPI library!"; return;
	}
	else {
		MPI_Comm_size(MPI_COMM_WORLD, &nranks);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	}

	double timer = MPI_Wtime();
	versioninfo(rank, nranks);
	if (rank == MASTER)
		cout << "Starting parallelized tracking of grains in " << Settings::Dimensionality << "d..." << endl;

	topoHdlP worker = new topoHdl;
	worker->set_MPICommunication(rank, nranks);
	worker->init_MPIDatatypes();

	//mode specific checks and reading of data before starting heavy I/O machinery
	if (Settings::AnalyzeTrajectoriesForward == true && Settings::AnalyzeTrajectoriesForwardMode == E_FORWARD_SELECTED) {
		if (worker->read_targets_fromlist() != true) {
			cerr << "ERROR::Worker::TrackingSeq " << worker->get_Rank() << " list of TargetGrainIDs was nonexistent or errornous!" << endl;
			delete worker; worker = NULL; MPI_Finalize();
			return;
		}
	}
	worker->myprofiler.logev("MPIStartup" + std::to_string(Settings::SimID), (double)(MPI_Wtime() - timer)); timer = MPI_Wtime();

//DATASET HANDLING
	//MK::as we work sequentially through the dataset we cannot easily check the boundary contact evolution of the grains
	//therefore we need first to probe the heavy dataset to generate a contact history table
	//the processes can then read from this table immediately at which snapshot a grain made contact if at all and hence
	//the analyze_elimbnd functions become much leaner and thus faster
	//MK::in general we always query the existence of the files but assume the heavy to remain unchanged

	if (worker->analyze_rvebnd_contacttime() == false) {
		delete worker; worker = NULL; MPI_Finalize(); return;
	}
	
	worker->myprofiler.logev("AnalyzeRVEBndContactTextureIO" + std::to_string(Settings::SimID), (double)(MPI_Wtime() - timer)); timer = MPI_Wtime();

	//now all nodes know at which point a grain if at all during its existence made contact with its boundary with its boundary
	//master sets up results file, into which the processes dump results for individual time steps,
	//so were running in parallel but are not required to keep all data at one time step in global memory!
	if (Settings::ProbeWorkPartMode == true) {
		cout << "...Worker " << worker->get_Rank() << " workpartitioning was probed, now terminating!" << endl;
		delete worker; worker = NULL; MPI_Finalize(); return;
	}

//PARALLELIZED DATASET PROCESSING BUT LOADING INDIVIDUAL FILES ONE AFTER ANOTHER
	//so far LocalDB entries to process the different analyses tasks were not established yet...
	//results files are initialized, processes can start performing individual analyses by reading files successively and processing
	//MK::however, we require to distinguish two categories of tasks depending whether they require prior filtering of grainIDs or not,
	//i.e. trivially executable analysis tasks with all IDs we have to distinguish from non-trivial filtering list dependent tasks
	//for the latter we observe that one requires in general different grainID lists, filtered for certain conditions::
	//list of surviving grains --> never touching a RVE boundary + those remaining in one reference dataset, which for BK tracking is Settings::SnapshotLast
	//list of matrix --> never ever touching a RVE boundary but may still include grains which eventually after some time have died
	//additionally combinations of these filter criteria may be necessary lists to consider via exclusion (Black) or inclusion (White) lists/criteria

	//thus, ALL processes have to load at least Settings::SnapshotLast by loading Texture_<Settings::SnapshotLast>.bin binary
	
	unsigned int localhealth = 1;
	bool prepare = worker->readSingleDataset(Settings::SnapshotLast);
	worker->myprofiler.logev("InitAnalysisReadSnapshotLast" + std::to_string(Settings::SimID), (double)(MPI_Wtime() - timer)); timer = MPI_Wtime();

	if (prepare == false)
		localhealth = 0;
	
	unsigned int globalhealth = 0;
	MPI_Allreduce(&localhealth, &globalhealth, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
	if (globalhealth != worker->get_nRanks() ) {
		cout << "ERROR::Worker::TrackingSeq " << worker->get_Rank() << " not all were able to load snapshot Settings::SnapshotLast!" << endl;
		delete worker; worker = NULL; MPI_Finalize(); return;
	}
	
	globalhealth = 0;
	if (worker->get_Rank() == MASTER) {
		if (worker->init_results_container() == true)
			globalhealth = 1;
	}
	MPI_Bcast(&globalhealth, 1, MPI_UNSIGNED, MASTER, MPI_COMM_WORLD);
	if (globalhealth != 1) {
		cout << "ERROR::Worker::TrackingSeq " << worker->get_Rank() << " MASTER was unable to setup results files!" << endl;
		delete worker; worker = NULL; MPI_Finalize(); return;
	}

	worker->myprofiler.logev("InitAnalysisMemSetup" + std::to_string(Settings::SimID), (double)(MPI_Wtime() - timer)); timer = MPI_Wtime();

	//not returned yet, then all processes have loaded the last timestep and can rely upon that all results container were initialized
//COOPERATIVE SUCCESSIVE TRACKING SEQUENTIALLY

	double tic = MPI_Wtime();

	//apply filter criteria
	std::vector<unsigned int>* survgrs = NULL;
	survgrs = worker->analyze_elimbnd_one_mpiioself("ADHOCSurvivors", NULL, Settings::SnapshotLast);

	std::vector<unsigned int>* matrixgrs = NULL;
	matrixgrs = worker->analyze_elimbnd_all_self("ADHOCMatrix", NULL, false);
		
	//##MK::invokes twice the reading of the last timestep
	std::vector<unsigned int>* rxgrs = NULL;
	rxgrs = worker->analyze_elimbnd_one_mpiioself("ADHOCRX", matrixgrs, Settings::SnapshotLast); 

	std::vector<unsigned int>* matrixgrs0survgrs = NULL;
	matrixgrs0survgrs = worker->analyze_elimbnd_all_self("ADHOCMatrixWithoutSurv", survgrs, false);

	std::vector<unsigned int>* fwtargets = NULL;
	if (Settings::AnalyzeTrajectoriesForward == true || Settings::AnalyzeTopologyDifferenceForward == true) {
		if (Settings::AnalyzeTrajectoriesForwardMode == E_FORWARD_ALL)
			fwtargets = worker->analyze_elimbnd_all_self("FWTargets", NULL, false);
		else if (Settings::AnalyzeTrajectoriesForwardMode == E_FORWARD_SELECTED)
			fwtargets = worker->analyze_elimbnd_all_self("FWTargets", NULL, true);
		else
			fwtargets = NULL; //remains NULL
	}
	//helper to accumulate local for some non-trival analysis tasks such as maxszgain_adhoc
	worker->init_helper();

	double toc = MPI_Wtime();
	worker->myprofiler.logev("SetupFilterLists" + std::to_string(Settings::SimID), (toc-tic)); tic = MPI_Wtime();

	//now we can process all analysis tasks, for this we load the heavy data per time step only, process, delete them, move to the next timestep
	//processes work in parallel on disjoint datasets...
	for (unsigned int fid = Settings::SnapshotFirst; fid <= Settings::SnapshotLast; fid += Settings::SnapshotOffset) {

		if (worker->workpartitioning_fid(fid) == true) {
			tic = MPI_Wtime();
				
			bool dbgood = worker->readDatasets(fid);
				
			toc = MPI_Wtime(); worker->ioprofiler.logev("LoadSnapshot." + std::to_string(Settings::SimID), (toc - tic)); tic = MPI_Wtime();

			//MK::load the heavy data for this via a database entry to enable reutilization of algorithms from the coorperatively parallel execution
			if (dbgood == true) {
				tic = MPI_Wtime();

				if (Settings::AnalyzeGrainSizeQuantiles == true) {
					worker->analyze_grainsize_quantiles_adhoc(fid);
					toc = MPI_Wtime(); worker->myprofiler.logev("AnalyzeGSDQuantiles." + std::to_string(Settings::SimID) + "." + std::to_string(fid), (toc - tic));  tic = MPI_Wtime();
				}
				if (Settings::AnalyzeSEE == true) {
					worker->analyze_see_adhoc(fid);
					toc = MPI_Wtime(); worker->myprofiler.logev("AnalyzeSEEDF." + std::to_string(Settings::SimID) + "." + std::to_string(fid), (toc - tic));  tic = MPI_Wtime();
				}
				if (Settings::AnalyzeMODF == true) {
					worker->analyze_modf_adhoc(fid);
					toc = MPI_Wtime(); worker->myprofiler.logev("AnalyzeMODF." + std::to_string(Settings::SimID) + "." + std::to_string(fid), (toc - tic));  tic = MPI_Wtime();
				}
				if (Settings::AnalyzeGSD == true) {
					worker->analyze_gsd_adhoc(fid);
					toc = MPI_Wtime(); worker->myprofiler.logev("AnalyzeGSD." + std::to_string(Settings::SimID) + "." + std::to_string(fid), (toc - tic));  tic = MPI_Wtime();
				}
				if (Settings::AnalyzeDrvForceSEE == true) {
					worker->analyze_drvforce_see_adhoc(fid);
					toc = MPI_Wtime(); worker->myprofiler.logev("AnalyzeDrvForceSEE." + std::to_string(Settings::SimID) + "." + std::to_string(fid), (toc - tic));  tic = MPI_Wtime();
				}
				if (Settings::AnalyzeMaxSizeGainForward == true) {
					worker->analyze_maxsizegain_fw_adhoc(fid);
					toc = MPI_Wtime(); worker->myprofiler.logev("AnalyzeMaxSzGainFW." + std::to_string(Settings::SimID) + "." + std::to_string(fid), (toc - tic));  tic = MPI_Wtime();
				}
				if (Settings::AnalyzeMeanDrvForceSEEForward == true) {
					worker->analyze_meandrvforcesee_adhoc(fid);
					toc = MPI_Wtime(); worker->myprofiler.logev("AnalyzeMeanDrvForceSEEFW." + std::to_string(Settings::SimID) + "." + std::to_string(fid), (toc - tic));  tic = MPI_Wtime();
				}

				if (Settings::AnalyzeApproxRXFraction == true) {
					worker->analyze_rxfraction_adhoc(fid, matrixgrs, rxgrs);
					toc = MPI_Wtime(); worker->myprofiler.logev("AnalyzeRXFraction." + std::to_string(Settings::SimID) + "." + std::to_string(fid), (toc - tic));  tic = MPI_Wtime();
				}
				if (Settings::AnalyzeTrajectoriesForward == true) {
					worker->analyze_trajectories_fw_adhoc(fid, fwtargets);
					toc = MPI_Wtime(); worker->myprofiler.logev("AnalyzeTrackingFW." + std::to_string(Settings::SimID) + "." + std::to_string(fid), (toc - tic));  tic = MPI_Wtime();
				}
				if (Settings::AnalyzeTrajectoriesBackward == true) {
					worker->analyze_trajectories_bk_adhoc(fid, survgrs);
					toc = MPI_Wtime(); worker->myprofiler.logev("AnalyzeTrackingBK." + std::to_string(Settings::SimID) + "." + std::to_string(fid), (toc - tic));  tic = MPI_Wtime();
				}
				if (Settings::AnalyzeAbnormalGrains == true) {
					worker->analyze_abnormalgraingrowth_adhoc(fid, survgrs, matrixgrs, matrixgrs0survgrs);
					toc = MPI_Wtime(); worker->myprofiler.logev("AnalyzeAGG." + std::to_string(Settings::SimID) + "." + std::to_string(fid), (toc - tic));  tic = MPI_Wtime();
				}
				if (Settings::AnalyzeSizeGainVsMatrixBackward == true) {
					worker->analyze_sizegain_vs_matrix_bk_adhoc(fid, survgrs, matrixgrs0survgrs); //MK::former functions passed matrixgrs but we seek to pass the matrix excluding survivors
					toc = MPI_Wtime(); worker->myprofiler.logev("AnalyzeSzVsMatrixBK." + std::to_string(Settings::SimID) + "." + std::to_string(fid), (toc - tic));  tic = MPI_Wtime();
				}
				if (Settings::AnalyzeTopologyDifferenceForward == true) {
					worker->analyze_topodiff_fw_adhoc(fid, fwtargets);
					toc = MPI_Wtime(); worker->myprofiler.logev("AnalyzeTopoDiffFW." + std::to_string(Settings::SimID) + "." + std::to_string(fid), (toc-tic));  tic = MPI_Wtime();
				}
				//MK::heavy data will be cleared prior to the loading of the next dataset or latest during destruction of worker class instance
			}
			else {
				cerr << "ERROR::Worker::TrackingSeq " << worker->get_Rank() << " dataset " << fid << " was not processable!" << endl;
			}
		}
		//else nothing, different process will take care of it!
	} //next snapshot

	//during the local analyzing of non-trivial tasks the worker have collected individual summary results which still need consolidation 
	//MK::the barrier is necessary, because we have to wait for all data to arrive
	double t1 = MPI_Wtime();

	MPI_Barrier(MPI_COMM_WORLD);
	worker->consolidate_helper();

	double t2 = MPI_Wtime();  worker->myprofiler.logev("MPIBarrierAndConsolidateAnalyses" + std::to_string(Settings::SimID), (t2-t1));
	worker->spit_profiling( "TrackSequential" );

	//were done with all analyses, so we can delete the filter lists
	if (survgrs != NULL) { delete survgrs; survgrs = NULL; }
	if (matrixgrs != NULL) { delete matrixgrs; matrixgrs = NULL; }
	if (rxgrs != NULL) { delete rxgrs; rxgrs = NULL; }
	if (matrixgrs0survgrs != NULL) { delete matrixgrs0survgrs; matrixgrs0survgrs = NULL; }
	if (fwtargets != NULL) { delete fwtargets; fwtargets = NULL; }
	

	delete worker; worker = NULL;
	MPI_Finalize();
}





void UtilityKNN( int pargc, char** pargv, unsigned int sid )
{
	int nranks, rank;
	MPI_Init( NULL, NULL );
	MPI_Comm_size(MPI_COMM_WORLD, &nranks);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	versioninfo( rank, nranks);
	if ( rank == MASTER) 
		cout << "Starting the KNN Utility tools..." << endl;

	spatQueryHdlP worker = new spatQueryHdl;
	worker->set_MPICommunication(rank, nranks);

	bool status = true;
	
	status = worker->init_spatquerybox();
	status = worker->read_microstructure_uds();
	status = worker->init_fillspatquerybox();
	status = worker->read_maxgrainsizes( "TopoTracer2D3D.SimID.109500.MaxSizeGainFW.F.10.O.1.L.9500.csv" );
	
	if (status == true) {
		worker->init_neighbordistances();

		for ( unsigned int fid = Settings::SnapshotFirst; fid <= Settings::SnapshotLast; fid += Settings::SnapshotOffset ) {
			unsigned int f = 1 + (fid - Settings::SnapshotFirst) / Settings::SnapshotOffset;
			if ( (f % worker->get_nRanks()) == worker->get_Rank() ) { //me working on this dataset

				if ( worker->read_longrange( fid ) == true ) {
					worker->arrivaltime_quantiles( fid );
				}

			}
			//other worker taking care of, all have loaded rawdata and computed with the same algorithm the neighbors
		}
	}
	
	delete worker;
	MPI_Finalize();

/*
	utilityHdlP worker = new utilityHdl;
	worker->set_MPICommunication(  rank, nranks );
	
	//worker->read_targetfile();
	worker->read_knnfile( Settings::MaxNumberOfKShells );

	worker->analyze_apriori_prediction();
	//worker->analyze_find_knn_for_targets();
	//worker->analyze_Gest_for_targets();

	delete worker;
	MPI_Finalize();
*/
}


void LongRangeEnvPersistenceAnalysis2D( int pargc, char** pargv )
{
	//MK::Init MPI capable of hybrid MPI/OpenMP even though the TopologyTracer currently at no point calls MPI commands within an OpenMP parallel region
	int supportlevel_desired = MPI_THREAD_FUNNELED;
	int supportlevel_provided;
	MPI_Init_thread(&pargc, &pargv, supportlevel_desired, &supportlevel_provided);

	int nranks, rank;
	if ( supportlevel_provided < supportlevel_desired ) {
		versioninfo( MASTER, SINGLE_PROCESS );
		cerr << "ERROR::LR2D Insufficient threading capabilities of the MPI library!"; return;
	}
	else {
		MPI_Comm_size(MPI_COMM_WORLD, &nranks);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	}
	versioninfo( rank, nranks );
	if ( rank == MASTER) 
		cout << "Starting computation of unbiased growth measures to higher-order neighbors analytically in 2d..." << endl;

//#ifdef USE_POLYTRI
	arealanaHdlP worker2d = NULL;
	worker2d = new arealanaHdl;
	worker2d->set_MPICommunication( rank, nranks );
	worker2d->init_MPIDatatypes();
	worker2d->init();
	worker2d->read_initialmicrostructure_metadata();

	//process all files
	//##MK::consider optimization of load balancing because as time progresses the population shrinks in total grain number hence formal complexity
	for ( unsigned int fid = Settings::SnapshotFirst; fid <= Settings::SnapshotLast; fid += Settings::SnapshotOffset ) {
		unsigned int f = 1 + (fid - Settings::SnapshotFirst) / Settings::SnapshotOffset;
		unsigned int nr = worker2d->get_nRanks();
		unsigned int me = worker2d->get_Rank();
		if ( (f % nr) == me ) { //##MK::naive maner of partitioning, because first files have higher processing complexity than following ones...

			if ( worker2d->read_snapshot_2d( fid ) == true ) {

				worker2d->modify_rawdata_gragles( fid );

				std::vector<unsigned int>* elim = NULL;
				worker2d->fdmeta = analysismetrics2d();

				elim = worker2d->triangularize_contours( fid ); //elim will be allocated if size == 0 all contours are fine!
				
				if ( elim != NULL && worker2d->fdmeta.TotalAreaEnclosedByContour > MINIMUM_REQUIRED_TOTALAREA ) { //at least a significant area of the total domain was triangularizable
					//triangularization was successful with at most only a few bad grains, proceed
					if ( elim->size() > 0 ) { //something to eliminate?
						worker2d->eliminate_non_triangularized_grains( fid, elim );
					}

					//MK::heuristic fraction of grains to accept eluded, proceed only when triangularization was successful!
					//and from this total not more than a total area of triangularizationSuccess[1] was covered by
					//not successfully triangularizable contours, then we perform the analysis at all
					if ( worker2d->fdmeta.TotalAreaNonTriangularizable < MAXIMUM_ADMISSIBLE_NONTRIANGULARIZED ) { //less than 1% grains
					//if ( elim->size() < 10 ) { 
						
						worker2d->envcharacterize_omp_naive( fid );
						worker2d->report( fid );
						
						//worker2d->metareport();
					}
				} //no, triangularization was unsuccessful
				else {
					cerr << endl << endl << "ERROR::Worker::Poly2Tri " << worker2d->get_Rank() << " area on unit square domain not successfully triangularizable " << worker2d->fdmeta.TotalAreaNonTriangularizable << endl;
 					cerr << "ERROR::Worker::Poly2Tri " << worker2d->get_Rank() << " Poly2Tri had too many problems triangularizing!" << endl << endl;
				}

				delete elim; elim = NULL;
				worker2d->reset_contour_memory( fid );

				cout << endl << endl;
			} //done with fid
		} //next fid
	}
	worker2d->spit_profiling();

	delete worker2d; worker2d = NULL;
//#endif

	MPI_Finalize();
}


void LongRangeEnvPersistenceAnalysis3D( int pargc, char** pargv )
{
	//discrete approximation of unbiased growth propensity measures for each grain answering the following question:
	//how do state variable values in a fixed reference volume about a grain evolve in time?

	//MK::initially my idea was a straightforward extension of the 2d code, i.e. grains build out of tetrahedra partially intersecting
	//inspection volume, however a quick trial showed this currently to be of too high formal complexity, ... therefore still discrete approximation
	//see further comments below

	//MK::Init MPI capable of hybrid MPI/OpenMP even though the TopologyTracer currently at no point calls MPI commands within an OpenMP parallel region
	int supportlevel_desired = MPI_THREAD_FUNNELED;
	int supportlevel_provided;
	MPI_Init_thread(&pargc, &pargv, supportlevel_desired, &supportlevel_provided);

	int nranks, rank;
	if ( supportlevel_provided < supportlevel_desired ) {
		versioninfo( MASTER, SINGLE_PROCESS ); 
		cerr << "ERROR::LR3D Insufficient threading capabilities of the MPI library!"; return;
	}
	else {
		MPI_Comm_size(MPI_COMM_WORLD, &nranks);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	}
	versioninfo( rank, nranks );
	if ( rank == MASTER) 
		cout << "Starting computation of unbiased growth measures to higher-order neighbors via MC in 3d..." << endl;

	volanaHdlP worker3d = NULL;
	worker3d = new volanaHdl;
	worker3d->set_MPICommunication( rank, nranks );
	worker3d->init_MPIDatatypes();
	worker3d->init_discrete_version();
	worker3d->read_initialmicrostructure_metadata();

	//process all files
	//##MK::consider optimization of load balancing because as time progresses the population shrinks in total grain number hence formal complexity
	for ( unsigned int fid = Settings::SnapshotFirst; fid <= Settings::SnapshotLast; fid += Settings::SnapshotOffset ) {
		unsigned int f = 1 + (fid - Settings::SnapshotFirst) / Settings::SnapshotOffset;
		unsigned int nr = worker3d->get_nRanks();
		unsigned int me = worker3d->get_Rank();
		if ( (f % nr) == me ) { //##MK::naive maner of partitioning, because first files have higher processing complexity than following
			worker3d->read_snapshot_3d( fid );
			worker3d->calc_grainvolume( fid );
			worker3d->envcharacterize_omp_naive( worker3d->AssgnSnapshotsMeta.at(worker3d->AssgnSnapshotsMeta.size()-1) );
			worker3d->report( worker3d->AssgnSnapshotsMeta.at(worker3d->AssgnSnapshotsMeta.size()-1) );
			worker3d->reset_xitable();
		}
	}

	worker3d->spit_profiling();

	delete worker3d; worker3d = NULL;
	MPI_Finalize();

	//further details to formal complexity of intersection computation
	//##MK::this code snippet implements a trial to compute the intersection volume analytically by utilizing triangularized contour hulls
	//##MK::It requires the TopologyTracer2D3D to link against the TetGen tetrahedralization library, however that one is not thread safe!
	/*
		//worker3d->test_tetgen();
		worker3d->init();
		worker3d->test_singlegrain_vtk( "GrainHull_1068Timestep_660.vtk", 1068 );
		worker3d->test_singlegrain_modifydata_GraGeLeS( 1068 );
		if ( worker3d->test_singlegrain_tetgen( 1068 ) == false ) { cerr << "ERROR::LR3D Tetrahedralization of hull " << 1068 << " was unsuccessful!" << endl; }
	*/
	//##MK::the trial with such a grain showed that given 500k of such guys each with approx. 2000 disjoint supporting points per hull, 
	//##MK::4000 tris, hence ~100KB/hull and an increase by a factor of 5-10 the number of triangles of the PLC results resulting in 
	//##MK::approx. 300GB only for tetrahedra  storage and one time step, tetrahedralization took approx. 2s, hence formal complexity is too high
	//##MK::if no multi-threaded tetrahedralization library available...
}


void ApproximateCurvatureFromContourData2D( int pargc, char** pargv )
{
	//MK::Init MPI capable of hybrid MPI/OpenMP even though the TopologyTracer currently at no point calls MPI commands within an OpenMP parallel region
	int supportlevel_desired = MPI_THREAD_FUNNELED;
	int supportlevel_provided;
	MPI_Init_thread(&pargc, &pargv, supportlevel_desired, &supportlevel_provided);

	int nranks, rank;
	if ( supportlevel_provided < supportlevel_desired ) {
		versioninfo ( MASTER, SINGLE_PROCESS );
		cerr << "ERROR::ApproxCurv Insufficient threading capabilities of the MPI library!";
		return;
	}
	else {
		//cout << "Initialize MPI_COMM_WORLD..." << endl;
		MPI_Comm_size(MPI_COMM_WORLD, &nranks);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		//cout << "\t\tMPI firing up rank " << rank << " within MPI_COMM_WORLD of " << nranks << endl;
	}
	versioninfo( rank, nranks );
	if ( rank == MASTER) 
		cout << "Starting approximation of integral mean curvature from 2D linearized contours..." << endl;

	curvapprxHdlP worker2d = new curvapprxHdl;
	worker2d->set_MPICommunication( rank, nranks );
	worker2d->init_MPIDatatypes();
	if ( worker2d->init() == true ) { 

		//##MK::in the future add support for what happens when one pointer on gcontour2d is NULL and implement
		//cross-check across MPI_Comm_world to assure that all workers are still there..., if not jump out prior processing any...

		for ( unsigned int fid = Settings::SnapshotFirst; fid <= Settings::SnapshotLast; fid = fid + Settings::SnapshotOffset ) {
			unsigned int f = 1 + (fid - Settings::SnapshotFirst) / Settings::SnapshotOffset;
			unsigned int nr = worker2d->get_nRanks();
			unsigned int me = worker2d->get_Rank();
			if ( (f % nr) == me ) { //##MK::solution so far not tested for mpiexec -n > 1
				//this work partitioning is very naive because the lower index fid datasets are orders of magnitude larger, therefore
				//better compute size of each dataset and distribute accordingly
				//MPIE::when working on the MPIE cluster though we utilize nRanks = 1 as the system is no lustre anyway and 
				//instead process each dataset with the maximum number of threads per process available...

				//MK::currently setup for mpiexec -n 1 i.e nRanks == 1 only! with export OMP_NUM_THREADS=8, export KMP_AFFINITY=verbose,scatter
				if ( worker2d->read_snapshot_contour( "GBContourPoints", fid ) == true ) {
					bool status = true;

					//make a container to store results local to the process in memory of main thread
					std::vector<wghtd_imcurv_res>* rbucket = NULL;
					rbucket = new vector<wghtd_imcurv_res>;
					worker2d->gresults_data2d.push_back( rbucket );
					worker2d->gresults_meta2d.push_back( imcurv_meta( fid, 0 ) );

					status = worker2d->integral_mean_curvature_approximation( fid );

					//implement results reporting
					if ( status == true ) {
						worker2d->report_wghtd_imcurv( pargv[SUFFIX], fid, true );
					}

					worker2d->reset_memory();
				}
			}
		}

	} //analyze at all

	delete worker2d;
	MPI_Finalize();
}


void FuseCapillaryActivity(  int pargc, char** pargv  )
{
//INITIALIZATION
	//MK::Init MPI capable of hybrid MPI/OpenMP even though the TopologyTracer currently at no point calls MPI commands within an OpenMP parallel region
	int supportlevel_desired = MPI_THREAD_FUNNELED;
	int supportlevel_provided;
	MPI_Init_thread(&pargc, &pargv, supportlevel_desired, &supportlevel_provided);

	int nranks, rank;
	if (supportlevel_provided < supportlevel_desired) {
		versioninfo(MASTER, SINGLE_PROCESS);
		cerr << "ERROR::CapillaryActivityTracking Insufficient threading capabilities of the MPI library!"; return;
	}
	else {
		MPI_Comm_size(MPI_COMM_WORLD, &nranks);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	}

	double timer = MPI_Wtime();
	versioninfo(rank, nranks);
	if (rank == MASTER)
		cout << "Starting auxiliary code to fuse capillary-activity data in " << Settings::Dimensionality << "d..." << endl;

	curvfuseHdlP worker = new curvfuseHdl;
	worker->set_MPICommunication(rank, nranks);
	worker->init_MPIDatatypes();

	//MPI_Barrier( MPI_COMM_WORLD );

//CHOOSE TARGETS for which at all to analyze
	if ( Settings::CapillaryTrackingTargets == E_CAP_ALL ) { //Settings are global to all processes
		if ( worker->init_targets() != true ) {
			cerr << "ERROR::Worker::CapillaryActivityTracker " << worker->get_Rank() << " target list was not compilable!" << endl;
			delete worker; worker = NULL; MPI_Finalize(); return;
		}
	}
	else {
		if (worker->read_targets_fromlist() != true) {
			cerr << "ERROR::Worker::CapillaryActivityTracker " << worker->get_Rank() << " list of TargetGrainIDs was nonexistent or errornous!" << endl;
			delete worker; worker = NULL; MPI_Finalize(); return;
		}
	}
	worker->myprofiler.logev("MPIStartup" + std::to_string(Settings::SimID), (double)(MPI_Wtime() - timer)); timer = MPI_Wtime();

//DATASET HANDLING
	//MK::work through capillary activity dump data by loading them, extracting relevant pieces of information and fusing into matrix with
	//rows as all targets and columns the timesteps

	unsigned int localhealth = 0;
	unsigned int globalhealth = 0;
	if ( Settings::CapillaryTrackingMode == E_CAP_HISTORY ) {
		if ( rank == MASTER ) //MK::only Master makes MPI I/O file objects existent on harddisks
			localhealth = worker->init_results_file();
		else
			localhealth = 0;
	}
	else { //...Mode == E_CAP_AVERAGE in which case all processes have to initialize their buffers
		localhealth = worker->init_localbuffer();
	}

	//meanwhile slaves wait for master prepping I/O files or prep local buffers in any case
	MPI_Allreduce(&localhealth, &globalhealth, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
	if ( Settings::CapillaryTrackingMode == E_CAP_HISTORY ) {
		if (globalhealth != 1) {
			cout << "ERROR::Worker::CapillaryActivityTracker " << worker->get_Rank() << " not all requested dump data exists!" << endl;
			//##MK::clear potential results file zombies generated by master
			delete worker; worker = NULL; MPI_Finalize(); return;
		}
	}
	else { //mode == E_CAP_AVERAGE
		if ( globalhealth != nranks ) {
			cout << "ERROR::Worker::CapillaryActivityTracker " << worker->get_Rank() << " not all processes have prepared their local buffer!" << endl;
			delete worker; worker = NULL; MPI_Finalize(); return;
		}
	}

	worker->myprofiler.logev("Initialization" + std::to_string(Settings::SimID), (double)(MPI_Wtime() - timer)); timer = MPI_Wtime();

//COOPERATIVE SUCCESSIVE TRACKING SEQUENTIALLY

	//apply filter criteria:: already done via reading list of target grainIDs
	//now we can process all analysis tasks, for this we load the heavy data per time step only, process, delete them, move to the next timestep
	//processes work in parallel on disjoint datasets, fusing everything in already existing file via writing to disjoint a priori known positions in main memory utilizing MPI I/O on MPI_COMM_SELF
	double tic = 0.0;
	double toc = 0.0;
	for (unsigned int fid = Settings::SnapshotFirst; fid <= Settings::SnapshotLast; fid += Settings::SnapshotOffset) {

		if (worker->workpartitioning_fid(fid) == true) {
			tic = MPI_Wtime();

			bool dbgood = worker->readDatasets(fid);

			toc = MPI_Wtime(); 
			worker->myprofiler.logev("LoadSnapshot." + std::to_string(fid), (toc - tic));

			//MK::load the heavy data for this via a database entry to enable reutilization of algorithms from the coorperatively parallel execution
			if (dbgood == true) {
				tic = MPI_Wtime();

				if ( Settings::CapillaryTrackingMode == E_CAP_HISTORY ) {
					worker->fuse_capillary_activity(fid);
				}

				if ( Settings::CapillaryTrackingMode == E_CAP_AVERAGE ) {
					worker->localaveraging_capillary_activity(fid);
				}

				toc = MPI_Wtime();
				worker->myprofiler.logev("FusedCapillaryActivity." + std::to_string(fid), (toc - tic));

				//if ( worker->get_Rank() == MASTER ) {
				std::cout << "...Worker " << worker->get_Rank() << " has successfully fused " << fid << std::endl;
				//}

				//MK::heavy data will be cleared prior to the loading of the next dataset or latest during destruction of worker class instance
			}
			else {
				cerr << "ERROR::Worker::CapillaryActivityTracking " << worker->get_Rank() << " dataset " << fid << " was not processable!" << endl;
			}
		}
		//else nothing, different process will take care of it!
	} //next snapshot

	//finally all results have been processed
	MPI_Barrier(MPI_COMM_WORLD);

	//only when coorperative analysis with multiple processes averaging requires communication via MPIReduce and master I/O spitting
	if ( Settings::CapillaryTrackingMode == E_CAP_AVERAGE ) {
		worker->spit_average_results();
	}

	worker->spit_profiling();

	delete worker; worker = NULL;
	MPI_Finalize();
}



int main(int argc, char** argv)
{
//PARAMETER, single_process I/O of XML file first...
	if ( argc < 2 ) { cerr << "ERROR::TopologyTracer2D3D cannot start without Parameter!" << endl; return 0; }
	unsigned int simid = atoi(argv[SIMID]);

	try { //to read user input
		if (argc > 1)	Settings::readXML(argv[PARAMETERFILE]);
		else			Settings::readXML();
	} 
	catch (exception& e) {
		cerr << "ERROR::TopologyTracer2D3D is unable to parse parameters file! Details:\n" << e.what() << endl << "Simulation will now halt." << endl; 
		return 0;
	}
	if ( Settings::checkUserInput() == false ) { cerr << "ERROR::TopologyTracer2D3D input failed the validity check!\n Please set realistic values into the XML file!" << endl; return 0; }
	Settings::SimID = simid;

//START THE ANALYSIS
	switch ( Settings::AnalysisMode )
	{
		case E_TRACKING_PARALLEL :
			TrackingInParallel( argc, argv, simid );
			break;
		case E_TRACKING_SEQUENTIAL :
			TrackingSequential( argc, argv, simid );
			break;
		case E_KNNCORRELATION :
			UtilityKNN( argc, argv, simid );
			break;
		case E_NUCSITES_LONGRANGE2D :
			LongRangeEnvPersistenceAnalysis2D( argc, argv );
			break;
		case E_NUCSITES_LONGRANGE3D :
			LongRangeEnvPersistenceAnalysis3D( argc, argv );
			break;
		case E_APPROX_CURVATURE_CONTOUR :
			ApproximateCurvatureFromContourData2D( argc, argv );
			break;
		case E_CAPILLARY_ACTIVITY :
			FuseCapillaryActivity( argc, argv );
			break;
		default : 
			cerr << "ERROR::TopologyTracer2D3D found no specific analysis mode chosen!" << endl;
			break;
	}
	return 0;
}
