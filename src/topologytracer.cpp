//MK::TopologyTracer is a software developed by Markus K{\u}hbach in 2015
//for questions and comments please contact markus.kuehbach at rwth-aachen.de, last major revision 20150720
//##MK::marks comments that where the code is in rather quick and dirty,i.e. parts of the code which have been validated operative as is but offer certainly room for further improvements


#include "topologytracer_kernel.h" 


using namespace std;

//input parameter
#define JOBID					1
#define NGRAINS_INITIALLY		2
#define ID_FIRSTDATASET			3
#define ID_OFFSETDATASET		4
#define ID_LASTDATASET			5
#define MEMORY_LIMIT_PERNODE	6


int main(int argc, char** argv) {

	topoHdlP worker = new topoHdl;

	MPI_Init (&argc, &argv); 
	MPI_Comm_rank (MPI_COMM_WORLD, &worker->myRank);
	MPI_Comm_size (MPI_COMM_WORLD, &worker->nRanks);
		worker->jobid = atoi( argv[JOBID] );
	worker->createMPIDataTypes();

	worker->N0 = atoi(argv[NGRAINS_INITIALLY]);
	worker->firstloaded = atoi(argv[ID_FIRSTDATASET]);
	worker->offsetloaded = atoi(argv[ID_OFFSETDATASET]);
	worker->lastloaded = atoi(argv[ID_LASTDATASET]);
	if ( worker->firstloaded < 0 || worker->lastloaded < 0 || worker->firstloaded > worker->lastloaded || worker->offsetloaded < 1 ) { cout << "Invalid input parameter!" << endl; MPI_Finalize(); delete worker; return 0; }
	worker->memlimitpernode = atof(argv[MEMORY_LIMIT_PERNODE]);
	if ( worker->memlimitpernode < MINIMUM_COMPLEXITY_PERNODE || worker->memlimitpernode > MAXIMUM_COMPLEXITY_PERNODE ) { cout << "Reconsider memory limits per node and number of nodes to perform the job!" << endl; MPI_Finalize(); delete worker; return 0; }

if ( worker->myRank == MASTER ) {
	cout << "TopologyTracer 2D/3D v0.1" << endl;
	cout << worker->myRank << "\t\t" << worker->N0 << ";" << worker->firstloaded << ";" << worker->offsetloaded << ";" << worker->lastloaded << endl;
}
/*
//Classical analysis with L.A. Barrales Mora 2D vertex code, loading data and plot evolution of surviving grains in area and shape space
//##MK::L.B. files at the moment only sequentially and format should be changed in the future
	if ( worker->nRanks > SEQUENTIAL ) { cout << "At the moment 2D vertex files only sequential!" << endl; MPI_Finalize(); delete worker; return 0; } 
	worker->read_rawdata_lbarrales_metadata_2d( 15000, "PopulationStatistics.dat" ); //##MK::15000 average number of grains during starting the simulation
	worker->read_rawdata_lbarrales_read_2d( "grainsMetrics.dat", 14029 ); //first/offset/lastloaded automatically!

	worker->analyse_grain2d_sizeevolution( 0, 1, 7706 );
	worker->analyse_grain2d_onegrain_vs_nearestnbors( 62458520, 0, 1, 7706 ); //largest grain of the population
	worker->analyse_grain2d_onegrain_vs_nearestnbors( 62491148, 0, 1, 7706 ); //initially growing then consumed due to competition
	worker->analyse_grain2d_onegrain_vs_nearestnbors( 62510972, 0, 1, 7706 ); //initially growing grain and then consumed

	MPI_Finalize();
	delete worker;
	return 0;
*/


//Classical analysis with 2D GraGrLeS - loading data and plot evolution of surviving grains in area and shape space
	worker->read_rawdata_cmiessen_read_2d( worker->N0, worker->firstloaded, worker->offsetloaded, worker->lastloaded );
	MPI_Barrier(MPI_COMM_WORLD);

	//worker->analyse_grain2d_sizeevolution( worker->firstloaded, worker->lastloaded );
	worker->analyse_grain2d_sizeevolution_nbormisori( worker->firstloaded, worker->lastloaded );

//MK::the following code serves as an example how the neighbors of a particular target grain can be identified
//take for instance worker->analyse_grain2d_onegrain_vs_nearestnbors( 25891, worker->firstloaded, worker->lastloaded );
//here grain ID 25891 is becoming tracked in the data interval [firstloaded, lastloaded]

	//the small
	MPI_Barrier(MPI_COMM_WORLD);	worker->analyse_grain2d_onegrain_vs_nearestnbors( 25891, worker->firstloaded, worker->lastloaded ); //CMiessen Grain 2608 becomes ID -1
	MPI_Barrier(MPI_COMM_WORLD);	worker->analyse_grain2d_onegrain_vs_nearestnbors( 26066, worker->firstloaded, worker->lastloaded ); //CMiessen Grain 2608 becomes ID -1
	MPI_Barrier(MPI_COMM_WORLD);	worker->analyse_grain2d_onegrain_vs_nearestnbors( 26150, worker->firstloaded, worker->lastloaded ); //CMiessen Grain 2608 becomes ID -1
	MPI_Barrier(MPI_COMM_WORLD);	worker->analyse_grain2d_onegrain_vs_nearestnbors( 27332, worker->firstloaded, worker->lastloaded ); //CMiessen Grain 2608 becomes ID -1
	MPI_Barrier(MPI_COMM_WORLD);	worker->analyse_grain2d_onegrain_vs_nearestnbors( 27189, worker->firstloaded, worker->lastloaded ); //CMiessen Grain 2608 becomes ID -1
	//the large
	MPI_Barrier(MPI_COMM_WORLD);	worker->analyse_grain2d_onegrain_vs_nearestnbors( 26161, worker->firstloaded, worker->lastloaded ); //CMiessen Grain 2608 becomes ID -1
	MPI_Barrier(MPI_COMM_WORLD);	worker->analyse_grain2d_onegrain_vs_nearestnbors( 25871, worker->firstloaded, worker->lastloaded ); //CMiessen Grain 2608 becomes ID -1
	MPI_Barrier(MPI_COMM_WORLD);	worker->analyse_grain2d_onegrain_vs_nearestnbors( 26103, worker->firstloaded, worker->lastloaded ); //CMiessen Grain 2608 becomes ID -1
	MPI_Barrier(MPI_COMM_WORLD);	worker->analyse_grain2d_onegrain_vs_nearestnbors( 26089, worker->firstloaded, worker->lastloaded ); //CMiessen Grain 2608 becomes ID -1
	MPI_Barrier(MPI_COMM_WORLD);	worker->analyse_grain2d_onegrain_vs_nearestnbors( 25843, worker->firstloaded, worker->lastloaded ); //CMiessen Grain 2608 becomes ID -1

//MK::safely shutdown the MPI
	MPI_Finalize();
	delete worker; 
	return 0;


//Classical analysis of Aboav-Weaire, Lewis law
	int DidEveryoneConstructTheKShells = worker->analyse_grain2d_construct_kshell(4, worker->firstloaded, worker->lastloaded);
cout << worker->myRank << " this rank has DidEveryoneConstructTheKShells." << endl;
	int Decision = 0;
	MPI_Allreduce ( &DidEveryoneConstructTheKShells, &Decision, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );

	if ( Decision == worker->nRanks ) { 
cout << worker->myRank << " all constructed their kshells, so analyse..." << endl;

		worker->analyse_grain2d_kshellevolution(4, worker->firstloaded, worker->lastloaded);
	}
	else { 
		cout << "KShell construction in all workers not successful." << endl; 
	}

	MPI_Finalize();
	delete worker;
	return 0;
}




/*
//##MK::DEBUG CODE

//how to get the file size in byte without opening, ##MK::mind when file exceeds 2^32 bytes
struct stat buf;
if ( stat( "Texture_0.ori", &buf ) != -1 ) {
	cout << "Int=" << buf.st_size << " Byte" << endl;
}
else {
	cout << "Not able to obtain the file size!" << endl;
}
return 0;

//how to call the orientation library with two Bunge Euler angle triplets
double eA[3], eB[3];
eA[0] = 0.0 / 180.0 * _PI_;
eA[1] = 0.0 / 180.0 * _PI_;
eA[2] = 0.0 / 180.0 * _PI_;
eB[0] = 56.9 / 180.0 * _PI_;
eB[1] = 30.7 / 180.0 * _PI_;
eB[2] = 59.1 / 180.0 * _PI_;
cout << worker->misorientationCubic( eA[0], eA[1], eA[2], eB[0], eB[1], eB[2] ) << endl;

//debug code to test the functionality of sending MPI_UNSIGNED
uint32_t a = 28000000;
uint32_t b = UINT32T_MAX;
uint32_t* buf = new uint32_t[2];
if ( worker->myRank == MASTER ) {
	buf[0] = a;
	buf[1] = b;
	MPI_Send( buf, 2, MPI_UNSIGNED, 1, 88, MPI_COMM_WORLD );
	cout << worker->myRank << " a = " << a << "\t\tb = " << b << endl;
}
if ( worker->myRank == 1 ) {
	MPI_Recv( buf, 2, MPI_UNSIGNED, MASTER, 88, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
	cout << worker->myRank << " a = " << buf[0] << "\t\tb = " << buf[1] << endl;
}
delete [] buf;
MPI_Finalize();
delete worker;
return 0;

*/