//MK::TopologyTracer is a software developed by Markus K{\u}hbach in 2015
//for questions and comments please contact markus.kuehbach at rwth-aachen.de, last major revision 20150720
//##MK::marks comments that where the code is in rather quick and dirty,i.e. parts of the code which have been validated operative as is but offer certainly room for further improvements


#include "topologytracer_kernel.h"


inline bool SortTimeStepAscending (const MPI_IO_GrainInfo2D &t1 , const MPI_IO_GrainInfo2D &t2) {
	return t1.fid < t2.fid;
}


topoHdl::topoHdl()
{
	jobid = MASTER;
	N0 = 0;
	nRanks = SEQUENTIAL;
	myRank = MASTER;

	firstloaded = 0;
	offsetloaded = 1;
	lastloaded = 0;
	memlimitpernode = 0.0;

	//default complexity model
	complexitymodel_pergrain = sizeof(grain2d) + ( AVERAGE_NEIGHBORS * sizeof(adjointb2d) );
	complexitymodel_a = 1;
	complexitymodel_t0 = 0;
	complexitymodel_q = 1;


	kshellmax = DEFAULT_KSHELL;
	nneighborsmax = DEFAULT_NNEIGHBORS;

	myMemGuard = 0.0;
	dimensionality = ANALYZE_IN_2D;
}


topoHdl::~topoHdl()
{
	//in all mydata free memory for grains
	for (uint32_t d = 0; d < mydata.size(); d++) {
		grain2dP grainbucket = mydata[d];

		for ( uint32_t gr = 0; gr < mydatasize[d]; gr++ ) {
			delete [] grainbucket[gr].neighbors2d;

			if ( grainbucket[gr].kshell != NULL ) {
				delete [] grainbucket[gr].kshell->n;
				delete [] grainbucket[gr].kshell->ids;
				//delete [] grainbucket[gr].kshell->disori;
			}
			delete [] grainbucket[gr].kshell;
		}

		delete [] grainbucket;

		mydata[d] = NULL;
		mydatasize[d] = 0;
	}

	for (uint32_t d = 0; d < gen_aboav_weaire_avfaces.size(); d++) {
		delete [] gen_aboav_weaire_avfaces[d];
	}

	for (uint32_t d = 0; d < gen_aboav_weaire_nnbors.size(); d++) {
		delete [] gen_aboav_weaire_nnbors[d];
	}

	for (uint32_t d = 0; d < this->disori_environment.size(); d++) {
		delete [] disori_environment[d];
	}
}


void topoHdl::init( void )
{

}


uint32_t topoHdl::get_closest_idealori( double * quat )
{
	//##MK::a dummy which is not yet implemented
	//##MK::I randomly oriented grains, i.e. those how are not close to any particular probed position in Euler space as RANDOM_ORIENTATION
	return 0;
}


uint32_t topoHdl::addOrientation ( double * bunge, uint32_t m, bool recategorize )
{
	//recategorize == true --> incoming orientations are checked for agreement with already known orientation
	//if not orientation is always add to the pool otherwise added
	int closestid = NOT_WITHIN_GRIDRESU; //MK::smaller than range of uint32_t so can be a valid marker
	double disori = 2 * _PI_;
	double closestdisori = RESOLUTION_SO3GRID;
	double qbunge[4], qcand[4];

	euler2quaternion( bunge, qbunge );

	if ( recategorize == true ) {
		//check disorientation to all other components in oripool
		for ( uint32_t cand = 0; cand < oripool.size(); cand++ ) {
			qcand[0] = oripool[cand].q0;
			qcand[1] = oripool[cand].q1;
			qcand[2] = oripool[cand].q2;
			qcand[3] = oripool[cand].q3;
			disori = misorientationCubicQxQ( qbunge[0], qbunge[1], qbunge[2], qbunge[3], qcand[0], qcand[1], qcand[2], qcand[3]);

			//closer to any other orientation that is already known?
			if (disori <= closestdisori) {
				closestdisori = disori;
				closestid = cand;
			}
			//MK::however still testing against all other to get the closest of the known components
		}

		if ( closestid != NOT_WITHIN_GRIDRESU ) {
			return closestid;
		}
	}

	//hasnt returned yet so either recategorization == false || all cand oris too distant apart in SO3 -> in either case add that orientation
	struct ori anori;
	anori.bunge1 = bunge[0];
	anori.bunge2 = bunge[1];
	anori.bunge3 = bunge[2];
	anori.q0 = qbunge[0];
	anori.q1 = qbunge[1];
	anori.q2 = qbunge[2];
	anori.q3 = qbunge[3];
	anori.id = m;
	anori.closest = get_closest_idealori( qbunge ); //##MK::needs implementation

	oripool.push_back( anori );

//cout << "ADD\t\t" << (myoripool.size() - 1) << "\t\t" << anori.bunge1 << "\t\t" << anori.bunge2 << "\t\t" << anori.bunge3 << "\t\t" << anori.q0 << "\t\t" << anori.q1 << "\t\t" << anori.q2 << "\t\t" << anori.q3 << endl;
	return (oripool.size() - 1);
}


inline double topoHdl::complexitymodel( uint32_t t )
{
	double elapsedtime = t;
	double comp = complexitymodel_pergrain * complexitymodel_a * pow( ( elapsedtime - complexitymodel_t0 ), complexitymodel_q);

	return comp;
}


inline double topoHdl::complexitymodel_linear ( uint32_t N0, uint32_t Nfin, uint32_t istart, uint32_t i, uint32_t iend )
{
	double comp = complexitymodel(i);

	//##MK::very simple load model where it is assumed that the number of grains decreases with constant rate over time
	if ( dimensionality == ANALYZE_IN_2D )
		comp = TYPICALSIZE_2D_GRAIN * ( Nfin + ((double) (N0 - Nfin)) * ( ((double) (i - istart)) / ((double) (iend - istart))) );
	if ( dimensionality == ANALYZE_IN_3D )
		comp = TYPICALSIZE_3D_GRAIN * ( Nfin + ((double) (N0 - Nfin)) * ( ((double) (i - istart)) / ((double) (iend - istart))) );

	return comp;
}


inline double topoHdl::complexitymodel_linear ( uint32_t ngr )
{
	double comp = 0.0;

	//##MK::very simple load model where it is assumed that the number of grains decreases with constant rate over time
	if ( dimensionality == ANALYZE_IN_2D )
		comp = TYPICALSIZE_2D_GRAIN * ngr;
	if ( dimensionality == ANALYZE_IN_3D )
		comp = TYPICALSIZE_3D_GRAIN * ngr;

	return comp;
}


void topoHdl::createMPIDataTypes( void )
{
	int elementCountsDomain1[2] = {2, 2};
	MPI_Aint displacementDomain1[2] = { 0, 8 };
	MPI_Datatype oldTypesDomain1[2] = { MPI_UNSIGNED, MPI_DOUBLE };
	MPI_Type_create_struct(2, elementCountsDomain1, displacementDomain1, oldTypesDomain1, &MPI_IO_GrainInfo2D_Type);

	MPI_Type_commit(&MPI_IO_GrainInfo2D_Type);


	int elementCountsDomain2[2] = {2, 1};
	MPI_Aint displacementDomain2[2] = { 0, 8 };
	MPI_Datatype oldTypesDomain2[2] = { MPI_UNSIGNED, MPI_DOUBLE };
	MPI_Type_create_struct(2, elementCountsDomain2, displacementDomain2, oldTypesDomain2, &MPI_IO_MetadataInfo_Type);

	MPI_Type_commit(&MPI_IO_MetadataInfo_Type);

	//##MK::add further here
}


bool topoHdl::read_rawdata_cmiessen_read_2d( uint32_t N0, uint32_t first, uint32_t offset, uint32_t last)
{
	//knowing now the complexity model of form N(t) = k*t^q the MPI nodes
	//fill all in parallel the assignment tableau

	//check reasonable constraints
	if (N0 <= 0) { cout << "ERR:Invalid number of initial grains!" << endl; return LOADING_DATASET_FAILED; }
	if ( (first % offset) != 0 || (last % offset) != 0 ) { cout << "ERR:Invalid first and last according to modulo offset scheme!" << endl; return LOADING_DATASET_FAILED; }
	if ( first < firstloaded || offset != offsetloaded || last > lastloaded || first > last || offset < 1 || last < 0 ) { cout << "Other error " << endl; return LOADING_DATASET_FAILED; }


	//all participate in job distribution
	//collect round robin complexity < MAXIMUM_COMPLEXITY_PERNODE, ##MK::CMiessen utilized fixed timestep logging at the moment
	stringstream dynfname;
	uint32_t howmanycollected = 0;
	int targetnode = MASTER;
	double cplexaccumulated = 0.0;
	//MK::in the initial version work was partitioned by fid % nRanks, this is not optimal because sequential time steps overwhich to extract data are spread on the nodes
	//MK::Hence now, fill available ranks up completely and test for message size

	for ( uint32_t f = first; f <= last; f += offset ) {

		uint32_t fid = (f/offsetloaded) - (firstloaded/offsetloaded);

		//int targetnode = fid % nRanks; //MK::works nRanks > 0

		//determine the expected computational complexity of the job
		double cplex = complexitymodel_linear ( N0, 0, (first/offset), (f/offset), (last/offset) );

		char status[2] = { TARGETNODE_REJECTS, DATA_INVALID };

		while ( status[0] == TARGETNODE_REJECTS && targetnode < nRanks ) {

			if ( myRank == targetnode ) {
				dynfname.str(std::string());
				dynfname << "Texture_" << f << ".ori";

				//system call, how large is the file?
				struct stat buf;
				if ( stat( dynfname.str().c_str() , &buf ) != -1 ) {
					cplex = buf.st_size;
					status[1] = DATA_VALID;
				}
				//else {status[1] = DATA_INVALID;} //already initialized invalidated

				//fitting still in my local memory?
				if ( cplexaccumulated < this->memlimitpernode ) { 
					status[0] = TARGETNODE_TAKES;
					cplexaccumulated += cplex;
				}

cout << "MyRank;File;Size/Bytes\t\t" << myRank << ";" << f << ";" << (myMemGuard + cplex) << endl;
			}
			//broadcast so that all can properly exit or participate in redistribution
			MPI_Bcast( status, 2, MPI_CHAR, targetnode, MPI_COMM_WORLD );

			if ( status[1] == DATA_INVALID ) { 
				cout << "Analyzing file " << f << " was unsuccessful!" << endl;
				return false;
			}


			if ( status[0] == TARGETNODE_TAKES ) 
				break;

			//TARGETNODE_REJECTED, proceed with next available node
			targetnode++;
		}


		if ( targetnode >= nRanks ) {
			cout << "Insufficient amount of node-local memory to read in all data!" << endl; return LOADING_DATASET_FAILED;
		}

		//now targetnode took it guaranteed!

		//MK::constructing the assignment array linear from 0, to tlast is essential
		struct metadata fi;
		fi.globalid = fid;
		fi.localid = UNKNOWN;
		if ( myRank == targetnode ) {
			fi.localid = howmanycollected;
			howmanycollected++;
		}
		fi.whichnode = targetnode;
		fi.ngrainsintotal = UNKNOWN;
		fi.complexity = UNKNOWN;

		assignment.push_back ( fi );
	}

cout << myRank << " finished distributing the work, now workers read in data in parallel ..." << endl;

	//now read in sequential all my files, workers find on what they have to work automatically
	for ( uint32_t f = first; f <= last; f += offset ) {

		uint32_t fid = (f/offsetloaded) - (firstloaded/offsetloaded);
		int targetnode = assignment[fid].whichnode;

		if ( myRank == targetnode ) {
			dynfname.str(std::string());
			dynfname << "Texture_" << f << ".ori";

			//uint32_t nelem[1] = { 0 };
			//double complex[1] = { 0.0 };
			bool success = false;

			success = read_rawdata_cmiessen_interpret_2d( dynfname.str().c_str(), fid ); //, nelem, complex ); //##DEBUG&assignment[fid] );

//cout << "MASTER scale = " << assignment[fid].localid << ";" << assignment[fid].ngrainsintotal << ";" << assignment[fid].complexity << "\n";

		}
	}

cout << myRank << " finished successfully reading its input files, now committing meta information" << endl;

	MPI_Barrier( MPI_COMM_WORLD );

	//redistribute information after reading
	for ( uint32_t f = first; f <= last; f += offset ) {

		uint32_t fid = (f/offsetloaded) - (firstloaded/offsetloaded);
		int targetnode = assignment[fid].whichnode;

		//stack array because MPI_Bcast requires pointer
		MPI_IO_MetadataInfo mdata[1];
		mdata[0].lid = NO_PROPERTY_ASSIGNED;
		mdata[0].ngrainsintotal = 0;
		mdata[0].complexity = 0.0;

		if ( myRank == targetnode ) {
			mdata[0].lid = assignment[fid].localid;
			mdata[0].ngrainsintotal = assignment[fid].ngrainsintotal;
			mdata[0].complexity = assignment[fid].complexity;
		}
		//have not all gone out?
		MPI_Bcast( &mdata, 1, MPI_IO_MetadataInfo_Type, targetnode, MPI_COMM_WORLD );

		assignment[fid].localid = mdata[0].lid;
		assignment[fid].ngrainsintotal = mdata[0].ngrainsintotal;
		assignment[fid].complexity = mdata[0].complexity;
//cout << "After BCast MASTER scale = " << assignment[fid].localid << ";" << assignment[fid].ngrainsintotal << ";" << assignment[fid].complexity << "\n";
	}

	if ( myRank == MASTER ) {
		stringstream logassgn_fn;
		logassgn_fn << "JobAssign.MASTER." << N0 << ".First." << first << ".Offset." << offset << ".Last." << last << ".csv";

		ofstream logassgn;
		logassgn.open( logassgn_fn.str().c_str() );
		logassgn << "fid;target;globalid;localid;ngrainsintotal;complexity(byte)\n";
		uint32_t fid;
		for (uint32_t f = first; f <= last; f += offset ) {
			fid = (f/offsetloaded) - (firstloaded/offsetloaded);
			logassgn << fid << ";" << assignment[fid].whichnode << ";" << assignment[fid].globalid << ";" << assignment[fid].localid << ";" << assignment[fid].ngrainsintotal << ";" << assignment[fid].complexity << "\n";
		}
		logassgn.flush();
		logassgn.close();
	}

cout << myRank << " information was redistributed successfully." << endl;

	return true;
}


bool topoHdl::read_rawdata_cmiessen_interpret_2d(const char* fname, uint32_t fid ) // uint32_t* nread, double* complexity) //##DEBUG, struct metadata * finfo)
{
	//get memory to store the data
	uint32_t grsize = CMIESSEN_ORIFILE_GRAINPRECACHE;
	grain2dP grbucket = NULL;
	grbucket = new grain2d[grsize];
	if ( grbucket == NULL ) { cout << "ERR::Memory allocation error in read_rawdata_cmiessen_interpret!" << endl; return false; }
	uint32_t grfilledin = 0;

	ifstream orifile;
	string oriline;
	istringstream line;
	string datapiece;

	long szofg2d = sizeof(grain2d);
	long szofadjointb2d = sizeof(adjointb2d);
//cout << szofg2d << "\t\t" << szofadjointb2d << endl;

	orifile.open( fname ); //##DEBUG"Texture_0.ori");
	if ( orifile.is_open() == true ) {
cout << myRank << "-th, " << fname << " opened ..." << endl;
		uint32_t nlines = 0; //CMiessen each line a disjoint grain
		double totalsize = 0.0;
		long s = 0;

		while ( orifile.good() == true ) {
			if ( nlines >= MAXIMUM_NUMBER_OFGRAINS ) {
				cout << "The dataset is too large but the reader requires all files to become read, so abort!" << endl; return false;
			}

			//CMiessen 2D GraGrLeS files have no header line
			getline( orifile, oriline );
 
			s = oriline.size(); //get true size complexity of the line

//cout << s << "\t\t" << nlines << endl;

			if ( s > 0 ) { //MK::okay there is content that is safe to be assumed a valid formatted line

				//enough memory in data containerleft?
				if (grfilledin >= grsize) { //add further memory
					grain2dP tmpgr = NULL;
					tmpgr = new grain2d[grfilledin + CMIESSEN_ORIFILE_GRAINPRECACHE]; //##MK::if the number of grains were known there is optimization potential here
					if ( tmpgr == NULL ) { cout << "ERR::Memory allocation error in read_rawdata_cmiessen_interpret!" << endl; return false; }

					for (uint32_t g = 0; g < grfilledin; g++) {
						tmpgr[g].globalid = grbucket[g].globalid;
						tmpgr[g].nfaces = grbucket[g].nfaces;
						tmpgr[g].boundarycontact = grbucket[g].boundarycontact;
						tmpgr[g].area = grbucket[g].area;
						tmpgr[g].areachange = grbucket[g].areachange;
						tmpgr[g].perimeter = grbucket[g].perimeter;
						tmpgr[g].totalboundaryenergy = grbucket[g].totalboundaryenergy;
						tmpgr[g].ori = grbucket[g].ori;
						tmpgr[g].neighbors2d = grbucket[g].neighbors2d;
						tmpgr[g].kshell = grbucket[g].kshell; //##DEBUGwasNULL;
					}

					delete [] grbucket;
					grbucket = tmpgr;
					grsize = grsize + CMIESSEN_ORIFILE_GRAINPRECACHE;
				}

				istringstream line( oriline );
				double oobunge[3];
				double markerid = 0; //dummy value

				getline( line, datapiece, '\t');		grbucket[nlines].globalid = (atoi( datapiece.c_str() )); // - 1); //because CMiessen indices from 1 to N rather than 0 to N-1
				getline( line, datapiece, '\t');		grbucket[nlines].nfaces = atoi( datapiece.c_str() );
				getline( line, datapiece, '\t');		grbucket[nlines].boundarycontact = atoi( datapiece.c_str() );
				getline( line, datapiece, '\t');		grbucket[nlines].area = atof( datapiece.c_str() );
				getline( line, datapiece, '\t');		grbucket[nlines].areachange = atof( datapiece.c_str() );
				getline( line, datapiece, '\t');		grbucket[nlines].perimeter = atof( datapiece.c_str() );
				getline( line, datapiece, '\t');		grbucket[nlines].totalboundaryenergy = atof( datapiece.c_str() );
				getline( line, datapiece, '\t');		oobunge[0] = atof( datapiece.c_str() );
				getline( line, datapiece, '\t');		oobunge[1] = atof( datapiece.c_str() );
				getline( line, datapiece, '\t');		oobunge[2] = atof( datapiece.c_str() );

				grbucket[nlines].ori = addOrientation( oobunge, markerid, true );

if ( myRank == MASTER && (nlines % 10000 == 0) ) cout << nlines << endl;
//if ( grbucket[nlines].nfaces < 3 ) cout << grbucket[nlines].nfaces << "\t" << grbucket[nlines].globalid << endl;

				grbucket[nlines].kshell = NULL;
				grbucket[nlines].neighbors2d = NULL;
				grbucket[nlines].neighbors2d = new adjointb2d[grbucket[nlines].nfaces];
				if (grbucket[nlines].neighbors2d == NULL) { cout << "ERR::Memory allocation error reading line " << nlines << "**cmiessen_interpret!" << endl; return false; }

				for (uint32_t nf = 0; nf < grbucket[nlines].nfaces; nf++) {
					getline( line, datapiece, '\t');	grbucket[nlines].neighbors2d[nf].globalid = atoi( datapiece.c_str() ); // - 1 ); //because CMiessen indices from 1 to N rather than 0 to N-1
					getline( line, datapiece, '\t');	grbucket[nlines].neighbors2d[nf].length = atof( datapiece.c_str() );
					//cout << "\t\t" << grbucket[nlines].neighbors2d[nf].globalid << " | " << grbucket[nlines].neighbors2d[nf].length << endl;
					getline( line, datapiece, '\t');
				}
				
//cout << grbucket[nlines].globalid << "\t" << grbucket[nlines].area << "\t" << grbucket[nlines].boundarycontact << endl;

				grfilledin++;
//cout << fid << ";" << totalsize << ";" << grbucket[nlines].nfaces << endl;
				totalsize = totalsize + (szofg2d + (grbucket[nlines].nfaces * szofadjointb2d));
				nlines++;
				
			}
		}

//cout << "Interpreting-------" << nlines << endl;
		//nread[0] = nlines;
		//complexity[0] = totalsize;

		//all data copied, thin out precached memory
		grain2dP fingr = NULL;
		fingr = new grain2d[nlines];
		if ( fingr == NULL ) { cout << "ERR::Memory allocation error in final interpret read_rawdata_cmiessen_interpret!" << endl; return false; }

		for (uint32_t g = 0; g < nlines; g++) {
			fingr[g].globalid = grbucket[g].globalid;
			fingr[g].nfaces = grbucket[g].nfaces;
			fingr[g].boundarycontact = grbucket[g].boundarycontact;
			fingr[g].area = grbucket[g].area;
			fingr[g].areachange = grbucket[g].areachange;
			fingr[g].perimeter = grbucket[g].perimeter;
			fingr[g].totalboundaryenergy = grbucket[g].totalboundaryenergy;
			fingr[g].ori = grbucket[g].ori;
			fingr[g].neighbors2d = grbucket[g].neighbors2d;
			fingr[g].kshell = grbucket[g].kshell;
		}

		delete [] grbucket;
		grbucket = NULL;
		grsize = 0;

		mydata.push_back( fingr );
		mydatasets.push_back( fid );
		mydatasize.push_back( nlines );
		myMemGuard = myMemGuard + totalsize;

		//fill in metadata that are only known to me, the local node when the file was successfully read
		assignment[fid].localid = mydatasets.size() - 1;
		assignment[fid].ngrainsintotal = nlines;
		assignment[fid].complexity = totalsize;

		orifile.close();

cout << fname << " processed and closed ..." << endl;

		return true;
	}
	else {
		cout << "ERROR::Rawdata file not found!" << endl; return false;
	}
}


bool topoHdl::read_rawdata_lbarrales_metadata_2d( uint32_t N0, const char* ensembleInfo )
{
	//read information how many grains to expect from population file that
	//first read the header
	ifstream ensinfofile;
	string ensinfoline;
	istringstream line;
	string datapiece;

	ensinfofile.open( ensembleInfo );
	if ( ensinfofile.is_open() == true ) {
cout << ensembleInfo << " opened ..." << endl;

		//kick one header line
		getline( ensinfofile, ensinfoline );

		uint32_t f = 0;
		double expectedComplexity = 0.0;

		while ( ensinfofile.good() == true ) {
			//read information and partition work, currently everything on MASTER

			getline( ensinfofile, ensinfoline );
			istringstream line( ensinfoline );

			if ( ensinfoline.size() < 6 )
				continue;

			//okay, seems reasonable snippet of data as I am expecting at least six pieces of information
			getline( line, datapiece, '\t');		double tt = atof( datapiece.c_str() );
			getline( line, datapiece, '\t');		uint32_t ngrains = atoi( datapiece.c_str() );

			//MK::ONLY ONE NODE
			struct metadata fi;
			fi.globalid = f;
			fi.localid = f;
			fi.whichnode = MASTER;
			fi.ngrainsintotal = ngrains;
			fi.complexity = complexitymodel_linear( ngrains );
			expectedComplexity = expectedComplexity + fi.complexity;

			if ( expectedComplexity >= this->memlimitpernode ) {
				cout << "It is very likely that there is not enough memory on the single node to execute the job so quit!" << endl; 
				return false;
			}

			assignment.push_back( fi );


			f = f + 1;
		}

		//##DEBUG::current
		this->firstloaded = 0;
		this->offsetloaded = 1;
		this->lastloaded = 7706; //f - 1;

		ensinfofile.close();

cout << ensembleInfo << " processed and closed ... now I have " << this->firstloaded << " with offset " << this->offsetloaded << " up to last " << this->lastloaded << endl;
		return true;
	}
	else {
		cout << "ERROR::Population information file not found!" << endl; 
		return false;
	}
}


bool topoHdl::read_rawdata_lbarrales_read_2d( const char* dfname, uint32_t ngrains )
{
	//all read the metrics file, and read the appropriate section if they host the data, 
	//##MK::clearly reading from individual files would be way better!
	//with the assignment we know the data topology ||headerline||noris|||headerline||assignment[0].ngrainsintotal lines||assignment[0+1]||...
	ifstream grmetricsfile;
	string grmetricsline;
	istringstream line;
	string datapiece;

	long szofg2d = sizeof(grain2d);
	long szofadjointb2d = sizeof(adjointb2d);


	grmetricsfile.open( dfname );
	if ( grmetricsfile.is_open() == true ) {
cout << dfname << " opened ..." << endl;

		//L.A. Barrales 2D vertex data have a first header line ot become kicked
		getline( grmetricsfile, grmetricsline );

		//MK::the output format is more complex, first all grains that ever participate are listed, then a continuous block of lines
		//separating the time steps without blank lines are listed...

		//read first all grains that somewhen participate
		double oobunge[3];
		uint32_t markerid;

		for ( uint32_t o = 0; o < ngrains; o++ ) {
			if ( grmetricsfile.good() == false ) { cout << "GrainMetrics file error during header reading of grain " << o << endl; return false; }

			getline( grmetricsfile, grmetricsline ); //format is ID, phi1, PHI, phi2
			istringstream line( grmetricsline );

			getline( line, datapiece, '\t');		markerid = atoi( datapiece.c_str() );
			getline( line, datapiece, '\t');		oobunge[0] = atof( datapiece.c_str() ) / 180.0 * _PI_;
			getline( line, datapiece, '\t');		oobunge[1] = atof( datapiece.c_str() ) / 180.0 * _PI_;
			getline( line, datapiece, '\t');		oobunge[2] = atof( datapiece.c_str() ) / 180.0 * _PI_;

			uint32_t addoid = addOrientation( oobunge, markerid, false ); //no categorization because Luis file has first all grain properties assuming them static!
		}
		cout << oripool.size() << " orientations imported successfully from dfname." << endl;

		//kick last header before the actual data
		getline ( grmetricsfile, grmetricsline );

		//now start interpreting the dataset
		//##MK::Barrier necessary if done cooperatively MPI_Barrier ( MPI_COMM_WORLD );

		for ( uint32_t f = 0; f < assignment.size(); f++) {
			//##MK::currently only working with one node so all data go into the MASTER
			unsigned long nlines = 0;
			long totalsize = 0;
			long s = 0;

			//assignment[f].ngrainsintotal = assignment[f].ngrainsintotal - 1; //##workaround 
			uint32_t ngr_expected = assignment[f].ngrainsintotal;

			//allocate memory for the grains
			grain2dP grbucket = NULL;
			grbucket = new grain2d[ngr_expected];
			if ( grbucket == NULL ) { cout << "ERR::Memory allocation error in read_rawdata_cmiessen_interpret!" << endl; return false; } //expected to become filled completely...

			for ( uint32_t g = 0; g < ngr_expected; g++) {
				if ( grmetricsfile.good() == false ) { cout << "File error in reading data at f/g = " << f << "\t" << g << endl; return false; }

				getline( grmetricsfile, grmetricsline);
 
				s = grmetricsline.size(); //get true size complexity of the line
 
				if ( s > 0 ) { //MK::okay there is content assume a valid formatted line
					//idx			time			area			perimeter		dv-v			Ngb	ListOfNeighbours	angle	kb_curvature	special

					istringstream line( grmetricsline );

					getline( line, datapiece, '\t');	grbucket[g].globalid = atoi( datapiece.c_str() );
					getline( line, datapiece, '\t');	double tt = atof( datapiece.c_str() );
					getline( line, datapiece, '\t');	grbucket[g].area = atof( datapiece.c_str() );
					getline( line, datapiece, '\t');	grbucket[g].perimeter = atof( datapiece.c_str() );
					getline( line, datapiece, '\t');	grbucket[g].areachange = atof( datapiece.c_str() );
					getline( line, datapiece, '\t');	grbucket[g].nfaces = atoi( datapiece.c_str() );
					//fill in other values
					grbucket[g].boundarycontact = NOTCHECKED;
					grbucket[g].totalboundaryenergy = NOTANALYZED;

//##DEBUG
if( g == 0) cout << f << ";" << tt << endl;
//if ( g == 0 ) cout << "START--f=" << f << "---tt=" << tt << endl;
//if ( g == (ngr_expected - 1) ) cout << "END--f=" << f << "---tt=" << tt << endl;

					uint32_t oripoolid = getOripoolID( grbucket[g].globalid );
					if ( oripoolid == UNKNOWN ) { cout << "Failed in f/g " << f << "\t" << g << " UNKNOWN ORI!" << endl; return false; }
					grbucket[g].ori = oripoolid;

					//analyze neighbors
					grbucket[g].kshell = NULL;
					grbucket[g].neighbors2d = NULL;
					grbucket[g].neighbors2d = new adjointb2d[grbucket[g].nfaces];
					if (grbucket[g].neighbors2d == NULL) { cout << "ERR::Memory allocation error reading line " << nlines << "**cmiessen_interpret!" << endl; return false; }

					//import neighbor ids
					for (uint32_t nf = 0; nf < grbucket[nlines].nfaces; nf++) {
						getline( line, datapiece, '\t');
						grbucket[g].neighbors2d[nf].globalid = atoi( datapiece.c_str() );
					}
				}

				totalsize = totalsize + (szofg2d + grbucket[g].nfaces * szofadjointb2d);
			} //next grain please

			mydata.push_back( grbucket );
			mydatasets.push_back( f );
			mydatasize.push_back( nlines );

			assignment[f].localid = mydatasets.size() - 1;
//##DEBUG::cout << "mydatasetsandSize" << mydatasets[f] << ";" << mydatasize[f] << endl;
			myMemGuard = myMemGuard + totalsize;

			//##DEBUGif ( f > 50 ) break;

//##DEBUGcout << "Assignment\t" << f << "\tanalyzed." << endl;
		} //check next assignment

		grmetricsfile.close();
		cout << dfname << " processed and closed ..." << endl;
		return true;
	}

	else {
		cout << "ERROR::Rawdata file not found!" << endl; return false; 
	}
}


uint32_t topoHdl::getOripoolID( uint32_t globalid )
{
	//##MK::scan for this data element as number might not be ordered and not contiguous, later ordering might be a helpful precomputation
	uint32_t id = UNKNOWN;
	uint32_t ncand = oripool.size();

	for ( uint32_t cand = 0; cand < ncand; cand++) {
		if ( oripool[cand].id == globalid ) {
			id = cand;
			break;
		}
	}
	return id;
}


uint32_t topoHdl::getGlobalID( uint32_t globalid, uint32_t whichmydataset)
{
	//##MK::scan for this data element as number might not be ordered and not contiguos, later ordering might be a helpful precomputation
	uint32_t id = UNKNOWN;
	uint32_t ncand = this->mydatasize[whichmydataset];

	grain2dP thegrains = mydata[whichmydataset];

	for ( uint32_t cand = 0; cand < ncand; cand++) {
		if ( thegrains[cand].globalid == globalid ) {
			id = thegrains[cand].globalid;
			break;
		}
	}
	return id;
}


uint32_t topoHdl::getNumberOfFaces( uint32_t globalid, uint32_t whichmydataset)
{
	//##MK::scan for this data element as number might not be ordered and not contiguos, later ordering might be a helpful precomputation
	uint32_t nfaces = UNKNOWN;
	uint32_t ncand = this->mydatasize[whichmydataset];

	grain2dP thegrains = mydata[whichmydataset];

	for ( uint32_t cand = 0; cand < ncand; cand++) {
		if ( thegrains[cand].globalid == globalid ) {
			nfaces = thegrains[cand].nfaces;
			break;
		}
	}

	return nfaces;
}


double topoHdl::getGrainVolume( uint32_t globalid, uint32_t whichmydataset)
{
	//##MK::scan for this data element as number might not be ordered and not contiguos, later ordering might be a helpful precomputation
	double volume = 0.0;
	uint32_t ncand = this->mydatasize[whichmydataset];

	for ( uint32_t cand = 0; cand < ncand; cand++) {
		if ( mydata[whichmydataset][cand].globalid == globalid ) {
			volume = mydata[whichmydataset][cand].area;
			break;
		}
	}

	return volume;
}


#define LAGB2HAGB			0.261799387
#define CUTOFF_ENERGY		0.017453292

inline double topoHdl::getMobilityTimesEnergy( double disori )
{
	//mobility model along the line of Humphreys, Rollett
	double mobility = 1.0 - ( 1.0 * exp ( - 5.0 * pow( (disori/LAGB2HAGB), 4.0 ) ));

	//energy model C. Miessen
	double gamma_hagb = 1.0;
	double energy = gamma_hagb;

	if ( disori <= LAGB2HAGB ) {
		//all energies of boundaries with disori <= 1.0/180.0*_PI_ are clamped to that value at 1.0deg
		if ( disori < CUTOFF_ENERGY ) disori = CUTOFF_ENERGY;

		energy = gamma_hagb * (disori/LAGB2HAGB) * (1.0 - log(disori/LAGB2HAGB));
	}
	//is there a clamp down at to disori --> 0.0?

	return ( mobility * energy );
}


double topoHdl::getAverageMobilityTimesEnergyPosition( uint32_t mydatapos, uint32_t whichmydataset )
{
	double value = 0.0;
	uint32_t ncand = this->mydatasize[whichmydataset];
	//MK::saves us a search as the position in whichmydataset is already known

	//how many faces and how long is the perimeter?
	uint32_t nf = mydata[whichmydataset][mydatapos].nfaces;
	uint32_t oid_target = mydata[whichmydataset][mydatapos].ori;
	adjointb2dP nbors = mydata[whichmydataset][mydatapos].neighbors2d;

	double len = 0.0;
	double disori = 0.0;
	double me = 0.0; //product mobility times energy
	double perimeter = 0.0;
	uint32_t nbor, nborpos, oid_nbor;
	for ( uint32_t f = 0; f < nf; f++ ) {
		len = nbors[f].length;

		//find this neighbor in local data
		nbor = nbors[f].globalid;

		for ( nborpos = 0; nborpos < ncand; nborpos++) { 
			if ( mydata[whichmydataset][nborpos].globalid == nbor ) break;
		}

		oid_nbor = mydata[whichmydataset][nborpos].ori;

		disori = this->misorientationCubicQxQ( oripool[oid_target].q0, oripool[oid_target].q1, oripool[oid_target].q2, oripool[oid_target].q3, oripool[oid_nbor].q0, oripool[oid_nbor].q1, oripool[oid_nbor].q2, oripool[oid_nbor].q3 );

		me = me + (getMobilityTimesEnergy( disori ) * len);

		perimeter = perimeter + len;

//cout << "\t\tpos;nf;f;nborid; len;peri;disori;me" << mydatapos << ";" << nf << ";" << f << ";" << nbor << ";" << len << ";" << perimeter << ";" << disori << ";" << me << endl;
	}

	if ( perimeter > DOUBLE_ACCURACY ) 
		value = me / perimeter;

	return value;
}



double topoHdl::getAverageMobilityTimesEnergyGlobalID( uint32_t globalid, uint32_t whichmydataset )
{
	double value = 0.0;
	uint32_t ncand = this->mydatasize[whichmydataset];
	uint32_t target;
	for ( target = 0; target < ncand; target++) {
		if ( mydata[whichmydataset][target].globalid == globalid ) 
			break;
	}

	//how many faces and how long is the perimeter?
	uint32_t nf = mydata[whichmydataset][target].nfaces;
	uint32_t oid_target = mydata[whichmydataset][target].ori;
	adjointb2dP nbors = mydata[whichmydataset][target].neighbors2d;

	double len = 0.0;
	double disori = 0.0;
	double me = 0.0; //product mobility times energy
	double perimeter = 0.0;
	uint32_t nbor, nborpos, oid_nbor;
	for ( uint32_t f = 0; f < nf; f++ ) {
		len = nbors[f].length;

		//find this neighbor in local data
		nbor = nbors[f].globalid;

		for ( nborpos = 0; nborpos < ncand; nborpos++) { 
			if ( mydata[whichmydataset][nborpos].globalid == nbor ) break;
		}

		oid_nbor = mydata[whichmydataset][nborpos].ori;

		disori = this->misorientationCubicQxQ( oripool[oid_target].q0, oripool[oid_target].q1, oripool[oid_target].q2, oripool[oid_target].q3, oripool[oid_nbor].q0, oripool[oid_nbor].q1, oripool[oid_nbor].q2, oripool[oid_nbor].q3 );

		me = me + (getMobilityTimesEnergy( disori ) * len);

		perimeter = perimeter + len;

//cout << "\t\tgid;nf;f;nborid; len;peri;disori;me" << globalid << ";" << nf << ";" << f << ";" << nbor << ";" << len << ";" << perimeter << ";" << disori << ";" << me << endl;
	}

	if ( perimeter > DOUBLE_ACCURACY ) 
		value = me / perimeter;

	return value;
}


void topoHdl::analyse_grain2d_sizeevolution( uint32_t first, uint32_t last)
{
	//mpiexec -n 200 topotrace_rzcluster 1000000 1400 3400  1.94s user 21.61s system 10% cpu 3:54.87 total

	//follow the size evolution for the successful grains
	//first find which node has last dataset
	if ( first < this->firstloaded || last > this->lastloaded || first > last || last < 0 || first < 0) { cout << "ERR::With the analysis limits." << endl; return; }
	//MK::so [first,last] is not exceeding beyond [firstloaded, lastloaded]

	uint32_t fid = (last / offsetloaded) - (firstloaded / offsetloaded);
	uint32_t last2assignmentid = UNKNOWN;
	for ( uint32_t a = 0; a < assignment.size(); a++ ) { //MK::simulation time increases for each dataset, so assignment is implicitly sorted ascending
		if ( assignment[a].globalid == fid ) {
			last2assignmentid = a;
			break;
		}
	}
	if ( last2assignmentid == UNKNOWN ) { cout << "It was not possible to identify the last timestep!" << endl; return; }

	int WhoKnowsLastTimestep = assignment[last2assignmentid].whichnode;
	uint32_t nGrainsWhichSurvived = assignment[last2assignmentid].ngrainsintotal; //e.g. if any
	if ( nGrainsWhichSurvived <= 0 ) { cout << "ERR::No grains survived or otherwise corrupted data in size evolution!" << endl; return; }


	//get a list of filter targets for all nodes to work on
	//WhoKnowsLastTimeStep compiles a list to bcast to all other workers
	uint32_t* GrainsWhichSurvived = NULL;
	GrainsWhichSurvived = new uint32_t[nGrainsWhichSurvived];
	if ( GrainsWhichSurvived == NULL ) { cout << "ERR::Allocating memory in analyse_grain2d" << endl; return; }

	if ( myRank == WhoKnowsLastTimestep ) {
		//output these grains the rowid and the globalid
		stringstream loglast_fn;
		loglast_fn << "Survivors.First." << first << ".Offset." << offsetloaded << ".Last." << last << ".csv";

		ofstream loglast;
		loglast.open( loglast_fn.str().c_str() );
		loglast << "RowID;GlobalID;NFaces;Area\n";


		//but which dataset
		uint32_t localid = UNKNOWN;
		for ( uint32_t i = 0; i < mydatasets.size(); i++ ) { //find which mydata container is it in my local data
			if ( mydatasets[i] == fid ) { localid = i; break; }
		}
		if ( localid == UNKNOWN ) { cout << myRank << " who supposedly kept the data for the last simstep is not able to find the data in his set!" << endl; return; }

		uint32_t largestgrain_gid = UNKNOWN;
		double largestgrain_area = 0.0;

		//copy all grains that survived to Bcast to the workers, so that WhoKnowsLastTimeStep can MPI_Bcast
		grain2dP thegrains = mydata[localid];
		for ( uint32_t gr = 0; gr < nGrainsWhichSurvived; gr++) {
			GrainsWhichSurvived[gr] = thegrains[gr].globalid;

			loglast << gr << ";" << thegrains[gr].globalid << ";" << thegrains[gr].nfaces << ";" << thegrains[gr].area << "\n";

			if ( mydata[localid][gr].area >= largestgrain_area ) {
				largestgrain_gid = thegrains[gr].globalid;
				largestgrain_area = thegrains[gr].area;
			}
		}

		loglast.flush();
		loglast.close();
cout << "THE LARGEST GRAIN IS " << largestgrain_gid << " with size " << setprecision(8) << largestgrain_area << endl;
	}

	MPI_Bcast( GrainsWhichSurvived, nGrainsWhichSurvived, MPI_UNSIGNED, WhoKnowsLastTimestep, MPI_COMM_WORLD );


	//open a binary file with collectively on all nodes
	MPI_File mpiiohdl_faces, mpiiohdl_vol, mpiiohdl_gbtype;
	MPI_Status mpiiosta_faces, mpiiosta_vol, mpiiosta_gbtype;

	//create C-consistent file name for MPI I/O
	ostringstream faces_fn, vol_fn, gbtype_fn;
	faces_fn << "topotrace.nfaces." << jobid << ".first." << first << ".offset." << offsetloaded << ".last." << last << ".bin";
	vol_fn << "topotrace.volume." << jobid << ".first." << first << ".offset." << offsetloaded << ".last." << last << ".bin";
	gbtype_fn << "topotrace.gbtype." << jobid << ".first." << first << ".offset." << offsetloaded << ".last." << last << ".bin";
	char* faces_fname = new char[faces_fn.str().size() + 1];		strcpy( faces_fname, faces_fn.str().c_str() );
	char* vol_fname = new char[vol_fn.str().size() + 1];			strcpy( vol_fname, vol_fn.str().c_str() );
	char* gbtype_fname = new char[gbtype_fn.str().size() + 1];		strcpy( gbtype_fname, gbtype_fn.str().c_str() );


	// open the file in create and write-only mode, ##only MASTER should finally write the file
	MPI_File_open(MPI_COMM_WORLD, faces_fname, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &mpiiohdl_faces);
	MPI_File_open(MPI_COMM_WORLD, vol_fname, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &mpiiohdl_vol);
	MPI_File_open(MPI_COMM_WORLD, gbtype_fname, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &mpiiohdl_gbtype);

	__int64 totalOffsetFaces = 0;
	__int64 totalOffsetVolume = 0;
	__int64 totalOffsetGBType = 0;

	//file format?
cout << "Size-evolution is output, the binary file has no header, format::GRAINS-IN-A-ROW-TIMESTEPS-IN-A-COLUMN" << endl;

	//create I/O buffer container to write timeslot after timeslot, only as many write ops as timeslots because one node has all grains
	double* volbuffer = NULL;
	volbuffer = new double[nGrainsWhichSurvived];
	if ( volbuffer == NULL ) { cout << "ERR::Allocating memory in analyse_grain2d" << endl; return; }
	uint32_t* nfacesbuffer = NULL;
	nfacesbuffer = new uint32_t[nGrainsWhichSurvived];
	if ( nfacesbuffer == NULL ) { cout << "ERR:Allocating memory in analyse_grain2d" << endl; return; }
	double* gbtypebuffer = NULL;
	gbtypebuffer = new double[nGrainsWhichSurvived];
	if ( gbtypebuffer == NULL ) { cout << "ERR::Allocating memory in analyze_grain2d" << endl; return; }

	//all nodes process in parallel, processes do not share data time-resolved data
	uint32_t myid;
	uint32_t howmanycolumns = 0;
	for ( uint32_t f = first; f <= last; f += offsetloaded ) {

		fid = (f/offsetloaded) - (firstloaded/offsetloaded);
		howmanycolumns++;

		if ( myRank == assignment[fid].whichnode ) { // I have this dataset, I analyze and write to file
//cout << "Size-evolution of f = " << f << " taken by node " << this->myRank << endl;
			myid = assignment[fid].localid;

			for ( uint32_t gr = 0; gr < nGrainsWhichSurvived; gr++) { //##MK::mind that the order in which WhoKnowsLastTimeStep delivers ids of surviving grains might not be preserved or strictly sorted
				nfacesbuffer[gr] = getNumberOfFaces ( GrainsWhichSurvived[gr], myid );
				volbuffer[gr] = getGrainVolume( GrainsWhichSurvived[gr], myid );
				gbtypebuffer[gr] = getAverageMobilityTimesEnergyGlobalID( GrainsWhichSurvived[gr], myid );
			}

			//calculate address in the file //globalid of the dataset * nGrainsWhichSurvived
			totalOffsetFaces = fid * nGrainsWhichSurvived * 4;
			totalOffsetVolume = fid * nGrainsWhichSurvived * 8;
			totalOffsetGBType = fid * nGrainsWhichSurvived * 8;

			MPI_File_write_at(mpiiohdl_faces, totalOffsetFaces, nfacesbuffer, nGrainsWhichSurvived, MPI_UNSIGNED, &mpiiosta_faces);
			MPI_File_write_at(mpiiohdl_vol, totalOffsetVolume, volbuffer, nGrainsWhichSurvived, MPI_DOUBLE, &mpiiosta_vol);
			MPI_File_write_at(mpiiohdl_gbtype, totalOffsetGBType, gbtypebuffer, nGrainsWhichSurvived, MPI_DOUBLE, &mpiiosta_gbtype);
		}
		//else, nothing because someone else will take care...

//cout << myRank << "\t" << nGrainsWhichSurvived << "\t" << mydatasets[mys] << endl;
	}

if ( myRank == MASTER ) { 
	cout << faces_fn.str().c_str() << " with  " << nGrainsWhichSurvived << " rows = grains survived, and " << howmanycolumns << " columns = timesteps in ascending order" << endl;
	cout << vol_fn.str().c_str() << " with " << nGrainsWhichSurvived << " rows = grains survived, and " << howmanycolumns << " columns = timesteps in ascending order" << endl;
	cout << gbtype_fn.str().c_str() << " with " << nGrainsWhichSurvived << " rows = grains survived, and " << howmanycolumns << " columns = timesteps in ascending order" << endl;
	cout << endl;
}

	delete [] nfacesbuffer;
	nfacesbuffer = NULL;
	delete [] volbuffer;
	volbuffer = NULL;
	delete [] gbtypebuffer;
	gbtypebuffer = NULL;
	delete [] GrainsWhichSurvived;
	GrainsWhichSurvived = NULL;
	delete [] faces_fname;
	faces_fname = NULL;
	delete [] vol_fname;
	vol_fname = NULL;
	delete [] gbtype_fname;
	gbtype_fname = NULL;

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_File_close(&mpiiohdl_faces);
	MPI_File_close(&mpiiohdl_vol);
	MPI_File_close(&mpiiohdl_gbtype);
}


uint32_t topoHdl::getOriIDOfNeighbor( uint32_t globalid, uint32_t thedataset )
{
	uint32_t ooid = UNKNOWN;
	uint32_t ncand = this->mydatasize[thedataset];
	grain2dP thegrains = mydata[thedataset];

	for ( uint32_t cand = 0; cand < ncand; cand++) {
		if ( thegrains[cand].globalid == globalid ) {
			ooid = thegrains[cand].ori;
			break;
		}
	}

	return ooid;
}


uint32_t topoHdl::checkCandidatesRiso( uint32_t datasetid, uint32_t whichmydataset )
{
	//is the candidate itself close to the cube orientation?
	uint32_t oid = mydata[whichmydataset][datasetid].ori;

	if ( this->misorientationCubic( 0.0, 0.0, 0.0, oripool[oid].bunge1, oripool[oid].bunge2, oripool[oid].bunge3 ) >= (20.0/180.0*_PI_) ) {
		return DO_NOT_INCLUDE;
	}

	//interesting grain, but also at least a certain neighborhood?
	uint32_t gid = DO_NOT_INCLUDE;
	uint32_t nnbors = this->mydata[whichmydataset][datasetid].nfaces;
	bool condition = false;

	for ( uint32_t nb = 0; nb < nnbors; nb++ ) {
		//##MK::DEBUG hard coded condition at least one S oriented neighbor close
		//orientation of the neighbor
		oid = getOriIDOfNeighbor( mydata[whichmydataset][datasetid].neighbors2d[nb].globalid, whichmydataset );

		if ( oid != UNKNOWN ) {
			if ( this->misorientationCubic( 59.0, 36.7, 63.4, oripool[oid].bunge1, oripool[oid].bunge2, oripool[oid].bunge3 ) <= (15.0/180.0*_PI_) ) {
				condition = true;
				//##MK::pull in gid assignment here...
				break;
			}
		}
	}

	if ( condition == true ) {
		gid = mydata[whichmydataset][datasetid].globalid;
	}

	return gid;
}


void topoHdl::analyse_grain2d_sizeevolution_nbormisori( uint32_t first, uint32_t last)
{
	//mpiexec -n 200 topotrace_rzcluster 1000000 1400 3400  1.94s user 21.61s system 10% cpu 3:54.87 total

	//MK::first identify grains located at the boundary of a bicrystal setup where the 
	//upper:S, 2deg
	//lower:Cube, 5deg
	//heuristic: find cube grains with at least one neighbor in S orientation, because at the moment no barycenter of the grain is provided

	//follow the size evolution for the successful grains first find which node has last dataset
	if ( first < this->firstloaded || last > this->lastloaded || first > last || last < 0 || first < 0) { cout << "ERR::With the analysis limits." << endl; return; }
	//MK::so [first,last] is not exceeding beyond [firstloaded, lastloaded]

	uint32_t fid = (last / offsetloaded) - (firstloaded / offsetloaded);
	uint32_t last2assignmentid = UNKNOWN;
	for ( uint32_t a = 0; a < assignment.size(); a++ ) { //MK::simulation time increases for each dataset, so assignment is implicitly sorted ascending
		if ( assignment[a].globalid == fid ) {
			last2assignmentid = a; break;
		}
	}
	if ( last2assignmentid == UNKNOWN ) { cout << "It was not possible to identify the last timestep!" << endl; return; }

	int WhoKnowsLastTimestep = assignment[last2assignmentid].whichnode;

	//node WhoKnowsLastTimeStep identifies all grains at the bicrystal boundary
	uint32_t nGrainsWhichSurvived = 0;
	uint32_t localid = UNKNOWN;
	uint32_t* GrainsWhichSurvived = NULL;
	uint32_t largestgrain_gid = UNKNOWN;
	double largestgrain_area = 0.0;
	if ( myRank == WhoKnowsLastTimestep ) { //which of my datasets is the from the specific timestep?
		for ( uint32_t i = 0; i < mydatasets.size(); i++ ) { //find which mydata container is it in my local data
			if ( mydatasets[i] == fid ) {
				localid = i; break;
			}
		}
		if ( localid == UNKNOWN ) { cout << myRank << " who supposedly kept the data for the last simstep is not able to find the data in his set!" << endl; return; } //##MK::change error handling here

		//output these grains
		stringstream loglast_fn;		loglast_fn << "GrainsAtBicrystalBoundary.First." << first << ".Offset." << offsetloaded << ".Last." << last << ".csv";
		ofstream loglast;				loglast.open( loglast_fn.str().c_str() );		loglast << "RowID;GlobalID;NFaces;Area\n";

		uint32_t ncand = mydatasize[localid];
		uint32_t* idbucket = NULL;
		idbucket = new uint32_t[ncand];
		if ( idbucket == NULL ) { cout << "ERR::Allocating memory in analyse_grain2d" << endl; return; }
		uint32_t targetid = DO_NOT_INCLUDE;
		for ( uint32_t cand = 0; cand < ncand; cand++ ) {
			targetid = checkCandidatesRiso( cand, localid ); //returns globalid

			if ( targetid != DO_NOT_INCLUDE ) {
				idbucket[nGrainsWhichSurvived] = targetid;
				nGrainsWhichSurvived++;

				loglast << cand << ";" << targetid << ";" << mydata[localid][cand].nfaces << ";" << mydata[localid][cand].area << "\n";
				if ( mydata[localid][cand].area >= largestgrain_area ) {
					largestgrain_gid = mydata[localid][cand].globalid;
					largestgrain_area = mydata[localid][cand].area;
				}
			}
		}
		cout << "THE LARGEST GRAIN IS " << largestgrain_gid << " with size " << setprecision(8) << largestgrain_area << endl;

		loglast.flush();
		loglast.close();

		//thin out idbucket to save MPI message length
		GrainsWhichSurvived = new uint32_t[nGrainsWhichSurvived];
		if ( GrainsWhichSurvived == NULL ) { cout << "ERR::Allocating memory in analyse_grain2d" << endl; return; }
		for ( uint32_t g = 0; g < nGrainsWhichSurvived; g++ ) {
			GrainsWhichSurvived[g] = idbucket[g];
		}

		delete [] idbucket;
	}

	MPI_Bcast( &nGrainsWhichSurvived, 1, MPI_UNSIGNED, WhoKnowsLastTimestep, MPI_COMM_WORLD );
	if ( nGrainsWhichSurvived <= 0 ) { cout << "ERR::No grains survived or otherwise corrupted data in size evolution!" << endl; return; }

	if ( myRank != WhoKnowsLastTimestep ) {
		GrainsWhichSurvived = new uint32_t[nGrainsWhichSurvived];
		if ( GrainsWhichSurvived == NULL ) { cout << "ERR::Allocating memory in analyse_grain2d" << endl; return; }
	}
	//WhoKnowsLastTimeStep compiles a list to bcast to all other workers
	MPI_Bcast( GrainsWhichSurvived, nGrainsWhichSurvived, MPI_UNSIGNED, WhoKnowsLastTimestep, MPI_COMM_WORLD );


	//open a binary file with collectively on all nodes
	MPI_File mpiiohdl_faces, mpiiohdl_vol, mpiiohdl_gbtype;
	MPI_Status mpiiosta_faces, mpiiosta_vol, mpiiosta_gbtype;

	//create C-consistent file name for MPI I/O
	ostringstream faces_fn, vol_fn, gbtype_fn;
	faces_fn << "topotrace.nfaces." << jobid << ".first." << first << ".offset." << offsetloaded << ".last." << last << ".bin";
	vol_fn << "topotrace.volume." << jobid << ".first." << first << ".offset." << offsetloaded << ".last." << last << ".bin";
	gbtype_fn << "topotrace.gbtype." << jobid << ".first." << first << ".offset." << offsetloaded << ".last." << last << ".bin";
	char* faces_fname = new char[faces_fn.str().size() + 1];		strcpy( faces_fname, faces_fn.str().c_str() );
	char* vol_fname = new char[vol_fn.str().size() + 1];			strcpy( vol_fname, vol_fn.str().c_str() );
	char* gbtype_fname = new char[gbtype_fn.str().size() + 1];		strcpy( gbtype_fname, gbtype_fn.str().c_str() );


	// open the file in create and write-only mode, ##only MASTER should finally write the file
	MPI_File_open(MPI_COMM_WORLD, faces_fname, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &mpiiohdl_faces);
	MPI_File_open(MPI_COMM_WORLD, vol_fname, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &mpiiohdl_vol);
	MPI_File_open(MPI_COMM_WORLD, gbtype_fname, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &mpiiohdl_gbtype);

	__int64 totalOffsetFaces = 0;
	__int64 totalOffsetVolume = 0;
	__int64 totalOffsetGBType = 0;

	//file format?
cout << "Size-evolution is output, the binary file has no header, format::GRAINS-IN-A-ROW-TIMESTEPS-IN-A-COLUMN" << endl;

	//create I/O buffer container to write timeslot after timeslot, only as many write ops as timeslots because one node has all grains
	double* volbuffer = NULL;
	volbuffer = new double[nGrainsWhichSurvived];
	if ( volbuffer == NULL ) { cout << "ERR::Allocating memory in analyse_grain2d" << endl; return; }
	uint32_t* nfacesbuffer = NULL;
	nfacesbuffer = new uint32_t[nGrainsWhichSurvived];
	if ( nfacesbuffer == NULL ) { cout << "ERR:Allocating memory in analyse_grain2d" << endl; return; }
	double* gbtypebuffer = NULL;
	gbtypebuffer = new double[nGrainsWhichSurvived];
	if ( gbtypebuffer == NULL ) { cout << "ERR::Allocating memory in analyze_grain2d" << endl; return; }

	//all nodes process in parallel, processes do not share data time-resolved data
	uint32_t myid;
	uint32_t howmanycolumns = 0;
	for ( uint32_t f = first; f <= last; f += offsetloaded ) {

		fid = (f/offsetloaded) - (firstloaded/offsetloaded);
		howmanycolumns++;

		if ( myRank == assignment[fid].whichnode ) { // I have this dataset, I analyze and write to file
//cout << "Size-evolution of f = " << f << " taken by node " << this->myRank << endl;
			myid = assignment[fid].localid;

			for ( uint32_t gr = 0; gr < nGrainsWhichSurvived; gr++) { //##MK::mind that the order in which WhoKnowsLastTimeStep delivers ids of surviving grains might not be preserved or strictly sorted
				nfacesbuffer[gr] = getNumberOfFaces ( GrainsWhichSurvived[gr], myid );
				volbuffer[gr] = getGrainVolume( GrainsWhichSurvived[gr], myid );
				gbtypebuffer[gr] = getAverageMobilityTimesEnergyGlobalID( GrainsWhichSurvived[gr], myid );
			}

			//calculate address in the file //globalid of the dataset * nGrainsWhichSurvived
			totalOffsetFaces = fid * nGrainsWhichSurvived * 4;
			totalOffsetVolume = fid * nGrainsWhichSurvived * 8;
			totalOffsetGBType = fid * nGrainsWhichSurvived * 8;

			MPI_File_write_at(mpiiohdl_faces, totalOffsetFaces, nfacesbuffer, nGrainsWhichSurvived, MPI_UNSIGNED, &mpiiosta_faces);
			MPI_File_write_at(mpiiohdl_vol, totalOffsetVolume, volbuffer, nGrainsWhichSurvived, MPI_DOUBLE, &mpiiosta_vol);
			MPI_File_write_at(mpiiohdl_gbtype, totalOffsetGBType, gbtypebuffer, nGrainsWhichSurvived, MPI_DOUBLE, &mpiiosta_gbtype);
		}
		//else, nothing because someone else will take care...

//cout << myRank << "\t" << nGrainsWhichSurvived << "\t" << mydatasets[mys] << endl;
	}

if ( myRank == MASTER ) { 
	cout << faces_fn.str().c_str() << " with  " << nGrainsWhichSurvived << " rows = grains survived, and " << howmanycolumns << " columns = timesteps in ascending order" << endl;
	cout << vol_fn.str().c_str() << " with " << nGrainsWhichSurvived << " rows = grains survived, and " << howmanycolumns << " columns = timesteps in ascending order" << endl;
	cout << gbtype_fn.str().c_str() << " with " << nGrainsWhichSurvived << " rows = grains survived, and " << howmanycolumns << " columns = timesteps in ascending order" << endl;
	cout << endl;
}

	delete [] nfacesbuffer;
	nfacesbuffer = NULL;
	delete [] volbuffer;
	volbuffer = NULL;
	delete [] gbtypebuffer;
	gbtypebuffer = NULL;
	delete [] GrainsWhichSurvived;
	GrainsWhichSurvived = NULL;
	delete [] faces_fname;
	faces_fname = NULL;
	delete [] vol_fname;
	vol_fname = NULL;
	delete [] gbtype_fname;
	gbtype_fname = NULL;

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_File_close(&mpiiohdl_faces);
	MPI_File_close(&mpiiohdl_vol);
	MPI_File_close(&mpiiohdl_gbtype);
}


uint32_t topoHdl::lastVitalitySign( uint32_t tgr )
{
	if ( mydata.size() != mydatasets.size() ) { cout << "Local mydata(sets) size inconsistency!" << endl; return NEVER_OBSERVED; }

	uint32_t fidmax = NEVER_OBSERVED;
	for ( uint32_t mys = 0; mys < mydata.size(); mys++) {
		bool found = false;
		uint32_t ngrains = mydatasize[mys];
		grain2dP thegrains = mydata[mys];

		for ( uint32_t g = 0; g < ngrains; g++) {
			if ( thegrains[g].globalid == tgr ) {
				found = true;
				break;
			}
		}

		if ( found == true ) {
			//which time step is this, get the fid = f/offsetloaded
			uint32_t fid = mydatasets[mys];

			if ( fid > fidmax ) {
				fidmax = fid;
			}
		}
		//all data have to be checked because nodes have not necessarily linear data-time-packages in their local memory
	}

	//assured accessible reference and the last time, the process identified the tgr
	return fidmax;
}


uint32_t* topoHdl::withMPIfindAllDissimilarNBors( uint32_t tgr, uint32_t lastfid, uint32_t& howmanyever )
{
	//##MK::several approaches possible, slightly tricky...
	//    ------        --       
	// ---      -----        ----
	//____________________________t-axis, obviously number of neighbors and time interval along which they are neighbors changes
	//naive: forall fid forall ranks MPISend/MPIRecv --> too many many too small messages, but short but blocking
	//naive: forall fid forall g --> too many comparisons in forced sequentialization, 

	//first all nodes collect in parallel all disjoint grain ids that somewhen became or were located next to tgr
	vector<uint32_t> mypart_nborgrids;
	vector<MPI_IO_GrainInfo2D> mypart_tgrevo;

	for ( uint32_t mys = 0; mys < mydata.size(); mys++) {
		//assuming that tgr has not vanished before lastfid + 1
		if ( mydatasets[mys] > lastfid ) continue;

		//find the position of tgr in mydata because globalids, e.g. grain identifiers are not necessarily running from 0 to mydatasize[mys]-1 and number of grains decreasing...
		uint32_t mygid = UNKNOWN;
		uint32_t ng = mydatasize[mys];
		grain2dP thegrains = mydata[mys];

		for ( uint32_t g = 0; g < ng; g++ ) {
			if ( thegrains[g].globalid == tgr ) {
				mygid = g;
				break;
			}
		}
		if ( mygid == UNKNOWN ) { cout << "ERR::Unable to locate the targetgrain in myRank/mydatasetid while searching dissimilar grain ids" << myRank << ";" << mys << endl; return NULL; }

		//I collect the information I have about the target grain
		MPI_IO_GrainInfo2D tg;
		tg.fid = mydatasets[mys];
		tg.nfaces = mydata[mys][mygid].nfaces;
		tg.area = mydata[mys][mygid].area;
		tg.mobilityenergy = getAverageMobilityTimesEnergyPosition( mygid, mys );
		mypart_tgrevo.push_back( tg );

//cout << myRank << " scanning through my dataset " << mys << " found tgr " << tgr << " as entry " << mygid << " with nfaces;area;" << this->mydata[mys][mygid].nfaces << ";" << this->mydata[mys][mygid].area << endl;

//collect the neighbors in this timestep
		for ( uint32_t ngr = 0; ngr < mydata[mys][mygid].nfaces; ngr++) {
			uint32_t nid = mydata[mys][mygid].neighbors2d[ngr].globalid;
			bool unknown = true;
			for ( uint32_t founds = 0; founds < mypart_nborgrids.size(); founds++ ) {
				if ( nid == mypart_nborgrids[founds] ) {
					unknown = false;
					break;
				}
			}
			if ( unknown == true ) {
				mypart_nborgrids.push_back( nid );
//cout << "TimeStep;" << mys << "\tnew neighbor\t" << nid << endl;
			}
		} //for all neighbors of mygid := tgr
	} //in all my datasets parallel in all nodes

	//wrap present for MPI to send all dissimilar neighbors I found
	uint32_t mylen = mypart_nborgrids.size();
 	uint32_t* myarr = NULL;
	myarr = new uint32_t[mylen];
	if ( myarr == NULL ) { cout << "ERR::Allocation error in withMPIFindAllDissimilarNBors!" << endl; return NULL; }
	for ( uint32_t i = 0; i < mylen; i++) {
		myarr[i] = mypart_nborgrids[i]; 
	}

	//now all nodes have performed this analysis, now one by one collect dissimilar once from the workers in master then broadcast result in all

MPI_Barrier( MPI_COMM_WORLD );

	//broadcast one by one the information of the target grain to the nodes
	for ( int rr = MASTER; rr < nRanks; rr++) {
		uint32_t ginfosize = mypart_tgrevo.size();

		MPI_Bcast( &ginfosize, 1, MPI_UNSIGNED, rr, MPI_COMM_WORLD);

		MPI_IO_GrainInfo2D* ginfobucket = NULL;
		ginfobucket = new MPI_IO_GrainInfo2D[ginfosize];
		if ( ginfobucket == NULL ) { cout << "Allocation error in MPI_Bcast ginfobucket in withMPI" << endl; return NULL; }
		//the rr node wraps its MPI present

		if ( myRank == rr ) {
			for ( uint32_t gg = 0; gg < mypart_tgrevo.size(); gg++ ) {
				ginfobucket[gg].fid = mypart_tgrevo[gg].fid;
				ginfobucket[gg].nfaces = mypart_tgrevo[gg].nfaces;
				ginfobucket[gg].area = mypart_tgrevo[gg].area;
				ginfobucket[gg].mobilityenergy = mypart_tgrevo[gg].mobilityenergy;
			}
		}

		MPI_Bcast( ginfobucket, ginfosize, MPI_IO_GrainInfo2D_Type, rr, MPI_COMM_WORLD);

		//incorporate information into own bucket
		for ( uint32_t t = 0; t < ginfosize; t++ ) {
			MPI_IO_GrainInfo2D tt;
			tt.fid = ginfobucket[t].fid;
			tt.nfaces = ginfobucket[t].nfaces;
			tt.area = ginfobucket[t].area;
			tt.mobilityenergy = ginfobucket[t].mobilityenergy;
			targetgrain_evolution.push_back(tt);
		}

		delete [] ginfobucket;
	} //for all nodes

	//sort evolution of the grains utilizing that fid are disjoint
	std::sort ( targetgrain_evolution.begin(), targetgrain_evolution.end(), SortTimeStepAscending );

//##DDEBUGMK::Dummy
//cout << endl; for (uint32_t ts = 0; ts < targetgrain_evolution.size(); ts++) { cout << "myRank;fid;nf;A;mgamma\t\t" << this->myRank << ";" << targetgrain_evolution[ts].fid << ";" << targetgrain_evolution[ts].nfaces << ";" << targetgrain_evolution[ts].area << ";" << targetgrain_evolution[ts].mobilityenergy << endl; } cout << endl;
//##DDEBUGMK::Dummy

	//aggregate dissimilar neighbors into list of disjoint ids that all nodes know
	vector<uint32_t> resid;

	//MASTER copies first its own
	if ( myRank == MASTER ) {
		for ( uint32_t j = 0; j < mylen; j++ )
			resid.push_back ( myarr[j] ); 
	}

	for (int r = (MASTER + 1); r < this->nRanks; r++) {
cout << "MASTER/WORKER data consolidation in withMPIDissimilar ones." << endl;
		//Master queries mylen from r
		uint32_t rsize = 0;
		uint32_t* rbucket = NULL;

		//send/recv array size
		if ( myRank == r )
			MPI_Send( &mylen, 1, MPI_UNSIGNED, MASTER, r, MPI_COMM_WORLD );
		if ( myRank == MASTER ) {
			MPI_Recv( &rsize, 1, MPI_UNSIGNED, r, r, MPI_COMM_WORLD, MPI_STATUS_IGNORE );

			rbucket = new uint32_t[rsize];
			if ( rbucket == NULL ) { cout << "ERR::Allocation erro in rbucket in withMPI..." << endl; return NULL; }
		}

		//send/receive data
		if ( myRank == r ) {
			MPI_Send( myarr, mylen, MPI_UNSIGNED, MASTER, r, MPI_COMM_WORLD ); //MPI standard requires tags > SHORT_MAX_RANGE
		}
		if ( myRank == MASTER ) {
			MPI_Recv( rbucket, rsize, MPI_UNSIGNED, r, r, MPI_COMM_WORLD, MPI_STATUS_IGNORE );

			//analyze data for disjoint ids
			bool unknown;
			uint32_t cand;
			for ( uint32_t rd = 0; rd < rsize; rd++ ) {
				unknown = true;
				cand = rbucket[rd];
				for ( uint32_t c = 0; c < resid.size(); c++ ) {
					if ( cand == resid[c] ) {
						unknown = false;
						break;
					}
				}
				if ( unknown == true ) {
					resid.push_back ( cand );
				}
			} //collect all disjoint ...
		}

		delete [] rbucket;

		//##MK::maybe obsolete because MPI_Send/Recv is blocking and MPI_Tags dijoint but better do so
		MPI_Barrier ( MPI_COMM_WORLD );

cout << "WithMPI collecting local container of neighbors from tgr = " << tgr << " from node " << r << endl;
	} //...from the next node

	//summarize results on the master
	uint32_t reslen = 0;
	uint32_t* master_res = NULL;

	if ( myRank == MASTER ) {
		reslen = resid.size();
		master_res = new uint32_t[reslen];
		if ( master_res == NULL ) { cout << myRank << "throws ERR::Allocation error for res IN THE MASTER in withMPIFind..." << endl; return NULL; }
		for ( uint32_t c = 0; c < reslen; c++) {
			master_res[c] = resid[c];
		}
	}
	//commit to all nodes the final size reslen
	MPI_Bcast( &reslen, 1, MPI_UNSIGNED, MASTER, MPI_COMM_WORLD );
	if ( reslen <= 0 ) { cout << "Invalid reslen message of size " << reslen << endl; }

	//commit to all nodes the disjoint neighbors
	if ( myRank != MASTER ) {
		master_res = new uint32_t[reslen];
		if ( master_res == NULL ) { cout << myRank << "throws ERR::Allocation error for res IN THE WORKER in withMPIFind..." << endl; return NULL; }
	}
	MPI_Bcast( master_res, reslen, MPI_UNSIGNED, MASTER, MPI_COMM_WORLD );

	//cleanup
	delete [] myarr;


//##DEBUG
//cout << "Found all neighbors ever;myRank;" << this->myRank << "--" << endl; for ( uint32_t k = 0; k < reslen; k++) { cout << master_res[k] << ";" << endl; } cout << endl;
//##DEBUG


	howmanyever = reslen;
	return master_res;
}


void topoHdl::analyse_grain2d_onegrain_vs_nearestnbors( uint32_t targetgr, uint32_t first, uint32_t last)
{
	//follow the evolution of a single grain in its environment of neighbors with respect to size and topology
	if (first < this->firstloaded || last > this->lastloaded || first > last || last < 0) {
		cout << "ERR:Invalid file descriptors for input data!" << endl; return;
	}

	uint32_t fidstart = (first / offsetloaded) - (firstloaded / offsetloaded);
	uint32_t fidend = (last / offsetloaded) - (firstloaded / offsetloaded);

	//up to which timestep did the grain survive, assuming that the grain is existent in the first place
	uint32_t LastTimeStepGrainExistedLocal = lastVitalitySign( targetgr );
	uint32_t LastTimeStepGrainExistedGlobal = NEVER_OBSERVED;
	MPI_Allreduce( &LastTimeStepGrainExistedLocal, &LastTimeStepGrainExistedGlobal, 1, MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD );
	if ( LastTimeStepGrainExistedGlobal < fidstart ) { cout << "ERR:Grain was not observed during [first, last]!" << endl; return; }

cout << this->myRank << " targetgrain " << targetgr << " was last seen in timestep " << LastTimeStepGrainExistedGlobal << endl;

	//in order to construct a matrix of all neighbors which became somewhen a neighbor to the targetgr we need to solve the problem that
	//the nodes do not know when grains became neighbor to the target grain, at the same time in order for 
	//utilization of MPI/IO all nodes need to know how many disjoint neighbors at any point [first, last] were there in total, such that fixed-size data stripes can be written in parallel

	MPI_Barrier ( MPI_COMM_WORLD );

	uint32_t HowManyNBorsEver = NO_NEIGHBORS_EVER;
	uint32_t* TheNBorsEver = NULL;

	TheNBorsEver = withMPIfindAllDissimilarNBors( targetgr, LastTimeStepGrainExistedGlobal, HowManyNBorsEver );

//##DEBUG
//cout << endl;for (int n = 0; n < HowManyNBorsEver; n++) { cout << TheNBorsEver[n] << ";"; }cout << endl;
//##DEBUG

	if ( HowManyNBorsEver == NO_NEIGHBORS_EVER || TheNBorsEver == NULL ) { cout << "ERR::Mistake in withMPIfindAllDissimilarNeighbors." << endl; return; }

	//as now the total number of neighbors ever in the dataset is known we can use MPI I/O in parallel
	//##MK::this is not the best solution as now each element from TheNBorsEver is queried even if it might not be present in the node
	MPI_File mpiiohdl_faces, mpiiohdl_vol, mpiiohdl_gb;
	MPI_Status mpiiosta_faces, mpiiosta_vol, mpiiosta_gb;

	ostringstream faces_fn, vol_fn, gb_fn;
	faces_fn << "followingtarget." << targetgr << ".nfaces." << jobid << ".first." << first << ".offset." << offsetloaded << ".last." << last << ".bin";
	vol_fn << "followingtarget." << targetgr << ".volume." << jobid << ".first." << first << ".offset." << offsetloaded << ".last." << last << ".bin";
	gb_fn << "followingtarget." << targetgr << ".gbtype." << jobid << ".first." << first << ".offset." << offsetloaded << ".last." << last << ".bin";
	char* faces_fname = new char[faces_fn.str().size() + 1];		strcpy( faces_fname, faces_fn.str().c_str() );
	char* vol_fname = new char[vol_fn.str().size() + 1];			strcpy( vol_fname, vol_fn.str().c_str() );
	char* gb_fname = new char[gb_fn.str().size() + 1];			strcpy( gb_fname, gb_fn.str().c_str() );


	// open the file in create and write-only mode
	MPI_File_open(MPI_COMM_WORLD, faces_fname, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &mpiiohdl_faces);
	MPI_File_open(MPI_COMM_WORLD, vol_fname, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &mpiiohdl_vol);
	MPI_File_open(MPI_COMM_WORLD, gb_fname, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &mpiiohdl_gb);
	__int64 toffset_faces = 0;
	__int64 toffset_vol = 0;
	__int64 toffset_gb = 0;

cout << "Following-Single-Grain-GRAINS-IN-ROWS-TIMESTEPS-IN-COLS" << endl;

	uint32_t* facebuffer = NULL;
	facebuffer = new uint32_t[HowManyNBorsEver];
	if ( facebuffer == NULL ) { cout << "Allocation error for the facebuffer in vs_nearestneighbors" << endl; }
	double* volbuffer = NULL;
	volbuffer = new double[HowManyNBorsEver];
	if ( volbuffer == NULL ) { cout << "Allocation error for the volbuffer in vs_nearestneighbors" << endl; }
	double* gbbuffer = NULL;
	gbbuffer = new double[HowManyNBorsEver];
	if ( gbbuffer == NULL ) { cout << "Allocation error for the gbbuffer in vs_nearestneighbors" << endl; }


	uint32_t fid;
	for ( uint32_t f = first; f <= last; f += offsetloaded ) {
		fid = ( f / offsetloaded) - (firstloaded / offsetloaded); //in relation to all timesteps

		//grain never observed in this timestep?
		if ( fid > LastTimeStepGrainExistedGlobal ) continue;

		//collect local data in buffers
		int targetnode = assignment[fid].whichnode;
		uint32_t mys = assignment[fid].localid;

		//set dummy state for all grains not of interest in the matrix
		for ( uint32_t nbever = 0; nbever < HowManyNBorsEver; nbever++ ) {
			facebuffer[nbever] = 0;
			volbuffer[nbever] = 0.0;
			gbbuffer[nbever] = 0.0;
		}

		//collect data
//cout << "f;myRank;targetnode;fid;ass[fid].which\t" << f << "\t" << myRank << "\t" << targetnode << "\t" << fid << "\t" << assignment[fid].whichnode << endl;
		if ( this->myRank == targetnode ) { //MK::globalids are disjoint
//cout << "In IF loop, mydatasize, localid" << "\t" << this->mydatasize[assignment[fid].localid] << "\t" << assignment[fid].localid << endl;
			//##MK::if volume from all neighbors desired...

			//##MK::parallel optimization potential here, for instance cache misses via mydata[mys]... can be avoided by loading a datastructure with the current neighbors
			//MK::here the choice was made to put lean memory over performance

			//search position of the targetgrain
			uint32_t positionOfTheTarget = UNKNOWN;
			uint32_t ngt = mydatasize[mys];
			grain2dP thegrains = mydata[mys];
			for ( uint32_t gt = 0; gt < ngt; gt++ ) {
				if ( thegrains[gt].globalid == targetgr ) {
					positionOfTheTarget = gt;
					break;
				}
			}
			if ( positionOfTheTarget == UNKNOWN ) { cout << "ERR::Unable to identify position of the target grain in f = " << fid << endl; return; }


			//analyze correlation properties of the neighbors to the targetgrain right now in fid
			//by finding from all ids that become at least once ever a neighbor to the target all those 
			//which are a neighbor right now in fid
			uint32_t gid;
			for ( uint32_t g = 0; g < ngt; g++) {
				gid = thegrains[g].globalid;

				//is this entry useful information for following a neighboring grain to targetgr?
				for ( uint32_t nbever = 0; nbever < HowManyNBorsEver; nbever++) {

					bool is_g_now_anbor_totarget = false;
					for ( uint32_t nn = 0; nn < thegrains[positionOfTheTarget].nfaces; nn++) { //is gid any of the current neighbors?
						if ( gid == mydata[mys][positionOfTheTarget].neighbors2d[nn].globalid ) {
							is_g_now_anbor_totarget = true;
							break; //stops testing neighborship, gid is a neighbor
						}
					}

					if ( gid == TheNBorsEver[nbever] && is_g_now_anbor_totarget == true ) { //then I am interest in the state of this neighbor with its faces and area
						facebuffer[nbever] = mydata[mys][g].nfaces;
						volbuffer[nbever] = mydata[mys][g].area;
						gbbuffer[nbever] = getAverageMobilityTimesEnergyGlobalID( gid, mys );

//cout << "nbever;gid;TheNBorsEver[nbever];nface;area;" << nbever << ";" << ";" << TheNBorsEver[nbever]; // << ";" << nface << ";" << setprecision(8) << area << endl;

						break; //as grains are disjoint g cannot be one of the following nbever ids!
					}
				} //for all nbors
//cout << endl;
			} //check all grains

			//MK::that first == firstloaded is not guaranteed!
			toffset_faces = ( (f/offsetloaded) - (first/offsetloaded) ) * HowManyNBorsEver * 4;
			toffset_vol =	( (f/offsetloaded) - (first/offsetloaded) ) * HowManyNBorsEver * 8;
			toffset_gb =	( (f/offsetloaded) - (first/offsetloaded) ) * HowManyNBorsEver * 8;


			MPI_File_write_at(mpiiohdl_faces, toffset_faces, facebuffer, HowManyNBorsEver, MPI_UNSIGNED, &mpiiosta_faces);
			MPI_File_write_at(mpiiohdl_vol, toffset_vol, volbuffer, HowManyNBorsEver, MPI_DOUBLE, &mpiiosta_vol);
			MPI_File_write_at(mpiiohdl_gb, toffset_gb, gbbuffer, HowManyNBorsEver, MPI_DOUBLE, &mpiiosta_gb);
		}
//cout << myRank << " writing f = " << f << " which is fid=" << fid << " Howmanynborsever;offsetFace;offsetVol " << HowManyNBorsEver << "\t" << (double) toffset_faces << "\t" << (double) toffset_vol << endl;
	} //all timesteps

	delete [] facebuffer;
	delete [] volbuffer;
	delete [] gbbuffer;
	delete [] TheNBorsEver;

	MPI_Barrier( MPI_COMM_WORLD );
	MPI_File_close(&mpiiohdl_faces);
	MPI_File_close(&mpiiohdl_vol);
	MPI_File_close(&mpiiohdl_gb);

	cout << faces_fn.str().c_str() << " with  " << HowManyNBorsEver << " rows = neighbors and " << (last/offsetloaded - first/offsetloaded + 1) << " columns = timesteps in ascending order" << endl;
	cout << vol_fn.str().c_str() << " with " << HowManyNBorsEver << " rows = neighbors and " << (last/offsetloaded - first/offsetloaded + 1) << " columns = timesteps in ascending order" << endl;
	cout << gb_fn.str().c_str() << " with " << HowManyNBorsEver << " rows = neighbors and " << (last/offsetloaded - first/offsetloaded + 1) << " columns = timesteps in ascending order" << endl;
	cout << endl;

	//output single grain evolution
	if ( myRank == MASTER ) {
		stringstream logsgr_fn;
		logsgr_fn << "followinggrain." << targetgr << ".Evolution." << this->jobid << ".csv";

		ofstream logsgr;
		logsgr.open( logsgr_fn.str().c_str() );
		logsgr << "TimeStep-" << targetgr << ";Fid;NFaces;Area;MobilityEnergyOverPerimeter\n";
		for ( uint32_t s = 0; s < targetgrain_evolution.size(); s++ ) {
			logsgr << s << ";" << targetgrain_evolution[s].fid << ";" << targetgrain_evolution[s].nfaces << ";" << targetgrain_evolution[s].area << ";" << targetgrain_evolution[s].mobilityenergy << "\n";
		}
		logsgr.flush();
		logsgr.close();
	}

	delete [] faces_fname;
	delete [] vol_fname;
	delete [] gb_fname;

	this->targetgrain_evolution.clear(); //MK::but this does not reduce the capacity!
}


uint32_t* topoHdl::cutTopoTail( uint32_t* arr, uint32_t CurrentSize, uint32_t WhichToCutFirst )
{
	//cuts an array of CurrentSize and copies the individual elements < WhichToCutFirst
	if ( CurrentSize < WhichToCutFirst ) { cout << "Array in cutTail can't cut arr farther apart than its size!" << endl; return NULL; }

	uint32_t* ennya_arr = NULL;
	ennya_arr = new uint32_t[WhichToCutFirst];
	if ( ennya_arr == NULL ) { cout << "Allocation error in cutTail." << endl; return NULL; }
	myMemGuard = myMemGuard + (WhichToCutFirst * sizeof(uint32_t));

	for (uint32_t j = 0; j < WhichToCutFirst; j++ ) {
		ennya_arr[j] = arr[j];
	}

	delete [] arr;
	myMemGuard = myMemGuard - (CurrentSize * sizeof(int));

	return ennya_arr;
}


uint32_t* topoHdl::getTopoMemory( uint32_t* arr, uint32_t oldsize, uint32_t wheretoplace, uint32_t& resized )
{
	//only called when remalloc necessary
	uint32_t* newarr = NULL;
	uint32_t newsize = 2 * wheretoplace;
	if ( newsize >= MAXIMUM_NUMBER_OFGRAINS ) { 
		cout << "Allocation error in getTopoMemory!" << endl;
		resized = 0;
		return newarr;
	}

	newarr = new uint32_t[newsize];
	if ( newarr == NULL ) { 
		cout << "ERROR::Allocation error in getTopoMemory!" << endl;
		resized = 0;
		return NULL; 
	}
	myMemGuard = myMemGuard + (newsize * sizeof(uint32_t));

	//copy old elements from arr into newarr
	for ( uint32_t j = 0; j < oldsize; j++ ) {
		newarr[j] = arr[j];
	}

	delete [] arr;
	myMemGuard = myMemGuard - (oldsize * sizeof(uint32_t));

	resized  = newsize;
	return newarr;
}

#define FIRST_FIRSTNN_NBOR		0
#define FIRST_SECONDNN_NBOR		1
#define FIRST_THIRDNN_NBOR		2
#define FIRST_FOURTHNN_NBOR		3


void topoHdl::find_neighbors( uint32_t whichset, uint32_t whichgrain, nbortopoP target )
{
	//currently implementation is hardcoded...
	uint32_t cutoff = 4; //DEFAULT_KSHELL;

	//get the grainid around which to construct the shell
	uint32_t cgrid = mydata[whichset][whichgrain].globalid;
	uint32_t cnfc = mydata[whichset][whichgrain].nfaces;
	uint32_t iisize = cnfc + pow( 6.0, (cutoff - 1.0) );

	//no safety implemented at the moment for the cutoff, get drop off memory for neighbors to start with
	uint32_t* nn = NULL;
	nn = new uint32_t[cutoff];
	if (nn == NULL) { cout << "Allocation error nn inside find_neighbors" << endl; return; }
	myMemGuard = myMemGuard + (cutoff*sizeof(uint32_t));
	for (uint32_t j = 0; j < cutoff; j++) nn[j] = 0;

//first all direct neighbors, nn contains how many neighbors in that k-shell
	nn[FIRST_FIRSTNN_NBOR] = 0; 

//construct nearest shell by analyzing the adjacent neighbors
	uint32_t* ii = NULL;
	ii = new uint32_t[iisize];
	myMemGuard = myMemGuard + (iisize * sizeof(uint32_t));
	if ( ii == NULL ) { cout << "Allocation error ii inside find_neighbors" << endl; return; }

	//counting disjoint neighbors that become identified during the analysis
	uint32_t disjoints = 0;

//begin with first layer, these ids are already known...
	for ( uint32_t fc = 0; fc < mydata[whichset][whichgrain].nfaces; fc++ ) {
		ii[fc] = mydata[whichset][whichgrain].neighbors2d[fc].globalid; //iisize >= nfaces guaranteed
		disjoints++;
	}
	nn[FIRST_SECONDNN_NBOR] = disjoints; //then the farther apart located neighbors
cout << "\t\t-->" << disjoints;

//proceed with second layer
	for (uint32_t j = nn[FIRST_FIRSTNN_NBOR]; j < nn[FIRST_SECONDNN_NBOR]; j++ ) {
		uint32_t jnfc = mydata[whichset][ii[j]].nfaces; //check the neighbors of the neighbors of cgrid

		for (uint32_t jnbors = 0; jnbors < jnfc; jnbors++) {
			bool unknown = true; //assume each neighbor's neighbor is a higher order neighbor not accounted yet and not the target grain itself...
			uint32_t jid = mydata[whichset][ii[j]].neighbors2d[jnbors].globalid;

			for ( uint32_t known = 0; known < disjoints; known++ ) {
				if ( jid == ii[known] || jid == cgrid ) {
					unknown = false;
					//break;
				}
			}
			if ( unknown == true ) { //unknown is very likely
//cout << "jnbors,unknown=" << jnbors << endl;
				uint32_t* iinew = ii;
				uint32_t currsize = iisize;
				if ( iisize <= disjoints ) { iinew = getTopoMemory( ii, iisize, disjoints, currsize ); } //i becomes deleted in array ii to be able to place at least disjoint elements
				ii = iinew;
				iisize = currsize;
				ii[disjoints] = jid;
				disjoints++;
	}	}	}
	nn[FIRST_THIRDNN_NBOR] = disjoints;
cout << "\t" << disjoints;

//proceed with third layer
	for (uint32_t j = nn[FIRST_SECONDNN_NBOR]; j < nn[FIRST_THIRDNN_NBOR]; j++ ) {
		uint32_t jnfc = mydata[whichset][ii[j]].nfaces; //check neighbors neighbors

		for (uint32_t jnbors = 0; jnbors < jnfc; jnbors++) {
			bool unknown = true; //assume each neighbor's neighbor is a higher order neighbor not accounted yet and not the target grain itself...
			uint32_t jid = mydata[whichset][ii[j]].neighbors2d[jnbors].globalid;

			for ( uint32_t known = 0; known < disjoints; known++ ) {
				if ( jid == ii[known] || jid == cgrid ) {
					unknown = false;
					//break;
				}
			}
			if ( unknown == true ) { //unknown is very likely
				uint32_t* iinew = ii;
				uint32_t currsize = iisize;
				if ( iisize <= disjoints ) { iinew = getTopoMemory( ii, iisize, disjoints, currsize ); }//in array ii to be able to place at least disjoint elements
				ii = iinew;
				iisize = currsize;
				ii[disjoints] = jid;
				disjoints++;
	}	}	}
	nn[FIRST_FOURTHNN_NBOR] = disjoints;
cout << "\t" << disjoints << endl;

	//free unnecessary memory trailing in ii
	ii = cutTopoTail( ii, iisize, disjoints );
	if ( ii == NULL ) { cout << "Invalid cutTail operation" << endl; return; }

	//carry over
	target->n = nn;
	target->ids = ii;
}


int topoHdl::analyse_grain2d_construct_kshell( uint32_t kmax, uint32_t first, uint32_t last )
{
	//for each completely disjoint data set and grain analyze disjointly the neighbor topology
	//MK::checking for memory is essential as, grain id environment in each data set for each grain is stored as many times as there are n-th order neighbors!
	//##CURRENTLY ALL TIME STEPS PARTICIPATION
	if (first < this->firstloaded || last > this->lastloaded || first > last || last < 0) {
		cout << "ERR:Invalid file descriptors for input data!" << endl; return 0;
	}


	for (uint32_t mys = 0; mys < mydatasets.size(); mys++) {
		//how many grains in this data set?
		uint32_t ngr = mydatasize[mys];

		for (uint32_t g = 0; g < ngr; g++) {
			nbortopoP nt = NULL;
			nt = new struct nbortopo[1];
			myMemGuard = myMemGuard + sizeof(nbortopo);
			if ( nt == NULL ) { cout << "Error in allocating memory at construct kshell mys/g = " << mys << "\t" << g << endl; return 0; }
			nt->n = NULL;
			nt->ids = NULL;
			//nt->disori = NULL;

//cout << mys << "\t" << g << "...." << endl;

			find_neighbors( mys, g, nt );

			if ( nt->n == NULL || nt->ids == NULL ) { 
				cout << "Error in allocating memory during find_neighbors at construct kshell mys/g = "<< mys << "\t" << g << endl;
				return 0;
			}

cout << mys << "\t" << g << "\t" << myMemGuard << "\t" << nt->n[0] << ";" << nt->n[1] << ";" <<nt->n[2] << ";" <<nt->n[3] << "||";
//##DEBUG
//for ( int kk = nt->n[0]; kk < nt->n[1]; kk++ ) { cout << nt->ids[kk] << ";"; } cout << "||";
//for ( int kk = nt->n[1]; kk < nt->n[2]; kk++ ) { cout << nt->ids[kk] << ";"; } cout << "||";
//for ( int kk = nt->n[2]; kk < nt->n[3]; kk++ ) { cout << nt->ids[kk] << ";"; } cout << endl;

			mydata[mys][g].kshell = nt;
		}
	}
	cout << this->myRank << " is at " << myMemGuard/1024/1024 << " MB when the kshell construction is now complete." << endl;

	return 1;
}


void topoHdl::analyse_grain2d_kshellevolution(uint32_t kmax, uint32_t first, uint32_t last)
{
	if (kmax > this->kshellmax ) { cout << "Invalid kshell max parameter!" << endl; return; }

	if (first < this->firstloaded || last > this->lastloaded || first > last || last < 0) {
		cout << "ERR:Invalid file descriptors for input data!" << endl; return;
	}

	//##MK::implement a check if I at all have the f % offset dataset...
//if ( this->nRanks > 1 ) { cout << "Currently executable on MASTER only." << endl; return; }

//##MK::this for (int f ... for loop can be changed to for (int mys = 0; mys < mydatasets.size(); mys++ ... 
//##MK::so then first the nodes perform locally the analysis, then either a) write to MPI_File directly, or b) MPISend to MASTER who aggregates a large matrix, i.e. kshellmax*nneighborsmax*((last-first)/offset +1), or c) MASTER hierarchy
//##currently method a)
	MPI_File mpiiohdl_nnbors, mpiiohdl_avfaces;
	MPI_Status mpiiosta_nnbors, mpiiosta_avfaces;

	//create C-consistent file name for MPI I/O
	ostringstream nnbors_fn, avfaces_fn;
	nnbors_fn << "topotrace.genaw.nnbors." << jobid << ".first." << first << ".offset." << offsetloaded << ".last." << last << ".bin";
	avfaces_fn << "topotrace.genaw.avfaces." << jobid << ".first." << first << ".offset." << offsetloaded << ".last." << last << ".bin";
	char* nnbors_fname = new char[nnbors_fn.str().size() + 1];		strcpy( nnbors_fname, nnbors_fn.str().c_str() );
	char* avfaces_fname = new char[avfaces_fn.str().size() + 1 ];	strcpy( avfaces_fname, avfaces_fn.str().c_str() );

	//open currently via method a)
	MPI_File_open(MPI_COMM_WORLD, nnbors_fname, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &mpiiohdl_nnbors);
	MPI_File_open(MPI_COMM_WORLD, avfaces_fname, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &mpiiohdl_avfaces);

	__int64 toffset_nnbors = 0;
	__int64 toffset_avfaces = 0;

	uint32_t arrlen = this->kshellmax * this->nneighborsmax;

cout << "IMPLICIT ADRESSING KSHELLMAX = " << this->kshellmax << " nneighborsmax = " << this->nneighborsmax << endl;

	for ( uint32_t f = first; f <= last; f += offsetloaded) {
		//who has the dataset makes the analysis and writes into the MPI file currently via method a)
		uint32_t fid = (f/offsetloaded) - (firstloaded/offsetloaded);

		if ( myRank == assignment[fid].whichnode ) { // I have this dataset, I analyze and write to the collectively shared file

			long* nnbors = NULL;
			nnbors = new long[arrlen];
			if ( nnbors == NULL ) { cout << "Allocation memory for nnbors in kshellevolution failed." << endl; return; }

			long* avfaces = NULL;
			avfaces = new long[arrlen];
			if ( avfaces == NULL ) { cout << "Allocation memory for avfaces in kshellevolution failed." << endl; return; }

			//long* disorienv = NULL;
			//disorienv = new long[arrlen];
			//if ( disorienv == NULL ) { cout << "Allocation memory in disorienv in kshellevolution failed." << endl; break; return; }

			for ( uint32_t i = 0; i < arrlen; i++ ) {
				nnbors[i] = 0;
				avfaces[i] = 0;
				//disorienv[i] = 0;
			}

			uint32_t topoclass;
			//shell = 0 filled first how many center grains of class?
			for (uint32_t g = 0; g < mydatasize[fid]; g++) {
				topoclass = mydata[fid][g].nfaces;

				nnbors[(0*nneighborsmax)+topoclass] += 1;
				avfaces[(0*nneighborsmax)+topoclass] += 1;
				//disorienv[(0*nneighborsmax)+topoclass] += 1;
			}

			for ( uint32_t shell = 1; shell < kshellmax; shell++ ) {
cout << "... processing fid " << fid << " ... shell " << shell << " ..." << endl;
				for ( uint32_t g = 0; g < mydatasize[fid]; g++ ) {
//cout << "\t\t\t" << g << endl;
					//uint32_t oid = mydata[fid][g].ori;
					//double ce1 = oripool[oid].bunge1;
					//double ce2 = oripool[oid].bunge2;
					//double ce3 = oripool[oid].bunge3;

					//get topological class of the g grain
					topoclass = mydata[fid][g].nfaces;
					nbortopoP ks = mydata[fid][g].kshell;

					//how many neighboring grains in the shell-th shell?
					nnbors[((shell)*this->nneighborsmax) + topoclass] += ks->n[shell];

					//how many faces do these neighbor add for the generalized Aboav Weaire?
					for ( uint32_t nb = ks->n[shell-1]; nb < ks->n[shell]; nb++ ) {
						uint32_t nid = ks->ids[nb];
						uint32_t nbtopoclass = mydata[fid][nid].nfaces;

						//uint32_t noid = mydata[fid][nid].ori;
						//double ne1 = oripool[noid].bunge1;
						//double ne2 = oripool[noid].bunge2;
						//double ne3 = oripool[noid].bunge3;

						//double disori = misorientationCubic ( ce1, ce2, ce3, ne1, ne2, ne3 );
						//if ( disori < LAGBtoHAGBTransition ) 
						//	disorienv[(shell*nneighborsmax)+topoclass] += 1;

						avfaces[(shell*nneighborsmax)+topoclass] += nbtopoclass; //mind that center grains are in topoclass not nbtopoclass!
					}
				} //all grains
			} //all shells

			double* nnbors_res = NULL;
			nnbors_res = new double[arrlen];
			if ( nnbors_res == NULL ) { cout << "Failed alloc in nnbors_res." << endl; return; }
			double* avfaces_res = NULL;
			avfaces_res = new double[arrlen];
			if ( avfaces_res == NULL ) { cout << "Failed alloc in avfaces_res." << endl; return; }
			//double* disorienv_res = NULL;
			//disorienv_res = new double[arrlen];
			//if ( disorienv_res == NULL ) { cout << "Failed alloc in disorienv_res." << endl; break; return; }
			for ( uint32_t i = 0; i < arrlen; i++) {
				nnbors_res[i] = 0.0;
				avfaces_res[i] = 0.0;
				//disorienv_res[i] = 0.0;
			}

			//postprocess results
			//shell = 0 is probability density of grains in class tc
			for ( uint32_t tc = 0; tc < this->nneighborsmax; tc++ ) {
cout << " there are " << mydatasize[fid] << " grains in total in fid " << fid << endl;
				if ( mydatasize[fid] > 0 ) {
					nnbors_res[(0*nneighborsmax)+tc] = ( (double) nnbors[(0*nneighborsmax)+tc] / (double) mydatasize[fid] );
					avfaces_res[(0*nneighborsmax)+tc] = ( (double) avfaces[(0*nneighborsmax)+tc] / (double) mydatasize[fid] );
					//disori...
				}
			}

			for ( uint32_t shell = 1; shell < kshellmax; shell++ ) {
cout << "... postprocessing fid " << fid << " ... shell " << shell << "..." << endl;
				for ( uint32_t tc = 0; tc < this->nneighborsmax; tc++ ) {
					//how many neighbors on average in the shell-th shell for center grains in topoclass tc? - Lewis law
					if ( nnbors[(0*nneighborsmax)+tc] > 0 ) {
						nnbors_res[(shell*nneighborsmax)+tc] = ( (double) nnbors[(shell*nneighborsmax)+tc] / (double) nnbors[(0*nneighborsmax)+tc] );
					}

					//how many faces on average for neighbors in the shell-th shell of center grain in topoclass tc for center grains of topological class tc? - Aboav Weaire law
					if ( nnbors[(shell*nneighborsmax)+tc] > 0 ) {
						avfaces_res[(shell*nneighborsmax)+tc] = ( (double) avfaces[(shell*nneighborsmax)+tc] / (double) nnbors[(shell*nneighborsmax)+tc] );
					}

					//how many LAGB neighbors in the shell-th shell of a center grain in class tc in relation to all neighbors?
					//##MK::incorporate in above
					//if ( nnbors[(shell*nneighborsmax)+tc] > 0 ) {
					//	disorienv_res[(shell*nneighborsmax)+tc] = ( (double) disorienv[(shell*nneighborsmax)+tc] / (double) nnbors[(shell*nneighborsmax)+tc] );
					//}
				}
			}

			//##MK::later simply overwrite to zero again to safe delete and reallocation effort
			delete [] nnbors;
			delete [] avfaces;
			//delete [] disorienv;


		//output results
/*##DEBUGcout << endl << endl << "TC;Shell;nnbors;avfaces\n"; //;disorienv\n";
for ( uint32_t i = 0; i < arrlen; i++) {
	uint32_t shell = i / nneighborsmax;		uint32_t tc = i % nneighborsmax;
	if ( i != ((shell*nneighborsmax)+tc) ) { cout << "Invalid indexing scheme!!" << endl; return; }
	cout << tc << ";" << shell << ";" << setprecision(8) << nnbors_res[i] << ";" << setprecision(8) << avfaces_res[i] << "\n"; //";" << disorienv_res[i] << "\n";
} cout << endl << endl;*/

			toffset_nnbors = (fid - (first/offsetloaded)) * arrlen * 8; //sizeof(double);
			toffset_avfaces = (fid - (first/offsetloaded)) * arrlen * 8; //sizeof(double);

			MPI_File_write_at( mpiiohdl_nnbors, toffset_nnbors, nnbors_res, arrlen, MPI_DOUBLE, &mpiiosta_nnbors );
			MPI_File_write_at( mpiiohdl_avfaces, toffset_avfaces, avfaces_res, arrlen, MPI_DOUBLE, &mpiiosta_avfaces );


			gen_aboav_weaire_avfaces.push_back( nnbors_res );
			gen_aboav_weaire_nnbors.push_back( avfaces_res );
			//disori_environment.push_back ( disorienv_res );
		}
		//else {}, nothing someone else meanwhile takes care of fid in parallel on the MPI files

	} //next data element, all are disjoint and can be issued trivial in parallel as only read no writes

	delete [] nnbors_fname;
	delete [] avfaces_fname;


	MPI_Barrier( MPI_COMM_WORLD );
	MPI_File_close(&mpiiohdl_nnbors);
	MPI_File_close(&mpiiohdl_avfaces);
}
