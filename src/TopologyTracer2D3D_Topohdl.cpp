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

//for autodetection of RVE extent of GraGLeS 3d data
#include <boost/filesystem.hpp>
using namespace boost::filesystem;


//#ifdef USE_POLYTRI
#include "thirdparty/PolyTri/PolyTri_Geometry.h"
//#endif

//C header for querying file size
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

//#ifdef __GNUC__
//typedef long long __int64;
//#endif


using namespace std;


//##MK::dummy replace with correct data when handling
#define OMP_GET_NUM_THREADS					1
#define OMP_GET_THREAD_NUM					0
//###OMP screen multithreaded the data to build local tree-like hash and reserve local memory for grains
//##OMP every thread does, but as regions are a complete non-overlapping coverage of space
//no grain will be assigned twice a thread
//ompindividually_scanGrainsAndHash( buf, it->ng );

inline bool SortDoubleAscending(const double &dbl1 , const double &dbl2) {
	return dbl1 < dbl2;
}

inline bool SortGIDAscending( const UDS_Grain3DInfo &aa1, const UDS_Grain3DInfo &aa2) {
	return aa1.gid < aa2.gid;
}

inline bool SortSQKeyDistAscending(const sqkey &aa1, const sqkey &aa2) {
	return aa1.dist < aa2.dist;
}

timeLogger::timeLogger(){
}

timeLogger::~timeLogger(){
	titles.clear();
	times.clear();
}

void timeLogger::logev(const string title, double time){
	titles.push_back(title);
	times.push_back(time);
}


unsigned int timeLogger::get_nentries( void ){
	if ( this->titles.size() != times.size() )
		return 0;
	else
		return titles.size();
}


MemRegion::MemRegion( void )
{
	mys = NULL;
	intID = NOT_YET_KNOWN;

	Geometry.xmi = 0.0;
	Geometry.xmx = 0.0;
	Geometry.ymi = 0.0;
	Geometry.ymx = 0.0;
	Geometry.zmi = 0.0;
	Geometry.zmx = 0.0;
}


MemRegion::~MemRegion( void )
{
	//GrainBucket cleared before
	
	//however a vector of objects grains is not standard such we need to delete pieces of information regarding faces manually
	//vector<Grain>::iterator g;
	//for ( g = GrainBucket.begin(); g != GrainBucket.end(); ++g ) {
	//    delete [] g->neighbors;
	//}
	unsigned int ngr = GrainBucket.size();
	for ( unsigned int g = 0; g < ngr; ++g ) {
		if (GrainBucket.at(g).neighbors != NULL) {
			delete [] GrainBucket[g].neighbors;
			GrainBucket[g].neighbors = NULL;
		}
	}
	//MK::one MUST NOT delete mys as it is a back-references!
}


bool MemRegion::writeLookupTable( unsigned int extGID, unsigned int mr, unsigned int lid, unsigned int location)
{
	//writes into the lookup table in bucket idrange lid that the grain with the external simulation id extGID it is located in MemRegion->GrainBucket[location]
	//this enables to reidentify with the LookupTable arbitrary and non-contiguous numeral names of grains without the necessity to scan the entire set of grains
	//to search desperately for only one ID, i.e. it 
	//eliminates --- once having the GrainHash information read --- all naive grain id search operations on the grainbucket
	unsigned int ncand = this->mys->LookupTable[lid].size;
	GrainHash* idbucket = this->mys->LookupTable[lid].bucket; //MK::Lookup table is shared across threads! 

	//therefore writing to idbucket can cause datarace!

	bool status = false;
	for ( unsigned int c = 0; c < ncand; c++ ) {
		if ( extGID == idbucket[c].gid ) { //found the grain
			idbucket[c].idx = location;
			status = true;
			break;
		}
	}

	return status;
}


unsigned int MemRegion::queryMemRegionViaLookupTable( unsigned int extGID, unsigned int lid )
{
	//MK::get the memory region in which grain with extGID should be stored but scan only a particular ID range to safe search time!
	unsigned int mres = NOT_YET_FOUND;
	unsigned int ncand = this->mys->LookupTable[lid].size;
	GrainHash* idbucket = this->mys->LookupTable[lid].bucket;
	for ( unsigned int c = 0; c < ncand; c++ ) {
		if ( extGID == idbucket[c].gid ) { //found the grain
			mres = idbucket[c].mrg_idx;
			break;
		}
	}

	return mres;
}



void MemRegion::addFaceData( double size, unsigned int refid, unsigned int nborid, unsigned int reflid, unsigned int nborlid )
{
	//refid refers to the grain to which grain nborid is a neighbor

	//scan only the grains in the ID bucket
	//unsigned int refm = this->intID;				//the memory region in which grain refid is stored
	unsigned int refidx = NOT_YET_FOUND;			//the index in the container
	unsigned int ngl = UNKNOWN_ID;
	GrainHash* thebucket = NULL;
	thebucket = this->mys->LookupTable[reflid].bucket;
	ngl = this->mys->LookupTable[reflid].size;
	for ( unsigned int gl = 0; gl < ngl; ++gl ) { //MK::was g++ and == refid logic
		if (thebucket[gl].gid != refid)
			continue;

		//not continued means item was found
		refidx = thebucket[gl].idx;
		break;
	}

	if ( refidx == NOT_YET_FOUND ) { cerr << "ERROR::Worker " << this->mys->mytopohdl->get_Rank() << " I cannot find a grain named refid in my local data! (refid/nborid/reflid/nborlid = )" << refid << ";" << nborid << ";" << reflid << ";" << nborlid << std::endl; return; }

	//location of grain refid which has this face was found, now instantiate in data container describing the reference grain
	unsigned int nextFreeSlot = this->GrainBucket[refidx].nfaces_identified;
	if ( nextFreeSlot >= this->GrainBucket[refidx].nfaces) { cerr << "ERROR::Worker " << this->mys->mytopohdl->get_Rank() << " inconsistent attempt detected to add more faces to a grain with name refid than there should be existing (refid/nborid/reflid/nborlid = )" << refid << ";" << nborid << ";" << reflid << ";" << nborlid << endl; return; }
	this->GrainBucket[refidx].neighbors[nextFreeSlot].size = size;
	this->GrainBucket[refidx].neighbors[nextFreeSlot].extNBID = nborid;

	//in which MemRegion and position can I find this neighbor for later and quick access?
	unsigned int nborm = NOT_YET_FOUND;
	unsigned int nboridx = NOT_YET_FOUND;
	thebucket = NULL;
	ngl = UNKNOWN_ID;

	thebucket = this->mys->LookupTable[nborlid].bucket;
	ngl = this->mys->LookupTable[nborlid].size;

	//##MK::DEBUBif ( refid == 48930 && nborid == 27028 ) cout << "DUMMY nborlidsize " <<  this->mys->LookupTable[nborlid].size << std::endl;
	//if ( refid == 48930 && nborid == 27028 ) std::cerr << "\t\tc/gid/mrgidx/idx " << c << ";" << thebucket[c].gid << ";" << thebucket[c].mrg_idx << ";" << thebucket[c].idx << endl;

	for ( unsigned int c = 0; c < ngl; ++c ) { //MK::was c++ and ==nborid logic
		if (thebucket[c].gid != nborid)
			continue;
		
		//not continued means found
		nborm = thebucket[c].mrg_idx;
		nboridx = thebucket[c].idx;
		break;
	}
	if ( nborm == NOT_YET_FOUND || nboridx == NOT_YET_FOUND ) { cerr << "ERROR::Worker " << this->mys->mytopohdl->get_Rank() << " I do not find the memregion and or the index of nborid again!\t\trefid/nborid/reflid/nborlid = " << refid << ";" << nborid << ";" << reflid << ";" << nborlid << std::endl; return; }

	this->GrainBucket[refidx].neighbors[nextFreeSlot].mrg_idx = nborm;
	this->GrainBucket[refidx].neighbors[nextFreeSlot].idx = nboridx;

	//add size
	this->GrainBucket[refidx].nfaces_identified = this->GrainBucket[refidx].nfaces_identified + 1;

//##MK::DUMMYif ( refid == 4329 || nborid == 4329 ) { cout << "\t\tNFacesIdentified " << this->GrainBucket[refidx].nfaces_identified << " refid/nborid = " << refid << ";" << nborid << std::endl; }
//cout << "\t\t\t-->refid/nborid/reflid/nborlid = " << refid << ";" << nborid << ";" << reflid << ";" << nborlid << endl;
}


unsigned int MemRegion::InWhichMemRegionIsTheGrain( unsigned int id, unsigned int lid )
{
	if ( this->mys->LookupTable[lid].size > 0 ) { //most likely case
		GrainHash* thebucket = this->mys->LookupTable[lid].bucket;
		unsigned int ngl = this->mys->LookupTable[lid].size;
		for (unsigned int gl = 0; gl < ngl; ++gl) {
			if (thebucket[gl].gid != id)
				continue;

			//not continued so found
			return thebucket[gl].mrg_idx;
		}

		//no return yet, but IDBucket had items is inconsistent!
		return NOT_YET_FOUND;
	}
	else {
		return NOT_YET_FOUND;
	}
}



bool MemRegion::omp_storeFaces( unsigned int this_mreg )
{
	//CALLED FROM WITHIN PARALLEL REGION
	//self-consistency
	//MK:: an essential convention is that Grain IDs run from 1;Ngrains a grain with extGID == 0 specifies the virtual grain
	//i.e. the domain, as such a face pair with one grain 0 is a face to the simulation domain boundary
	//nfaces of the grains accounts for all faces!

	if ( this_mreg != this->intID ) {
		cerr << "ERROR::Worker " << this->mys->mytopohdl->get_Rank() << " a work-sharing organizational inconsistency was detected in that this_mreg and this->intID are different in thread " << OMP_GET_THREAD_NUM << endl; return false;
	}

	//unsigned int mtarget = FIRST_REGION;
	unsigned int luidA = 0;
	unsigned int luidB = 0;
	unsigned int nftotal = this->mys->mytopohdl->fcbufsize;
	//MK::in the binary file, faces between grains are always stored only once!
	//in the TopologyTracer, however we would like to have them for each grain to simplify perimeter-neighbor-property correlation computations!

	MPI_3DFaceIO* thedata = this->mys->mytopohdl->fcbuf;

	for (unsigned int f = 0; f < nftotal; f++ ) {
		//check for the first grain
		luidA = thedata[f].gA / Settings::LookupMaxGrainIDRange; //gA > gB by definition
		luidB = thedata[f].gB / Settings::LookupMaxGrainIDRange;

//cout << "f/gA/gB = " << f << ";" << thedata[f].gA << ";" << thedata[f].gB << "---" << InWhichMemRegionIsTheGrain( thedata[f].gA, luidA ) << ";" << InWhichMemRegionIsTheGrain( thedata[f].gB, luidB ) << endl;
		//we do not addFaces the zero grain (THE_DOMAIN_ITSELF) has with grains A, B both != 0 and >=1
		if ( thedata[f].gA != THE_DOMAIN_ITSELF ) {
			if ( InWhichMemRegionIsTheGrain( thedata[f].gA, luidA ) == this_mreg ) { //I host a grain with this ID in my data
				this->addFaceData( thedata[f].size, thedata[f].gA, thedata[f].gB, luidA, luidB ); //ID of the neighbor, not of grain gA
				//if ( thedata[f].gA == 4329 || thedata[f].gB == 4329 ) cout << "\t\tAdd/gA/gB/" << thedata[f].gA << ";" << thedata[f].gB << std::endl;
			}
		}
		if ( thedata[f].gB != THE_DOMAIN_ITSELF ) {
			if ( InWhichMemRegionIsTheGrain( thedata[f].gB, luidB ) == this_mreg ) { //I host a grain with this ID in my data
				this->addFaceData( thedata[f].size, thedata[f].gB, thedata[f].gA, luidB, luidA ); //ID vice versa grain gA, not gB is the neighbor
				//if ( thedata[f].gA == 4329 || thedata[f].gB == 4329 ) cout << "\t\tAdd/gB/gA/" << thedata[f].gB << ";" << thedata[f].gA << std::endl;
			}
		}
		//MK::it is possible that one of the grain of a pair is not in the same MemRegion!

		//##MK::occasional promptif ( f % 1000 == 0 ) { cout << "\t\t\tFace " << f << " written in the database" << endl; }

	} //all faces from MPIIO file processed

	return checkFaceReassignmentConsistency();
}


bool MemRegion::omp_storeGrainData( void )
{
	//CALLED FROM WITHIN THREAD PARALLEL REGION
	//register the "zero grain" which encodes THE_DOMAIN_ITSELF i.e. enables that grains have defined faces with boundaries
	unsigned int mtarget = FIRST_REGION;
	unsigned int luid = THE_DOMAIN_ITSELF / Settings::LookupMaxGrainIDRange;
	if ( mtarget == intID ) { 	//is being executed in only one thread for each snapshot, namely thread FIRST_REGION
			struct Grain domainbnd;
			domainbnd.extGID = THE_DOMAIN_ITSELF;
			domainbnd.nfaces = 0; //the neighbors of the domain, I the boundary grains are not included
			domainbnd.boundary = BOUNDARY_CONTACT;
			domainbnd.size = 0.0;
			domainbnd.surface = 0.0;
			domainbnd.phi1 = 0.0;
			domainbnd.Phi = 0.0;
			domainbnd.phi2 = 0.0;
			domainbnd.x = 0.0;
			domainbnd.y = 0.0;
			domainbnd.z = 0.0;
			domainbnd.see = 0.0;
			domainbnd.nfaces_identified = 0;
			this->GrainBucket.push_back( domainbnd );
			this->GrainBucket[this->GrainBucket.size()-1].neighbors = NULL; //no tracking of faces for domain
			if ( writeLookupTable( domainbnd.extGID, mtarget, luid, (this->GrainBucket.size()-1) ) == false ) { cerr << "ERROR::TrackingParallel::OmpStoreGrain write back error in LookupTable for ZERO GRAIN!" << endl; return false; }
	}

	//process through all the grains and pick those necessary...
	unsigned int ngr = this->mys->mytopohdl->grbufsize;
	MPI_3DGrainIO* thedata = this->mys->mytopohdl->grbuf;

	for ( unsigned int gr = 0; gr < ngr; gr++ ) {
		luid = thedata[gr].id / Settings::LookupMaxGrainIDRange; //##MK::performance problem?
		mtarget = queryMemRegionViaLookupTable( thedata[gr].id, luid );
		if ( mtarget == NOT_YET_FOUND ) { 
			cerr << "Data organization inconsistency in omp_storeGrainData for grain " << gr << " with extID = " << thedata[gr].id << endl;
			return false;
		}

		if ( mtarget == intID ) { //I take care about the grain so then data transfer
			struct Grain agr;
			agr.extGID = thedata[gr].id;
			agr.nfaces = thedata[gr].NeighbourCount;
			agr.boundary = thedata[gr].intersectsBoundaryGrain;
			agr.size = thedata[gr].size;
			agr.surface = thedata[gr].surfaceArea;
			agr.phi1 = thedata[gr].phi1;
			agr.Phi = thedata[gr].Phi;
			agr.phi2 = thedata[gr].phi2;
			agr.x = thedata[gr].x;
			agr.y = thedata[gr].y;
			agr.z = thedata[gr].z;
			agr.see = thedata[gr].BulkEnergy;
			agr.nfaces_identified = 0; //##MK::no faces were linked to the grain yet successfully

			this->GrainBucket.push_back( agr );

//##MK::occasional promptif ( agr.extGID % 10000 == 0 ) { cerr << "\t\t\tGrain " << agr.extGID << " written in the database, check phi1 = " << this->GrainBucket[this->GrainBucket.size()-1].phi1 << endl; }

			//get memory to store the faces
			Face* nb = NULL;
			nb = (Face*) new Face[agr.nfaces];
			this->GrainBucket[this->GrainBucket.size()-1].neighbors = NULL;
			this->GrainBucket[this->GrainBucket.size()-1].neighbors = nb;

			if ( writeLookupTable( agr.extGID, mtarget, luid, (this->GrainBucket.size()-1) ) == false ) { 
				cerr << "Write back error in LookupTable" << endl; return false;
			}
		}

	} //proceed with next possible candidate in much MemRegion

	return true;
}


bool MemRegion::omp_allocGrainMemory( unsigned int ngr_per_mr )
{
	//CALLED FORM WITHIN PARALLEL REGION
	this->GrainBucket.reserve( ngr_per_mr );
	//##exception handling

	return true;
}


unsigned int* MemRegion::nborid_malloc( unsigned int size )
{
	if ( size < 1 ) return NULL;

	unsigned int* nb = NULL;
	nb = new unsigned int[size];
	if ( nb == NULL ) {
		cerr << "ERROR::TrackingParallel::Nborid_malloc allocation error during nborid_malloc!" << endl;
		return NULL;
	}
	for ( unsigned int i = 0; i < size; ++i ) {
		nb[i] = REASONABLE_CLEARING_VALUE;
	}

	return nb;
}


void MemRegion::nborid_delete( unsigned int* nb, ksmetaP smeta )
{
	if ( nb != NULL ) {
		delete [] nb;
		nb = NULL;
	}

	if ( smeta != NULL ) {
		delete [] smeta;
		smeta = NULL;
	}
}


unsigned int* MemRegion::nborid_expand( unsigned int* nb_old, unsigned int size_old, unsigned int size_new )
{
	if ( size_new <= size_old ) { //catch no increase and contraction attempts
		return NULL;
	}

	//allocate a larger buffer
	unsigned int* nb_new = NULL;
	nb_new = new unsigned int[size_new];
	if ( nb_new == NULL ) {
		cerr << "ERROR::TrackingParallel::Nborid_expand allocation error during nborid_malloc!" << endl;
		return NULL;
	}

	//copy content from nb_old over
	for ( unsigned int i = 0; i < size_old; ++i ) {
		nb_new[i] = nb_old[i];
	}

	//fill remaining with default values
	for ( unsigned int i = size_old; i < size_new; ++i ) {
		nb_new[i] = REASONABLE_CLEARING_VALUE;
	}

	//delete nb_old properly
	delete [] nb_old;
	nb_old = NULL;

	return nb_new;
}


unsigned int* MemRegion::nborid_resize( unsigned int* nb_old, unsigned int size_old, unsigned int copylimit )
{
	if ( copylimit-1 > size_old ) { //nb_old[copylimit] will be accessed due to inclusive reads
		cerr << "ERROR::TrackingParallel::Nborid_resize allocation error during nborid_resize!" << endl;
		return NULL;
	}

	//copys content from nb_old en block [0,copylimit-1] into buffer with size of exactly copylimit
	unsigned int size_new = copylimit;
	unsigned int* nb_new = NULL;
	nb_new = (unsigned int*) new unsigned int[size_new];
	if ( nb_new == NULL ) { 
		cerr << "ERROR::TrackingParallel::Nborid_resize allocation error during nborid_resize nb_new!" << endl;
		return NULL;
	}

	for ( unsigned int i = 0; i <= copylimit; ++i ) { //will be optimized by compiler
		nb_new[i] = nb_old[i];
	}

	return nb_new;
}


bool MemRegion::nborid_add_noduplicates( unsigned int* nbid, unsigned int last_to_compare, unsigned int gid )
{
	//returns true when the grain id gid is not in nbid and places the gid then
	//else returns false as gid was a duplicate without the necessity of a placement

	//MK::caller has assured that reading up to last_to_compare stays within allocated range of nbid and there is
	//space to add another item
	//bool foundduplicate = false;
	for ( unsigned int i = 0; i <= last_to_compare; i++ ) { //50:50 duplicates expected
		if ( nbid[i] != gid ) continue;

		//a duplicate in the range [0, last_to_compare]
		return false;
	}

	//no duplicate found to returned, so not a duplicate
	nbid[last_to_compare+1] = gid;
	return true;
}


void MemRegion::find_gid( unsigned int ggid, unsigned int* wmr, unsigned int* wid )
{
	if ( ggid == THE_DOMAIN_ITSELF )
		return;

	unsigned int cand_luid = ggid / Settings::LookupMaxGrainIDRange;
	GrainHash* thebucket = this->mys->LookupTable[cand_luid].bucket;

	unsigned int cand_lusize = this->mys->LookupTable[cand_luid].size;
	for ( unsigned int g = 0; g < cand_lusize; g++ ) { //scan only the most interesting candidates
		if ( thebucket[g].gid != ggid )
			continue;

		//grain id found!
		*wmr = thebucket[g].mrg_idx;
		*wid = thebucket[g].idx;
		return;
	}
}


bool MemRegion::omp_knn( unsigned int kkmax )
{
	//##MK::improve exception handling!!
	//loop over the grains in my bucket
	//MK::the tracking of n-th nearest neighbors in large networks requires eventually the tracking of thousands of IDs per grain, hence this function
	//processes one grain after another, while other OpenMP threads work on another queue of grains in threadlocal memory and the same timestep
	//nbIDs carries an implicit array of all neighbors in the 0-th, 1-th, 2-th, 3-th, ..., n-th shell
	//nbShell[k] identifies all nbIDs indices (inclusive reading) between which IDs of only the k-th shell are found, nbcnt assures the filling of the buffer
	bool stillgood = true;
	ksmetaP nbShell = NULL; //guides the program by detailing which IDs on nbIDs belong to which (k-th) shell
	//unsigned int* nbIDs = NULL; //keeps these grain IDs
	vector<unsigned int> nbIDs;
	unsigned int nbIDs_size = 0;
	unsigned int nbcnt = 0;

	vector<Grain>::iterator itG = this->GrainBucket.begin(); //filter zero grain
	while ( itG->extGID == 0 )
		itG++;
	for (  ; itG != this->GrainBucket.end(); itG++ ) {
		double timer = MPI_Wtime();
		nbcnt = 0;
		nbShell = (ksmetaP) new ksmeta[1+kkmax];
		if ( nbShell == NULL ) { std::cerr << "Allocation error prior 0-th shell (nbShell)!" << std::endl; stillgood = false; break; }
		//nbIDs = nborid_malloc( DEFAULT_KNN_BUFFERSIZE );
		nbIDs.clear();
		//if ( nbIDs == NULL ) { std::cerr << "Allocation error prior 0-th shell (nbIDs)!" << std::endl; stillgood = false; delete [] nbShell; break; }
		//nbIDs_size = DEFAULT_KNN_BUFFERSIZE;

		//add the target itself in the 0-th shell
		nbShell[0].first = 0;

		//nbIDs[nbcnt] = itG->extGID;
		nbIDs.push_back(itG->extGID);
		nbcnt++;

		nbShell[0].last = nbIDs.size() - 1; //nbcnt - 1; //MK::-1 because [first,last] inclusive reading of nbIDs
		//if ( (nbcnt-1) != nbShell[0].last ) cerr << "Logical error 0" << std::endl; //##MK::dummy
		

///////////////////////////////////////////THIS SECTION CAN BE ELIMINATED WHEN SETTING k = 1 IN THE ANALYSIS OF THE n-th neighbors!
		//add direct neighbors of the faces in the 1-th shell
		if ( kkmax >= 1 ) {
			if ( itG->nfaces > 0 ) {
				/*if ( (1 + itG->nfaces + 2) > nbIDs_size ) { //assure properly sized array before writing
					unsigned int* switchID = nborid_expand( nbIDs, nbIDs_size, (nbIDs_size + 2*itG->nfaces) ); //copies over and deletes nbIDs
					if ( switchID == NULL ) { cerr << "Allocation error prior 1-th shell!" << std::endl; stillgood = false; delete [] nbShell; delete [] nbIDs; break; }
					nbIDs = switchID;
					nbIDs_size = nbIDs_size + 2*itG->nfaces;
				}*/

				nbShell[1].first = nbShell[0].last + 1;

				for ( unsigned int n = 0; n < itG->nfaces; ++n ) { 
					//nbIDs[nbcnt] = itG->neighbors[n].extNBID;
					nbIDs.push_back( itG->neighbors[n].extNBID );
					//nbcnt++;
				}

				nbShell[1].last = nbIDs.size() - 1; //nbcnt - 1;
				//if ( (nbcnt-1) != nbShell[1].last ) cerr << "Logical error 1" << std::endl; //##MK::dummy
			}
		}

//////////////////////////////////////////

		//add n-th neighbors iteratively
		for ( unsigned int k = 2; k <= kkmax; k++ ) {
			nbShell[k].first = nbShell[k-1].last + 1; //##MK::potentially unsafe when allocation in subsequent for loop fails! when no one is found

			for ( unsigned int cand = nbShell[k-1].first; cand <= nbShell[k-1].last; cand++ ) {	//read inclusive k-1-th shell ids, i.e. [first,last] accept unknowns
				unsigned int whichmr = NOT_ASSIGNED_YET;
				unsigned int whichidx = NOT_ASSIGNED_YET;

				this->find_gid( nbIDs[cand], &whichmr, &whichidx );

				//##MK::dummy, cannot incur if database has been initialized correctly
				if ( whichmr == NOT_ASSIGNED_YET || whichidx == NOT_ASSIGNED_YET ) { cerr << "ERROR::LookupTableMissingEntry gid/whichmr/whichidx " << nbIDs[cand] << ";" << whichmr << ";" << whichidx << std::endl; continue; } //skip this neighbor...

				//now get the 1-order neighbors of the candidate and inspect for whom of them to add
				Face* TestTheseNeighbors = this->mys->mrg[whichmr]->GrainBucket[whichidx].neighbors; //for sure will invoke remote accesses...
				unsigned int HowManyOfThese = this->mys->mrg[whichmr]->GrainBucket[whichidx].nfaces; //most neighbors are located in thread-local memory!

				if ( HowManyOfThese > 0 ) {
					/*unsigned int ExpansionTest = nbcnt + HowManyOfThese + 2;
					if ( ExpansionTest > nbIDs_size ) { //expansion is necessary
						
						if ( ExpansionTest < (nbIDs_size + DEFAULT_KNN_RECACHE) )	//test whether better to get directly a larger expansion
							ExpansionTest = nbIDs_size + DEFAULT_KNN_RECACHE;
						else														//than an incremental one...
							ExpansionTest = nbIDs_size + 2*HowManyOfThese;

						unsigned int* switchID = nborid_expand( nbIDs, nbIDs_size, ExpansionTest ); //#MK::improve here by reallocating more!
						if ( switchID == NULL ) { cerr << "Allocation error prior " << k << "-th shell!" << std::endl; stillgood = false; break; }
						nbIDs = switchID;
						nbIDs_size = ExpansionTest;
					}*/
					for ( unsigned int thisnbor = 0; thisnbor < HowManyOfThese; thisnbor++ ) {
						//if ( nborid_add_noduplicates( nbIDs, nbShell[k-1].last, TestTheseNeighbors[thisnbor].extNBID ) == true )
						//	nbcnt++; //###in no duplicates nborid it must be make a comment why last works
						bool found = false;
						for ( unsigned int i = 0; i <= nbIDs.size(); i++ ) { //50:50 duplicates expectednbShell[k-1].last
							if ( nbIDs[i] == TestTheseNeighbors[thisnbor].extNBID || TestTheseNeighbors[thisnbor].extNBID == THE_DOMAIN_ITSELF ) {
								found = true;
								break;
							}
						}
						if ( found == false ) 
							nbIDs.push_back( TestTheseNeighbors[thisnbor].extNBID );
					}
				}
			}

			nbShell[k].last = nbIDs.size() - 1; //nbcnt - 1;
			//if ( (nbcnt-1) != nbShell[k].last ) cerr << "Logical error " << k << std::endl; //##MK::dummy
		}

		//unsigned int* theShell = nborid_resize( nbIDs, nbIDs_size, nbShell[kkmax].last );
		//delete [] nbIDs;
		//nbIDs = NULL;

		//##do something physically with these IDs ...
		cout << "kNN for grain/nfaces/nb_size/kkmaxlast/time " << itG->extGID << ";" << itG->nfaces << ";" << nbIDs_size << ";" << nbShell[kkmax].last << ";" << (double) MPI_Wtime() - timer << "\t\t"; //<< std::endl;
		for ( unsigned int kk = 0; kk <= kkmax; kk++ )	cout << nbShell[kk].first << "__" << nbShell[kk].last << ";";
		cout << endl;
		for ( unsigned int nn = 0; nn < nbIDs.size(); nn++) cout << nbIDs[nn] << ";";
		cout << endl;
		//##this is the dummy at the moment...


		//nborid_delete( theShell, nbShell );
		nbIDs.clear();
		if ( nbShell != NULL) { delete [] nbShell; nbShell = NULL; }
	}

	return stillgood;
}

#define NOT_TOUCHING	0x00
#define TOUCHING		0x01

bool MemRegion::omp_knn_elimbnd( unsigned int kkmax, std::ofstream &outputfile )
{
	//MK::the tracking of n-th nearest neighbors in large networks requires eventually the tracking of thousands of IDs per grain, hence this function
	//processes one grain after another, while other OpenMP threads work on another queue of grains in threadlocal memory and the same timestep
	//nbIDs carries an implicit array of all neighbors in the 0-th, 1-th, 2-th, 3-th, ..., n-th shell
	//nbShell[k] identifies all nbIDs indices (inclusive reading) between which IDs of only the k-th shell are found, nbcnt assures the filling of the buffer
	//this version stops the identification of a k-shell to a particular grain as soon as a grain with boundary contact is found

	char* TouchingDomainBnd = NULL;
	unsigned int nn = 1 + Settings::LargestGrainID;
	TouchingDomainBnd = (char*) new char[nn]; for (unsigned int c = 0; c < nn; ++c ) { TouchingDomainBnd[c] = NOT_TOUCHING; }
	if ( TouchingDomainBnd == NULL ) { cerr << "Unable to inspect boundary contact on thread " << this->intID << std::endl; return false; }
	//detect in local portion of dataset whether a grain at some point touches the boundary
	for ( unsigned int mr = 0; mr < this->mys->mrg.size(); mr++ ) {
		unsigned int mrngr = this->mys->mrg[mr]->GrainBucket.size();
		for ( unsigned int g = 0; g < mrngr; g++ ) {
			if ( this->mys->mrg[mr]->GrainBucket[g].boundary != BOUNDARY_CONTACT ) continue;
			TouchingDomainBnd[this->mys->mrg[mr]->GrainBucket[g].extGID] = TOUCHING;
		}
	}


	bool stillgood = true;
	bool bndcontact = false;
	ksmetaP nbShell = NULL; //guides the program by detailing which IDs on nbIDs belong to which (k-th) shell
	vector<unsigned int> nbIDs;
	unsigned int Statistics[3] = {0, 0, 0};

	vector<Grain>::iterator itG = this->GrainBucket.begin(); //filter zero grain
	while ( itG->extGID == 0 )
		itG++;
	for (  ; itG != this->GrainBucket.end(); itG++ ) {
		//double timer = MPI_Wtime();
		Statistics[KNN_STATS_TOTAL] += 1;
		bndcontact = false;

		if ( itG->boundary == BOUNDARY_CONTACT ) { Statistics[KNN_STATS_EXPELLED] += 1; bndcontact = true; continue; }

		nbShell = NULL;
		nbShell = (ksmetaP) new ksmeta[1+kkmax];
		if ( nbShell == NULL ) { cerr << "Allocation error prior 0-th shell (nbShell)!" << std::endl; stillgood = false; break; }
		nbIDs.clear();

		//add the target itself in the 0-th shell
		nbShell[0].first = 0;
		nbIDs.push_back(itG->extGID);
		nbShell[0].last = nbIDs.size() - 1;
		
		//add k-th neighbors iteratively (k = 1, 2, ..., kkmax)
		for ( unsigned int k = 1; k <= kkmax; k++ ) {
			if ( bndcontact == true ) { bndcontact = true; break; } //identify only neighbors without boundary contact

			nbShell[k].first = nbShell[k-1].last + 1; //##MK::potentially unsafe when allocation in subsequent for loop fails! when no one is found

			for ( unsigned int cand = nbShell[k-1].first; cand <= nbShell[k-1].last; cand++ ) {	//read inclusive k-1-th shell ids, i.e. [first,last] accept unknowns
				unsigned int whichmr = NOT_ASSIGNED_YET;
				unsigned int whichidx = NOT_ASSIGNED_YET;

				this->find_gid( nbIDs[cand], &whichmr, &whichidx );

				//##MK::dummy, cannot incur if database has been initialized correctly
				if ( whichmr == NOT_ASSIGNED_YET || whichidx == NOT_ASSIGNED_YET ) { cerr << "ERROR::LookupTableMissingEntry gid/whichmr/whichidx " << nbIDs[cand] << ";" << whichmr << ";" << whichidx << std::endl; continue; } //skip this neighbor...

				if ( this->mys->mrg[whichmr]->GrainBucket[whichidx].boundary == BOUNDARY_CONTACT ) { bndcontact = true; break; }

				//now get the 1-order neighbors of the candidate and inspect whom of them to add as a disjoint neighbor into nbIDs
				Face* TestTheseNeighbors = this->mys->mrg[whichmr]->GrainBucket[whichidx].neighbors; //for sure will invoke remote accesses...
				unsigned int HowManyOfThese = this->mys->mrg[whichmr]->GrainBucket[whichidx].nfaces; //most neighbors are located in thread-local memory!

//if ( itG->extGID == 6454463 ) { for (unsigned int tz = 0; tz < HowManyOfThese; tz++ ) { cout << "KShell " << k << ";" << tz << ";" << TestTheseNeighbors[tz].extNBID << std::endl; } }

				if ( HowManyOfThese > 0 ) {
					for ( unsigned int thisnbor = 0; thisnbor < HowManyOfThese; thisnbor++ ) {
						if ( bndcontact == true ) { break; }
						bool found = false;
						for ( unsigned int i = 0; i < nbIDs.size(); i++ ) { //was <= nbIDs.size() 50:50 duplicates expectednbShell[k-1].last

							//how to improve::either run the outmost loop only to k < kkmax and then run k <= kkmax during which testing for boundary contact of each candidate --> expensive when most grains are candidates...
							//allow in general the accounting of a neighborhood to the domain_itself --> currently this is disallowed!, however purely via faces file this may be inaccurate, therefore better add for all grains with boundary contact a grain the DOMAIN_ITSELF
							//test in each iteration for all candidates --> enormously expensive to test neighborhood
							//add neighborhood flag in neighbors to grains --> add memory and bools causes disadvantages padding

							if ( k == kkmax ) { //##MK::compromise solution to prevent excessive cache misses during neighborhood detection
								if ( TouchingDomainBnd[TestTheseNeighbors[thisnbor].extNBID] == TOUCHING ) {
//if  ( itG->extGID == 6454463 ) cout << TestTheseNeighbors[thisnbor].extNBID << " kicks the analysis because bnd!" << std::endl;
									bndcontact = true; break;
								}
							}
							if ( nbIDs[i] == TestTheseNeighbors[thisnbor].extNBID ) {
//if  ( itG->extGID == 6454463 ) cout << TestTheseNeighbors[thisnbor].extNBID << " kicks the analysis because nbIDs[i] known" << std::endl;
								found = true; break;
							}
							if ( TestTheseNeighbors[thisnbor].extNBID == THE_DOMAIN_ITSELF ) { //##MK::this does not test whether thisnbor touches the boundary
								bndcontact = true; break;
								//found = true; break;
							}
						}
						if ( found == false )
							nbIDs.push_back( TestTheseNeighbors[thisnbor].extNBID );
					}
				}
			}

			nbShell[k].last = nbIDs.size() - 1;
		}

		if ( bndcontact == true ) {
			Statistics[KNN_STATS_EXPELLED] += 1;
			nbIDs.clear();
			if ( nbShell != NULL) { delete [] nbShell; nbShell = NULL; }
			continue;
//cout << "EXPELLED::kNN for grain/nfaces/nb_size/kkmaxlast/time " << itG->extGID << ";" << itG->nfaces << ";" << nbIDs.size() << ";" << nbShell[kkmax].last << ";" << (double) MPI_Wtime() - timer << std::endl;
		}

		//do something physically with these IDs ...

//cout << "ACCEPTED::kNN for grain/nfaces/nb_size/kkmaxlast/time " << itG->extGID << ";" << itG->nfaces << ";" << nbIDs.size() << ";" << nbShell[kkmax].last << ";" << (double) MPI_Wtime() - timer << "\t\t" << itG->x << ";" << itG->y << std::endl;

		//double* kshell_perimeter = new double[kkmax];
		//for (unsigned int j = 0; j < kkmax; j++ ) { kshell_perimeter[j] = INVALID_LENGTH; } //j-th shell e.g. if kkmax=4 j=0,1,2,3
		kshell_props* kshell_meta = new kshell_props[kkmax];
		for ( unsigned int kj = 0; kj < kkmax; kj++ ) {
			kshell_meta[kj].ref_phi1 = itG->phi1;
			kshell_meta[kj].ref_PHI = itG->Phi;
			kshell_meta[kj].ref_phi2 = itG->phi2;
			kshell_meta[kj].ref_see = itG->see;
		}

		this->omp_knn_kshell_perimeter( nbIDs, nbShell, kshell_meta, kkmax ); //kshell_perimeter, kkmax );

		//##MK::when OMP parallelized this would require the storage of the thrad local data, otherwise concurrency for shared resource logknn!

		//intelligent output, do not output grains without any numerical resolvable neighbors, these are those shortly before becoming deleted!
		bool suspicious = false;
		for ( unsigned int kk = 1; kk <= kkmax; kk++ ) {
			if ( (int) ((nbShell[kk].last - nbShell[kk].first + 1) - (nbShell[kk-1].last - nbShell[kk-1].first + 1)) <= (int) 0 ) { suspicious = true; }
//if ( itG->extGID == 1399812 ) cout << "Grai1399812 " << (int) ((nbShell[kk].last - nbShell[kk].first + 1) - (nbShell[kk-1].last - nbShell[kk-1].first + 1)) << "\t\t" << (int) suspicious << std::endl;
		}
		bool NoNeighbors = false;
		if (kkmax >= 1 && ( (int) (nbShell[kkmax].last - nbShell[kkmax].first + 1) == (int) 0 ) ) { NoNeighbors = true; }

		if ( NoNeighbors != true ) {
			outputfile << itG->extGID << ";" << setprecision(8) << itG->x << ";" << setprecision(8) << itG->y << ";";
			for ( unsigned int kj = 0; kj <= kkmax; kj++ )	outputfile << (int) (nbShell[kj].last - nbShell[kj].first + 1) << ";";
			for ( unsigned int kj = 0; kj <= kkmax; kj++ )	outputfile << setprecision(8) << kshell_meta[kj].area_ksgrains << ";";
			for ( unsigned int kj = 0; kj < kkmax; kj++ )	outputfile << setprecision(8) << kshell_meta[kj].ksbnd_peri << ";"; //MK::mind only to kj < kkmax because regionboundaries not beyond kk-th shell!
			for ( unsigned int kj = 0; kj < kkmax; kj++ )	outputfile << setprecision(8) << kshell_meta[kj].area_hagbthresholded << ";";
			for ( unsigned int kj = 0; kj < kkmax; kj++ )	outputfile << setprecision(8) << kshell_meta[kj].area_velocity_w << ";";
			//for ( unsigned int kj = 0; kj < kkmax; kj++ )	outputfile << setprecision(8) << kshell_meta[kj].area_disori2ref_w << ";";
			//for ( unsigned int kj = 0; kj < kkmax; kj++ )	outputfile << setprecision(8) << kshell_meta[kj].area_dsee2ref_w << ";";
			outputfile << endl;
		}

/*cout << "\t\tKShellPerimeter\t\t";
for (unsigned int j = 0; j < kkmax; j++ ) { cout << kshell_meta[j].ksbnd_peri << "---" << kshell_meta[j].area_disori2ref_w << "---" << kshell_meta[j].area_dsee2ref_w << std::endl; } //kshell_perimeter[j] << ";" << std::endl;}
cout << std::endl;*/

		//delete [] kshell_perimeter;
		delete [] kshell_meta;

		//MK::DEBUG ONLY
		//if ( itG->extGID == 3466945 ) {
		//output only stats for those who loose faces
			//cout << itG->extGID << ";" << itG->size << endl;

			cout << itG->extGID << ";" << itG->size << ";" << setprecision(8) << itG->x << ";" << setprecision(8) << itG->y << "\t\t";
			for ( unsigned int kk = 0; kk <= kkmax; kk++ )	cout << (nbShell[kk].last - nbShell[kk].first + 1) << ";";
			for ( unsigned int kk = 0; kk <= kkmax; kk++ ) cout << nbShell[kk].first << "__" << nbShell[kk].last << ";";
			cout << endl;

			//intelligence, if in any case a reduction of number with increasing size report
			if ( NoNeighbors == true ) {
				std::cerr << "NO-NEIGHBORS;" << itG->extGID << ";" << setprecision(8) << itG->size << ";" << setprecision(8) << itG->x << ";" << setprecision(8) << itG->y << ";" << endl;
			}

			if ( suspicious == true ) {
				Statistics[KNN_STATS_SUSPICIOUS] += 1;
				if ( NoNeighbors != true ) {
					//std::cerr << itG->extGID << ";" << itG->size << ";" << itG->x << ";" << itG->y << endl;
					std::cerr << itG->extGID << ";" << setprecision(8) << itG->size << ";" << setprecision(8) << itG->x << ";" << setprecision(8) << itG->y << ";";
					for ( unsigned int kk = 0; kk <= kkmax; kk++ ) std::cerr << nbShell[kk].first << "__" << nbShell[kk].last << ";";
					for ( unsigned int kk = 0; kk <= kkmax; kk++ ) std::cerr << (nbShell[kk].last - nbShell[kk].first + 1) << ";";
					/*std::cerr << endl;
					for ( unsigned int kk = 0; kk <= kkmax; kk++ ) {
						for ( unsigned int nn = nbShell[kk].first; nn <= nbShell[kk].last; nn++ ) {
							std::cerr << nbIDs[nn] << ";" << kk << ";" << setprecision(8) << this->omp_getArea( nbIDs[nn] ) << ";" << setprecision(8) << this->omp_getBaryCenterX( nbIDs[nn]) << ";" << setprecision(8) << this->omp_getBaryCenterY( nbIDs[nn]) << endl;
						}
					}*/
					std::cerr << endl;
				}
			}
		//}

		nbIDs.clear();
		if ( nbShell != NULL) { delete [] nbShell; nbShell = NULL; }
	}

	delete [] TouchingDomainBnd; TouchingDomainBnd = NULL;

cout << "kNN with eliminating boundaries completed Total/Accepted/Expelled/Suspicious(ButKept) = " << (unsigned int) Statistics[KNN_STATS_TOTAL] << ";" << (unsigned int) (Statistics[KNN_STATS_TOTAL] - Statistics[KNN_STATS_EXPELLED]) << ";" << Statistics[KNN_STATS_EXPELLED] << ";" << Statistics[KNN_STATS_SUSPICIOUS] << endl;

	return stillgood;
}


inline double MemRegion::disori2mobility( double theta )
{
	double m = 1.0 - (0.99 * exp( -5.0 * (pow(	theta / DEG2RAD(15.0) , 9.0))));
	return m;
}


void MemRegion::omp_knn_kshell_perimeter( std::vector<unsigned int>& kids, ksmetaP meta, kshell_props* res, unsigned int nk )
{
	//this function calculates properties along those boundary faces of a target grain which enclose its k-th shell of neighbors
	//e.g. the 0-th kshell perimeter is the perimeter/grain boundary surface area of the target itself
	//the 1-th kshell perimeter is the accumulated boundary length between the 1-order neighbors to the target and their neighbors in the 2-th kshell
	//so in general the k-th kshell perimeter is the accumulated boundary length between the k-order neighbors to a target and their neighbors in the k+1-th kshell of the target
	//hence nk identifies how many such perimeter to calculate starting at and including the 0-th perimeter
	//MK::assumes that boundary grains have been filtered out, and that we do not track the zero grain!

	for ( unsigned int i = 0; i < nk; i++ ) { //MK:: < nk necessary because accessing elements[nk]
		//double wperi_see = 0.0;
		//double wperi_dis = 0.0;
		double peri = 0.0;
		double perithresholded = 0.0; //perimeter with specific properties, ##MK::at the moment a classical engineer's threshold to scan the HAGB-fraction of the grain
		double wperi_mp = 0.0;
		double ksarea = 0.0; //area of all the grains in this shell

		for ( unsigned int kgr = meta[i].first; kgr <= meta[i].last; kgr++ ) { //get boundary of k-th shell grain to its k+1-th neighbors
			unsigned int kgr_wmr = NOT_ASSIGNED_YET;
			unsigned int kgr_wid = NOT_ASSIGNED_YET;

			this->find_gid( kids[kgr], &kgr_wmr, &kgr_wid );
			if ( kgr_wmr == NOT_ASSIGNED_YET || kgr_wid == NOT_ASSIGNED_YET ) { cerr << "ERROR::LookupTableMissingEntry gid/whichmr/whichidx " << kids[kgr] << ";" << kgr_wmr << ";" << kgr_wid << std::endl; break; }

			//get list of neighbors of this grain
			/*double kgr_euler[3]; //##MK::utilize the quaternion directly to safe euler2quaternion conversions...
			kgr_euler[0] = this->mys->mrg[kgr_wmr]->GrainBucket[kgr_wid].phi1;
			kgr_euler[1] = this->mys->mrg[kgr_wmr]->GrainBucket[kgr_wid].Phi;
			kgr_euler[2] = this->mys->mrg[kgr_wmr]->GrainBucket[kgr_wid].phi2;*/

			ksarea += this->mys->mrg[kgr_wmr]->GrainBucket[kgr_wid].size;
			Face* KGrNeighbors = this->mys->mrg[kgr_wmr]->GrainBucket[kgr_wid].neighbors;
			unsigned int KGrHowManyOfThese = this->mys->mrg[kgr_wmr]->GrainBucket[kgr_wid].nfaces;

			for ( unsigned int kgrn = 0; kgrn < KGrHowManyOfThese; kgrn++ ) {
				//check which of the k-th grains KGrNeighbors[kgrn] is included in the following id interval [ kids[meta[i+1].first], kids[meta[i+1].last] ] of k+1-th shell ids
				unsigned int kplus1gr = NOT_ASSIGNED_YET;
				for ( unsigned int c = meta[i+1].first; c <= meta[i+1].last; c++ ) {
					if ( kids[c] != KGrNeighbors[kgrn].extNBID ) continue;
					//found, so store and can go out because ids are unique
					kplus1gr = KGrNeighbors[kgrn].extNBID;
					break;
				}
				if ( kplus1gr == NOT_ASSIGNED_YET ) { //verified that kplus1gr is no grain in the i+1-th kshell of the target
					continue;
				}

				//verified that kplus1gr is not only in the i+1-th kshell to the target but also a neighbor to KGrNeighbors[kgrn]
				//and -- it was found in the database, hence it forms a boundary segment contributing to the perimeter of the i-th kshell
				double dperi = KGrNeighbors[kgrn].size;
				peri = peri + dperi;

				//see difference nbor - ref
				//wperi_see += ( this->mys->mrg[KGrNeighbors[kgrn].mrg_idx]->GrainBucket[KGrNeighbors[kgrn].idx].see - res[i].ref_see ) * KGrNeighbors[kgrn].size;

				//disorientation angle along perimeter
				double kgrn_euler[3];
				kgrn_euler[0] = this->mys->mrg[KGrNeighbors[kgrn].mrg_idx]->GrainBucket[KGrNeighbors[kgrn].idx].phi1;
				kgrn_euler[1] = this->mys->mrg[KGrNeighbors[kgrn].mrg_idx]->GrainBucket[KGrNeighbors[kgrn].idx].Phi;
				kgrn_euler[2] = this->mys->mrg[KGrNeighbors[kgrn].mrg_idx]->GrainBucket[KGrNeighbors[kgrn].idx].phi2;

				//disorientation between k and k+1 grain
				//double disori = this->mys->mytopohdl->misorientationCubic( kgr_euler[0], kgr_euler[1], kgr_euler[2], kgrn_euler[0], kgrn_euler[1], kgrn_euler[2] );
				//disorientation between ref and k+1 grain
				double disori = this->mys->mytopohdl->misorientationCubic( res[i].ref_phi1, res[i].ref_PHI, res[i].ref_phi2,   kgrn_euler[0], kgrn_euler[1], kgrn_euler[2] );
				//wperi_dis += disori * KGrNeighbors[kgrn].size;

				//m*p with m \in [0, 1]
				wperi_mp = wperi_mp + ( this->disori2mobility(disori) * ( this->mys->mrg[KGrNeighbors[kgrn].mrg_idx]->GrainBucket[KGrNeighbors[kgrn].idx].see - res[i].ref_see ) * dperi );

				if ( disori >= Settings::HAGBDetectionThreshold ) 
					perithresholded = perithresholded + dperi;

				//cout << "\t\t\t" << RAD2DEG((double) disori / PI) << std::endl;
			}

		} //inspect how the next member of the k-shell contributes to peri

		res[i].ksbnd_peri = peri;
		if ( peri > DBL_EPSILON ) {
			res[i].area_hagbthresholded = perithresholded / peri;
			res[i].area_velocity_w = wperi_mp / peri;
			//res[i].area_disori2ref_w = wperi_dis / peri;
			//res[i].area_dsee2ref_w = wperi_see / peri;
		}
		res[i].area_ksgrains = ksarea;
	} //process the next shell


	//analyze the highest-order kshell for grain-based data separately, as the region bnd between the nk and nk+1 shell cannot be computed as only nk-order shells are tracked
	unsigned int i = nk;
	double ksarea = 0.0; //area of all the grains in this shell

	for ( unsigned int kgr = meta[i].first; kgr <= meta[i].last; kgr++ ) { //get boundary of k-th shell grain to its k+1-th neighbors
		unsigned int kgr_wmr = NOT_ASSIGNED_YET;
		unsigned int kgr_wid = NOT_ASSIGNED_YET;

		this->find_gid( kids[kgr], &kgr_wmr, &kgr_wid );
		if ( kgr_wmr == NOT_ASSIGNED_YET || kgr_wid == NOT_ASSIGNED_YET ) { cerr << "ERROR::LookupTableMissingEntry gid/whichmr/whichidx " << kids[kgr] << ";" << kgr_wmr << ";" << kgr_wid << std::endl; break; }

		ksarea += this->mys->mrg[kgr_wmr]->GrainBucket[kgr_wid].size;
	}
	res[i].area_ksgrains = ksarea;
}


double MemRegion::omp_getArea(unsigned int gid )
{
	//##MK::dummy in the local dataset find the grain if found get its barycenter x coordinate
	double res = -1.0;

	unsigned int luid = gid / Settings::LookupMaxGrainIDRange;

	//find grain with the local hashtable
	unsigned int lu_size = this->mys->LookupTable[luid].size;
	GrainHash* thebucket = this->mys->LookupTable[luid].bucket;
	unsigned int whichmr = NOT_ASSIGNED_YET;
	unsigned int whichidx = NOT_ASSIGNED_YET;

	for ( unsigned int g = 0; g < lu_size; g++ ) { //scan most interesting candidates
		if ( thebucket[g].gid == gid ) {
			whichmr = thebucket[g].mrg_idx;
			whichidx = thebucket[g].idx;
			break;
		}
	}

	if ( whichmr == NOT_ASSIGNED_YET || whichidx == NOT_ASSIGNED_YET ) {
//cerr << "ERROR detected during determination of number of faces per grain " << gid << endl;
		return res;
	}

	res = this->mys->mrg[whichmr]->GrainBucket[whichidx].size;
	return res;
}



double MemRegion::omp_getBaryCenterX(unsigned int gid )
{
	//##MK::dummy in the local dataset find the grain if found get its barycenter x coordinate
	double res = -1.0;

	unsigned int luid = gid / Settings::LookupMaxGrainIDRange;

	//find grain with the local hashtable
	unsigned int lu_size = this->mys->LookupTable[luid].size;
	GrainHash* thebucket = this->mys->LookupTable[luid].bucket;
	unsigned int whichmr = NOT_ASSIGNED_YET;
	unsigned int whichidx = NOT_ASSIGNED_YET;

	for ( unsigned int g = 0; g < lu_size; g++ ) { //scan most interesting candidates
		if ( thebucket[g].gid == gid ) {
			whichmr = thebucket[g].mrg_idx;
			whichidx = thebucket[g].idx;
			break;
		}
	}

	if ( whichmr == NOT_ASSIGNED_YET || whichidx == NOT_ASSIGNED_YET ) {
//cerr << "ERROR detected during determination of number of faces per grain " << gid << endl;
		return res;
	}

	res = this->mys->mrg[whichmr]->GrainBucket[whichidx].x;
	return res;
}


double MemRegion::omp_getBaryCenterY(unsigned int gid )
{
	//##MK::dummy in the local dataset find the grain if found get its barycenter x coordinate
	double res = -1.0;

	unsigned int luid = gid / Settings::LookupMaxGrainIDRange;

	//find grain with the local hashtable
	unsigned int lu_size = this->mys->LookupTable[luid].size;
	GrainHash* thebucket = this->mys->LookupTable[luid].bucket;
	unsigned int whichmr = NOT_ASSIGNED_YET;
	unsigned int whichidx = NOT_ASSIGNED_YET;

	for ( unsigned int g = 0; g < lu_size; g++ ) { //scan most interesting candidates
		if ( thebucket[g].gid == gid ) {
			whichmr = thebucket[g].mrg_idx;
			whichidx = thebucket[g].idx;
			break;
		}
	}

	if ( whichmr == NOT_ASSIGNED_YET || whichidx == NOT_ASSIGNED_YET ) {
//cerr << "ERROR detected during determination of number of faces per grain " << gid << endl;
		return res;
	}

	res = this->mys->mrg[whichmr]->GrainBucket[whichidx].y;
	return res;
}


void MemRegion::memhello( void )
{
	cout << "The memregion " << this->intID << " from the " << this->mys->mytopohdl->get_Rank() << "-th worker says hello!, my address is = " << this << endl;
}


void MemRegion::dummyCheckGrainBucketConsistency( void )
{
	unsigned int dn = 0;
	cout << "\nChecking consistency of face detection algorithm!" << endl;
	std::vector<Grain>::iterator g;
	long acctex = 0;
	long accfac = 0;
	for ( g = this->GrainBucket.begin(); g != this->GrainBucket.end(); g++ ) {
		acctex = acctex + g->nfaces;
		accfac = accfac + g->nfaces_identified;

		if ( g->nfaces != g->nfaces_identified ) {
			dn = g->nfaces - g->nfaces_identified;
			cerr << "Worker " << this->mys->mytopohdl->get_Rank() << " finds for grain " << g->extGID << " bnd = " << g->boundary << " an inconsistency in the number of faces nfaces/identified (" << g->nfaces << "---" << g->nfaces_identified << ") delta = " << dn << endl;
		}
	}

	cout << "ACCTEX = " << acctex << endl;
	cout << "ACCFAC = " << accfac << endl;
	cout << "NGR = " << this->GrainBucket.size()-1 << endl; //no zero grain
}


bool MemRegion::checkFaceReassignmentConsistency( void )
{
	std::vector<Grain>::iterator g;
	for ( g = this->GrainBucket.begin(); g != this->GrainBucket.end(); g++ ) {

//##MK::dummyif ( g->extGID == 6454463 ) {  cout << "Grain 6454463 has nfaces/nfaces_identified = " << g->nfaces << ";" << g->nfaces << "\t\t"; for ( unsigned int nb = 0; nb < g->nfaces; nb++ ) cout << g->neighbors[nb].extNBID << ";"; cout << endl; }

		if ( g->nfaces == g->nfaces_identified ) 
			continue;

		//not continued means inconsistency
cout << "ERROR::Worker " << this->mys->mytopohdl->get_Rank() << " CheckFaceReassgn\t\t" << g->nfaces << "\t\t" << g->nfaces_identified << "\t\t" << g->extGID << endl;
		return false;
	}
	return true;
}


Snapshot::Snapshot( void )
{
	meta.extID = NOT_YET_KNOWN;

	meta.nMemRegionsXYZ = NOT_ASSIGNED_YET;
	meta.nMemRegionsX = NOT_ASSIGNED_YET;
	meta.nMemRegionsY = NOT_ASSIGNED_YET;
	meta.nMemRegionsZ = NOT_ASSIGNED_YET;

	meta.ng = 0;
	meta.nf = 0;

	meta.RID = NOT_ASSIGNED_YET;
	meta.size_ng = 0.0;
	meta.size_nf = 0.0;
	
	mytopohdl = NULL; //MK::must not be deleted because is only back-referencing the owner!
}


Snapshot::~Snapshot( void )
{
	//mrg and LookupTable have own clean-up routines
	for ( unsigned int m = 0; m < mrg.size(); m++ ) {
		if (mrg.at(m) != NULL) {
			delete mrg.at(m);
			mrg.at(m) = NULL;
		}
	}
	
	//vector<int*> STL does not know that the pointed to object needs deallocation
	for ( unsigned int lu = 0; lu < LookupTable.size(); lu++ ) {
		if (LookupTable.at(lu).bucket != NULL) {
			delete[] LookupTable.at(lu).bucket;
			LookupTable.at(lu).bucket = NULL;
		}
	}
	//mytopohdl = NULL; //do not delete only back-reference!
}


bool Snapshot::readFacesForEachGrain( void )
{
	//a consistent amount of memory to store data for faces for each grain has already been set up in 
	//MemRegion::GrainBucket[].face
//##OMP pragma omp parallel for
	//MK::idea is simple all threads scan the array check if for each face they host one of the two grain IDS gA or gB and if so add for that grain the respective face as the a neighbor in the preallocated space
	MemRegion* am = NULL;
	bool mysuccess = true;
	for ( unsigned int tid = 0; tid < OMP_GET_NUM_THREADS; tid++ ) {
		if ( tid == OMP_GET_THREAD_NUM ) {
			am = this->mrg[tid];
			mysuccess = am->omp_storeFaces( tid );
		}
	}
	//##MK::exception handling

	//end of parallel region

	if ( mysuccess == true )
		cout << "...Worker " << this->mytopohdl->get_Rank() << " wrote Faces_" << this->meta.extID << " data successfully into database" << endl;
	else 
		cerr << "ERROR::Worker " << this->mytopohdl->get_Rank() << " read facial information for snapshot " << this->meta.extID << " was unsuccessful!" << endl;

	return mysuccess;
}


unsigned int Snapshot::baryxyz2mid( double bcx, double bcy, double bcz )
{
	//reads grain barycenter coordinates and determines in which memory region the barycenter --- i.e. the grain --- is located
	unsigned int mid = NOT_ASSIGNED_YET;
	MemRegion* am = NULL;

	bool found = false;
	for ( unsigned int mit = 0; mit < this->mrg.size(); ++mit ) {
		am = this->mrg[mit];
		if ( bcx < am->Geometry.xmi ) continue;
		if ( bcx > am->Geometry.xmx ) continue;
		if ( bcy < am->Geometry.ymi ) continue;
		if ( bcy > am->Geometry.ymx ) continue;
		if ( bcz < am->Geometry.zmi ) continue;
		if ( bcz > am->Geometry.zmx ) continue;

		//obviously in the mit-th MemRegion
		mid = mit;
		found = true;
		break; 
		//the breaking of testing a potential location in other regions also avoids the problem of numerical
		//unstable position exactly at the domain boundary!
	}
	//what happens if no region is found?
	if (found == false) {
		mid = FIRST_REGION;
	}

	return mid;
}


bool Snapshot::planDataPartitioning( void )
{
	//implements the mapping of the grains in response to their spatial position on the memory regions
	//MK::a zero grain, namely THE_DOMAIN_ITSELF must be initialized first!
	
	//grainIDs are stored in ID buckets to speed up finding of IDs during nearest neighbor computations
	unsigned int mrg_size = this->meta.nMemRegionsXYZ;
	unsigned int lu_size = ( Settings::LargestGrainID / Settings::LookupMaxGrainIDRange ) + 1; //MK::in case of downcasting

	//determine mapping of grain to memRegionID																							   //##MK::size of vectors may be problematic as it may exceed stack size given that .size() = sizeof(unsigned int)*ngrains per snapshot...
	vector<unsigned int>* gr2mr = NULL;
	try { gr2mr = new vector<unsigned int>; }
	catch (std::bad_alloc &exc) {
		cerr << "ERROR::Worker " << this->mytopohdl->get_Rank() << " allocation error for gr2mr" << endl; return false;
	}
	vector<unsigned int>* gr2luid = NULL;
	try { gr2luid = new vector<unsigned int>;  }
	catch (std::bad_alloc &exc) {
		cerr << "ERROR::Worker " << this->mytopohdl->get_Rank() << " allocation error for gr2luid" << endl; 
		if (gr2mr != NULL) delete gr2mr; gr2mr = NULL; return false;
	}

	unsigned int mtarget = NOT_ASSIGNED_YET;
	unsigned int luid = NOT_ASSIGNED_YET;

	//MK::initialize THE_DOMAIN_ITSELF, i.e. the zero grain, it lives by definition in the FIRST_REGION always
	mtarget = FIRST_REGION;
	luid = THE_DOMAIN_ITSELF / Settings::LookupMaxGrainIDRange;
	gr2mr->push_back( mtarget );
	gr2luid->push_back( luid );

	//##MK::identification process of grain location relative position in domain could be parallelized!
	unsigned int ngr = this->meta.ng;
	for ( unsigned int g = 0; g < ngr; ++g ) { //assign a grain a MemoryRegion
		mtarget = NOT_ASSIGNED_YET;
		mtarget = baryxyz2mid( this->mytopohdl->grbuf[g].x, this->mytopohdl->grbuf[g].y, this->mytopohdl->grbuf[g].z );

		luid = NOT_ASSIGNED_YET;
		luid = this->mytopohdl->grbuf[g].id / Settings::LookupMaxGrainIDRange;

//##DEBUGif ( this->mytopohdl->grbuf[g].id == 27027 ) cout << "DUMMY27027 bufid/mtarget/luid " <<  this->mytopohdl->grbuf[g].id << ";" << mtarget << ";" << luid << std::endl;

		//MK::employ not the best branch prediction for the benefit of code clarity and diagnostic capability...
		if (mtarget != NOT_ASSIGNED_YET && luid != NOT_ASSIGNED_YET && luid < lu_size) { //most likely case, a valid memory region was found
			gr2mr->push_back(mtarget); //MK::push_back possible because we walk in ascending gID order
			gr2luid->push_back(luid);
			continue;
		}

		//no valid was found, so its worthwhile to get some diagnostics
		if ( mtarget == NOT_ASSIGNED_YET ) { cerr << "ERROR::Worker " << this->mytopohdl->get_Rank() << " invalid MemRegion identified " << g << endl; }
		else if ( luid == NOT_ASSIGNED_YET ) { cerr << "ERROR::Worker " << this->mytopohdl->get_Rank() << " invalid, namely no Lookup ID identified " << g << endl;  }
		else if ( luid >= lu_size ) { cerr << "ERROR::Worker " << this->mytopohdl->get_Rank() << " invalid, namely too high Lookup ID identified " << g << endl; }
		else { cerr << "ERROR::Worker " << this->mytopohdl->get_Rank() << " invalid, unknown mistake " << g << endl; }
		if (gr2mr != NULL) { delete gr2mr; gr2mr = NULL; }
		if (gr2luid != NULL) { delete gr2luid; gr2luid = NULL; }
		return false;
	}

	//consistent?
	if ( gr2mr->size() != (ngr + ONE_FOR_THE_DOMAIN) || gr2luid->size() != (ngr + ONE_FOR_THE_DOMAIN) ) {
		cerr << "ERROR::Worker " << this->mytopohdl->get_Rank() << " memRegion cnts " << gr2mr->size() << " lookup cnts " << gr2luid->size() << " expected " << (ngr + ONE_FOR_THE_DOMAIN) << " inconsistent!" << endl;
		if (gr2mr != NULL) { delete gr2mr; gr2mr = NULL; }
		if (gr2luid != NULL) { delete gr2luid; gr2luid = NULL; }
		return false;
	}

	//accumulate counts for memory regions
	unsigned int* mcnt = NULL;
	try { mcnt = new unsigned int[mrg_size]; }
	catch (std::bad_alloc &exc) {
		cerr << "ERROR::Worker " << this->mytopohdl->get_Rank() << " unable to allocate counts per memory region!" << endl;
		if (gr2mr != NULL) { delete gr2mr; gr2mr = NULL; }
		if (gr2luid != NULL) { delete gr2luid; gr2luid = NULL; }
		return false;
	}
	for ( unsigned int mi = 0; mi < mrg_size; ++mi ) { mcnt[mi] = 0; }

	//accumulate counts for how many grains with gID in a particular numeral range \in \mathbb{N}
	unsigned int* lucnt = NULL;
	try { lucnt = new unsigned int[lu_size]; }
	catch (std::bad_alloc &exc) {
		cerr << "ERROR::Worker " << this->mytopohdl->get_Rank() << " unable to allocate counts per lookup table bucket!" << endl;
		if (gr2mr != NULL) { delete gr2mr; gr2mr = NULL; }
		if (gr2luid != NULL) { delete gr2luid; gr2luid = NULL; }
		if (mcnt != NULL) { delete[] mcnt; mcnt = NULL; }
		return false;
	}
	for (unsigned int li = 0; li < lu_size; ++li ) { lucnt[li] = 0; }

	//accumulate counts for the memregions and the gIDranges over all grains
	for (unsigned int mi = 0; mi < (ngr + ONE_FOR_THE_DOMAIN); mi++) { mcnt[gr2mr->at(mi)] += 1; }
	for (unsigned int li = 0; li < (ngr + ONE_FOR_THE_DOMAIN); li++) { lucnt[gr2luid->at(li)] += 1; }

	//##OMP pragma omp parallel
	MemRegion* am = NULL;
	bool mysuccess = true;
	for ( unsigned int tid = 0; tid < OMP_GET_NUM_THREADS; tid++ ) {
		if ( tid == OMP_GET_THREAD_NUM ) {
			am = this->mrg.at(tid);
			mysuccess = am->omp_allocGrainMemory( mcnt[tid] );
		}
	}
	//##exception handling, not as smart nMemRegions may have been set differently then omp_get_num_threads evaluates...

	//##OMP potentially lay out lookup table in thread 0 memory
	//##MK::consider packing in extra functions
	
	//create lookup table
	this->LookupTable.reserve( lu_size );

	for ( unsigned int lu = 0; lu < lu_size; lu++ ) {
		GrainHash* ahash = NULL;
		if (lucnt[lu] > 0) {
			try { ahash = new GrainHash[lucnt[lu]]; }
			catch (std::bad_alloc &exc) {
				cerr << "ERROR::Worker " << this->mytopohdl->get_Rank() << " error during memory allocation ahash for lu " << lu << endl;
				if (gr2mr != NULL) { delete gr2mr; gr2mr = NULL; }
				if (gr2luid != NULL) { delete gr2luid; gr2luid = NULL; }
				if (mcnt != NULL) { delete[] mcnt; mcnt = NULL; }
				if (lucnt != NULL) { delete[] lucnt; lucnt = NULL; }
				return false;
			}
		}
		this->LookupTable.push_back ( IDRangeBucket(0, lucnt[lu], ahash) );
	}
	//##MK::here definately end of parallel region

	//register THE_DOMAIN_ITSELF
	unsigned int entryid = UNKNOWN_ID;
	mtarget = gr2mr->at(THE_DOMAIN_ITSELF);
	luid = gr2luid->at(THE_DOMAIN_ITSELF);
	entryid = this->LookupTable[luid].next;
	this->LookupTable[luid].bucket[entryid].gid = THE_DOMAIN_ITSELF;
	this->LookupTable[luid].bucket[entryid].mrg_idx = mtarget;
	this->LookupTable[luid].bucket[entryid].idx = NOT_ASSIGNED_YET; //thread fills this in
	this->LookupTable[luid].next += 1; //##MK::this->LookupTable[luid].next + 1;

	if (this->meta.ng != this->mytopohdl->grbufsize) {
		cerr << "ERROR::Worker " << this->mytopohdl->get_Rank() << " MetaNG and GrBufSize dont match up!" << endl;
		if (gr2mr != NULL) { delete gr2mr; gr2mr = NULL; }
		if (gr2luid != NULL) { delete gr2luid; gr2luid = NULL; }
		if (mcnt != NULL) { delete[] mcnt; mcnt = NULL; }
		if (lucnt != NULL) { delete[] lucnt; lucnt = NULL; }
		return false;
	}

	//read out grain buffer and locate grains in buffer	
	for ( unsigned int g = 0; g < this->mytopohdl->grbufsize; g++ ) {
		mtarget = gr2mr->at(1+g);
		luid = gr2luid->at(1+g);
		//##MK::if planning to safe memory in the future do so for the luids because MemRegions are more costly to compute

		entryid = this->LookupTable[luid].next;
		//MK::here "next" means the next free slot in the LookupTable bucket for the IDRange luid
		if ( entryid >= this->LookupTable[luid].size ) { 
			cerr << "ERROR::Worker " << this->mytopohdl->get_Rank() << " attempting to fill more grains into LookupTable bucket than planned, there must be a logical mistake!" << endl; 
			if (gr2mr != NULL) { delete gr2mr; gr2mr = NULL; }
			if (gr2luid != NULL) { delete gr2luid; gr2luid = NULL; }
			if (mcnt != NULL) { delete[] mcnt; mcnt = NULL; }
			if (lucnt != NULL) { delete[] lucnt; lucnt = NULL; }
			return false;
		}

		this->LookupTable[luid].bucket[entryid].gid = this->mytopohdl->grbuf[g].id;
		this->LookupTable[luid].bucket[entryid].mrg_idx = mtarget; //the threads will check against this assignment
		this->LookupTable[luid].bucket[entryid].idx = NOT_ASSIGNED_YET; //the threads have to fill this in because only they know the idx --- the position in the thread-local MemRegion::GrainBucket
		this->LookupTable[luid].next = this->LookupTable[luid].next + 1;
	}

	//MK::Debug consistency checks, LookupTable[luid].next == LookupTable[luid].size for all luid
	for ( unsigned int li = 0; li < lu_size; li++ ) {
		if ( this->LookupTable[li].next != this->LookupTable[li].size ) {
			cerr << "ERROR::Worker " << this->mytopohdl->get_Rank() << " LookupTable_ " << this->meta.extID << " is inconsistent!" << endl; 
			if (gr2mr != NULL) { delete gr2mr; gr2mr = NULL; }
			if (gr2luid != NULL) { delete gr2luid; gr2luid = NULL; }
			if (mcnt != NULL) { delete[] mcnt; mcnt = NULL; }
			if (lucnt != NULL) { delete[] lucnt; lucnt = NULL; }
			return false;
		}
	}
	//finally =) I know how many grains in which MemRegion and in the Lookup table

	//MK::now the threads can read out the grains in parallel again for fill their properties into local database entries
	//##OMP parallel region
	mysuccess = true;
	for ( unsigned int tid = 0; tid < OMP_GET_NUM_THREADS; tid++ ) {
		if ( tid == OMP_GET_THREAD_NUM ) { 
			am = this->mrg[tid];
			mysuccess = am->omp_storeGrainData();
			if (mysuccess == false) {
				cerr << "ERROR::Worker " << this->mytopohdl->get_Rank() << " MemRegionHdl " << tid << " was unable to store the grains!" << endl;
				if (gr2mr != NULL) { delete gr2mr; gr2mr = NULL; }
				if (gr2luid != NULL) { delete gr2luid; gr2luid = NULL; }
				if (mcnt != NULL) { delete[] mcnt; mcnt = NULL; }
				if (lucnt != NULL) { delete[] lucnt; lucnt = NULL; }
				return false;
			}
		}
	}
	//##MK::though at the cost of a little bit cacheline wrestle for the LookupTable ...
	//##MK::exception handling
	if (gr2mr != NULL) { delete gr2mr; gr2mr = NULL; }
	if (gr2luid != NULL) { delete gr2luid; gr2luid = NULL; }
	if (mcnt != NULL) { delete [] mcnt; mcnt = NULL; }
	if (lucnt != NULL) { delete [] lucnt; lucnt = NULL; }

cout << "...Worker " << this->mytopohdl->get_Rank() << " wrote Texture_" << this->meta.extID << " successfully into database" << endl;

	return true;
}


bool Snapshot::omp_initOneMemRegion( unsigned int xx, unsigned yy, unsigned zz )
{
	//MK::CALLED FROM WITHIN OMP PRAGMA PARALLEL REFGION!
	//thread local allocation of Snapshot::MemRegion memory

	//generate a new memRegion, as function encapsulates, memory is ccNUMA local to thread, provided that memory allocator is smart...
	MemRegionP am = NULL;
	try {
		am = new MemRegion;
	}
	catch (std::bad_alloc &exc) {
		cerr << "ERROR::Worker " << this->mytopohdl->get_Rank() << " unable to allocate memory to init a memRegion!" << endl;
		return false;
	}

	am->mys = this; //back-reference to the calling snapshot class object that from now on is the owner of the memory region!
	am->intID = this->mcxyz2intID( xx, yy, zz );
	am->Geometry.xmi = ((double) xx ) * ( Settings::InitialDomainEdgeLength / ((double) Settings::MemRegionsX) ); //##MK::numerically range checking may be difficult...
	am->Geometry.xmx = ((double) (xx+1) ) * (Settings::InitialDomainEdgeLength / ((double) Settings::MemRegionsX) );
	am->Geometry.ymi = ((double) yy ) * (Settings::InitialDomainEdgeLength / ((double) Settings::MemRegionsY) );
	am->Geometry.ymx = ((double) (yy+1)) * (Settings::InitialDomainEdgeLength / ((double) Settings::MemRegionsY) );
	am->Geometry.zmi = ((double) zz ) * (Settings::InitialDomainEdgeLength / ((double) Settings::MemRegionsZ) );
	am->Geometry.zmx = ((double) (zz+1)) * (Settings::InitialDomainEdgeLength / ((double) Settings::MemRegionsZ) );

	//register the memory region in the snapshot class object
	this->mrg.push_back( NULL );
	this->mrg[this->mrg.size()-1] = am;

	return true;

//cout << "OMPMemRegionInitialized = " << setprecision(18) << this->mrg[this->mrg.size()-1]->Geometry.xmi << "--vs--" << am->Geometry.xmi << "\t\t" << setprecision(18) << this->mrg[this->mrg.size()-1]->Geometry.xmx << "--vs--" << am->Geometry.xmx << "\t\t" << setprecision(18) << am->Geometry.ymi << setprecision(18) << am->Geometry.zmi << endl;
}


inline unsigned int Snapshot::mcxyz2intID ( unsigned int mcx, unsigned int mcy, unsigned int mcz )
{
	//return implicit id of MemRegion in nx stacking in ny stacked in nz
	unsigned int id = mcx + ( mcy * this->meta.nMemRegionsX ) + ( mcz * this->meta.nMemRegionsX * this->meta.nMemRegionsY );
	return id;
}


bool Snapshot::initMemRegions( void )
{
	//building a 3d aggregate of cuboidal static memory regions that cover non-overlapping the unit domain in 3d
	unsigned int npx = this->meta.nMemRegionsX;
	unsigned int npy = this->meta.nMemRegionsY;
	unsigned int npz = this->meta.nMemRegionsZ;

	bool status = true;
	unsigned int mrid = FIRST_REGION; //memory region ID
	for ( unsigned int z = 0; z < npz; z++ ) { //coordinate system is x pointing to the right, y pointing in the plane, z pointing upwards
		for ( unsigned int y = 0; y < npy; y++ ) {
			for ( unsigned int x = 0; x < npx; x++ ) {
				mrid = this->mcxyz2intID( x, y, z );
				if ( OMP_GET_THREAD_NUM == mrid ) { //MK::each thread initializes only a single region
					status = this->omp_initOneMemRegion(x, y, z);
				}
			}
		}
	}
	
	return status;
}


bool Snapshot::fcbuf2database( void )
{
	bool filegood = true;
	
	filegood = readFacesForEachGrain();
	if ( filegood == false ) { cerr << "ERROR::Worker " << this->mytopohdl->get_Rank() << " was unable to readOut Grain Boundary Face information for " << this->meta.extID << endl; return false; }

	return filegood;
}


bool Snapshot::grbuf2database( void )
{
	//snapshot class implements the spatial partitioning of the dataset for machining off large dataset via threading
	bool filegood = true;
	
	filegood = initMemRegions();
	if ( filegood == false ) { 
		cerr << "ERROR::Worker " << this->mytopohdl->get_Rank() << " snapshot class object was unable to initialize MemRegions for " << this->meta.extID << endl; 
		return false;
	}
	//by now a NUMA-aware container for the spatial partitioning of the grain data are available

	filegood = planDataPartitioning();
	if ( filegood == false ) { 
		cerr << "ERROR::Worker " << this->mytopohdl->get_Rank() << " was unable to build lookup table and assign the grains to the MemRegions " << this->meta.extID << endl; return false; }
	//by now a lookup table is available to find where each grain is spatially located and to find grain IDs quickly

	return filegood;
}


void Snapshot::snhello( void )
{
	cout << "Hello from snapshot " << this->meta.extID << endl;
}


bool Snapshot::analyze_knn( unsigned int kkmax, std::ofstream &ofstr )
{
	//identification of k-nearest neighbors to all targets with their neighbors unless the target or its nbors has boundary contact
	unsigned int npx = this->meta.nMemRegionsX;
	unsigned int npy = this->meta.nMemRegionsY;
	unsigned int npz = this->meta.nMemRegionsZ;

//#pragma omp parallel
	bool mystatus = false;
	unsigned int m = 0;
	for ( unsigned int z = 0; z < npz; z++ ) {
		for ( unsigned int y = 0; y < npy; y++ ) {
			for ( unsigned int x = 0; x < npx; x++ ) {
				m = this->mcxyz2intID( x, y, z );
				if ( OMP_GET_THREAD_NUM == m ) { //evaluates only once to true in each thread
					MemRegion* am = this->mrg[m];
					mystatus = am->omp_knn_elimbnd( kkmax, ofstr );
				}
			}
		}
	}
	//##MK::exception handling
//pragma omp parallel

	return mystatus;
}


inline double Snapshot::calc_quantile( std::vector<double>* thecdf, double p )
{
	//Definition 3 of R. J. Hyndman and Y. Fan, Sample Quantiles in Statistical Packages, The American Statistician, Vol 50, 1996, 361-365
	//j = |_pn+m_| largest integer not greater than pn+m
	double n = thecdf->size(); //guaranteed to be >= 1
	double m = -1.0/2.0;
	double arg = p*n + m;
	long j = 0;
	if ( arg > QUANTILES_EPSILONLIMIT ) { //arg is positive so j at least 0
		j = arg;

		double jj = j;
		while ( jj > arg ) {
			j--;
			jj = j;
		}
		if ( j < 0 )
			return thecdf->at(0);
		
		//now j is largest integer not greater than arg
		if ( j > (thecdf->size() - 2) )
			return thecdf->at(thecdf->size()-1);

		double g = arg - jj;
		double gamma = 1.0; 
		if ( g > 0.0 && g <= QUANTILES_EPSILONLIMIT && (j % 2 == 0) ) //zero but not negative and j even
			gamma = 0.0;

		return ( (1.0 - gamma) * thecdf->at(j) + gamma * thecdf->at(j+1) );
	}
	//arg < EPSILONLIMIT
	return thecdf->at(0);
}


void Snapshot::analyze_vol_quants( MPI_VOLSTATS* meta, double* res )
{
	//###process through all my memory regions and their grains
	std::vector<double>* cdf = NULL;
	try { cdf = new std::vector<double>; }
	catch (std::bad_alloc &exc) {
		cerr << "ERROR::Worker::AnalyzeGrainSizeQuantiles " << this->mytopohdl->get_Rank() << " unable to allocate cdf!" << endl;
		return;
	}
	unsigned int ngr = 0;
	unsigned int ngr_elimbnd = 0;

	//MK::in the level-set approach the volume is normalized in the gridsize
	double normsize2physical = 1.0;
	if ( Settings::Dimensionality == TWODIMENSIONAL ) {
		normsize2physical = SQR( METER2MICRON(Settings::PhysicalDomainSize) );
	}
	if ( Settings::Dimensionality == THREEDIMENSIONAL ) {
		normsize2physical = CUBE( METER2MICRON(Settings::PhysicalDomainSize) );
	}

	for ( unsigned int mr = 0; mr < this->mrg.size(); mr++ ) { //collect target grains
		for ( unsigned int g = 0; g < this->mrg[mr]->GrainBucket.size(); g++ ) {
			if (  this->mrg[mr]->GrainBucket[g].extGID != THE_DOMAIN_ITSELF ) {
				ngr++;
				if ( this->mrg[mr]->GrainBucket[g].boundary != BOUNDARY_CONTACT ) {
					ngr_elimbnd++;
					cdf->push_back( this->mrg[mr]->GrainBucket[g].size * normsize2physical );
				}
			}
		} //next memory region please
	}

	//file in any case an analysis result
	meta->ng = ngr; //total grains
	meta->ng_elimbnd = ngr_elimbnd; //grains in the filet
	meta->totalsize_ng_elimbnd = 0.0;
	meta->meansize_ng_elimbnd = 0.0;
	meta->varsize_ng_elimbnd = 0.0;

	if ( cdf->size() < 1 )  { //no grain found or all not considered because of boundary contact?, MK::
		cout << "WARNING::Worker " << this->mytopohdl->get_Rank() << " no population found at all!" << endl;
		delete cdf; cdf = NULL; return;
	}
	if (ngr_elimbnd < MINIMUM_SAMPLE_SIZE_FOR_VARIANCE) {
		cout << "WARNING::Worker " << this->mytopohdl->get_Rank() << " population is too small for achieving sufficient convergence of t to normal distribution!" << endl;
		delete cdf; cdf = NULL; return;
	}

	//most likely though grains exists so calculate descriptive statistics of these
	double sum = 0.0;
	for ( unsigned int cg = 0; cg < cdf->size(); ++cg ) 
		sum += cdf->at(cg);

	double mu = sum / ((double) ngr_elimbnd);
	double var = 0.0;
	for ( unsigned int cg = 0; cg < cdf->size(); ++cg )
		var = var + SQR(cdf->at(cg) - mu);

	var = var / ((double) ( ngr_elimbnd - 1));

	meta->totalsize_ng_elimbnd = sum; //population size sufficient for t Distribution significantly converged to normal distribution...
	meta->meansize_ng_elimbnd = mu;
	meta->varsize_ng_elimbnd = var;

	//quantiles
	std::sort( cdf->begin(), cdf->end(), SortDoubleAscending );

	for ( unsigned int q = 0; q < STATISTICS_VOL_NQUANTILES; ++q ) //fill first with dummy values
		res[q] = 0.0;

	//then compute quantiles
	for ( unsigned int q = 1; q < STATISTICS_VOL_NQUANTILES; q++ ) {
		double p = ((double) q) / ((double) STATISTICS_VOL_NQUANTILES);
		res[q] = calc_quantile( cdf, p );

//cout << "\t\t" << q << ";" << p << ";" << res[q] << endl;
	}

	delete cdf; cdf = NULL;
}


unsigned int Snapshot::gsd_binning(double normsize)
{
	unsigned int nbins = Settings::get_GSDBinCount();
	unsigned int binend = Settings::get_GSDBinEnd();

	double b = (normsize - Settings::GSDBinMin) / Settings::GSDBinWidth;

	if (b >= 0.0) { //so at least Settings::SeeBinMin, the most likely case
		unsigned int bi = std::ceil(b);
		if (bi < (nbins - 1)) { //most likely case
			return bi;
		}
		//implicit else, map to bin for grains above highest value, i.e. > Settings::SeeBinMax
		return (1 + binend);
	}
	//implicit else, map to below smallest bin value, i.e. <= Settings::SeeBinMin
	return 0;
}


void Snapshot::analyze_gsd_histcnt( MPI_VOLSTATS* meta, double* res ) 
{
	//###process through all my memory regions and their grains
	std::vector<double>* val = NULL;
	try { val = new std::vector<double>; }
	catch (std::bad_alloc &exc) {
		cerr << "ERROR::Worker::AnalyzeGrainSizeQuantiles " << this->mytopohdl->get_Rank() << " unable to allocate val!" << endl;
		return;
	}
	unsigned int ngr = 0;
	unsigned int ngr_elimbnd = 0;

	//MK::in the level-set approach the volume is normalized in the gridsize
	double normsize2physical = 1.0;
	if ( Settings::Dimensionality == TWODIMENSIONAL ) {
		normsize2physical = SQR( METER2MICRON(Settings::PhysicalDomainSize) );
	}
	if ( Settings::Dimensionality == THREEDIMENSIONAL ) {
		normsize2physical = CUBE( METER2MICRON(Settings::PhysicalDomainSize) );
	}

	for ( unsigned int mr = 0; mr < this->mrg.size(); mr++ ) { //collect target grains
		for ( unsigned int g = 0; g < this->mrg[mr]->GrainBucket.size(); g++ ) {
			if (  this->mrg[mr]->GrainBucket[g].extGID != THE_DOMAIN_ITSELF ) {
				ngr++;
				if ( this->mrg[mr]->GrainBucket[g].boundary != BOUNDARY_CONTACT ) {
					ngr_elimbnd++;
					val->push_back( this->mrg[mr]->GrainBucket[g].size * normsize2physical );
				}
			}
		} //next memory region please
	}

	//file in any case an analysis result
	meta->ng = ngr; //total grains
	meta->ng_elimbnd = ngr_elimbnd; //grains in the filet
	meta->totalsize_ng_elimbnd = 0.0;
	meta->meansize_ng_elimbnd = 0.0;
	meta->varsize_ng_elimbnd = 0.0;

	if ( val->size() < 1 )  { //no grain found or all not considered because of boundary contact?, MK::
		cout << "WARNING::Worker " << this->mytopohdl->get_Rank() << " no population found at all!" << endl;
		if ( val != NULL ) { delete val; val = NULL; }
		return;
	}
	if (ngr_elimbnd < MINIMUM_SAMPLE_SIZE_FOR_VARIANCE) {
		cout << "WARNING::Worker " << this->mytopohdl->get_Rank() << " population is too small for achieving sufficient convergence of t to normal distribution!" << endl;
		if (val != NULL ) { delete val; val = NULL; }
		return;
	}

	//most likely though grains exists so calculate descriptive statistics of these
	double sum = 0.0;
	for ( unsigned int cg = 0; cg < val->size(); ++cg ) 
		sum += val->at(cg);

	double mu = sum / ((double) ngr_elimbnd);
	double _mu = 1.0 / mu;
	double var = 0.0;
	//variance computation, normalization, and binning
	for ( unsigned int cg = 0; cg < val->size(); ++cg ) {
		var = var + SQR(val->at(cg) - mu);
		double tmp = val->at(cg) * _mu;
		unsigned int bin = gsd_binning( tmp );
		res[bin] += 1.0;
	}

	var = var / ((double) ( ngr_elimbnd - 1));

	meta->totalsize_ng_elimbnd = sum; //population size sufficient for t Distribution significantly converged to normal distribution...
	meta->meansize_ng_elimbnd = mu;
	meta->varsize_ng_elimbnd = var;

	if ( val != NULL ) { delete val; val = NULL; }
}


topoHdl::topoHdl( void )
{
	myRank = MASTER;
	nRanks = SINGLE_PROCESS;

	grbufsize = 0;
	grbuf = NULL;

	fcbufsize = 0;
	fcbuf = NULL;

	bndbufsize = 0;
	bndbuf = NULL;

	helper_maxszgain_fw = NULL;
	helper_meandrvfsee_fw = NULL;
}


topoHdl::~topoHdl( void )
{
	for ( unsigned int db = 0; db < LocalDB.size(); ++db ) {
		if (LocalDB.at(db) != NULL) {
			delete LocalDB.at(db); 
			LocalDB.at(db) = NULL;
		}
	}

	//GrainIDWhiteList clears up after itself

	if ( grbuf != NULL ) {
		delete[] grbuf; grbuf = NULL;
		grbufsize = 0;
	}

	if (fcbuf != NULL) {
		delete[] fcbuf; fcbuf = NULL;
		fcbufsize = 0;
	}

	if (bndbuf != NULL) {
		delete[] bndbuf; bndbuf = NULL;
		bndbufsize = 0;
	}

	if (helper_maxszgain_fw != NULL) {
		delete[] helper_maxszgain_fw;
		helper_maxszgain_fw = NULL;
	}

	if (helper_meandrvfsee_fw != NULL) {
		delete[] helper_meandrvfsee_fw;
		helper_meandrvfsee_fw = NULL;
	}
}


void topoHdl::set_MPICommunication( int r, int nr ) {
	myRank = r;
	nRanks = nr;
}


int topoHdl::get_Rank( void ) {
	return myRank;
}


int topoHdl::get_nRanks( void ) {
	return nRanks;
}

void topoHdl::hellow( void )
{
	cout << "Hello world from process " << myRank << endl;
}


void topoHdl::init_MPIDatatypes( void ) 
{
	/*initializes user-defined MPI datatypes*/
	//MPI_Datatype MPI_VOLSTATS_Type
	MPI_Type_contiguous(2, MPI_UNSIGNED, &MPI_GrainBndContactIO_Type);

	//MPI_Datatype MPI_3DFaceIO_Type;
	int elementCounts3[2] = {1, 2};
	MPI_Aint displacements3[2] = {0, 1 * MPIIO_OFFSET_INCR_DOUBLE};
	MPI_Datatype oldTypes3[2] = {MPI_DOUBLE, MPI_UNSIGNED};
	MPI_Type_create_struct(2, elementCounts3, displacements3, oldTypes3, &MPI_3DFaceIO_Type);


	//MPI_Datatype MPI_2DGrainIO_Type
	int elementCounts2[2] = {9, 4 };
	MPI_Aint displacements2[2] = {0, 9 * MPIIO_OFFSET_INCR_DOUBLE};
	MPI_Datatype oldTypes2[2] = {MPI_DOUBLE, MPI_UNSIGNED};
	MPI_Type_create_struct(2, elementCounts2, displacements2, oldTypes2, &MPI_2DGrainIO_Type);

	//MPI_Datatype MPI_3DGrainIO_Type
	int elementCounts4[2] = {10, 4 };
	MPI_Aint displacements4[2] = {0, 10 * MPIIO_OFFSET_INCR_DOUBLE};
	MPI_Datatype oldTypes4[2] = {MPI_DOUBLE, MPI_UNSIGNED};
	MPI_Type_create_struct(2, elementCounts4, displacements4, oldTypes4, &MPI_3DGrainIO_Type);

	//MPI_Datatype MPI_SMDIO_Type
	int elementCounts5[2] = {2, 8 };
	MPI_Aint displacements5[2] = {0, 2 * MPIIO_OFFSET_INCR_DOUBLE};
	MPI_Datatype oldTypes5[2] = {MPI_DOUBLE, MPI_UNSIGNED};
	MPI_Type_create_struct(2, elementCounts5, displacements5, oldTypes5, &MPI_SMDIO_Type);

	//MPI_Datatype MPI_SEEAV_Type
	int elementCounts6[2] = {2, 2};
	MPI_Aint displacements6[2] = {0, 2 * MPIIO_OFFSET_INCR_DOUBLE};
	MPI_Datatype oldTypes6[2] = {MPI_DOUBLE, MPI_UNSIGNED};
	MPI_Type_create_struct(2, elementCounts6, displacements6, oldTypes6, &MPI_SEEAV_Type);

	//MPI_Datatype MPI_VOLSTATS_Type
	MPI_Type_contiguous( 5, MPI_DOUBLE, &MPI_VOLSTATS_Type );

	//MPI_Datatype MPI_RXSTATS_Type
	MPI_Type_contiguous( 5, MPI_DOUBLE, &MPI_RXSTATS_Type );

	//MPI_QUANTILES_Type;
	MPI_Type_contiguous( STATISTICS_VOL_NQUANTILES, MPI_DOUBLE, &MPI_QUANTILES_Type);

	//MPI_NUCBAIHIR_Type;
	//MPI_Type_contiguous( 1+4+2+2+2, MPI_DOUBLE, &MPI_NUCBAIHIR_Type);

	//MPI_AGGSTATS_Type
	MPI_Type_contiguous( (3+3+3+4), MPI_DOUBLE, &MPI_AGGSTATS_Type);

	MPI_Type_commit(&MPI_GrainBndContactIO_Type);
	MPI_Type_commit(&MPI_3DFaceIO_Type);
	MPI_Type_commit(&MPI_2DGrainIO_Type);
	MPI_Type_commit(&MPI_3DGrainIO_Type);
	MPI_Type_commit(&MPI_SMDIO_Type);
	MPI_Type_commit(&MPI_SEEAV_Type);
	MPI_Type_commit(&MPI_VOLSTATS_Type);
	MPI_Type_commit(&MPI_RXSTATS_Type);
	MPI_Type_commit(&MPI_QUANTILES_Type);
	//MPI_Type_commit(&MPI_NUCBAIHIR_Type);
	MPI_Type_commit(&MPI_AGGSTATS_Type);


	//cout << "\t\tMPI IO Types committed ..." << endl;
}


bool topoHdl::checkIfDatasetIsComplete( void )
{
	/*checks whether or not there are empty files or non-existent files*/
	bool accept = true;
	for ( vector<SnapshotMetaData>::iterator s = DatasetInfo.begin(); s != DatasetInfo.end(); ++s ) {
		if ( s->ng == NO_GRAINS || s->nf == NO_FACES ) {
			accept = false;
			break;
		}
	}

	if ( accept == true ) 
		cout << "...Worker " << myRank << " dataset complete" << endl;
	else 
		cerr << "ERROR::Worker " << myRank << " dataset incomplete or faulty!" << endl;

	return accept;
}


bool topoHdl::queryFileSizes( void )
{
	/*analyzes the file size of Texture_*.bin and Faces_*.bin*/
	string fn;
	double filesize[2] = {0.0, 0.0};
	double nobjects[2] = {0, 0};
	struct stat buf;

	cout << "...Worker " << myRank << " querying file sizes." << endl;
	for (unsigned int s = Settings::SnapshotFirst; s <= Settings::SnapshotLast; s = s + Settings::SnapshotOffset) { //check size of Texture_<s>.bin and Faces_<s>.bin by system call
		fn = "Texture_" + std::to_string(s) + ".bin";
		if (stat(fn.c_str(), &buf) != -1)
			filesize[TEXTUREFILE] = buf.st_size;
		else
			filesize[TEXTUREFILE] = 0.0;

		fn = "Faces_" + std::to_string(s) + ".bin";
		if (stat(fn.c_str(), &buf) != -1)
			filesize[FACEFILE] = buf.st_size;
		else
			filesize[FACEFILE] = 0.0;

		if (Settings::Dimensionality == THREEDIMENSIONAL) {
			nobjects[TEXTUREFILE] = filesize[TEXTUREFILE] / sizeof(MPI_3DGrainIO);
		}
		else if (Settings::Dimensionality == TWODIMENSIONAL) {
			nobjects[TEXTUREFILE] = filesize[TEXTUREFILE] / sizeof(MPI_2DGrainIO);
		}
		nobjects[FACEFILE] = filesize[FACEFILE] / sizeof(MPI_3DFaceIO);

		struct SnapshotMetaData smd = SnapshotMetaData(s, (Settings::MemRegionsX*Settings::MemRegionsY*Settings::MemRegionsZ),
			Settings::MemRegionsX, Settings::MemRegionsY, Settings::MemRegionsZ,
			nobjects[TEXTUREFILE], nobjects[FACEFILE], NOT_ASSIGNED_YET, filesize[TEXTUREFILE], filesize[FACEFILE]);

		DatasetInfo.push_back(smd);

		if (this->get_Rank() == MASTER) {
			cout << "\t\t" << "s;extSID;ng;nf;sng;snf = " << s << "---" << DatasetInfo[DatasetInfo.size() - 1].extID << ";" << DatasetInfo[DatasetInfo.size() - 1].ng << ";" << DatasetInfo[DatasetInfo.size() - 1].nf << ";" << setprecision(18) << DatasetInfo[DatasetInfo.size() - 1].size_ng << ";" << DatasetInfo[DatasetInfo.size() - 1].size_nf << endl;
		}
	}

	//##MK::list at the moment only sorted because of the aforementioned for loop!

	return checkIfDatasetIsComplete();
}


bool topoHdl::partitionWorkload(void)
{
	//needs only one process to execute, implements model to partition dataset as efficiently as possible on MPI processes available
	//MK::current model is to assign each process a certain block of contiguous snapshot (ids) and to fill up with raw data of size up to Settings::LocalDatabaseMaximumSize*/
	cout << "...Worker " << myRank << " partitioning workload" << endl;

	//first of all check whether entire dataset fits into the available processes at all
	vector<SnapshotMetaData>::iterator s;
	double MaxDatabaseVolume = nRanks * Settings::LocalDatabaseMaximumSize; //in Bytes
	double AggregatedDatasetVolume = 0.0;
	for (s = DatasetInfo.begin(); s != DatasetInfo.end(); ++s) {
		AggregatedDatasetVolume = AggregatedDatasetVolume + (s->size_ng + s->size_nf);
	}

	if (AggregatedDatasetVolume > (EMPIRICAL_OVERHEAD_FACTOR * MaxDatabaseVolume)) {
		cerr << "ERROR::Worker " << myRank << " dataset does not fit into aggregated memory of all processes" << endl;
		cerr << "ERROR::Worker " << myRank << " increase LocalDatabaseMaximumSize or utilize more MPI processes on more cores" << endl;
		return false;
	}

	//enough memory in total on all ranks should be available so attempting a partitioning is worthwhile...
	bool success = true;
	int rid = MASTER;
	double svol = 0.0;	//how much the pair Texture_/Faces_<fid>.bin will occupy for snapshot s
	double pvol = 0.0;	//how much the node already has...
	for (s = DatasetInfo.begin(); s != DatasetInfo.end(); ++s) {
		svol = s->size_ng + s->size_nf;

		if ((pvol + svol) <= Settings::LocalDatabaseMaximumSize) { //rank can still take this rawdata set pair
			s->RID = rid;
			pvol = pvol + svol;
			continue;
		}

		//else, need to choose the next rank in line, if any...
		rid++;
		pvol = 0.0;

		if (rid >= nRanks) {
			success = false; break;
		} //...fail, no further nodes available!

		//success, still a node available
		s->RID = rid;
		pvol = pvol + svol;
	}

	if (success == true) {
		writeWorkloadPartitioning();
		cout << "...Worker " << myRank << " the entire rawdata set was successfully partitioned on the computing nodes" << endl;
	}

	return success;
}


void topoHdl::writeWorkloadPartitioning(void)
{
	string logassgn_fn;
	logassgn_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".TrackingWorkloadDistro.csv";

	ofstream logassgn;
	logassgn.open(logassgn_fn.c_str(), std::ofstream::out | std::ofstream::trunc);
	if (logassgn.is_open() == true) {
		logassgn << "SnapshotID;SizeNG;SizeNF;NG;NF;AssignedToRank\n";
		logassgn << "Unit;Byte;Byte;1;1;MPIRankID\n";
		//region partitioning not outputted as it is execution manner dependent, in particular on OMP_NUM_THREADS!

		for (unsigned int s = 0; s < DatasetInfo.size(); s++) {
			logassgn << DatasetInfo.at(s).extID << ";" << setprecision(18) << DatasetInfo.at(s).size_ng << ";" << DatasetInfo.at(s).size_nf << ";" << DatasetInfo.at(s).ng << ";" << DatasetInfo.at(s).nf << ";" << DatasetInfo.at(s).RID << "\n";
		}
		logassgn.flush();
		logassgn.close();
	}
	else {
		cerr << "ERROR::Worker::WriteWorkloadPartitioning " << this->get_Rank() << " unable to write" << endl;
	}
}


void topoHdl::dummycontrol(void)
{
	//MK::quick and dirty dummy output
	cout << endl;
	for (int r = MASTER; r < nRanks; r++) {
		if (myRank == r) {
			vector<SnapshotMetaData>::iterator s;
			cout << endl << "Rank " << r << " DatasetInfo" << endl;
			for (s = DatasetInfo.begin(); s != DatasetInfo.end(); ++s) {
				cout << myRank << "\t\t" << s->extID << ";" << setprecision(18) << s->size_ng << ";" << setprecision(18) << s->size_nf << ";" << s->ng << ";" << s->nf << ";" << s->RID << endl;
			}
			cout << endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
}


void topoHdl::dummy_ascii_faces_binary(unsigned int extid)
{
	if (get_Rank() == MASTER) {
		string dummyfaces_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".CurrentFaces.Rank." + std::to_string(myRank)  + ".ID." + std::to_string(extid) + ".csv";
		ofstream dummyfaces;
		dummyfaces.open(dummyfaces_fn.c_str(), std::ofstream::out | std::ofstream::trunc);
		if (dummyfaces.is_open() == true ) {
			dummyfaces << "size;gA;gB\n";

			for (unsigned int f = 0; f < this->fcbufsize; f++) {
				dummyfaces << setprecision(8) << this->fcbuf[f].size << ";" << (unsigned int) this->fcbuf[f].gA << ";" << (unsigned int) this->fcbuf[f].gB << endl;
			}
			dummyfaces.flush();
			dummyfaces.close();
		}
		else {
			cerr << "ERROR::Worker::DummyASCIIFaces " << this->get_Rank() << " unable to write" << endl;
		}
	}
}


void topoHdl::dummy_ascii_grains_binary_3d(unsigned int extid)
{
	if (get_Rank() == MASTER) {
		string dummygrains_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".CurrentGrains.Rank." + std::to_string(myRank) + ".ID." + std::to_string(extid) + ".csv";
		ofstream dummygrains;
		dummygrains.open(dummygrains_fn.c_str(), std::ofstream::out | std::ofstream::trunc);
		if (dummygrains.is_open() == true) {
			dummygrains << "ID;NeighbourCount;Boundary;Vol;Surface;GBEnergy;BulkEnergy;phi1;Phi;phi2;x;y;z\n";

			for (unsigned int g = 0; g < this->grbufsize; g++) {
				//cout << "g/rho" << g << "\t\t" << this->grbuf[g].BulkEnergy << ";" << setprecision(8) << this->grbuf[g].BulkEnergy << endl;
				dummygrains << (unsigned int) this->grbuf[g].id << ";" << this->grbuf[g].NeighbourCount << ";" << this->grbuf[g].intersectsBoundaryGrain << ";" << setprecision(8) << this->grbuf[g].size << ";" << setprecision(8) << this->grbuf[g].surfaceArea << ";" << setprecision(8) << this->grbuf[g].GBEnergy << ";" << setprecision(8) << this->grbuf[g].BulkEnergy << ";" << setprecision(8) << this->grbuf[g].phi1 << ";" << setprecision(8) << this->grbuf[g].Phi << ";" << setprecision(8) << this->grbuf[g].phi2 << ";" << setprecision(8) << this->grbuf[g].x << ";" << setprecision(8) << this->grbuf[g].y << ";" << setprecision(8) << this->grbuf[g].z << endl;
			}
			dummygrains.flush();
			dummygrains.close();
		}
		else {
			cerr << "ERROR::Worker::DummyASCIIGrains " << this->get_Rank() << " unable to write" << endl;
		}
	}
}


bool topoHdl::determineWorkload(void)
{
	/*queries file size to plan the workload distribution*/
	if (myRank == MASTER) cout << "Identifying workload ..." << endl;

	if (queryFileSizes() == true) { //catch if files are non-existent or faulty!
		if (partitionWorkload() == true)
			return true;
	}
	return false;
}


void topoHdl::commitWorkload(void)
{
	//MK::assure to call by all processes! synchronizes SnapshotMetaData across them
	if (myRank == MASTER) {
		cout << "...Worker " << myRank << " committing workload" << "\n";
	}
	int nsnapshots[1] = { 0 };
	nsnapshots[0] = (int)DatasetInfo.size();		//##MK::usually all known grain coarsening dataset have substantially less than 2^32 - 1 timesteps!

	MPI_Bcast(&nsnapshots, 1, MPI_INT, MASTER, MPI_COMM_WORLD);

	//allocate bucket
	MPI_SMDIO* smdbuf = NULL;
	smdbuf = new MPI_SMDIO[nsnapshots[0]];

	if (myRank == MASTER) { //master fills bucket others wait
		for (int s = 0; s < nsnapshots[0]; s++) {
			smdbuf[s].ng = DatasetInfo[s].ng;
			smdbuf[s].nf = DatasetInfo[s].nf;
			smdbuf[s].size_ng = DatasetInfo[s].size_ng;
			smdbuf[s].size_nf = DatasetInfo[s].size_nf;
			smdbuf[s].extSID = DatasetInfo[s].extID;
			smdbuf[s].RID = DatasetInfo[s].RID;
			smdbuf[s].nmrx = DatasetInfo[s].nMemRegionsX;
			smdbuf[s].nmry = DatasetInfo[s].nMemRegionsY;
			smdbuf[s].nmrz = DatasetInfo[s].nMemRegionsZ;
			smdbuf[s].nmrxyz = DatasetInfo[s].nMemRegionsXYZ;
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(smdbuf, nsnapshots[0], MPI_SMDIO_Type, MASTER, MPI_COMM_WORLD); //MK::what happens when DatasetInfo.size() exceeds high value?

	//slaves read out
	if (myRank != MASTER) {
		for (int s = 0; s < nsnapshots[0]; s++) {
			DatasetInfo.push_back(SnapshotMetaData(smdbuf[s].extSID, smdbuf[s].nmrxyz, smdbuf[s].nmrx, smdbuf[s].nmry, smdbuf[s].nmrz,
				smdbuf[s].ng, smdbuf[s].nf, smdbuf[s].RID, smdbuf[s].size_ng, smdbuf[s].size_nf));
		}
	}

	delete[] smdbuf; smdbuf = NULL;

	cout << "...Worker " << myRank << " committed to his workload" << endl;
}



bool topoHdl::initMemoryForOneSnapshot(std::vector<SnapshotMetaData>::iterator sii, bool cleardb)
{
	//process generates private data container as snapshot class objects
	//processes know details from all data via DataInfo construct

	if (cleardb == true) { //MK::because we seek to initialize only one database object instance for TrackingSequentially
		for (unsigned int db = 0; db < LocalDB.size(); ++db) {
			if (LocalDB.at(db) != NULL) {
				delete LocalDB.at(db);
				LocalDB.at(db) = NULL;
			}
		}
		this->LocalDB.clear();
	}
	
	//generate single class object completely anew (cleardb==true) or (cleardb==false) keep accumulating class objects pointed to by thereafter by this->LocalDB
	SnapshotP asnap = NULL;
	try { asnap = new Snapshot; }
	catch (std::bad_alloc &exc) {
		cerr << "ERROR::Worker " << this->get_Rank() << " allocation error for memory of Snapshot class object!" << endl; return false;
	}

	asnap->mytopohdl = this;
	asnap->meta = SnapshotMetaData(sii->extID, sii->nMemRegionsXYZ, sii->nMemRegionsX, sii->nMemRegionsY, sii->nMemRegionsZ, sii->ng, sii->nf, this->myRank, sii->size_ng, sii->size_nf);

	//carry over pointer for the future
	this->LocalDB.push_back(asnap);
	
	return true;
}


bool topoHdl::nfaces_fcbuf2grbuf( std::vector<SnapshotMetaData>::iterator sii )
{
	//MK::this function assures that a) only faces between grains with id existent in the Texture_sii->extID.bin file exist and
	//b) resolves counting inconsistencies which in seldom cases cause --- because of numerics --- that some grains of grain pairs not to know each other being neighbors...

	//which grain IDs exist?
	unsigned int ngr = 1 + Settings::LargestGrainID;

	vector<bool>* gidExistsInTextureFile = NULL;
	gidExistsInTextureFile = new vector<bool>; //hash of existence
	for ( unsigned int b = 0; b < ngr; ++b ) //assume no grain exists
		gidExistsInTextureFile->push_back(false);

	for ( unsigned int g = 0; g < this->grbufsize; g++ ) //set explicitly IDs which do exist
		(*gidExistsInTextureFile)[this->grbuf[g].id] = true;
	
	//create a hashtable of face counts
	unsigned int* cnts = NULL;
	try { cnts = new unsigned int[ngr]; } //cnts is also an ID hash
	catch (std::bad_alloc &exc) {
		cerr << "ERROR::Worker " << this->get_Rank() << " allocation error to identify face counts in nfaces_fcbuf2grbuf" << endl;
		return false;
	}
	for ( unsigned int gid = THE_DOMAIN_ITSELF; gid <= Settings::LargestGrainID; ++gid ) { 
		cnts[gid] = 0;
	}

	//perform the counting by filling the hashtable and catch data inconsistencies
	for ( unsigned int f = 0; f < this->fcbufsize; f++ ) {
		//unexpected large IDs
		if ( fcbuf[f].gA > Settings::LargestGrainID || fcbuf[f].gB > Settings::LargestGrainID ) { //very unlikely case
			cerr << "ERROR::Worker " << this->get_Rank() << " data inconsistence for face " << f << " grainIDs larger as expected" << endl;
			if (gidExistsInTextureFile != NULL) { delete gidExistsInTextureFile; gidExistsInTextureFile = NULL; }
			if (cnts != NULL) { delete[] cnts; cnts = NULL; }
			return false;
		} //fools speculative execution for dummy purpose...

		//MK::do not exclude at that stage cases in which any of the two gA or gB equals THE_DOMAIN_ITSELF, as this only indicates
		//---at least in GraGLeS boundary contact for simulations without explicit periodic boundary conditions!

		//faces between grains with ids not included in the Texture file
		if ( (*gidExistsInTextureFile)[fcbuf[f].gA] == false ) {
			fcbuf[f].gA = THE_DOMAIN_ITSELF;
			fcbuf[f].gB = THE_DOMAIN_ITSELF;
		}
		if ( (*gidExistsInTextureFile)[fcbuf[f].gB] == false ) {
			fcbuf[f].gB = THE_DOMAIN_ITSELF;
			fcbuf[f].gA = THE_DOMAIN_ITSELF;
		}

		//account the face for a grain it it is not shared with the computational domain boundary (in case of open-boundary conditions)
		//MK::in case of periodic boundary conditions also in GraGLeS there should not be any face neighboring the "zero grain" <=> the computational domain ITSELF...
		if ( fcbuf[f].gA != THE_DOMAIN_ITSELF )
			cnts[fcbuf[f].gA] += 1;
		if ( fcbuf[f].gB != THE_DOMAIN_ITSELF )
			cnts[fcbuf[f].gB] += 1;
	}

	//update actual consistent number of faces shared with only the domain
	//MK::write back to grbuf by utilizing cnts as a hashtable at the expense of cache-misses for unsorted input data or ID fields with large gaps...
	unsigned int gid = NOT_ASSIGNED_YET;
	for ( unsigned int g = 0; g < this->grbufsize; ++g ) { //THE_DOMAIN_ITSELF has not account of its number of faces!
		gid = grbuf[g].id;
		
		if ( this->grbuf[g].NeighbourCount != cnts[gid] ) {
			this->grbuf[g].NeighbourCount = cnts[gid]; //reset neighbor count to detected number of faces
		}
	}

	//now the number of faces identified in grbuf is consistent with the fcbuf data
	//MK::this procedure is necessary as programs which consider grain boundary faces as local members may because of numerical inaccuracies 
	//not contain or have detected in every possible topological and geometrical configuration sufficiently resolvable faces 
	//and hence recognize all neighbors. In effect, a face between grain A and B may have been detected by grain A but not by B
	
	delete gidExistsInTextureFile; gidExistsInTextureFile = NULL;
	delete [] cnts; cnts = NULL;
	
	return true;
}


bool topoHdl::read_mpiio_faces_binary_2d3d( std::vector<SnapshotMetaData>::iterator sii )
{
	double rtimer = MPI_Wtime();

	//reject too large files that cannot read with MPI with the current explicit unsigned int formalism
	double tmp = sii->size_nf / (double) sizeof(MPI_3DFaceIO);
	if ((tmp - 2.0) >= (double)std::numeric_limits<unsigned int>::max()) {
		cerr << "ERROR::Worker " << this->get_Rank() << " Faces file " << sii->extID << " is too large to be read with the current data structure!" << endl;
		return false;
	}

	unsigned int TotalBlocksToRead = ( sii->size_nf / Settings::MPIReadBlockLength) + 1;
	unsigned int ElementsPerBlock = Settings::MPIReadBlockLength / sizeof(MPI_3DFaceIO); //MK::casts down, i.e. implicit round off
	unsigned int elementsTotal = sii->nf;
	unsigned int elementsRead = 0;
	unsigned int elementsNow = 0;

	//get a large enough buffer
	if (this->fcbuf != NULL) { delete[] this->fcbuf; }
	this->fcbuf = NULL;
	this->fcbufsize = 0;
	try { this->fcbuf = new MPI_3DFaceIO[elementsTotal]; }
	catch (std::bad_alloc &exc) { 
		cerr << "ERROR::Worker " << this->get_Rank() << " MPI_3DFaceIO buffer allocation error!" << endl; return false;
	}
	this->fcbufsize = elementsTotal;

	MPI_File ioReadFileHdl;
	MPI_Status ioReadFileStatus;

	string fn = "Faces_" + std::to_string(sii->extID) + ".bin";
	MPI_File_open(MPI_COMM_SELF, fn.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &ioReadFileHdl);
	MPI_File_seek(ioReadFileHdl, 0, MPI_SEEK_SET);

	//cout << "\t\tReading file in " << TotalBlocksToRead << " blocks with " << ElementsPerBlock << " elements each" << endl;
	for ( unsigned int b = 0; b < TotalBlocksToRead; b++ ) {
		elementsNow = ElementsPerBlock;
		if ( (elementsTotal - elementsRead) < elementsNow ) elementsNow = elementsTotal - elementsRead;
		if ( elementsNow > 0 ) {
			MPI_3DFaceIO* rbuf = NULL;
			try { rbuf = new MPI_3DFaceIO[elementsNow]; }
			catch (std::bad_alloc &exc) {
				cerr << "ERROR::Worker " << this->get_Rank() << " unable to allocate buffer in read_binary_faces!" << endl;
				MPI_File_close(&ioReadFileHdl); return false;
			}

			MPI_File_read( ioReadFileHdl, rbuf, elementsNow, MPI_3DFaceIO_Type, &ioReadFileStatus);
			//implicit advance of fp

			for ( unsigned int e = 0; e < elementsNow; ++e) { //explicit data transfer
				fcbuf[elementsRead+e].size = rbuf[e].size;
				fcbuf[elementsRead+e].gA = rbuf[e].gA;
				fcbuf[elementsRead+e].gB = rbuf[e].gB;
			}

			elementsRead = elementsRead + elementsNow;

//cout << "\t\tBlockID/elementsRead/elementsNow/elementsTotal--time = " << b << "/" << elementsRead << "/" << elementsNow << "/" << elementsTotal << "\t\t\t" << (MPI_Wtime() - rtimer) << "\t\tseconds" << endl;
			delete [] rbuf; rbuf = NULL;
		}
	} //data read
	MPI_File_close(&ioReadFileHdl);

	//##MK::add a back-reference to the hosting topoHdl

	cout << "...Worker " << this->get_Rank() << " read Faces_" << sii->extID << " in " << setprecision(8) << (MPI_Wtime() - rtimer) << " seconds" << endl;

	//MK::transcode the binary content into an ASCII
	if (Settings::DeveloperMode == true) {
		this->dummy_ascii_faces_binary(sii->extID);
	}

	return true;
}


void topoHdl::modelspecific_rawdata_modification( std::vector<SnapshotMetaData>::iterator siii )
{
	double rtimer = MPI_Wtime();
	//this function serves to implement user-defined rescaling operations on the raw data
	if ( Settings::SimModel == E_LEVELSET ) {
		//MK::GraGLeS outputs a scaled BulkEnergy
		double scaling = Settings::HAGBEnergy / (Settings::PhysicalDomainSize * Settings::DislocEnPerM);

		for (unsigned int g = 0; g < this->grbufsize; ++g ) {
			this->grbuf[g].BulkEnergy = this->grbuf[g].BulkEnergy * scaling;
		}
	}
	cout << "...Worker " << this->get_Rank() << " performed user-defined unit scaling operations on raw data " << siii->extID << " in " << setprecision(8) << (MPI_Wtime() - rtimer) << " seconds" << endl;
}


bool topoHdl::read_mpiio_grains_binary_2d( std::vector<SnapshotMetaData>::iterator sii )
{
	double rtimer = MPI_Wtime();

	//reject too large files that cannot be read with MPI with the current explicit unsigned int formalism
	double tmp = sii->size_ng / (double) sizeof(MPI_2DGrainIO);
	if ((tmp - 2.0) >= (double)std::numeric_limits<unsigned int>::max()) {
		cerr << "ERROR::Worker " << this->get_Rank() << " Texture file " << sii->extID << " is too large to be read with the current data structure!" << endl;
		return false;
	}

	//plan reading process...
	unsigned int TotalBlocksToRead = ( sii->size_ng / Settings::MPIReadBlockLength) + 1;
	unsigned int ElementsPerBlock = Settings::MPIReadBlockLength / sizeof(MPI_2DGrainIO); //##MK::we read objects of a 2D structure!
	unsigned int elementsTotal = sii->ng;
	unsigned int elementsRead = 0;
	unsigned int elementsNow = 0;

	//get large enough buffer
	//MK::we use a MPI_3DGrainIO bucket because internally the grains are always 3D!
	if (this->grbuf != NULL) {
		delete[] this->grbuf;
	}
	this->grbuf = NULL;
	this->grbufsize = 0;
	try { this->grbuf = new MPI_3DGrainIO[elementsTotal]; }
	catch (std::bad_alloc &exc) {
		cerr << "ERROR::Worker " << this->get_Rank() << " MPI_3DGrainIO allocation error in fucntion read_mpiio_grains_binary_2d!" << endl;
		return false;
	}
	this->grbufsize = elementsTotal;

	MPI_File ioReadFileHdl;
	MPI_Status ioReadFileStatus;

	string fn = "Texture_"  + std::to_string(sii->extID) + ".bin";
	MPI_File_open(MPI_COMM_SELF, fn.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &ioReadFileHdl);
	MPI_File_seek(ioReadFileHdl, 0, MPI_SEEK_SET);

	//cout << "\t\tReading file in " << TotalBlocksToRead << " blocks with " << ElementsPerBlock << " elements each" << endl;
	for ( unsigned int b = 0; b < TotalBlocksToRead; b++ ) {
		elementsNow = ElementsPerBlock;
		if ( (elementsTotal - elementsRead) < elementsNow ) elementsNow = elementsTotal - elementsRead;
		if ( elementsNow > 0 ) {
			MPI_2DGrainIO* rbuf = NULL;
			try { rbuf = new MPI_2DGrainIO[elementsNow]; }
			catch (std::bad_alloc &exc) {
				cerr << "ERROR::Worker " << this->get_Rank() << " MPI_2DGrainIO buffer allocation error!" << endl; 
				MPI_File_close(&ioReadFileHdl); return false;
			}

			MPI_File_read( ioReadFileHdl, rbuf, elementsNow, MPI_2DGrainIO_Type, &ioReadFileStatus);

			//explicit datatransfer
			for ( unsigned int e = 0; e < elementsNow; ++e ) {
				grbuf[elementsRead+e].size = rbuf[e].volume;
				grbuf[elementsRead+e].surfaceArea = rbuf[e].perimeter;
				grbuf[elementsRead+e].GBEnergy = rbuf[e].GBEnergy;
				grbuf[elementsRead+e].BulkEnergy = rbuf[e].BulkEnergy;
				grbuf[elementsRead+e].phi1 = rbuf[e].phi1;
				grbuf[elementsRead+e].Phi = rbuf[e].Phi;
				grbuf[elementsRead+e].phi2 = rbuf[e].phi2;
				grbuf[elementsRead+e].x = rbuf[e].x;
				grbuf[elementsRead+e].y = rbuf[e].y;
				grbuf[elementsRead+e].z = FLATTEN_THIRD_DIMENSION;
				grbuf[elementsRead+e].id = rbuf[e].id;
				grbuf[elementsRead+e].NeighbourCount = rbuf[e].NeighbourCount;
				grbuf[elementsRead+e].intersectsBoundaryGrain = rbuf[e].intersectsBoundaryGrain;
				grbuf[elementsRead+e].pad = rbuf[e].pad;
			}

			elementsRead = elementsRead + elementsNow;

//cout << "\t\tBlockID/elementsRead/elementsNow/elementsTotal--time = " << b << "/" << elementsRead << "/" << elementsNow << "/" << elementsTotal << "\t\t\t" << (MPI_Wtime() - rtimer) << "\t\tseconds" << endl;

			delete [] rbuf; rbuf = NULL;
		}
	} //data read
	MPI_File_close(&ioReadFileHdl);

	cout << "...Worker " << this->get_Rank() << " read 2D-Texture_" << sii->extID << " in " << setprecision(8) << (MPI_Wtime() - rtimer) << " seconds" << endl;

	this->modelspecific_rawdata_modification( sii );

	//transcode binary to ASCII
	if ( Settings::DeveloperMode == true ) 
		this->dummy_ascii_grains_binary_3d( sii->extID ); //ASCII format does not distinguish between 2d and 3d but if in doubt FLATTENS missing dimension

	return true;
}


bool topoHdl::read_mpiio_grains_binary_3d( std::vector<SnapshotMetaData>::iterator sii )
{
	double rtimer = MPI_Wtime();

	//reject too large files that cannot be read with MPI with the current explicit unsigned int formalism
	double tmp = sii->size_ng / (double) sizeof(MPI_3DGrainIO);
	if ((tmp - 2.0) >= (double)std::numeric_limits<unsigned int>::max()) {
		cerr << "ERROR::Worker " << this->get_Rank() << " Texture file " << sii->extID << " is too large to be read with the current data structure!" << endl;
		return false;
	}

	//plan reading process...
	unsigned int TotalBlocksToRead = ( sii->size_ng / Settings::MPIReadBlockLength) + 1;
	unsigned int ElementsPerBlock = Settings::MPIReadBlockLength / sizeof(MPI_3DGrainIO);
	unsigned int elementsTotal = sii->ng;
	unsigned int elementsRead = 0;
	unsigned int elementsNow = 0;

	//get large enough buffer
	if (this->grbuf != NULL) {
		delete[] this->grbuf;
	}
	this->grbuf = NULL;
	this->grbufsize = 0;
	try { this->grbuf = new MPI_3DGrainIO[elementsTotal]; }
	catch (std::bad_alloc &exc) {
		cerr << "ERROR::Worker " << this->get_Rank() << " MPI_3DGrainIO buffer allocation error!" << endl; 
		return false;
	}
	this->grbufsize = elementsTotal;

	MPI_File ioReadFileHdl;
	MPI_Status ioReadFileStatus;

	string fn = "Texture_" + std::to_string(sii->extID) + ".bin";
	MPI_File_open(MPI_COMM_SELF, fn.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &ioReadFileHdl);
	MPI_File_seek(ioReadFileHdl, 0, MPI_SEEK_SET);

	//cout << "\t\tReading file in " << TotalBlocksToRead << " blocks with " << ElementsPerBlock << " elements each" << endl;
	for ( unsigned int b = 0; b < TotalBlocksToRead; b++ ) {
		elementsNow = ElementsPerBlock;
		if ( (elementsTotal - elementsRead) < elementsNow ) elementsNow = elementsTotal - elementsRead;
		if ( elementsNow > 0 ) {
			MPI_3DGrainIO* rbuf = NULL;
			try { rbuf = new MPI_3DGrainIO[elementsNow]; }
			catch (std::bad_alloc &exc) {
				cerr << "ERROR::Worker " << this->get_Rank() << " unable to allocate buffer in read_binary_grains3d!" << endl;
				MPI_File_close(&ioReadFileHdl); return false;
			}

			MPI_File_read( ioReadFileHdl, rbuf, elementsNow, MPI_3DGrainIO_Type, &ioReadFileStatus);
						
			for ( unsigned int e = 0; e < elementsNow; ++e ) { //explicit datatransfer
				grbuf[elementsRead+e].size = rbuf[e].size;
				grbuf[elementsRead+e].surfaceArea = rbuf[e].surfaceArea;
				grbuf[elementsRead+e].GBEnergy = rbuf[e].GBEnergy;
				grbuf[elementsRead+e].BulkEnergy = rbuf[e].BulkEnergy;
				grbuf[elementsRead+e].phi1 = rbuf[e].phi1;
				grbuf[elementsRead+e].Phi = rbuf[e].Phi;
				grbuf[elementsRead+e].phi2 = rbuf[e].phi2;
				grbuf[elementsRead+e].x = rbuf[e].x;
				grbuf[elementsRead+e].y = rbuf[e].y;
				grbuf[elementsRead+e].z = rbuf[e].z;
				grbuf[elementsRead+e].id = rbuf[e].id;
				grbuf[elementsRead+e].NeighbourCount = rbuf[e].NeighbourCount;
				grbuf[elementsRead+e].intersectsBoundaryGrain = rbuf[e].intersectsBoundaryGrain;
				grbuf[elementsRead+e].pad = rbuf[e].pad;
			}

			elementsRead = elementsRead + elementsNow;

//cout << "\t\tBlockID/elementsRead/elementsNow/elementsTotal--time = " << b << "/" << elementsRead << "/" << elementsNow << "/" << elementsTotal << "\t\t\t" << (MPI_Wtime() - rtimer) << "\t\tseconds" << endl;

			delete [] rbuf; rbuf = NULL;
		}
	} //data read
	MPI_File_close(&ioReadFileHdl);

	cout << "Worker " << this->get_Rank() << " read 3D-Texture_" << sii->extID << " in " << setprecision(8) << (MPI_Wtime() - rtimer) << " seconds" << endl;

	this->modelspecific_rawdata_modification( sii );

	//transcode binary file to ASCII
	if ( Settings::DeveloperMode == true ) 
		this->dummy_ascii_grains_binary_3d( sii->extID );

	return true;
}


bool topoHdl::read_targets_fromlist( void ) {
	double timer = MPI_Wtime();

	GrainIDWhiteList.clear();

	ifstream trgfile;
	string trgline;
	istringstream line;
	string datapiece;

	trgfile.open( Settings::TargetGIDFromFilename );
	if ( trgfile.is_open() == true ) {
		//format of the file is only IDs, headerless, single column
		//single header
		getline( trgfile, trgline );
		
		while ( trgfile.good() == true ) {
			getline( trgfile, trgline );
			unsigned int s = trgline.length();

			if ( s > 0 ) { //only a target id 
				long tgid = std::stol( trgline ); //accept only gIDs on (THE_DOMAIN_ITSELF,Settings::LargestGrainID]
				if ( tgid > THE_DOMAIN_ITSELF && tgid <= Settings::LargestGrainID )
					GrainIDWhiteList.push_back( tgid );
				else {
					cerr << "ERROR::Attempting to feed an invalid grainID via TargetGrainIDs supplementary file!" << endl; trgfile.close(); return false;
				}
			}
			//ignore blank lines
		}
		trgfile.close();
		cout << "...Worker " << this->get_Rank() << " loaded " << GrainIDWhiteList.size() << " valid user-selected grainIDs to work on selectively in " << (double) (MPI_Wtime() - timer) << " seconds" << endl;
		return true;
	}
	return false;
}


bool topoHdl::readSnapshot( std::vector<SnapshotMetaData>::iterator si, bool ClearDatabase )
{
	double t1 = 0.0; //profiling
	double t2 = 0.0;
	double timer1 = MPI_Wtime();

	//process reads in an entire snapshot with externalID, when returning always true we are still good
	if ( initMemoryForOneSnapshot( si, ClearDatabase ) == false ) { 
		cerr << "ERROR::Worker " << this->get_Rank() << " detected a memory error during creation of SnapshotClass object for " << si->extID << endl; 
		return false;
	}

	t1 = MPI_Wtime();
	if ( this->read_mpiio_faces_binary_2d3d ( si ) == false ) {
		cerr << "ERROR::Worker " << this->get_Rank() << " detected file reading or memory error during snapshot instance attempting Faces_" << si->extID << " file access!" << endl; return false;
	}
	t2 = MPI_Wtime();
	ioprofiler.logev("MPIIOFaces" + std::to_string(si->extID), (t2 - t1));

	t1 = MPI_Wtime();
	if ( Settings::Dimensionality == THREEDIMENSIONAL ) {
		if ( this->read_mpiio_grains_binary_3d ( si ) == false ) { 
			cerr << "ERROR::Worker " << this->get_Rank() << " detected file reading or memory error during snapshot instance attempting Texture_" << si->extID << " file access!" << endl;
			return false;
		}
	} 
	else if ( Settings::Dimensionality == TWODIMENSIONAL ) {
		if ( this->read_mpiio_grains_binary_2d ( si ) == false ) { 
			cerr << "ERROR::Worker " << this->get_Rank() << " detected file reading or memory error during snapshot instance attempting Texture_" << si->extID << " file access!" << endl;
			return false;
		}
	}
	else {
		cerr << "ERROR::Worker " << this->get_Rank() << " invalid dimensionality during readSnapshot!" << endl;
		return false;
	}
	t2 = MPI_Wtime();
	ioprofiler.logev("MPIIOTexture" + std::to_string(si->extID), (t2 - t1));

	t1 = MPI_Wtime();
	//count faces per grain on fcbuf and update nfaces for all grains in grbuf
	if ( this->nfaces_fcbuf2grbuf( si ) == false ) {
		cerr << "ERROR::Worker " << this->get_Rank() << " counting of exact number of faces failed!" << endl;
		return false;
	}

	//MK::mind that the "adhoc" analyzer routines --- in contrast to the "db" equivalents demand that there is only one class object instance in the LocalDB!
	//register the grains
	if (ClearDatabase == true) {
		if (this->LocalDB.size() != 1) {
			cerr << "ERROR::Worker " << this->get_Rank() << " desiring one after another tracking but having more than one database object!" << endl; return false;
		}
		if (this->LocalDB.at(0)->grbuf2database() == false) {
			cerr << "ERROR::Worker " << this->get_Rank() << " the transfer of the grain-based data failed!" << endl;
			return false;
		}
	}
	else {
		//register grains, now only just as much memory to store faces such that nfaces == nfaces_identified
		if (this->LocalDB.at(this->LocalDB.size() - 1)->grbuf2database() == false) {
			cerr << "ERROR::Worker " << this->get_Rank() << " the transfer of the grain-based data failed!" << endl;
			return false;
		}
	}

	//either case we now no longer need the heavy data I/O buffer for Texture files
	if (this->grbuf != NULL) { delete[] this->grbuf; }
	this->grbuf = NULL; this->grbufsize = 0;
	
	//register faces
	if (ClearDatabase == true) {
		if (this->LocalDB.at(0)->fcbuf2database() == false) {
			cerr << "ERROR::Worker " << this->get_Rank() << " the transfer of the face-based data failed!" << endl;
			return false;
		}
	}
	else {
		if (this->LocalDB[this->LocalDB.size() - 1]->fcbuf2database() == false) {
			cerr << "ERROR::Worker " << this->get_Rank() << " the transfer of the face-based data failed!" << endl;
			return false;
		}
	}

	//either case we now no longer need the heavy data I/O buffer for Faces files
	if (this->fcbuf != NULL) { delete[] this->fcbuf; }
	this->fcbuf = NULL; this->fcbufsize = 0;

	t2 = MPI_Wtime();
	ioprofiler.logev("DBHandling" + std::to_string(si->extID), (t2 - t1));

	double timer2 = MPI_Wtime();
	ioprofiler.logev("ReadSnapshot" + std::to_string(si->extID), (timer2 - timer1));
	cout << "...Worker " << this->get_Rank() << " read Snapshot_" << si->extID << " successfully in " << (timer2-timer1) << " seconds" << endl;
	//Texture_<id> and Faces_<id> now successful read and registered in the database
	return true;
}


bool topoHdl::readDatasets( void )
{
	//constructs local database of Grains and Faces by reading sequential all my files 
	//with MPI I/O while other nodes meanwhile do the same in parallel with MPI_COMM_SELF
	bool stillgood = true;
	vector<SnapshotMetaData>::iterator s;
	for ( s = DatasetInfo.begin(); s != DatasetInfo.end(); ++s ) {
		if ( s->RID == myRank ) { //I take care of the dataset

			//snapshot class reads in the filehand over to the S
			stillgood = readSnapshot( s, false ); //MK::we do not clear the database but accumulate all belonging to me for the operation mode TrackingParallel

			if ( stillgood == false ) { 
				cerr << "ERROR::Worker::ReadDatasets " << this->get_Rank() << " detected a read error for snapshot " << s->extID << endl; 
				return false;
			}
		}
		//success remains true
	}
	return stillgood;
}


bool topoHdl::readDatasets(unsigned int fd)
{
	//constructs local database of Grains and Faces by reading only the file with extID == fd 
	//with MPI I/O while other nodes meanwhile do the same in parallel with MPI_COMM_SELF
	bool stillgood = false;
	vector<SnapshotMetaData>::iterator s;
	for (s = DatasetInfo.begin(); s != DatasetInfo.end(); ++s) {
		if (s->extID != fd)
			continue;

		//not continued means meta data to heavy data <fd> were found
		stillgood = readSnapshot(s, true);
		return stillgood;
	}

	return stillgood;
}


unsigned int topoHdl::WhichOfMyDatasets(unsigned int sid)
{
	//returns if for MPI process has a snapshot with a global external ID of sid if not returns NOT_ONEOFMINE
	for ( unsigned int s = 0; s < this->LocalDB.size(); ++s ) { //MK::was s++ and == sid logic
		if (this->LocalDB[s]->meta.extID != sid)
			continue;

		//not continued means found
		return s;
	}
	return NOT_ONEOFMINE;
}


struct target_props topoHdl::getTargetProperties( unsigned int gid, unsigned int lsid, bool getarea, bool getnfaces, bool gethagbfrac, bool getmobdsee, bool getdsee, bool getnfdiff )
{
	//MK::the function finds the grain only once and reads and computes all possible properties right away to 1-th order neighbors to improve data reutilization
	unsigned int luid = gid / Settings::LookupMaxGrainIDRange;
	unsigned int lu_size = this->LocalDB[lsid]->LookupTable[luid].size;
	GrainHash* thebucket = this->LocalDB[lsid]->LookupTable[luid].bucket;
	unsigned int whichmr = NOT_ASSIGNED_YET;
	unsigned int whichidx = NOT_ASSIGNED_YET;

	for ( unsigned int g = 0; g < lu_size; g++ ) { //scan most interesting candidates
		if (thebucket[g].gid != gid)
			continue;
		//not continued, .gid == gid so candidate found
		whichmr = thebucket[g].mrg_idx;
		whichidx = thebucket[g].idx;
		break;
	}

	struct target_props tp(GRAINVOLUME_ON_ERROR, HAGBFRACTION_ON_ERROR, MOBDSEE_ON_ERROR, DSEE_ON_ERROR, NFACES_ON_ERROR, NFACESDIFF_ON_ERROR);
	if ( whichmr == NOT_ASSIGNED_YET || whichidx == NOT_ASSIGNED_YET ) 
		return tp;

	//read the grain and the relevant properties at once
	Grain* ag = &(this->LocalDB[lsid]->mrg[whichmr]->GrainBucket[whichidx]);

	if ( getarea == true )
		tp.area = ag->size;

	if ( getnfaces == true )
		tp.nfacesk1 = ag->nfaces_identified;

	if ( gethagbfrac == true || getmobdsee == true ) { //disorientation required
		double peri = 0.0;
		double perithreshold = 0.0;
		double wperi_mp = 0.0;
	
		for ( unsigned int nb = 0; nb < ag->nfaces_identified; nb++ ) {
			if ( ag->neighbors[nb].extNBID != THE_DOMAIN_ITSELF ) {
				whichmr = ag->neighbors[nb].mrg_idx; //once overwritten we cannot reidentify in this function the grain gid...
				whichidx = ag->neighbors[nb].idx;

				double nbsee = this->LocalDB[lsid]->mrg[whichmr]->GrainBucket[whichidx].see;
				double disori = misorientationCubic( ag->phi1, ag->Phi, ag->phi2,  this->LocalDB[lsid]->mrg[whichmr]->GrainBucket[whichidx].phi1, this->LocalDB[lsid]->mrg[whichmr]->GrainBucket[whichidx].Phi, this->LocalDB[lsid]->mrg[whichmr]->GrainBucket[whichidx].phi2 );

				//total grain boundary segment length (2d), area (3d) not to boundary grains
				peri = peri + ag->neighbors[nb].size;

				//compute HAGB fraction
				if ( disori >= Settings::HAGBDetectionThreshold )
					perithreshold = perithreshold + ag->neighbors[nb].size;

				//compute Mob2SEE
				//segment length averaged stored elastic energy induced effective migration speed
				wperi_mp = wperi_mp + ( Settings::HAGBMobility * dis2mob(disori) * (Settings::DislocEnPerM*(nbsee - ag->see) * ag->neighbors[nb].size) ); //relative mobility N/m^^2*m^2/m^2=N/m^2
			}
		}

		if ( peri > DBL_EPSILON ) {
			tp.hagbfraction = perithreshold / peri;
			tp.mob_dsee = wperi_mp / peri;
		}
	}

	if (getdsee == true) {
		double peri = 0.0;
		double wperi_p = 0.0;

		for (unsigned int nb = 0; nb < ag->nfaces_identified; nb++) {
			if (ag->neighbors[nb].extNBID != THE_DOMAIN_ITSELF) {
				whichmr = ag->neighbors[nb].mrg_idx; //once overwritten we cannot reidentify in this function the grain gid...
				whichidx = ag->neighbors[nb].idx;

				double nbsee = this->LocalDB[lsid]->mrg[whichmr]->GrainBucket[whichidx].see;
				//total grain boundary segment length (2d), area (3d) not to boundary grains
				peri = peri + ag->neighbors[nb].size;
				//segment length averaged stored elastic energy induced effective driving force
				wperi_p = wperi_p + ( Settings::DislocEnPerM * (nbsee - ag->see) * ag->neighbors[nb].size ); //relative mobility N/m^2*m^2/m^2=N/m^2=Pa
			}
		}
		if (peri > DBL_EPSILON) {
			tp.dsee = wperi_p / peri;
		}
	}

	if (getnfdiff == true ) {
		int nf_diff = 0;
		int me_nf = 0;
		int nb_nf = 0;
		if ( ag->nfaces_identified < static_cast<unsigned int>(std::numeric_limits<int>::max()) ) { //most likely 
			me_nf = ag->nfaces_identified;

			for (unsigned int nb = 0; nb < ag->nfaces_identified; nb++) {
				if (ag->neighbors[nb].extNBID != THE_DOMAIN_ITSELF) {
					whichmr = ag->neighbors[nb].mrg_idx; //once overwritten we cannot reidentify in this function the grain gid...
					whichidx = ag->neighbors[nb].idx;

					if ( this->LocalDB[lsid]->mrg[whichmr]->GrainBucket[whichidx].nfaces_identified < static_cast<unsigned int>(std::numeric_limits<int>::max()) ) {
						nb_nf = this->LocalDB[lsid]->mrg[whichmr]->GrainBucket[whichidx].nfaces_identified;
						nf_diff += (me_nf - nb_nf);
					}
					else {
						std::cout << "WARNING::Grain with gID " << this->LocalDB[lsid]->mrg[whichmr]->GrainBucket[whichidx].extGID << " has more faces than describable with datatype int!" << std::endl;
					}
				}
			}
		}
		else {
			std::cout << "WARNING::Grain with gID " << ag->extGID << " has more faces than describable with datatype int!" << std::endl;
		}

		tp.nfacesk1diff = nf_diff;
	}

	return tp;
}


struct rxstats topoHdl::calcRXFraction( unsigned int lsid, std::vector<unsigned int>* dgr, std::vector<unsigned int>* rxg ) {
	//first of all find all grains in the list dgr and accumulate their volume
	struct rxstats rs;
	rs.ndefg = dgr->size();
	rs.nrxg = rxg->size();

	unsigned int ndg = dgr->size(); //compute domain area covered by matrix first including the rxg
	for ( unsigned int dg = 0; dg < ndg; dg++  ) {
		//find the grain in my data set
		unsigned int luid = dgr->at(dg) / Settings::LookupMaxGrainIDRange;
		unsigned int lu_size = this->LocalDB[lsid]->LookupTable[luid].size;
		GrainHash* thebucket = this->LocalDB[lsid]->LookupTable[luid].bucket;
		unsigned int whichmr = NOT_ASSIGNED_YET;
		unsigned int whichidx = NOT_ASSIGNED_YET;

		for ( unsigned int g = 0; g < lu_size; g++ ) { //scan only interesting candidates
			if ( thebucket[g].gid != dgr->at(dg) ) 
				continue;
			//not continued, item found
			whichmr = thebucket[g].mrg_idx;
			whichidx = thebucket[g].idx;
			break;
		}
		if ( whichmr == NOT_ASSIGNED_YET || whichidx == NOT_ASSIGNED_YET ) { continue; } //##MK::is comment is wrong cout << "WARNING during RXEVO matrix grain " << dgr->at(dg) << " not found!" << endl; break; }

		rs.totalsize_elimbnd += this->LocalDB[lsid]->mrg[whichmr]->GrainBucket[whichidx].size;
	}

	unsigned nrg = rxg->size();
	for ( unsigned int rg = 0; rg < nrg; rg++ ) { //find only the grains that are considered as to be recrystallized
		unsigned int luid = rxg->at(rg) / Settings::LookupMaxGrainIDRange;
		unsigned int lu_size = this->LocalDB[lsid]->LookupTable[luid].size;
		GrainHash* thebucket = this->LocalDB[lsid]->LookupTable[luid].bucket;
		unsigned int whichmr = NOT_ASSIGNED_YET;
		unsigned int whichidx = NOT_ASSIGNED_YET;

		for ( unsigned int g = 0; g < lu_size; g++ ) { //scan most interesting candidates
			if ( thebucket[g].gid != rxg->at(rg) ) 
				continue;

			whichmr = thebucket[g].mrg_idx;
			whichidx = thebucket[g].idx;
			break;
		}
		if ( whichmr == NOT_ASSIGNED_YET || whichidx == NOT_ASSIGNED_YET ) { continue; } //##MK::this comment is wrong cout << "WARNING during RXEVO rxed grain " << rxg->at(rg) << " not found!" << endl; break; }

		rs.totalsize_targets += this->LocalDB[lsid]->mrg[whichmr]->GrainBucket[whichidx].size;
	}

	//compute RX fraction
	if ( rs.totalsize_elimbnd > RXEVO_EPSILON ) { 
		rs.X = rs.totalsize_targets / rs.totalsize_elimbnd; //how much of the total area
	}

	return rs;
}


struct nucmodel_bh topoHdl::analyzeNuc_BaileyHirsch( unsigned int gid, unsigned int lsid ) {
	struct nucmodel_bh res;
	res.grainfound = false;

	//find the grain and his neighbors
	unsigned int luid = gid / Settings::LookupMaxGrainIDRange;
	unsigned int lu_size = this->LocalDB[lsid]->LookupTable[luid].size;
	GrainHash* thebucket = this->LocalDB[lsid]->LookupTable[luid].bucket;
	unsigned int whichmr = NOT_ASSIGNED_YET;
	unsigned int whichidx = NOT_ASSIGNED_YET;

	for ( unsigned int g = 0; g < lu_size; g++ ) { //scan most interesting candidates
		if ( thebucket[g].gid != gid ) continue;

		whichmr = thebucket[g].mrg_idx;
		whichidx = thebucket[g].idx;
		break;
	}

	if ( whichmr == NOT_ASSIGNED_YET || whichidx == NOT_ASSIGNED_YET ) { return res; }

	//grain found, so read all relevant properties at once
	res.grainfound = true;

	Grain* ag = &(this->LocalDB[lsid]->mrg[whichmr]->GrainBucket[whichidx]);

	double peri = 0.0;
	double perithreshold = 0.0;
	double hv_hagb_drho_bailey = 0.0;
	double hv_all_drho_bailey = 0.0;
	double hv_hagb_drho_bate = 0.0;
	double hv_all_drho_bate = 0.0;

	//HAGB fraction
	for ( unsigned int nb = 0; nb < ag->nfaces_identified; nb++ ) {
		if ( ag->neighbors[nb].extNBID != THE_DOMAIN_ITSELF ) {
			whichmr = ag->neighbors[nb].mrg_idx;
			whichidx = ag->neighbors[nb].idx;	 //##MK::memory access check?

			double nbsee = this->LocalDB[lsid]->mrg[whichmr]->GrainBucket[whichidx].see;
			double disori = misorientationCubic( ag->phi1, ag->Phi, ag->phi2,  this->LocalDB[lsid]->mrg[whichmr]->GrainBucket[whichidx].phi1, this->LocalDB[lsid]->mrg[whichmr]->GrainBucket[whichidx].Phi, this->LocalDB[lsid]->mrg[whichmr]->GrainBucket[whichidx].phi2 );

			peri = peri + ag->neighbors[nb].size; //boundary segment length

			if ( disori >= Settings::HAGBDetectionThreshold ) {
				perithreshold = perithreshold + ag->neighbors[nb].size;
				hv_hagb_drho_bailey = hv_hagb_drho_bailey + (nbsee * Heaviside1Minusf(nbsee, ag->see) * ag->neighbors[nb].size );
				hv_hagb_drho_bate = hv_hagb_drho_bate + (nbsee * pow(Heaviside1Minusf(nbsee, ag->see), 0.5) * ag->neighbors[nb].size );
			}

			hv_all_drho_bailey = hv_all_drho_bailey + (nbsee * Heaviside1Minusf(nbsee, ag->see) * ag->neighbors[nb].size );
			hv_all_drho_bate = hv_all_drho_bate + (nbsee * pow(Heaviside1Minusf(nbsee, ag->see), 0.5) * ag->neighbors[nb].size );

//cout << "\t\t\t\t nbsee/agsee/eval/nbsize = " << nbsee << ";" << ag->see << ";" << Heaviside1Minusf( nbsee, ag->see ) << ";" << ag->neighbors[nb].size << endl;
		}
	}

	if ( peri > DBL_EPSILON ) {
		res.hagbfraction = perithreshold / peri;

		hv_all_drho_bailey = hv_all_drho_bailey / peri;
		hv_all_drho_bate = hv_all_drho_bate / peri;
	}

	if ( perithreshold > DBL_EPSILON ) {
		hv_hagb_drho_bailey = hv_hagb_drho_bailey / perithreshold;
		hv_hagb_drho_bate = hv_hagb_drho_bate / perithreshold;
	}

	//double mult = 1.0 / PI;
	//if ( Settings::Dimensionality == TWODIMENSIONAL )		mult = mult * SQR(Settings::PhysicalDomainSize );
	//if ( Settings::Dimensionality == THREEDIMENSIONAL )		mult = mult * CUBE(Settings::PhysicalDomainSize );

	//compute radii and critical nuclei sizes according to Bailey Hirsch and Bate Hutchinson/Bailey Hirsch model
	if ( Settings::Dimensionality == TWODIMENSIONAL ) { 
		res.r = pow ( (ag->size * SQR(Settings::PhysicalDomainSize) / PI), 0.5 );
		if ( res.hagbfraction >= 0.25 ) res.r_hagb25 = res.r; //setting of this value means, this is a sub-grain with at least 25% high-angle grain boundary
		if ( res.hagbfraction >= 0.50 ) res.r_hagb50 = res.r;
		if ( res.hagbfraction >= 0.75 ) res.r_hagb75 = res.r;
	}
	if ( Settings::Dimensionality == THREEDIMENSIONAL ) {
		res.r = pow ( (ag->size * CUBE(Settings::PhysicalDomainSize) / ((4.0/3.0)*PI) ), (1.0/3.0) );
		if ( res.hagbfraction >= 0.25 ) res.r_hagb25 = res.r; //setting of this value means, this is a sub-grain with at least 25% high-angle grain boundary
		if ( res.hagbfraction >= 0.50 ) res.r_hagb50 = res.r;
		if ( res.hagbfraction >= 0.75 ) res.r_hagb75 = res.r;
	}

	if ( hv_hagb_drho_bailey >= BH_MINIMUM_RHO ) res.rc_baihir_hagb = ( 4.0 * Settings::HAGBEnergy ) / ( Settings::DislocEnPerM * hv_hagb_drho_bailey );
	if ( hv_all_drho_bailey >= BH_MINIMUM_RHO ) res.rc_baihir_all = ( 4.0 * Settings::HAGBEnergy ) / ( Settings::DislocEnPerM * hv_all_drho_bailey );
	if ( hv_hagb_drho_bate >= BH_MINIMUM_RHO ) res.rc_bathut_hagb = ( 4.0 * Settings::HAGBEnergy ) / ( Settings::DislocEnPerM * hv_hagb_drho_bate );
	if ( hv_all_drho_bate >= BH_MINIMUM_RHO ) res.rc_bathut_all = ( 4.0 * Settings::HAGBEnergy ) / ( Settings::DislocEnPerM * hv_all_drho_bate );

//cout << "\t\tlsid/gid/found/r//hagb/25/50/75//rcritBaiH/ABatH/A//peri/perithreshold/hv values " << lsid << ";" << gid << ";" << (int) res.grainfound << "\t\t" << res.r << ";" << res.hagbfraction << ";" << res.r_hagb25 << ";" << res.r_hagb50 << ";" << res.r_hagb75 << "\t\t" << res.rc_baihir_hagb << ";" << res.rc_baihir_all << ";" << res.rc_bathut_hagb << ";" << res.rc_bathut_all << "\t\t" << peri << ";" << perithreshold << ";" << hv_hagb_drho_bailey << ";" << hv_all_drho_bailey << ";" << hv_hagb_drho_bate << ";" << hv_all_drho_bate << endl; 

	return res;
}
	

std::vector<unsigned int>* topoHdl::analyze_elimbnd_all_self(std::string logfn_method, std::vector<unsigned int>* BlackList, bool UseGrainIDWhiteList)
{
	//MK::utilizes only process-locally the bndbuf --- if existent --- without MPI communication to identify
	//a list of all grains in all time steps which never touched the boundary and are not in the Blacklist
	//the order of the list values is the same for each process because the order in this->bndbuf is the same for all processes
	if (this->bndbuf == NULL) {
		return NULL;
	}

	std::vector<unsigned int>* targets = NULL;
	try { targets = new std::vector<unsigned int>; }
	catch (std::bad_alloc &exc) {
		cerr << "ERROR::Worker::AnalyzeElimbndAllSelf " << this->get_Rank() << " I cannot build a TargetList" << endl;
		return NULL;
	}
	
	//utilize a BlackListInclusionHash to identify quicker than naive if a grain is in a black list
	std::vector<bool>* BlackListInclusionHash = NULL;
	if (BlackList != NULL) {
		try { BlackListInclusionHash = new vector<bool>; }
		catch (std::bad_alloc &exc) {
			cerr << "ERROR::Worker " << this->get_Rank() << " I cannot build a BlackListInclusionHash" << endl; //fall back to imperformant solution
		}
		if (BlackListInclusionHash != NULL) {
			unsigned int ngr = 1 + Settings::LargestGrainID;
			BlackListInclusionHash->reserve(ngr);
			for (unsigned int gid = 0; gid < ngr; ++gid) {
				BlackListInclusionHash->push_back(false);
			}
			for (unsigned int tgid = 0; tgid < BlackList->size(); tgid++) {
				BlackListInclusionHash->at(BlackList->at(tgid)) = true;
			}
		}
	}

	//generate the filter list
	//is there no white list to consider? whitelist means take these grains
	if (UseGrainIDWhiteList == false) {
		//is there a blacklist or not
		if (BlackList == NULL) {
			for (unsigned int c = 0; c < this->bndbufsize; ++c) {
				if ( this->bndbuf[c].whenhitboundary > Settings::SnapshotLast ) { //grain never made boundary contact
					targets->push_back(c);
				}
				//else we do not want to consider it!
			}
		}
		else { //consider blacklist as well
			for (unsigned int c = 0; c < this->bndbufsize; ++c) { //start at 1 to exclude THE_DOMAIN_ITSELF,i.e. the zero grain explicitly
				if (this->bndbuf[c].whenhitboundary > Settings::SnapshotLast ) {
					if (BlackListInclusionHash != NULL) {
						if (BlackListInclusionHash->at(c) == false) { //not in blacklist
							targets->push_back(c);
						}
						//else gID is in BlackList and shall not be considered...
					}
					else { //not Hash was generated, so fallback to imperformant checking of whether gID c is in BlackList or not
						bool inblacklist = false;
						for (unsigned int black = 0; black < BlackList->size(); ++black) {
							if ( c != BlackList->at(black))
								continue;
							inblacklist = true;
						} //all have to items of the blacklist have to be tested
						if (inblacklist == false) {
							targets->push_back(c);
						}
					}
				}
			} //next grain
		}
	}
	else { //UseGrainIDWhiteList == true
		for (unsigned int whiteid = 0; whiteid < this->GrainIDWhiteList.size(); ++whiteid) {
			unsigned int c = GrainIDWhiteList[whiteid];
			if (BlackList == NULL) {
				if ( this->bndbuf[c].whenhitboundary > Settings::SnapshotLast) {
					targets->push_back(c);
				}
			}
			else { //consider blacklist as well
				if (this->bndbuf[c].whenhitboundary > Settings::SnapshotLast) {
					if (BlackListInclusionHash != NULL) {
						if (BlackListInclusionHash->at(c) == false) { //not in blacklist
							targets->push_back(c);
						}
						//else gID is in BlackList and shall not be considered...
					}
					else { //no Hash was generated, so fallback to imperformant checking of whether gID c is in BlackList or not
						bool inblacklist = false;
						for (unsigned int black = 0; black < BlackList->size(); ++black) {
							if (c != BlackList->at(black))
								continue;
							inblacklist = true;
						}
						if (inblacklist == false) {
							targets->push_back(c);
						}
					}
				}
			}
		} //next whiteid
	}

	//Hash no longer needed
	if (BlackListInclusionHash != NULL) { 
		delete BlackListInclusionHash; 
		BlackListInclusionHash = NULL;
	}


	if (targets->size() < 1) {
		cerr << "ERROR::Worker::AnalyzeElimbndAllSelf " << this->get_Rank() << " I build a target list but it contains no IDs!" << endl;
		if (targets != NULL) { delete targets; targets = NULL; }
		return NULL;
	}

	//no we are sure targets exists so store them to file
	string log_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + "." + logfn_method + ".All.GrainsElimbnd.Rank." + std::to_string(this->myRank) + ".csv";
	ofstream loglast;
	loglast.open(log_fn.c_str(), std::ofstream::out | std::ofstream::trunc);
	if (loglast.is_open() == true) {
		loglast << "GrainID\n";
		for (unsigned int g = 0; g < targets->size(); ++g) 
			loglast << targets->at(g) << "\n";	
		loglast.flush();
		loglast.close();
	}
	else {
		cerr << "ERROR::Worker::AnalyzeElimBndAll " << this->get_Rank() << " unable to open logfile " << log_fn << endl;
	}

	return targets;
	//MK::DO NOT FORGET TO "delete targets" OUTSIDE THIS FUNCTION!
}


std::vector<unsigned int>* topoHdl::analyze_elimbnd_all_db( std::string logfn_method, std::vector<unsigned int>* BlackList, bool UseGrainIDWhiteList )
{
	//MK::utilizes coorperative the process-local databases and MPI communication to identify
	//cooperatively over all processes compiles a list of all grains in all time steps which never touched the boundary and are not in the Blacklist
	//the order of the list values is the same for each process
	//which do not touch the boundary, one rank distributes to all other
	//MK::size of the boundary bucket has to be 1+Settings::LargestGrainID because it is used as a hash with extGID...

	//find the max(GID) over all datasets
	unsigned int ngr = 0;
	for ( unsigned int s = 0; s < this->LocalDB.size(); s++ ) {
		for ( unsigned int mr = 0; mr < this->LocalDB[s]->mrg.size(); mr++ ) {
			unsigned int mrngr = this->LocalDB[s]->mrg[mr]->GrainBucket.size();
			for ( unsigned int g = 0; g < mrngr; ++g ) {
				if (this->LocalDB[s]->mrg[mr]->GrainBucket[g].extGID < ngr)
					continue;
				//not continued >=
				ngr = this->LocalDB[s]->mrg[mr]->GrainBucket[g].extGID;
//if ( this->LocalDB[s]->mrg[mr]->GrainBucket[g].extGID >= ngr ) //automatically the zero grain is expelled as it has ID = 0
//	ngr = this->LocalDB[s]->mrg[mr]->GrainBucket[g].extGID;
	}}}
	unsigned int worldngr = 0;
	if (this->get_nRanks() > SINGLE_PROCESS) {
		MPI_Allreduce(&ngr, &worldngr, 1, MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);
	}
	else { worldngr = ngr; }
	//cout << "Worker " << myRank << " my maximum ID is " << ngr << " world maximum ID is " << worldngr << endl;

		//1+worldngr because accounting for zero grain
		//recvbuffer and sendbuf are utilized as bitflag fields which indicate BOUNDARY_CONTACT or not
		//##MK::further optimatization pack 8bits at once...
	char* recvbuf = NULL;
	try { recvbuf = new char[1 + worldngr]; }
	catch (std::bad_alloc &exc) {
		cerr << "ERROR::Worker::AnalyzeElimBndAll " << this->get_Rank() << " allocation of recvbuf failed!" << endl; return NULL;
	}
	char* sendbuf = NULL;
	try { sendbuf = new char[1 + worldngr]; }
	catch (std::bad_alloc &exc) {
		cerr << "ERROR::Worker::AnalyzeElimBndAll " << this->get_Rank() << " allocation of recvbuf failed!" << endl;
		if (recvbuf != NULL) { delete[] recvbuf; recvbuf = NULL; }
		return NULL;
	}

	//##MK::minimalistic it is fact that only MASTER needs a recvbuffer but not a sendbuf as he places the data directly in the recvbuffer
	//how we initialize a recv and sendbuf in allow all ranks to work on a recvbuf to compile the list of targets on their own while all others do the same
	//these are ID state fields
	for (unsigned int c = 0; c < 1 + worldngr; ++c) { recvbuf[c] = DO_ANALYZE; }
	for (unsigned int c = 0; c < 1 + worldngr; ++c) { sendbuf[c] = DO_ANALYZE; }

	//detect in local portion of dataset whether a grain at some point touches the boundary
	unsigned int nelim = 0;
	if (myRank == MASTER) { //master operates inplace on recvbuffer
		for (unsigned int s = 0; s < this->LocalDB.size(); s++) {
			for (unsigned int mr = 0; mr < this->LocalDB[s]->mrg.size(); mr++) {
				unsigned int mrngr = this->LocalDB[s]->mrg[mr]->GrainBucket.size();
				for (unsigned int g = 0; g < mrngr; ++g) {
					if (this->LocalDB[s]->mrg[mr]->GrainBucket[g].boundary != BOUNDARY_CONTACT)
						continue;
					if (recvbuf[this->LocalDB[s]->mrg[mr]->GrainBucket[g].extGID] != DO_NOT_ANALYZE)
						nelim++; //prevent multiple counting
					recvbuf[this->LocalDB[s]->mrg[mr]->GrainBucket[g].extGID] = DO_NOT_ANALYZE;
				}
			}
		}
	}
	else { //slaves operate on sendbuf
		for (unsigned int s = 0; s < this->LocalDB.size(); s++) {
			for (unsigned int mr = 0; mr < this->LocalDB[s]->mrg.size(); mr++) {
				unsigned int mrngr = this->LocalDB[s]->mrg[mr]->GrainBucket.size();
				for (unsigned int g = 0; g < mrngr; ++g) {
					if (this->LocalDB[s]->mrg[mr]->GrainBucket[g].boundary != BOUNDARY_CONTACT)
						continue;
					if (sendbuf[this->LocalDB[s]->mrg[mr]->GrainBucket[g].extGID] != DO_NOT_ANALYZE)
						nelim++; //to prevent double counting
					sendbuf[this->LocalDB[s]->mrg[mr]->GrainBucket[g].extGID] = DO_NOT_ANALYZE;
				}
			}
		}
	}
	//separation into master and not the master allows to use the MPI_IN_PLACE option during the MPI_Reduce process

//cout << "Worker " << this->get_Rank() << " eliminated " << nelim << " grains." << endl;

	if (this->get_nRanks() > SINGLE_PROCESS) {
		if (myRank == MASTER) {
			MPI_Reduce(MPI_IN_PLACE, recvbuf, 1 + worldngr, MPI_CHAR, MPI_BAND, MASTER, MPI_COMM_WORLD);
			//BAND in combination with 0x00 for DO_NOT_ANALYZE and 0x01 for DO_ANALYZE sets all bits to 0 and therefore the entire chars in sendbuf to 0x00 if not all bits in all chars for one grain in the processes are set to 0x01
			//i.e. any process with the buffer[] entry being set to DO_NOT_ANALYZE will cause the expelling of the grain
		}
		else {
			MPI_Reduce(sendbuf, NULL, 1 + worldngr, MPI_CHAR, MPI_BAND, MASTER, MPI_COMM_WORLD);
		}
		//now the MASTER knows in recvbuf which targets to consider and which not

		MPI_Barrier(MPI_COMM_WORLD); //##MK::may be obsolete
		MPI_Bcast(recvbuf, 1 + worldngr, MPI_CHAR, MASTER, MPI_COMM_WORLD);
	}
	else {
		//MASTER already has all pieces of information inplace in recvbuffer
	}

	//analyze which ids whom to consider
	std::vector<unsigned int>* targets = NULL;
	targets = new std::vector<unsigned int>;

	//utilize a BlackListInclusionHash to identify quicker than naive if a grain is in a black list
	std::vector<bool>* BlackListInclusionHash = NULL;
	if (BlackList != NULL) {
		try { BlackListInclusionHash = new vector<bool>; }
		catch (std::bad_alloc &exc) {
			cerr << "ERROR::Worker " << this->get_Rank() << " I cannot build a BlackListInclusionHash" << endl; //fall back to imperformant solution
		}
		if (BlackListInclusionHash != NULL) {
			unsigned int ngr = 1 + Settings::LargestGrainID;
			BlackListInclusionHash->reserve(ngr);
			for (unsigned int gid = 0; gid < ngr; ++gid)
				BlackListInclusionHash->push_back(false);
			for (unsigned int tgid = 0; tgid < BlackList->size(); tgid++)
				BlackListInclusionHash->at(BlackList->at(tgid)) = true;
		}
	}

	//is there no white list to consider? whitelist means take these grains
	if (UseGrainIDWhiteList == false) {
		//is there a blacklist or not
		if (BlackList == NULL) {
			for (unsigned int c = 1; c < (1 + worldngr); ++c) { //start at 1 to exclude THE_DOMAIN_ITSELF,i.e. the zero grain explicitly
				if (recvbuf[c] == DO_ANALYZE) //##MK::earlier version silently read c <= 1+worldngr which may cause memory leak...
					targets->push_back(c);
			}
		}
		else { //consider blacklist as well
			for (unsigned int c = 1; c < (1 + worldngr); ++c) { //start at 1 to exclude THE_DOMAIN_ITSELF,i.e. the zero grain explicitly
				if (recvbuf[c] == DO_ANALYZE) {
					if (BlackListInclusionHash != NULL) {
						if (BlackListInclusionHash->at(c) == false) { //not in blacklist
							targets->push_back(c);
						}
						//else gID is in BlackList and shall not be considered...
					}
					else { //not Hash was generated, so fallback to imperformant checking of whether gID c is in BlackList or not
						bool inblacklist = false;
						for (unsigned int black = 0; black < BlackList->size(); ++black) {
							if (c != BlackList->at(black))
								continue;
							inblacklist = true;
						} //all have to items of the blacklist have to be tested
						if (inblacklist == false) {
							targets->push_back(c);
						}
					}
				}
			} //next grain
		}
	}
	else { //UseGrainIDWhiteList == true
		for (unsigned int whiteid = 0; whiteid < this->GrainIDWhiteList.size(); ++whiteid) {
			unsigned int c = GrainIDWhiteList[whiteid];
			if (BlackList == NULL) {
				if (recvbuf[c] == DO_ANALYZE) {
					targets->push_back(c);
				}
			}
			else { //consider blacklist as well
				if (recvbuf[c] == DO_ANALYZE) {
					if (BlackListInclusionHash != NULL) {
						if (BlackListInclusionHash->at(c) == false) { //not in blacklist
							targets->push_back(c);
						}
						//else gID is in BlackList and shall not be considered...
					}
					else { //no Hash was generated, so fallback to imperformant checking of whether gID c is in BlackList or not
						bool inblacklist = false;
						for (unsigned int black = 0; black < BlackList->size(); ++black) {
							if (c != BlackList->at(black))
								continue;
							inblacklist = true;
						}
						if (inblacklist == false) {
							targets->push_back(c);
						}
					}
				}
			}
		} //next whiteid
	}

	if (BlackListInclusionHash != NULL) {
		delete BlackListInclusionHash;
		BlackListInclusionHash = NULL;
	}

	if (targets->size() < 1) {
		cerr << "ERROR::Worker::AnalyzeElimBndAllDB " << this->get_Rank() << " I build a target list but it contains no IDs!" << endl;
		if (targets != NULL) { delete targets; targets = NULL; }
		if (recvbuf != NULL) { delete[] recvbuf; recvbuf = NULL; }
		if (sendbuf != NULL) { delete[] sendbuf; sendbuf = NULL; }
		return NULL;
	}

	string log_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + "." + logfn_method + ".All.GrainsElimbnd.Rank." + std::to_string(this->myRank) + ".csv";
	ofstream loglast;
	loglast.open(log_fn.c_str(), std::ofstream::out | std::ofstream::trunc);
	if (loglast.is_open() == true) {
		loglast << "GrainID\n";
		for (unsigned int g = 0; g < targets->size(); ++g) { loglast << targets->at(g) << "\n"; }
		loglast.flush();
		loglast.close();
	}
	else {
		cerr << "ERROR::Worker::AnalyzeElimBndAll " << this->get_Rank() << " unable to open logfile " << log_fn << endl;
	}
	
	if (recvbuf != NULL) { delete[] recvbuf; recvbuf = NULL; }
	if (sendbuf != NULL) { delete[] sendbuf; sendbuf = NULL; }
	
	//MK::may be obsolete
	MPI_Barrier ( MPI_COMM_WORLD );

	return targets;
	//MK::DO NOT FORGET TO "delete targets" OUTSIDE THIS FUNCTION!
}


std::vector<unsigned int>* topoHdl::analyze_elimbnd_one_mpiioself(std::string logfn_method, std::vector<unsigned int>* OnlyTheseTargets, unsigned int extID)
{
	//MK::requires this->grbuf already loaded with Texture_<extID>.bin !
	if (this->grbuf == NULL) {
		cerr << "ERROR::Worker::AnalyzeElimOneMPIIOSELF " << this->get_Rank() << " this->grbuf was not loaded!" << endl;
		return NULL;
	}

	//allocate memory for target grainIDs
	std::vector<unsigned int>* targets = NULL;
	try { targets = new vector<unsigned int>; }
	catch (std::bad_alloc &exc) {
		cerr << "ERROR::Worker::AnalyzeElimOneMPIIOSELF " << this->get_Rank() << " unable to allocate memory for targets in elimbnd_one_mpiioself" << endl;
		return NULL;
	}

	//build WhiteListHash for OnlyTheseTargets
	vector<bool>* WhiteListHash = NULL;
	if (OnlyTheseTargets != NULL) {
		try { WhiteListHash = new vector<bool>; }
		catch (std::bad_alloc &exc) {
			cerr << "ERROR::Worker::AnalyzeElimOneMPIIOSELF " << this->get_Rank() << " unable to build WhiteListHash!" << endl;
			//WhiteListHash remains NULL so algo will fall back to imperformant linear filtering of targets...
		}
		if (WhiteListHash != NULL) {
			unsigned int ngr = 1 + Settings::LargestGrainID;
			WhiteListHash->reserve(ngr);
			for (unsigned int gid = 0; gid < ngr; gid++) 
				WhiteListHash->push_back(false); //fill with false
			for (unsigned int tgid = 0; tgid < OnlyTheseTargets->size(); tgid++) 
				WhiteListHash->at(OnlyTheseTargets->at(tgid)) = true; //reset desired to true 
		}
	}

	//filter targets based potentially OnlyTheseTargets
	for (unsigned int g = 0; g < this->grbufsize; g++) {
		if (this->grbuf[g].id != THE_DOMAIN_ITSELF && this->grbuf[g].intersectsBoundaryGrain != BOUNDARY_CONTACT) {
			unsigned int gid = this->grbuf[g].id; //also one of those to be included?

			if (OnlyTheseTargets == NULL) { //no restriction of the target gID
				targets->push_back(gid);
			}
			else { //a constraining list of target gID was passed
				if (WhiteListHash != NULL) { //and we have an hash to speed up the identification process via the WhiteListHash
					if (WhiteListHash->at(gid) == true) {
						targets->push_back(gid);
					}
				}
				else { //fall back to imperformant sequential scanning of the WhiteList for finding whether to consider or not
					for (unsigned int ot = 0; ot < OnlyTheseTargets->size(); ++ot) {  //scan whether the rxgrain is include in the list of grains which never made boundary contact
						if (gid != (*OnlyTheseTargets).at(ot))
							continue;
						//not continued means gid was found in the list OnlyTheseTargets and is therefore considered
						targets->push_back(gid);
						break;
					}
				}
			} 
		} //filtering of g done
	}
	
	//WhiteListHash if existent is now no longer needed
	if (WhiteListHash != NULL) {
		delete WhiteListHash; WhiteListHash = NULL;
	}

	//avoid passing list without any targets
	if ( targets->size() < 1 ) {
		cerr << "ERROR::Worker::AnalyzeElimBndOneMPIIOSELF " << this->get_Rank() << " I build a target list but it contains no IDs!" << endl;
		if (targets != NULL) { delete targets; targets = NULL; }
		return NULL;
	}

	//DEBUG write results
	string log_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + "." + logfn_method + ".F." + std::to_string(Settings::SnapshotFirst) + ".O." + std::to_string(Settings::SnapshotOffset) + ".L." + std::to_string(Settings::SnapshotLast) + ".GrainElimbnd.Rank." + std::to_string(this->myRank) + ".csv";
	ofstream loglast;
	loglast.open(log_fn.c_str(), std::ofstream::out | std::ofstream::trunc);
	if (loglast.is_open() == true) {
		loglast << "GrainID\n"; //##MK::obviously old flaw Area(allPropsTimeStep = " << extID << "); SEE; Phi1; PHI; Phi2; x; y; z\n";
		for (unsigned int g = 0; g < targets->size(); ++g)
			loglast << targets->at(g) << "\n";
		loglast.flush();
		loglast.close();
	}
	else {
		cerr << "ERROR::Worker::AnalyzeElimBndOne " << this->get_Rank() << " unable to write " << log_fn << endl;
	}

	return targets;
}


std::vector<unsigned int>* topoHdl::analyze_elimbnd_one_db( std::string logfn_method, std::vector<unsigned int>* OnlyTheseTargets, unsigned int extID ) {
	//cooperatively over all processes compiles a list of all grains in time step extID do not touch the boundary
	//utilizes process local databases instead of performing I/O on heavy data
	//and are additionally included in the gID list OnlyTheseTargets, if such list is passed, i.e. OnlyTheseTargets != NULL 
	//one rank distributes to all other, list values are consecutive for all processes
	unsigned int whichrank = snapshot2rank( extID );
	if ( whichrank == UNABLE_TO_IDENTIFY_RANK ) { return NULL; }

	std::vector<unsigned int>* targets = NULL;
	unsigned int ntargets = NOT_YET_FOUND;
	unsigned int starget = NOT_YET_FOUND;

	if ( myRank == whichrank ) { //rank who has the data gets number and GIDs of targets
		for ( unsigned int s = 0; s < this->LocalDB.size(); s++ ) {
			if ( LocalDB[s]->meta.extID == extID ) {
				starget = s;
				break;
			}
		}
		MPI_Bcast( &starget, 1, MPI_UNSIGNED, whichrank, MPI_COMM_WORLD );
	}
	else {
		MPI_Bcast( &starget, 1, MPI_UNSIGNED, whichrank, MPI_COMM_WORLD );
	}
	if ( starget == NOT_YET_FOUND ) { return NULL; } //##MK:: was stargets, all check whether the dataset was found

	//scans for targets and outputs help file
	if ( myRank == whichrank ) {
		cout << "...Worker " << this->get_Rank() << " identifies targets on LocalDB " << starget << endl;

		//build WhiteListHash for OnlyTheseTargets, ##MKfor the moment
		vector<bool>* WhiteListHash = NULL;
		if (OnlyTheseTargets != NULL) {
			try { WhiteListHash = new vector<bool>; }
			catch (std::bad_alloc &exc) {
				cerr << "ERROR::Worker " << this->get_Rank() << " unable to build WhiteListHash!" << endl;
				//WhiteListHash remains NULL so algo will fall back to imperformant linear filtering of targets...
			}
			if (WhiteListHash != NULL) {
				unsigned int ngr = 1 + Settings::LargestGrainID;
				WhiteListHash->reserve(ngr);
				for (unsigned int gid = 0; gid < ngr; gid++) //fill with false
					WhiteListHash->push_back(false);
				for (unsigned int tgid = 0; tgid < OnlyTheseTargets->size(); tgid++) 
					WhiteListHash->at(OnlyTheseTargets->at(tgid)) = true;
			}
		}
		
		//filter targets based potentially OnlyTheseTargets
		targets = new std::vector<unsigned int>;

		for ( unsigned int mr = 0; mr < this->LocalDB[starget]->mrg.size(); mr++ ) {
			MemRegion* themr = this->LocalDB[starget]->mrg[mr];
			for ( unsigned int g = 0; g < themr->GrainBucket.size(); g++ ) {
				if ( themr->GrainBucket[g].extGID != THE_DOMAIN_ITSELF && themr->GrainBucket[g].boundary != BOUNDARY_CONTACT ) {
					//also one of those to be included?
					unsigned int gid = themr->GrainBucket[g].extGID;
					if (OnlyTheseTargets == NULL) { //no restriction of the target gID
						targets->push_back(gid);
					}
					else { //a constraining list of target gID was passed
						if (WhiteListHash != NULL) {
							if (WhiteListHash->at(gid) == true) {
								targets->push_back(themr->GrainBucket[g].extGID);
							}
						}
						else { //fall back to imperformant list scanning solution
							for (unsigned int ot = 0; ot < OnlyTheseTargets->size(); ++ot) {  //scan whether the rxgrain is include in the list of grains which never made boundary contact
								if (gid != (*OnlyTheseTargets).at(ot))
									continue;
								//not continued means gid was found in the list OnlyTheseTargets and is therefore considered
								targets->push_back(themr->GrainBucket[g].extGID);
								break;
							}
						}
					} //filtering of targets done
				}
			} //next grain from bucket
		} //next memregion

		//WhiteListHash no longer needed
		if (WhiteListHash != NULL) {
			delete WhiteListHash; WhiteListHash = NULL;
		}

		if ( targets != NULL ) 
			ntargets = targets->size();

		MPI_Bcast( &ntargets, 1, MPI_UNSIGNED, whichrank, MPI_COMM_WORLD );
	}
	else { //slaves wait for whichrank to send
		MPI_Bcast( &ntargets, 1, MPI_UNSIGNED, whichrank, MPI_COMM_WORLD );
	}

	if ( ntargets == NOT_YET_FOUND || ntargets < 1 ) { //no targets or other errors
		cerr << "ERROR::Worker::AnalyzeElimBndOne " << this->get_Rank() << " no target identified" << endl;
		if (myRank == whichrank) {
			if (targets != NULL) { delete targets; targets = NULL; }
		}
		return NULL;
	}

	//synchronize list of targets to all
	MPI_Barrier( MPI_COMM_WORLD ); //##MK::may not be necessary!

	unsigned int* gidtargets = NULL;
	try { gidtargets = new unsigned int[ntargets]; }
	catch (std::bad_alloc &exc) {
		cerr << "ERROR::Worker::AnalyzeElimBndOne " << this->get_Rank() << " allocation of gidtargets buffer failed!" << endl;
		if (targets != NULL) { //check necessary!
			delete targets; targets = NULL;
		}
		return NULL;
	}
	
	if ( myRank == whichrank ) { //whichrank broadcasts to slaves
		for (unsigned int gs = 0; gs < ntargets; ++gs ) 
			gidtargets[gs] = targets->at(gs);

		MPI_Bcast( gidtargets, ntargets, MPI_UNSIGNED, whichrank, MPI_COMM_WORLD );
	}
	else { 
		MPI_Bcast( gidtargets, ntargets, MPI_UNSIGNED, whichrank, MPI_COMM_WORLD );	//slaves listen to whichrank

		//now also the slaves initialize their copy of the targets
		targets = new std::vector<unsigned int>;
		for ( unsigned int gs = 0; gs < ntargets; ++gs ) { targets->push_back( gidtargets[gs] ); }
	}

	if (gidtargets != NULL) {
		delete [] gidtargets; gidtargets = NULL;
	}
	
	//master writes results
	if ( myRank == whichrank ) {
		string log_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + "." + logfn_method + ".OneBL.GrainElimbnd.Rank." + std::to_string(this->myRank) + ".csv";
		ofstream loglast;
		loglast.open( log_fn.c_str(), std::ofstream::out | std::ofstream::trunc );
		if (loglast.is_open() == true) {
			loglast << "GrainID\n"; //##MK::obviously old flaw Area(allPropsTimeStep = " << extID << "); SEE; Phi1; PHI; Phi2; x; y; z\n";

			for (unsigned int g = 0; g < targets->size(); ++g)
				loglast << targets->at(g) << "\n";

			loglast.flush();
			loglast.close();
		}
		else {
			cerr << "ERROR::Worker::AnalyzeElimBndOne " << this->get_Rank() << " unable to write " << log_fn << endl;
		}
	}

	MPI_Barrier( MPI_COMM_WORLD ); //##MK::may not be necessary!

	return targets;
	//MK::DO NOT FORGET TO "delete targets" OUTSIDE THIS FUNCTION!
}


void topoHdl::analyze_trajectories_fw_db( void )
{
	//targets are only those grains which in their lifetime never touch the boundary
	//hence first a cooperative analysis is necessary to identify which grains did not touch the boundary
	//for this task the processes inspect first their own datasets and identify for these grains whether they still exist and if so then touching the boundary
	//second they synchronize
	double timer = MPI_Wtime();
	
	std::vector<unsigned int>* thetargets = NULL;
	if ( Settings::AnalyzeTrajectoriesForwardMode == E_FORWARD_ALL ) {
		thetargets = this->analyze_elimbnd_all_db( "FW", NULL, false ); //everybody has to call this, dont forget to delete the pointer!
	}
	else if ( Settings::AnalyzeTrajectoriesForwardMode == E_FORWARD_SELECTED ) {
		thetargets = this->analyze_elimbnd_all_db( "FW", NULL, true );
	}
	else {
		cerr << "ERROR::Analyze trajectories forward impossible, invalid mode!" << endl;
		if (thetargets != NULL) { delete thetargets; thetargets = NULL; }
		return;
	}

	//MK::strictly here should be MPI communication to insure process consistency!
	if ( thetargets->size() < 1 ) { 
		cerr << "ERROR::Worker " << this->get_Rank() << " Analyze trajectories forward impossible, no grain survived" << endl; 
		if (thetargets != NULL) { delete thetargets; thetargets = NULL; }
		return;
	}

	//as all ranks derived the targets collectively with the same algorithm and calling arguments the vector thetargets has the same ID order in each rank

	int ncol = ((Settings::SnapshotLast - Settings::SnapshotFirst) / Settings::SnapshotOffset + 1);
	int nrow = thetargets->size();
	cout << "Writing rows (ngrains) = " << nrow << " columns (time) = " << ncol << endl;

	//create I/O buffer container to write timeslot after timeslot, only as many write ops as timeslots because one node has all grains
	unsigned int* nfacesbuffer = NULL;
	try { nfacesbuffer = new unsigned int[nrow]; }
	catch (std::bad_alloc &exc) {
		cerr << "ERROR::Worker::AnalyzeForwardTracking " << this->get_Rank() << " allocating memory for nfacesbuffer failed!" << endl;
		if (thetargets != NULL) { delete thetargets; thetargets = NULL; }
	}
	double* volbuffer = NULL;
	try { volbuffer = new double[nrow]; }
	catch (std::bad_alloc &exc) {
		cerr << "ERROR::Worker::AnalyzeForwardTracking " << this->get_Rank() << " allocating memory for volbuffer failed!" << endl;
		if (thetargets != NULL) { delete thetargets; thetargets = NULL; }
		if (nfacesbuffer != NULL) { delete [] nfacesbuffer; nfacesbuffer = NULL; }
	}
	double* hagbfracbuffer = NULL;
	try { hagbfracbuffer = new double[nrow]; }
	catch (std::bad_alloc &exc) {
		cerr << "ERROR::Worker::AnalyzeForwardTracking " << this->get_Rank() << " allocating memory for hagbfracbuffer failed!" << endl;
		if (thetargets != NULL) { delete thetargets; thetargets = NULL; }
		if (nfacesbuffer != NULL) { delete[] nfacesbuffer; nfacesbuffer = NULL; }
		if (volbuffer != NULL) { delete[] volbuffer; volbuffer = NULL; }
	}
	double* mobdseebuffer = NULL;
	try { mobdseebuffer = new double[nrow]; }
	catch (std::bad_alloc &exc) {
		cerr << "ERROR::Worker::AnalyzeForwardTracking " << this->get_Rank() << " allocating memory for mobdseebuffer failed!" << endl;
		if (thetargets != NULL) { delete thetargets; thetargets = NULL; }
		if (nfacesbuffer != NULL) { delete[]  nfacesbuffer; nfacesbuffer = NULL; }
		if (volbuffer != NULL) { delete[]  volbuffer; volbuffer = NULL; }
		if (hagbfracbuffer != NULL) { delete[]  hagbfracbuffer; hagbfracbuffer = NULL; }
	}
	
	//open three files volume and nfaces and see
	MPI_File ioHdlFAC, ioHdlVOL, ioHdlHAGB, ioHdlMOBDSEE;
	MPI_Status ioStaFAC, ioStaVOL, ioStaHAGB, ioStaMOBDSEE;

	string suffix =  "F." + std::to_string(Settings::SnapshotFirst) + ".O." + std::to_string(Settings::SnapshotOffset) + ".L." + std::to_string(Settings::SnapshotLast) + ".NC." + std::to_string(ncol) + ".NR." + std::to_string(nrow) + ".bin";
	string faces_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".FW.NF." + suffix;
	string vol_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".FW.VOL." + suffix; 
	string hagbfrac_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".FW.HAGB." + suffix;
	string mobdsee_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".FW.MOBDSEE." + suffix;
	
	// open the file in create and write-only mode, ##only MASTER should finally write the file
	MPI_File_open(MPI_COMM_WORLD, faces_fn.c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlFAC);
	MPI_File_open(MPI_COMM_WORLD, vol_fn.c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlVOL);
	MPI_File_open(MPI_COMM_WORLD, hagbfrac_fn.c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlHAGB);
	MPI_File_open(MPI_COMM_WORLD, mobdsee_fn.c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlMOBDSEE);

	long long totalOffsetFaces = 0;
	long long totalOffsetVolume = 0;
	long long totalOffsetHAGB = 0;
	long long totalOffsetMOBDSEE = 0;
	
	unsigned int fid = UNKNOWN_ID;
	unsigned int localid = UNKNOWN_ID;

	//all nodes process in parallel, processes do not share time-resolved data
	for ( unsigned int f = Settings::SnapshotFirst; f <= Settings::SnapshotLast; f += Settings::SnapshotOffset ) {
		fid = (f/Settings::SnapshotOffset) - ( Settings::SnapshotFirst / Settings::SnapshotOffset );

		localid = WhichOfMyDatasets( f );

//cout << "My local ID is = " << localid << endl;

		if ( localid != NOT_ONEOFMINE ) {
			for ( unsigned int t = 0; t < thetargets->size(); t++ ) { //equivalent to gid < nrow
				struct target_props targetp;
				targetp = this->getTargetProperties( thetargets->at(t), localid, true, true, true, true, false, false );
				//MK::for this to work correctly the order in thetargets has to be the same in every rank!
				nfacesbuffer[t] = targetp.nfacesk1;
				volbuffer[t] = targetp.area;
				hagbfracbuffer[t] = targetp.hagbfraction;
				mobdseebuffer[t] = targetp.mob_dsee;
			}

			//calculate address in the file //globalid of the dataset * nrow
			totalOffsetFaces = fid * nrow * 4;
			totalOffsetVolume = fid * nrow * 8;
			totalOffsetHAGB = fid * nrow * 8;
			totalOffsetMOBDSEE = fid * nrow * 8;

			MPI_File_write_at( ioHdlFAC, totalOffsetFaces, nfacesbuffer, nrow, MPI_UNSIGNED, &ioStaFAC);
			MPI_File_write_at( ioHdlVOL, totalOffsetVolume, volbuffer, nrow, MPI_DOUBLE, &ioStaVOL);
			MPI_File_write_at( ioHdlHAGB, totalOffsetHAGB, hagbfracbuffer, nrow, MPI_DOUBLE, &ioStaHAGB);
			MPI_File_write_at( ioHdlMOBDSEE, totalOffsetMOBDSEE, mobdseebuffer, nrow, MPI_DOUBLE, &ioStaMOBDSEE);
		}
		//else, nothing because different node will take care of it...
	}

	if (thetargets != NULL) { delete thetargets; thetargets = NULL; }
	if (nfacesbuffer != NULL) { delete[]  nfacesbuffer; nfacesbuffer = NULL; }
	if (volbuffer != NULL) { delete[]  volbuffer; volbuffer = NULL; }
	if (hagbfracbuffer != NULL) { delete[]  hagbfracbuffer; hagbfracbuffer = NULL; }
	if (mobdseebuffer != NULL) { delete[] mobdseebuffer; mobdseebuffer = NULL; }

	MPI_Barrier(MPI_COMM_WORLD);

	MPI_File_close(&ioHdlFAC);
	MPI_File_close(&ioHdlVOL);
	MPI_File_close(&ioHdlHAGB);
	MPI_File_close(&ioHdlMOBDSEE);

	cout << "...Worker " << this->get_Rank() << " analyzing trajectories forward = " << (MPI_Wtime() - timer) << " seconds" << endl;
}


inline unsigned int topoHdl::snapshot2rank( unsigned int extID )
{
	for ( unsigned int s = 0; s < this->DatasetInfo.size(); s++ ) {
		if ( DatasetInfo[s].extID != extID ) 
			continue;

		return DatasetInfo[s].RID;
	}
	cerr << "ERROR::Unable to identify rank in snapshot2rank for extID " << extID << endl;
	return UNABLE_TO_IDENTIFY_RANK;
}


unsigned int topoHdl::ngr_elimbnd( unsigned int extID )
{
	unsigned int starget = NOT_YET_FOUND;
	for ( unsigned int s = 0; s < this->LocalDB.size(); s++ ) {
		if (LocalDB[s]->meta.extID != extID)
			continue;
		//not continued, found
		starget = s;
		break;
	}
	if ( starget == NOT_YET_FOUND ) return NOT_YET_FOUND;

	unsigned int ngr = 0;
	for ( unsigned int mr = 0; mr < this->LocalDB[starget]->mrg.size(); mr++ ) {
		MemRegion* themr = this->LocalDB[starget]->mrg[mr];
		for ( unsigned int g = 0; g < themr->GrainBucket.size(); g++ ) {
			if ( themr->GrainBucket[g].extGID != THE_DOMAIN_ITSELF && themr->GrainBucket[g].boundary != BOUNDARY_CONTACT ) {
				ngr++;
			}
		}
	}
	return ngr;
}


void topoHdl::gid_elimbnd( unsigned int extID, unsigned int* res )
{
	//collect Grains in Snapshot extID that have no boundary contact
	unsigned int starget = NOT_YET_FOUND;
	for ( unsigned int s = 0; s < this->LocalDB.size(); s++ ) {
		if (LocalDB[s]->meta.extID != extID)
			continue;
		//not continued so found
		starget = s;
		break;
	}
	if ( starget == NOT_YET_FOUND ) return; //leaving all elements in res invalid

cout << "...Worker::AnalyzeBackwardsTracking " << this->get_Rank() << " identifies targets on LocalDB " << starget << endl;
	
	//collect the candidates and write the meta data directly to reutilize cache content
	string log_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".BK.GrainElimbnd.Rank." + std::to_string(this->get_Rank()) + ".csv";
	ofstream loglast;
	loglast.open(log_fn.c_str(), std::ofstream::out | std::ofstream::trunc);
	if (loglast.is_open() == true) {
		loglast << "GrainID;Area;SEE;Phi1;PHI;Phi2;x;y;z\n";
		loglast << ";1;1/m^2;rad;rad;rad;1;1;1\n";

		unsigned int ngr = 0;
		for (unsigned int mr = 0; mr < this->LocalDB[starget]->mrg.size(); mr++) {
			MemRegion* themr = this->LocalDB[starget]->mrg[mr];
			for (unsigned int g = 0; g < themr->GrainBucket.size(); g++) {
				if (themr->GrainBucket[g].extGID != THE_DOMAIN_ITSELF && themr->GrainBucket[g].boundary != BOUNDARY_CONTACT) {
					res[ngr] = themr->GrainBucket[g].extGID;
					ngr++;
					loglast << themr->GrainBucket[g].extGID << ";" << themr->GrainBucket[g].size << ";" << themr->GrainBucket[g].see << ";" << themr->GrainBucket[g].phi1 << ";" << themr->GrainBucket[g].Phi << ";" << themr->GrainBucket[g].phi2 << ";" << themr->GrainBucket[g].x << ";" << themr->GrainBucket[g].y << ";" << themr->GrainBucket[g].z << "\n";
				}
			}
		}
		loglast.flush();
		loglast.close();
	}
}


void topoHdl::analyze_trajectories_bk_db( void )
{
	double timer = MPI_Wtime();

	unsigned int rankSnapshotLast = snapshot2rank( Settings::SnapshotLast );
	if ( rankSnapshotLast == UNABLE_TO_IDENTIFY_RANK ) { return; }

	unsigned int nGrainsSurvived = NOT_YET_FOUND; //MK::value of this constant sets hard limit for maximum eligible number of grains 2billion...
	if ( myRank == rankSnapshotLast ) //MK::ranks do not share raw data other than meta
		nGrainsSurvived = ngr_elimbnd( Settings::SnapshotLast );

	MPI_Bcast( &nGrainsSurvived, 1, MPI_UNSIGNED, rankSnapshotLast, MPI_COMM_WORLD );
	if ( nGrainsSurvived == NOT_YET_FOUND ) { cerr << "ERROR::Worker::AnalyzeBackwardTracking " << this->get_Rank() << " analyze trajectories backward failed, data inconsistency!" << endl; return; }
	if ( nGrainsSurvived < 1 ) { cerr << "ERROR::Worker::AnalyzeBackwardTracking " << this->get_Rank() << " analyze trajectories backward impossible, no grain survived" << endl; return; }

	//now all ranks know how many grains to consider
	unsigned int localstatus = 1;
	unsigned int* gidGrainsSurvived = NULL;
	//all ranks allocate memory to store the IDs of the survivors
	try { gidGrainsSurvived = new unsigned int[nGrainsSurvived]; }
	catch (std::bad_alloc &exc) { localstatus = 0; }
	unsigned int worldstatus = 0;
	MPI_Allreduce(&localstatus, &worldstatus, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
	if (worldstatus != this->get_nRanks()) { //if not all have the memory available
		if (gidGrainsSurvived != NULL) { delete[] gidGrainsSurvived; gidGrainsSurvived = NULL; }
		return;
	}
	//MK::consistence check
	for (unsigned int gs = 0; gs < nGrainsSurvived; ++gs ) 
		gidGrainsSurvived[gs] = INVALID;

	//all ranks have memory to store the IDs
	if ( myRank == rankSnapshotLast ) {
		gid_elimbnd ( Settings::SnapshotLast, gidGrainsSurvived );
		MPI_Bcast( gidGrainsSurvived, nGrainsSurvived, MPI_UNSIGNED, rankSnapshotLast, MPI_COMM_WORLD );
	}
	else {
		MPI_Bcast( gidGrainsSurvived, nGrainsSurvived, MPI_UNSIGNED, rankSnapshotLast, MPI_COMM_WORLD );
	}

	//MK::the Bcast of the specifically arranged list of grain IDs from the MASTER to the slaves assures that all processes
	//have the same view which gID is at the first, the second, the third position, and so forth...

	//##MK::debug, additional check may be obsolete
	unsigned int localinvalid = 0;
	for (unsigned int target = 0; target < nGrainsSurvived; ++target) {
		if (gidGrainsSurvived[target] != INVALID)
			continue;
		//not continued, so INVALID
		localinvalid++;
		//cerr << "ERROR::Worker::AnalyzeBackwardTrajectory " << this->get_Rank() << " analyze trajectories backward impossible because of data inconsistency for targets" << endl;
	}
	unsigned int worldinvalid = 0;
	MPI_Allreduce(&localinvalid, &worldinvalid, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
	if (worldinvalid != 0) {
		cerr << "ERROR::Worker::AnalyzeBackwardTrajectory " << this->get_Rank() << " analyze trajectories backward impossible because of data inconsistency for targets" << endl;
		if (gidGrainsSurvived != NULL) { delete[] gidGrainsSurvived; gidGrainsSurvived = NULL; }
		return;
	}

	//now all ranks know all targets with their gIDs in an array in the same order, so lets start processing
	int ncol = ((Settings::SnapshotLast - Settings::SnapshotFirst) / Settings::SnapshotOffset + 1);
	int nrow = nGrainsSurvived;
	cout << "...Worker::AnalyzeBackwardsTracking " << this->get_Rank() << " writing rows (ngrains) = " << nrow << " columns (time) = " << ncol << endl;
	
	//create I/O buffer container to write timeslot after timeslot, only as many write ops as timeslots because one node has all grains
	unsigned int* nfacesbuffer = NULL;
	try { nfacesbuffer = new unsigned int[nrow]; }
	catch (std::bad_alloc &exc) {
		cerr << "ERROR::Worker::AnalyzeBackwardsTracking " << this->get_Rank() << " allocating memory for nfacesbuffer failed!" << endl;
		if (gidGrainsSurvived != NULL) { delete[] gidGrainsSurvived; gidGrainsSurvived = NULL; }
	}
	double* volbuffer = NULL;
	try { volbuffer = new double[nrow]; }
	catch (std::bad_alloc &exc) {
		cerr << "ERROR::Worker::AnalyzeBackwardsTracking " << this->get_Rank() << " allocating memory for volbuffer failed!" << endl;
		if (gidGrainsSurvived != NULL) { delete[] gidGrainsSurvived; gidGrainsSurvived = NULL; }
		if (nfacesbuffer != NULL) { delete[] nfacesbuffer; nfacesbuffer = NULL; }
	}
	double* hagbfracbuffer = NULL;
	try { hagbfracbuffer = new double[nrow]; }
	catch (std::bad_alloc &exc) {
		cerr << "ERROR::Worker::AnalyzeBackwardsTracking " << this->get_Rank() << " allocating memory for hagbfracbuffer failed!" << endl;
		if (gidGrainsSurvived != NULL) { delete[] gidGrainsSurvived; gidGrainsSurvived = NULL; }
		if (nfacesbuffer != NULL) { delete[] nfacesbuffer; nfacesbuffer = NULL; }
		if (volbuffer != NULL) { delete[] volbuffer; volbuffer = NULL; }
	}
	double* mobdseebuffer = NULL;
	try { mobdseebuffer = new double[nrow]; }
	catch (std::bad_alloc &exc) {
		cerr << "ERROR::Worker::AnalyzeBackwardsTracking " << this->get_Rank() << " allocating memory for mobdseebuffer failed!" << endl;
		if (gidGrainsSurvived != NULL) { delete[] gidGrainsSurvived; gidGrainsSurvived = NULL; }
		if (nfacesbuffer != NULL) { delete[]  nfacesbuffer; nfacesbuffer = NULL; }
		if (volbuffer != NULL) { delete[]  volbuffer; volbuffer = NULL; }
		if (hagbfracbuffer != NULL) { delete[]  hagbfracbuffer; hagbfracbuffer = NULL; }
	}


	//open three files volume and nfaces and see
	MPI_File ioHdlFAC, ioHdlVOL, ioHdlHAGB, ioHdlMOBDSEE;
	MPI_Status ioStaFAC, ioStaVOL, ioStaHAGB, ioStaMOBDSEE;

	string suffix = "F." + std::to_string(Settings::SnapshotFirst) + ".O." + std::to_string(Settings::SnapshotOffset) + ".L." + std::to_string(Settings::SnapshotLast) + ".NC." + std::to_string(ncol) + ".NR." + std::to_string(nrow) + ".bin";
	string faces_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".BK.NF." + suffix;
	string vol_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".BK.VOL." + suffix;
	string hagbfrac_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".BK.HAGB." + suffix;
	string mobdsee_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".BK.MOBDSEE." + suffix;
	
	// open the file in create and write-only mode, ##only MASTER should finally write the file
	MPI_File_open(MPI_COMM_WORLD, faces_fn.c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlFAC);
	MPI_File_open(MPI_COMM_WORLD, vol_fn.c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlVOL);
	MPI_File_open(MPI_COMM_WORLD, hagbfrac_fn.c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlHAGB);
	MPI_File_open(MPI_COMM_WORLD, mobdsee_fn.c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlMOBDSEE);

	long long totalOffsetFaces = 0;
	long long totalOffsetVolume = 0;
	long long totalOffsetHAGB = 0;
	long long totalOffsetMOBDSEE = 0;

	unsigned int fid = UNKNOWN_ID;
	unsigned int localid = UNKNOWN_ID;
	unsigned int gid = NOT_YET_KNOWN;

	//all nodes process in parallel, processes do not share time-resolved data
	for ( unsigned int f = Settings::SnapshotFirst; f <= Settings::SnapshotLast; f += Settings::SnapshotOffset ) {
		fid = (f/Settings::SnapshotOffset) - ( Settings::SnapshotFirst / Settings::SnapshotOffset );

		localid = WhichOfMyDatasets( f );

//cout << "My local ID is = " << localid << endl;

		if ( localid != NOT_ONEOFMINE ) { //so my responsibility
			for ( unsigned int target = 0; target < nGrainsSurvived; target++ ) {
				unsigned int gid = gidGrainsSurvived[target]; //assured to be a grain not at the boundary and not THE_DOMAIN_ITSELF
				struct target_props targetp;
				targetp = this->getTargetProperties( gid, localid, true, true, true, true, false, false );
				nfacesbuffer[target] = targetp.nfacesk1;
				volbuffer[target] = targetp.area;
				hagbfracbuffer[target] = targetp.hagbfraction;
				mobdseebuffer[target] = targetp.mob_dsee;
				//##MK::further optimization potential write single function that requires only once to identify the grain and its neighbors
			}

			//calculate address in the file //globalid of the dataset * nrow
			totalOffsetFaces = fid * nrow * 4;
			totalOffsetVolume = fid * nrow * 8;
			totalOffsetHAGB = fid * nrow * 8;
			totalOffsetMOBDSEE = fid * nrow * 8;

			MPI_File_write_at( ioHdlFAC, totalOffsetFaces, nfacesbuffer, nrow, MPI_UNSIGNED, &ioStaFAC);
			MPI_File_write_at( ioHdlVOL, totalOffsetVolume, volbuffer, nrow, MPI_DOUBLE, &ioStaVOL);
			MPI_File_write_at( ioHdlHAGB, totalOffsetHAGB, hagbfracbuffer, nrow, MPI_DOUBLE, &ioStaHAGB);
			MPI_File_write_at( ioHdlMOBDSEE, totalOffsetMOBDSEE, mobdseebuffer, nrow, MPI_DOUBLE, &ioStaMOBDSEE);
		}
		//else, nothing because different node will take care of it...
	}

	if (gidGrainsSurvived != NULL) { delete[] gidGrainsSurvived; gidGrainsSurvived = NULL; }
	if (nfacesbuffer != NULL) { delete[]  nfacesbuffer; nfacesbuffer = NULL; }
	if (volbuffer != NULL) { delete[]  volbuffer; volbuffer = NULL; }
	if (hagbfracbuffer != NULL) { delete[]  hagbfracbuffer; hagbfracbuffer = NULL; }
	if (mobdseebuffer != NULL) { delete[] mobdseebuffer; mobdseebuffer = NULL; }
	
	MPI_Barrier(MPI_COMM_WORLD);

	MPI_File_close(&ioHdlFAC);
	MPI_File_close(&ioHdlVOL);
	MPI_File_close(&ioHdlHAGB);
	MPI_File_close(&ioHdlMOBDSEE);

	cout << "...Worker::AnalyzeTrackingBackwardss " << this->get_Rank() << " analyzing trajectories backwards = " << (MPI_Wtime() - timer) << " seconds" << endl;
}


struct descr_stats topoHdl::getDescrStats( std::vector<unsigned int>* cand, unsigned int lsid, bool getarea )
{
	struct descr_stats res;
	for ( unsigned int c = 0; c < cand->size(); c++ ) {
		//find grain cgid
		unsigned int cgid = cand->at(c);

		unsigned int luid = cgid / Settings::LookupMaxGrainIDRange;
		unsigned int lu_size = this->LocalDB[lsid]->LookupTable[luid].size;
		GrainHash* thebucket = this->LocalDB[lsid]->LookupTable[luid].bucket;
		unsigned int whichmr = NOT_ASSIGNED_YET;
		unsigned int whichidx = NOT_ASSIGNED_YET;

		for ( unsigned int g = 0; g < lu_size; g++ ) { //scan most interesting candidates
			if ( thebucket[g].gid != cgid ) 
				continue;

			whichmr = thebucket[g].mrg_idx;
			whichidx = thebucket[g].idx;
			break;
		}
		if ( whichmr == NOT_ASSIGNED_YET || whichidx == NOT_ASSIGNED_YET ) //next one if grain not there
			continue;

		//add value of existing grains only cgid to
		res.n++;
		res.sum += this->LocalDB[lsid]->mrg[whichmr]->GrainBucket[whichidx].size;
	}

	if (res.n > 0)
		res.mean = res.sum / res.n;
	
	return res;
}


#define MATRIX_GRAINS_REMAINING		0
#define MATRIX_GRAINS_MEAN			1

void topoHdl::analyze_sizegain_vs_matrix_bk_db( void )
{
	//for all survivors identify in each time step how their size (area/volume) is relative to the mean of the matrix, 
	//the matrix are all grains that do not touch the boundary and are not included in the list of survivors
	//survivors - all the grains without boundary contact in the final data set, most likely recrystallizing
	//matrix - all grains ever considered with never boundary contact but excluding the survivors
	double timer = MPI_Wtime();
	
	std::vector<unsigned int>* survivors = NULL;
	survivors = analyze_elimbnd_one_db("SizeGainVsMatrixSurvivors", NULL, Settings::SnapshotLast);

	unsigned int localstatus = 1;
	if ( survivors == NULL ) { localstatus = 0; }
	unsigned int worldstatus = 0;
	MPI_Allreduce( &localstatus, &worldstatus, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD );
	if (worldstatus != nRanks) {
		cerr << "ERROR::Worker::SzGainVsMatrix " << this->get_Rank() << " analysis SizeGain unable to determine survivors!" << endl;
		if (survivors != NULL) { delete survivors; survivors = NULL; }
		return;
	}

	std::vector<unsigned int>* matrix = NULL;
	matrix = analyze_elimbnd_all_db( "SizeGainVsMatrixSurvivors", survivors, false );

	localstatus = 1;
	if ( matrix == NULL ) { localstatus = 0; }
	worldstatus = 0;
	MPI_Allreduce( &localstatus, &worldstatus, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD );
	if (worldstatus != nRanks) {
		cerr << "ERROR::Worker::SzGainVsMatrix " << this->get_Rank() << " analysis SizeGain unable to determine matrix!" << endl;
		if (survivors != NULL) { delete survivors; survivors = NULL; }
		if (matrix != NULL) { delete matrix; matrix = NULL; }
		return;
	}

	//now all ranks know how many survivors and matrix grains there are and their ids, so lets start processing
	int ncol = ((Settings::SnapshotLast - Settings::SnapshotFirst) / Settings::SnapshotOffset + 1);
	int nrow = 2 + survivors->size(); //+2 for the number of grains in the mean value and the mean value itself
	cout << "...Worker::SzGainVsMatrix " << this->get_Rank() << " writing rows (ngrains) = " << nrow << " columns (time) = " << ncol << endl;

	//create I/O buffer container to write timeslot after timeslot, only as many write ops as timeslots because one node has all grains
	double* relsizebuffer = NULL;
	try { relsizebuffer = new double[nrow]; }
	catch (std::bad_alloc &exc) { 
		cerr << "ERROR::Worker::SzGainVsMatrix " << this->get_Rank() << " allocating memory in analyse_sizegain relsize" << endl;
		if (survivors != NULL) { delete survivors; survivors = NULL; }
		if (matrix != NULL) { delete matrix; matrix = NULL; }
		return;
	}

	//open three files volume and nfaces and see
	MPI_File ioHdlREL;
	MPI_Status ioStaREL;

	string relsize_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".SIZEGAINBK.F." + std::to_string(Settings::SnapshotFirst) + ".O." + std::to_string(Settings::SnapshotOffset) + ".L." + std::to_string(Settings::SnapshotLast) + ".NC." + std::to_string(ncol) + ".NR." + std::to_string(nrow) + ".bin";
	MPI_File_open(MPI_COMM_WORLD, relsize_fn.c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlREL);

	long long totalOffsetREL = 0;
	unsigned int fid = UNKNOWN_ID;
	unsigned int localid = UNKNOWN_ID;

	//all nodes process in parallel, processes do not share time-resolved data
	for ( unsigned int f = Settings::SnapshotFirst; f <= Settings::SnapshotLast; f += Settings::SnapshotOffset ) {
		fid = (f/Settings::SnapshotOffset) - ( Settings::SnapshotFirst / Settings::SnapshotOffset );

		localid = WhichOfMyDatasets( f );

		if ( localid != NOT_ONEOFMINE ) { //so my responsibility
			for ( unsigned int s = 0; s < survivors->size(); s++ ) { relsizebuffer[s] = 0.0; }

			struct descr_stats matrixstats;
			matrixstats = this->getDescrStats( matrix, localid, true );

			relsizebuffer[MATRIX_GRAINS_REMAINING] = matrixstats.n;
			relsizebuffer[MATRIX_GRAINS_MEAN] = matrixstats.mean;

			if ( matrixstats.mean > DBL_EPSILON ) {
				for ( unsigned int surv = 0; surv < survivors->size(); surv++ ) {
					unsigned int gid = survivors->at(surv);
					struct target_props targetp;
					targetp = this->getTargetProperties( gid, localid, true, false, false, false, false, false );
					relsizebuffer[2+surv] = targetp.area / matrixstats.mean;
				}
			}

			totalOffsetREL = fid * nrow * 8;
			MPI_File_write_at( ioHdlREL, totalOffsetREL, relsizebuffer, nrow, MPI_DOUBLE, &ioStaREL);
		}
		//else, nothing because different node will take care of it...
	}

	if (survivors != NULL) { delete survivors; survivors = NULL; }
	if (matrix != NULL) { delete matrix; matrix = NULL; }
	if (relsizebuffer != NULL) { delete[] relsizebuffer; relsizebuffer = NULL; }

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_File_close(&ioHdlREL);

	cout << "...Worker " << this->get_Rank() << " analyzing SizeGainVsMatrix backwards = " << (MPI_Wtime() - timer) << " seconds" << endl;
}


void topoHdl::analyze_maxsizegain_fw_db( void )
{
	//set up intermediate array of 0.0 on [1:Settings::LargestGrainID]
	unsigned int localhealth = 1;
	unsigned int ngr = 1 + Settings::LargestGrainID;
	double* local_max_size = NULL;
	try { local_max_size = new double[ngr]; } 
	catch (std::bad_alloc& exc) {
		localhealth = 0;
		cerr << "ERROR::Worker::AnalyzeMaxSizeGainFW " << this->get_Rank() << " unable to get memory analyze_maxsizegain_forward!" << endl; return;
	}
	for ( unsigned int gid = 0; gid < ngr; ++gid ) { //initialize local max size to zero
		local_max_size[gid] = 0.0; 
	}
	
	unsigned int globalhealth = 0;
	if (nRanks > SINGLE_PROCESS) {
		MPI_Allreduce(&localhealth, &globalhealth, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
		if (globalhealth != 1 * this->get_nRanks() ) {
			cerr << "ERROR::WorkerAnalyzeMaxSizeGainFW " << this->get_Rank() << " ranks were unable to allocate memory for measuring maximum size!" << endl;
			if (local_max_size != NULL) { delete[] local_max_size; local_max_size = NULL; return; }
		}
	} //SINGLE_PROCESS
	else {
		if (localhealth != 1) {
			cerr << "ERROR::WorkerAnalyzeMaxSizeGainFW " << this->get_Rank() << " MASTER unable to allocate memory for measuring maximum size!" << endl;
			if (local_max_size != NULL) { delete[] local_max_size; local_max_size = NULL; return; }
		}
	}

	//all nodes process through local data in parallel, processes do not share time-resolved data
	unsigned int localid = UNKNOWN_ID;
	for (unsigned int f = Settings::SnapshotFirst; f <= Settings::SnapshotLast; f += Settings::SnapshotOffset) {
		localid = WhichOfMyDatasets(f);

		if (localid != NOT_ONEOFMINE) {

			calc_localmaxsize(localid, local_max_size);

		} //done with local snapshot fid
		//else, nothing because different node will take care of it...
	}

	//wait for all to obtain the local maximum values
	MPI_Barrier(MPI_COMM_WORLD); 

	if (this->get_Rank() == MASTER) {
		MPI_Reduce(MPI_IN_PLACE, local_max_size, ngr, MPI_DOUBLE, MPI_MAX, MASTER, MPI_COMM_WORLD);
	}
	else {
		MPI_Reduce(local_max_size, local_max_size, ngr, MPI_DOUBLE, MPI_MAX, MASTER, MPI_COMM_WORLD);
	}

	//master outputs csv list
	if (this->get_Rank() == MASTER ) {
		string log_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".MaxSizeGainFW.F." + std::to_string(Settings::SnapshotFirst) + ".O." + std::to_string(Settings::SnapshotOffset) + ".L." + std::to_string(Settings::SnapshotLast) + ".csv";
		ofstream logmxsz;
		logmxsz.open(log_fn.c_str(), std::ofstream::out | std::ofstream::trunc);
		if (logmxsz.is_open() == true) {
			//header
			logmxsz << "GrainID;MaxSize\n";
			if(Settings::Dimensionality == TWO_DIMENSIONS)
				logmxsz << ";micron^2\n";
			if(Settings::Dimensionality == THREE_DIMENSIONS)
				logmxsz << ";micron^3\n";

			//normalized area on [0,1]^d to physical size in micron^d
			double norm2real = SQR(METER2MICRON(Settings::PhysicalDomainSize)); //TWO_DIMENSIONS
			if (Settings::Dimensionality == THREE_DIMENSIONS)
				norm2real = CUBE(METER2MICRON(Settings::PhysicalDomainSize));
			
			//data writing
			for (unsigned int gid = 0; gid < ngr; ++gid) {
				if (local_max_size[gid] > DOUBLE_EPSILON) //only grains != BOUNDARY_CONTACT are accounted for their maximum size
					logmxsz << gid << ";" << setprecision(18) << local_max_size[gid] * norm2real << "\n";
			}
			logmxsz.flush();
			logmxsz.close();
		}
		else { cerr << "ERROR::Worker::AnalyzeMaxSizeGainFW " << this->get_Rank() << " is unable to open results file " << log_fn << endl; }
	}
	
	if (local_max_size != NULL) {
		delete[] local_max_size; local_max_size = NULL;
	}

	//##MK::master outputs csv list or binary double string
}


void topoHdl::analyze_knn_naive_db( unsigned int kmax )
{
	//identifies shells of higher order neighbors to target provided that neither these nor the target has boundary contact, i.e. is neighboring THE_DOMAIN_ITSELF
	//##MK::naive in this function is that only one process identifies the neighbors and may only utilize OpenMP on the nodes

	double timer = MPI_Wtime();
	cout << "...Worker::AnalyzeKNN " << this->get_Rank() << " performing k-nearest neighbor analyses for " << kmax << " shells" << endl;

	unsigned int localid = UNKNOWN_ID;
	for (unsigned int f = Settings::SnapshotFirst; f <= Settings::SnapshotFirst; f += Settings::SnapshotOffset) { //##MK::change f <= Settings::SnapshotLast
		localid = WhichOfMyDatasets(f);

		if (localid != NOT_ONEOFMINE) {
			double ltimer = MPI_Wtime();
			//initialize individual output file for each timestep
			string logknn_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".KNN.FID." + std::to_string(f) + ".csv";

			ofstream logknn;
			logknn.open( logknn_fn.c_str() );
			if (logknn.is_open() == true) {
				logknn << "GrainID;x;y;";
				for (unsigned int kj = 0; kj <= kmax; kj++) logknn << "GrainsIn-" << kj << "-thShell" << ";";
				for (unsigned int kj = 0; kj <= kmax; kj++) logknn << "SumArea-" << kj << "-thShell" << ";";
				for (unsigned int kj = 0; kj < kmax; kj++) logknn << "BndLength-" << kj << "-thShell" << ";";
				for (unsigned int kj = 0; kj < kmax; kj++) logknn << "HAGBFraction-" << kj << "-thShell" << ";";
				for (unsigned int kj = 0; kj < kmax; kj++) logknn << "MdSEEwBndLength-" << kj << "-thShell" << ";";
				logknn << endl;

				//timesteps are never shared among processes, hence, only one process will perform a knn
				bool decision = false;

				decision = this->LocalDB[localid]->analyze_knn(kmax, logknn);

				logknn.flush();
				logknn.close();
				cout << "...Worker::AnalyzeKNN " << this->get_Rank() << " kNN analysis was successful in " << (double)(MPI_Wtime() - ltimer) << " seconds" << endl;
			}
		}
		//else, nothing because different node will take care of it...
	} //next timestep
	cout << "...Worker::AnalyzeKNN " << this->get_Rank() << " all kNN analyses were successful in " << (double) (MPI_Wtime() - timer) << " seconds" << endl;
}


unsigned int topoHdl::modf_binning(double disori) {
	unsigned int nbins = Settings::get_MODFBinCount();
	unsigned int binend = Settings::get_MODFBinEnd();

	double b = (disori - Settings::DisoriAngleBinMin) / Settings::DisoriAngleBinWidth;

	if (b >= 0.0) { //so at least Settings::DisoriAngleBinMin, the most likely case
		unsigned int bi = std::ceil(b);
		if (bi < (nbins - 1)) { //most likely case
			return bi;
		}
		//implicit else, map to above heighest value, i.e. > Settings::DisoriAngleBinMax
		return (1 + binend);
	}
	//implicit else, map to below smallest bin value, i.e. <= Settings::DisoriAngleBinMin
	return 0;
}


unsigned int topoHdl::see_binning(double see) {
	unsigned int nbins = Settings::get_SEEBinCount();
	unsigned int binend = Settings::get_SEEBinEnd();

	double b = (see - Settings::SeeBinMin) / Settings::SeeBinWidth;

	if (b >= 0.0) { //so at least Settings::SeeBinMin, the most likely case
		unsigned int bi = std::ceil(b);
		if (bi < (nbins - 1)) { //most likely case
			return bi;
		}
		//implicit else, map to bin for grains above highest value, i.e. > Settings::SeeBinMax
		return (1 + binend);
	}
	//implicit else, map to below smallest bin value, i.e. <= Settings::SeeBinMin
	return 0;
}


void topoHdl::calc_modf( unsigned int lsid, double* binnedmodf ) {

	unsigned int nbins = Settings::get_MODFBinCount();
	
	//clear results buffer
	for (unsigned int b = 0; b < nbins; ++b) {
		binnedmodf[b] = 0.0;
	}

	//##MK::here a simpler version that is not yet #### OpenMP parallel is utilized, better however is to aggregate individual counts over the memory regions
	double totalbndlen = 0.0;
	for ( unsigned int mr = 0; mr < this->LocalDB[lsid]->mrg.size(); mr++ ) {
		unsigned int ngmr = this->LocalDB[lsid]->mrg[mr]->GrainBucket.size();
		for ( unsigned int g = 0; g < ngmr; g++ ) {
			//track only boundary segments between grain A and B when both grains do not make contact with the boundary
			if ( this->LocalDB[lsid]->mrg[mr]->GrainBucket[g].boundary != BOUNDARY_CONTACT ) {

				unsigned int gnb = this->LocalDB[lsid]->mrg[mr]->GrainBucket[g].nfaces_identified;
				double geuler[3] = { this->LocalDB[lsid]->mrg[mr]->GrainBucket[g].phi1, this->LocalDB[lsid]->mrg[mr]->GrainBucket[g].Phi, this->LocalDB[lsid]->mrg[mr]->GrainBucket[g].phi2 };

				//MK::access to the zero grain is safe because the zero grains has nfaces_identified == 0
				for ( unsigned int nb = 0; nb < gnb; nb++ ) { //find neighbors to grain g
					unsigned int whichmr = NOT_ASSIGNED_YET;
					unsigned int whichidx = NOT_ASSIGNED_YET;
					double dbndlen = 0.0;

					whichmr = this->LocalDB[lsid]->mrg[mr]->GrainBucket[g].neighbors[nb].mrg_idx;
					whichidx = this->LocalDB[lsid]->mrg[mr]->GrainBucket[g].neighbors[nb].idx;
					dbndlen = this->LocalDB[lsid]->mrg[mr]->GrainBucket[g].neighbors[nb].size * 0.5; //account for only half the length as pairs A and B know each other and have the same numerical value for size

					if ( whichmr == NOT_ASSIGNED_YET || whichidx == NOT_ASSIGNED_YET ) { 
						cerr << "ERROR::Worker " << this->get_Rank() << " have not found a neighbor that should be there!" << endl; 
						continue; //MK::rather implement a fatal
					}

					//boundary contact? only if the neighbor is not at the boundary!
					if ( this->LocalDB[lsid]->mrg[whichmr]->GrainBucket[whichidx].boundary != BOUNDARY_CONTACT ) { 
						totalbndlen += dbndlen;

						double disori =  misorientationCubic( geuler[0], geuler[1], geuler[2],  this->LocalDB[lsid]->mrg[whichmr]->GrainBucket[whichidx].phi1, this->LocalDB[lsid]->mrg[whichmr]->GrainBucket[whichidx].Phi, this->LocalDB[lsid]->mrg[whichmr]->GrainBucket[whichidx].phi2 );
						unsigned int bin = modf_binning(disori);
						binnedmodf[bin] += dbndlen;
					}
				} //next potential boundary to a neighbor of g
			}
		} //next grain from this memory region
	} //next memory region

	//translate into CDF
	//normalize distribution
	if ( totalbndlen > MODF_EPSILON ) {
		double _totalbndlen = 1.0 / totalbndlen;
		double accumulator = 0.0;
		for ( unsigned int bin = 0; bin < nbins; bin++ ) {
			accumulator += binnedmodf[bin];
			binnedmodf[bin] = accumulator * _totalbndlen;
		}
	}
	else { cerr << "ERROR::Worker " << this->get_Rank() << " normalization of CDF in MODF was impossible!" << endl; }

	cout << "...Worker " << this->get_Rank() << " localID " << lsid << " so many bins " << nbins << " and so much boundary in total " << setprecision(8) << totalbndlen << endl;
}


void topoHdl::calc_seedf( unsigned int lsid, double* binnedseedf ) {

	unsigned int nbins = Settings::get_SEEBinCount();
	
	//clear results buffer
	for ( unsigned int b = 0; b < nbins; ++b ) {
		binnedseedf[b] = 0.0;
	}

	double totalsize = 0.0;
	//##MK::here a simpler version that is not yet #### OpenMP parallel is utilized, better however is to aggregate individual counts over the memory regions
	for ( unsigned int mr = 0; mr < this->LocalDB[lsid]->mrg.size(); mr++ ) {
		unsigned int ngmr = this->LocalDB[lsid]->mrg[mr]->GrainBucket.size();
		for ( unsigned int g = 0; g < ngmr; g++ ) {
			//track, i.e. do not exclude boundary contact grains, because we are interested in the evolution of the total SEE content in the domain
			double localsize = this->LocalDB[lsid]->mrg[mr]->GrainBucket[g].size;
			totalsize += localsize;

			double see = this->LocalDB[lsid]->mrg[mr]->GrainBucket[g].see;
			unsigned int bin = see_binning( see );
			binnedseedf[bin] += localsize;
		} //next grain from this memory region
	} //next memory region

	//translate into CDF
	//normalize distribution
	if ( totalsize > SEE_EPSILON ) {
		double _totalsize = 1.0 / totalsize;
		double accumulator = 0.0;
		for ( unsigned int bin = 0; bin < nbins; bin++ ) {
			accumulator += binnedseedf[bin];
			binnedseedf[bin] = accumulator * _totalsize;
		}
	}
	else { 
		cerr << "ERROR::Worker " << this->get_Rank() << " normalization of CDF in SEE was impossible!" << endl;
	}

	cout << "Worker " << this->get_Rank() << " localID " << lsid << " so many bins " << nbins << " and so much size in total " << setprecision(16) << totalsize << endl;
}


void topoHdl::calc_localmaxsize(unsigned int lsid, double* localbuffer) {
	//MK::do not clear the local buffer as we keep accumulating higher values over the entire snapshot interval [SnapshotFirst;SnapshotLast]
	
	for (unsigned int mr = 0; mr < this->LocalDB[lsid]->mrg.size(); mr++) {
		unsigned int ngmr = this->LocalDB[lsid]->mrg[mr]->GrainBucket.size();
		for (unsigned int g = 0; g < ngmr; g++) {
			unsigned int gid = this->LocalDB[lsid]->mrg[mr]->GrainBucket[g].extGID;
			double localsize = this->LocalDB[lsid]->mrg[mr]->GrainBucket[g].size;

			if (this->LocalDB[lsid]->mrg[mr]->GrainBucket[g].boundary != BOUNDARY_CONTACT) { //most likely case
				if ( localsize > localbuffer[gid] ) {
					localbuffer[gid] = localsize;
				}
			} //else volume remains 0.0 to indicate grain is not considered
		} //next grain from this memory region
	} //next memory region

	//cout << "Worker " << this->get_Rank() << " localID " << lsid << " accumulated local maximum size " << endl;
}


void topoHdl::calc_localavdrvfrc(unsigned int lsid, seeav* localbuffer)
{
	//MK::do not clear the local buffer as we keep accumulating higher values over the entire snapshot interval [SnapshotFirst;SnapshotLast]
	unsigned int llsid = lsid;

	for (unsigned int mr = 0; mr < this->LocalDB[lsid]->mrg.size(); mr++) {
		unsigned int ngmr = this->LocalDB[lsid]->mrg[mr]->GrainBucket.size();
		for (unsigned int g = 0; g < ngmr; g++) {
			unsigned int ggid = this->LocalDB[lsid]->mrg[mr]->GrainBucket[g].extGID;

			if (this->LocalDB[lsid]->mrg[mr]->GrainBucket[g].boundary != BOUNDARY_CONTACT) { //most likely case
				struct target_props tmp = this->getTargetProperties( ggid, llsid, false, false, false, true, true, false );

				localbuffer[ggid].add( tmp.mob_dsee, tmp.dsee );
			}
		} //next grain from this memory region
	} //next memory region

	//cout << "Worker " << this->get_Rank() << " localID " << lsid << " accumulated local maximum size " << endl;
}


void topoHdl::analyze_modf_db( void ) {
	//processes work through their files populate a stripe array with the binned data
	//MK::computing the MODF excluding boundary grains for all time steps on [Settings::SnapshotFirst, Settings::SnapshotLast]
	double timer = MPI_Wtime();

	int ncol = 1 + ((Settings::SnapshotLast - Settings::SnapshotFirst) / Settings::SnapshotOffset + 1); //1 + timesteps because binends added as well
	int nrow = Settings::get_MODFBinCount();

	cout << "...Worker " << this->get_Rank() << " writing rows (MODF bins) = " << nrow << " columns (First col are bin ends, rest are time CDF-data) = " << ncol << endl;

	//create I/O buffer container to write timeslot after timeslot, only as many write ops as timeslots because one node has all grains
	double* modfbuffer = NULL;
	try { modfbuffer = new double[nrow]; }
	catch (std::bad_alloc &exc) {
		cerr << "ERROR::Worker " << this->get_Rank() << " allocating memory in analyze_modf!" << endl; return;
	}
	for ( unsigned int b = 0; b < nrow; ++b ) { 
		modfbuffer[b] = 0.0;
	}

	//open MPI files
	MPI_File ioHdlMODF;
	MPI_Status ioStaMODF;

	string modf_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".MODF.F." + std::to_string(Settings::SnapshotFirst) + ".O." + std::to_string(Settings::SnapshotOffset) + ".L." + std::to_string(Settings::SnapshotLast) + ".NC." + std::to_string(ncol) + ".NR." + std::to_string(nrow) + ".bin";
	
	// open the file in create and write-only mode, ##only MASTER should finally write the file
	MPI_File_open(MPI_COMM_WORLD, modf_fn.c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlMODF);

	long long totalOffsetMODF = 0;
	unsigned int fid = UNKNOWN_ID;
	unsigned int localid = UNKNOWN_ID;

	//MASTER fills in bin end values
	if ( myRank == MASTER ) {
		for ( unsigned int b = 0; b < (nrow-1); b++ ) { 
			modfbuffer[b] = Settings::DisoriAngleBinMin + RAD2DEG(((double) b) * Settings::DisoriAngleBinWidth); //user input is in degree
		}
		modfbuffer[nrow-1] = RAD2DEG(MODF_INFTY); //back in degrees for the user

		//limit last bin to 1.0, //if ( (((double) (nrow-1+1)) * Settings::DisoriAngleBinWidth) >= 1.0 ) modfbuffer[nrow-1] = 1.0;

		//write first columns which are the binends
		totalOffsetMODF = 0 * nrow * 8;
		MPI_File_write_at( ioHdlMODF, totalOffsetMODF, modfbuffer, nrow, MPI_DOUBLE, &ioStaMODF);
	}

	//all nodes process in parallel, processes do not share time-resolved data
	for( unsigned int f = Settings::SnapshotFirst; f <= Settings::SnapshotLast; f += Settings::SnapshotOffset ) {
		fid = (f / Settings::SnapshotOffset) - ( Settings::SnapshotFirst / Settings::SnapshotOffset );

		localid = WhichOfMyDatasets( f );

//cout << "My local ID is = " << localid << endl;

		if ( localid != NOT_ONEOFMINE ) { //so my responsibility
			
			calc_modf( localid, modfbuffer );

			totalOffsetMODF = (1 * nrow * 8) + (fid * nrow * 8); //MK::offset nrow*8 :because binends preceed

			MPI_File_write_at( ioHdlMODF, totalOffsetMODF, modfbuffer, nrow, MPI_DOUBLE, &ioStaMODF); //results are in radiant
		}
		//else, nothing because different node will take care of it...
	}

	if (modfbuffer != NULL) {
		delete[] modfbuffer; modfbuffer = NULL;
	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_File_close(&ioHdlMODF);

	cout << "Worker " << this->get_Rank() << " analyzing MODF = " << (MPI_Wtime() - timer) << " seconds" << endl;
}


void topoHdl::analyze_gsd_db( void )
{
	//####
}


void topoHdl::analyze_drvforce_see_db( void )
{
	//for all grains in a timestep which never touch the boundary compute driving force difference over perimeter
	double timer = MPI_Wtime();

	unsigned int localid = UNKNOWN_ID;
	//all nodes process in parallel, processes do not share time-resolved data
	for (unsigned int f = Settings::SnapshotFirst; f <= Settings::SnapshotLast; f += Settings::SnapshotOffset) {
		localid = WhichOfMyDatasets(f);
		
		if (localid != NOT_ONEOFMINE) {
			string log_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".DrvForceSEE.FID." + std::to_string(f) + ".csv";
			ofstream logsee;
			logsee.open(log_fn.c_str(), std::ofstream::out | std::ofstream::trunc);
			if (logsee.is_open() == true) {
				logsee << "GrainID;SegmentLengthAveragedSEEDrvForce\n";
				logsee << ";Pa\n";

				//process each grain
				for (unsigned int mr = 0; mr < this->LocalDB[localid]->mrg.size(); mr++) {
					unsigned int mrngr = this->LocalDB[localid]->mrg[mr]->GrainBucket.size();
					for (unsigned int g = 0; g < mrngr; ++g) {
						if (this->LocalDB[localid]->mrg[mr]->GrainBucket[g].boundary != BOUNDARY_CONTACT) {
							unsigned int gid = this->LocalDB[localid]->mrg[mr]->GrainBucket[g].extGID;
							struct target_props targetp;
							targetp = this->getTargetProperties( gid, localid, false, false, false, false, true, false );
							logsee << gid << ";" << setprecision(18) << targetp.dsee << "\n";
						}
					}
				}
				logsee.flush();
				logsee.close();
			}
			else {
				cerr << "ERROR::Worker::AnalyzeDrvForceSEE " << this->get_Rank() << " unable to open file " << log_fn.c_str() << endl;
			}			
		}
		//else, nothing because different node will take care of it...
	}

	cout << "...Worker " << this->get_Rank() << " analyzing driving force forward = " << (MPI_Wtime() - timer) << " seconds" << endl;
}


void topoHdl::analyze_see_db( void ) {
	//compute the distribution of stored elastic energy values of sub-grains for the entire domain, including boundary grains!
	//we employ a fixed binning and report a CDF to avoid "overbinning" effects
	//processes work through their files populate a stripe array with the binned data for all time steps in [Settings::SnapshotFirst, Settings::SnapshotLast]
	double timer = MPI_Wtime();

	int ncol = 1 + ((Settings::SnapshotLast - Settings::SnapshotFirst) / Settings::SnapshotOffset + 1); //1 + timesteps because binends added as well
	int nrow = Settings::get_SEEBinCount();
	if ( this->get_Rank() == MASTER ) {
		cout << "...Worker " <<  this->get_Rank() << " Writing rows (SEE bins) = " << nrow << " columns (First col are bin ends, rest are time CDF-data) = " << ncol << endl;
	}

	//create I/O buffer container to write timeslot after timeslot, only as many write ops as timeslots because one node has all grains
	double* seedfbuffer = NULL;
	try { seedfbuffer = new double[nrow]; }
	catch (std::bad_alloc &exc) {
		cerr << "ERROR::Worker " << this->get_Rank() << " allocating memory in analyze_see" << endl; return;
	}
	for ( unsigned int b = 0; b < nrow; ++b ) { seedfbuffer[b] = 0.0; }

	//open MPI files
	MPI_File ioHdlSEEDF;
	MPI_Status ioStaSEEDF;

	string seedf_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".SEEDF.F." + std::to_string(Settings::SnapshotFirst) + ".O." + std::to_string(Settings::SnapshotOffset) + ".L." + std::to_string(Settings::SnapshotLast) + ".NC." + std::to_string(ncol) + ".NR." + std::to_string(nrow) + ".bin";

	// open the file in create and write-only mode, ##only MASTER should finally write the file
	MPI_File_open(MPI_COMM_WORLD, seedf_fn.c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlSEEDF);

	long long totalOffsetSEEDF = 0;
	unsigned int fid = UNKNOWN_ID;
	unsigned int localid = UNKNOWN_ID;

	//MASTER fills in bin end values first
	if ( this->get_Rank() == MASTER ) {
		for ( unsigned int b = 0; b < (nrow-1); ++b ) { //CDFs so we seedfbuffer stores binend values
			seedfbuffer[b] = Settings::SeeBinMin + (((double) b) * Settings::SeeBinWidth);
//cout << "BinID/BinEnd = " << b << "\t\t" << setprecision(16) << seedfbuffer[b] << endl;
		}
		seedfbuffer[nrow-1] = SEE_INFTY;
//cout << "BinID/BinEnd = " << (nrow-1) << "\t\t" << setprecision(16) << seedfbuffer[nrow-1] << endl;

		//write first columns which are the binends
		totalOffsetSEEDF = 0 * nrow * 8;
		MPI_File_write_at( ioHdlSEEDF, totalOffsetSEEDF, seedfbuffer, nrow, MPI_DOUBLE, &ioStaSEEDF); //results are in 1/m^2
	}

	//all nodes process in parallel, processes do not share time-resolved data
	for ( unsigned int f = Settings::SnapshotFirst; f <= Settings::SnapshotLast; f += Settings::SnapshotOffset ) {
		fid = (f / Settings::SnapshotOffset) - ( Settings::SnapshotFirst / Settings::SnapshotOffset );

		localid = WhichOfMyDatasets( f );

//cout << "My local ID is = " << localid << endl;

		if ( localid != NOT_ONEOFMINE ) { //so my responsibility
			
			calc_seedf( localid, seedfbuffer );

			totalOffsetSEEDF = (1 * nrow * 8) + (fid * nrow * 8); //MK::offset nrow*8 :because binends preceed

			MPI_File_write_at( ioHdlSEEDF, totalOffsetSEEDF, seedfbuffer, nrow, MPI_DOUBLE, &ioStaSEEDF); //results are in radiant
		}
		//else, nothing because different node will take care of it...
	}

	if (seedfbuffer != NULL) {
		delete[] seedfbuffer; seedfbuffer = NULL;
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_File_close(&ioHdlSEEDF);

	cout << "...Worker " << this->get_Rank() << " analyzing SEEDF = " << (MPI_Wtime() - timer) << " seconds" << endl;
}


void topoHdl::analyze_grainsize_quantiles_db(void) {
	//compiles cooperatively a table of quantile values eliminating grains at boundaries
	//as the structure coarsens and therefore relatively the more and more grains touch the boundary the counting statistics progressively become weaker
	double timer = MPI_Wtime();
	cout << "...Worker " << this->get_Rank() << " calculating cooperatively quantile values ..." << endl;

	//get memory to store for all my time steps the MPI_Quantiles and the MPI_VOLSTATS
	MPI_VOLSTATS* myvol = NULL;
	try { myvol = new MPI_VOLSTATS[this->LocalDB.size()]; }
	catch (std::bad_alloc &exc) {
		cerr << "ERROR::Worker::AnalyzeGrainSizeQuantiles " << this->get_Rank() << " unable to get memory for storing results in MPI_VOLSTATS" << endl; return;
	}
	double** myquant = NULL;
	try { myquant = new double*[this->LocalDB.size()]; }
	catch (std::bad_alloc &exc) {
		cerr << "ERROR::Worker::AnalyzeGrainSizeQuantiles " << this->get_Rank() << " unable to get memory for storing results in myquant" << endl;
		if (myvol != NULL) { delete[] myvol; myvol = NULL; }
		return;
	}
	bool status = true;
	for (unsigned int s = 0; s < this->LocalDB.size(); s++) {
		myquant[s] = NULL;
		try { myquant[s] = new double[STATISTICS_VOL_NQUANTILES]; }
		catch (std::bad_alloc &exc) {
			cerr << "ERROR::Worker::AnalyzeGrainSizeQuantiles " << this->get_Rank() << " unable to get memory for storing results in myquant" << endl;
			status = false;
		}
		if (status == true) { //fill memory in case memory was allocatable...
			for (unsigned int i = 0; i < STATISTICS_VOL_NQUANTILES; ++i) {
				myquant[s][i] = 0.0;
			}
		}
		else {  //...error handling in case if not
			cerr << "ERROR::Worker::AnalyzeGrainSizeQuantiles " << this->get_Rank() << " unable to get memory for vector quantiles" << endl;
			if (myvol != NULL) { delete[] myvol; myvol = NULL; }
			if (myquant != NULL) {
				for (unsigned int s = 0; s < this->LocalDB.size(); s++) {
					if (myquant[s] != NULL) { delete[] myquant[s]; myquant[s] = NULL; }
				}
				delete[] myquant; myquant = NULL; return;
			}
		}
	}

	//process data
	for (unsigned int s = 0; s < this->LocalDB.size(); s++) {
		this->LocalDB[s]->analyze_vol_quants(&myvol[s], myquant[s]);
		//##MK::Debug
		//for ( unsigned int i = 0; i < STATISTICS_VOL_NQUANTILES; i++ ) cout << "\t\t" << i << ";" << myquant[s][i] << endl;
	}

	//cooperative write to file
	int ncol = ((Settings::SnapshotLast - Settings::SnapshotFirst) / Settings::SnapshotOffset + 1); //##MK:add 1+ to enable in the future the adding of quantile arguments as the first row
	int nrow = STATISTICS_VOL_NQUANTILES; //the same for all worker
	if (this->get_Rank() == MASTER) {
		cout << "...Worker " << this->get_Rank() << " Writing rows (quantiles) = " << nrow << " columns (time) = " << ncol << endl;
	}

	//open two files quants and meta
	MPI_File ioHdlQUANT, ioHdlMETA;
	MPI_Status ioStaQUANT, ioStaMETA;

	string quant_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".VolQuantiles.F." + std::to_string(Settings::SnapshotFirst) + ".O." + std::to_string(Settings::SnapshotOffset) + ".L." + std::to_string(Settings::SnapshotLast) + ".NC." + std::to_string(ncol) + ".NR." + std::to_string(nrow) + ".bin";
	string meta_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".VolMeta.F." + std::to_string(Settings::SnapshotFirst) + ".O." + std::to_string(Settings::SnapshotOffset) + ".L." + std::to_string(Settings::SnapshotLast) + ".NC." + std::to_string(ncol) + ".NR.1.bin";

	// open the file in create and write-only mode, ##only MASTER should finally write the file
	MPI_File_open(MPI_COMM_WORLD, quant_fn.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlQUANT);
	MPI_File_open(MPI_COMM_WORLD, meta_fn.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlMETA);

	long long totalOffsetQuant = 0;
	long long totalOffsetMeta = 0;

	unsigned int fid = UNKNOWN_ID;
	unsigned int localid = UNKNOWN_ID;

	//all nodes process in parallel, processes do not share time-resolved data
	for (unsigned int f = Settings::SnapshotFirst; f <= Settings::SnapshotLast; f += Settings::SnapshotOffset) {
		fid = (f / Settings::SnapshotOffset) - (Settings::SnapshotFirst / Settings::SnapshotOffset);

		localid = WhichOfMyDatasets(f);

		if (localid != NOT_ONEOFMINE) {
			//calculate address in the file //globalid of the dataset * nrow
			totalOffsetQuant = fid * nrow * 8;
			totalOffsetMeta = fid * 1 * (5 * 8);

			//cout << "fid/nrow " << fid << ";" << nrow << ";" << totalOffsetQuant << ";" << totalOffsetMeta << endl;

			MPI_File_write_at(ioHdlQUANT, totalOffsetQuant, myquant[localid], nrow, MPI_DOUBLE, &ioStaQUANT);
			MPI_File_write_at(ioHdlMETA, totalOffsetMeta, &(myvol[localid]), 1, MPI_VOLSTATS_Type, &ioStaMETA);
		}
		//else, nothing because different node will take care of it...
	}

	MPI_Barrier(MPI_COMM_WORLD); //##MK::maybe unnecessary...

	MPI_File_close(&ioHdlQUANT);
	MPI_File_close(&ioHdlMETA);

	//clear memory
	for (unsigned int s = 0; s < this->LocalDB.size(); s++) {
		if (myquant[s] != NULL) {
			delete[] myquant[s]; myquant[s] = NULL;
		}
	}
	if (myquant != NULL) {
		delete[] myquant; myquant = NULL;
	}
	if (myvol != NULL) {
		delete[] myvol; myvol = NULL;
	}

	cout << "...Worker " << this->get_Rank() << " quantile value computation was successful in " << (double)(MPI_Wtime() - timer) << " on rank " << this->get_Rank() << endl;
}


void topoHdl::analyze_classical_nucmodels_db( void ) {
	//for each time step track for all the grains in the network whether they comply with the classical nucleation models of Bailey/Hirsch et. al.
	//keep in mind that when tracking with Offset > 1 grains may apparently disappear from the data because they died in between two sampled integration steps!

	double timer = MPI_Wtime();

	//Bailey et Hirsch in Phil Mag  1960, p833 and follows state:
	//a pre-existent (high-angle) grain boundary segment between two grains A and B with rhoA >= rhoB may bulge when
	//the radius of the (spherical grain) >= 4gamma/(0.5Gb^2*(rhoA*(1.0-f))) with f = rhoB/rhoA
	//Bate and Hutchinson in Script. Mat 36, 2, 195-198, 1997 correct this statement by stating
	//>= 4gamma/(0.5Gb^2*(rhoA*(1.0-f)^0.5))
	//from the context of Bailey's treatment the grain boundary segment has to be a HAGB one
	//Bailey/Hirsch's treatment motivated many people in the SEM/EBSD community to pragmatically identify the necessity for a 
	//"threshold HAGB fraction over the perimeter for a nucleus to be identified as a successful one" and thereby segment an EBSD map 
	//into rxed and deformed grains

	//hence, this function identifies for each grain in the simulation incrementally when it first complies 
	//with either the original or the modified Bailey Hirsch criterion, this requires:
	//1. get targets for the analysis, the array has gid which are consecutive and the same in every rank!
	std::vector<unsigned int>* targets = NULL;
	targets = this->analyze_elimbnd_all_db("BAILEYHIRSCHTargets", NULL, false);
	//every grain can be a potential target as long as it does not touch the boundary, if the sub-grain disappears at some point no problem
	//except for then, thats it with being a successful nucleus =)

	std::vector<unsigned int>* survivors = NULL;
	survivors = analyze_elimbnd_one_db("BAILEYHIRSCHSurvivors", NULL, Settings::SnapshotLast);

	//hanlding of potential inconsistence
	unsigned int localhealth = 0;
	if (targets != NULL && survivors != NULL) localhealth = 2;
	if (nRanks > 1) {
		unsigned int globalhealth = 0;
		MPI_Allreduce(&localhealth, &globalhealth, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
		if (globalhealth != 2 * nRanks) {
			cerr << "ERROR::Worker::AnalyzeBaileyHirsch " << this->get_Rank() << " unable to identify survivors and/or targets for Bailey Hirsch analysis, will stop!" << endl;
			if (targets != NULL) { delete targets; targets = NULL; }
			if (survivors != NULL) { delete survivors; survivors = NULL; }
			return;
		}
	}
	else {
		if (localhealth != 2) {
			cerr << "ERROR::Worker::AnalyzeBaileyHirsch " << this->get_Rank() << " unable to identify survivors and/or targets for Bailey Hirsch analysis, will stop!" << endl;
			if (targets != NULL) { delete targets; targets = NULL; }
			if (survivors != NULL) { delete survivors; survivors = NULL; }
			return;
		}
	}

	//2. create process local hash value that states when for the first time a criterion in the local dataset is met default value is maximum value Settings::SnapshotLast
	if (Settings::LargestGrainID > MAXIMUM_NUMBER_OF_TARGETS) {
		cerr << "Too many targets for the analysis, will stop!" << endl;
		if (targets != NULL) { delete targets; targets = NULL; }
		if (survivors != NULL) { delete survivors; survivors = NULL; }
		return;
	}

	//MK::this hash keeps track of at which integration timestep (if at all) a sub-grain
	//complies with the classical BaileyHirsch criteria and modifications thereof
	unsigned int ngr = 1 + Settings::LargestGrainID;
	unsigned int* NucStatus = NULL;
	try { NucStatus = new unsigned int[CRITERION_METHODS*ngr]; }
	catch (std::bad_alloc &exc) {
		cerr << "ERROR::Worker::AnalyzeBaileyHirsch " << this->get_Rank() << " unable to allocate memory for NucStatus!" << endl;
		if (targets != NULL) { delete targets; targets = NULL; }
		if (survivors != NULL) { delete survivors; survivors = NULL; }
		return;
	}

	//because CRITERION_METHODS different critical radius models to satisfy (1+ to use the array as a hash)
	for (unsigned int hash = 0; hash < ngr; ++hash) { //starts at 0 though it is a hash to overwrite random memory
		NucStatus[CRITERION_METHODS*hash + R_RC_BAIHIR_HAGB] = Settings::SnapshotLast + 1; //setting moment when grain becomes a nucleus beyond Settings::SnapshotLast indicates the grain never became a nucleus, MK::set like this in every node to allow for MPI_Reduce with via MPI_MIN op
		NucStatus[CRITERION_METHODS*hash + R_RC_BAIHIR_ALL] = Settings::SnapshotLast + 1;
		NucStatus[CRITERION_METHODS*hash + R_RC_BATHUT_HAGB] = Settings::SnapshotLast + 1;
		NucStatus[CRITERION_METHODS*hash + R_RC_BATHUT_ALL] = Settings::SnapshotLast + 1;
		NucStatus[CRITERION_METHODS*hash + R50_RC_BAIHIR_HAGB] = Settings::SnapshotLast + 1;
		NucStatus[CRITERION_METHODS*hash + R50_RC_BAIHIR_ALL] = Settings::SnapshotLast + 1;
		NucStatus[CRITERION_METHODS*hash + R50_RC_BATHUT_HAGB] = Settings::SnapshotLast + 1;
		NucStatus[CRITERION_METHODS*hash + R50_RC_BATHUT_ALL] = Settings::SnapshotLast + 1;
	}

	if (myRank == MASTER) {
		cout << "Starting to perform Bailey-Hirsch analysis" << endl;
	}

	unsigned int localid = NOT_YET_KNOWN;
	for (unsigned int f = Settings::SnapshotFirst; f <= Settings::SnapshotLast; f += Settings::SnapshotOffset) {
		localid = WhichOfMyDatasets(f);

		if (localid != NOT_ONEOFMINE) { //so my responsibility

			for (unsigned int t = 0; t < targets->size(); t++) { // analyze for single targets
				unsigned int tgid = targets->at(t);

				struct nucmodel_bh magic;

				magic = analyzeNuc_BaileyHirsch(tgid, localid);

				if (magic.grainfound == true) { //evaluate actual values against criterion values, if critically sized, reset "incubation time" in NucStatus array for sub-grain tgid
					if (magic.r >= magic.rc_baihir_hagb && f <= NucStatus[CRITERION_METHODS*tgid + R_RC_BAIHIR_HAGB])
						NucStatus[CRITERION_METHODS*tgid + R_RC_BAIHIR_HAGB] = f;
					if (magic.r >= magic.rc_baihir_all && f <= NucStatus[CRITERION_METHODS*tgid + R_RC_BAIHIR_ALL])
						NucStatus[CRITERION_METHODS*tgid + R_RC_BAIHIR_ALL] = f;
					if (magic.r >= magic.rc_bathut_hagb && f <= NucStatus[CRITERION_METHODS*tgid + R_RC_BATHUT_HAGB])
						NucStatus[CRITERION_METHODS*tgid + R_RC_BATHUT_HAGB] = f;
					if (magic.r >= magic.rc_bathut_all && f <= NucStatus[CRITERION_METHODS*tgid + R_RC_BATHUT_ALL])
						NucStatus[CRITERION_METHODS*tgid + R_RC_BATHUT_ALL] = f;

					if (magic.r_hagb50 >= magic.rc_baihir_hagb && f <= NucStatus[CRITERION_METHODS*tgid + R50_RC_BAIHIR_HAGB])
						NucStatus[CRITERION_METHODS*tgid + R50_RC_BAIHIR_HAGB] = f;
					if (magic.r_hagb50 >= magic.rc_baihir_all && f <= NucStatus[CRITERION_METHODS*tgid + R50_RC_BAIHIR_ALL])
						NucStatus[CRITERION_METHODS*tgid + R50_RC_BAIHIR_ALL] = f;
					if (magic.r_hagb50 >= magic.rc_bathut_hagb && f <= NucStatus[CRITERION_METHODS*tgid + R50_RC_BATHUT_HAGB])
						NucStatus[CRITERION_METHODS*tgid + R50_RC_BATHUT_HAGB] = f;
					if (magic.r_hagb50 >= magic.rc_bathut_all && f <= NucStatus[CRITERION_METHODS*tgid + R50_RC_BATHUT_ALL])
						NucStatus[CRITERION_METHODS*tgid + R50_RC_BATHUT_ALL] = f;
				}
			} //next target
			  //cout << myRank << " Bailey-Hirsch performed on f = " << f << endl;
		}
		//else, nothing because different node will take care of it...
	} //next timestep

	  //3. MPI_Reduce with MPI_MIN to find for all relevant grains when they nucleated
	  //##MK::here a very conservative parallelization strategy with many barriers is utilized to assure correctness
	for (unsigned int CR = 0; CR < CRITERION_METHODS; CR++) { //MK::criterion-wise necessary because the entire NucStatus array may exceed the capacity of the MPI message transport buffers!
		unsigned int* NucStatusBuffer = NULL;
		try { NucStatusBuffer = new unsigned int[ngr]; }
		catch (std::bad_alloc &exc) {
			cerr << "ERROR::Worker::AnalyzeBaileyHirsch " << this->get_Rank() << " unable to allocate memory for NucStatusBuffer!" << endl;
			if (targets != NULL) { delete targets; targets = NULL; }
			if (survivors != NULL) { delete survivors; survivors = NULL; }
			if (NucStatus != NULL) { delete[] NucStatus; NucStatus = NULL; }
			return;
		}

		//ranks fill in parallel local buffer
		for (unsigned int hash = 0; hash < ngr; ++hash) {
			NucStatusBuffer[hash] = NucStatus[CRITERION_METHODS*hash + CR];
		}

		MPI_Barrier(MPI_COMM_WORLD);

		if (nRanks > 1) {
			if (myRank == MASTER)
				MPI_Reduce(MPI_IN_PLACE, NucStatusBuffer, ngr, MPI_UNSIGNED, MPI_MIN, MASTER, MPI_COMM_WORLD);
			else
				MPI_Reduce(NucStatusBuffer, NULL, ngr, MPI_UNSIGNED, MPI_MIN, MASTER, MPI_COMM_WORLD);
		}
		else {
			//master has already all pieces of information
		}

		if (myRank == MASTER) {
			for (unsigned int hash = 0; hash <= Settings::LargestGrainID; ++hash) {
				NucStatus[CRITERION_METHODS*hash + CR] = NucStatusBuffer[hash];
			}
			//cout << "Master aggregate criterion data for the " << CR << "-th criterion" << endl;
		}

		MPI_Barrier(MPI_COMM_WORLD);

		if (NucStatusBuffer != NULL) { delete[] NucStatusBuffer; NucStatusBuffer = NULL; }
	} //aggregate results for next criterion

	  //cout << "Worker " << myRank << " aggregated all criterion data" << endl;

	  //4. write log file
	if (myRank == MASTER) {
		double ftimer = MPI_Wtime();
		string log_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".BaileyHirschAnalysis.F." + std::to_string(Settings::SnapshotFirst) + ".O." + std::to_string(Settings::SnapshotOffset) + ".L." + std::to_string(Settings::SnapshotLast) + ".csv";
		ofstream loglast;
		loglast.open(log_fn.c_str());
		if (loglast.is_open() == true) {
			loglast << "GrainID;R_RC_BAIHIR_HAGB;R_RC_BAIHIR_ALL;R_RC_BATHUT_HAGB;R_RC_BATHUT_ALL;R50_RC_BAIHIR_HAGB;R50_RC_BAIHIR_ALL;R50_RC_BATHUT_HAGB;R50_RC_BATHUT_ALL;Survivor?\n";

			stringstream line;
			for (unsigned int t = 0; t < targets->size(); t++) {
				unsigned int gid = targets->at(t);
				//included in survivors?
				bool asurvivor = false;
				for (unsigned int s = 0; s < survivors->size(); s++) {
					if (gid != survivors->at(s))
						continue;
					//not continued found
					asurvivor = true; break;
				}

				unsigned int entry = CRITERION_METHODS*targets->at(t);

				//speculative caching
				line.str(std::string());
				line << targets->at(t) << ";" << NucStatus[entry + R_RC_BAIHIR_HAGB] << ";" << NucStatus[entry + R_RC_BAIHIR_ALL] << ";" << NucStatus[entry + R_RC_BATHUT_HAGB] << ";" << NucStatus[entry + R_RC_BATHUT_ALL] << ";";
				line << NucStatus[entry + R50_RC_BAIHIR_HAGB] << ";" << NucStatus[entry + R50_RC_BAIHIR_ALL] << ";" << NucStatus[entry + R50_RC_BATHUT_HAGB] << ";" << NucStatus[entry + R50_RC_BATHUT_ALL] << ";" << (int)asurvivor;

				if (NucStatus[entry + R_RC_BAIHIR_HAGB] < (Settings::SnapshotLast + 1)) { loglast << line.str() << endl; continue; }
				if (NucStatus[entry + R_RC_BAIHIR_ALL] < (Settings::SnapshotLast + 1)) { loglast << line.str() << endl; continue; }
				if (NucStatus[entry + R_RC_BATHUT_HAGB] < (Settings::SnapshotLast + 1)) { loglast << line.str() << endl; continue; }
				if (NucStatus[entry + R_RC_BAIHIR_ALL] < (Settings::SnapshotLast + 1)) { loglast << line.str() << endl; continue; }
				if (NucStatus[entry + R50_RC_BAIHIR_HAGB] < (Settings::SnapshotLast + 1)) { loglast << line.str() << endl; continue; }
				if (NucStatus[entry + R50_RC_BAIHIR_ALL] < (Settings::SnapshotLast + 1)) { loglast << line.str() << endl; continue; }
				if (NucStatus[entry + R50_RC_BATHUT_HAGB] < (Settings::SnapshotLast + 1)) { loglast << line.str() << endl; continue; }
				if (NucStatus[entry + R50_RC_BATHUT_ALL] < (Settings::SnapshotLast + 1)) { loglast << line.str() << endl; continue; }
			}

			loglast.flush();
			loglast.close();
		}

		cout << "Worker::AnalyzeBaileyHirsch " << this->get_Rank() << " wrote BaileyHirsch to file in " << (MPI_Wtime() - ftimer) << " seconds" << endl;
	}


	if (targets != NULL) { delete targets; targets = NULL; }
	if (survivors != NULL) { delete survivors; survivors = NULL; }
	if (NucStatus != NULL) { delete[] NucStatus; NucStatus = NULL; }

	cout << "Worker::AnalyzeBaileyHirsch " << this->get_Rank() << " analyzed Bailey Hirsch criterion in " << (MPI_Wtime() - timer) << " seconds" << endl;
}


unsigned int topoHdl::IsIncluded( unsigned int ggid, std::vector<unsigned int>* cand ) {
	for ( unsigned int c = 0; c < cand->size(); c++ ) {
		if ( ggid != cand->at(c) ) {
			continue;
		}
		//included!
		return 1;
	}

	//not found, so not included!
	return 0;
}


struct aggstats topoHdl::inspectForAbnormalGrains( unsigned int lsid, std::vector<unsigned int>* surv, std::vector<unsigned int>* matr, std::vector<unsigned int>* targ )
{
	//computes the instantaneous domain area covered by the three grain populations
	struct aggstats res;

//analyze for the survivors
	unsigned int gid = 0;
	struct descr_stats surv_stats;
	for ( unsigned int s = 0; s < surv->size(); s++ ) {
		gid = surv->at(s);
		unsigned int luid = gid / Settings::LookupMaxGrainIDRange;
		unsigned int lu_size = this->LocalDB[lsid]->LookupTable[luid].size;
		GrainHash* thebucket = this->LocalDB[lsid]->LookupTable[luid].bucket;
		unsigned int whichmr = NOT_ASSIGNED_YET;
		unsigned int whichidx = NOT_ASSIGNED_YET;

		for ( unsigned int g = 0; g < lu_size; g++ ) {
			if ( thebucket[g].gid != gid ) 
				continue;
			//not continued so found
			whichmr = thebucket[g].mrg_idx;
			whichidx = thebucket[g].idx;
			break;
		}
		if (whichmr == NOT_ASSIGNED_YET || whichidx == NOT_ASSIGNED_YET) {
			cerr << "ERROR::Worker::AnalyzeHolmRollett " << this->get_Rank() << " inspecting for AGG survivor grain " << gid << " not found!" << endl; 
			continue;
		}
		
		surv_stats.n++;
		surv_stats.sum += this->LocalDB[lsid]->mrg[whichmr]->GrainBucket[whichidx].size;
	}

	if ( surv_stats.n > 0 ) 
		surv_stats.mean = surv_stats.sum / ((double) surv_stats.n);	//##MK::no student t distribution...


//analyze for the matrix
	struct descr_stats matr_stats;
	for ( unsigned int m = 0; m < matr->size(); m++ ) {
		gid = matr->at(m);
		unsigned int luid = gid / Settings::LookupMaxGrainIDRange;
		unsigned int lu_size = this->LocalDB[lsid]->LookupTable[luid].size;
		GrainHash* thebucket = this->LocalDB[lsid]->LookupTable[luid].bucket;
		unsigned int whichmr = NOT_ASSIGNED_YET;
		unsigned int whichidx = NOT_ASSIGNED_YET;

		for ( unsigned int g = 0; g < lu_size; g++ ) {
			if ( thebucket[g].gid != gid ) continue;

			whichmr = thebucket[g].mrg_idx;
			whichidx = thebucket[g].idx;
			break;
		}
		if ( whichmr == NOT_ASSIGNED_YET || whichidx == NOT_ASSIGNED_YET ) { 
			//cerr << "ERROR::Worker::AnalyzeHolmRollett " << this->get_Rank() << " inspecting for AGG matrix grain " << gid << " not found!" << endl;
			continue;
		}
		matr_stats.n++;
		matr_stats.sum += this->LocalDB[lsid]->mrg[whichmr]->GrainBucket[whichidx].size;
	}

	if ( matr_stats.n > 0 ) 
		matr_stats.mean = matr_stats.sum / ((double) matr_stats.n);

//analyze for targets
	struct descr_stats targ_stats;
	for ( unsigned int t = 0; t < targ->size(); t++ ) {
		gid = targ->at(t);
		unsigned int luid = gid / Settings::LookupMaxGrainIDRange;
		unsigned int lu_size = this->LocalDB[lsid]->LookupTable[luid].size;
		GrainHash* thebucket = this->LocalDB[lsid]->LookupTable[luid].bucket;
		unsigned int whichmr = NOT_ASSIGNED_YET;
		unsigned int whichidx = NOT_ASSIGNED_YET;
		
		for ( unsigned int g = 0; g < lu_size; g++ ) {
			if ( thebucket[g].gid != gid ) continue;

			whichmr = thebucket[g].mrg_idx;
			whichidx = thebucket[g].idx;
			break;
		}
		if ( whichmr == NOT_ASSIGNED_YET || whichidx == NOT_ASSIGNED_YET ) { 
			//cerr << "ERROR::Worker::AnalyzeHolmRollett " << this->get_Rank() << " inspecting for AGG target grain " << gid << " not found!" << endl;
			continue;
		}

		targ_stats.n++;
		targ_stats.sum += this->LocalDB[lsid]->mrg[whichmr]->GrainBucket[whichidx].size;
	}

	if ( targ_stats.n > 0 ) 
		targ_stats.mean = targ_stats.sum / ((double) targ_stats.n);

	//interpret the generic descriptive stats
	res.totalsize_survivors = surv_stats.sum; //normalized area/volume on [0,1]^Settings::Dimensionality unit domain
	res.totalsize_matrix = matr_stats.sum;
	res.totalsize_targets = targ_stats.sum;

	if ( Settings::Dimensionality == TWODIMENSIONAL ) { //EHRs
		res.rmean_survivors = pow( (surv_stats.mean*SQR(METER2MICRON(Settings::PhysicalDomainSize)))/PI, 0.5 ); // / METER2MICRON))/PI , 0.5 ); //radii in micron
		res.rmean_matrix = pow( (matr_stats.mean*SQR(METER2MICRON(Settings::PhysicalDomainSize)))/PI , 0.5 );
		res.rmean_targets = pow( (targ_stats.mean*SQR(METER2MICRON(Settings::PhysicalDomainSize)))/PI , 0.5 );
	}
	if ( Settings::Dimensionality == THREEDIMENSIONAL ) {
		res.rmean_survivors = pow( (surv_stats.mean*CUBE(METER2MICRON(Settings::PhysicalDomainSize)))/((4.0/3.0)*PI) , (1.0/3.0) ); //radii in micron
		res.rmean_matrix = pow( (matr_stats.mean*CUBE(METER2MICRON(Settings::PhysicalDomainSize)))/((4.0/3.0)*PI) , (1.0/3.0) );
		res.rmean_targets = pow( (targ_stats.mean*CUBE(METER2MICRON(Settings::PhysicalDomainSize)))/((4.0/3.0)*PI) , (1.0/3.0) );
	}

	res.n_survivors = surv_stats.n;
	res.n_matrix = matr_stats.n;
	res.n_targets = targ_stats.n;

//determine abnormal grains from the list of targets
	if ( matr_stats.mean > DBL_EPSILON && targ_stats.mean > DBL_EPSILON ) {
		for ( unsigned int t = 0; t < targ->size(); t++ ) {
			gid = targ->at(t);
			unsigned int luid = gid / Settings::LookupMaxGrainIDRange;
			unsigned int lu_size = this->LocalDB[lsid]->LookupTable[luid].size;
			GrainHash* thebucket = this->LocalDB[lsid]->LookupTable[luid].bucket;
			unsigned int whichmr = NOT_ASSIGNED_YET;
			unsigned int whichidx = NOT_ASSIGNED_YET;
		
			for ( unsigned int g = 0; g < lu_size; g++ ) {
				if ( thebucket[g].gid != gid ) continue;

				whichmr = thebucket[g].mrg_idx;
				whichidx = thebucket[g].idx;
				break;
			}
			if ( whichmr == NOT_ASSIGNED_YET || whichidx == NOT_ASSIGNED_YET ) { 
				//cerr << "ERROR::Worker::AnalyzeHolmRollett " << this->get_Rank() << " inspecting for AGG target grain " << gid << " not found!" << endl;
				continue;
			}

			//calc EHR for this target grain
			double tsz = this->LocalDB[lsid]->mrg[whichmr]->GrainBucket[whichidx].size;
			if ( Settings::Dimensionality == TWODIMENSIONAL ) 
				tsz = pow( (tsz*SQR(METER2MICRON(Settings::PhysicalDomainSize)))/PI , 0.5 );
			if ( Settings::Dimensionality == THREEDIMENSIONAL )
				tsz = pow( (tsz*CUBE(METER2MICRON(Settings::PhysicalDomainSize)))/((4.0/3.0)*PI) , (1.0/3.0) );

			//considerable as abnormal with respect to matrix grains?
			if ( tsz / res.rmean_matrix >= CRITICAL_RADIUS_RATIO ) {
				res.n_agg_m++;
				res.n_agg_m_insurv = res.n_agg_m_insurv + IsIncluded( gid, surv ); //##MK::interal optimization of IsIncluded possible by hash
			}
			if ( tsz / res.rmean_targets >= CRITICAL_RADIUS_RATIO ) {
				res.n_agg_t++;
				res.n_agg_t_insurv = res.n_agg_t_insurv + IsIncluded( gid, surv );
			}
		} //next target
	}

	return res;
}


void topoHdl::analyze_abnormalgraingrowth_db( void ) {
	//Holm, Miodownik and Rollett in Acta Mat 2003, 51, 2701 analyzed abnormal sub-grain coarsening
	double timer = MPI_Wtime();
	//this function computes the instantaneous descriptive means of the matrix grains and counts how many grains are considered as to grow at the moment abnormally
	//we define three lists of grain ids
	//survivors - all the grains without boundary contact in the final data set, most likely AGG grains and remaining matrix
	//targets - all grains ever considered with never boundary contact
	//matrix - the targets but excluding the survivors

	std::vector<unsigned int>* survivors = NULL;
	survivors = analyze_elimbnd_one_db("AGGHolmRollSurvivors", NULL, Settings::SnapshotLast);
	if (survivors == NULL) {
		cerr << "ERROR::Worker::AnalyzeHolmRollett " << this->get_Rank() << " unable to determine survivors!" << endl;
		return;
	}

	//cout << myRank << "-th rank has " << survivors->size() << " survivors determined!" << endl;

	std::vector<unsigned int>* targets = NULL;
	targets = analyze_elimbnd_all_db("AGGHolmRollTargets", NULL, false);
	if (targets == NULL) {
		cerr << "ERROR::Worker::AnalyzeHolmRollett " << this->get_Rank() << " unable to determine targets!" << endl;
		if (survivors != NULL) { delete survivors; survivors = NULL; }
		return;
	}

	//cout << myRank << "-th rank has " << targets->size() << " targets determined!" << endl;

	std::vector<unsigned int>* matrix = NULL;
	matrix = analyze_elimbnd_all_db("AGGHolmRollMatrix", survivors, false);
	if (matrix == NULL) {
		cerr << "ERROR::Worker::AnalyzeHolmRollett " << this->get_Rank() << " unable to determine matrix" << endl;
		if (survivors != NULL) { delete survivors; survivors = NULL; }
		if (targets != NULL) { delete targets; targets = NULL; }
		return;
	}

	//cout << myRank << "-th rank has " << matrix->size() << " matrix determined!" << endl;

	//now all know the matrixgrains and the rxgrains so they can compute rxgrains
	int ncol = ((Settings::SnapshotLast - Settings::SnapshotFirst) / Settings::SnapshotOffset + 1);
	int nrow = 1;
	cout << "...Worker " << this->get_Rank() << " writing rows (only one-dimension of AGG_STATS) = " << nrow << " columns (time) = " << ncol << endl;

	//create I/O buffer container to write timeslot after timeslot, only as many write ops as timeslots because one node has all grains
	MPI_AGGSTATS* aggbuffer = NULL;
	try { aggbuffer = new MPI_AGGSTATS[1]; }
	catch (std::bad_alloc &exc) {
		cerr << "ERROR::Worker::AnalyzeHolmRollett " << this->get_Rank() << " unable to allocate result memory" << endl;
		if (survivors != NULL) { delete survivors; survivors = NULL; }
		if (targets != NULL) { delete targets; targets = NULL; }
		if (matrix != NULL) { delete matrix; matrix = NULL; }
		return;
	}

	MPI_File ioHdlAGG;
	MPI_Status ioStaAGG;

	string agg_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".AGGHolmRoll.F." + std::to_string(Settings::SnapshotFirst) + ".O." + std::to_string(Settings::SnapshotOffset) + ".L." + std::to_string(Settings::SnapshotLast) + ".NC." + std::to_string(ncol) + ".NR." + std::to_string(nrow) + ".bin";
	MPI_File_open(MPI_COMM_WORLD, agg_fn.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlAGG);

	long long totalOffsetAGG = 0;
	unsigned int fid = UNKNOWN_ID;
	unsigned int localid = UNKNOWN_ID;

	//all nodes process in parallel, processes do not share time-resolved data
	for (unsigned int f = Settings::SnapshotFirst; f <= Settings::SnapshotLast; f += Settings::SnapshotOffset) {
		fid = (f / Settings::SnapshotOffset) - (Settings::SnapshotFirst / Settings::SnapshotOffset);

		localid = WhichOfMyDatasets(f);

		if (localid != NOT_ONEOFMINE) {
			struct aggstats agr;

			agr = inspectForAbnormalGrains(localid, survivors, matrix, targets);

			cout << agr.totalsize_survivors << ";" << agr.totalsize_matrix << ";" << agr.totalsize_targets << ";" << agr.rmean_survivors << ";" << agr.rmean_matrix << ";" << agr.rmean_targets << agr.n_survivors << ";" << agr.n_matrix << ";" << agr.n_targets << endl;

			aggbuffer[0].totalsize_survivors = agr.totalsize_survivors;
			aggbuffer[0].totalsize_matrix = agr.totalsize_matrix;
			aggbuffer[0].totalsize_targets = agr.totalsize_targets;

			aggbuffer[0].rmean_survivors = agr.rmean_survivors;
			aggbuffer[0].rmean_matrix = agr.rmean_matrix;
			aggbuffer[0].rmean_targets = agr.rmean_targets;

			aggbuffer[0].n_survivors = agr.n_survivors;
			aggbuffer[0].n_matrix = agr.n_matrix;
			aggbuffer[0].n_targets = agr.n_targets;
			aggbuffer[0].n_agg_m = agr.n_agg_m;
			aggbuffer[0].n_agg_t = agr.n_agg_t;
			aggbuffer[0].n_agg_m_insurv = agr.n_agg_m_insurv;
			aggbuffer[0].n_agg_t_insurv = agr.n_agg_t_insurv;

			//calculate address in the file //globalid of the dataset * nrow
			totalOffsetAGG = fid * nrow * ((3 + 3 + 3 + 4) * 8);

			MPI_File_write_at(ioHdlAGG, totalOffsetAGG, aggbuffer, 1, MPI_AGGSTATS_Type, &ioStaAGG);
		}
		//else, nothing because different node will take care of it...
	}

	if (aggbuffer != NULL) {
		delete[] aggbuffer; aggbuffer = NULL;
	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_File_close(&ioHdlAGG);

	if (survivors != NULL) { delete survivors; survivors = NULL; }
	if (targets != NULL) { delete targets; targets = NULL; }
	if (matrix != NULL) { delete matrix; matrix = NULL; }

	cout << "Worker " << this->get_Rank() << " analyzed AGG Holm Rollett = " << (MPI_Wtime() - timer) << " seconds" << endl;
}


void topoHdl::analyze_rxfraction_db( unsigned int extid ) {
	//identifies all remaining grains that never touched the boundary
	//##MK::here we assume that all grains that in time step extid do not touch the boundary also before never touched the boundary
	//computes thereafter the area for these grains in relation to the total area for each timestep of all grains which do not touch the boundary
	double timer = MPI_Wtime();

	std::vector<unsigned int>* matrixgrains = NULL;
	matrixgrains = this->analyze_elimbnd_all_db( "RXEVO", NULL, false ); //everybody has to call this, dont forget to delete the pointer!

	if ( matrixgrains != NULL ) {
		if ( matrixgrains->size() < 1 ) { 
			cerr << "ERROR::Worker::AnalyzeApproxRX " << this->get_Rank() << " analyze approximate recrystallized fraction impossible, no grain survived" << endl;
			if (matrixgrains != NULL) { delete matrixgrains; matrixgrains = NULL; }
			return;
		}
	}
	else {
		cerr << "ERROR::Worker::AnalyzeApproxRX " << this->get_Rank() << " analyze approximate recrystallized fraction impossible, impossible to identify survivors" << endl;
		return;
	}

	std::vector<unsigned int>* rxgrains = NULL;
	rxgrains = this->analyze_elimbnd_one_db("RXEVO", matrixgrains, extid); //##MK::better work with bitfield to check inclusion in whitelist

	if (rxgrains != NULL) {
		if (rxgrains->size() < 1) {
			cerr << "ERROR::Worker:AnalyzeApproxRX " << this->get_Rank() << " analyze approximate recrystallized fraction impossible, no grain survived" << endl;
			if (matrixgrains != NULL) { delete matrixgrains; matrixgrains = NULL; }
			if (rxgrains != NULL) { delete rxgrains; rxgrains = NULL; }
			return;
		}
	}
	else {
		cerr << "ERROR::Worker:AnalyzeApproxRX " << this->get_Rank() << " analyze approximate recrystallized fraction impossible, impossible to identify survivors" << endl;
		if (matrixgrains != NULL) { delete matrixgrains; matrixgrains = NULL; }
		return;
	}

	//now all know the matrixgrains and the rxgrains so they can compute rxgrains
	int ncol = ((Settings::SnapshotLast - Settings::SnapshotFirst) / Settings::SnapshotOffset + 1);
	int nrow = 1;
	if (this->get_Rank() == MASTER) {
		cout << "...Worker " << this->get_Rank() << " writing rows (only one-dimension of RX_STATS) = " << nrow << " columns (time) = " << ncol << endl;
	}

	//create I/O buffer container to write timeslot after timeslot, only as many write ops as timeslots because one node has all grains
	MPI_RXSTATS* rxbuffer = NULL;
	try { rxbuffer = new MPI_RXSTATS[1]; }
	catch (std::bad_alloc &exc) {
		cerr << "ERROR::Worker::AnalyzeApproxRX " << this->get_Rank() << " unable to allocate memory for MPI_RXSTATS" << endl;
		if (matrixgrains != NULL) { delete matrixgrains; matrixgrains = NULL; }
		if (rxgrains != NULL) { delete rxgrains; rxgrains = NULL; }
		return;
	}

	MPI_File ioHdlRX;
	MPI_Status ioStaRX;

	string rx_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".RXEVO.F." + std::to_string(Settings::SnapshotFirst) + ".O." + std::to_string(Settings::SnapshotOffset) + ".L." + std::to_string(Settings::SnapshotLast) + ".NC." + std::to_string(ncol) + ".NR.1.bin";

	MPI_File_open(MPI_COMM_WORLD, rx_fn.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlRX);

	long long totalOffsetRX = 0;
	unsigned int fid = UNKNOWN_ID;
	unsigned int localid = UNKNOWN_ID;

	//all nodes process in parallel, processes do not share time-resolved data
	for (unsigned int f = Settings::SnapshotFirst; f <= extid; f += Settings::SnapshotOffset) {
		fid = (f / Settings::SnapshotOffset) - (Settings::SnapshotFirst / Settings::SnapshotOffset);

		localid = WhichOfMyDatasets(f);

		if (localid != NOT_ONEOFMINE) {
			struct rxstats rxs;
			rxs = calcRXFraction(localid, matrixgrains, rxgrains);

			//cout << rxs.totalsize_elimbnd << ";" << rxs.totalsize_targets << ";" << rxs.X << ";" << rxs.ndefg << ";" << rxs.nrxg << endl;

			rxbuffer[0].totalsize_elimbnd = rxs.totalsize_elimbnd;
			rxbuffer[0].totalsize_targets = rxs.totalsize_targets;
			rxbuffer[0].X = rxs.X;
			rxbuffer[0].ndefg = rxs.ndefg;
			rxbuffer[0].nrxg = rxs.nrxg;

			//calculate address in the file //globalid of the dataset * nrow
			totalOffsetRX = fid * nrow * (5 * 8);

			MPI_File_write_at(ioHdlRX, totalOffsetRX, rxbuffer, 1, MPI_RXSTATS_Type, &ioStaRX);
		}
		//else, nothing because different node will take care of it...
	}

	if (rxbuffer != NULL) { delete[] rxbuffer; rxbuffer = NULL; }

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_File_close(&ioHdlRX);

	if (matrixgrains != NULL) { delete matrixgrains; matrixgrains = NULL; }
	if (rxgrains != NULL) { delete rxgrains; rxgrains = NULL; }

	cout << "...Worker " << this->get_Rank() << " analyzed recrystallized fraction in " << (MPI_Wtime() - timer) << " seconds" << endl;
}


void topoHdl::spit_profiling( string mode )
{
	//MYProfiler
	//first overview of profiling for single process execution...
	string mylogfname;
	mylogfname = "TopologyTracer." + mode + ".SimID." + std::to_string(Settings::SimID) + ".Rank." + std::to_string(this->get_Rank()) + ".MyProfiling.csv";
	ofstream mylogfile;
	mylogfile.open(mylogfname.c_str(), std::ofstream::out | std::ofstream::trunc);
	if (mylogfile.is_open() == true) {
		//header
		mylogfile << "What;WallClock\n";
		mylogfile << ";s\n";
		mylogfile << ";MPI_Wtime\n";

		for (unsigned int l = 0; l < myprofiler.get_nentries(); l++) {
			mylogfile << myprofiler.titles[l] << ";" << setprecision(8) << myprofiler.times[l] << "\n";
		}
		mylogfile.flush();
		mylogfile.close();
	}
	else {
		cerr << "ERROR::Worker::TopoHdl " << this->get_Rank() << " unable to write myprofiling" << endl;
	}

	//IOProfiler
	//first overview of profiling for single process execution...
	string iologfname;
	iologfname = "TopologyTracer.TrackParallel.SimID." + std::to_string(Settings::SimID) + ".Rank." + std::to_string(this->get_Rank()) + ".IoProfiling.csv";
	ofstream iologfile;
	iologfile.open(iologfname.c_str(), std::ofstream::out | std::ofstream::trunc);
	if (iologfile.is_open() == true) {
		//header
		iologfile << "What;WallClock\n";
		iologfile << ";s\n";
		iologfile << ";MPI_Wtime\n";

		for (unsigned int l = 0; l < ioprofiler.get_nentries(); l++) {
			iologfile << ioprofiler.titles[l] << ";" << setprecision(8) << ioprofiler.times[l] << "\n";
		}
		iologfile.flush();
		iologfile.close();
	}
	else {
		cerr << "ERROR::Worker::TopoHdl " << this->get_Rank() << " unable to write i/o-profiling" << endl;
	}
}


bool topoHdl::readSingleDataset(unsigned int fd)
{
	bool filegood = false;
	vector<SnapshotMetaData>::iterator s;
	for (s = DatasetInfo.begin(); s != DatasetInfo.end(); ++s) { //all ranks know all
		if (s->extID != fd)
			continue;
		//not continued so found
		filegood = true;
		break;
	}

	if (filegood == true) {
		filegood = readSingleDataset(s, true, false);
	}
	return filegood;			
}



bool topoHdl::readSingleDataset( std::vector<SnapshotMetaData>::iterator si, bool loadgrains, bool loadfaces )
{
	bool status = true;
	if (loadgrains == true && loadfaces == true) {
		//##MK::full load for processing data in TrackingSequentially
		//#######################
		//#######################
		return status;
	}

	if (loadgrains == true && loadfaces == false) { //only for probing for instance bnd contact
		double t1 = MPI_Wtime();
		if (this->grbuf != NULL) { delete[] this->grbuf; this->grbuf = NULL; this->grbufsize = 0; }

		//single load of Texture_si->extID.bin without writing into database...
		if (Settings::Dimensionality == THREEDIMENSIONAL) {
			status = this->read_mpiio_grains_binary_3d(si);
		}
		else if (Settings::Dimensionality == TWODIMENSIONAL) { 
			status = this->read_mpiio_grains_binary_2d(si);
		}
		else { return false; } //cout << "ERROR::Worker " << this->get_Rank() << " invalid dimensionality during readSingleDataset!" << endl; 
			
		double t2 = MPI_Wtime();
		ioprofiler.logev("MPIIOTexture" + std::to_string(si->extID), (t2 - t1));
		cout << "...Worker " << this->get_Rank() << " loaded Texture_" << si->extID << ".bin in " << (t2 - t1) << " seconds" << endl;
		return status;
	}

	//##MK::all other cases not implemented
	return false;
}


bool topoHdl::analyze_rvebnd_contacttime(void)
{
	double timer = MPI_Wtime();
	//probe all files for existence with all workers
	if (queryFileSizes() == false) {
		cerr << "ERROR::Worker::TrackingSequentiallyRVEContact " << this->get_Rank() << " dataset not complete or otherwise I/O faulty!" << endl;
		return false;
	}
	//not returned, then DatasetInfo is populated in all ranks, so we start processing

	if (Settings::ProbeBoundaryContact == false) { //if we know already via a file which grains made boundary contact, we proceed and load the file directly

		if ( load_rvebnd_contacttime() == false )
			return false;
		else
			return true; //now we know the boundary contacts
	}
	//else hence not returned so we have to generate first such a file by loading in parallel the Texture_<fid>.bin heavy data and process the contacts for all grains...
	unsigned int ngr = 1 + Settings::LargestGrainID;
	
	unsigned int* when = NULL;
	try { when = new unsigned int[ngr]; }
	catch (std::bad_alloc &exc) {
		cerr << "ERROR::Worker::TrackingSequentiallyRVEContact " << this->get_Rank() << " unable to allocate results array!" << endl;
		return false;
	}
	
	//fill local container of potential boundary contact time
	unsigned int NeverDuringExistence = std::numeric_limits<unsigned int>::max();
	for (unsigned int gid = 0; gid < ngr; ++gid) { when[gid] = NeverDuringExistence; }

	//identify boundary existence with the heavy data
	unsigned int localgidmax = 0;
	unsigned int localhealth = 1;
	for (unsigned int fid = Settings::SnapshotFirst; fid <= Settings::SnapshotLast; fid += Settings::SnapshotOffset) {
		//round robin work partitioning
		unsigned int f = 1 + fid / Settings::SnapshotOffset - Settings::SnapshotFirst / Settings::SnapshotOffset;
		if (f % nRanks == myRank) {
			vector<SnapshotMetaData>::iterator s;
			for (s = DatasetInfo.begin(); s != DatasetInfo.end(); ++s) { //all ranks know all
				if (s->extID != fid)
					continue;
				//not continued so found
				break;
			}
			
			bool filegood = readSingleDataset( s, true, false );
			
			if (filegood == true) { //process through
				unsigned int gid = 0;
				for (unsigned int g = 0; g < this->grbufsize; ++g) {
					gid = this->grbuf[g].id;
					//check for potential boundary contact
					if (gid != THE_DOMAIN_ITSELF && this->grbuf[g].intersectsBoundaryGrain != BOUNDARY_CONTACT) { //no change in NeverDuringExistence necessary because grain does not make boundary contact
						continue;
					}
					//not continued, so either zero grain or boundary contact, most likely within ID user-defined maximum ID range
					if (gid < ngr) {
						if (fid <= when[gid]) {
							when[gid] = fid;
						}
					}
					else { //if not user input and data set inconsistent!
						cerr << "ERROR::Worker::TrackingSequentiallyRVEContact " << this->get_Rank() << " found grain ID larger than expected!" << endl;
						localhealth = 0;
						break; //no longer analyzing any						
					}
					//determine maximum GrainID in dataset
					if (gid > localgidmax) { localgidmax = gid; }
				}
			}
			else {
				cerr << "ERROR::Worker::TrackingSequentiallyRVEContact " << this->get_Rank() << " dataset " << "Texture_" << fid << ".bin non-existent!" << endl;
				localhealth = 0;
				break;			
			}			
		} //probe next workload for me...
		//else other one will take care of it
	}

	unsigned int globalhealth = 0;
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Allreduce(&localhealth, &globalhealth, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
	if (globalhealth != nRanks) {
		cerr << "ERROR::Worker::TrackingSequentiallyRVEContact " << this->get_Rank() << " either a dataset is missing or we found a grain ID larger as than expected!" << endl;
		return false;
	}

	unsigned int globalgidmax = 0;
	MPI_Allreduce(&localgidmax, &globalgidmax, 1, MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);
	cout << "...Worker " << this->get_Rank() << " found maximum ID to be " << globalgidmax << endl;
	
	//aggregate contact times data on master
	if (nRanks > SINGLE_PROCESS) {
		if (myRank == MASTER) {
			MPI_Reduce(MPI_IN_PLACE, when, ngr, MPI_UNSIGNED, MPI_MIN, MASTER, MPI_COMM_WORLD);
		}
		else {
			MPI_Reduce(when, NULL, ngr, MPI_UNSIGNED, MPI_MIN, MASTER, MPI_COMM_WORLD);
		}
	}
	//else no aggregation not necessary because master has data already "inplace" a hash in ascending ID order


	//output results for future use
	if (myRank == MASTER) {
		//clear zero grain
		when[THE_DOMAIN_ITSELF] = 0;
		//and clear all grains with IDs > globalgidmax
		for (unsigned int gid = globalgidmax + 1; gid < ngr; gid++) {
			when[gid] = 0;
		}
		//MK::otherwise we would later bookkeep non-existent IDs in filter lists

		string log_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".GrainBoundaryContact.F." + std::to_string(Settings::SnapshotFirst) + ".O." + std::to_string(Settings::SnapshotOffset) + ".L." + std::to_string(Settings::SnapshotLast) + ".csv";
		ofstream logbnd;
		logbnd.open(log_fn.c_str(), std::ofstream::out | std::ofstream::trunc);
		if (logbnd.is_open() == true) {
			//... with a header
			logbnd << "GrainID;WhenMakingContact(Timestep)\n";
			for (unsigned int gid = 0; gid < ngr; ++gid) {
				logbnd << gid << ";" << when[gid] << "\n";
			}
			logbnd.flush();
			logbnd.close();
		}
		else {
			cerr << "ERROR::Worker::TrackingSequentiallyRVEContact " << this->get_Rank() << " unable to output to results file!" << endl;
			return false;
		}
	}

	//transfer data to this->bndbuf
	if (this->bndbuf != NULL) { delete[] this->bndbuf; this->bndbuf = NULL; this->bndbufsize = 0; }
	try { this->bndbuf = new MPI_GrainBndContactIO[ngr]; }
	catch (std::bad_alloc &exc) {
		cout << "...Worker " << this->get_Rank() << " analyze boundary contact in " << (MPI_Wtime() - timer) << " seconds BUT was unable to allocate bndbuf!" << endl;
		cerr << "ERROR::Worker::AnalyzeRVEBndContacttime" << this->get_Rank() << " bndbuf unallocatable!" << endl;
		if (when != NULL) { delete[] when; when = NULL; }
		return false;
	}
	this->bndbufsize = ngr;
	for (unsigned int gid = 0; gid < ngr; ++gid) {
		this->bndbuf[gid].extGID = gid;
		this->bndbuf[gid].whenhitboundary = when[gid];
	}

	if (when != NULL) { delete[] when; when = NULL; }
	cout << "...Worker " << this->get_Rank() << " analyze boundary contact in " << (MPI_Wtime() - timer) << " seconds" << endl;
	return true;
}


bool topoHdl::load_rvebnd_contacttime(void)
{
	double t1 = MPI_Wtime();
	if (this->bndbuf != NULL) { delete[] this->bndbuf; this->bndbuf = NULL; this->bndbufsize = 0; }

	this->bndbuf = NULL;
	unsigned int ngr = 1 + Settings::LargestGrainID;
	try { this->bndbuf = new MPI_GrainBndContactIO[ngr]; }
	catch (std::bad_alloc &exc) {
		cerr << "ERROR::Worker::LoadRVEBndContacttime" << this->get_Rank() << " file inaccessible!" << endl;
		return false;
	}
	this->bndbufsize = ngr;
	
	ifstream bndfile;
	string bndline;
	istringstream line;
	string datapiece;

	bndfile.open(Settings::GrainBndContactFilename);
	if (bndfile.is_open() == true) {
		//kick header
		getline(bndfile, bndline);

		//read data
		while (bndfile.good() == true) {
			getline(bndfile, bndline);
			if (bndline.size() > 0) {
				istringstream line(bndline);

				getline(line, datapiece, ';');	unsigned int gid = atoi(datapiece.c_str());
				getline(line, datapiece, ';');  unsigned int when = atoi(datapiece.c_str());
				if (gid < ngr) {
					this->bndbuf[gid].extGID = gid;
					this->bndbuf[gid].whenhitboundary = when;
				}
				else {
					cerr << "ERROR::Worker::LoadRVEBndContacttime" << this->get_Rank() << " found invalid ID in GrainBoundaryContactTime file!" << endl;
					bndfile.close();
					return false;
				}
			}
		}
		//done
		bndfile.close();
	}
	else {
		cerr << "ERROR::Worker::LoadRVEBndContacttime" << this->get_Rank() << " unable to open GrainBoundaryContactTime file!" << endl;
		bndfile.close();
		return false;
	}

	double t2 = MPI_Wtime();
	ioprofiler.logev("GrainBndContactLoad", (t2 - t1));
	cout << "...Worker " << this->get_Rank() << " GrainBoundaryContactTime loaded in " << (t2-t1) << " seconds" << endl;
	return true;
}


bool topoHdl::master_mpiio_init_grainsize_quantiles(void)
{
	//matrix of fixed number of quantiles for all timesteps
	int ncol = ((Settings::SnapshotLast - Settings::SnapshotFirst) / Settings::SnapshotOffset + 1); //##MK:: in the future1 + timesteps because add quantile arguments as well
	int nrow = STATISTICS_VOL_NQUANTILES; //the same for all worker
	cout << "...MASTER " << this->get_Rank() << " Writing rows (quantiles) = " << nrow << " columns (time) = " << ncol << endl;
	
	/*double* quantbuffer = NULL;
	try { quantbuffer = new double[nrow]; }
	catch (std::bad_alloc &exc) {
		cerr << "ERROR::MASTER " << this->get_Rank() << " allocating memory in prepping analyze_gsd" << endl; return false;
	}
	//MASTER fills in bin end values first
	for (unsigned int b = 0; b < nrow; ++b) { //CDFs so we seedfbuffer stores binend values
		quantbuffer[b] = ((double)(b+1)) / ((double)STATISTICS_VOL_NQUANTILES);
	}*/

	//open two files quants and meta
	MPI_File ioHdlQUANT, ioHdlMETA;
	MPI_Status ioStaQUANT, ioStaMETA;

	string quant_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".VolQuantiles.F." + std::to_string(Settings::SnapshotFirst) + ".O." + std::to_string(Settings::SnapshotOffset) + ".L." + std::to_string(Settings::SnapshotLast) + ".NC." + std::to_string(ncol) + ".NR." + std::to_string(nrow) + ".bin";
	string meta_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".VolMeta.F." + std::to_string(Settings::SnapshotFirst) + ".O." + std::to_string(Settings::SnapshotOffset) + ".L." + std::to_string(Settings::SnapshotLast) + ".NC." + std::to_string(ncol) + ".NR.1.bin";

	// open the file in create and write-only mode but do not write data
	MPI_File_open(MPI_COMM_SELF, quant_fn.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlQUANT);
	/*long long totalOffsetQUANT = 0 * nrow * 8; //write first columns which are the quantile arguments
	MPI_File_write_at(ioHdlQUANT, totalOffsetQUANT, quantbuffer, nrow, MPI_DOUBLE, &ioStaQUANT); //fractional*/
	MPI_File_close(&ioHdlQUANT);

	MPI_File_open(MPI_COMM_SELF, meta_fn.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlMETA);
	MPI_File_close(&ioHdlMETA);
	
	/*if ( quantbuffer != NULL ) { delete [] quantbuffer; quantbuffer = NULL; }*/
	return true;
}


bool topoHdl::master_mpiio_init_see(void)
{
	int ncol = 1 + ((Settings::SnapshotLast - Settings::SnapshotFirst) / Settings::SnapshotOffset + 1); //1 + timesteps because binends added as well
	int nrow = Settings::get_SEEBinCount();
	cout << "...MASTER " << this->get_Rank() << " Writing rows (SEE bins) = " << nrow << " columns (First col are bin ends, rest are time CDF-data) = " << ncol << endl;
	
	//create I/O buffer container to write timeslot after timeslot, only as many write ops as timeslots because one node has all grains
	double* seedfbuffer = NULL;
	try { seedfbuffer = new double[nrow]; }
	catch (std::bad_alloc &exc) {
		cerr << "ERROR::MASTER " << this->get_Rank() << " allocating memory in prepping analyze_see" << endl; return false;
	}
	
	//MASTER fills in bin end values first
	for (unsigned int b = 0; b < (nrow - 1); ++b) { //CDFs so we seedfbuffer stores binend values
		seedfbuffer[b] = Settings::SeeBinMin + (((double)b) * Settings::SeeBinWidth);
	}
	seedfbuffer[nrow - 1] = SEE_INFTY;

	//open MPI files
	MPI_File ioHdlSEEDF;
	MPI_Status ioStaSEEDF;

	string seedf_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".SEEDF.F." + std::to_string(Settings::SnapshotFirst) + ".O." + std::to_string(Settings::SnapshotOffset) + ".L." + std::to_string(Settings::SnapshotLast) + ".NC." + std::to_string(ncol) + ".NR." + std::to_string(nrow) + ".bin";
	MPI_File_open(MPI_COMM_SELF, seedf_fn.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlSEEDF);
	long long totalOffsetSEEDF = 0 * nrow * 8; //write first columns which are the binends
	MPI_File_write_at(ioHdlSEEDF, totalOffsetSEEDF, seedfbuffer, nrow, MPI_DOUBLE, &ioStaSEEDF); //results are in 1/m^2
	MPI_File_close(&ioHdlSEEDF);

	if (seedfbuffer != NULL) { delete[] seedfbuffer; seedfbuffer = NULL; }
	return true;
}


bool topoHdl::master_mpiio_init_modf(void)
{
	int ncol = 1 + ((Settings::SnapshotLast - Settings::SnapshotFirst) / Settings::SnapshotOffset + 1); //1 + timesteps because binends added as well
	int nrow = Settings::get_MODFBinCount();
	cout << "...MASTER " << this->get_Rank() << " writing rows (MODF bins) = " << nrow << " columns (First col are bin ends, rest are time CDF-data) = " << ncol << endl;

	//create I/O buffer container to write timeslot after timeslot, only as many write ops as timeslots because one node has all grains
	double* modfbuffer = NULL;
	try { modfbuffer = new double[nrow]; }
	catch (std::bad_alloc &exc) {
		cerr << "ERROR::MASTER " << this->get_Rank() << " allocating memory in prepping analyze_modf!" << endl; return false;
	}
	
	for (unsigned int b = 0; b < (nrow - 1); b++) {
		modfbuffer[b] = Settings::DisoriAngleBinMin + RAD2DEG(((double)b) * Settings::DisoriAngleBinWidth); //user input is in degree
	}
	modfbuffer[nrow - 1] = RAD2DEG(MODF_INFTY); //back in degrees for the user

	//open MPI files
	MPI_File ioHdlMODF;
	MPI_Status ioStaMODF;

	string modf_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".MODF.F." + std::to_string(Settings::SnapshotFirst) + ".O." + std::to_string(Settings::SnapshotOffset) + ".L." + std::to_string(Settings::SnapshotLast) + ".NC." + std::to_string(ncol) + ".NR." + std::to_string(nrow) + ".bin";
	MPI_File_open(MPI_COMM_SELF, modf_fn.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlMODF);
	long long totalOffsetMODF = 0 * nrow * 8;
	MPI_File_write_at(ioHdlMODF, totalOffsetMODF, modfbuffer, nrow, MPI_DOUBLE, &ioStaMODF);
	MPI_File_close(&ioHdlMODF);
	
	if (modfbuffer != NULL) { delete[] modfbuffer; modfbuffer = NULL; }
	return true;
}


bool topoHdl::master_mpiio_init_gsd(void)
{
	int ncol = 1 + ((Settings::SnapshotLast - Settings::SnapshotFirst) / Settings::SnapshotOffset + 1); //1 + timesteps because binends added as well
	int nrow = Settings::get_GSDBinCount();
	cout << "...MASTER " << this->get_Rank() << " Writing rows (GSD bins) = " << nrow << " columns (First col are bin ends, rest are time CDF-data) = " << ncol << endl;

	//create I/O buffer container to write timeslot after timeslot, only as many write ops as timeslots because one node has all grains
	MPI_VOLSTATS dummy;
	double* gsdbuffer = NULL;
	try { gsdbuffer = new double[nrow]; }
	catch (std::bad_alloc &exc) {
		cerr << "ERROR::MASTER " << this->get_Rank() << " allocating memory in prepping analyze_gsd binbuffer" << endl; return false;
	}

	dummy.totalsize_ng_elimbnd = 0.0;
	dummy.meansize_ng_elimbnd = 0.0;
	dummy.varsize_ng_elimbnd = 0.0;
	dummy.ng = 0.0;
	dummy.ng_elimbnd = 0.0;

	//MASTER fills in bin end values first
	for (unsigned int b = 0; b < (nrow - 1); ++b) { //CDFs so we seedfbuffer stores binend values
		gsdbuffer[b] = Settings::GSDBinMin + (((double)b) * Settings::GSDBinWidth);
	}
	gsdbuffer[nrow - 1] = GSD_INFTY;

	//open MPI files
	MPI_File ioHdlGSD, ioHdlMETA;
	MPI_Status ioStaGSD, ioStaMETA;

	string meta_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".GSDMeta.F." + std::to_string(Settings::SnapshotFirst) + ".O." + std::to_string(Settings::SnapshotOffset) + ".L." + std::to_string(Settings::SnapshotLast) + ".NC." + std::to_string(ncol) + ".NR.1.bin";
	string gsd_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".GSDHist.F." + std::to_string(Settings::SnapshotFirst) + ".O." + std::to_string(Settings::SnapshotOffset) + ".L." + std::to_string(Settings::SnapshotLast) + ".NC." + std::to_string(ncol) + ".NR." + std::to_string(nrow) + ".bin";

	MPI_File_open(MPI_COMM_SELF, meta_fn.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlMETA);
	long long totalOffsetMETA = 0;
	MPI_File_write_at(ioHdlMETA, totalOffsetMETA, &dummy, 1, MPI_VOLSTATS_Type, &ioStaMETA);
	MPI_File_close(&ioHdlMETA);


	MPI_File_open(MPI_COMM_SELF, gsd_fn.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlGSD);
	long long totalOffsetGSD = 0 * nrow * 8; //write first columns which are the binends
	MPI_File_write_at(ioHdlGSD, totalOffsetGSD, gsdbuffer, nrow, MPI_DOUBLE, &ioStaGSD);
	MPI_File_close(&ioHdlGSD);

	if (gsdbuffer != NULL) { delete[] gsdbuffer; gsdbuffer = NULL; }
	return true;
}


bool topoHdl::master_mpiio_init_drvfrc_see(void)
{
	//MK::one-file per timestep, individual process ASCII I/O therefore no prepping necessary!
	return true;
}


bool topoHdl::master_mpiio_init_maxszgain_fw(void)
{
	//MK::master will accumulate results into single csv file, therefore no prepping necessary, results file instead
	//will be generated in main function
	return true;
}


bool topoHdl::master_mpiio_init_meandrvfrcsee_fw(void)
{
	//MK::master will accumulate results into single csv file, therefore no prepping necessary, results instead in main
	return true;
}


bool topoHdl::master_mpiio_init_rxfrac(void)
{
	//now all know the matrixgrains and the rxgrains so they can compute rxgrains
	int ncol = ((Settings::SnapshotLast - Settings::SnapshotFirst) / Settings::SnapshotOffset + 1);
	int nrow = 1;
	cout << "...MASTER " << this->get_Rank() << " writing rows (only one-dimension of RX_STATS) = " << nrow << " columns (time) = " << ncol << endl;
	
	MPI_File ioHdlRX;
	MPI_Status ioStaRX;

	string rx_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".RXEVO.F." + std::to_string(Settings::SnapshotFirst) + ".O." + std::to_string(Settings::SnapshotOffset) + ".L." + std::to_string(Settings::SnapshotLast) + ".NC." + std::to_string(ncol) + ".NR.1.bin";
	MPI_File_open(MPI_COMM_SELF, rx_fn.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlRX);
	MPI_File_close(&ioHdlRX);
	return true;
}


bool topoHdl::master_mpiio_init_fw(void) 
{
	//targets are only those grains which in their lifetime never touch the boundary
	std::vector<unsigned int>* thetargets = NULL;
	if (Settings::AnalyzeTrajectoriesForwardMode == E_FORWARD_ALL) {
		thetargets = this->analyze_elimbnd_all_self("FW", NULL, false);
	}
	else if (Settings::AnalyzeTrajectoriesForwardMode == E_FORWARD_SELECTED) {
		thetargets = this->analyze_elimbnd_all_self("FW", NULL, true);
	}
	else {
		cerr << "ERROR::MASTER " << this->get_Rank() << " analyze trajectories forward impossible, invalid mode!" << endl;
		return false;
	}

	//catch case when targets vector allocated but filled with nothing
	if (thetargets->size() > 0) {
		int ncol = ((Settings::SnapshotLast - Settings::SnapshotFirst) / Settings::SnapshotOffset + 1);
		int nrow = thetargets->size();
		cout << "...MASTER " << this->get_Rank() << " writing rows (ngrains) = " << nrow << " columns (time) = " << ncol << endl;

		//generate files for processes to fill them later individually
		MPI_File ioHdlFAC, ioHdlVOL, ioHdlHAGB, ioHdlMOBDSEE;
		MPI_Status ioStaFAC, ioStaVOL, ioStaHAGB, ioStaMOBDSEE;

		string suffix = "F." + std::to_string(Settings::SnapshotFirst) + ".O." + std::to_string(Settings::SnapshotOffset) + ".L." + std::to_string(Settings::SnapshotLast) + ".NC." + std::to_string(ncol) + ".NR." + std::to_string(nrow) + ".bin";
		string faces_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".FW.NF." + suffix;
		string vol_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".FW.VOL." + suffix;
		string hagbfrac_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".FW.HAGB." + suffix;
		string mobdsee_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".FW.MOBDSEE." + suffix;

		// open the file in create and write-only mode, ##only MASTER should finally write the file
		MPI_File_open(MPI_COMM_SELF, faces_fn.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlFAC);
		MPI_File_close(&ioHdlFAC);

		MPI_File_open(MPI_COMM_SELF, vol_fn.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlVOL);
		MPI_File_close(&ioHdlVOL);

		MPI_File_open(MPI_COMM_SELF, hagbfrac_fn.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlHAGB);
		MPI_File_close(&ioHdlHAGB);

		MPI_File_open(MPI_COMM_SELF, mobdsee_fn.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlMOBDSEE);
		MPI_File_close(&ioHdlMOBDSEE);
	}
	else {
		cerr << "ERROR::MASTER " << this->get_Rank() << " no targets found!" << endl;
		if (thetargets != NULL) { delete thetargets; thetargets = NULL; }
		return false;
	}

	if (thetargets != NULL) { delete thetargets; thetargets = NULL; }
	return true;
}


bool topoHdl::master_mpiio_init_bk(void)
{
	//necessity to reload Settings::SnapshotLast
	std::vector<unsigned int>* survivors = NULL;
	survivors = analyze_elimbnd_one_mpiioself("TrackingBK", NULL, Settings::SnapshotLast);
	if (survivors == NULL) {
		cerr << "ERROR::MASTER " << this->get_Rank() << " no survivors found!" << endl;
		return false;
	}

	if (survivors->size() > 0) {
		int ncol = ((Settings::SnapshotLast - Settings::SnapshotFirst) / Settings::SnapshotOffset + 1);
		int nrow = survivors->size();
		cout << "...MASTER " << this->get_Rank() << " writing rows (ngrains) = " << nrow << " columns (time) = " << ncol << endl;

		MPI_File ioHdlFAC, ioHdlVOL, ioHdlHAGB, ioHdlMOBDSEE;
		MPI_Status ioStaFAC, ioStaVOL, ioStaHAGB, ioStaMOBDSEE;

		string suffix = "F." + std::to_string(Settings::SnapshotFirst) + ".O." + std::to_string(Settings::SnapshotOffset) + ".L." + std::to_string(Settings::SnapshotLast) + ".NC." + std::to_string(ncol) + ".NR." + std::to_string(nrow) + ".bin";
		string faces_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".BK.NF." + suffix;
		string vol_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".BK.VOL." + suffix;
		string hagbfrac_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".BK.HAGB." + suffix;
		string mobdsee_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".BK.MOBDSEE." + suffix;

		MPI_File_open(MPI_COMM_SELF, faces_fn.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlFAC);
		MPI_File_close(&ioHdlFAC);

		MPI_File_open(MPI_COMM_SELF, vol_fn.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlVOL);
		MPI_File_close(&ioHdlVOL);

		MPI_File_open(MPI_COMM_SELF, hagbfrac_fn.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlHAGB);
		MPI_File_close(&ioHdlHAGB);

		MPI_File_open(MPI_COMM_SELF, mobdsee_fn.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlMOBDSEE);
		MPI_File_close(&ioHdlMOBDSEE);
	}
	else {
		cerr << "ERROR::MASTER " << this->get_Rank() << " no survivors found!" << endl;
		if (survivors != NULL) { delete survivors; survivors = NULL; }
		return false;
	}

	if (survivors != NULL) { delete survivors; survivors = NULL; }
	return true;
}


bool topoHdl::master_mpiio_init_agg(void)
{
	int ncol = ((Settings::SnapshotLast - Settings::SnapshotFirst) / Settings::SnapshotOffset + 1);
	int nrow = 1;
	cout << "...MASTER " << this->get_Rank() << " writing rows (only one-dimension of AGG_STATS) = " << nrow << " columns (time) = " << ncol << endl;

	MPI_File ioHdlAGG;
	MPI_Status ioStaAGG;

	string agg_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".AGGHolmRoll.F." + std::to_string(Settings::SnapshotFirst) + ".O." + std::to_string(Settings::SnapshotOffset) + ".L." + std::to_string(Settings::SnapshotLast) + ".NC." + std::to_string(ncol) + ".NR." + std::to_string(nrow) + ".bin";
	
	MPI_File_open(MPI_COMM_SELF, agg_fn.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlAGG);
	MPI_File_close(&ioHdlAGG);

	return true;
}


bool topoHdl::master_mpiio_init_szgain_bk(void)
{
	std::vector<unsigned int>* survivors = NULL;
	survivors = analyze_elimbnd_one_mpiioself("SzGainVsMatrixSurv", NULL, Settings::SnapshotLast);
	if (survivors == NULL) {
		cerr << "ERROR::MASTER " << this->get_Rank() << " no survivors found!" << endl;
		return false;
	}

	if ( survivors->size() > 0 ) {
		int ncol = ((Settings::SnapshotLast - Settings::SnapshotFirst) / Settings::SnapshotOffset + 1);
		int nrow = 2 + survivors->size(); //+2 for the number of grains in the mean value and the mean value itself
		cout << "...MASTER " << this->get_Rank() << " writing rows (ngrains) = " << nrow << " columns (time) = " << ncol << endl;

		MPI_File ioHdlREL;
		MPI_Status ioStaREL;

		string relsize_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".SIZEGAINBK.F." + std::to_string(Settings::SnapshotFirst) + ".O." + std::to_string(Settings::SnapshotOffset) + ".L." + std::to_string(Settings::SnapshotLast) + ".NC." + std::to_string(ncol) + ".NR." + std::to_string(nrow) + ".bin";
		MPI_File_open(MPI_COMM_SELF, relsize_fn.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlREL);
		MPI_File_close(&ioHdlREL);
	}
	else {
		cerr << "ERROR::MASTER " << this->get_Rank() << " no survivors found!" << endl;
		if (survivors != NULL) { delete survivors; survivors = NULL; }
		return false;
	}

	if (survivors != NULL) { delete survivors; survivors = NULL; }
	return true;
}


bool topoHdl::master_mpiio_init_topodiff_fw(void)
{
	//targets are only those grains which in their lifetime never touch the boundary
	std::vector<unsigned int>* thetargets = NULL;
	if (Settings::AnalyzeTrajectoriesForwardMode == E_FORWARD_ALL) {
		thetargets = this->analyze_elimbnd_all_self("FW", NULL, false);
	}
	else if (Settings::AnalyzeTrajectoriesForwardMode == E_FORWARD_SELECTED) {
		thetargets = this->analyze_elimbnd_all_self("FW", NULL, true);
	}
	else {
		cerr << "ERROR::MASTER " << this->get_Rank() << " analyze trajectories forward impossible, invalid mode!" << endl;
		return false;
	}

	//catch case when targets vector allocated but filled with nothing
	if (thetargets->size() > 0) {
		int ncol = ((Settings::SnapshotLast - Settings::SnapshotFirst) / Settings::SnapshotOffset + 1);
		int nrow = thetargets->size();
		cout << "...MASTER " << this->get_Rank() << " writing rows (ngrains) = " << nrow << " columns (time) = " << ncol << endl;

		//generate files for processes to fill them later individually
		MPI_File ioHdlNFDiff;
		MPI_Status ioStaNFDiff;

		string suffix = "F." + std::to_string(Settings::SnapshotFirst) + ".O." + std::to_string(Settings::SnapshotOffset) + ".L." + std::to_string(Settings::SnapshotLast) + ".NC." + std::to_string(ncol) + ".NR." + std::to_string(nrow) + ".bin";
		string nfdiff_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".FW.NFDIFF." + suffix;

		// open the file in create and write-only mode, ##only MASTER should finally write the file
		MPI_File_open(MPI_COMM_SELF, nfdiff_fn.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlNFDiff);
		MPI_File_close(&ioHdlNFDiff);
	}
	else {
		cerr << "ERROR::MASTER " << this->get_Rank() << " no targets found!" << endl;
		if (thetargets != NULL) { delete thetargets; thetargets = NULL; }
		return false;
	}

	if (thetargets != NULL) { delete thetargets; thetargets = NULL; }
	return true;
}

//bool topoHdl::master_mpiio_init_knn(void) {}
//bool topoHdl::master_mpiio_init_classicalnuc(void) {}


bool topoHdl::init_results_container(void)
{
	//##MK::functions which are commented out are not implemented!
	if (Settings::AnalyzeGrainSizeQuantiles == true) {
		if (master_mpiio_init_grainsize_quantiles() == false)	return false; //MK::exemplarily now for all, if not returned, initialized and alive...
	}

	if (Settings::AnalyzeSEE == true) {
		if (master_mpiio_init_see() == false)					return false;
	}

	if (Settings::AnalyzeMODF == true) {
		if (master_mpiio_init_modf() == false)					return false;
	}

	if (Settings::AnalyzeGSD == true) {
		if (master_mpiio_init_gsd() == false)					return false;
	}

	if (Settings::AnalyzeDrvForceSEE == true) {
		if (master_mpiio_init_drvfrc_see() == false)			return false;
	}

	if (Settings::AnalyzeMaxSizeGainForward == true) {
		if (master_mpiio_init_maxszgain_fw() == false)			return false;
	}

	if (Settings::AnalyzeMeanDrvForceSEEForward == true) {
		if (master_mpiio_init_meandrvfrcsee_fw() == false)		return false;
	}

	if (Settings::AnalyzeApproxRXFraction == true) {
		if (master_mpiio_init_rxfrac() == false)				return false;
	}
	
	if (Settings::AnalyzeTrajectoriesForward == true) {
		if (master_mpiio_init_fw() == false)					return false;
	}
	
	if (Settings::AnalyzeTrajectoriesBackward == true) {
		if (master_mpiio_init_bk() == false)					return false;
	}
	
	if (Settings::AnalyzeAbnormalGrains == true) {
		if (master_mpiio_init_agg() == false)					return false;
	}
	
	if (Settings::AnalyzeSizeGainVsMatrixBackward == true) {
		if (master_mpiio_init_szgain_bk() == false)				return false;
	}

	if (Settings::AnalyzeTopologyDifferenceForward == true) {
		if (master_mpiio_init_topodiff_fw() == false)			return false;
	}

	//if (Settings::AnalyzeKNN == true) {
	//	if (master_mpiio_init_knn() == false) return false;
	//}

	//if (Settings:AnalyzeClassicalNucModels == true) {
	//	if (master_mpiio_init_classicalnuc() == false) return false;
	//}

	return true;
}


bool topoHdl::workpartitioning_fid(unsigned int fid)
{
	//decide if this->get_Rank() processes or not
	unsigned int f = 1 + (fid / Settings::SnapshotOffset) - (Settings::SnapshotFirst / Settings::SnapshotOffset);
	if (f % this->get_nRanks() == this->get_Rank())
		return true;
	else
		return false;
}


void topoHdl::analyze_grainsize_quantiles_adhoc(unsigned int fd)
{
	//get grain size quantile values for snapshot fd
	//worker does so completely independently writing shared MPI I/O to well-defined position in already existing file

	double timer = MPI_Wtime();
	cout << "...Worker " << this->get_Rank() << " calculating cooperatively quantile values ..." << endl;

	//get memory to store for all my time steps the MPI_Quantiles and the MPI_VOLSTATS
	MPI_VOLSTATS myvol;
	double* myquant = NULL;
	try { myquant = new double[STATISTICS_VOL_NQUANTILES]; }
	catch (std::bad_alloc &exc) {
		cerr << "ERROR::Worker::AnalyzeGrainSizeQuantilesADHOC " << this->get_Rank() << " unable to get memory for storing results in myquant" << endl;
		return;
	}
	for (unsigned int i = 0; i < STATISTICS_VOL_NQUANTILES; ++i) 
		myquant[i] = 0.0;

	//process data
	this->LocalDB.at(0)->analyze_vol_quants(&myvol, myquant);
	
	int ncol = ((Settings::SnapshotLast - Settings::SnapshotFirst) / Settings::SnapshotOffset + 1); //##MK::add 1+ to enable in the future the adding of quantile arguments as the first row
	int nrow = STATISTICS_VOL_NQUANTILES;
	
	//open two files quants and meta
	MPI_File ioHdlQUANT, ioHdlMETA;
	MPI_Status ioStaQUANT, ioStaMETA;

	string quant_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".VolQuantiles.F." + std::to_string(Settings::SnapshotFirst) + ".O." + std::to_string(Settings::SnapshotOffset) + ".L." + std::to_string(Settings::SnapshotLast) + ".NC." + std::to_string(ncol) + ".NR." + std::to_string(nrow) + ".bin";
	string meta_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".VolMeta.F." + std::to_string(Settings::SnapshotFirst) + ".O." + std::to_string(Settings::SnapshotOffset) + ".L." + std::to_string(Settings::SnapshotLast) + ".NC." + std::to_string(ncol) + ".NR.1.bin";

	// open the file in create and write-only mode, ##only MASTER should finally write the file
	unsigned int f = (fd / Settings::SnapshotOffset) - (Settings::SnapshotFirst / Settings::SnapshotOffset);
	long long totalOffsetQuant = f * nrow * 8;
	long long totalOffsetMeta = f * 1 * (5 * 8);
	
	MPI_File_open(MPI_COMM_SELF, quant_fn.c_str(), MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlQUANT);
	MPI_File_write_at(ioHdlQUANT, totalOffsetQuant, myquant, nrow, MPI_DOUBLE, &ioStaQUANT);
	MPI_File_close(&ioHdlQUANT);

	MPI_File_open(MPI_COMM_SELF, meta_fn.c_str(), MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlMETA);
	MPI_File_write_at(ioHdlMETA, totalOffsetMeta, &myvol, 1, MPI_VOLSTATS_Type, &ioStaMETA);
	MPI_File_close(&ioHdlMETA);

	if (myquant != NULL) { delete[] myquant; myquant = NULL; }
	cout << "...Worker " << this->get_Rank() << " quantile value computation was successful in " << (double)(MPI_Wtime() - timer) << endl;
}


void topoHdl::analyze_see_adhoc(unsigned int fd)
{
	//get see cdf for snapshot fd
	//worker does so completely independently writing shared MPI I/O to well-defined position in already existing file
	double timer = MPI_Wtime();

	int ncol = 1 + ((Settings::SnapshotLast - Settings::SnapshotFirst) / Settings::SnapshotOffset + 1); //1 + timesteps because binends added as well
	int nrow = Settings::get_SEEBinCount();

	//create I/O buffer container to write timeslot after timeslot, only as many write ops as timeslots because one node has all grains
	double* seedfbuffer = NULL;
	try { seedfbuffer = new double[nrow]; }
	catch (std::bad_alloc &exc) {
		cerr << "ERROR::Worker::AnalyzeSEEADHOC " << this->get_Rank() << " allocating memory in analyze_see" << endl; 
		return;
	}
	
	calc_seedf( 0, seedfbuffer);

	//open MPI files
	string seedf_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".SEEDF.F." + std::to_string(Settings::SnapshotFirst) + ".O." + std::to_string(Settings::SnapshotOffset) + ".L." + std::to_string(Settings::SnapshotLast) + ".NC." + std::to_string(ncol) + ".NR." + std::to_string(nrow) + ".bin";
	MPI_File ioHdlSEEDF;
	MPI_Status ioStaSEEDF;

	unsigned int f = (fd / Settings::SnapshotOffset) - (Settings::SnapshotFirst / Settings::SnapshotOffset);
	long long totalOffsetSEEDF = (1 * nrow * 8) + (f * nrow * 8);
	
	MPI_File_open(MPI_COMM_SELF, seedf_fn.c_str(), MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlSEEDF);
	MPI_File_write_at(ioHdlSEEDF, totalOffsetSEEDF, seedfbuffer, nrow, MPI_DOUBLE, &ioStaSEEDF); //results are in radiant
	MPI_File_close(&ioHdlSEEDF);
	
	if (seedfbuffer != NULL) { delete[] seedfbuffer; seedfbuffer = NULL; }

	cout << "...Worker " << this->get_Rank() << " analyzing SEEDF = " << (MPI_Wtime() - timer) << " seconds" << endl;
}


void topoHdl::analyze_modf_adhoc(unsigned int fd)
{
	//get see cdf for snapshot fd
	//worker does so completely independently writing shared MPI I/O to well-defined position in already existing file
	double timer = MPI_Wtime();

	int ncol = 1 + ((Settings::SnapshotLast - Settings::SnapshotFirst) / Settings::SnapshotOffset + 1); //1 + timesteps because binends added as well
	int nrow = Settings::get_MODFBinCount();

	//create I/O buffer container to write timeslot after timeslot, only as many write ops as timeslots because one node has all grains
	double* modfbuffer = NULL;
	try { modfbuffer = new double[nrow]; }
	catch (std::bad_alloc &exc) {
		cerr << "ERROR::Worker::AnalyzeMODFADHOC " << this->get_Rank() << " allocating memory in analyze_modf!" << endl; return;
	}
	
	calc_modf( 0, modfbuffer);

	string modf_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".MODF.F." + std::to_string(Settings::SnapshotFirst) + ".O." + std::to_string(Settings::SnapshotOffset) + ".L." + std::to_string(Settings::SnapshotLast) + ".NC." + std::to_string(ncol) + ".NR." + std::to_string(nrow) + ".bin";

	//open MPI files
	MPI_File ioHdlMODF;
	MPI_Status ioStaMODF;
	
	unsigned int f = (fd / Settings::SnapshotOffset) - (Settings::SnapshotFirst / Settings::SnapshotOffset);
	long long totalOffsetMODF = (1 * nrow * 8) + (f * nrow * 8);

	// open the file in create and write-only mode, ##only MASTER should finally write the file
	MPI_File_open(MPI_COMM_SELF, modf_fn.c_str(), MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlMODF);
	MPI_File_write_at(ioHdlMODF, totalOffsetMODF, modfbuffer, nrow, MPI_DOUBLE, &ioStaMODF); //results are in radiant
	MPI_File_close(&ioHdlMODF);	

	if (modfbuffer != NULL) { delete[] modfbuffer; modfbuffer = NULL; }
	cout << "Worker " << this->get_Rank() << " analyzing MODF = " << (MPI_Wtime() - timer) << " seconds" << endl;
}


void topoHdl::analyze_gsd_adhoc(unsigned int fd)
{
	//get arithmetic average grain size normalized grain sizes histogram for snapshot fd
	//worker does so completely independently writing shared MPI I/O to well-defined position in already existing file

	double timer = MPI_Wtime();
	cout << "...Worker " << this->get_Rank() << " calculating cooperatively grain size distribution histogram ..." << endl;

	int ncol = 1 + ((Settings::SnapshotLast - Settings::SnapshotFirst) / Settings::SnapshotOffset + 1); //1 + timesteps because binends added as well
	int nrow = Settings::get_GSDBinCount();

	//get memory to buffer histogram cnts
	MPI_VOLSTATS myvol;
	double* histcnt = NULL;
	try { histcnt = new double[nrow]; }
	catch (std::bad_alloc &exc) {
		cerr << "ERROR::Worker::AnalyzeGSDADHOC " << this->get_Rank() << " unable to get memory for storing results in histcnt" << endl;
		return;
	}
	for (unsigned int i = 0; i < nrow; ++i) 
		histcnt[i] = 0.0;

	//process data
	this->LocalDB.at(0)->analyze_gsd_histcnt(&myvol, histcnt);

	//open two files quants and meta
	MPI_File ioHdlGSD, ioHdlMETA;
	MPI_Status ioStaGSD, ioStaMETA;

	string meta_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".GSDMeta.F." + std::to_string(Settings::SnapshotFirst) + ".O." + std::to_string(Settings::SnapshotOffset) + ".L." + std::to_string(Settings::SnapshotLast) + ".NC." + std::to_string(ncol) + ".NR.1.bin";
	string gsd_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".GSDHist.F." + std::to_string(Settings::SnapshotFirst) + ".O." + std::to_string(Settings::SnapshotOffset) + ".L." + std::to_string(Settings::SnapshotLast) + ".NC." + std::to_string(ncol) + ".NR." + std::to_string(nrow) + ".bin";

	// open the file in create and write-only mode, ##only MASTER should finally write the file
	unsigned int f = (fd / Settings::SnapshotOffset) - (Settings::SnapshotFirst / Settings::SnapshotOffset);
	long long totalOffsetMeta = (5 * 8) + f * 1 * (5 * 8);
	long long totalOffsetGSD = (1 * nrow * 8) + f * nrow * 8;

	MPI_File_open(MPI_COMM_SELF, meta_fn.c_str(), MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlMETA);
	MPI_File_write_at(ioHdlMETA, totalOffsetMeta, &myvol, 1, MPI_VOLSTATS_Type, &ioStaMETA);
	MPI_File_close(&ioHdlMETA);

	MPI_File_open(MPI_COMM_SELF, gsd_fn.c_str(), MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlGSD);
	MPI_File_write_at(ioHdlGSD, totalOffsetGSD, &(histcnt[0]), nrow, MPI_DOUBLE, &ioStaGSD);
	MPI_File_close(&ioHdlGSD);

	if (histcnt != NULL) { delete [] histcnt; histcnt = NULL; }
	cout << "...Worker " << this->get_Rank() << " GSD histogram computation was successful in " << (double)(MPI_Wtime() - timer) << endl;
}


void topoHdl::analyze_drvforce_see_adhoc(unsigned int fd)
{
	//for all grains in a timestep which never touch the boundary compute driving force difference over perimeter
	double timer = MPI_Wtime();

	unsigned int f = (fd / Settings::SnapshotOffset) - (Settings::SnapshotFirst / Settings::SnapshotOffset);
	unsigned int localid = 0; //MK::because each process has only one database object at a time!

	string log_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".DrvForceSEE.FID." + std::to_string(fd) + ".csv";
	ofstream logsee;
	logsee.open(log_fn.c_str(), std::ofstream::out | std::ofstream::trunc);
	if (logsee.is_open() == true) {
		logsee << "GrainID;SegmentLengthAveragedSEEDrvForce\n";
		logsee << ";Pa\n";

		//process each grain
		for (unsigned int mr = 0; mr < this->LocalDB[localid]->mrg.size(); mr++) {
			unsigned int mrngr = this->LocalDB[localid]->mrg[mr]->GrainBucket.size();
			for (unsigned int g = 0; g < mrngr; ++g) {
				if (this->LocalDB[localid]->mrg[mr]->GrainBucket[g].boundary != BOUNDARY_CONTACT) {
					unsigned int gid = this->LocalDB[localid]->mrg[mr]->GrainBucket[g].extGID;
					struct target_props targetp;
					targetp = this->getTargetProperties(gid, localid, false, false, false, false, true, false);
					logsee << gid << ";" << setprecision(18) << targetp.dsee << "\n";
				}
			}
		}
		logsee.flush();
		logsee.close();
	}
	else {
		cerr << "ERROR::Worker::AnalyzeDrvForceSEEADHOC " << this->get_Rank() << " unable to open file " << log_fn.c_str() << endl;
	}
	
	cout << "...Worker " << this->get_Rank() << " analyzing driving force forward = " << (MPI_Wtime() - timer) << " seconds" << endl;
}


void topoHdl::init_helper(void)
{
	unsigned int ngr = 1 + Settings::LargestGrainID;

	if ( Settings::AnalyzeMaxSizeGainForward == true ) {
		this->helper_maxszgain_fw = NULL;
		try { this->helper_maxszgain_fw = new double[ngr]; }
		catch (std::bad_alloc &exc) {
			cerr << "ERROR::Worker::TrackingSeqInitHelper " << this->get_Rank() << " unable to allocate helper_maxszgain_fw!" << endl;
			return;
		}

		for (unsigned int g = 0; g < ngr; ++g) {
			helper_maxszgain_fw[g] = 0.0;
		}
	}

	if ( Settings::AnalyzeMeanDrvForceSEEForward == true ) {
		this->helper_meandrvfsee_fw = NULL;
		try { this->helper_meandrvfsee_fw = new seeav[ngr]; }
		catch (std::bad_alloc &exc) {
			cerr << "ERROR::Worker::TrackingSeqInitHelper " << this->get_Rank() << " unable to allocate helper_meandrvfsee_fw!" << endl;
			return;
		}

		for (unsigned int g = 0; g < ngr; ++g) {
			helper_meandrvfsee_fw[g].gid = g;
		}
	}
}


void topoHdl::consolidate_helper(void)
{
	//so far the processes have just measured the maximum size a grain ever obtained in their own local data
	//but we still have to synchronized globally, via collective MPI
	double timer = MPI_Wtime();

	//MK::settings the same within all processes in MPI_COMM_WORLD
	if ( Settings::AnalyzeMaxSizeGainForward == true ) {

		//MPI_Barrier(MPI_COMM_WORLD);
		unsigned int ngr = 1 + Settings::LargestGrainID;

		if (this->myRank == MASTER) {
			MPI_Reduce( MPI_IN_PLACE, this->helper_maxszgain_fw, ngr, MPI_DOUBLE, MPI_MAX, MASTER, MPI_COMM_WORLD);
		}
		else {
			MPI_Reduce( this->helper_maxszgain_fw, this->helper_maxszgain_fw, ngr, MPI_DOUBLE, MPI_MAX, MASTER, MPI_COMM_WORLD);
		}

		//##MK::necessary, to assure that all data to output were reduced in the MASTER
		MPI_Barrier(MPI_COMM_WORLD);

		if (this->get_Rank() == MASTER) {
			string log_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".MaxSizeGainFW.F." + std::to_string(Settings::SnapshotFirst) + ".O." + std::to_string(Settings::SnapshotOffset) + ".L." + std::to_string(Settings::SnapshotLast) + ".csv";
			ofstream logmxsz;
			logmxsz.open(log_fn.c_str(), std::ofstream::out | std::ofstream::trunc);
			if (logmxsz.is_open() == true) {
				//header
				logmxsz << "GrainID;MaxSize\n";
				if (Settings::Dimensionality == TWO_DIMENSIONS)
					logmxsz << ";micron^2\n";
				if (Settings::Dimensionality == THREE_DIMENSIONS)
					logmxsz << ";micron^3\n";

				//normalized area on [0,1]^d to physical size in micron^d
				double norm2real = SQR(METER2MICRON(Settings::PhysicalDomainSize)); //TWO_DIMENSIONS
				if (Settings::Dimensionality == THREE_DIMENSIONS)
					norm2real = CUBE(METER2MICRON(Settings::PhysicalDomainSize));

				//data writing
				for (unsigned int gid = 0; gid < ngr; ++gid) {
					if ( this->helper_maxszgain_fw[gid] > DOUBLE_EPSILON) //only grains != BOUNDARY_CONTACT are accounted for their maximum size
						logmxsz << gid << ";" << setprecision(18) << this->helper_maxszgain_fw[gid] * norm2real << "\n";
				}
				logmxsz.flush();
				logmxsz.close();
	
				cout << "...Worker " << this->get_Rank() << " consolidated process local intermediate data in " << (double)(MPI_Wtime() - timer) << endl;
			}
			else { cerr << "ERROR::Worker::ConsolidateMaxSizeGainFW " << this->get_Rank() << " is unable to open results file " << log_fn << endl; }
		}

	}

	if ( Settings::AnalyzeMeanDrvForceSEEForward == true && this->nRanks == 1 ) { //MK::hardcoded safety as algorithm currently only for single process

		unsigned int ngr = 1 + Settings::LargestGrainID;

		//MK::in this special case master has already all pieces of information so write results
		if (this->get_Rank() == MASTER) {
			string log_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".AvDrvFrcSEEFW.F." + std::to_string(Settings::SnapshotFirst) + ".O." + std::to_string(Settings::SnapshotOffset) + ".L." + std::to_string(Settings::SnapshotLast) + ".csv";
			ofstream logavsee;
			logavsee.open(log_fn.c_str(), std::ofstream::out | std::ofstream::trunc);
			if (logavsee.is_open() == true) {
				//header
				logavsee << "GrainID;AvMigrationSpeed;AvDrivingForce;N\n";
				logavsee << ";m/s;Pa;1\n";

				//normalized area on [0,1]^d to physical size in micron^d
				double norm2real = SQR(METER2MICRON(Settings::PhysicalDomainSize)); //TWO_DIMENSIONS
				if (Settings::Dimensionality == THREE_DIMENSIONS)
					norm2real = CUBE(METER2MICRON(Settings::PhysicalDomainSize));

				//data writing
				for (unsigned int gid = 0; gid < ngr; ++gid) {
					if ( this->helper_meandrvfsee_fw[gid].n > 0) { //only grains that existed
						double _n = 1.0 / ((double) this->helper_meandrvfsee_fw[gid].n);
						logavsee << gid << ";" << setprecision(18) << this->helper_meandrvfsee_fw[gid].psee * _n << ";" << this->helper_meandrvfsee_fw[gid].vsee << ";" << this->helper_meandrvfsee_fw[gid].n << "\n";
					}
				}
				logavsee.flush();
				logavsee.close();
	
				cout << "...Worker " << this->get_Rank() << " consolidated process local intermediate data in " << (double)(MPI_Wtime() - timer) << endl;
			}
			else { cerr << "ERROR::Worker::ConsolidateMeanDrvForceSEEFW " << this->get_Rank() << " is unable to open results file " << log_fn << endl; }
		}

	}

}


void topoHdl::analyze_maxsizegain_fw_adhoc(unsigned int fd)
{
	unsigned int localid = 0;

	calc_localmaxsize( localid, this->helper_maxszgain_fw );
}


void topoHdl::analyze_meandrvforcesee_adhoc(unsigned int fd)
{
	unsigned int localid = 0;

	calc_localavdrvfrc( localid, this->helper_meandrvfsee_fw );
}


void topoHdl::analyze_rxfraction_adhoc(unsigned int fd, std::vector<unsigned int>* matrix, std::vector<unsigned int>* rx)
{
	//identifies all remaining grains that never touched the boundary
	//##MK::here we assume that all grains that in time step extid do not touch the boundary also before never touched the boundary
	//computes thereafter the area for these grains in relation to the total area for each timestep of all grains which do not touch the boundary
	double timer = MPI_Wtime();
	if (matrix == NULL || rx == NULL) {
		cerr << "ERROR::Worker::AnalyzeApproxRXADHOC " << this->get_Rank() << " filter input is non-existent!" << endl;
		return;
	}

	//create I/O buffer container to write timeslot after timeslot, only as many write ops as timeslots because one node has all grains
	MPI_RXSTATS* rxbuffer = NULL;
	try { rxbuffer = new MPI_RXSTATS[1]; }
	catch (std::bad_alloc &exc) {
		cerr << "ERROR::Worker::AnalyzeApproxRXADHOC " << this->get_Rank() << " results for MPI_RXSTATS not allocatable!" << endl;
		return;
	}

	unsigned int localid = 0;
	struct rxstats rxs;
	rxs = calcRXFraction(localid, matrix, rx);

	rxbuffer[0].totalsize_elimbnd = rxs.totalsize_elimbnd;
	rxbuffer[0].totalsize_targets = rxs.totalsize_targets;
	rxbuffer[0].X = rxs.X;
	rxbuffer[0].ndefg = rxs.ndefg;
	rxbuffer[0].nrxg = rxs.nrxg;

	//write in correct position in already generated file, independently
	MPI_File ioHdlRX;
	MPI_Status ioStaRX;

	unsigned int f = (fd / Settings::SnapshotOffset) - (Settings::SnapshotFirst / Settings::SnapshotOffset);
	int ncol = ((Settings::SnapshotLast - Settings::SnapshotFirst) / Settings::SnapshotOffset + 1);
	int nrow = 1;
	long long totalOffsetRX = f * nrow * (5 * 8); //because nrow=1;

	string rx_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".RXEVO.F." + std::to_string(Settings::SnapshotFirst) + ".O." + std::to_string(Settings::SnapshotOffset) + ".L." + std::to_string(Settings::SnapshotLast) + ".NC." + std::to_string(ncol) + ".NR.1.bin";
	
	MPI_File_open(MPI_COMM_SELF, rx_fn.c_str(), MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlRX);
	MPI_File_write_at(ioHdlRX, totalOffsetRX, rxbuffer, 1, MPI_RXSTATS_Type, &ioStaRX);
	MPI_File_close(&ioHdlRX);

	if (rxbuffer != NULL) {
		delete[] rxbuffer; rxbuffer = NULL;
	}

	cout << "...Worker::AnalyzeRXFractionADHOC " << this->get_Rank() << " analyzed approximate recrystallized volume fraction = " << (MPI_Wtime() - timer) << " seconds" << endl;
}


void topoHdl::analyze_trajectories_fw_adhoc(unsigned int fd, std::vector<unsigned int>* thetargets)
{
	double timer = MPI_Wtime();

	if (thetargets == NULL) {
		cerr << "ERROR::Worker::AnalyzeTrackingFWADHOC " << this->get_Rank() << " filter input is non-existent!" << endl;
		return;
	}

	//all workers have the same identically arranged list of target grain IDs because they performed the identical analysis routine beforehand
	int ncol = ((Settings::SnapshotLast - Settings::SnapshotFirst) / Settings::SnapshotOffset + 1);
	int nrow = thetargets->size();
	cout << "...Worker::AnalyzeForwardTracking " << this->get_Rank() << " writing rows (ngrains) = " << nrow << " columns (time) = " << ncol << endl;

	//create I/O buffer container to write timeslot after timeslot, only as many write ops as timeslots because one node has all grains
	unsigned int* nfacesbuffer = NULL;
	try { nfacesbuffer = new unsigned int[nrow]; }
	catch (std::bad_alloc &exc) {
		cerr << "ERROR::Worker::AnalyzeTrackingFWADHOC " << this->get_Rank() << " allocating memory for nfacesbuffer failed!" << endl;
	}
	double* volbuffer = NULL;
	try { volbuffer = new double[nrow]; }
	catch (std::bad_alloc &exc) {
		cerr << "ERROR::Worker::AnalyzeTrackingFWADHOC " << this->get_Rank() << " allocating memory for volbuffer failed!" << endl;
		if (nfacesbuffer != NULL) { delete[] nfacesbuffer; nfacesbuffer = NULL; }
	}
	double* hagbfracbuffer = NULL;
	try { hagbfracbuffer = new double[nrow]; }
	catch (std::bad_alloc &exc) {
		cerr << "ERROR::Worker::AnalyzeTrackingFWADHOC " << this->get_Rank() << " allocating memory for hagbfracbuffer failed!" << endl;
		if (nfacesbuffer != NULL) { delete[] nfacesbuffer; nfacesbuffer = NULL; }
		if (volbuffer != NULL) { delete[] volbuffer; volbuffer = NULL; }
	}
	double* mobdseebuffer = NULL;
	try { mobdseebuffer = new double[nrow]; }
	catch (std::bad_alloc &exc) {
		cerr << "ERROR::Worker::AnalyzeTrackingFWADHOC " << this->get_Rank() << " allocating memory for mobdseebuffer failed!" << endl;
		if (nfacesbuffer != NULL) { delete[]  nfacesbuffer; nfacesbuffer = NULL; }
		if (volbuffer != NULL) { delete[]  volbuffer; volbuffer = NULL; }
		if (hagbfracbuffer != NULL) { delete[]  hagbfracbuffer; hagbfracbuffer = NULL; }
	}

	unsigned int f = (fd / Settings::SnapshotOffset) - (Settings::SnapshotFirst / Settings::SnapshotOffset);
	unsigned int localid = 0;

	for (unsigned int target = 0; target < thetargets->size(); target++) {
		unsigned int gid = thetargets->at(target); //assured to be a grain not at the boundary and not THE_DOMAIN_ITSELF

		struct target_props targetp;
		targetp = this->getTargetProperties(gid, localid, true, true, true, true, false, false);

		nfacesbuffer[target] = targetp.nfacesk1;
		volbuffer[target] = targetp.area;
		hagbfracbuffer[target] = targetp.hagbfraction;
		mobdseebuffer[target] = targetp.mob_dsee;
	}

	//calculate address in already existent file where to store the results
	long long totalOffsetFaces = f * nrow * 4;
	long long totalOffsetVolume = f * nrow * 8;
	long long totalOffsetHAGB = f * nrow * 8;
	long long totalOffsetMOBDSEE = f * nrow * 8;

	MPI_File ioHdlFAC, ioHdlVOL, ioHdlHAGB, ioHdlMOBDSEE;
	MPI_Status ioStaFAC, ioStaVOL, ioStaHAGB, ioStaMOBDSEE;

	string suffix = "F." + std::to_string(Settings::SnapshotFirst) + ".O." + std::to_string(Settings::SnapshotOffset) + ".L." + std::to_string(Settings::SnapshotLast) + ".NC." + std::to_string(ncol) + ".NR." + std::to_string(nrow) + ".bin";
	string faces_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".FW.NF." + suffix;
	string vol_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".FW.VOL." + suffix;
	string hagbfrac_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".FW.HAGB." + suffix;
	string mobdsee_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".FW.MOBDSEE." + suffix;

	//write to well-defined position in already existent file
	MPI_File_open(MPI_COMM_SELF, faces_fn.c_str(), MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlFAC);
	MPI_File_open(MPI_COMM_SELF, vol_fn.c_str(), MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlVOL);
	MPI_File_open(MPI_COMM_SELF, hagbfrac_fn.c_str(), MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlHAGB);
	MPI_File_open(MPI_COMM_SELF, mobdsee_fn.c_str(), MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlMOBDSEE);

	MPI_File_write_at(ioHdlFAC, totalOffsetFaces, nfacesbuffer, nrow, MPI_UNSIGNED, &ioStaFAC);
	MPI_File_write_at(ioHdlVOL, totalOffsetVolume, volbuffer, nrow, MPI_DOUBLE, &ioStaVOL);
	MPI_File_write_at(ioHdlHAGB, totalOffsetHAGB, hagbfracbuffer, nrow, MPI_DOUBLE, &ioStaHAGB);
	MPI_File_write_at(ioHdlMOBDSEE, totalOffsetMOBDSEE, mobdseebuffer, nrow, MPI_DOUBLE, &ioStaMOBDSEE);
	
	MPI_File_close(&ioHdlFAC);
	MPI_File_close(&ioHdlVOL);
	MPI_File_close(&ioHdlHAGB);
	MPI_File_close(&ioHdlMOBDSEE);

	if (nfacesbuffer != NULL) { delete[]  nfacesbuffer; nfacesbuffer = NULL; }
	if (volbuffer != NULL) { delete[]  volbuffer; volbuffer = NULL; }
	if (hagbfracbuffer != NULL) { delete[]  hagbfracbuffer; hagbfracbuffer = NULL; }
	if (mobdseebuffer != NULL) { delete[] mobdseebuffer; mobdseebuffer = NULL; }

	cout << "...Worker::AnalyzeTrackingFWADHOC " << this->get_Rank() << " analyzing trajectories forward = " << (MPI_Wtime() - timer) << " seconds" << endl;
}


void topoHdl::analyze_trajectories_bk_adhoc(unsigned int fd , std::vector<unsigned int>* survivors)
{
	double timer = MPI_Wtime();

	if (survivors == NULL ) {
		cerr << "ERROR::Worker::AnalyzeTrackingBWADHOC " << this->get_Rank() << " filter input is non-existent!" << endl;
		return;
	}
	
	unsigned int nGrainsSurvived = survivors->size();

	//all workers have the same identically arranged list of survivor grain IDs because they performed the identical survivor analysis routine beforehand
	int ncol = ((Settings::SnapshotLast - Settings::SnapshotFirst) / Settings::SnapshotOffset + 1);
	int nrow = nGrainsSurvived;
	cout << "...Worker::AnalyzeBackwardsTracking " << this->get_Rank() << " writing rows (ngrains) = " << nrow << " columns (time) = " << ncol << endl;

	//create I/O buffer container to write timeslot after timeslot, only as many write ops as timeslots because one node has all grains
	unsigned int* nfacesbuffer = NULL;
	try { nfacesbuffer = new unsigned int[nrow]; }
	catch (std::bad_alloc &exc) {
		cerr << "ERROR::Worker::AnalyzeTrackingBWADHOC " << this->get_Rank() << " allocating memory for nfacesbuffer failed!" << endl;
	}
	double* volbuffer = NULL;
	try { volbuffer = new double[nrow]; }
	catch (std::bad_alloc &exc) {
		cerr << "ERROR::Worker::AnalyzeTrackingBWADHOC " << this->get_Rank() << " allocating memory for volbuffer failed!" << endl;
		if (nfacesbuffer != NULL) { delete[] nfacesbuffer; nfacesbuffer = NULL; }
	}
	double* hagbfracbuffer = NULL;
	try { hagbfracbuffer = new double[nrow]; }
	catch (std::bad_alloc &exc) {
		cerr << "ERROR::Worker::AnalyzeTrackingBWADHOC " << this->get_Rank() << " allocating memory for hagbfracbuffer failed!" << endl;
		if (nfacesbuffer != NULL) { delete[] nfacesbuffer; nfacesbuffer = NULL; }
		if (volbuffer != NULL) { delete[] volbuffer; volbuffer = NULL; }
	}
	double* mobdseebuffer = NULL;
	try { mobdseebuffer = new double[nrow]; }
	catch (std::bad_alloc &exc) {
		cerr << "ERROR::Worker::AnalyzeTrackingBWADHOC " << this->get_Rank() << " allocating memory for mobdseebuffer failed!" << endl;
		if (nfacesbuffer != NULL) { delete[]  nfacesbuffer; nfacesbuffer = NULL; }
		if (volbuffer != NULL) { delete[]  volbuffer; volbuffer = NULL; }
		if (hagbfracbuffer != NULL) { delete[]  hagbfracbuffer; hagbfracbuffer = NULL; }
	}

	unsigned int f = (fd / Settings::SnapshotOffset) - (Settings::SnapshotFirst / Settings::SnapshotOffset);
	unsigned int localid = 0;
	
	for (unsigned int target = 0; target < nGrainsSurvived; target++) {
		unsigned int gid = survivors->at(target); //assured to be a grain not at the boundary and not THE_DOMAIN_ITSELF

		struct target_props targetp;
		targetp = this->getTargetProperties(gid, localid, true, true, true, true, false, false);

		nfacesbuffer[target] = targetp.nfacesk1;
		volbuffer[target] = targetp.area;
		hagbfracbuffer[target] = targetp.hagbfraction;
		mobdseebuffer[target] = targetp.mob_dsee;
	}
	
	//calculate address in already existent file where to store the results
	long long totalOffsetFaces = f * nrow * 4;
	long long totalOffsetVolume = f * nrow * 8;
	long long totalOffsetHAGB = f * nrow * 8;
	long long totalOffsetMOBDSEE = f * nrow * 8;

	MPI_File ioHdlFAC, ioHdlVOL, ioHdlHAGB, ioHdlMOBDSEE;
	MPI_Status ioStaFAC, ioStaVOL, ioStaHAGB, ioStaMOBDSEE;

	string suffix = "F." + std::to_string(Settings::SnapshotFirst) + ".O." + std::to_string(Settings::SnapshotOffset) + ".L." + std::to_string(Settings::SnapshotLast) + ".NC." + std::to_string(ncol) + ".NR." + std::to_string(nrow) + ".bin";
	string faces_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".BK.NF." + suffix;
	string vol_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".BK.VOL." + suffix;
	string hagbfrac_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".BK.HAGB." + suffix;
	string mobdsee_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".BK.MOBDSEE." + suffix;
	
	//write to well-defined position in already existent file
	MPI_File_open(MPI_COMM_SELF, faces_fn.c_str(), MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlFAC);
	MPI_File_open(MPI_COMM_SELF, vol_fn.c_str(), MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlVOL);
	MPI_File_open(MPI_COMM_SELF, hagbfrac_fn.c_str(), MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlHAGB);
	MPI_File_open(MPI_COMM_SELF, mobdsee_fn.c_str(), MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlMOBDSEE);

	MPI_File_write_at(ioHdlFAC, totalOffsetFaces, nfacesbuffer, nrow, MPI_UNSIGNED, &ioStaFAC);
	MPI_File_write_at(ioHdlVOL, totalOffsetVolume, volbuffer, nrow, MPI_DOUBLE, &ioStaVOL);
	MPI_File_write_at(ioHdlHAGB, totalOffsetHAGB, hagbfracbuffer, nrow, MPI_DOUBLE, &ioStaHAGB);
	MPI_File_write_at(ioHdlMOBDSEE, totalOffsetMOBDSEE, mobdseebuffer, nrow, MPI_DOUBLE, &ioStaMOBDSEE);
	
	MPI_File_close(&ioHdlFAC);
	MPI_File_close(&ioHdlVOL);
	MPI_File_close(&ioHdlHAGB);
	MPI_File_close(&ioHdlMOBDSEE);

	if (nfacesbuffer != NULL) { delete[]  nfacesbuffer; nfacesbuffer = NULL; }
	if (volbuffer != NULL) { delete[]  volbuffer; volbuffer = NULL; }
	if (hagbfracbuffer != NULL) { delete[]  hagbfracbuffer; hagbfracbuffer = NULL; }
	if (mobdseebuffer != NULL) { delete[] mobdseebuffer; mobdseebuffer = NULL; }
	
	cout << "...Worker::AnalyzeTrackingBKADHOC " << this->get_Rank() << " analyzing trajectories backwards = " << (MPI_Wtime() - timer) << " seconds" << endl;
}


void topoHdl::analyze_abnormalgraingrowth_adhoc(unsigned int fd, std::vector<unsigned int>* survivors, std::vector<unsigned int>* matrix, std::vector<unsigned int>* matrix0surv )
{
	//Holm, Miodownik and Rollett in Acta Mat 2003, 51, 2701 analyzed abnormal sub-grain coarsening
	double timer = MPI_Wtime();

	if (survivors == NULL || matrix == NULL || matrix0surv == NULL) {
		cerr << "ERROR::Worker::AnalyzeAbnormalGrainGrowthADHOC " << this->get_Rank() << " filter input is non-existent or incomplete!" << endl;
		return;
	}

	//this function computes the instantaneous descriptive means of the matrix grains and counts how many grains are considered as to grow at the moment abnormally
	//we define three lists of grain ids
	//survivors - all the grains without boundary contact in the final data set, most likely AGG grains and remaining matrix
	//matrix - all grains ever considered with never boundary contact (can include also survivors)
	//matrix0surv - excludes the survivors

	//now all know the matrixgrains and the rxgrains so they can compute rxgrains
	int ncol = ((Settings::SnapshotLast - Settings::SnapshotFirst) / Settings::SnapshotOffset + 1);
	int nrow = 1;
	
	//create I/O buffer container to write timeslot after timeslot, only as many write ops as timeslots because one node has all grains
	MPI_AGGSTATS* aggbuffer = NULL;
	try { aggbuffer = new MPI_AGGSTATS[1]; }
	catch (std::bad_alloc &exc) {
		cerr << "ERROR::Worker::AnalyzeAbnormalGrainGrowthADHOC " << this->get_Rank() << " unable to allocate result memory" << endl;
		return;
	}

	struct aggstats agr;
	unsigned int localid = 0;

	agr = inspectForAbnormalGrains(localid, survivors, matrix, matrix0surv);

	//cout << agr.totalsize_survivors << ";" << agr.totalsize_matrix << ";" << agr.totalsize_targets << ";" << agr.rmean_survivors << ";" << agr.rmean_matrix << ";" << agr.rmean_targets << agr.n_survivors << ";" << agr.n_matrix << ";" << agr.n_targets << endl;
	aggbuffer[0].totalsize_survivors = agr.totalsize_survivors;
	aggbuffer[0].totalsize_matrix = agr.totalsize_matrix;
	aggbuffer[0].totalsize_targets = agr.totalsize_targets;

	aggbuffer[0].rmean_survivors = agr.rmean_survivors;
	aggbuffer[0].rmean_matrix = agr.rmean_matrix;
	aggbuffer[0].rmean_targets = agr.rmean_targets;

	aggbuffer[0].n_survivors = agr.n_survivors;
	aggbuffer[0].n_matrix = agr.n_matrix;
	aggbuffer[0].n_targets = agr.n_targets;
	aggbuffer[0].n_agg_m = agr.n_agg_m;
	aggbuffer[0].n_agg_t = agr.n_agg_t;
	aggbuffer[0].n_agg_m_insurv = agr.n_agg_m_insurv;
	aggbuffer[0].n_agg_t_insurv = agr.n_agg_t_insurv;

	string agg_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".AGGHolmRoll.F." + std::to_string(Settings::SnapshotFirst) + ".O." + std::to_string(Settings::SnapshotOffset) + ".L." + std::to_string(Settings::SnapshotLast) + ".NC." + std::to_string(ncol) + ".NR." + std::to_string(nrow) + ".bin";
	unsigned int f = (fd / Settings::SnapshotOffset) - (Settings::SnapshotFirst / Settings::SnapshotOffset);
	long long totalOffsetAGG = f * nrow * ((3 + 3 + 3 + 4) * 8);
	
	MPI_File ioHdlAGG;
	MPI_Status ioStaAGG;
	
	MPI_File_open(MPI_COMM_SELF, agg_fn.c_str(), MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlAGG);
	MPI_File_write_at(ioHdlAGG, totalOffsetAGG, aggbuffer, 1, MPI_AGGSTATS_Type, &ioStaAGG);
	MPI_File_close(&ioHdlAGG);

	if (aggbuffer != NULL) { 
		delete[] aggbuffer; aggbuffer = NULL; 
	}

	cout << "...Worker " << this->get_Rank() << " analyzed AGG Holm Rollett = " << (MPI_Wtime() - timer) << " seconds" << endl;
}


void topoHdl::analyze_sizegain_vs_matrix_bk_adhoc(unsigned int fd, std::vector<unsigned int>* survivors, std::vector<unsigned int>* matrix)
{
	//for all survivors identify in each time step how their size (area/volume) is relative to the mean of the matrix, 
	//the matrix are all grains that do not touch the boundary and are not included in the list of survivors
	//survivors - all the grains without boundary contact in the final data set, most likely recrystallizing
	//matrix - all grains ever considered with never boundary contact but excluding the survivors
	double timer = MPI_Wtime();

	if (survivors == NULL || matrix == NULL ) {
		cerr << "ERROR::Worker::AnalyzeSzGainVsMatrixBKADHOC " << this->get_Rank() << " filter input is non-existent or incomplete!" << endl;
		return;
	}

	//write well-defined sized dataset at known position cooperatively
	int ncol = ((Settings::SnapshotLast - Settings::SnapshotFirst) / Settings::SnapshotOffset + 1);
	int nrow = 2 + survivors->size(); //+2 for the number of grains in the mean value and the mean value itself

	double* relsizebuffer = NULL;
	try { relsizebuffer = new double[nrow]; }
	catch (std::bad_alloc &exc) {
		cerr << "ERROR::Worker::AnalyzeSzGainVsMatrixBKADHOC " << this->get_Rank() << " allocating memory relsizebuffer failed!" << endl;
		return;
	}

	for ( unsigned int s = 0; s < nrow; ++s ) 
		relsizebuffer[s] = 0.0;

	//process average grain size of matrix
	unsigned int localid = 0;
	struct descr_stats matrixstats;
	matrixstats = this->getDescrStats(matrix, localid, true);

	relsizebuffer[MATRIX_GRAINS_REMAINING] = matrixstats.n;
	relsizebuffer[MATRIX_GRAINS_MEAN] = matrixstats.mean;

	//process size gain in comparison to this matrix average
	if (matrixstats.mean > DBL_EPSILON) {
		for (unsigned int surv = 0; surv < survivors->size(); surv++) {
			unsigned int gid = survivors->at(surv);
			struct target_props targetp;
			targetp = this->getTargetProperties(gid, localid, true, false, false, false, false, false);

			relsizebuffer[2 + surv] = targetp.area / matrixstats.mean;
		}
	}

	string relsize_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".SIZEGAINBK.F." + std::to_string(Settings::SnapshotFirst) + ".O." + std::to_string(Settings::SnapshotOffset) + ".L." + std::to_string(Settings::SnapshotLast) + ".NC." + std::to_string(ncol) + ".NR." + std::to_string(nrow) + ".bin";
	unsigned int f = (fd / Settings::SnapshotOffset) - (Settings::SnapshotFirst / Settings::SnapshotOffset);
	long long totalOffsetREL = f * nrow * 8;
	
	MPI_File ioHdlREL;
	MPI_Status ioStaREL;

	MPI_File_open(MPI_COMM_SELF, relsize_fn.c_str(), MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlREL);
	MPI_File_write_at(ioHdlREL, totalOffsetREL, relsizebuffer, nrow, MPI_DOUBLE, &ioStaREL);
	MPI_File_close(&ioHdlREL);

	if (relsizebuffer != NULL) { delete[] relsizebuffer; relsizebuffer = NULL; }
	
	cout << "...Worker::AnalyzeSzGainVsMatrixBKADHOC " << this->get_Rank() << " analyzing SizeGainVsMatrix backwards = " << (MPI_Wtime() - timer) << " seconds" << endl;
}


void topoHdl::analyze_topodiff_fw_adhoc(unsigned int fd, std::vector<unsigned int>* thetargets)
{
	double timer = MPI_Wtime();

	if (thetargets == NULL) {
		cerr << "ERROR::Worker::AnalyzeTrackingTopologyDiffFWADHOC " << this->get_Rank() << " filter input is non-existent!" << endl;
		return;
	}

	//all workers have the same identically arranged list of target grain IDs because they performed the identical analysis routine beforehand
	int ncol = ((Settings::SnapshotLast - Settings::SnapshotFirst) / Settings::SnapshotOffset + 1);
	int nrow = thetargets->size();
	cout << "...Worker::AnalyzeTopologyDifferenceForward " << this->get_Rank() << " writing rows (ngrains) = " << nrow << " columns (time) = " << ncol << endl;

	//create I/O buffer container to write timeslot after timeslot, only as many write ops as timeslots because one node has all grains
	int* nfdiffbuffer = NULL;
	try { nfdiffbuffer = new int[nrow]; }
	catch (std::bad_alloc &exc) {
		cerr << "ERROR::Worker::AnalyzeTopologyDiffFWADHOC " << this->get_Rank() << " allocating memory for nfdiffbuffer failed!" << endl;
	}

	unsigned int f = (fd / Settings::SnapshotOffset) - (Settings::SnapshotFirst / Settings::SnapshotOffset);
	unsigned int localid = 0;

	for (unsigned int target = 0; target < thetargets->size(); target++) {
		unsigned int gid = thetargets->at(target); //assured to be a grain not at the boundary and not THE_DOMAIN_ITSELF

		struct target_props targetp;
		targetp = this->getTargetProperties(gid, localid, false, false, false, false, false, true);

		nfdiffbuffer[target] = targetp.nfacesk1diff;
	}

	//calculate address in already existent file where to store the results
	long long totalOffsetDiffFaces = f * nrow * 4;

	MPI_File ioHdlNFDiff;
	MPI_Status ioStaNFDiff;

	string suffix = "F." + std::to_string(Settings::SnapshotFirst) + ".O." + std::to_string(Settings::SnapshotOffset) + ".L." + std::to_string(Settings::SnapshotLast) + ".NC." + std::to_string(ncol) + ".NR." + std::to_string(nrow) + ".bin";
	string nfdiff_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".FW.NFDIFF." + suffix;

	//write to well-defined position in already existent file
	MPI_File_open(MPI_COMM_SELF, nfdiff_fn.c_str(), MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlNFDiff);
	MPI_File_write_at(ioHdlNFDiff, totalOffsetDiffFaces, nfdiffbuffer, nrow, MPI_INT, &ioStaNFDiff);
	MPI_File_close(&ioHdlNFDiff);


	if (nfdiffbuffer != NULL) { delete[]  nfdiffbuffer; nfdiffbuffer = NULL; }

	cout << "...Worker::AnalyzeTopologyDiffFWADHOC " << this->get_Rank() << " analyzing trajectories forward = " << (MPI_Wtime() - timer) << " seconds" << endl;
}


//void topoHdl::analyze_knn_naive_adhoc(unsigned int fd) {}
//void topoHdl::analyze_classical_nucmodels_adhoc(unsigned int fd) {}


spatQueryHdl::spatQueryHdl(void) 
{
	unsigned int ngr = 1 + Settings::LargestGrainID;

	ppp = NULL;
	try { ppp = new point3d[ngr]; }
	catch (std::bad_alloc &exc) {
		//##
	}
	for (unsigned int gid = 0; gid < ngr; ++gid) {
		ppp[gid].x = 0.0;
		ppp[gid].y = 0.0;
		ppp[gid].z = 0.0;
	}

	chi = NULL;
	try { chi = new double[ngr]; }
	catch (std::bad_alloc &exc) {
		//##
	}
	for (unsigned int gid = 0; gid < ngr; ++gid) {
		chi[gid] = std::numeric_limits<double>::lowest(); //##error value
	}

	mxsz = NULL;
	try { mxsz = new double[ngr]; }
	catch (std::bad_alloc &exc) {
		//##
	}
	for (unsigned int gid = 0; gid < ngr; ++gid) {
		mxsz[gid] = std::numeric_limits<double>::lowest(); //##error value
	}

	process = NULL;
	try { process = new bool[ngr]; }
	catch (std::bad_alloc &exc) {
		//##
	}
	for (unsigned int gid = 0; gid < ngr; ++gid) {
		process[gid] = false;
	}
	
	myRank = MASTER;
	nRanks = SINGLE_PROCESS;	
}


spatQueryHdl::~spatQueryHdl(void)
{
	for (unsigned int b = 0; b < pp3.size(); ++b) {
		if (pp3.at(b) != NULL) {
			delete pp3[b];
			//pp3[b] = NULL;
		}
	}
	pp3.clear();

	unsigned int ngr = 1 + Settings::LargestGrainID;
	for (unsigned int gid = 0; gid < ngr; ++gid ) {
		if ( nbors.at(gid) != NULL ) {
			delete nbors[gid]; 
			//nbors[gid] = NULL;
		}
	}
	nbors.clear();	

	if (ppp != NULL) {
		delete[] ppp;
		ppp = NULL;
	}
	if (chi != NULL) {
		delete[] chi;
		chi = NULL;
	}
	if (mxsz != NULL) {
		delete[] mxsz;
		mxsz = NULL;
	}
	if (process != NULL) {
		delete[] process;
		process = NULL;
	}
}


bool spatQueryHdl::init_spatquerybox(void)
{
	double ndbl = (1.0 / Settings::LongRangeRadiusMax) + 1.0 + DBL_EPSILON;
	rve = sqb(ndbl, ndbl, 1, 0); //unit square slab partitioned in cubic subdomains

	for (unsigned int b = 0; b < rve.nxyz; ++b) {
		std::vector<mpoint3d>* bucket = NULL;
		try { bucket = new std::vector<mpoint3d>; }
		catch (std::bad_alloc &exc) {
			cerr << "ERROR::Worker::InitQueryBox " << this->get_Rank() << " unable to allocate memory!" << endl;
			//##MK::dealloc memory for buckets
			return false;
		}
		pp3.push_back(bucket);
	}

	cout << "...Worker " << this->get_Rank() << " queryBox with nx " << rve.nx << " and nxyz " << pp3.size() << " build" << endl;
	return true;
}


inline unsigned int spatQueryHdl::pos2bucket(double x, double y, double z)
{
	double tmpd;
	unsigned int tmpui;
	unsigned int xx;
	unsigned int yy;
	unsigned int zz;

	//distribute on unit cube
	tmpd = rve.nx;
	tmpd *= x;
	tmpui = tmpd;
	xx = (tmpui != rve.nx) ? tmpui : tmpui - 1;

	tmpd = rve.ny;
	tmpd *= y;
	tmpui = tmpd;
	yy = (tmpui != rve.ny) ? tmpui : tmpui - 1;

	tmpd = rve.nz;
	tmpd *= z;
	tmpui = tmpd;
	zz = (tmpui != rve.nz) ? tmpui : tmpui - 1;

	return (xx + yy*rve.nx + zz*rve.nxyz);
}


bool spatQueryHdl::read_microstructure_uds(void)
{
	//reading Initial UDS File Microstructure.uds for identification of initial positions of sub-grains
	double tic = MPI_Wtime();
	ifstream udsfile;
	string udsline;
	istringstream line;
	string datapiece;

	udsfile.open(Settings::UDSDataFromFilename);
	if (udsfile.is_open() == false) {
		cerr << "ERROR::Worker " << this->get_Rank() << " unable to load file " << Settings::UDSDataFromFilename << endl;
		return false;
	}
	else { //jump over header first
		unsigned int nheaderlines = 1 + 1 + 2 + 2 + 1 + 1;
		unsigned int j = 0;
		while (udsfile.good() == true && j < nheaderlines) {
			getline(udsfile, udsline);
			j++;
		}

		unsigned int gid = 0;
		double x0 = 0.0;
		double y0 = 0.0;
		double z0 = 0.0;
		while (udsfile.good() == true) { //begin reading metadata
			getline(udsfile, udsline);

			if (udsline.size() > 0) { //not an empty line
				istringstream line(udsline);

				getline(line, datapiece, '\t');	gid = atoi(datapiece.c_str());
				getline(line, datapiece, '\t'); x0 = atof(datapiece.c_str()); //promotion of unsigned int to double no problemB
				getline(line, datapiece, '\t'); y0 = atof(datapiece.c_str());

				getline(line, datapiece); //discard the rest, implicitly reading up to '\n'

//cout << gid << "\t\t" << x0 << "\t\t" << y0 << endl;

				x0 = x0 / Settings::InitialDomainEdgeLength; //scale to unit cube
				y0 = y0 / Settings::InitialDomainEdgeLength;
				//z0 = z0 / Settings::InitialDomainEdgeLength;

				ppp[gid] = point3d( x0, y0, z0 );
				process[gid] = true;
			} //next line...
		} //...if still good
		
		udsfile.close();
	}
	double toc = MPI_Wtime();
	cout << "...Worker " << this->get_Rank() << " " << Settings::UDSDataFromFilename << " grain meta data were loaded successfully in " << (toc-tic) << " seconds" << endl;
	return true;
}

#define ARTIFICIAL_DRIVINGFORCE (1.0e6)

bool spatQueryHdl::read_longrange( unsigned int fd )
{
	//reading results of the longrange analysis into chi
	double tic = MPI_Wtime();
	
	string fn = "TopoTracer2D3D.SimID.509500.LR2D.FID." + std::to_string(fd) + ".RCLASS.0.csv"; 
	
	unsigned int ngr = 1 + Settings::LargestGrainID;
	for (unsigned int gid = 0; gid < ngr; ++gid) {
		chi[gid] = std::numeric_limits<double>::lowest(); //##error value
	}

	ifstream lrfile;
	string lrline;
	istringstream line;
	string datapiece;

	lrfile.open( fn.c_str() );
	if (lrfile.is_open() == false) {
		cerr << "ERROR::Worker " << this->get_Rank() << " unable to load file " << fn << endl;
		return false;
	}
	else { //jump over header first
		unsigned int nheaderlines = 15;
		unsigned int j = 0;
		while (lrfile.good() == true && j < nheaderlines) {
			getline(lrfile, lrline);
			j++;
		}

		unsigned int gid = 0;
		double mobenv = 0.0;

		while (lrfile.good() == true) { //begin reading metadata
			getline(lrfile, lrline);

			if (lrline.size() > 0) { //not an empty line
				istringstream line(lrline);

				getline(line, datapiece, ';');	gid = atoi(datapiece.c_str());
				getline(line, datapiece, ';');
				getline(line, datapiece, ';');

				getline(line, datapiece, ';');
				getline(line, datapiece, ';');
				getline(line, datapiece, ';');

				getline(line, datapiece, ';');
				getline(line, datapiece, ';');
				getline(line, datapiece, ';');

				getline(line, datapiece, ';');	mobenv = atof(datapiece.c_str());

				//cout << gid << "\t\t" << mobenv << endl;

				if (gid > THE_DOMAIN_ITSELF && gid < ngr && mobenv > 0.0) { //valid grain ID and valid value because >> std::numeric_limits<double>::lowest()
					chi[gid] = mobenv;
				}
				else {
					cerr << "ERROR::Worker " << this->get_Rank() << " invalid grain ID " << fn << endl;
					lrfile.close();
					return false;
				}
			} //next line...
		} //...if still good
		lrfile.close();
	}
	double toc = MPI_Wtime();
	cout << "...Worker " << this->get_Rank() << " " << fn.c_str() << " grain meta data were loaded successfully in " << (toc-tic) << " seconds" << endl;
	return true;
}


bool spatQueryHdl::init_fillspatquerybox(void)
{
	double tic = MPI_Wtime();

	//minus sampling method to avoid bias because of inspection areas protruding beyond nonperiodic RVE
	double minBound = 2.0*Settings::LongRangeRadiusMax; //METER2MICRON((2.0*Settings::LongRangeRadiusMax) * Settings::PhysicalDomainSize);
	double maxBound = 1.0 - 2.0*Settings::LongRangeRadiusMax; //METER2MICRON((1.0 - (2.0*Settings:LongRangeRadiusMax)) * Settings::PhysicalDomainSize);

	unsigned int ngr = 1 + Settings::LargestGrainID;
	for (unsigned int gid = 0; gid < ngr; ++gid) {
		if (process[gid] == true) {
			unsigned int b = pos2bucket(ppp[gid].x, ppp[gid].y, ppp[gid].z);

			bool consider = true;
			if ( ppp[gid].x > minBound && ppp[gid].x < maxBound && ppp[gid].y > minBound && ppp[gid].y < maxBound) { //##MK::z not probed
				pp3.at(b)->push_back(mpoint3d( ppp[gid].x, ppp[gid].y, ppp[gid].z, gid, false ));
			}
			else {
				pp3.at(b)->push_back(mpoint3d( ppp[gid].x, ppp[gid].y, ppp[gid].z, gid, true ));
			}			
		}
	}

	double toc = MPI_Wtime();
	cout << "...Worker " << this->get_Rank() << " pp3 filled in " << (toc - tic) << " seconds" << endl;
	return true;
}


bool spatQueryHdl::read_maxgrainsizes(string fn)
{
	//reading TopoTracer2D3D sim output for identification of maximum size of grain before ceasing
	double tic = MPI_Wtime();
	ifstream csvfile;
	string csvline;
	istringstream line;
	string datapiece;

	csvfile.open( fn.c_str() );
	if (csvfile.is_open() == false) {
		cerr << "ERROR::Worker " << this->get_Rank() << " unable to load file " << fn.c_str() << endl;
		return false;
	}
	else { //jump over header first
		unsigned int nheaderlines = 2;
		unsigned int j = 0;
		while (csvfile.good() == true && j < nheaderlines) {
			getline(csvfile, csvline);
			j++;
		}

		unsigned int gid = 0;
		double val = 0.0;
		while (csvfile.good() == true) { //begin reading metadata
			getline(csvfile, csvline);

			if (csvline.size() > 0) { //not an empty line
				istringstream line(csvline);

				getline(line, datapiece, ';');	gid = atoi(datapiece.c_str());
				getline(line, datapiece, ';');	val = atof(datapiece.c_str());
			
				//cout << gid << "\t\t" << val << endl;

				mxsz[gid] = val;
			} //next line...
		} //...if still good
		csvfile.close();
	}

	double toc = MPI_Wtime();
	cout << "...Worker " << this->get_Rank() << " " << fn.c_str() << " grain meta data were loaded successfully in " << (toc-tic) << " seconds" << endl;
	return true;
}


std::vector<sqkey>* spatQueryHdl::rangesearch(unsigned int gID, point3d p, double r)
{
	std::vector<sqkey>* res = NULL;
	try { res = new std::vector<sqkey>; }
	catch (std::bad_alloc &exc) {
		return NULL;
	}

	//normalized distances on unit cube
	double wxmi = p.x - r;
	if (wxmi < 0.0) wxmi = 0.0;
	double wxmx = p.x + r;
	if (wxmx > (1.0-DBL_EPSILON) ) wxmx = 1.0 - DBL_EPSILON; //MK::DBL_EPSILON must be large enough for unsigned int cast of 9.9999 to 9 in all cases!

	double wymi = p.y - r;
	if (wymi < 0.0) wymi = 0.0;
	double wymx = p.y + r;
	if (wymx > (1.0-DBL_EPSILON) ) wymx = 1.0 - DBL_EPSILON;

	double wzmi = p.z - r;
	if (wzmi < 0.0) wzmi = 0.0;
	double wzmx = p.z + r;
	if (wzmx > (1.0-DBL_EPSILON) ) wzmx = 1.0 - DBL_EPSILON;

	double scale_x = rve.nx;	wxmi *= scale_x;	wxmx *= scale_x;	//implict promotion of unsinged int to double is safe
	double scale_y = rve.ny;	wymi *= scale_y;	wymx *= scale_y;
	double scale_z = rve.nz;	wzmi *= scale_z;	wzmx *= scale_z;
	
	//query only point grids cells into which the sphere protrudes
	for (unsigned int zz = wzmi; zz <= wzmx; zz++) {
		unsigned int zoff = zz*rve.nxy;
		for (unsigned int yy = wymi; yy <= wymx; yy++) {
			unsigned int yzoff = yy*rve.nx + zoff;
			for (unsigned int xx = wxmi; xx <= wxmx; xx++) {
				unsigned int b = xx + yzoff;
				std::vector<mpoint3d>* bucket = pp3.at(b);
				unsigned int n = bucket->size();

				for (unsigned int i = 0; i < n; ++i) {
					//double uu = bucket->at(i).x - p.x;
					//double vv = bucket->at(i).y - p.y;
					//double ww = bucket->at(i).z - p.z;
					double dSQR = SQR(bucket->at(i).x-p.x) + SQR(bucket->at(i).y-p.y) + SQR(bucket->at(i).z-p.z);

					if ( dSQR <= SQR(r) ) { //inside spherical region!
						res->push_back( sqkey( sqrt(dSQR), bucket->at(i).gid, bucket->at(i).exclude ) );
					}
				} //next point bucket
			} //next xline of cells
		} //next xyslab of cells
	} //done checking intruded cells

	return res;
}

#define USE_MYOMP

void spatQueryHdl::init_neighbordistances(void)
{
	double tic = MPI_Wtime();
	//compute for all points all neighbors in circular/spherical region of radius R, memorize them in  all neighbroconduct analysis for grains for which we have a value
	unsigned int ngr = 1 + Settings::LargestGrainID;
	for (unsigned int gid = 0; gid < ngr; ++gid ) { //gid 2 all neighbors
		nbors.push_back( NULL );
	}

	unsigned long N = 0;

#ifdef USE_MYOMP
	#pragma omp parallel
	{
		unsigned int tid = omp_get_thread_num();
		unsigned int ntid = omp_get_num_threads();
		unsigned int ngr = 1 + Settings::LargestGrainID;
#else
		unsigned int tid = 0;
		unsigned int ntid = 1;
#endif

	
		double SearchRadius = Settings::LongRangeRadiusMax; //METER2MICRON * Settings::PhysicalDomainSize);
	
		unsigned int kk = 0;
		unsigned long elem = 0;

		for ( unsigned int b = 0; b < pp3.size(); ++b ) {
			if ( ((b+1) % ntid) == tid ) {
				std::vector<mpoint3d>* bucket = pp3.at(b);
				unsigned int n = bucket->size();
				for ( unsigned int i = 0; i < n; ++i ) {
					if ( bucket->at(i).exclude == false ) {
						//get all neighbors in circular region about point i employing exact euclidean distance in algorithm but using tree to reduce total checks
						std::vector<sqkey>* cand = NULL;
						//2-3ms
						cand = rangesearch( bucket->at(i).gid, point3d(bucket->at(i).x, bucket->at(i).y, bucket->at(i).z), SearchRadius);
				
						std::sort(cand->begin(), cand->end(), SortSQKeyDistAscending);

						//##MK::now comparable to Matlab2015a rangesearch with prebuild kdtree, 
						//cout << bucket->at(i).gid << "\t\t" << cand->size() << endl;
						//kk++; if ( kk % 1000 != 0 ) continue; cout << kk << endl;
						//elem += cand->size();

						if ( cand->size() > 0 ) {
#ifdef USE_MYOMP
							#pragma omp critical
							{
#endif
								nbors.at( bucket->at(i).gid ) = cand; //MK::threads write to disjoint positions in nbors!
#ifdef USE_MYOMP
							}
#endif
						}
						else { //allocated but no neighbor found
							if (cand != NULL) { delete cand; cand = NULL; }
						}
				
					} //next point i
				} //next bucket b
			}
		}

#ifdef USE_MYOMP
		#pragma omp critical
		{
			N += elem;
		}
	}
#endif

	cout << N << endl;

	double toc = MPI_Wtime();
	cout << "...Worker " << this->get_Rank() << " envcharacterize in " << (toc - tic) << " seconds" << endl;
} 


double spatQueryHdl::quantile_asc_sorted_dbl( std::vector<double>* dat, double qnt )
{
	//computes qnt's quantile by utilizing linear interpolation with ascendingly sorted list of size > 0
	//MK::only finite number of distance values and interpolation makes no sense because there are no grains in between!
	/*
	double qn = qnt * (double) dat->size();
	unsigned int whom = qn;
	if ( whom != dat->size() ) 
		return dat->at(whom);
	else
		return dat->at(dat->size()-1);
	*/
		
	//quick check for extrema
	if ( qnt < DBL_EPSILON ) //min 
		return dat->at(0);
	if ( qnt > (1.0-DBL_EPSILON) ) //and max
		return dat->at(dat->size()-1);
	
	//exhaustive using linear interpolation
	double n = dat->size();
	double dn = 1.0 / ((double) n);
	double cdf = 0.0;
	double xl, xr, yl, yr;
	double res = std::numeric_limits<double>::lowest();

	for ( unsigned int i = 0; i < n; i++ ) {
		cdf += dn;
		if ( cdf < qnt )
			continue;

		//not continued so cdf >= qnt
		//stop and interpolate linearly
		unsigned int il = (i != 0) ? i-1 : 0;
		unsigned int ir = i;
		xl = dat->at(il);
		xr = dat->at(ir);
		yl = 1+il;			yl *= dn;
		yr = 1+ir;			yr *= dn;

		res = xl + (qnt-yl)/(yr-yl) * (xr-xl);
		break;
	}
	return res;
}

#define QUANT_MIN		1
#define QUANT_MAX		99
#define QUANT_WIDTH		99-1+1

void spatQueryHdl::arrivaltime_quantiles( unsigned int fd )
{
	double tic = MPI_Wtime();

	//for all neighbors with measured neighbors compute arrival time distribution
	std::vector<arrivaltime_results>* globalbuffer = NULL;
	globalbuffer = new std::vector<arrivaltime_results>;
	
	unsigned int ngr = 1 + Settings::LargestGrainID;
	for ( unsigned int gid = 0; gid < ngr; ++gid ) {
		globalbuffer->push_back( arrivaltime_results() );
	}

	#pragma omp parallel shared(ngr)
	{
		unsigned int tid = omp_get_thread_num();
		unsigned int ntid = omp_get_num_threads();

		std::vector<arrivaltime_results>* localbuffer = NULL;
		localbuffer = new std::vector<arrivaltime_results>;

		std::vector<double>* tau = NULL;
		tau = new vector<double>;

		std::vector<sqkey>* nb = NULL;
		
		for ( unsigned int gid = 1; gid < ngr; ++gid ) {
			if ( gid % ntid == tid ) { //my work
				nb = nbors.at(gid);
				if ( nb != NULL && process[gid] == true ) { //only for grains for which neighbors were found neighbors found
					tau->clear();
			
					double gid_speed = ARTIFICIAL_DRIVINGFORCE * chi[gid];
					if ( gid_speed > 0.0 ) {
						double nb_speed = std::numeric_limits<double>::lowest();
			
						unsigned int nc = nb->size();
						for ( unsigned int c = 1; c < nc; ++c ) { //start at one to jump over myself
							nb_speed = ARTIFICIAL_DRIVINGFORCE * chi[nb->at(c).pid];
							if ( nb_speed > 0.0 ) { //valid value
								//model 1 choose difference
								//double arrdiff1 = (nb->at(c).dist / gid_speed) - (nb->at(c).dist / nb_speed); //difference in arrival time
						
								//model 2 choose quotient
								double arrdiff2 = (nb->at(c).dist / gid_speed) / (nb->at(c).dist / nb_speed);
								//negative means gid has a growth advantage
						
								tau->push_back( arrdiff2 );
							}
							//invalid chi means most likely grain died between the previous and this timestep
						} //next candidate

						if ( tau->size() > 0 ) { //at least any values found?
							std::sort( tau->begin(), tau->end() );

							//integration "Matano-Boltzmann style"
							unsigned int nt = tau->size();
							double dy = 1.0/((double) nt);
							double integralval = 0.0;
							for ( unsigned int t = 0; t < nt; ++t ) {
								integralval += dy * tau->at(t);
							}
								
							localbuffer->push_back( arrivaltime_results( gid, tau->size(), integralval, mxsz[gid] ) );
							//logfile << gid << ";" << tau->size() << ";" << integralval << ";" << mxsz[gid] << "\n";
						}
					} //next grain to process
				}
			} //next grain for me
		}
		
		#pragma omp critical
		{
			for ( unsigned int i = 0; i < localbuffer->size(); ++i ) {
				unsigned int gID = localbuffer->at(i).gid;
				globalbuffer->at(gID) = localbuffer->at(i);
			}
		}

		if ( localbuffer != NULL ) { delete localbuffer; localbuffer = NULL; }
		if ( tau != NULL ) { delete tau; tau = NULL; }
		nb = NULL;
	} //end of OMP parallel region
	
	string logfname = "TopoTracer2D3D.SimID.509500.Competitiveness.FID." + std::to_string(fd) + ".csv";
	ofstream logfile;
	logfile.open( logfname.c_str(), std::ofstream::out | std::ofstream::trunc );
	if ( logfile.is_open() == false ) {
		cerr << "ERROR::Worker " << this->get_Rank() << " unable to open " << logfname.c_str() << endl;
		return;
	}
	logfile << "GrainID;NeighborsProbed;IntegralValue;MaxSizeObtainedEverMicron^Dim\n";
	for ( unsigned int i = 0; i < globalbuffer->size(); ++i ) {
		if ( globalbuffer->at(i).gid != THE_DOMAIN_ITSELF ) {
			logfile << globalbuffer->at(i).gid << ";" << globalbuffer->at(i).nt << ";" << globalbuffer->at(i).intval << ";" << globalbuffer->at(i).mxsize << "\n";
		}
	}
	logfile.flush();
	logfile.close();

	if ( globalbuffer != NULL ) { delete globalbuffer; globalbuffer = NULL; }
	
	double toc = MPI_Wtime();
	cout << "...Worker " << this->get_Rank() << " characterize arrival time in " << (toc - tic) << " seconds" << endl;
}


utilityHdl::utilityHdl( void ) {
	//kkmin implicitly 0
	kkmax = Settings::MaxNumberOfKShells;
	nRanks = SINGLE_PROCESS;
	myRank = MASTER;
}


utilityHdl::~utilityHdl( void ) {
	for ( unsigned int luid = 0; luid < KNNDatabase.size(); luid++  ) {
		if ( KNNDatabase[luid] != NULL ) { //when IDs in this range where there
			for ( std::vector<knn_meta>::iterator it = KNNDatabase[luid]->begin(); it != KNNDatabase[luid]->end(); it++ ) {
				if (it->GrainsInShell != NULL) { delete[] it->GrainsInShell;	it->GrainsInShell = NULL; }
				if (it->AreaOfKShell != NULL) { delete[] it->AreaOfKShell;		it->AreaOfKShell = NULL; }
				if (it->BndLen != NULL) { delete[] it->BndLen;			it->BndLen = NULL; }
				if (it->HAGB != NULL) { delete[] it->HAGB;				it->HAGB = NULL; }
				if (it->MdSEE != NULL) { delete[] it->MdSEE;			it->MdSEE = NULL; }
			}
			delete KNNDatabase[luid];
			KNNDatabase[luid] = NULL;
		}
	}
}


void utilityHdl::set_MPICommunication( int r, int nr ) {
	myRank = r;
	nRanks = nr;
}


int utilityHdl::get_Rank( void ) {
	return myRank;
}


int utilityHdl::get_nRanks( void ) {
	return nRanks;
}


bool utilityHdl::init_KNNDatabase( void ) {

	unsigned int lu_size = 1 + (Settings::LargestGrainID / Settings::LookupMaxGrainIDRange);
	//initialize KNNDatabase with NULL first
	for ( unsigned int luid = 0; luid < lu_size; luid++ ) {
		KNNDatabase.push_back( NULL );
	}

	for ( unsigned int luid = 0; luid < lu_size; luid++ ) {
		std::vector<knn_meta>* knnmeta_bucket = NULL;
		try { knnmeta_bucket = new std::vector<knn_meta>; }
		catch (std::bad_alloc &exc) {
			cerr << "ERROR::Initializing KNNDatabase!" << endl;
			for (unsigned int j = 0; j < luid; ++j) {
				if (KNNDatabase[j] != NULL) {
					delete KNNDatabase[j]; KNNDatabase[j] = NULL;
				}
			}
			return false;
		}
		
		KNNDatabase[luid] = knnmeta_bucket; //store pointer to database entries
	}

	return true;
}


bool utilityHdl::read_targetfile( void ) {
	double timer = MPI_Wtime();

	TargetGIDs.clear();

	ifstream trgfile;
	string trgline;
	istringstream line;
	string datapiece;

	trgfile.open( Settings::TargetGIDFromFilename );
	if ( trgfile.is_open() == true ) {
		//has no header
		//getline( trgfile, trgline );

		while ( trgfile.good() == true ) {
			getline( trgfile, trgline );
			unsigned int s = trgline.size();

			if ( s > 0 ) { //only a target id 
				unsigned int tgid = std::stoul( trgline ); //potential overflow
				TargetGIDs.push_back( tgid );
//cout << tgid << "\t\t" << TargetGIDs[TargetGIDs.size()-1] << endl;
			}
		}
		trgfile.close();
	}
	
	cout << TargetGIDs.size() << " targets were loaded successfully in " << (double) (MPI_Wtime() - timer) << " seconds" << endl;

	return true;
}


bool utilityHdl::read_knnfile( unsigned int kkmax ) {
	double timer = MPI_Wtime();

	ifstream knnfile;
	string knnline;
	istringstream line;
	string datapiece;

	knnfile.open( Settings::KNNDataFromFilename );
	if ( knnfile.is_open() == true ) {
		//kick header
		getline( knnfile, knnline );

		if ( init_KNNDatabase() == false ) {
			return false;
		}

		//##MK::analyze header structure to catch if it matches with Settings::maximum i.e. kkmax

		//read data
		bool memorygood = true;
		while ( knnfile.good() == true ) {
			getline( knnfile, knnline );
			if ( knnline.size() > 0 ) {
				istringstream line( knnline );

				struct knn_meta ak;
				
				getline(line, datapiece, ';');	ak.gid = atoi( datapiece.c_str() );
				getline(line, datapiece, ';');	ak.x = atof( datapiece.c_str() );
				getline(line, datapiece, ';');	ak.y = atof( datapiece.c_str() );

				ak.kshellbased = kkmax + 1;
				ak.peribased = kkmax; //##MK

				try { ak.GrainsInShell = new unsigned int[ak.kshellbased]; }
				catch (std::bad_alloc &exc) { memorygood = false; break; }
				try { ak.AreaOfKShell = new double[ak.kshellbased]; }
				catch (std::bad_alloc &exc) { memorygood = false; break; }
				try { ak.BndLen = new double[ak.peribased]; }
				catch (std::bad_alloc &exc) { memorygood = false; break; }
				try { ak.HAGB = new double[ak.peribased]; }
				catch (std::bad_alloc &exc) { memorygood = false; break; }
				try { ak.MdSEE = new double[ak.peribased]; }
				catch (std::bad_alloc &exc) { memorygood = false; break; }

				//memory guaranteed accessible
				for ( unsigned int kk = 0; kk <= kkmax; kk++ ) {
					getline(line, datapiece, ';');		ak.GrainsInShell[kk] = atoi( datapiece.c_str() );
				}
				for ( unsigned int kk = 0; kk <= kkmax; kk++ ) {
					getline(line, datapiece, ';');		ak.AreaOfKShell[kk] = atof( datapiece.c_str() );
				}
				for ( unsigned int kk = 0; kk < kkmax; kk++ ) {
					getline(line, datapiece, ';');		ak.BndLen[kk] = atof( datapiece.c_str() );
				}
				for ( unsigned int kk = 0; kk < kkmax; kk++ ) {
					getline(line, datapiece, ';');		ak.HAGB[kk] = atof( datapiece.c_str() );
				}
				for ( unsigned int kk = 0; kk < kkmax; kk++ ) {
					getline(line, datapiece, ';');		ak.MdSEE[kk] = atof( datapiece.c_str() );
				}

				//in which luid
				unsigned int luid = ak.gid / Settings::LookupMaxGrainIDRange;

				KNNDatabase[luid]->push_back( ak );
//cout << KNNDatabase[KNNDatabase.size()-1].gid << ";" << KNNDatabase[KNNDatabase.size()-1].GrainsInShell[kkmax] << ";" << KNNDatabase[KNNDatabase.size()-1].AreaOfKShell[kkmax] << ";" << KNNDatabase[KNNDatabase.size()-1].BndLen[kkmax-1] << ";" << KNNDatabase[KNNDatabase.size()-1].HAGB[kkmax-1] << ";" << KNNDatabase[KNNDatabase.size()-1].MdSEE[kkmax-1] << endl;
			}
		}

		//if (memorygood == false) destructor of utilityHdl deletes potentially already allocated memory

		knnfile.close();
	}

	//##MK::consistency check//for ( unsigned int luid = 0; luid < KNNDatabase.size(); luid++ ) { cout << luid << ";" << KNNDatabase[luid]->size(); }

cout << "The KNNDatabase with " << KNNDatabase.size() << " elements was loaded in " << (double) (MPI_Wtime() - timer) << " seconds" << endl;

	return true;
}


void utilityHdl::analyze_apriori_prediction( void ) {
	//MK::naive, only one process, compute from KNNDatabase the sum(HAGB)/Nk and sum(mdeltarho)/Nk
	//to allow the correlation of the actual size gain of a grain to the initial long-range environment
	if ( kkmax < 1 ) { cerr << "Set kkmax at least = 1 !" << endl; return; }

double timer = MPI_Wtime();

	string log_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".APrioriLongRangeEnvironment.csv";
	ofstream logpredict;
	logpredict.open ( log_fn.c_str() );
	if (logpredict.is_open() == true) {

		logpredict << "GrainID;Sum(HAGB)/Nk;Sum(MOBDSEE)/Nk" << "\n";

		double _kkmax = 1.0 / ((double)kkmax);

		for (unsigned int i = 0; i < KNNDatabase.size(); i++) {
			if (KNNDatabase.at(i) != NULL) {
				for (unsigned int j = 0; j < KNNDatabase.at(i)->size(); j++) {
					double sumhagb = 0.0;
					for (unsigned int kk = 0; kk < kkmax; kk++) 
						sumhagb += KNNDatabase.at(i)->at(j).HAGB[kk];

					double summob = 0.0;
					for (unsigned int kk = 0; kk < kkmax; kk++) 
						summob += KNNDatabase.at(i)->at(j).MdSEE[kk];

					sumhagb = sumhagb * _kkmax;
					summob = summob * _kkmax;
					logpredict << KNNDatabase.at(i)->at(j).gid << ";" << sumhagb << ";" << summob << "\n";
				}
			}
		}
		logpredict.flush();
		logpredict.close();
	}
	
	cout << "A priori predictions were madein " << (double) (MPI_Wtime() - timer) << " seconds" << endl;
}


void utilityHdl::analyze_find_knn_for_targets( void ) {

	double timer = MPI_Wtime();

	if ( this->TargetGIDs.size() > 0 ) {
		//generate an outputfile ...
		string log_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".DataAssgnKNNforTargets.csv";
		ofstream logknntargets;
		logknntargets.open ( log_fn.c_str(), std::ofstream::out | std::ofstream::trunc );
		if (logknntargets.is_open() == true) {
			//... with a header
			logknntargets << "GrainID;x(AtKNNTimeStep);y(AtKNNTimeStep);";
			for (unsigned int kj = 0; kj <= kkmax; kj++) logknntargets << "GrainsIn-" << kj << "-thShell" << ";";
			for (unsigned int kj = 0; kj <= kkmax; kj++) logknntargets << "SumArea-" << kj << "-thShell" << ";";
			for (unsigned int kj = 0; kj < kkmax; kj++) logknntargets << "BndLength-" << kj << "-thShell" << ";";
			for (unsigned int kj = 0; kj < kkmax; kj++) logknntargets << "HAGBFraction-" << kj << "-thShell" << ";";
			for (unsigned int kj = 0; kj < kkmax; kj++) logknntargets << "MdSEEwBndLength-" << kj << "-thShell" << ";";
			logknntargets << "\n";

			//##MK::find the targets, relatively stupidly
			bool found = false;
			for (unsigned int t = 0; t < this->TargetGIDs.size(); t++) {
				unsigned int tgrid = TargetGIDs[t];
				unsigned int luidt = tgrid / Settings::LookupMaxGrainIDRange;
				found = false;
				for (unsigned int i = 0; i < KNNDatabase[luidt]->size(); i++) {
					if (tgrid != KNNDatabase[luidt]->at(i).gid) {
						continue;
					}
					//found because not continued, so IO
					found = true;
					logknntargets << KNNDatabase[luidt]->at(i).gid << ";" << KNNDatabase[luidt]->at(i).x << ";" << KNNDatabase[luidt]->at(i).y << ";";
					for (unsigned int kj = 0; kj <= kkmax; kj++) logknntargets << KNNDatabase[luidt]->at(i).GrainsInShell[kj] << ";";
					for (unsigned int kj = 0; kj <= kkmax; kj++) logknntargets << KNNDatabase[luidt]->at(i).AreaOfKShell[kj] << ";";
					for (unsigned int kj = 0; kj < kkmax; kj++) logknntargets << KNNDatabase[luidt]->at(i).BndLen[kj] << ";";
					for (unsigned int kj = 0; kj < kkmax; kj++) logknntargets << KNNDatabase[luidt]->at(i).HAGB[kj] << ";";
					for (unsigned int kj = 0; kj < kkmax; kj++) logknntargets << KNNDatabase[luidt]->at(i).MdSEE[kj] << ";";
					logknntargets << "\n";
					//cout << "tgrid/gid = " << tgrid << ";" << KNNDatabase[luidt]->at(i).gid << endl;
					break;
				}
				if (found == false) {
					cerr << "Target " << tgrid << " was not found probably because the KNN-shell extends long and hence many grains close to the domain boundary" << endl;
					cerr << " were untrackable as some of the neighbors in a shell touched the boundary. In the coarsening sim still though such a target may consume only" << endl;
					cerr << " a few shells prior to the point in time at which a BK tracking analysis was made. Hence the error!" << endl;
				} //annotate error but continue to find other grains...
			} //for all targets

			logknntargets.flush();
			logknntargets.close();
		}
		else {
			cerr << "ERROR::Worker::AnalyzeFindKNN " << this->get_Rank() << " unable to open logknntargets file!" << endl;
		}
	}

cout << "The targets were tracked in " << (double) (MPI_Wtime() - timer) << " seconds" << endl;
}


void utilityHdl::analyze_Gest_for_targets( void ) {
	double timer = MPI_Wtime();

	unsigned int nt = TargetGIDs.size();

	if (nt == 0) 
		return;
		
	//build size hash
	std::vector<point3d> tpp;
	for ( unsigned int t = 0; t < nt; t++ ) {
		unsigned int tgrid = TargetGIDs[t];
		unsigned int luidt = tgrid / Settings::LookupMaxGrainIDRange;
		bool found = false;
		for ( unsigned int i = 0; i < KNNDatabase[luidt]->size(); ++i ) {
			if ( tgrid != KNNDatabase[luidt]->at(i).gid ) {
				continue;
			}
			//not continued so found
			found = true;
			//##MK::fill z coordinate with life in 3D implementation
			tpp.push_back( point3d(KNNDatabase[luidt]->at(i).x, KNNDatabase[luidt]->at(i).y, 0.0) );
			break;
		}
		if ( found == false ) {
			cerr << "Target " << tgrid << " was not found probably because the KNN-shell extends long and hence many grains close to the domain boundary" << endl;
			cerr << " were untrackable as some of the neighbors in a shell touched the boundary. In the coarsening sim still though such a target may consume only" << endl; 
			cerr << " a few shells prior to the point in time at which a BK tracking analysis was made. Hence the error!" << endl; 
		} //claim that the list of KNN valid targets for the G2 analysis is incomplete
	}

	if ( tpp.size() == TargetGIDs.size() )
		cout << "The KNN environment for all targets was found estimating now with " << tpp.size() << " targets." << endl;
	else {
		cerr << "The KNN environment was NOT FOUND for all targets, hence stopping" << endl;
		return;
	}

	//generate an outputfile ...
	string log_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".NearestNeighborAnalysis.csv";
	ofstream logg2targets;
	logg2targets.open ( log_fn.c_str(), std::ofstream::out | std::ofstream::trunc );

	if (logg2targets.is_open() == true) {
		logg2targets << "TargetGID;x(AtKNNTimeStep);y(AtKNNTimeStep);z(AtKNNTimeStep);NearestNBorDist;NearestNBorGID\n";

		for (unsigned int t = 0; t < nt; t++) {
			double tx = tpp[t].x;
			double ty = tpp[t].y;
			double tz = tpp[t].z;
			double mindist = 2.0 * SQR(pow(3.0, 0.5) * 1.0); //twice room SQR of room diagonal in unit cube //##MK CUBE(2 * 1.0 * pow(3.0, 0.5));
			unsigned int minid = nt + 1;

			for (unsigned int tl = 0; tl < t; tl++) { //relatively naive check all neighbors in domain, but split in two loops each to avoid comparison against oneselv
				double dist = SQR(tpp[tl].x - tx) + SQR(tpp[tl].y - ty); //##MK::add when truly in 3d, i.e. tz != 0.0 +SQR(tpp[tl].z - tz);
				if (dist > mindist) //most likely farther away
					continue;
				//<= mindist
				mindist = dist;
				minid = tl;
			}
			for (unsigned int tu = t + 1; tu < nt; tu++) { //two for loops to avoid finding myself as being closest to myself
				double dist = SQR(tpp[tu].x - tx) + SQR(tpp[tu].y - ty); //##MK::see above +SQR(tpp[tu].z - tz);
				if (dist > mindist) 
					continue;
				//<= mindist
				mindist = dist;
				minid = tu;
			}

			logg2targets << TargetGIDs[t] << ";" << tx << ";" << ty << ";" << tz << ";" << mindist << ";" << TargetGIDs[minid] << "\n";
		}
		logg2targets.flush();
		logg2targets.close();
	}

	cout << "The nearest neighbor distribution was found in " << (double) (MPI_Wtime() - timer) << " seconds" << endl;
}


//#ifdef USE_POLYTRI

arealanaHdl::arealanaHdl( void ) 
{
}


arealanaHdl::~arealanaHdl( void )
{
	for ( unsigned int c = 0; c < gcontour2d.size(); c++  ) {
		if ( gcontour2d[c] != NULL ) { //when IDs in this range where there
			delete gcontour2d[c];
			gcontour2d[c] = NULL;
		}

	}

	//##MK::memory leak?
	for ( unsigned int c = 0; c < gsegments2d.size(); c++ ) {
		if ( gsegments2d[c] != NULL ) {
			delete gsegments2d[c];
			gsegments2d[c] = NULL;
		}
	}

#ifdef COMPILE_IN_IF_DESIRED
	for ( unsigned int t = 0; t < trial.size(); t++ ) {
		if ( trial[t] != NULL ) {
			delete trial[t];
			trial[t] = NULL;
		}
	}
#endif

	for ( unsigned int rclass = 0; rclass < ChiTable.size(); rclass++ ) {
		if ( ChiTable[rclass] != NULL ) {
			delete ChiTable[rclass];
			ChiTable[rclass] = NULL;
		}
	}
}


void arealanaHdl::set_MPICommunication( int r, int nr ) {
	myRank = r;
	nRanks = nr;
}


void arealanaHdl::init_MPIDatatypes( void )
{
	//MPI_Datatype MPI_GBContourPoint_Type;
	int elementCounts1[2] = {4, 2};
	MPI_Aint displacements1[2] = {0, 4 * MPIIO_OFFSET_INCR_DOUBLE};
	MPI_Datatype oldTypes1[2] = {MPI_DOUBLE, MPI_UNSIGNED};
	MPI_Type_create_struct(2, elementCounts1, displacements1, oldTypes1, &MPI_GBContourPoint_Type);


	MPI_Type_commit(&MPI_GBContourPoint_Type);
}


int arealanaHdl::get_Rank( void ) {
	return myRank;
}


int arealanaHdl::get_nRanks( void ) {
	return nRanks;
}


bool arealanaHdl::read_microstructure_uds_2d( void ) {
	//reading Initial UDS File Microstructure.uds for computation of unbiased growth prognosis, we need this file to get to know 
	//all meta data to the grains whihc the GBContourPoints_<fid>.bin files do not contain
	double timer = MPI_Wtime();
	ifstream udsfile;
	string udsline;
	istringstream line;
	string datapiece;

	struct isectionmeta dat;
	
	udsfile.open( Settings::UDSDataFromFilename );
	if ( udsfile.is_open() == true ) { //read in file line-by-line
		//ump over header first
		unsigned int nheaderlines = 1 + 1 + 2 + 2 + 1 + 1;
		unsigned int j = 0;
		while ( udsfile.good() == true && j < nheaderlines ) {
			getline( udsfile, udsline);
			j++;
		}

		while ( udsfile.good() == true ) { //begin reading metadata
			getline( udsfile, udsline );

			if ( udsline.size() > 0 ) { //not an empty line
				istringstream line( udsline );
				//interpret
				getline(line, datapiece, '\t');	dat.gid = atoi( datapiece.c_str() );
				getline(line, datapiece, '\t'); dat.x0 = atoi( datapiece.c_str() ); //initial reference position as center of AABB
				getline(line, datapiece, '\t'); dat.y0 = atoi( datapiece.c_str() );
				
				double bunge[3];
				getline(line, datapiece, '\t'); bunge[0] = atof( datapiece.c_str() );
				getline(line, datapiece, '\t'); bunge[1] = atof( datapiece.c_str() );
				getline(line, datapiece, '\t'); bunge[2] = atof( datapiece.c_str() );
				double quat[4];
				euler2quaternion( bunge, quat );
				dat.q0 = quat[0];		dat.q1 = quat[1];		dat.q2 = quat[2];		dat.q3 = quat[3];

				getline(line, datapiece, '\t'); //xmin, xmax, ymin, ymax
				getline(line, datapiece, '\t');
				getline(line, datapiece, '\t');
				getline(line, datapiece, '\t');
				getline(line, datapiece, '\t'); //vol
				getline(line, datapiece, '\t'); dat.see = atof( datapiece.c_str() ); //the 2d version of the GraGLeS solver outputs directly dislocation in unit 1/m^2
				
				/*
				//##MK::when potentially reading from a MicrostructureDiagnostics.uds use this code to eliminate the trailing columns of additional grain-related details
				//eliminate the rest 10 columns, rather crappy design, i know ##MK
				getline(line, datapiece, '\t');
				getline(line, datapiece, '\t');
				getline(line, datapiece, '\t');
				getline(line, datapiece, '\t');
				getline(line, datapiece, '\t');
				getline(line, datapiece, '\t');
				getline(line, datapiece, '\t');
				getline(line, datapiece, '\t');
				getline(line, datapiece, '\t');
				getline(line, datapiece, '\t');
				*/
				
				//cout << dat.gid << "\t" << bunge[0] << ";" << bunge[1] << ";" << bunge[2] << "\t" << dat.q0 << ";" << dat.q1 << ";" << dat.q2 << ";" << dat.q3 << ";" << dat.see << endl;
				//cout << dat.gid << ";" << dat.see << endl;
				//invokes copyconstructor for stl types
				if ( gmetadata.at(dat.gid).gid == THE_DOMAIN_ITSELF ) { //most likely case
					gmetadata.at(dat.gid) = dat;
				} 
				else { 
					cerr << "ERROR::Worker " << this->get_Rank() << " inconsistency namely the same ID twice!" << endl;
					udsfile.close(); return false;
				}
			} //next line...
		} //...if still good
		udsfile.close();
	}
	else { 
		cerr << "ERROR::Worker " << this->get_Rank() << " unable to load file " << Settings::UDSDataFromFilename << endl;
		return false;
	}

	myprofiler.logev( "ReadMetadataInitially", (double) (MPI_Wtime() - timer) );
	cout << "...Worker " << this->get_Rank() << " " << Settings::UDSDataFromFilename << " grain meta data were loaded successfully in " << (double) (MPI_Wtime() - timer) << " seconds" << endl;
	return true;
}


void arealanaHdl::read_initialmicrostructure_metadata( void )
{
	//MK::reads a Microstructure_<fd>.uds ASCII file to get the mapping of grainIDs to grain properties
	if ( read_microstructure_uds_2d() == false ) { cerr << "ERROR::Unable to locate Microstructure.uds file" << endl; return; }
}


bool arealanaHdl::read_snapshot_2d( unsigned int fd )
{
	//read MPI I/O contours and distribute on available grains via gcontour2d
	double timer = MPI_Wtime();

	///MPI_COMM_SELF read of an implicit binary of GBContourPoint structs
	string fn = "GBContourPoints_" + std::to_string(fd)  + ".bin";
	
	//probe file size of fn to check if consistent
	struct stat buf;
	double filesize[1] = {0.0};
	unsigned long fsz = 0;
	unsigned long szcp = sizeof(struct GBContourPoint);
	if ( stat( fn.c_str(), &buf ) != -1 ) {
		filesize[0] = buf.st_size;

		fsz = filesize[0];
		if ( fsz % szcp != 0 ) { //file size is not a multiple of sizeof GBContourPoint, indicates file is potentially corrupted!
			cerr << "ERROR::Worker " << this->get_Rank() << " File " << fn << " potentially corrupted!" << endl;
			return false;
		}
	}
	else { cerr << "ERROR::Worker " << this->get_Rank() << " file size of " << fn.c_str() << " not detectable!" << endl; 
		return false;
	}

	//file exists, size is known so read it in as contiguous blocks
	//MK::MPI library only read blocks of maximum size at once but contour files can be larger than 4GB
	//no global buffering of the file content required, MK::extremely beneficial because reducing memory overhead...
	MPI_File ioReadFileHdl;
	MPI_Status ioReadFileStatus;

	MPI_File_open(MPI_COMM_SELF, fn.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &ioReadFileHdl);
	//native view, MPI_BYTE
	MPI_File_seek(ioReadFileHdl, 0, MPI_SEEK_SET); //set implicit fp to beginning of file

	bool stillgood = true;
	unsigned long elementsRead = 0; //MK::not unsigned int because file potentially more than UINT32MAX elements
	unsigned long elementsTotal = fsz / szcp; //MK::safe because we have assured fsz to be an integer multiple of szcp
	unsigned int elementsNow = 0; //read in blocks of struct GBContourPoint begin at offset=0, MPI_CHAR native interpretation
	unsigned int elementsPerBlock = Settings::MPIReadBlockLength / sizeof(struct GBContourPoint);

	while ( (elementsRead < elementsTotal) && stillgood == true ) {
		//formulate a single block read via MPI i/O
		elementsNow = elementsPerBlock; //MK::not unsigned long , compatibility with MPI library

		 //handle remaining potentially shorter rest of the data
		if ( (elementsTotal - elementsRead) < ((unsigned long) elementsNow) ) {
			unsigned long diff = elementsTotal - elementsRead;
			if ( diff > std::numeric_limits<unsigned int>::max() ) {
				cerr << "ERROR::Worker " << this->get_Rank() << " remaining block handling on " << fn.c_str() << " failed!" << endl;
				stillgood = false;
				break; //MK::do not return false, because then ioReradFileHdl is not closed!
			}
			elementsNow = (unsigned int) diff; //MK::downcasting from unsigned long to unsigned int is now safe
		}

		if ( stillgood == true ) {
			//set up a buffer to store MPI I/O result
			GBContourPoint* rbuf = NULL;
			try { rbuf = new GBContourPoint[elementsNow]; }
			catch (std::bad_alloc &exc) { 
				cerr << "ERROR::Worker " << this->get_Rank() << " allocation error in read_snapshot_2d_mpi_gragles!" << endl; 
				MPI_File_close(&ioReadFileHdl); return false;
			}

			MPI_File_read( ioReadFileHdl, rbuf, elementsNow, MPI_GBContourPoint_Type, &ioReadFileStatus);
			//implicitly advancing fp

			//interpret data immediately to gcontour2d
			unsigned gid = 0;
			for ( unsigned int e = 0; e < elementsNow; ++e ) {
				gid = rbuf[e].myid;
				if ( gid != THE_DOMAIN_ITSELF && gid <= Settings::LargestGrainID ) { //##MK::debug safety to catch inconsistencies
					//this loading imprints an explicit order on the grain-local container preserving the order of GraGLeS
					if ( gcontour2d[gid] != NULL ) { //not the first time to initialize contourpoints for this grain
						gcontour2d[gid]->push_back( point2d(rbuf[e].x, rbuf[e].y) );
						continue;
					}
					//least likely case first time to add a contour point
					std::vector<point2d>* contourpoints = NULL;
					try { contourpoints = new vector<point2d>; }
					catch (std::bad_alloc &exc) { 
						cerr << "ERROR::Worker " << this->get_Rank() << " allocation error during adding of contour points!" << endl;
						delete [] rbuf; MPI_File_close(&ioReadFileHdl); return false;
					}
					
					gcontour2d[gid] = contourpoints;
					gcontour2d[gid]->push_back( point2d(rbuf[e].x, rbuf[e].y) );
					continue;
					//MK::mind that the algorithm requires the point data to be sorted either clock- or anti-clockwise!
					//GraGLeS 2d currently outputs clockwise into the GBContourPoints file
					//and outputs as the last point to each contour/grain a copy of the first
				}
				//problems with the grainID and the contour point occurred
				cerr << "ERROR::Worker " << this->get_Rank() << " inconsistent contour point found!" << endl;
				stillgood = false; break;
			}
			//buffer processed completely
			elementsRead = elementsRead + elementsNow;
			delete [] rbuf; rbuf = NULL;
		}
	} //proceed with reading next block

	MPI_File_close(&ioReadFileHdl);

	myprofiler.logev( "Read" + fn, (double) (MPI_Wtime() - timer) );
	cout << "...Worker " << this->get_Rank() << " " << fn.c_str() << " grain contour data were loaded successfully in " << (double) (MPI_Wtime() - timer) << " seconds" << endl;
	
	return stillgood;
}


struct aabb2d arealanaHdl::calc_polygon_aabb( std::vector<point2d>* thecontour, unsigned int gID ) {
	struct aabb2d res;
	for ( unsigned int i = 0; i < thecontour->size(); i++ ) { //was -1 ####MK::
		if ( thecontour->at(i).x <= res.xmi ) res.xmi = thecontour->at(i).x;
		if ( thecontour->at(i).y <= res.ymi ) res.ymi = thecontour->at(i).y;
		if ( thecontour->at(i).x >= res.xmx ) res.xmx = thecontour->at(i).x;
		if ( thecontour->at(i).y >= res.ymx ) res.ymx = thecontour->at(i).y;
	}
	res.gid = gID;
//cout << "Grain " << res.gid << " xmi/ymi/xmx/ymx = " << res.xmi << ";" << res.ymi << ";" << res.xmx << ";" << res.ymx << endl;

	return res;
}


double arealanaHdl::calc_polygon_area( std::vector<point2d>* thecontour ) {
	if ( thecontour == NULL ) return 0.0;
	if ( thecontour->size() == 0 ) return 0.0;

	double A = 0.0;
	for ( unsigned int i = 0; i < thecontour->size() - 1; i++ ) {
		A = A + (thecontour->at(i).x * thecontour->at(i+1).y - thecontour->at(i+1).x * thecontour->at(i).y);
	}
	return 0.5*A;
}


void arealanaHdl::init( void ) {
	//init contour buckets, as ID hashs
	double timer = MPI_Wtime();

	struct isectionmeta dummy;
	unsigned int ngr = 1 + Settings::LargestGrainID;
	for ( unsigned int gid = 0; gid < ngr; gid++ ) {
		gcontour2darea.push_back( 0.0 ); //area not yet known
		gcontour2d.push_back( NULL );
		gsegments2d.push_back( NULL );
		//set up meta data grain ID hash, 
		//MK::a compromise approach, namely if grainIDs are sparse, the hash provokes frequent far field fetches for which searching through a list would be faster
		//however, if gmetadata is populated dense, searching would always have to be executed for all to avoid one far field miss
		gmetadata.push_back( dummy );
	}

	myprofiler.logev( "InitializeMemory", (double) (MPI_Wtime() - timer) );
}


void arealanaHdl::reverse_contours( void ) {
	//changes arrangement of contour line points from C. Mießen's GraGeLeS outputs ######MKclock-clockwise we need clock-wise!
	for ( unsigned int gc = 0; gc < gcontour2d.size(); gc++ ) {
		if ( gcontour2d.at(gc) != NULL ) 
			std::reverse( gcontour2d.at(gc)->begin(), gcontour2d.at(gc)->end() );
	}
}


void arealanaHdl::remove_core_duplicates( void ) {
	//FOR GRAGLES MUST BE EXECUTED AFTER REVERSAL OF CONTOUR POINT ORDER!
	//MK::numerical issues arising in GraGLeS as detailed here::
	/*0.477081	0.00421709	54698	0.00419541 //-->that the first and the last are the same coordinates is expected for a closed contour
	0.476898	0.00434444	54698	0.00419541
	0.476871	0.00440044	54698	0.0300042
	0.476819	0.00458379	54698	0.0300042
	0.476819	0.00476714	54698	0.0124919
	0.476898	0.00481005	54698	0.0333015
	0.477023	0.00476714	54698	0.0333015
	0.477081	0.0047468	54698	0.0333015
	0.477251	0.00458379	54698	0.0333015
	0.477264	0.00454116	54698	0.00728969
	0.477304	0.00440044	54698	0.00728969
	0.477264	0.00434915	54698	0.0755061
	0.477081	0.00421709	54698	0.0755061 -->##but these numerical duplicates must be either rounding or numerics issues, so we have to delete these points
	0.477081	0.00421705	54698	0.0755061
	0.477081	0.00421709	54698	0.00419541*/

	//cout << "\t\tRemoving core duplicates (including the last point!)..." << endl;
	vector<double> inspectorx, inspectory;
	for ( unsigned int gc = 0; gc < gcontour2d.size(); gc++ ) {
		if ( gcontour2d.at(gc) != NULL ) { //inspect for duplicates
			//if ( gc == 54698 ) double stopper = 54698.0; //dummy

			inspectorx.clear();
			inspectory.clear();
			inspectorx.push_back( gcontour2d[gc]->at(0).x ); //hardcoded access is safe because for loop does not execute if contour is non existent or of size 0
			inspectory.push_back( gcontour2d[gc]->at(0).y );

			for ( unsigned int cp = 1; cp < gcontour2d[gc]->size(); cp++ ) { //O(n^2) problem....
				bool dupl = false;
				double testx = gcontour2d[gc]->at(cp).x;
				double testy = gcontour2d[gc]->at(cp).y;
				//inspector xy grow by adding non-duplicate points
				for ( unsigned int i = 0; i < inspectorx.size(); i++ ) { //inspectorx and y same length..
					if ( testx != inspectorx[i] ) continue; //if the same tests for y
					if ( testy != inspectory[i] ) continue; //if the same remains in loop
					//not continued so duplicate found, which we do not add! (thereby implicitly removing last duplicate point which is the first of the contour
					dupl = true;
				}
				if ( dupl == false ) { //most likely case
					inspectorx.push_back( testx );
					inspectory.push_back( testy );
				}
			}

			//update gcontour2d points by removing the old with potential duplicates and adding all unique ones
			gcontour2d.at(gc)->clear();
			unsigned int unique = inspectorx.size();
			if ( unique != inspectory.size() ) { cerr << "ERROR::Removing core duplicates for grainID " << gc << endl; }

			for ( unsigned int i = 0; i < inspectorx.size(); i++ ) {
				struct point2d ap;
				ap.x = inspectorx[i];
				ap.y = inspectory[i];
				gcontour2d.at(gc)->push_back( ap );
			}
		} //next grain
	}
}


void arealanaHdl::remove_lastpoints( void ) {
	//removes last point in each contour as C. Mießen's GraGeLeS exports both start and end point (for GNU closed contour plotting purposes)
	for ( unsigned int gc = 0; gc < gcontour2d.size(); gc++ ) {
		if ( gcontour2d.at(gc) != NULL )
			gcontour2d.at(gc)->pop_back();
	}
}


void arealanaHdl::remove_degenerated_grains( void ) {
	//very seldomly it may occurr that the modified contour line contains only two numerically different points, thus not forming at least one triangle,
	//we need to exclude such a grain from the analysis, as otherwise the Poly2Tri triangularization will fail and stop the entire workflow
	//a good example is 
	/*
	0.94256	0.106324	1075735	0.0417861
	0.94256	0.106325	1075735	0.0313236
	0.94256	0.106324	1075735	0.0313236
	0.94256	0.106324	1075735	0.0417861
	0.94256	0.106324	1075735	0.0417861
	*/ //we see immediately that within output accuracy all x are the same and there are only two disjoint y coordinate values,
	//so either increase output accuracy, however numerically there can always be a very small triangle, so better exclude the grain
	//it'll very likely become consumed by its neighbors anyway in the next time steps...

	//MK::for the "Quantifying Necessity vs. Sufficiency" paper this test was not included because the Contour data were double precision!

	for ( unsigned int gc = 0; gc < gcontour2d.size(); gc++ ) {
		if ( gcontour2d.at(gc) != NULL ) { 
			//prior call to remove_insufficiently_supported_contours() assured gcontour2d.at(gc)->size() >= 3
			bool collinear = true; //disprove
			double x = gcontour2d[gc]->at(0).x; //access to 0 possible because size() >=3
			double y = gcontour2d[gc]->at(0).y;
			for ( unsigned int p = 1; p < gcontour2d.at(gc)->size(); p++ ) { 
				//if at least one point not on the same line, no collinearity
				if ( fabs(gcontour2d.at(gc)->at(p).x - x) > COLLINEARITY_ACCURACY && fabs(gcontour2d.at(gc)->at(p).y - y) > COLLINEARITY_ACCURACY ) {
					collinear = false; break;
				}
				//if never breaked all points on the line assumption of collinear == true still holds, (very unlikely case)
			} //collinearity test done...
			if ( collinear == false )
				continue;
					
			//##MK::this function does not check collinearity for all against all!

			//less than 3 points, or collinear == true and therefore not yet continued, i.e. degenerated polygon!
			//erase data of this grain but not the metadata
			delete gcontour2d.at(gc); gcontour2d.at(gc) = NULL;

			for ( unsigned int g = 0; g < contours_bary_quick.size(); g++ ) {
				if ( contours_bary_quick[g].gid != gc ) 
					continue;
				//gc found!
				contours_bary_quick.erase( contours_bary_quick.begin() + g ); 
				break;
			}
			for ( unsigned int g = 0; g < contours_aabb_quick.size(); g++ ) {
				if ( contours_aabb_quick[g].gid != gc )
					continue;
				contours_aabb_quick.erase( contours_aabb_quick.begin() + g ); 
				break;
			}
			//if ( collinear == true )
			cerr << "WARNING::Worker " << this->get_Rank() << " grain " << gc << " was removed because of collinearity, same point, or less than 3 points" << endl;
		} //next gc
	}
}


void arealanaHdl::calc_contourarea( void )
{
	//compute prior to any elimination of grains on sufficiently supported contours (i.e. those of size > 3)
	for ( unsigned int gc = 0; gc < gcontour2d.size(); gc++ ) {
		if ( gcontour2d.at(gc) != NULL ) {
			double A = 0.0;
			unsigned int np = gcontour2d.at(gc)->size() - 1;
			
			for ( unsigned int p = 0; p < np; p++) { //Gaussian Trapez Formula
				A = A + (gcontour2d.at(gc)->at(p).y + gcontour2d.at(gc)->at(p+1).y) 
					* (gcontour2d.at(gc)->at(p).x - gcontour2d.at(gc)->at(p+1).x);
			}
			gcontour2darea.at(gc) = fabs(0.5*A);
		}
	}
	
	/*double current_volume = 0;
	double h = m_owningGrain->get_h();

	for (unsigned int i = 0; i < m_grainBoundary.size() - 1; i++) {
		// Gaussian Trapez Formula:
		current_volume += (m_grainBoundary[i].y + m_grainBoundary[i + 1].y)
				* (m_grainBoundary[i].x - m_grainBoundary[i + 1].x);
	}*/
}


void arealanaHdl::remove_insufficiently_supported_contours( void )
{
	for ( unsigned int gc = 0; gc < gcontour2d.size(); gc++ ) {
		if ( gcontour2d.at(gc) != NULL ) {
			if ( gcontour2d[gc]->size() > 3 ) 
				continue;

			//else, insufficient number of supporting points for this contour
			delete gcontour2d.at(gc);
			gcontour2d.at(gc) = NULL; //MK::flag contour implicitly as do not analyze because none existent!
		}
	}
}


void arealanaHdl::modify_rawdata_gragles( unsigned int fd ) {
	//implement model output specific preprocessing operations here...
	cout << "...Worker " << this->get_Rank() << " checks contour data for sufficient number of supporting points, collinearity ... " << endl;
	
	double timer1 = MPI_Wtime();
	calc_contourarea();
	double timer2 = MPI_Wtime();
	cout << "...Worker " << this->get_Rank() << " calculated contour-based area of the grain in " << timer2 - timer1 << " seconds" << endl; 
	myprofiler.logev(  "ContourBasedArea" + std::to_string(fd), (double) (timer2 - timer1) ); 

	timer1 = MPI_Wtime();
	remove_insufficiently_supported_contours();
	timer2 = MPI_Wtime();
	cout << "...Worker " << this->get_Rank() << " sufficiency of contour discretization checked in " << timer2 - timer1 << " seconds" << endl; 
	myprofiler.logev(  "CheckSuppPointSufficiency" + std::to_string(fd), (double) (timer2 - timer1) ); 
	
	timer1 = timer2;
	calc_barycenter();
	timer2 = MPI_Wtime();
	cout << "...Worker " << this->get_Rank() << " calculated barycenter2d in " << timer2 - timer1 << " seconds" << endl;
	myprofiler.logev( "CalcBarycentre" + std::to_string(fd), (double) (timer2 - timer1) );
	
	timer1 = timer2;
	calc_aabb();
	timer2 = MPI_Wtime(); 
	cout << "...Worker " << this->get_Rank() << " calculated aabb2d in " << timer2 - timer1 << " seconds" << endl;
	myprofiler.logev( "CalcAABB" + std::to_string(fd), (double) (timer2 - timer1) );

	timer1 = timer2;
	reverse_contours(); //##MK::not necessary GraGeLeS outputs gnu files clockwise, as PolyTri requires?
	timer2 = MPI_Wtime();
	//cout << "...Worker " << this->get_Rank() << " reversed contours in " << timer2 - timer1 << " seconds" << endl;
	myprofiler.logev( "ReverseContours" + std::to_string(fd), (double) (timer2 - timer1) );
	
	//timer1 = timer2;
	//remove_core_duplicates();  //required when importing contour points with limited precision, like for classical gnu file, or floats translated into ASCII
	//timer2 = MPI_Wtime();
	//cout << "...Worker " << this->get_Rank() << " removed core duplicates including the lastpoint " << timer2 - timer1 << " seconds" << endl;
	//myprofiler.logev( "RemoveCoreDupl" + std::to_string(fd), (double) (timer2 - timer1) );
	
	timer1 = timer2;
	remove_lastpoints(); //MK::DO NOT EXECUTE GBContourPoints_<fid>.bin outputs the last point the duplicate of the first point
	timer2 = MPI_Wtime();
	cout << "...Worker " << this->get_Rank() << " removed last points in " << timer2 - timer1 << " seconds" << endl;
	myprofiler.logev( "RemoveLastContourPoint" + std::to_string(fd), (double) (timer2 - timer1) );

	//MK::for the "Quantifying Necessity vs. Sufficiency" paper this test was not included because the Contour data were double precision!
	//timer1 = timer2;
	//remove_degenerated_grains();
	//timer2 = MPI_Wtime();
	//cout << "...Worker " << this->get_Rank() << " removed degenerated grains " << timer2 - timer1 << " seconds" << endl;
	//myprofiler.logev( "RemoveDegenerated" + std::to_string(fd), (double) (timer2 - timer1) );

	cout << "...Worker " << this->get_Rank() << " adjusted GraGeLeS2D-specific output format successfully!" << endl;
}


std::vector<unsigned int>* arealanaHdl::triangularize_contours( unsigned int fd) {
	//##MK::the triangularization is at the moment not straightforward to parallelize as the PolyTri library
	//utilizes currently static constructs for the SplayTree which are not guaranteed to be threadsafe...
	//triangularize and report which gIDs had problems, memorized them in the blacklist to enable eliminating them in the next step
	double timer = MPI_Wtime();

	//metrics
	double DomainAreaInTotal = 0.0; //to account quantitatively the significance of eliminating contours
	double DomainAreaUnpartitioned = 0.0;
	unsigned int npolygons = 0;
	unsigned long ntriangles = 0;

	std::vector<unsigned int>* blacklist = NULL;
	try { blacklist = new std::vector<unsigned int>; }
	catch (std::bad_alloc &exc) { 
		cerr << "ERROR::Worker " << this->get_Rank() << " unable to allocate memory for the black list!" << endl; 
		return NULL;
	}
	//this->test_polytri();

	cout << "...Worker " << this->get_Rank() << " triangularizing via the PolyTri library..." << endl;
	for ( unsigned int gc = 0; gc < gcontour2d.size(); gc++ ) {
		//MK::clock-wise definition of points! no last point duplicate!
		if ( gcontour2d.at(gc) != NULL ) {			
			double contourarea = gcontour2darea.at(gc); //##MK::far field fetch
			DomainAreaInTotal += contourarea;

			Polygon* contour = new Polygon( gcontour2d.at(gc) ); //parse constructor, utilize the PolyTri library code to do the core stuff
			
			//MK::in case of being interest in a single grain
			/*if ( gc == 716218 ) {
				for ( unsigned int p = 0; p < gcontour2d.at(gc)->size(); p++ ) {
					cout << setprecision(32) << gcontour2d.at(gc)->at(p).x << ";" << gcontour2d.at(gc)->at(p).y << endl;  
				}
			}*/
			
			if ( contour->is_healthy() == false ) {
				cerr << "ERROR::Worker " << this->get_Rank() << " Polygon " << gc << " is not healthy!" << endl; //MK::catastrophic
				DomainAreaUnpartitioned += contourarea;

				blacklist->push_back( gc );
				delete contour; contour = NULL; continue;
				//MK::one may alos jump out because this usually indicates a significant problem with the contour line geometry return NULL;
			}

			if ( contour->triangulation() == false ) { //accept some of such guys, if not too many... MK::however they reduce the total grian population
				//and eventually cause the analysis of some targets to become biased because a grain was punched out from the inspection area...
				cerr << "ERROR::Worker " << this->get_Rank() << " unable to triangularize contour " << gc << endl;
				DomainAreaUnpartitioned += contourarea;

				blacklist->push_back( gc );
				delete contour; contour = NULL; continue;
			}

			//triangularization was successful
			gsegments2d.at(gc) = contour->giveTriangles();

			if ( gsegments2d.at(gc) != NULL ) { //passing of triangle collection successful

#ifdef PERSISTENCE_VERBOSE
				cout << "\t\tGrain " << gc << " has " << gsegments2d[gc]->size() << " triangles." << endl;
#endif
				//delete gsegments2d.at(gc); //##MK::DUMMY
				npolygons++;
				ntriangles += gsegments2d[gc]->size();
			} 
			else { //passing was unsuccessful, catastrophic
				cerr << "ERROR::Worker " << this->get_Rank() << " grain " << gc << " has no triangles." << endl;
				delete contour; contour = NULL;
				delete blacklist; blacklist = NULL;
				fdmeta.TotalAreaEnclosedByContour = 0.0;
				fdmeta.TotalAreaTriangularizable = 0.0; //all other struct elements set already by constructor
				fdmeta.TotalAreaNonTriangularizable = 1.0; //GraGLeS working on [0,1]^2 domain
				return NULL;
			}

			delete contour; contour = NULL;
		} //done with processing existent contour for grain gc
	} //next grain <=> polygon

	myprofiler.logev( "Triangularize" + std::to_string(fd), ((double) MPI_Wtime() - timer) );
	cout << "...Worker " << this->get_Rank() << " area of all contours adds up to " << DomainAreaInTotal << " for which " << DomainAreaUnpartitioned << " were not triangularizable" << endl;
	cout << "...Worker " << this->get_Rank() << " triangularized all contours in " << (MPI_Wtime() - timer) << " seconds" << endl;
	cout << "...Worker " << this->get_Rank() << " a total of " << blacklist->size() << " grains are suspious and require elimination!" << endl;
	cout << "...Worker " << this->get_Rank() << " a total of " << npolygons << " grains triangularized into a total of " << ntriangles << " triangles" << endl;

	fdmeta.TotalAreaEnclosedByContour = DomainAreaInTotal;
	fdmeta.TotalAreaTriangularizable = DomainAreaInTotal - DomainAreaUnpartitioned;
	fdmeta.TotalAreaNonTriangularizable = DomainAreaUnpartitioned;
	fdmeta.TotalTriangleCnts = ntriangles;
	fdmeta.TotalPolygonsTriangularized = npolygons;
	fdmeta.TotalPolygonsNonTriangularizable = blacklist->size();

	return blacklist;
}


void arealanaHdl::calc_barycenter( void ) {
	for ( unsigned int gc = 0; gc < gcontour2d.size(); gc++ ) {
		if ( gcontour2d.at(gc) != NULL ) {
			unsigned int gID = gc; //MK::utilize that gcontour2d.size == 1 + Settings::LargestGrainID

			//reference position is always initial location of the sub-grain x0, y0
			//MK::we have to promote the unsigned int to doubles
			double x0 = gmetadata.at(gID).x0;
			double y0 = gmetadata.at(gID).y0;
			x0 /= Settings::InitialDomainEdgeLength; //living on square domain
			y0 /= Settings::InitialDomainEdgeLength;

			contours_bary_quick.push_back( barycenter2d(x0, y0, gID)  );
		}
	}
}


void arealanaHdl::calc_aabb( void ) {
	//init axis-aligned bounding box in 2D
	for ( unsigned int gc = 0; gc < gcontour2d.size(); gc++ ) { 
		if ( gcontour2d.at(gc) != NULL ) {
			contours_aabb_quick.push_back( calc_polygon_aabb( gcontour2d.at(gc) , gc ) );
		}
	}
}


struct overlap2d arealanaHdl::calc_overlap_polygon_circle( unsigned int gidx, double xcir, double ycir, double r, bool render )
{
	unsigned int ntri = gsegments2d.at(gidx)->size();
	double res = 0.0;
	double leakage = 0.0;

	for ( unsigned int tri = 0; tri < ntri; tri++ ) {
		double x[3] = { gsegments2d.at(gidx)->at(tri).x1, gsegments2d.at(gidx)->at(tri).x2, gsegments2d.at(gidx)->at(tri).x3 };
		double y[3] = { gsegments2d.at(gidx)->at(tri).y1, gsegments2d.at(gidx)->at(tri).y2, gsegments2d.at(gidx)->at(tri).y3 };

		if (render == true) {
			cout << gidx << ";" << setprecision(32) << x[0] << ";" << y[0] << ";" << x[1] << ";" << y[1] << ";" << x[2] << ";" << y[2] << endl;
		}

		//compute individual triangle circle overlap analytically, either works, then add to res, otherwise numerical fails add to leakage
		double A_analytical = intersection_area_triangle_circle( x[0], y[0], x[1], y[1], x[2], y[2], xcir, ycir, r );

		//most likely computation is correct
		if ( A_analytical >= 0.0 ) {
			res = res + A_analytical;
			continue;
		}

		//not continued, so numerical issues during overlap computation
		double A_triangle = triangle_area( x[1]-x[0], y[1]-y[0], x[2]-x[0], y[2]-y[0] );
		//heuristically account for reference area leakages
		leakage += A_triangle;
		cerr << "WARNING::Leakage gidx/tri/Atri/Aana/xy012xcycr;" << gidx << ";" << tri << ";" << setprecision(18) << A_triangle << ";" << setprecision(18) << A_analytical << ";" << setprecision(18) << x[0] << ";" << setprecision(18) << y[0] << ";" << setprecision(18) << x[1] << ";" << setprecision(18) << y[1] << ";" << setprecision(18) << x[2] << ";" << setprecision(18) << y[2] << ";" << setprecision(18) << xcir << ";" <<  setprecision(18) << ycir << ";" << setprecision(18) << r << endl;
		//we consider this triangle as leaking but it does not completely invalidate the overlap computation for the entire polygon!

		//do not add to total because overlap considered as leakage
		//res = res + 0.0;

#ifdef TRIANGLE_AREA_VIA_RASTERIZING
		//##MK::check correct size of Settings::InitialDomainEdgeLength
		double A_triangle = triangle_area( x[1]-x[0], y[1]-y[0], x[2]-x[0], y[2]-y[0] );
		double A_discrete = intersection_area_triangle_circle_mc( x[0], y[0], x[1], y[1], x[2], y[2], xcir, ycir, r, Settings::InitialDomainEdgeLength );
		//consistency
		if ( A_analytical > A_triangle || fabs(1.0 - (A_analytical/A_discrete) ) > 0.02 ) {
			//cout << "Wrong analytical solution;Atri/Aanalytic/Adiscr;" <<
			cout << setprecision(18) << A_triangle << ";" <<  setprecision(18) << A_analytical << ";" <<  setprecision(18) << A_discrete;
			cout << ";" <<  setprecision(18) << x[0] << ";" <<  setprecision(18) << y[0] << ";" <<  setprecision(18) << x[1] << ";" <<  setprecision(18) << y[1] << ";" <<  setprecision(18) << x[2] << ";" <<  setprecision(18) << y[2];
			cout << ";" <<  setprecision(18) << xcir << ";" <<  setprecision(18) << ycir << ";" <<  setprecision(18) << r << endl;
		}
		//else { cout << "Atri/Aanalytic/Adiscr;" << A_triangle << ";" <<  setprecision(18) << A_analytical << ";" <<  setprecision(18) << A_discrete << endl; }
#endif
	}

	return overlap2d( res, leakage );
}


void arealanaHdl::eliminate_non_triangularized_grains( unsigned int fd, std::vector<unsigned int>* gid_to_kick )
{
	//eliminate grainIDs from the analysis for which the triangularization was unsuccessful...
	double timer = MPI_Wtime();
	
	//I am safe to access and have sth to eliminate
	cout << "...Worker " << this->get_Rank() << " eliminating faulty polygons..." << endl;
	for ( unsigned int gk = 0; gk < gid_to_kick->size(); gk++ ) {
		unsigned int killgid = gid_to_kick->at(gk);
cout << "...Worker " << this->get_Rank() << " killing contour for grainID " << killgid << endl;

		if ( gcontour2d.at(killgid) != NULL ) {
			delete gcontour2d.at(killgid); //contour exists but was unsuccessful to triangularize
			gcontour2d.at(killgid) = NULL;
		}

		if ( gsegments2d.at(killgid) != NULL ) { //##MK::access was done via at(gk) very likely faulty!
			delete gsegments2d.at(killgid); //segments as well
			gsegments2d.at(killgid) = NULL;
		}

		for ( unsigned int g = 0; g < contours_bary_quick.size(); g++ ) { //##MK::hash-based algorithms would avoid this O(N^2) part...
			if ( contours_bary_quick[g].gid != killgid ) 
				continue;
			//gk found!
			contours_bary_quick.erase( contours_bary_quick.begin() + g ); 
			break;
		}
		for ( unsigned int g = 0; g < contours_aabb_quick.size(); g++ ) {
			if ( contours_aabb_quick[g].gid != killgid )
				continue;
			contours_aabb_quick.erase( contours_aabb_quick.begin() + g ); 
			break;
		}
	} //next grainID to eliminate
	
	myprofiler.logev( "ElimNonTriangularizable" + std::to_string(fd), (double) (MPI_Wtime() - timer) );
	cout << "...Worker " << this->get_Rank() << " eliminated a total of " << gid_to_kick->size() << " grains that were non-triangularizable!" << endl;
}


//MK::MUST NOT DEFINE PLOTGID if in doubt what the source code of arealanaHdl::envcharacterize_omp_naive does exactly!
//the define is meant only to visualize for an isolated grain its triangularization and of all the inspected neighbors however we are in an parallel region!
//#define PLOTGID

void arealanaHdl::envcharacterize_omp_naive( unsigned int fd ) {
	//MK::we work on the unit surface [0,1] \in {\mathbb{R}}^2
	//double halfGbSQR = Settings::DislocEnPerM; //0.5*Gb^2
	//OMP-parallel version without considering explicit thread local data placement..., therefore naive OMP parallelism, lets see how far we get with this...
	double timer = MPI_Wtime();
	cout << "...Worker " << this->get_Rank() << " starting LongRangeAnalysis Overlap Computations on " << fd << "..." << endl;

	double LeakageArea = 0.0;
	unsigned long LeakageCnts = 0;
	unsigned long ntriangles_inspected_total = 0;
	
	//master initiates process-local ChiTable
	unsigned int rc = 0;
	for ( double r = Settings::LongRangeRadiusMin; r <= Settings::LongRangeRadiusMax; r = r + Settings::LongRangeRadiusIncr ) { //implicit order
		std::vector<chi>* achi = NULL;
		achi = new vector<chi>;
		ChiTable.push_back( achi );
		rc++;
	}

	#pragma omp parallel
	{
		//shared by default:: ntriangles_inspected_total, timer, and ChiTable
		
		double localleakage_area = 0.0;
		unsigned long localleakage_cnts = 0;
		//threads allocate thread local results memory
		std::vector<std::vector<chi>*> LocalChiTable; //container of pointer only
		
		//all what follows is in a parallel region and thus private to the thread
		unsigned int ng = contours_bary_quick.size(); //MK::not an ID hash as its size is not necessarily 1+Settings::LargestGrainID!
		vector<unsigned int>* candidates = NULL;
		try { candidates = new vector<unsigned int>; }
		catch (std::bad_alloc &exc) { cerr << "ERROR::Allocation error for thread local chi memory in envcharacterize_longrange_core!" << endl; }

		//prepare LocalChiTable
		unsigned int rclass = 0;
		for ( double r = Settings::LongRangeRadiusMin; r <= Settings::LongRangeRadiusMax; r = r + Settings::LongRangeRadiusIncr ) { //implicit order
			std::vector<chi>* achi = NULL;
			try { achi = new vector<chi>; }
			catch (std::bad_alloc &exc) { cerr << "ERROR::Allocation error for thread local chi memory in envcharacterize_longrange_core!" << endl; }
			LocalChiTable.push_back( achi );
			rclass++;
		}

		unsigned long ntriangles_inspected_local = 0;
		unsigned int nthreads = omp_get_num_threads();
		unsigned int mythread = omp_get_thread_num();
		/*#pragma omp master
		{
			cout << "...Worker " << this->get_Rank() << " performing analysis with " << nthreads << " OpenMP threads!" << endl;
		}*/

		//local copies of result arrays (in threadlocal memory) #pragma omp barrier
		
		unsigned int plotgid = 772878; //which grain ID to plot

		for ( unsigned int g = 0; g < ng; g++ ) { //loop over all still existing grains, BE CAREFUL g reads not as a grainID but an index!

#ifdef PLOTGID
			//cout << g << "\t\t" << contours_bary_quick[g].gid << "\t\t" << plotgid << "\t\t" << WorkPartitioning(g, nthreads, mythread) << endl;
			if (contours_bary_quick[g].gid == plotgid) {
#endif
				if (WorkPartitioning(g, nthreads, mythread) == true) { //if I take care of this grain

					double gtimer = omp_get_wtime();

					//quick screening to include only grains whose circular investigation area falls completely inside the domain
					double rmax = Settings::LongRangeRadiusMax;
					unsigned int gid = contours_bary_quick[g].gid;
					double xinv = contours_bary_quick[g].x;
					double yinv = contours_bary_quick[g].y;
					double xinv_mr = xinv - rmax; //living on [0,1]^2 \in {\mathcal{R}}^2
					double xinv_pr = xinv + rmax;
					double yinv_mr = yinv - rmax;
					double yinv_pr = yinv + rmax;

#ifdef PLOTGID
					cout << setprecision(32) << xinv_mr << ";" << xinv_pr << "\t\t\t" << yinv_mr << ";" << yinv_pr << endl;
#endif

					if (xinv_mr < 0.0) continue; //MK::causes frequent branch misprediction but at an early stage and at most ng times!
					if (xinv_pr > 1.0) continue;
					if (yinv_mr < 0.0) continue;
					if (yinv_pr > 1.0) continue;

					//not continued?, so screening passed, so work on this grain, now multiple circle computations would always require to compute the neighbors for each radius
					//MK::therefore: build a hashtable of all possible candidate neighbors for the largest circle instead of naive testing all against all r times...

					//MK::DEBUG to inspect algorithm for some exemplary guys...
					//if ( gid == 152205 || gid == 401356 || gid == 513695 ) cout << "Lost guys still in game!" << endl;

					//eliminate all neighboring candidate grains whose AABB do not overlap at all with the circle
					candidates->clear();
					for (unsigned int candg = 0; candg < g; candg++) { //##MK::now a compromise between lean tree-based find of close objects and naive all to all
						if (contours_aabb_quick[candg].xmx < xinv_mr) continue; //too much left
						if (contours_aabb_quick[candg].xmi > xinv_pr) continue; //too much right
						if (contours_aabb_quick[candg].ymx < yinv_mr) continue; //too much low
						if (contours_aabb_quick[candg].ymi > yinv_pr) continue; //too much up 
						//not continued-->overlapping of AABB from candg with target, so add candidate id to the list of candidates
						candidates->push_back(candg); //MK::again these are indices on contours_bary and aabb NOT grainIDs!
					}
					//MK::twice such a for loop two avoid checking the only once occurring case of a selfintersection...
					for (unsigned int candg = g + 1; candg < ng; candg++) { //two loops to avoid self-reference
						if (contours_aabb_quick[candg].xmx < xinv_mr) continue;
						if (contours_aabb_quick[candg].xmi > xinv_pr) continue;
						if (contours_aabb_quick[candg].ymx < yinv_mr) continue;
						if (contours_aabb_quick[candg].ymi > yinv_pr) continue;
						candidates->push_back(candg);
					}

					//core analysis, check for all candidates whether any of there triangles partially intrude the inspection circle
					unsigned long ntriangles_inspected = 0;
					if (candidates->size() > 0) { //compute overlap
						double qme[4] = { gmetadata[gid].q0, gmetadata[gid].q1, gmetadata[gid].q2, gmetadata[gid].q3 };
						double seeme = gmetadata[gid].see;

						rclass = 0;
						for (double r = Settings::LongRangeRadiusMin; r <= Settings::LongRangeRadiusMax; r = r + Settings::LongRangeRadiusIncr) { //loop over radii
							//identify all triangles of grains inside the circular region of radius r about (xinv, yinv), characterize environment about grain g
							struct chi mp;
							mp.x = xinv;
							mp.y = yinv;
							mp.r = r;
							mp.gid = gid;
							ntriangles_inspected = ntriangles_inspected + gsegments2d.at(gid)->size();

							struct overlap2d ares;
#ifdef PLOTGID
							cout << "plotgid;" << plotgid << endl;
							cout << "r/x/y;" << setprecision(32) << mp.r << "; " << mp.x << "; " << mp.y << endl;

							ares = calc_overlap_polygon_circle(mp.gid, mp.x, mp.y, mp.r, true);

							cout << endl << endl; //inspection triangles follow...
#else
							ares = calc_overlap_polygon_circle(mp.gid, mp.x, mp.y, mp.r, false);
#endif

							mp.Aself = ares.area;
							if (ares.leakage > 0.0) {
								localleakage_area += ares.leakage; localleakage_cnts++;
							}

							double chi_drv = 0.0;
							double chi_spd = 0.0;
							double chi_dis = 0.0;
							double FractionalAreaInCircle = 0.0;
							double FractionalAreaInCircleTotal = 0.0;

							for (unsigned int c = 0; c < candidates->size(); c++) { //MK::benefit::candidates->size() is now much shorter then ng !
								unsigned int thatguy = candidates->at(c);
								if (contours_aabb_quick[thatguy].xmx < (mp.x - r)) continue; //##MK::cache misses
								if (contours_aabb_quick[thatguy].xmi > (mp.x + r)) continue;
								if (contours_aabb_quick[thatguy].ymx < (mp.y - r)) continue;
								if (contours_aabb_quick[thatguy].ymi > (mp.y + r)) continue;

								//not continued therefore AABB of candidate candg has overlap with circle so there may be triangles overlapping the circle
								unsigned int gidx = contours_bary_quick[thatguy].gid; //far field fetch ... MUST NOT BE candg as ng is usually much smaller than gmetadata->size()

								ntriangles_inspected = ntriangles_inspected + gsegments2d.at(gidx)->size();

								struct overlap2d ar1;
#ifdef PLOTGID
								ar1 = calc_overlap_polygon_circle(gidx, mp.x, mp.y, mp.r, true);
#else
								ar1 = calc_overlap_polygon_circle(gidx, mp.x, mp.y, mp.r, false);
#endif

								FractionalAreaInCircle = ar1.area;
								if (ar1.leakage > 0.0) {
									localleakage_area += ar1.leakage; localleakage_cnts++;
								}

								double qcand[4] = { gmetadata[gidx].q0, gmetadata[gidx].q1, gmetadata[gidx].q2, gmetadata[gidx].q3 };
								double seecand = gmetadata[gidx].see; //##MK::unfortunate far field fetch on hash-vector...!
								double mob = dis2mob(misorientationCubic(qme, qcand));
				//cout << "seecand/seeme = " << seecand << ";" << seeme << endl;
								chi_drv = chi_drv + (Settings::DislocEnPerM*(seecand - seeme) * FractionalAreaInCircle);
								chi_spd = chi_spd + (Settings::HAGBMobility * mob *  Settings::DislocEnPerM*(seecand - seeme) * FractionalAreaInCircle);
								chi_dis = chi_dis + (Settings::HAGBMobility * mob * FractionalAreaInCircle);
								FractionalAreaInCircleTotal = FractionalAreaInCircleTotal + FractionalAreaInCircle;
							} //tested all candidates for radius r

							mp.Anbors = FractionalAreaInCircleTotal;
							mp.Adiff = PI*SQR(mp.r) - mp.Aself - mp.Anbors;
							if (FractionalAreaInCircleTotal > 0.0) {
								chi_drv = chi_drv / FractionalAreaInCircleTotal;
								chi_spd = chi_spd / FractionalAreaInCircleTotal;
								chi_dis = chi_dis / FractionalAreaInCircleTotal;
								mp.Pdrv = chi_drv;
								mp.Pspd = chi_spd;
								mp.Pdis = chi_dis;
							}

							LocalChiTable.at(rclass)->push_back(mp);
							rclass++;
							//cout << "---->" << mp.gid << ";" << mp.x << ";" << mp.y << ";" << mp.r << ";" << mp.Aself << ";" << mp.Anbors << ";"<< mp.Adiff << ";" << mp.P << "\t" << FractionalAreaInCircleTotal << "\t took " << (double) (MPI_Wtime() - gtimer ) << " seconds " << endl;
							//##MK::reutilize already computed areas, but limited because circle enlarges...
						} //next radius
					} //done processing existent candidates
					else { //MK:: when candidates is empty that means all neighbors do not overlap with the AABB of the circle, hence the grain has grown at least to the size of the circle LongRangeRadiusMax size!
						//still though we have to keep track of this as otherwise we would loose an ID for the time series analysis
						rclass = 0;
						for (double r = Settings::LongRangeRadiusMin; r <= Settings::LongRangeRadiusMax; r = r + Settings::LongRangeRadiusIncr) { //loop over radii
							struct chi mp;
							mp.gid = gid;
							mp.x = xinv;
							mp.y = yinv;
							mp.r = r;

							ntriangles_inspected = ntriangles_inspected + gsegments2d.at(gid)->size();

							struct overlap2d ar2;
							ar2 = calc_overlap_polygon_circle(mp.gid, mp.x, mp.y, mp.r, false);

							mp.Aself = ar2.area;
							if (ar2.leakage > 0.0) {
								localleakage_area += ar2.leakage; localleakage_cnts++;
							}

							mp.Anbors = 0.0;
							mp.Adiff = PI*SQR(mp.r) - mp.Aself - mp.Anbors;
							mp.Pdrv = 0.0; //overwriting error value std::numeric_limits<double>::lowest()
							mp.Pspd = 0.0;
							mp.Pdis = 0.0;

							LocalChiTable.at(rclass)->push_back(mp);
							rclass++;
						} //next radius
					}
					ntriangles_inspected_local = ntriangles_inspected_local + ntriangles_inspected;
#ifdef PERSISTENCE_VERBOSE
#pragma omp critical
					{
						cout << "---->" << gid << "\t\t" << ntriangles_inspected << "\t\t" << setprecision(8) << (double)(omp_get_wtime() - gtimer) << " seconds " << endl;
					}
#endif
				} //executing thread done with grain g
#ifdef PLOTGID
			}
#endif
		} //end of processing grains

		#pragma omp critical
		{
			//threads reduce local ntriangle counts
			ntriangles_inspected_total = ntriangles_inspected_total + ntriangles_inspected_local;
			LeakageArea = LeakageArea + localleakage_area;
			LeakageCnts = LeakageCnts + localleakage_cnts;

			//threads copy local results
			unsigned int rcc = 0;
			for ( double r = Settings::LongRangeRadiusMin; r <= Settings::LongRangeRadiusMax; r = r + Settings::LongRangeRadiusIncr ) { //implicit order
for (unsigned int localg = 0; localg < LocalChiTable.at(rcc)->size(); localg++) {
	ChiTable.at(rcc)->push_back(LocalChiTable.at(rcc)->at(localg));
}

//threads delete local data
delete LocalChiTable.at(rcc); LocalChiTable.at(rcc) = NULL;
rcc++;
			}
		} //critical syncing of results done

	} //end of parallel region

	fdmeta.OverlapLeakageTotal = LeakageArea;
	fdmeta.OverlapLeakageCnts = LeakageCnts;
	fdmeta.TotalTriangleVisitsOverlap = ntriangles_inspected_total;

	myprofiler.logev("EnvCharacterizeOMP" + std::to_string(fd), (double)(MPI_Wtime() - timer));
	cout << "...Worker " << this->get_Rank() << " LongRangePersistency analysis with a total of " << fdmeta.TotalTriangleVisitsOverlap << " triangle inspections performed in " << setprecision(8) << (double)(MPI_Wtime() - timer) << " seconds" << endl;
	cout << "...Worker " << this->get_Rank() << " TotalLeakage " << setprecision(18) << fdmeta.OverlapLeakageTotal << " area in " << fdmeta.OverlapLeakageCnts << " incidences" << endl;
	//MK::the TotalLeakage is the sum of triangles thus NOT the area of the not accounted for intersection area!
}


void arealanaHdl::report(unsigned int fd) {
	double timer = MPI_Wtime();

	cout << "...Worker " << this->get_Rank() << " writing results" << endl;
	string prefix = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".LR2D.FID." + std::to_string(fd) + ".RCLASS.";
	string postfix = ".csv";
	string logfname;
	ofstream logfile;

	unsigned int rclass = 0;
	for (double r = Settings::LongRangeRadiusMin; r <= Settings::LongRangeRadiusMax; r = r + Settings::LongRangeRadiusIncr) { //loop over radii classes
		logfname = prefix + std::to_string(rclass) + postfix;
		logfile.open(logfname.c_str());
		logfile << "SnapshotID;" << fd << "\n";
		logfile << "r;" << r << "\n";
		logfile << "InspectionArea;" << setprecision(32) << SQR(METER2MICRON(r*Settings::PhysicalDomainSize))*PI << ";micron\n";
		logfile << "TotalAreaEnclosedByContours;" << setprecision(32) << fdmeta.TotalAreaEnclosedByContour << "\n";
		logfile << "TotalAreaEnclosedButNotTriangularizable;" << setprecision(32) << fdmeta.TotalAreaNonTriangularizable << "\n";
		logfile << "TotalTriangleCnts;" << setprecision(32) << fdmeta.TotalTriangleCnts << "\n";
		logfile << "TotalPolygonsTriangularized;" << setprecision(32) << fdmeta.TotalPolygonsTriangularized << "\n";
		logfile << "TotalPolygonsNonTriangularizable;" << setprecision(32) << fdmeta.TotalAreaNonTriangularizable << "\n";
		logfile << "OverlapLeakageAreaTotal;" << setprecision(32) << fdmeta.OverlapLeakageTotal << "\n";
		logfile << "OverlapLeakageNumberIncidences;" << setprecision(32) << fdmeta.OverlapLeakageCnts << "\n";
		logfile << "OverlapComputationTotalTrianglesVisits;" << setprecision(32) << fdmeta.TotalTriangleVisitsOverlap << "\n";
		logfile << "\n";
		logfile << "Grain;x;y;r;Aself;Anbors;Adiff;AreaAvEffDrivForceToGrowth;AreaAvEffMigrationSpeedToGrowth;AreaAvEffMobilityEnvironment\n";
		logfile << "Grain;1;1;1;1;1;1;J/m^3;m/s;m^4/Js\n";
		logfile << "RClass=" << METER2MICRON(r*Settings::PhysicalDomainSize) << " micron;x;y;Aself;Anbors;Adiff;AreaAvEffDrivForceToGrowth;AreaAvEffMigrationSpeedToGrowth;AreaAvEffMobilityEnvironment\n";

		std::vector<chi>* ths = ChiTable.at(rclass);
		for (unsigned long i = 0; i < ths->size(); i++) {
			logfile << ths->at(i).gid << ";" << setprecision(16) << ths->at(i).x << ";" << ths->at(i).y << ";" << ths->at(i).r << ";";
			logfile << ths->at(i).Aself << ";" << ths->at(i).Anbors << ";" << ths->at(i).Adiff << ";";
			logfile << ths->at(i).Pdrv << ";" << ths->at(i).Pspd << ";" << ths->at(i).Pdis << "\n";
		}
		logfile << endl;

		logfile.flush();
		logfile.close();
		rclass++;
	}

	myprofiler.logev("Reporting" + std::to_string(fd), (double)(MPI_Wtime() - timer));
	cout << "...Worker " << this->get_Rank() << " wrote results for fid " << fd << " in " << (MPI_Wtime() - timer) << " seconds" << endl;
}


void arealanaHdl::metareport(void) {
	string logfname;
	ofstream logfile;
	logfname = "LRAnalysis2D.SimID." + to_string(Settings::SimID) + ".DummyMetaReport.csv";

	logfile.open(logfname.c_str(), std::ofstream::out | std::ofstream::trunc);
	if (logfile.is_open() == true) {
		logfile << "Grain;x0;y0;see\n";
		logfile << "Grain;1;1;Pa\n";
		logfile << "Grain;x0;y0;see\n";
		for (unsigned int g = 0; g < contours_bary_quick.size(); g++) { //not an ID hash
			unsigned int gd = contours_bary_quick[g].gid;
			logfile << gmetadata[gd].gid << ";" << setprecision(16) << gmetadata[gd].x0 << ";" << setprecision(16) << gmetadata[gd].y0 << ";" << setprecision(16) << gmetadata[gd].see << "\n";
		}
		logfile.flush();
		logfile.close();
	}
	else {
		cerr << "ERROR::Worker::ArealAnaHdl " << this->get_Rank() << " unable to write metareport!" << endl;
	}
}


void arealanaHdl::spit_profiling(void)
{
	//first overview of profiling for sinlge process execution...
	string logfname;
	logfname = "TopologyTracer.LR2D.SimID." + std::to_string(Settings::SimID) + ".Rank." + std::to_string(this->get_Rank()) + ".Profiling.csv";
	ofstream logfile;
	logfile.open(logfname.c_str(), std::ofstream::out | std::ofstream::trunc);
	if (logfile.is_open() == true) {
		//header
		logfile << "What;WallClock\n";
		logfile << ";s\n";
		logfile << ";MPI_Wtime\n";

		for (unsigned int l = 0; l < myprofiler.get_nentries(); l++) {
			logfile << myprofiler.titles[l] << ";" << setprecision(8) << myprofiler.times[l] << "\n";
		}
		logfile.flush();
		logfile.close();
	}
	else {
		cerr << "ERROR::Worker::ArealAnaHdl " << this->get_Rank() << " unable to write profiling data!" << endl;
	}
}


void arealanaHdl::reset_contour_memory( unsigned int fd ) 
{
	double timer = MPI_Wtime();

	//delete all data in temporary data structures to allow the rank do load another set without requiring meta data re-reading
	for ( unsigned int gc = 0; gc < gcontour2d.size(); gc++  ) {
		if ( gcontour2d.at(gc) != NULL ) { //when IDs in this range where there
			delete gcontour2d[gc];
			gcontour2d[gc] = NULL;
		}
	}
	//do not clear the entire gcontour2d!

	//##MK::memory leak?
	for ( unsigned int gc = 0; gc < gsegments2d.size(); gc++ ) {
		if ( gsegments2d.at(gc) != NULL ) {
			delete gsegments2d[gc];
			gsegments2d[gc] = NULL;
		}
	}
	//do not clear the entire gsegments2d

#ifdef COMPILE_IN_IF_DESIRED
	for ( unsigned int t = 0; t < trial.size(); t++ ) {
		if ( trial[t] != NULL ) {
			delete trial[t];
			trial[t] = NULL;
		}
	}
#endif

	contours_bary_quick.clear();
	contours_aabb_quick.clear();
	for ( unsigned int rclass = 0; rclass < ChiTable.size(); rclass++ ) {
		if ( ChiTable.at(rclass) != NULL ) {
			delete ChiTable[rclass];
			ChiTable[rclass] = NULL;
		}
	}
	ChiTable.clear();

	myprofiler.logev( "ResetMemory" + std::to_string(fd), (double) (MPI_Wtime() - timer) );
	cout << "...Worker " << this->get_Rank() << " cleared local memory from processing snapshot fid " << fd << endl;
}

//#endif


volanaHdl::volanaHdl( void ) {}


volanaHdl::~volanaHdl( void ) {


	for ( unsigned int snp = 0; snp < AssgnSnapshots.size(); snp++ ) {
		if ( AssgnSnapshots.at(snp) != NULL ) {
			delete AssgnSnapshots[snp];
			AssgnSnapshots[snp] = NULL;
		}
	}

	for ( unsigned int snp = 0; snp < MetaSnapshots.size(); snp++ ) {
		if ( MetaSnapshots.at(snp) != NULL ) {
			delete MetaSnapshots[snp];
			MetaSnapshots[snp] = NULL;
		}
	}

	for ( unsigned int snp = 0; snp < XiTable.size(); snp++ ) {
		if ( XiTable.at(snp) != NULL ) {
			delete XiTable[snp]; //delete on vector calls vector destructor
			XiTable[snp] = NULL;
		}
	}
}


void volanaHdl::set_MPICommunication( int r, int nr ) {
	myRank = r;
	nRanks = nr;
}


void volanaHdl::init_MPIDatatypes( void ) {}


int volanaHdl::get_Rank( void ) {
	return myRank;
}


int volanaHdl::get_nRanks( void ) {
	return nRanks;
}


bool volanaHdl::read_snapshot_volume_mpi_gragles( std::string fn, unsigned int extx, unsigned int exty, unsigned int extz ) {
	//MPI_COMM_SELF read of an implicit raw of size extx*exty*extz in a right-handed coordinate system with x||RD, y||TD, z||ND
	//probe file size of fn to check if consistent
	struct stat buf;
	double filesize[1] = {0.0};
	if ( stat( fn.c_str(), &buf ) != -1 ) {
		filesize[0] = buf.st_size;
		//double maximum_accepted_size = (double) (sizeof(unsigned int) * std::numeric_limits<unsigned int>:: max()) - 2.0;
		double maximum_accepted_size = (double) (std::numeric_limits<unsigned int>:: max()) - 1;
		if ( filesize[0] >= maximum_accepted_size  ) { //file too large
			cerr << "ERROR::File is larger than what current MPI I/O routine is capable of reading in safely!" << endl;
			return false;
		}
		//##MK::rewrite in the future by exchanging nxyz_is to unsigned long and then changing also the blocklength management to unsigned long

		unsigned int nxyz_is = filesize[0] / sizeof(unsigned int); //MK::safe only because of previous check!
		unsigned int nxyz_exp = extx * exty * extz;
		if ( nxyz_is != nxyz_exp ) { cerr << "ERROR::Expected file size and actual size of file " << fn << " do not match!" << endl; return false; }
	}
	else { cerr << "ERROR::File " << fn << " not accessible!" << endl; return false; }

	//MK::file exists so read it in as contiguous blocks total number of data elements < uint32max
	unsigned int TotalBlocksToRead = ( filesize[0] / Settings::MPIReadBlockLength) + 1;
	unsigned int ElementsPerBlock = Settings::MPIReadBlockLength / sizeof(unsigned int);
	unsigned int elementsTotal = filesize[0] / sizeof(unsigned int);
	unsigned int elementsRead = 0;
	unsigned int elementsNow = 0;

	std::vector<unsigned int>* idbuf = NULL; //get large enough results buffer on heap
	try { idbuf = new std::vector<unsigned int>; } 
	catch (std::bad_alloc &exc) { cerr << "ERROR::Allocation error read_snapshot_volume_mpi_gragles!" << endl; 
		return false;
	}

	MPI_File ioReadFileHdl;
	MPI_Status ioReadFileStatus;

	MPI_File_open(MPI_COMM_SELF, fn.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &ioReadFileHdl);
	MPI_File_seek(ioReadFileHdl, 0, MPI_SEEK_SET);

	for ( unsigned int b = 0; b < TotalBlocksToRead; b++ ) {
		elementsNow = ElementsPerBlock;
		if ( (elementsTotal - elementsRead) < elementsNow ) //distribute remaining elements to avoid reading past end of the file
			elementsNow = elementsTotal - elementsRead;
		if ( elementsNow > 0 ) {
			unsigned int* rbuf = NULL;
			try { 
				rbuf = (unsigned int*) new unsigned int[elementsNow];
			} 
			catch (std::bad_alloc &exc) {
				cerr << "ERROR::Allocation error in read_snapshot_volume_mpi_gragles!" << endl;  delete idbuf; return false;
			}

			MPI_File_read( ioReadFileHdl, rbuf, elementsNow, MPI_UNSIGNED, &ioReadFileStatus);

			//explicit datatransfer
			for ( unsigned int e = 0; e < elementsNow; ++e ) {
				idbuf->push_back( rbuf[e] );
			}

			elementsRead = elementsRead + elementsNow;

//cout << "\t\tBlockID/elementsRead/elementsNow/elementsTotal--time = " << b << "/" << elementsRead << "/" << elementsNow << "/" << elementsTotal << "\t\t\t" << (MPI_Wtime() - rtimer) << "\t\tseconds" << endl;

			delete [] rbuf; rbuf = NULL;
		}
	} //data read
	MPI_File_close(&ioReadFileHdl);

	this->AssgnSnapshots.push_back( idbuf );
	return true;
}


bool volanaHdl::read_snapshot_volume( unsigned int fd, unsigned int size ) {
	double timer1 = MPI_Wtime();

	string fname;
	fname = "Container_" + std::to_string(fd) + "_size_" + std::to_string(size-1) + ".raw"; //size-1 because in GraGLeS file name referes to last accessible array index first index 0

	//MPI read of implicit container
	if ( read_snapshot_volume_mpi_gragles( fname, size, size, size ) == false ) {
		cerr << "ERROR::MPI I/O on " << fname << " failed!" << endl;
		return false;
	}

	double timer2 = MPI_Wtime();
	myprofiler.logev( "ReadVolumeData" + std::to_string(fd), (double) (timer2 - timer1) );
	cout << "...Worker " << this->get_Rank() << " " << fname << " volume snapshot MPI I/O successfully in " << (double) (timer2 - timer1) << " seconds" << endl;

	//##MK::add OpenMP thread-local duplicating of the dataset

	return true;
}


bool volanaHdl::read_snapshot_meta_txt_gragles( std::string fn ) {
	//MK::format of Microstructure .uds file
	double _scaler = 1.0;
	if ( Settings::SimModel == E_LEVELSET ) {
		_scaler = Settings::HAGBEnergy / (Settings::PhysicalDomainSize * Settings::DislocEnPerM);
	}
	else { cerr << "ERROR::Scaling is set to default, better reconsider!" << endl; return false; }
	
	ifstream txtfile;
	string txtline;
	istringstream line;
	string datapiece;

	struct UDS_Grain3DInfo ag;
	std::vector<UDS_Grain3DInfo>* mdat = NULL;
	try {
		mdat = new std::vector<UDS_Grain3DInfo>;
	} 
	catch (std::bad_alloc &exc) {
		cerr << "ERROR::Unable to allocate data in read_snapshot_meta_txt_gragles!" << endl; return false;
	}

	double bunge[3] = {0.0, 0.0, 0.0};
	double quat[4] = {1.0, 0.0, 0.0, 0.0};

	txtfile.open( fn.c_str() );
	if ( txtfile.is_open() == true ) { //read in file
		//GraGLeS Microstructure_<fid>.uds has no header!

		while ( txtfile.good() == true ) {
			getline( txtfile, txtline );

			if ( txtline.size() > 0 ) { //not empty
				istringstream line( txtline );

				getline( line, datapiece, '\t' );	ag.gid = atoi( datapiece.c_str() );
				getline( line, datapiece, '\t' ); //##MK::where do the Bunge angle start?
				getline( line, datapiece, '\t' );
				getline( line, datapiece, '\t' );
				getline( line, datapiece, '\t' );	bunge[0] = atof( datapiece.c_str() ); //already in radiant
				getline( line, datapiece, '\t' );	bunge[1] = atof( datapiece.c_str() );
				getline( line, datapiece, '\t' );	bunge[2] = atof( datapiece.c_str() );

				euler2quaternion( bunge, quat );
				ag.q0 = quat[0];
				ag.q1 = quat[1];
				ag.q2 = quat[2];
				ag.q3 = quat[3];

				getline( line, datapiece, '\t' ); //xmin, xmax, ymin, ymax, zmin, zmax
				getline( line, datapiece, '\t' );
				getline( line, datapiece, '\t' );
				getline( line, datapiece, '\t' );
				getline( line, datapiece, '\t' );
				getline( line, datapiece, '\t' );
				getline( line, datapiece, '\t' ); //volume
				getline( line, datapiece, '\t' );	ag.see = atof( datapiece.c_str() ) * _scaler; //MK::GraGLeS outputs scaled dislocation density

				mdat->push_back( ag );
			}
		}
		txtfile.close();

		//MK::IDs are required in ascended order
		std::sort( mdat->begin(), mdat->end(), SortGIDAscending);

		//register dataset
		MetaSnapshots.push_back ( mdat );
		return true;
	}
	else { 
		cerr << "ERROR::Unable to load file " << fn << endl; delete mdat; 
		return false;
	}
}


bool volanaHdl::read_snapshot_meta( unsigned int fd ) {
	double timer1 = MPI_Wtime();

	string fname;
	fname = "Microstructure_" + std::to_string(fd) + ".uds";

	//sequential read of txt file
	if ( read_snapshot_meta_txt_gragles( fname ) == false ) {
		cerr << "ERROR::ASCII I/O on " << fname << " failed!" << endl;
		return false;
	}
	double timer2 = MPI_Wtime();
	myprofiler.logev( "ReadMetaData" + std::to_string(fd), (double) (timer2 - timer1) );
	cout << fname << " metadata snapshot MPI I/O successfully in " << (double) (timer2 - timer1) << " seconds" << endl;

	return true;
}


unsigned int volanaHdl::read_snapshot_rve_extent( unsigned int fd )
{
	//MK::GraGLeS 3d Container_<fid>_size_<extent>.raw <extent> is RVE-size-1
	
	//check existence of candidate file, should only be one file!
	std::string fncand ( "Container_" + std::to_string(fd) + "_size_" );
	std::vector<std::string> dircontent;
	path p( boost::filesystem::current_path() );
	for (auto i = directory_iterator(p); i != directory_iterator(); i++) { //we eliminate directories
		if (!is_directory(i->path()))
			dircontent.push_back( i->path().filename().string() );
		else
			continue;
	}
	//nothing found at all?
	if( dircontent.size() == 0 ) { cerr << "ERROR::No file with prefix " << fncand << " found in current directory!" << endl; return 0; }
	
	//more than one candidate?
	unsigned int candidates = 0;
	unsigned int thatcandidate = 0;
	for (unsigned int fn = 0; fn < dircontent.size(); fn++) {
		std::string str1 ( dircontent.at(fn) );
		if ( dircontent.at(fn).compare(0, fncand.size(), fncand) == 0) { //20
			thatcandidate = fn;
			candidates++;
		}
	}
	//there should not be more than one file with the same prefix fncand but different remaining string <size>.raw!
	if ( candidates > 1 ) { cerr << "ERROR::Target snapshot volume file with prefix " << fncand << " potentially not unique!" << endl; return 0; }

	//okay, only one so read size knowing format of snapshot voolume filenames is fncand + std::to_string(<unsigned int>) + ".raw"
	string szstr = dircontent.at(thatcandidate).substr( fncand.size() );
	szstr.pop_back(); szstr.pop_back(); szstr.pop_back(); szstr.pop_back(); //pop out ".raw"
	return std::stoul( szstr ) + 1;
}


bool volanaHdl::read_snapshot_3d( unsigned int fd )
{
	unsigned int rve_extent = read_snapshot_rve_extent(fd);
	if ( rve_extent < E_LEVELSET_RVE_EXTENT_MIN || rve_extent > E_LEVELSET_RVE_EXTENT_MAX ) {
		cerr << "ERROR::Invalid autodetected RVE extent " << rve_extent << endl; return false;
	}
	
	if ( read_snapshot_volume( fd, rve_extent ) == false ) {
		cerr << "ERROR::Reading of volume data failed!" << endl; return false;
	}
	if ( read_snapshot_meta( fd ) == false ) {
		cerr << "ERROR::Reading of specific time step meta data failed!" << endl; return false;
	}

	//remember this rank works on this fid snapshot
	this->AssgnSnapshotsMeta.push_back( snapshotmeta(fd, rve_extent) );

	return true;
}


void volanaHdl::calc_grainvolume( unsigned int fd )
{
	unsigned int snapshotid = std::numeric_limits<unsigned int>::max();
	for ( unsigned int snp = 0; snp < this->AssgnSnapshotsMeta.size(); snp++ ) {
		if ( this->AssgnSnapshotsMeta.at(snp).fid == fd ) {
			snapshotid = snp;
			break;
		}
	}
	if ( snapshotid == std::numeric_limits<unsigned int>::max() ) { 
		cerr << "ERROR::Worker " << this->get_Rank() << " has not found the local dataset!" << endl; return;
	}

	//clear temporary where we collect the volume for only one fid, namely the one this rank is currently working on
	ghullvolume.clear();
	unsigned int ngr = 1 + Settings::LargestGrainID;
	for ( unsigned int gid = 0; gid < ngr; ++gid ) {
		this->ghullvolume.push_back(0);
	}

	std::vector<unsigned int>* gIDfield = AssgnSnapshots.at(snapshotid);
	unsigned int fieldsize = CUBE(AssgnSnapshotsMeta.at(snapshotid).size);
	for ( unsigned px = 0; px < fieldsize; ++px ) {
		unsigned int gid = gIDfield->at(px);
		ghullvolume.at(gid)++;
	}

	//check how many grainIDs have not at all a single voxel in the container, this is possible, because the discrete sampling of 
	//a LevelSet domain with grains smaller than the volume of a single voxel may cause that grains in the lower tail of the distribution have no voxel
	unsigned int nGrainsWithVolume = 0;
	unsigned int nGrainsWithoutVolume = 0;
	for ( unsigned int gid = (THE_DOMAIN_ITSELF + 1); gid < ngr; ++gid ) {
		if ( ghullvolume.at(gid) > 0 ) { //most likely
			nGrainsWithVolume++;
			continue;
		}
		nGrainsWithoutVolume++;
		//cout << "Grain " << gid << " has no voxel!" << endl;
	}
	cout << "...Worker " << this->get_Rank() << " A total of " << nGrainsWithoutVolume << " out of " << (nGrainsWithVolume+nGrainsWithoutVolume) << " grains remaining were too small of at least one voxel volume!" << endl;
}


void volanaHdl::write_envsphere( unsigned int refgid, unsigned int fid, int bx, int by, int bz, int rd, int xmin, int xmax, int ymin, int ymax, int zmin, int zmax, int size )
{
#define DELTASEE			0x01
#define MOBILITY			0x02
#define MARKOUTSIDE			(1.0e30)
	//output propertiy differences in environment of refgid
	//mode controls which quantity to output
	//mode == DELTASEE, everything not in sphere setting to 1.0e30 1/m^2
	//mode == DISORI, everything not in sphere setting to 1.0e30 rad to threshold later with Paraview
	//write MPI binary file, implicit order with [imin,imax] x lines along +y stacking xy layers along +z
	double qref[4] = { gmetadata.at(refgid).q0, gmetadata.at(refgid).q1, gmetadata.at(refgid).q2, gmetadata.at(refgid).q3 };
	double seeref = gmetadata.at(refgid).see;

	int lim = size;
	int extent[3] = {(xmax-xmin+1),(ymax-ymin+1),(zmax-zmin+1)};
	int nxyz = extent[0]*extent[1]*extent[2]; //<<INT32MAX
	int rdiscr2 = SQR(rd);

	//allocate container
	double * wbuf1 = NULL;
	wbuf1 = new double[nxyz];
	double * wbuf2 = NULL;
	wbuf2 = new double[nxyz];
	double * wbuf3 = NULL;
	wbuf3 = new double[nxyz];

	for ( unsigned int i = 0; i < nxyz; i++ ) {
		wbuf1[i] = MARKOUTSIDE; //i; //1.0e30;
		wbuf2[i] = MARKOUTSIDE; //i; //1.0e30;
		wbuf3[i] = MARKOUTSIDE;
	}

	unsigned int c = 0; //fill implicitly the write buffer wbuf
	int zoff, yzoff, z2, y2z2, x2y2z2;
	for ( int z = zmin; z <= zmax; z++ ) { //checking all cells in axis-aligned cube
		zoff = z*lim*lim;
		z2 = SQR(z - bz);
		for ( int y = ymin; y <= ymax; y++ ) {
			yzoff = y*lim + zoff;
			y2z2 = SQR(y - by) + z2;
			for ( int x = xmin; x <= xmax; x++ ) {
				x2y2z2 = SQR(x - bx) + y2z2;
				if ( x2y2z2 <= rdiscr2 ) { //52% of the case in sphere
					unsigned int there = x+yzoff; //uint to int okay because both x and yzoff > 0
					unsigned int who = AssgnSnapshots.at(0)->at(there);
					//if ( who == refgid ) double halter = 0.0; //##MK::Debug purpose only
					double qwho[4] = { gmetadata.at(who).q0, gmetadata.at(who).q1, gmetadata.at(who).q2, gmetadata.at(who).q3 };
					double seewho = gmetadata.at(who).see;
					//if ( what == DELTASEE ) 
					double seediff = seewho - seeref;
					wbuf1[c] = seediff;
					//if ( what == MOBILITY )
					double theta = misorientationCubic( qref, qwho );
					wbuf2[c] = RAD2DEG(theta); //disorientation angle in degree
					double mob = dis2mob ( theta );
					wbuf3[c] = mob;

//cout << x << ";" << y << ";" << z << "\t\t\t" << who << "\t\t\t" << seeref << ";" << seediff << ";" << theta << ";" << mob << endl;
				} //checked in sphere test
				//else case already handled during wbuf initialization
				//advance implicit file positioner
				c++;
			} // x done
		} // y done
	} // z done

	//write file with MPI
	MPI_File ioHdl1, ioHdl2, ioHdl3;
	MPI_Status ioSta;

	string fn1 = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".Grain." + std::to_string(refgid) + ".FID." + std::to_string(fid) + ".DELTASEE.SIZE." + std::to_string(extent[0]) + ".raw";
	string fn2 = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".Grain." + std::to_string(refgid) + ".FID." + std::to_string(fid) + ".DISORI.SIZE." + std::to_string(extent[0]) + ".raw";
	string fn3 = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".Grain." + std::to_string(refgid) + ".FID." + std::to_string(fid) + ".MOBILITY.SIZE." + std::to_string(extent[0]) + ".raw";
	MPI_File_open(MPI_COMM_SELF, fn1.c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdl1); // open the file in create and write-only mode, only caller writes to file
	MPI_File_open(MPI_COMM_SELF, fn2.c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdl2);
	MPI_File_open(MPI_COMM_SELF, fn3.c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdl3);

	long long totalOffset1 = 0;
	long long totalOffset2 = 0;
	long long totalOffset3 = 0;
	
	MPI_File_write_at( ioHdl1, totalOffset1, wbuf1, nxyz, MPI_DOUBLE, &ioSta);
	MPI_File_write_at( ioHdl2, totalOffset2, wbuf2, nxyz, MPI_DOUBLE, &ioSta);
	MPI_File_write_at( ioHdl3, totalOffset3, wbuf3, nxyz, MPI_DOUBLE, &ioSta);


	MPI_Barrier(MPI_COMM_WORLD); //##MK::not vital
	MPI_File_close(&ioHdl1);
	MPI_File_close(&ioHdl2);
	MPI_File_close(&ioHdl3);

	delete [] wbuf1;
	delete [] wbuf2;
	delete [] wbuf3;
}


void volanaHdl::envcharacterize_omp_naive( snapshotmeta thisone ) {
	//MonteCarlo based approach to compute xi as in the 2D case
	//MK::we work on the unit cube [0,1] \in {\mathbb{R}}^3
	//double halfGbSQR = Settings::DislocEnPerM; //0.5*Gb^2
	//##MK::OMP-parallel version without considering explicit thread local data placement..., therefore naive OMP parallelism...
	double timer = MPI_Wtime();

	unsigned int snapshotid = std::numeric_limits<unsigned int>::max();
	for ( unsigned int snp = 0; snp < this->AssgnSnapshotsMeta.size(); snp++ ) {
		if ( this->AssgnSnapshotsMeta.at(snp).fid == thisone.fid ) {
			snapshotid = snp;
			break;
		}
	}
	if ( snapshotid == std::numeric_limits<unsigned int>::max() ) { cerr << "ERROR::Worker " << this->get_Rank() << " has not found the local dataset!" << endl; return; }

	//obtain radius of inspection sphere, which is fixed for each grain and time step
	//MK::because of grid coarsening in GraGLeS the discrete representation of this physical radius changes though, now take this into account...
	double GridCoarsedSize = (double) thisone.size;
	int lim = thisone.size; //unsigned to int, the pixel coordinate lim is not accessible [0,lim)
	//by translating relative extent of inspection sphere about each grain into current size radius of the discrete sphere
	double scc = Settings::LongRangeRadiusMax * GridCoarsedSize; //MK:: the parameter LongRangeRadius is \in (0,1) \in \mathbb{R} i.e. it is a relative quantity, relative to Settings::PhysicalDomainLength;, so even if the discretization in GraGLeS reduces the smaller the number of grains we inspect the same physical environment
	int rmax = scc; //how many pixel is the radius of the largest inspection sphere?
	if ( 2*rmax >= (lim-1) ) { cerr << "ERROR::Analysis is definately impossible for all grains as LongRangeRadiusMax*DomainEdgeLength is too large!" << endl; return; }
	scc = Settings::LongRangeRadiusMin * GridCoarsedSize;
	int rmin = scc;
	if( rmin < 1 ) { cerr << "ERROR::Analysis is definately impossible as LongRangeRadiusMin*DomainEdgeLength is too small!" << endl; return; }
	
	//even though the analysis is now in principle possible we do not track currently over periodic boundaries and thus have to accept that
	//the total number of grains for which the inspection volume does not protrude outside the RVE may be very small...

	//initiates global XiTable
	unsigned int rc = 0;
	for ( double r = Settings::LongRangeRadiusMin; r <= Settings::LongRangeRadiusMax; r = r + Settings::LongRangeRadiusIncr ) { //implicit order of results array xi
		std::vector<xi>* axi = NULL;
		try { 
			axi = new vector<xi>;
			XiTable.push_back( axi );
			rc++;
		} 
		catch (std::bad_alloc &exc) {
			cerr << "ERROR::Allocation error for XiTable in envcharacterize_omp_naive!" << endl;
			for ( unsigned int i = 0; i < XiTable.size(); i++ ) { delete XiTable.at(i); } 
			return;
		}		
	}

	#pragma omp parallel //without further specifying snpashotid, GirdCoarsedSize, lim, rmax, rc are shared
	{
		unsigned int nthreads = omp_get_num_threads();
		unsigned int mythread = omp_get_thread_num();

		//threads allocate local memory for results in local memory
		std::vector<std::vector<xi>*> LocalXiTable; //container of pointer only

		//they prepare their LocalXiTable
		unsigned int rclass = 0;
		for ( double r = Settings::LongRangeRadiusMin; r <= Settings::LongRangeRadiusMax; r = r + Settings::LongRangeRadiusIncr ) {
			std::vector<xi>* axi = NULL;
			try { 
				axi = new vector<xi>;
				LocalXiTable.push_back( axi );
				rclass++;		
			}
			catch (std::bad_alloc &exc) { cerr << "ERROR::Allocation error for thread local xi memory in envcharacterize_longrange_core!" << endl; }
		}

		//#pragma omp barrier
		//##MK::barrier not necessary because LocalXiTable is threadlocal
		
		std::vector<unsigned int>* this_snapshot = this->AssgnSnapshots.at(snapshotid);
		unsigned int ngr = 1 + Settings::LargestGrainID;
		for ( unsigned int gid = 1; gid < ngr; gid++ ) { //loop over all grain ids, process existing
			if ( gid != THE_DOMAIN_ITSELF ) { //##MK::filter out a particular grain && gd == 113590 ) {
				if ( WorkPartitioning(gid, nthreads, mythread) == true ) { //if I take care of this grain
					if ( this->ghullvolume.at(gid) > 0 ) { //grain has still measureable volume
					
						//double gtimer = omp_get_wtime();
						//quick screening to include only grains whose spherical investigation area falls completely inside the RVE domain
						//map initial AABB center <=> initial sub-grain center from resolution InitialDomainEdgeLength in pixel on GridCoarsedStructure
						double mapping[3] = { (double) gmetadata.at(gid).x0, (double) gmetadata.at(gid).y0, (double) gmetadata.at(gid).z0 }; //MK::gmetadata works as a hash as it was fed with the initial structure and hence contains all grains!
						mapping[0] /= Settings::InitialDomainEdgeLength; //absolute pixel coordinate initial structure into relative coordinate in initial domain, MK::must not be thisone.size!
						mapping[1] /= Settings::InitialDomainEdgeLength;
						mapping[2] /= Settings::InitialDomainEdgeLength;
						mapping[0] *= GridCoarsedSize; //relative position times pixel results again in pixel in double representation
						mapping[1] *= GridCoarsedSize;
						mapping[2] *= GridCoarsedSize;
					
						int bx0 = mapping[0]; //casting required to get pixel, ##MK::coordinate wise casting from uint to int okay, because values << INT32_MAX stop always accessing this as it causes remote memory access..
						int by0 = mapping[1];
						int bz0 = mapping[2];

						if ( bx0 < rmax )			continue; //position \pm radius protruding outside domain, next gid!
						if ( bx0 > (lim-rmax-1) )	continue;
						if ( by0 < rmax )			continue;
						if ( by0 > (lim-rmax-1) )	continue;
						if ( bz0 < rmax )			continue;
						if ( bz0 > (lim-rmax-1) )	continue;

						//okay, process this grain
						double qme[4] = { gmetadata.at(gid).q0, gmetadata.at(gid).q1, gmetadata.at(gid).q2, gmetadata.at(gid).q3 };
						double seeme = gmetadata.at(gid).see;

						rclass = 0; //##MK::that can be optimized

						/*////DEBUG
						//which coordinates is the occupied by AssgnSnapshot->at(0)
						for ( unsigned int zz = 0; zz < lim; zz++ ) {
							for ( unsigned int yy = 0; yy < lim; yy++ ) {
								for ( unsigned int xx = 0; xx < lim; xx++ ) {
									if ( AssgnSnapshots.at(0)->at(xx+(yy*lim)+(zz*lim*lim)) == gd )
										cout << gd << "--" << xx << "\t\t" << yy << "\t\t" << zz << endl;
						}}}
						////DEBUG*/

						for ( double r = Settings::LongRangeRadiusMin; r <= Settings::LongRangeRadiusMax; r = r + Settings::LongRangeRadiusIncr ) { //loop over radii
							//double gtimer = omp_get_wtime();
							struct xi mp;
							mp.gid = gid;
							mp.x = bx0;		mp.y = by0;		mp.z = bz0;		mp.r = r;
							unsigned int ct = 0; //total pixel cnt in inspection sphere, therefrom covered
							unsigned int cs = 0; //how many of these are covered by the target itself
							unsigned int cn = 0; //cnt nbors

							double sc = r * GridCoarsedSize;
							int rdiscr = sc; //##MK::practically now 3D simulation domain for grain growth extxi >= INT32_RANGE
							int rdiscr2 = SQR(rdiscr); //MK::to avoid multiplications and sqrt //##MK::safely less than INT32_MAX
							int xmi, xmx, ymi, ymx, zmi, zmx;
							xmi = bx0 - rdiscr; //global pixel coordinates on GridCoarsedDomain
							xmx = bx0 + rdiscr;
							ymi = by0 - rdiscr;
							ymx = by0 + rdiscr;
							zmi = bz0 - rdiscr;
							zmx = bz0 + rdiscr; //all these ints are >= 0 !

							//MK::DEBUG write_envsphere( gid, thisone.fid, bx0, by0, bz0, rdiscr, xmi, xmx, ymi, ymx, zmi, zmx, lim );

							//std::vector<mind>* cache = NULL;
							//cache = new std::vector<mind>; //accounting for disjoint grains found in inspection volume
							//std::vector<unsigned int>* grid = NULL; //reference do not delete
							//grid = &(this->AssgnSnapshots.at(0));
							std::vector<mind> cache;

							unsigned int who;
							bool found = false; //status flag indicating if we have found the grain at position x,y,z
							int zoff, yzoff, z2, y2z2;
							for ( int z = zmi; z <= zmx; z++ ) { //checking all cells in axis-aligned cube about interspection sphere
								zoff = z*lim*lim;
								z2 = SQR(z - bz0);
								for ( int y = ymi; y <= ymx; y++ ) {
									yzoff = y*lim + zoff;
									y2z2 = SQR(y - by0) + z2;
									for ( int x = xmi; x <= xmx; x++ ) {
										if ( (SQR(x - bx0) + y2z2) <= rdiscr2 ) { //52% of the cases in sphere
											ct++;
											//who = grid[x+yzoff];
											who = this_snapshot->at(x+yzoff);
											if ( who != gid ) { //not me, most likely
												found = false;
												for ( unsigned int i = 0; i < cache.size(); i++ ) {
													if ( cache[i].gid != who ) //not found the correct grain to cache most likely the cnt?
														continue;

													//not continued? cached hid-id who already
													found = true;
													cache[i].cnt += 1;
													break; //do not search any longer
												}
												if ( found == false ) { cache.push_back( mind( who, 1) ); }//a previously not encountered neighbor id
											} //no, in fact me!
											else { cs++; }
										} //done in sphere test
									} //x line done
								} // y slab done
							} // z done, entire inspection done


							if ( cache.size() > 0 ) {
								for (unsigned int j = 0; j < cache.size(); j++) { cn += cache[j].cnt; }
								if ( (cn + cs) != ct ) { cerr << "ERROR::Worker " << this->get_Rank() << " analysis for grain " << gid << " inconsistent!" << endl; }
							}

							mp.Vtotal = ct;		mp.Vself = cs;		mp.Vnbors = cn;
							//report only for gid grains that still exists, i.e. have at least a single cell
							if ( cs > 0 ) {
								if ( cache.size() > 0 ) { //accumulate results from neighbors neighbors
									//MK::the key time saving is here that we calculate the disori to each neighbor only once instead of as many times as it shares cells in the inspection sphere
									double xi_drv = 0.0;
									double xi_spd = 0.0;
									double xi_dis = 0.0;
									double FractionalApproxVolNeighbors = 0.0;
									for ( unsigned int j = 0; j < cache.size(); j++ ) { //found somebody, perform computation of xi value
										unsigned int hd = cache[j].gid;
										double FractionalApproxVolInSphere = cache[j].cnt; //implicit uint32 to double no problem
										FractionalApproxVolNeighbors += FractionalApproxVolInSphere;
										double qcand[4] = { gmetadata.at(hd).q0, gmetadata.at(hd).q1, gmetadata.at(hd).q2, gmetadata.at(hd).q3 };
										double seecand = gmetadata.at(hd).see;
										double mob = dis2mob(misorientationCubic(qme, qcand));
										//MK::(seecnad-seeme) positive means SEE-driven growth, negative SEE-induced rather the neighbor consumes me
										xi_drv = xi_drv + (Settings::DislocEnPerM * (seecand - seeme) * FractionalApproxVolInSphere);
										xi_spd = xi_spd + (Settings::HAGBMobility * mob * Settings::DislocEnPerM*(seecand - seeme) * FractionalApproxVolInSphere);
										xi_dis = xi_dis + (Settings::HAGBMobility * mob * FractionalApproxVolInSphere);
									}
									if ( FractionalApproxVolNeighbors > 0.1 ) { //##MK::Counts must be at least 1 cell...
										xi_drv /= FractionalApproxVolNeighbors;
										xi_spd /= FractionalApproxVolNeighbors;
										xi_dis /= FractionalApproxVolNeighbors;
									}
									mp.Pdrv = xi_drv;
									mp.Pspd = xi_spd;
									mp.Pdis = xi_dis;
								}
								else { //myself larger than the discrete sphere
									mp.Pdrv = 0.0;
									mp.Pspd = 0.0;
									mp.Pdis = 0.0;
								}

								LocalXiTable.at(rclass)->push_back( mp );
							} //end of reporting existent gd grain

							cache.clear();
							rclass++;

							/*#pragma omp critical
							{
								cout << "---->" << mp.gid << " grain volume in Container " << this->ghullvolume.at(mp.gid) << ";" << mp.x << ";" << mp.y << ";" << mp.z << ";" << mp.r << ";" << mp.Vtotal << ";" << mp.Vself << ";"<< mp.Vnbors << ";" << mp.r << "\t\t" << setprecision(16) << mp.Pdrv << "\t\t" << mp.Pspd << " took " << setprecision(8) << (double) (omp_get_wtime() - gtimer ) << " seconds " << endl;
							}*/

						} //next radius
						//##MK::one could gain significant improvements by checking the inclusion only once for the largest radius, write into a buffer of included cells and successively test only those
					} //done with this grain
				} //grain with volume > 0
			} //executing thread done with grain g
		} //end of processing grains

		#pragma omp barrier
		//MK::necessary, no coordinated writing of results when not all are done, there is still improvement potential here...

		#pragma omp critical
		{
			//threads copy local results
			unsigned int rcc = 0;
			for ( double r = Settings::LongRangeRadiusMin; r <= Settings::LongRangeRadiusMax; r = r + Settings::LongRangeRadiusIncr ) { //implicit order
				for ( unsigned int localg = 0; localg < LocalXiTable.at(rcc)->size(); localg++ ) {
					//only grains with volume > 0 are collected
					//MK::the order of the grainIDs in the list is non-deterministic between groups of grains from different threads to avoid sorting as the point in time when the threads enter the critical region may vary
					XiTable.at(rcc)->push_back( LocalXiTable.at(rcc)->at(localg) );
				}

				//threads delete their local data
				delete LocalXiTable.at(rcc); LocalXiTable.at(rcc) = NULL;
				rcc++;
			}
		}
		LocalXiTable.clear();

	} //end of parallel region

	string tname = "EnvCharacterizeOMP" + std::to_string(thisone.fid);
	myprofiler.logev( tname, (double) (MPI_Wtime() - timer) );

	cout << "...Worker " << this->get_Rank() << " analysis on " << thisone.fid << " took " << setprecision(8) << (double) (MPI_Wtime() - timer) << " seconds" << endl;
}


void volanaHdl::init_discrete_version( void ) {
	//init contour buckets, as ID hashs
	double timer = MPI_Wtime();

	struct isectionmeta dummy;
	unsigned int ngr = 1 + Settings::LargestGrainID;
	for ( unsigned int gid = 0; gid < ngr; gid++ ) {
		gmetadata.push_back( dummy );
	}

	myprofiler.logev( "InitializeMemory", (double) (MPI_Wtime() - timer) );
}


bool volanaHdl::read_microstructure_uds_3d( void ) {
	//reading an UDS File Microstructure.uds for computation unbiased chi values by adding the orientation and SEE of the grains that which is not part of the GNU file
	//ASCII file with a header comprising 6 lines in 2D and 8 lines in 3D utilizing only either empty lines (LF or CRLF) or tab separated strings or a single string
	if ( Settings::Dimensionality == TWO_DIMENSIONS ) {
		cerr << "ERROR::3D UDS Microstructure.uds required!" << endl;
		return false;
	}

	ifstream udsfile;
	string udsline;
	istringstream line;
	string datapiece;

	struct isectionmeta dat;
	
	udsfile.open( Settings::UDSDataFromFilename );
	if ( udsfile.is_open() == true ) { //read in file, which has no header, line-by-line
		//identify header first
		unsigned int nheaderlines = 1 + 1 + 3 + 3 + 1 + 1; //MK::THREE_DIMENSIONS

//jump over header
unsigned int j = 0;
while (udsfile.good() == true && j < nheaderlines) {
	getline(udsfile, udsline);
	j++;
}

//begin reading metadata
while (udsfile.good() == true) {
	getline(udsfile, udsline);

	if (udsline.size() > 0) { //not an empty line
		istringstream line(udsline);
		//interpret
		getline(line, datapiece, '\t');	dat.gid = atoi(datapiece.c_str());
		getline(line, datapiece, '\t'); dat.x0 = atoi(datapiece.c_str()); //initial reference position as center of AABB
		getline(line, datapiece, '\t'); dat.y0 = atoi(datapiece.c_str()); //these coordinates are pixel and have to be translated into a relative coordinate by dividing with Settings::InitialDomainEdgeLength
		getline(line, datapiece, '\t'); dat.z0 = atoi(datapiece.c_str());

		double bunge[3];
		getline(line, datapiece, '\t'); bunge[0] = atof(datapiece.c_str());
		getline(line, datapiece, '\t'); bunge[1] = atof(datapiece.c_str());
		getline(line, datapiece, '\t'); bunge[2] = atof(datapiece.c_str());
		double quat[4];
		euler2quaternion(bunge, quat);
		dat.q0 = quat[0];		dat.q1 = quat[1];		dat.q2 = quat[2];		dat.q3 = quat[3];

		getline(line, datapiece, '\t'); //xmin, xmax, ymin, ymax
		getline(line, datapiece, '\t');
		getline(line, datapiece, '\t');
		getline(line, datapiece, '\t'); //zmin, zmax
		getline(line, datapiece, '\t');
		getline(line, datapiece, '\t');
		getline(line, datapiece, '\t'); //vol
		getline(line, datapiece, '\t'); dat.see = atof(datapiece.c_str()); //MK::SEE with no scaling

		/* ##MK::add these to parse rest of the line of a MicrostructureDiagnostics.uds file
		//only for MicrostructureDiagnostics.uds to eliminate the trailing detailed infos
		//eliminate the rest 10 columns, rather crappily
		getline(line, datapiece, '\t');
		getline(line, datapiece, '\t');
		getline(line, datapiece, '\t');
		getline(line, datapiece, '\t');
		getline(line, datapiece, '\t');
		getline(line, datapiece, '\t');
		getline(line, datapiece, '\t');
		getline(line, datapiece, '\t');
		getline(line, datapiece, '\t');
		getline(line, datapiece, '\t');
		*/

		//cout << dat.gid << "\t" << bunge[0] << ";" << bunge[1] << ";" << bunge[2] << "\t" << dat.q0 << ";" << dat.q1 << ";" << dat.q2 << ";" << dat.q3 << ";" << dat.see << endl;
		//cout << dat.gid << ";" << dat.see << endl;
		//invokes copyconstructor for stl types
		if (gmetadata.at(dat.gid).gid == THE_DOMAIN_ITSELF) { //replace dummy data in existing hashtable
			gmetadata.at(dat.gid) = dat;
		}
		else {
			cerr << "ERROR::Inconsistency namely the same ID twice! " << dat.gid << endl;
			udsfile.close();
			return false;
		}
	} //next line...
} //...if still good
udsfile.close();
	}
	else {
		cerr << "ERROR::Unable to load file " << Settings::UDSDataFromFilename << endl; return false;
	}

	return true;
}


void volanaHdl::read_initialmicrostructure_metadata(void) {
	//MK::reads a Microstructure_<fd>.uds ASCII file to get the mapping of grainIDs to grain properties
	double timer = MPI_Wtime();

	if (read_microstructure_uds_3d() == false) { cerr << "ERROR::Unable to locate Microstructure.uds file" << endl; return; }

	myprofiler.logev("ReadMetadataInitially", (double)(MPI_Wtime() - timer));
	cout << "...Worker " << this->get_Rank() << " " << Settings::UDSDataFromFilename << " grain meta data were loaded successfully in " << (double)(MPI_Wtime() - timer) << " seconds" << endl;
}


void volanaHdl::report(snapshotmeta thisone)
{
	double timer = MPI_Wtime();

	string prefix;
	string postfix;
	string logfname;
	ofstream logfile;
	prefix = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".LR3D.FID." + std::to_string(thisone.fid) + ".RCLASS.";
	postfix = ".csv";

	unsigned int rclass = 0;
	for (double r = Settings::LongRangeRadiusMin; r <= Settings::LongRangeRadiusMax; r = r + Settings::LongRangeRadiusIncr) { //loop over radii classes
		logfname = prefix + std::to_string(rclass) + postfix;
		logfile.open(logfname.c_str(), std::ofstream::out | std::ofstream::trunc);
		if (logfile.is_open() == true) {
			logfile << "GrainID;x;y;z;r;Vtotal;Vself;Vnbors;VolAvEffSignedDrivingForcePrognosis;VolAvEffSignedMigrationSpeedPrognosis;VolAvEffMobilityEnv\n";
			logfile << "Grain;px;px;px;px;px;px;px;J/m^3;m/s;m^4/Js\n";
			logfile << "RClass=" << r << ";x;y;z;r;Vtotal;Vself;Vnbors;VolAvEffSignedDrivingForcePrognosis;VolAvEffSignedMigrationSpeedPrognosis;VolAvEffMobilityEnv\n";

			std::vector<xi>* ths = XiTable.at(rclass);
			for (unsigned long i = 0; i < ths->size(); i++) {
				if (ths->at(i).Vnbors > 0) {
					logfile << ths->at(i).gid << ";" << setprecision(16) << ths->at(i).x << ";" << ths->at(i).y << ";" << ths->at(i).z << ";" << ths->at(i).r << ";";
					logfile << ths->at(i).Vtotal << ";" << ths->at(i).Vself << ";" << ths->at(i).Vnbors << ";";
					logfile << setprecision(18) << ths->at(i).Pdrv << ";" << setprecision(18) << ths->at(i).Pspd << ";" << ths->at(i).Pdis <<"\n";
				}
				//else, grain consumed entire inspection volume already so Pdrv = Pspd = 0.0 by definition!
			}
			logfile << endl;
			logfile.flush();
			logfile.close();
		}
		else {
			cerr << "ERROR::Worker::VolAnaHdl " << this->get_Rank() << " unable to write report" << endl;
		}
		rclass++;
	}
	cout << "...Worker " << this->get_Rank() << " reported the results successfully in " << (double) (MPI_Wtime() - timer) << endl;

	string tname = "Report" + std::to_string(thisone.fid);
	myprofiler.logev( tname, (double) (MPI_Wtime() - timer) );
}


void volanaHdl::metareport( void ) {
	string logfname;
	ofstream logfile;
	logfname = "LRAnalysis3D.SimID." + to_string(Settings::SimID) + ".DummyMetaReport.csv";
	logfile.open ( logfname.c_str(), std::ofstream::out | std::ofstream::trunc );

	if (logfile.is_open() == true) {
		logfile << "Grain;x0;y0;see\n";
		logfile << "Grain;1;1;Pa\n";
		logfile << "Grain;x0;y0;see\n";
		for (unsigned int g = 0; g < gmetadata.size(); g++) {
			logfile << endl;
		}
		logfile.flush();
		logfile.close();
	}
	else {
		cerr << "ERROR:Worker::VolAnaHdl " << this->get_Rank() << " unable to write metareport" << endl;
	}
}


void volanaHdl::reset_xitable( void )
{
	for ( unsigned int rcc = 0; rcc < XiTable.size(); rcc++ ) {
		if ( XiTable.at(rcc) != NULL ) {
			XiTable.at(rcc)->clear();
		}
	}
}


void volanaHdl::spit_profiling( void )
{
	//first overview of profiling for sinlge process execution...
	string logfname;
	logfname = "TopologyTracer.LR3D.SimID." + std::to_string(Settings::SimID) + ".Rank." + std::to_string(this->get_Rank()) + ".Profiling.csv";
	ofstream logfile;
	logfile.open( logfname.c_str(), std::ofstream::out | std::ofstream::trunc );
	if (logfile.is_open() == true) {
		//header
		logfile << "WhatTask;WallClock\n";
		logfile << ";s\n";
		logfile << ";MPI_Wtime\n";

		for (unsigned int l = 0; l < myprofiler.get_nentries(); l++) {
			logfile << myprofiler.titles[l] << ";" << setprecision(8) << myprofiler.times[l] << "\n";
		}
		logfile.flush();
		logfile.close();
	}
	else {
		cerr << "ERROR::Worker::VolAnaHdl " << this->get_Rank() << " unable to write profiling data" << endl;
	}
}





#include "thirdparty/Eigen/Dense"
using namespace Eigen;


curvapprxHdl::curvapprxHdl(){}


curvapprxHdl::~curvapprxHdl()
{
	unsigned int ng = 1 + Settings::LargestGrainID;
	for ( unsigned int g = 0; g < ng; g++ ) {
		if ( gcontour2d[g] != NULL ) {
			delete gcontour2d[g];
		}
	}
	gcontour2d.clear();

	//upon destruction of analyzer instance we can also clean the main data
	unsigned int ni = gresults_data2d.size();
	for ( unsigned int i = 0; i < ni; i++ ) {
		if ( gresults_data2d[i] != NULL ) {
			delete gresults_data2d[i];
		}
	}
	gresults_data2d.clear();
	gresults_meta2d.clear();
}


void curvapprxHdl::set_MPICommunication( int r, int nr ) {
	myRank = r;
	nRanks = nr;
}


void curvapprxHdl::init_MPIDatatypes( void )
{
	//MPI_Datatype MPI_GBContourPoint_Type;
	int elementCounts1[2] = {4, 2};
	MPI_Aint displacements1[2] = {0, 4 * MPIIO_OFFSET_INCR_DOUBLE};
	MPI_Datatype oldTypes1[2] = {MPI_DOUBLE, MPI_UNSIGNED};
	MPI_Type_create_struct(2, elementCounts1, displacements1, oldTypes1, &MPI_GBContourPoint_Type);

	//MPI_Datatype MPI_GBJunctionPoint_Type;
	int elementCounts2[2] = {2, 2};
	MPI_Aint displacements2[2] = {0, 2 * MPIIO_OFFSET_INCR_DOUBLE};
	MPI_Datatype oldTypes2[2] = {MPI_DOUBLE, MPI_UNSIGNED};
	MPI_Type_create_struct(2, elementCounts2, displacements2, oldTypes2, &MPI_GBJunctionPoint_Type);

	MPI_Type_commit(&MPI_GBContourPoint_Type);
	MPI_Type_commit(&MPI_GBJunctionPoint_Type);
}


int curvapprxHdl::get_Rank( void ) {
	return myRank;
}

int curvapprxHdl::get_nRanks( void ) {
	return nRanks;
}


bool curvapprxHdl::init( void )
{
	//init contour buckets, as ID hashs
	double timer = MPI_Wtime();

	unsigned int n = 1 + Settings::LargestGrainID;
	for ( unsigned int c = 0; c < n; c++ ) {
		gcontour2d.push_back( NULL );
	}

	for ( unsigned int c = 0; c < n; c++ ) {
		//replace NULL pointer with pointer to array of contour points per grain
		std::vector<GBContourPoint>* cpbucket = NULL;
		cpbucket = new vector<GBContourPoint>;
		if ( cpbucket == NULL ) {
			cerr << "ERROR::Worker " << this->get_Rank() << " was unable to initialize working array!" << endl;
			return false;
		}

		gcontour2d.at(c) = cpbucket;
	}

	myprofiler.logev( "InitializeMemory", (double) (MPI_Wtime() - timer) );
	return true;
}


bool curvapprxHdl::read_snapshot_contour_mpi_gragles( std::string fn ) 
{
	//MPI_COMM_SELF read of an implicit binary of GBContourPoint structs
	//probe file size of fn to check if consistent
	struct stat buf;
	double filesize[1] = {0.0};
	unsigned long fsz = 0;
	unsigned long szcp = sizeof(struct GBContourPoint);
	if ( stat( fn.c_str(), &buf ) != -1 ) {
		filesize[0] = buf.st_size;

		fsz = filesize[0];
		if ( fsz % szcp != 0 ) { //file size is not a multiple of sizeof GBContourPoint, indicates file is potentially corrupted!
			cerr << "ERROR::File potentially corrupted!" << endl; 
			return false;
		}
	}
	else { cerr << "ERROR::File size of " << fn.c_str() << " not detectable!" << endl; return false; } //filesize not detectable

	//file exists, size is known so read it in as contiguous blocks
	//MK::MPI library only read blocks of maximum size at once but contour files can be larger than 4GB
	//no global buffering of the file content required, MK::extremely beneficial because reducing memory overhead...
	MPI_File ioReadFileHdl;
	MPI_Status ioReadFileStatus;

	MPI_File_open(MPI_COMM_SELF, fn.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &ioReadFileHdl);
	MPI_File_seek(ioReadFileHdl, 0, MPI_SEEK_SET); //set implicit fp to beginning of file

	bool stillgood = true;
	unsigned long elementsRead = 0; //MK::not unsigned int because file potentially more than UINT32MAX elements
	unsigned long elementsTotal = fsz / szcp;
	unsigned int elementsNow = 0; //read in blocks of struct GBContourPoint beigin at offset=0, MPI_CHAR native interpretation
	unsigned int elementsPerBlock = Settings::MPIReadBlockLength / sizeof(struct GBContourPoint);

	while ( (elementsRead < elementsTotal) && stillgood == true ) {
		//formulate a single block read via MPI i/O
		elementsNow = elementsPerBlock; //MK::not unsigned long , compatibility with MPI library

		 //handle remaining potentially shorter rest of the data
		if ( (elementsTotal - elementsRead) < ((unsigned long) elementsNow) ) {
			unsigned long diff = elementsTotal - elementsRead;
			if ( diff > std::numeric_limits<unsigned int>::max() ) {
				cerr << "ERROR::Remaining block handling on " << fn.c_str() << " failed!" << endl;
				stillgood = false;
				break; //MK::do not return false, because then ioReradFileHdl is not closed!
			}
			elementsNow = (unsigned int) diff; //MK::downcasting from unsigned long to unsigned int is now safe
		}

		if ( stillgood == true ) {
			//set up a buffer to store MPI i/O result
			GBContourPoint* rbuf = NULL;
			rbuf = new GBContourPoint[elementsNow];
			if ( rbuf == NULL ) { cerr << "ERROR::Allocation error in read_snapshot_contour_mpi_gragles!" << endl; 
				stillgood = false; break;
			}

			MPI_File_read( ioReadFileHdl, rbuf, elementsNow, MPI_GBContourPoint_Type, &ioReadFileStatus);
			//implicitly advancing fp

			//interpret data immediately to gcontour2d
			unsigned gid = 0;
			for ( unsigned int e = 0; e < elementsNow; ++e ) {
				gid = rbuf[e].myid;
				if ( gid != THE_DOMAIN_ITSELF && gid <= Settings::LargestGrainID ) { //##MK::debug safety to catch inconsistencies
					//this loading imprints an explicit order on the grain-local container preserving the order of GraGLeS

					//MK::mind that the algorithm requires the point data to be sorted either clock- or anti-clockwise!
					gcontour2d.at(gid)->push_back( GBContourPoint( rbuf[e].x, rbuf[e].y, rbuf[e].energy, rbuf[e].mobility, rbuf[e].myid, rbuf[e].nborid ) );
				}
			}//buffer processed

			elementsRead = elementsRead + elementsNow;

			delete [] rbuf; rbuf = NULL;
		}
	}

	MPI_File_close(&ioReadFileHdl);
	return stillgood;
}


bool curvapprxHdl::read_snapshot_contour( std::string prefix, unsigned int fid )
{
	double timer1 = MPI_Wtime();
	string fname;
	fname = prefix + "_" + std::to_string(fid) + ".bin";

	//MPI read of implicit container
	if ( read_snapshot_contour_mpi_gragles( fname ) == false ) {
		cerr << "ERROR::MPI I/O on " << fname << " failed!" << endl;
		return false;
	}
	double timer2 = MPI_Wtime();
	myprofiler.logev( "ReadContourPointData" + std::to_string(fid), (double) (timer2 - timer1) );
	cout << "...Worker " << this->get_Rank() << " loaded contour snapshot " << fid << " via MPI I/O successfully in " << (double) (timer2 - timer1) << " seconds" << endl;

	if ( Settings::TranslateBinary2GNU == true ) {
		write_snapshot_contour_gnu_gragles( fid );
	}

	modelspecific_rawdata_modification();

	return true;
}


void curvapprxHdl::write_snapshot_contour_gnu_gragles( unsigned int fd ) 
{
	//MK::Original GraGLeS GNU File format requires for each grain as the last point a duplicate of the first point to be written followed by a blank line
	double timer1 = MPI_Wtime();
	string gnu_fn;
	gnu_fn = "Network_Timestep_" + std::to_string(fd) + ".gnu";

	ofstream gnu;
	gnu.open( gnu_fn.c_str(), std::ofstream::out | std::ofstream::trunc );
	if (gnu.is_open() == true) {
		unsigned int ngr = 1 + Settings::LargestGrainID;
		unsigned int ncp = 0;
		std::vector<GBContourPoint>* these = NULL;
		for (unsigned int gid = 1; gid < ngr; gid++) {
			these = gcontour2d.at(gid);
			ncp = these->size();
			if (ncp > 0) {
				for (unsigned int p = 0; p < ncp; p++) { //x, y, //relative mobility
					gnu << setprecision(8) << these->at(p).x << "\t" << these->at(p).y << "\t" << setprecision(4) << these->at(p).mobility << "\n";
				}
				gnu << "\n"; //blank line
			}
		} //next grain

		gnu.flush();
		gnu.close();
	}
	else {
		cerr << "ERROR::Worker::CurvApproxHdl " << this->get_Rank() << " unable to write to gnufile" << endl;
	}

	double timer2 = MPI_Wtime();
	myprofiler.logev( "ReportingResults" + std::to_string(fd), (double) (timer2 - timer1) );
	cout << "...Worker " << this->get_Rank() << " performed translation of " << fd << " *.bin to *.gnu successfully in " << (double) (timer2 - timer1) << " seconds" << endl;
}


void curvapprxHdl::modelspecific_rawdata_modification( void )
{
	//MK::this function serves to implement user-defined rescaling operations on the raw data
	double rtimer = MPI_Wtime();

	if ( Settings::SimModel == E_LEVELSET ) {
		//MK::modify contour data by eliminating last point
		unsigned int ngr = 1 + Settings::LargestGrainID;
		unsigned int ncp = 0;

		for ( unsigned int gid = 1; gid < ngr; gid++ ) {
			ncp = this->gcontour2d.at(gid)->size();
			if ( ncp > 4 ) { //most likely case, most grains have at least three triple junction supporting points + duplicate of the last point
				//MK::assume no change of output order from GraGLeS was introduced by reading the intput
				//then we can eliminate the last point as it is a duplicate
				this->gcontour2d.at(gid)->pop_back();
			}
			//cout << "Popped back was/now = " << ncp << "\t\t" << this->gcontour2d.at(gid)->size() << endl;
		}

		//MK::map all relative point coordinates into absolute physical coordinates
		for ( unsigned int gid = 1; gid < ngr; gid++ ) {
			ncp = this->gcontour2d.at(gid)->size();
			if ( ncp > 0 ) {
				std::vector<GBContourPoint>* thesepoints = this->gcontour2d.at(gid);
				for ( unsigned int p = 0; p < ncp; p++ ) {
					thesepoints->at(p).x *= METER2MICRON(Settings::PhysicalDomainSize);
					thesepoints->at(p).y *= METER2MICRON(Settings::PhysicalDomainSize);					
				}
			}
		}
	}
	cout << "...Worker " << this->get_Rank() << " performed user-defined unit scaling operations on raw data within " << setprecision(6) << (MPI_Wtime() - rtimer) << " seconds" << endl;
}


#define WITH_TJP		0
#define WITHOUT_TJP		1
struct wghtd_imcurv_res curvapprxHdl::ompcrit_imean_curv_core( unsigned int g )
{
	//MUST BE CALLED FROM WITHIN PARALLEL REGION BUT ONLY BY A SINGLE THREAD AT A TIME!
	//compute integral mean curvature of grain with contour this->gcontour2d.at(g)->
	vector<GBContourPoint>* cp = this->gcontour2d.at(g);

	//bool curve_orientation_is_clockwise = true; //MK::GraGLeS outputs contour points in clockwise order, verify by det[(1,xa,ya,1,xb,yb,1,xc,yc)^T] < 0.0 i.e. negative

	//MK::assumes contour points are unique! and in particualr the last point duplicate was popped back already, i.e. no duplicates!
	unsigned int ni = cp->size();
	unsigned int pre; //contour point indices: preceeding, ...
	unsigned int curr;
	unsigned int post;
	unsigned int pre_nbid; //utilize change in neighbor id to detect proximity to triple point and therefore exclusion from the dataset
	unsigned int curr_nbid; 
	unsigned int post_nbid;

	//point not next to triple line
	unsigned int nallvertices = 0;
	unsigned int nvirtualvertices = 0;
	double clen_curr = 0.0; //current segment length
	double wlen_total[2] = {0.0, 0.0};
	double drvf_total[2] = {0.0, 0.0};
	double spd_total[2] = {0.0, 0.0};
	double imcurv[2] = {0.0, 0.0}; //only integral mean curvature for the grain kappa * clen_curr

	//follow the path of contour points and approximate integral mean curvature by pointwise fitting of circle to cp not at junctions
	bool success = true; //can we compute values for all points or get other inconsistencies?
	for ( unsigned int i = 0; i < ni; i++ ) { //access to ni is faulty!
		pre = (i > 0) ? i-1 : ni-1;
		curr = i;
		post = (i < ni-1) ? i+1 : 0;

		//get local conditions at points 0,1,2 follow along contour
		double p0x = cp->at(pre).x;
		double p0y = cp->at(pre).y;
		pre_nbid = cp->at(pre).nborid;
		double p1x = cp->at(curr).x;
		double p1y = cp->at(curr).y;
		curr_nbid = cp->at(curr).nborid;
		double p2x = cp->at(post).x;
		double p2y = cp->at(post).y;
		post_nbid = cp->at(post).nborid;
		
		//compute length of bisector pre-curr and curr_post
		double pre_curr = 0.5 * sqrt(SQR(p1x-p0x) + SQR(p1y-p0y));  //##MK::optimization for sure possible, aka reutilization of values, etc...
		double curr_post = 0.5 * sqrt(SQR(p2x-p1x) + SQR(p2y-p1y));

		//effective segment length on which local curvature is assumed acting
		clen_curr = pre_curr + curr_post;
		
		//compute local curvature by fitting circle to p0,p1,p2 to obtain principal curvature radius approximation, solve equation system with Eigen
		double kappa = 0.0;
		Matrix3d A;
		A << p0x, p0y, 1.0, p1x, p1y, 1.0, p2x, p2y, 1.0;
		Vector3d b;
		b << -1.0*( SQR(p0x)+SQR(p0y) ), -1.0*( SQR(p1x)+SQR(p1y) ), -1.0*( SQR(p2x)+SQR(p2y) );
		Vector3d x;
		x << A.jacobiSvd(ComputeFullU|ComputeFullV).solve(b);

		vector<double> v2;
		v2.resize(x.size());
		VectorXd::Map(&v2[0],x.size()) = x;

		//only take into account if computation of integral mean curvature is possible
		if ( std::isnan(v2[0]) == false && std::isnan(v2[1]) == false && std::isnan(v2[2]) == false ) { //all values are not NaN, most likely
			double test = (SQR(v2[0])+SQR(v2[1]))/4.0 - v2[2];

			if ( test > 1.0e-14 ) { //sqrt(test) numerically safe, if small sqrt(test)--> 0.0 anyway so curvature diverge to INFTY
				kappa = fabs(1.0 / sqrt(test)); //always positive
				
				//identify local concavity via sign of detO
				double detO = +1.0*(p1x*p2y - p2x*p1y) -1.0*(p0x*p2y - p2x*p0y) +1.0*(p0x*p1y - p1x*p0y);
					
				//MK::GraGLeS clockwise, i.e. if det[(1,xf,yf,xg,yg,xh,yh)^T] < 0.0 segment at xg,yg is convex, if == 0, collinear, if > 0.0 concave
				//we consider capillary driving forces that assuming them positive when they cause the grain to shrink
				//i.e. approximating a circle as an n-fold regular polyhrdron with opening angles n/(2pi) the circle 
				//experiences a positive driving force causing its shrinkage, hence for the clockwise defined n-fold polyhedron
				//at all points we are convex, i.e. det < 0.0, thus, if detO < 0.0, as GraGLeS runs clockwise
				//we count the segment to be convex experiencing positive driving force
				//instead if detO is > 0.0 driving force is counted for as positive
				if ( detO > 0.0 ) { //opposite local concavity as clockwise counting (detO < 0.0, convex, positive force)
					//so segment locally concave --> driving force negative
					kappa *= -1.0;
				}
			}
			else { //##MK::technically curvature INFINTE, do not consider and dont worry about concav convex...
				//leave kappa zero, as initialized
				kappa = 0.0;
			}
		
			//we now have a local curvature so we can add it but have to distinguish whether or not considering the vertex to be close to tjp or not
			nallvertices++;
			wlen_total[WITH_TJP] = wlen_total[WITH_TJP] + clen_curr;
			drvf_total[WITH_TJP] = drvf_total[WITH_TJP] + (cp->at(curr).energy * kappa * clen_curr);
			spd_total[WITH_TJP] = spd_total[WITH_TJP] + (cp->at(curr).mobility * cp->at(curr).energy * kappa * clen_curr);
			imcurv[WITH_TJP] = imcurv[WITH_TJP] + (kappa * clen_curr); //accumulated local curvature values times approximate length on which they act into integral mean curvature

			if ( pre_nbid == curr_nbid && curr_nbid == post_nbid ) { //contour point is a "virtual vertex" discretizing curvature sufficiently far from junction
				nvirtualvertices++;
				wlen_total[WITHOUT_TJP] = wlen_total[WITHOUT_TJP] + clen_curr;
				drvf_total[WITHOUT_TJP] = drvf_total[WITHOUT_TJP] + (cp->at(curr).energy * kappa * clen_curr); //[spd] = 1*1/micron*micron
				spd_total[WITHOUT_TJP] = spd_total[WITHOUT_TJP] + (cp->at(curr).mobility * cp->at(curr).energy * kappa * clen_curr ); //[drvf] = 1*1*1/micron*micron
				imcurv[WITHOUT_TJP] = imcurv[WITHOUT_TJP] + (kappa * clen_curr); //[kappa] = 1/micron, [clen_curr] = micron, i.e. 
			}	
		}
		else { //potential flaw in the computation of the curvature, hence, better send error value
			nallvertices = 0;				nvirtualvertices = 0;		clen_curr = 0.0;
			wlen_total[WITH_TJP] = 0.0;		wlen_total[WITHOUT_TJP] = 0.0;
			drvf_total[WITH_TJP] = 0.0;		drvf_total[WITHOUT_TJP] = 0.0;
			spd_total[WITH_TJP] = 0.0;		spd_total[WITHOUT_TJP] = 0.0;
			imcurv[WITH_TJP] = 0.0;			imcurv[WITHOUT_TJP] = 0.0;
			success = false;
			break;
		} //eventually stop computation because of numerical flaws
	} //analyze next contour pointnext virtual vertex

	//are values for both WITH_TJP and WITHOUT_TJP computable at all?
	if ( nallvertices == 0 || wlen_total[WITH_TJP] <= DBL_EPSILON )
		success = false;
	if ( nvirtualvertices == 0 || wlen_total[WITHOUT_TJP] <= DBL_EPSILON )
		success = false;

	wghtd_imcurv_res ar; //prepare output set error values proformer
	ar.imcurv_with_tjp = std::numeric_limits<double>::max();
	ar.imcurv_without_tjp = std::numeric_limits<double>::max();
	ar.effcapindspeed_with_tjp = std::numeric_limits<double>::max();
	ar.effcapindspeed_without_tjp = std::numeric_limits<double>::max();
	ar.effcapdrvforce_with_tjp = std::numeric_limits<double>::max();
	ar.effcapdrvforce_without_tjp = std::numeric_limits<double>::max();
	ar.cp_with_tjp = nallvertices;
	ar.cp_without_tjp = nvirtualvertices;
	ar.gID = g;
	ar.status = success;

	if ( success == true ) { //do we have at all a trustworthy account of data?
		spd_total[WITH_TJP] = spd_total[WITH_TJP] / wlen_total[WITH_TJP];
		drvf_total[WITH_TJP] = drvf_total[WITH_TJP] / wlen_total[WITH_TJP];
		//MK::mind that cp->at(curr).energy is scaled to Settings::HAGBEnergy and cp->at(curr).mobility is scale to Settings::HAGBMobility
		ar.effcapindspeed_with_tjp = spd_total[WITH_TJP] * 1.0e6 * Settings::HAGBMobility * Settings::HAGBEnergy; //kappa was in 1/micron
		ar.effcapdrvforce_with_tjp = drvf_total[WITH_TJP] * 1.0e6 * Settings::HAGBEnergy;
		ar.imcurv_with_tjp = imcurv[WITH_TJP];
				
		spd_total[WITHOUT_TJP] = spd_total[WITHOUT_TJP] / wlen_total[WITHOUT_TJP];
		drvf_total[WITHOUT_TJP] = drvf_total[WITHOUT_TJP] / wlen_total[WITHOUT_TJP];
		ar.effcapindspeed_without_tjp = spd_total[WITHOUT_TJP] * 1.0e6 * Settings::HAGBMobility * Settings::HAGBEnergy;
		ar.effcapdrvforce_without_tjp = drvf_total[WITHOUT_TJP] * 1.0e6 * Settings::HAGBEnergy;
		ar.imcurv_without_tjp = imcurv[WITHOUT_TJP];
	}

	return ar;
}


bool curvapprxHdl::integral_mean_curvature_approximation( unsigned int fd )
{
	Eigen::setNbThreads(1); //disable potential threading of Eigen for subsequent calls
	double timer1 = MPI_Wtime();
	unsigned int ngr = 1 + Settings::LargestGrainID;

	//MK::find correct position on gresults_data2d as ranks partition processing of fid on [SnapshotFirst::Offset::SnapshotLast] potentially arbitrarily
	unsigned int thiscandidate = UNKNOWN_CANDIDATE;
	for ( unsigned int cand = 0; cand < this->gresults_meta2d.size(); cand++ ) {
		if ( this->gresults_meta2d[cand].fid != fd ) { continue; } 
		//implicit else
		thiscandidate = cand;
	}
	if ( thiscandidate == UNKNOWN_CANDIDATE ) { cerr << "ERROR::Worker " << this->get_Rank() << " unable to find local gresults_data2d index!" << endl; return false; }

	//determine correct link to results container
	std::vector<wghtd_imcurv_res>* rbucket = this->gresults_data2d.at(thiscandidate);
	
	//prepare processing by assigning dummy results to all grains
	double dblerrval = std::numeric_limits<double>::max();
	unsigned int uinterrval = std::numeric_limits<unsigned int>::max();
	for ( unsigned int gid = 0; gid < (1+Settings::LargestGrainID); gid++ ) { //master generates result
		rbucket->push_back( wghtd_imcurv_res(dblerrval, dblerrval, dblerrval, dblerrval, dblerrval, dblerrval, uinterrval, uinterrval, gid, false) );
	} //now for each grain a result exists but the values are flagged as still errorneous/unprocessed

	//process letting threads process curvature values in parallel and writeback on shared link rbucket
	#pragma omp parallel shared(ngr, rbucket)
	{
		unsigned int mythreadid = omp_get_thread_num();
		unsigned int nthreads = omp_get_num_threads();
		//std::vector<wghtd_imcurv_res>* rbucket = this->gresults_data2d.at(this->gresults_data2d.size()-1);

		for ( unsigned int gid = 0; gid < ngr; gid++ ) {
			if ( gid % nthreads == mythreadid ) { //I process
				if ( gcontour2d.at(gid)->size() > 4 ) { //contour data for grain gid exist at all
					//wghtd_imcurv_res tmp = ompcrit_imean_curv_core( gid );
					rbucket->at(gid) = ompcrit_imean_curv_core( gid ); //no race because threads write disjoint items!
					//cout << "Grain " << gid << " sum i=1^N m_i*gamma_i*L_i*kappa_i / sum i=1^N m_i*gamma_i*L_i m/s = " << tmp.cp_without_tjp << "\t\t" << setprecision(32) << tmp.effcapindspeed_without_tjp << endl;
					//rbucket->push_back( tmp );
				}
			} //next mygrain
		} //next grain
	} //parallel region end

	double timer2 = MPI_Wtime();
	myprofiler.logev( "ProcessOMPContourData" + std::to_string(fd), (double) (timer2 - timer1) );
	cout << "...Worker " << this->get_Rank() << " performed OMP-parallelized contour data processing for " << fd << " successfully in " << (double) (timer2 - timer1) << " seconds" << endl;
	return true;
}


void curvapprxHdl::report_wghtd_imcurv( std::string prefix, unsigned int fd, bool clear_rbucket )
{
	double timer1 = MPI_Wtime();
	//find correct output container to read from
	std::vector<wghtd_imcurv_res>* rbucket = NULL;
	unsigned int thiscandidate = UNKNOWN_CANDIDATE;
	for ( unsigned int cand = 0; cand < this->gresults_meta2d.size(); cand++ ) {
		if ( this->gresults_meta2d[cand].fid != fd ) { continue; } 
		//implicit else
		thiscandidate = cand;
	}
	if ( thiscandidate == UNKNOWN_CANDIDATE ) { cerr << "ERROR::Worker " << this->get_Rank() << " unable to find local gresults_data2d index!" << endl; return; }
	rbucket = this->gresults_data2d.at(thiscandidate);

	string curvlog_fn;
	curvlog_fn = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + "." + prefix + ".FID." + std::to_string(fd)  + ".CurvApprx.csv";

	ofstream curvlog;
	curvlog.open( curvlog_fn.c_str(), std::ofstream::out | std::ofstream::trunc );
	if (curvlog.is_open() == true) {
		curvlog << "GrainID;effCapillaryIndDrvForceWithoutTJP;effCapIndDrvForceWithTJP;effSpeedWithoutTJP;effSpeedWithTJP;intMeanCurvWithoutTJP;intMeanCurvWithTJP;nVirtualVertices;nAllVertices\n";
		curvlog << "1;Pa=N/m^2=Nm/m^3=J/m^3;Pa=J/m^3;m/s;m/s;micron*1/micron;micron*1/micron;1;1\n";

		for (unsigned int i = 0; i < rbucket->size(); i++) { //only output grains with valid computation result
			if (rbucket->at(i).status == true) { //no error occurred
				curvlog << rbucket->at(i).gID << ";" << setprecision(8) << rbucket->at(i).effcapdrvforce_without_tjp << ";" << rbucket->at(i).effcapdrvforce_with_tjp << ";";
				curvlog << rbucket->at(i).effcapindspeed_without_tjp << ";" << rbucket->at(i).effcapindspeed_with_tjp << ";";
				curvlog << rbucket->at(i).imcurv_without_tjp << ";" << rbucket->at(i).imcurv_with_tjp << ";";
				curvlog << rbucket->at(i).cp_without_tjp << ";" << rbucket->at(i).cp_with_tjp << "\n";
			}
		}
		curvlog.flush();
		curvlog.close();
	}
	else {
		cerr << "ERROR::Worker::CurvApproxHdl " << this->get_Rank() << " unable to write report_wght_imcurv" << endl;
	}

	if ( clear_rbucket == true ) {
		rbucket->clear();
	}

	double timer2 = MPI_Wtime();
	myprofiler.logev( "ReportingResults" + std::to_string(fd), (double) (timer2 - timer1) );
	cout << "...Worker " << this->get_Rank() << " performed result I/O for " << fd << " successfully in " << (double) (timer2 - timer1) << " seconds" << endl;
}


void curvapprxHdl::reset_memory( void )
{
	unsigned int ng = 1 + Settings::LargestGrainID;
	for ( unsigned int g = 0; g < ng; g++ ) {
		if ( gcontour2d[g] != NULL ) { //do not delete the vector but clear content of existing vector
			gcontour2d[g]->clear();
		}
	}

	//##MK::do not clear the results to allow filling potentially large results table...
}




curvfuseHdl::curvfuseHdl( void )
{
	//generate hash of dumpbuffer data sufficient to fastly index grains
	for ( unsigned int gid = 0; gid <= Settings::LargestGrainID; ++gid ) {
		dumpbuffer.push_back( capdump() );
	}
}


curvfuseHdl::~curvfuseHdl()
{
	//GrainIDWhiteList cleared automatically
	//dumpbuffer cleared automatically
}


void curvfuseHdl::set_MPICommunication( int r, int nr ) {
	myRank = r;
	nRanks = nr;
}


void curvfuseHdl::init_MPIDatatypes( void )
{
	//MPI_Datatype MPI_CapillaryDump_Type;
	int elementCounts1[2] = {3, 4};
	MPI_Aint displacements1[2] = {0, 3 * MPIIO_OFFSET_INCR_DOUBLE};
	MPI_Datatype oldTypes1[2] = {MPI_DOUBLE, MPI_UNSIGNED};
	MPI_Type_create_struct(2, elementCounts1, displacements1, oldTypes1, &MPI_CapillaryDump_Type);

	MPI_Type_commit(&MPI_CapillaryDump_Type);
}


int curvfuseHdl::get_Rank( void ) {
	return myRank;
}

int curvfuseHdl::get_nRanks( void ) {
	return nRanks;
}


bool curvfuseHdl::read_targets_fromlist( void )
{
	double timer = MPI_Wtime();

	GrainIDWhiteList.clear();

	ifstream trgfile;
	string trgline;
	istringstream line;
	string datapiece;

	trgfile.open( Settings::TargetGIDFromFilename );
	if ( trgfile.is_open() == true ) {
		//format of the file is only IDs, headerless, single column
		//single header
		getline( trgfile, trgline );
		
		while ( trgfile.good() == true ) {
			getline( trgfile, trgline );
			unsigned int s = trgline.length();

			if ( s > 0 ) { //only a target id 
				long tgid = std::stol( trgline ); //accept only gIDs on (THE_DOMAIN_ITSELF,Settings::LargestGrainID]
				if ( tgid > THE_DOMAIN_ITSELF && tgid <= Settings::LargestGrainID )
					GrainIDWhiteList.push_back( tgid );
				else {
					cerr << "ERROR::Attempting to feed an invalid grainID via TargetGrainIDs supplementary file!" << endl; trgfile.close(); return false;
				}
			}
			//ignore blank lines
		}
		trgfile.close();
		cout << "...Worker " << this->get_Rank() << " loaded " << GrainIDWhiteList.size() << " valid user-selected grainIDs to work on selectively in " << (double) (MPI_Wtime() - timer) << " seconds" << endl;
		return true;
	}
	return false;
}

bool curvfuseHdl::init_targets( void )
{
	//double timer = MPI_Wtime();

	GrainIDWhiteList.clear();
	for ( unsigned int tgid = THE_DOMAIN_ITSELF; tgid <= Settings::LargestGrainID; ++tgid ) {
		GrainIDWhiteList.push_back( tgid );
	}
	return true;
}


/*
unsigned int curvfuseHdl::check_dumpdata_existence( void )
{
	//##MK:: implement at last, for now assume all files exist!

	//1 on success, 0 otherwise
	return 1;
}
*/


void curvfuseHdl::reset_dumpdata_todefaults( void )
{
	size_t ngid = dumpbuffer.size();
	for ( size_t gid = 0; gid < ngid; ++gid ) {
			dumpbuffer.at(gid) = capdump();
	}
}


unsigned int curvfuseHdl::init_localbuffer( void )
{
	//##MK::change to map in the future
	unsigned int ni = GrainIDWhiteList.size();
	for ( unsigned int i = 0; i < ni; ++i ) {
		localacc.push_back( capav(GrainIDWhiteList.at(i)) );
	}

	if ( localacc.size() != GrainIDWhiteList.size() )
		return 0;
	else 
		return 1;
}


unsigned int curvfuseHdl::init_results_file( void )
{
	//matrix of fixed number of quantiles for all timesteps
	int ncol = ((Settings::SnapshotLast - Settings::SnapshotFirst) / Settings::SnapshotOffset + 1); //##MK:: in the future1 + timesteps because add quantile arguments as well
	int nrow = GrainIDWhiteList.size(); //as many as there are targets
	cout << "...MASTER " << this->get_Rank() << " Writing rows (number of targets) = " << nrow << " columns (time) = " << ncol << endl;

	MPI_File ioHdlICURV, ioHdlPCURV, ioHdlVCURV, ioHdlTSUPP, ioHdlVSUPP;
	//MPI_Status ioStaICURV, ioStaPCURV, ioStaVCURV, ioHdlTSUPP, ioHdlVSUPP;

	string prefix = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".CapTracking"; 
	string suffix = std::to_string(Settings::SnapshotFirst) + ".O." + std::to_string(Settings::SnapshotOffset) + ".L." + std::to_string(Settings::SnapshotLast) + ".NC." + std::to_string(ncol) + ".NR." + std::to_string(nrow) + ".bin";
	string icurv_fn = prefix + ".ICURV." + suffix;
	string pcurv_fn = prefix + ".PCURV." + suffix;
	string vcurv_fn = prefix + ".VCURV." + suffix;
	string tsupp_fn = prefix + ".TSUPP." + suffix;
	string vsupp_fn = prefix + ".VSUPP." + suffix;

	// open the file in create and write-only mode but do not write data
	MPI_File_open(MPI_COMM_SELF, icurv_fn.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlICURV);
	MPI_File_close(&ioHdlICURV);

	MPI_File_open(MPI_COMM_SELF, pcurv_fn.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlPCURV);
	MPI_File_close(&ioHdlPCURV);

	MPI_File_open(MPI_COMM_SELF, vcurv_fn.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlVCURV);
	MPI_File_close(&ioHdlVCURV);

	MPI_File_open(MPI_COMM_SELF, tsupp_fn.c_str(),  MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlTSUPP);
	MPI_File_close(&ioHdlTSUPP);

	MPI_File_open(MPI_COMM_SELF, vsupp_fn.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlVSUPP);
	MPI_File_close(&ioHdlVSUPP);

	//1 on success, 0 otherwise
	return 1;
}


bool curvfuseHdl::readDatasets( unsigned int fid )
{
	//MPI_COMM_SELF read of a binary GBCurvatureApproximation_<fid>.bin struct files
	//these are headerless, a container of (3*double,4*unsigned int), little endian
	string fn = "GBCurvatureApproximation_" + std::to_string(fid) + ".bin";

	//probe file size of fn to check if consistent
	struct stat buf;
	double filesize[1] = {0.0};
	unsigned long fsz = 0;
	unsigned long szcp = (3*8) + (4*4);
	if ( stat( fn.c_str(), &buf ) != -1 ) {
		filesize[0] = buf.st_size;

		fsz = filesize[0];
		if ( fsz % szcp != 0 ) { //file size is not a multiple of expected struct size, indicates file is potentially corrupted!
			cerr << "ERROR::File potentially corrupted!" << endl; 
			return false;
		}
	}
	else { cerr << "ERROR::File size of " << fn.c_str() << " not detectable!" << endl; return false; } //filesize not detectable

	//sets all hash values to sensible defaults such that no longer existent grains can be filtered out easily..
	reset_dumpdata_todefaults();

	//file exists, size is known so read it in as contiguous blocks
	//MK::MPI library only read blocks of maximum size at once but contour files can be larger than 4GB
	//no global buffering of the file content required, MK::extremely beneficial because reducing memory overhead...
	MPI_File ioReadFileHdl;
	MPI_Status ioReadFileStatus;

	MPI_File_open(MPI_COMM_SELF, fn.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &ioReadFileHdl);
	MPI_File_seek(ioReadFileHdl, 0, MPI_SEEK_SET); //set implicit fp to beginning of file

	bool stillgood = true;
	unsigned long elementsRead = 0; //MK::not unsigned int because file potentially more than UINT32MAX elements
	unsigned long elementsTotal = fsz / szcp;
	unsigned int elementsNow = 0; //read in blocks of struct GBContourPoint beigin at offset=0, MPI_CHAR native interpretation
	unsigned int elementsPerBlock = Settings::MPIReadBlockLength / ((3*8) + (4*4));

	while ( (elementsRead < elementsTotal) && stillgood == true ) {
		//formulate a single block read via MPI i/O
		elementsNow = elementsPerBlock; //MK::not unsigned long , compatibility with MPI library

		 //handle remaining potentially shorter rest of the data
		if ( (elementsTotal - elementsRead) < ((unsigned long) elementsNow) ) {
			unsigned long diff = elementsTotal - elementsRead;
			if ( diff > std::numeric_limits<unsigned int>::max() ) {
				cerr << "ERROR::Remaining block handling on " << fn.c_str() << " failed!" << endl;
				stillgood = false;
				break; //MK::do not return false, because then ioReradFileHdl is not closed!
			}
			elementsNow = (unsigned int) diff; //MK::downcasting from unsigned long to unsigned int is now safe
		}

		if ( stillgood == true ) {
			//set up a buffer to store MPI i/O result
			MPI_CapillaryDump* rbuf = NULL;
			rbuf = new MPI_CapillaryDump[elementsNow];
			if ( rbuf == NULL ) { cerr << "ERROR::Allocation error in curvfuseHdl::readDatasets!" << endl; 
				stillgood = false; break;
			}

			MPI_File_read( ioReadFileHdl, rbuf, elementsNow, MPI_CapillaryDump_Type, &ioReadFileStatus);
			//implicitly advancing fp

			//interpret data immediately to gcontour2d
			unsigned gid = 0;
			for ( unsigned int e = 0; e < elementsNow; ++e ) {
				gid = rbuf[e].gid;
				if ( gid != THE_DOMAIN_ITSELF && gid <= Settings::LargestGrainID ) { //##MK::debug safety to catch inconsistencies
					//this loading imprints an explicit order on the grain-local container preserving the order of GraGLeS

					//MK::mind that the algorithm requires the point data to be sorted either clock- or anti-clockwise!
					dumpbuffer.at(gid) = capdump(rbuf[e].imcurv, rbuf[e].pcurv, rbuf[e].vcurv, rbuf[e].nvt, rbuf[e].nvv);
				}
			}//buffer processed

			elementsRead = elementsRead + elementsNow;

			delete [] rbuf; rbuf = NULL;
		}
	}

	MPI_File_close(&ioReadFileHdl);
	return stillgood;
}


bool curvfuseHdl::workpartitioning_fid(unsigned int fid)
{
	//decide if this->get_Rank() processes or not
	unsigned int f = 1 + (fid / Settings::SnapshotOffset) - (Settings::SnapshotFirst / Settings::SnapshotOffset);
	if (f % this->get_nRanks() == this->get_Rank())
		return true;
	else
		return false;
}



void curvfuseHdl::fuse_capillary_activity( unsigned int fid)
{
	//matrix of fixed number of quantiles for all timesteps and processes
	int ncol = ((Settings::SnapshotLast - Settings::SnapshotFirst) / Settings::SnapshotOffset + 1); //##MK:: in the future1 + timesteps because add quantile arguments as well
	int nrow = GrainIDWhiteList.size(); //as many as there are targets

	//generate output buffer
	//##MK::quick-and-dirty only one single write of entire buffer
	//##MK::better segment into pieces as in readDatasets...
	double* ic = NULL;
	try { ic = new double[nrow]; }
	catch (std::bad_alloc &exc) {
		std::cout << "ERROR::CapillaryTracking allocation error when attempting writing results for snapshot " << fid << " with ic buffer!" << endl;
	}

	double* pc = NULL;
	try { pc = new double[nrow]; }
	catch (std::bad_alloc &exc) {
		std::cout << "ERROR::CapillaryTracking allocation error when attempting writing results for snapshot " << fid << " with pc buffer!" << endl;
		if ( ic != NULL ) { delete [] ic; ic = NULL; }
		return;
	}

	double* vc = NULL;
	try { vc = new double[nrow]; }
	catch (std::bad_alloc &exc) {
		std::cout << "ERROR::CapillaryTracking allocation error when attempting writing results for snapshot " << fid << " with vc buffer!" << endl;
		if ( ic != NULL ) { delete [] ic; ic = NULL; }
		if ( pc != NULL ) { delete [] pc; pc = NULL; }
		return;
	}

	unsigned int* tsupp = NULL;
	try { tsupp = new unsigned int[nrow]; }
	catch (std::bad_alloc &exc) {
		std::cout << "ERROR::CapillaryTracking allocation error when attempting writing results for snapshot " << fid << " with tsupp buffer!" << endl;
		if ( ic != NULL ) { delete [] ic; ic = NULL; }
		if ( pc != NULL ) { delete [] pc; pc = NULL; }
		if ( vc != NULL ) { delete [] vc; vc = NULL; }
		return;
	}

	unsigned int* vsupp = NULL;
	try { vsupp = new unsigned int[nrow]; }
	catch (std::bad_alloc &exc) {
		std::cout << "ERROR::CapillaryTracking allocation error when attempting writing results for snapshot " << fid << " with vsupp buffer!" << endl;
		if ( ic != NULL ) { delete [] ic; ic = NULL; }
		if ( pc != NULL ) { delete [] pc; pc = NULL; }
		if ( vc != NULL ) { delete [] vc; vc = NULL; }
		if ( tsupp != NULL ) { delete [] tsupp; tsupp = NULL; }
		return;
	}

	//MK::set sensible defaults to avoid writing undefined data for died out grains
	double dblerror = std::numeric_limits<double>::max();
	for ( unsigned int i = 0; i < nrow; ++i ) { 
		ic[i] = dblerror;
		pc[i] = dblerror;
		vc[i] = dblerror;
		tsupp[i] = 0;
		vsupp[i] = 0;
	}

	//fill buffers with targets values
	unsigned int ntargets = this->GrainIDWhiteList.size();
	for ( unsigned int tgr = 0; tgr < ntargets; ++tgr ) {
		unsigned int which = GrainIDWhiteList.at(tgr);
		capdump these = dumpbuffer.at( which );
		if ( these.nvt != 0 ) { //MK::nvt is upon default construction initialized to 0 serving as a flag
			ic[tgr] = these.imcurv;
			pc[tgr] = these.pcurv;
			vc[tgr] = these.vcurv;
			tsupp[tgr] = these.nvt;
			vsupp[tgr] = these.nvv;
		}
	}

	//write results to file via MPI_COMM_SELF parallel I/O at fixed a priori known position in results matrix
	MPI_File ioHdlICURV, ioHdlPCURV, ioHdlVCURV, ioHdlTSUPP, ioHdlVSUPP;
	MPI_Status ioStaICURV, ioStaPCURV, ioStaVCURV, ioStaTSUPP, ioStaVSUPP;

	string prefix = "TopoTracer2D3D.SimID." + std::to_string(Settings::SimID) + ".CapTracking"; 
	string suffix = std::to_string(Settings::SnapshotFirst) + ".O." + std::to_string(Settings::SnapshotOffset) + ".L." + std::to_string(Settings::SnapshotLast) + ".NC." + std::to_string(ncol) + ".NR." + std::to_string(nrow) + ".bin";
	string icurv_fn = prefix + ".ICURV." + suffix;
	string pcurv_fn = prefix + ".PCURV." + suffix;
	string vcurv_fn = prefix + ".VCURV." + suffix;
	string tsupp_fn = prefix + ".TSUPP." + suffix;
	string vsupp_fn = prefix + ".VSUPP." + suffix;

	unsigned int f = (fid / Settings::SnapshotOffset) - (Settings::SnapshotFirst / Settings::SnapshotOffset);
	long long totalOffset = 0;

	// open the file in write-only mode
	MPI_File_open(MPI_COMM_SELF, icurv_fn.c_str(), MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlICURV);
	totalOffset = f * nrow * 8;
	MPI_File_write_at(ioHdlICURV, totalOffset, ic, nrow, MPI_DOUBLE, &ioStaICURV);
	MPI_File_close(&ioHdlICURV);

	MPI_File_open(MPI_COMM_SELF, pcurv_fn.c_str(), MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlPCURV);
	totalOffset = f * nrow * 8;
	MPI_File_write_at(ioHdlPCURV, totalOffset, pc, nrow, MPI_DOUBLE, &ioStaPCURV);
	MPI_File_close(&ioHdlPCURV);

	MPI_File_open(MPI_COMM_SELF, vcurv_fn.c_str(), MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlVCURV);
	totalOffset = f * nrow * 8;
	MPI_File_write_at(ioHdlVCURV, totalOffset, vc, nrow, MPI_DOUBLE, &ioStaVCURV);
	MPI_File_close(&ioHdlVCURV);

	MPI_File_open(MPI_COMM_SELF, tsupp_fn.c_str(), MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlTSUPP);
	totalOffset = f * nrow * 4;
	MPI_File_write_at(ioHdlTSUPP, totalOffset, tsupp, nrow, MPI_UNSIGNED, &ioStaTSUPP);
	MPI_File_close(&ioHdlTSUPP);

	MPI_File_open(MPI_COMM_SELF, vsupp_fn.c_str(), MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdlVSUPP);
	totalOffset = f * nrow * 4;
	MPI_File_write_at(ioHdlVSUPP, totalOffset, vsupp, nrow, MPI_UNSIGNED, &ioStaVSUPP);
	MPI_File_close(&ioHdlVSUPP);

	if ( ic != NULL ) { delete [] ic; ic = NULL; }
	if ( pc != NULL ) { delete [] pc; pc = NULL; }
	if ( vc != NULL ) { delete [] vc; vc = NULL; }
	if ( tsupp != NULL ) { delete [] tsupp; tsupp = NULL; }
	if ( vsupp != NULL ) { delete [] vsupp; vsupp = NULL; }
}


void curvfuseHdl::localaveraging_capillary_activity( unsigned int fid) 
{
	//accumulate values in localbuffer
	unsigned int ntargets = this->GrainIDWhiteList.size();
	for ( unsigned int tgr = 0; tgr < ntargets; ++tgr ) { //localacc[tgr] accumulates for grain ID ..WhiteList.at(tgr)
		unsigned int gid = GrainIDWhiteList.at(tgr);
		capdump these = dumpbuffer.at( gid ); //a hash table
		if ( these.nvt != 0 ) { //MK::nvt is upon default construction initialized to 0 serving as a flag whether grain exists or not
			localacc[tgr].add( these );
		}
	}
}


void curvfuseHdl::spit_average_results(void)
{
	//we know every process has vector of imcurv
	//##MK::quick and dirty, on runtime error checks..
	MPI_CapillaryDump* sndbuf = NULL;
	MPI_CapillaryDump* rcvbuf = NULL;

	unsigned int ntargets = GrainIDWhiteList.size();
	if ( this->get_nRanks() > 1 ) {
		if ( this->get_Rank() == MASTER ) {
			rcvbuf = new MPI_CapillaryDump[ntargets];
			for ( unsigned int trg = 0; trg < ntargets; ++trg ) { //master clears buffer to avoid numerical data corruption
				rcvbuf[trg].imcurv = 0.0;
				rcvbuf[trg].pcurv = 0.0;
				rcvbuf[trg].vcurv = 0.0;
				rcvbuf[trg].nvt = 0;
				rcvbuf[trg].nvv = 0;
				rcvbuf[trg].gid = 0;
				rcvbuf[trg].padding = 0;
			}
		}
		else {
			sndbuf = new MPI_CapillaryDump[ntargets];
			for ( unsigned int trg = 0; trg < ntargets; ++trg ) { //slaves also fill buffer with local values
				sndbuf[trg].imcurv = localacc.at(trg).imcurv;
				sndbuf[trg].pcurv = localacc.at(trg).pcurv;
				sndbuf[trg].vcurv = localacc.at(trg).vcurv;
				sndbuf[trg].nvt = localacc.at(trg).nvt;
				sndbuf[trg].nvv = localacc.at(trg).nvv;
				sndbuf[trg].gid = GrainIDWhiteList.at(trg);
				sndbuf[trg].padding = localacc.at(trg).n;
			}
		}
		//wait for buffers to be ready, potentially overlapping possible because only master rcv buffer and r first snd have to be allocated...

		MPI_Barrier( MPI_COMM_WORLD );

		//aggregation only necessary if more than one process, successively aggregate results on buffers of same length in master
		for ( unsigned int r = 1; r < this->get_nRanks(); ++r ) { //start at one because master process has already his local results
			if ( this->get_Rank() == MASTER ) {
				MPI_Recv(rcvbuf, ntargets, MPI_CapillaryDump_Type, r, r, MPI_COMM_WORLD, MPI_STATUS_IGNORE );

				//MPI_Recv does not return until message arrived because unsent messages can result in deadlock in this blocking communication
				for ( unsigned int trg = 0; trg < ntargets; ++trg ) { //master clears buffer to avoid numerical data corruption
					localacc.at(trg).imcurv += rcvbuf[trg].imcurv;		rcvbuf[trg].imcurv = 0.0;
					localacc.at(trg).pcurv += rcvbuf[trg].pcurv;		rcvbuf[trg].pcurv = 0.0;
					localacc.at(trg).vcurv += rcvbuf[trg].vcurv;		rcvbuf[trg].vcurv = 0.0;
					localacc.at(trg).nvt += rcvbuf[trg].nvt;			rcvbuf[trg].nvt = 0;
					localacc.at(trg).nvv += rcvbuf[trg].nvv;			rcvbuf[trg].nvv = 0;
																		rcvbuf[trg].gid = 0;
					localacc.at(trg).n += rcvbuf[trg].padding;			rcvbuf[trg].padding = 0;
				}
			}
			else if ( this->get_Rank() == r ) {
				MPI_Send(sndbuf, ntargets, MPI_CapillaryDump_Type, MASTER, r, MPI_COMM_WORLD);
			}
			//else only MASTER and r involved in point-to-point communication
			MPI_Barrier(MPI_COMM_WORLD);
		}
	}

	if ( sndbuf != NULL ) { delete [] sndbuf; sndbuf = NULL; }
	if ( rcvbuf != NULL ) { delete [] rcvbuf; rcvbuf = NULL; }

	MPI_Barrier( MPI_COMM_WORLD );

	//ausgabe in textfile
	if ( this->get_Rank() == MASTER ) {
		string logfname;
		logfname = "TopologyTracer.CapTracking2D.SimID." + std::to_string(Settings::SimID) + ".Average.csv";
		ofstream logfile;
		logfile.open( logfname.c_str(), std::ofstream::out | std::ofstream::trunc );
		if (logfile.is_open() == true) {
			//header
			logfile << setprecision(18) << "GrainID;ICURV/n;PCURV/n;VCURV/n;nvt/n;nvv/n;n\n";

			for (unsigned int trg = 0; trg < ntargets; ++trg) {
				unsigned int gID = GrainIDWhiteList.at(trg);
				if ( localacc[trg].n > 0 ) { //only output existent grains for which sufficient capillary tracking data are available...
					double _cnts = 1.0 / ((double) localacc[trg].n);
					logfile << gID << ";" << localacc[trg].imcurv * _cnts << ";" << localacc[trg].pcurv * _cnts << ";" <<  localacc[trg].vcurv * _cnts;
					logfile << ";" <<  ((double) localacc[trg].nvt) * _cnts << ";" <<  ((double)  localacc[trg].nvv * _cnts) << ";" << localacc[trg].n << "\n";
				}
			}
			logfile.flush();
			logfile.close();
		}
		else {
			cerr << "ERROR::Worker::CapillaryTracking " << this->get_Rank() << " unable to write average results" << endl;
		}
	}
}


void curvfuseHdl::spit_profiling( void ) 
{
	//first overview of profiling for single process execution...
	string logfname;
	logfname = "TopologyTracer.CapTracking2D.SimID." + std::to_string(Settings::SimID) + ".Rank." + std::to_string(this->get_Rank()) + ".Profiling.csv";
	ofstream logfile;
	logfile.open( logfname.c_str(), std::ofstream::out | std::ofstream::trunc );
	if (logfile.is_open() == true) {
		//header
		logfile << "WhatTask;WallClock\n";
		logfile << ";s\n";
		logfile << ";MPI_Wtime\n";

		for (unsigned int l = 0; l < myprofiler.get_nentries(); l++) {
			logfile << myprofiler.titles[l] << ";" << setprecision(8) << myprofiler.times[l] << "\n";
		}
		logfile.flush();
		logfile.close();
	}
	else {
		cerr << "ERROR::Worker::CapillaryTracking " << this->get_Rank() << " unable to write profiling data" << endl;
	}
}


