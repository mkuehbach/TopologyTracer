//MK::TopologyTracer is a software developed by Markus K{\u}hbach in 2015
//for questions and comments please contact markus.kuehbach at rwth-aachen.de, last major revision 20150720
//##MK::marks comments that where the code is in rather quick and dirty,i.e. parts of the code which have been validated operative as is but offer certainly room for further improvements

#ifndef _TOPOLOGYTRACER_KERNEL_H_
#define	_TOPOLOGYTRACER_KERNEL_H_

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
//##DEBUG#include <math.h>
#include <iomanip>
#include <algorithm>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>


#include <mpi.h>

//#DEBUG#include "topologytracer_defs.h"
#include "topologytracer_math.h"

using namespace std;


//switch on for GNU, ##MK::prevents safe conversion to 32bit machine...
typedef long __int64;


struct ori
{
	double bunge1;			//Bunge (3,1,3) convention Euler angles
	double bunge2;
	double bunge3;
	double q0;				//unit quaternion representing the orientation
	double q1;
	double q2;
	double q3;
	uint32_t id;
	uint32_t closest;		//an ID to refer to this particular orientation
};
typedef struct ori * oriP;


struct adjointb2d
{
	uint32_t globalid;
	double length;
};
typedef adjointb2d * adjointb2dP;


typedef struct
{
	uint32_t fid;		//timestep
	uint32_t nfaces;	//the number of faces of this grain
	double area;		//the area of this grain
	double mobilityenergy; //mobility*energy weighted along segments along the perimeter
} MPI_IO_GrainInfo2D;


typedef struct
{
	uint32_t lid;
	uint32_t ngrainsintotal;
	double complexity;
} MPI_IO_MetadataInfo;


struct nbortopo
{
	uint32_t* n;
	uint32_t* ids; //[0,...,kshell_n[0]) the first layer, [kshell_n[0]...kshell_n[1]) second layer, and so forth...
	//double* disori; //disorientation angle among the grains
};
typedef nbortopo * nbortopoP;


struct grain2d
{
	uint32_t globalid;			//pieces of information at which timestep?
	uint32_t nfaces;			//how many faces does the grain have?
	uint32_t ori;				//an index referencing orientations from oripool
	bool boundarycontact;
	//##MK::three further unsigned chars or bools possible, or combine into bitfield

	double area;				//the current size of the grain
	double areachange;			//##MK::a simulation/experimental analysis dependent change of the area
	double perimeter;			//the length of the perimeter
	double totalboundaryenergy;	//how much energy is stored \Sum i 1 nfaces = gamma_i * length(i)

	adjointb2dP neighbors2d;	//the neighbors of the grain
	nbortopoP kshell;			//shells of local neighbors
};
typedef grain2d * grain2dP;


//3D grains
struct face
{
	uint32_t gi;				//ids of disjoint grains forming the face
	uint32_t gj;
	double facearea;
	//##MK::add further scalar properties

	uint32_t nbindingtplines;	//how many triple lines limiting the area?
	uint32_t* bindingtplines;	//the ids of all these triple lines
};
typedef face * faceP;


struct tpline
{
	uint32_t gi;				//the three grains meeting at a triple line
	uint32_t gj;
	uint32_t gk;
	uint32_t qp1;				//the limiting quadruple points of the triple line
	uint32_t qp2;
	double length;				//some measure of length of the triple line
};
typedef tpline * tplineP;


struct quadruple
{
	double x;					//the position of the quadruple point
	double y;
	double z;
	double someprop;			//##MK::a property value of the quadruple point
};
typedef quadruple * quadrupleP;


struct grain3d
{
	uint32_t timestep;			//pieces of information at which timestep?
	uint32_t ori;				//an index referencing orientation from oripool

	uint32_t nfaces;			//how many disjoint faces has the grains?
	uint32_t ntplines;			//how many disjoint triplelines are limiting the faces of the grain?

	uint32_t nqpoints;			//how many quadruple points does the grain link to?
	bool boundarycontact;		//is there boundary contact?

	//##MK::fill in padding if desired here

	double volume;
	faceP faces;
	tplineP tplines;
	quadrupleP qpoints;
};
typedef grain3d * grain3dP;



struct metadata
{
	uint32_t globalid;
	uint32_t localid;
	uint32_t ngrainsintotal;	//how many structural elements are expected?
	int whichnode;				//##MK::stays int, in the next years probably not an MPI universe of 2^31 processes...
	double complexity;			//how computationally expensive is working on this dataset expected? ##MK:: a function of the number of grains, their topology...
};
typedef metadata * metadataP;



class topoHdl : public mathMethods
{
public:
	topoHdl();
	~topoHdl();

	

//MPI-related intracomm identification
	int nRanks;
	int myRank;

	MPI_Datatype MPI_IO_GrainInfo2D_Type;
	MPI_Datatype MPI_IO_MetadataInfo_Type;

//job file bookkeeping
	uint32_t jobid;
	uint32_t N0;
	uint32_t firstloaded;					//MK::datasets are consecutively numbered starting from 0 + offset
	uint32_t offsetloaded;
	uint32_t lastloaded;
	double memlimitpernode;


//global bookkeeping
	std::vector<metadata> assignment;		//which process works on which datasets
	double complexitymodel_pergrain;
	double complexitymodel_a;
	double complexitymodel_t0;
	double complexitymodel_q;
	//##MK::constants in a simple load complexity model which assumes only the total number of grains decisive
	//evoluting in the form of N(t) = a*(iter - iter0)^q


//prototypes, make user-defined readers later private and public the toplevel function
	//##DEBUGvoid read_rawdata_cmiessen_distribute(int first, int offset, int last);
	void createMPIDataTypes( void );

//importer functions to read simulation data
	bool read_rawdata_cmiessen_read_2d( uint32_t N0, uint32_t first, uint32_t offset, uint32_t last);	//C. Miessen 2D GraGrLeS
	bool read_rawdata_lbarrales_metadata_2d( uint32_t N0, const char * ensembleInfo );
	bool read_rawdata_lbarrales_read_2d( const char* dfname , uint32_t ngrains);											//L.A. Barrales-Mora 2D Vertex code

	uint32_t getOripoolID( uint32_t globalid );
	uint32_t getGlobalID( uint32_t globalid, uint32_t whichmydataset);
	uint32_t getNumberOfFaces( uint32_t globalid, uint32_t whichmydataset);
	double getGrainVolume( uint32_t globalid, uint32_t whichmydataset);
	inline double getMobilityTimesEnergy( double disori );
	double getAverageMobilityTimesEnergyPosition( uint32_t mydatapos, uint32_t whichmydataset );
	double getAverageMobilityTimesEnergyGlobalID( uint32_t globalid, uint32_t whichmydataset );
	uint32_t getOriIDOfNeighbor( uint32_t globalid, uint32_t thedataset );
	uint32_t checkCandidatesRiso( uint32_t datasetid, uint32_t whichmydataset );

//currently implemented analysis routines
	void analyse_grain2d_sizeevolution(uint32_t first, uint32_t last);
	void analyse_grain2d_sizeevolution_nbormisori(uint32_t first, uint32_t last);

	uint32_t lastVitalitySign( uint32_t tgr );
	uint32_t* withMPIfindAllDissimilarNBors( uint32_t tgr, uint32_t lastfid, uint32_t& howmanyever );
	void analyse_grain2d_onegrain_vs_nearestnbors( uint32_t targetgr, uint32_t first, uint32_t last);
	int analyse_grain2d_construct_kshell(uint32_t kmax, uint32_t first, uint32_t last);
	void analyse_grain2d_kshellevolution(uint32_t kmax, uint32_t first, uint32_t last);

private: 
//prototypes
	void init( void );
	//##DEBUGvoid read_rawdata_cmiessen_finfo( const char* fname, metadata * finfo );
	bool read_rawdata_cmiessen_interpret_2d ( const char* fname, uint32_t fid ); // uint32_t* nread, double* complexity );//, struct metadata * finfo );

	void find_neighbors( uint32_t whichset, uint32_t whichgrain, nbortopoP target );
	uint32_t get_closest_idealori( double * quat );
	uint32_t addOrientation( double * bunge, uint32_t m, bool recategorization );		//returns either the closest match of already known ones or adds to pool and then returns the ID

	inline double complexitymodel( uint32_t t );
	inline double complexitymodel_linear ( uint32_t N0, uint32_t Nfin, uint32_t istart, uint32_t i, uint32_t iend );
	inline double complexitymodel_linear ( uint32_t ngr );

	uint32_t* cutTopoTail( uint32_t* arr, uint32_t CurrentSize, uint32_t WhichToCutFirst );
	uint32_t* getTopoMemory( uint32_t* arr, uint32_t oldsize, uint32_t wheretoplace, uint32_t& resized );

//jobrelated

//database
	std::vector<uint32_t> mydatasets;							//global ID of the dataset
	std::vector<uint32_t> mydatasize;							//how many grains in this dataset
	std::vector<grain2dP> mydata;								//the datasets as starting with a pointer to a container of grain structs, avoiding vector of vector of grains
	std::vector<ori> oripool;									//the orientations utilized
	std::vector<double*> gen_aboav_weaire_avfaces;				//generalized Aboav Weaire law, how many faces in shell-th shell of a center grain in class k
	std::vector<double*> gen_aboav_weaire_nnbors;				//how many neighbors in shell-th shell of a center grain in class k
	std::vector<double*> disori_environment;					//how many of the neighbors of a grain in class shell-th are low-angle misori?
	std::vector<MPI_IO_GrainInfo2D> targetgrain_evolution;		//a collection
	uint32_t kshellmax;
	uint32_t nneighborsmax;

	double myMemGuard;						//how much the node is currently loaded

	//flags and operation mode
	unsigned char dimensionality;
};
typedef class topoHdl * topoHdlP;


#endif	/* TOPOLOGYTRACER_KERNEL_H */
