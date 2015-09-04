//MK::TopologyTracer is a software developed by Markus K{\u}hbach in 2015
//for questions and comments please contact markus.kuehbach at rwth-aachen.de, last major revision 20150720
//##MK::marks comments that where the code is in rather quick and dirty,i.e. parts of the code which have been validated operative as is but offer certainly room for further improvements


#ifndef TOPOLOGYTRACER_DEFS_H
#define	TOPOLOGYTRACER_DEFS_H

//MPI related
#define MASTER									0
#define SEQUENTIAL								1
#define MINIMUM_COMPLEXITY_PERNODE				((10.0)*(1024.0)*(1024.0))		//in byte
#define MAXIMUM_COMPLEXITY_PERNODE				((4096.0)*(1024.0)*(1024.0))	//in byte

//user-defined dataset related
#define CMIESSEN_ORIFILE_MAXIMUMNEIGHBORS		100
#define CMIESSEN_ORIFILE_BUFFERLENGTH			2500
#define CMIESSEN_ORIFILE_GRAINPRECACHE			1000

#define NOTASSIGNEDYET							-1
#define NOTCHECKED								false
#define NOTANALYZED								(0.0)
#define AVERAGE_NEIGHBORS						6						//Mullins 2D
#define COMPLEXITY_PER_NEIGHBOR					50						//byte
#define	IMPOSSIBLE_NFACES_AS_DUMMY				0						//there is not any grain with such number of faces

#define DEFAULT_KSHELL							4
#define	DEFAULT_NNEIGHBORS						25
#define DEFAULT_EXPECTATION_DIFF_NBORS			1000


//definition of ideal texture component / ideal / standardlagen
#define RANDOM_ORIENTATION						0
#define MAX_DISORI_FCC							(1.099)						//MacKenzie m-3m, 1..., 62.8/180*_PI_
#define RESOLUTION_SO3GRID						(0.00174532925199433)		//0.1 degree raster
#define DISORICACHE_MEMLIMIT					(200.0)						//maximum size of lower-triangle matrix implemented as implicitly addressed array of types double
#define MAXDISORI_TO_40DEG111					(0.174532925199433)			//10/180 * _PI_;
#define MAXDISORI_LAGB2HAGB						(0.261799387799149)			//15/180 * _PI_;
#define NOT_WITHIN_GRIDRESU						-1


//complexity modeling
#define TYPICALSIZE_2D_GRAIN					112
#define TYPICALSIZE_3D_GRAIN					1106


//operation modes
#define ANALYZE_IN_2D							0x02
#define ANALYZE_IN_3D							0x03

#define LOADING_DATASET_FAILED					false
#define LOADING_DATASET_SUCCESS					true


//numerical bounds
#define UINT32T_MAX								4294967294 //2^32-1
#define MAXIMUM_NUMBER_OFGRAINS					4096000000 //I would like to see this simulation exceeding it, just let me know ...
#define NO_PROPERTY_ASSIGNED					4294967295
#define DATA_VALID								0x02
#define DATA_INVALID							0x01
#define TARGETNODE_TAKES						0x02
#define TARGETNODE_REJECTS						0x01
#define UNKNOWN									4294967294
#define NEVER_OBSERVED							0
#define NO_NEIGHBORS_EVER						0

//decisions
#define DO_NOT_INCLUDE							4294967294 //2^32-2

#endif

