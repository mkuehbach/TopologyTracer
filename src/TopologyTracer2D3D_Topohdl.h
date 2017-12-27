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


#ifndef __TOPOLOGYTRACER2D3D_TOPOHDL_H__
#define __TOPOLOGYTRACER2D3D_TOPOHDL_H__

#include "TopologyTracer2D3D_Math.h"



#ifndef TOPOLOGYTRACER2D3D_TOPOHDL_H_
#define TOPOLOGYTRACER2D3D_TOPOHDL_H_

#define SINGLE_REGION							1


//#define USE_POLYTRI

class timeLogger {

public:
	timeLogger();
	~timeLogger();

	void logev(const string title, double time);
	unsigned int get_nentries( void );
	vector<double> times;
	vector<string> titles;
};


class MemRegion;
class topoHdl;
class mathMethods;

class Snapshot
{
	//friend class topoHdl;
public:
	Snapshot( void );
	~Snapshot( void );
		
	bool readFacesForEachGrain( void );
	unsigned int baryxyz2mid( double bcx, double bcy, double bcz );
	bool planDataPartitioning( void );

	bool omp_initOneMemRegion( unsigned int x, unsigned y, unsigned z );
	
	inline unsigned int mcxyz2intID ( unsigned int mcx, unsigned int mcy, unsigned int mcz );
	bool initMemRegions( void );
	bool fcbuf2database( void );
	bool grbuf2database( void );
	void snhello( void );
	bool analyze_knn( unsigned int kkmax, std::ofstream &ofstr );
	inline double calc_quantile( std::vector<double>* thecdf, double p );
	void analyze_vol_quants( MPI_VOLSTATS* meta, double* res );
	unsigned int gsd_binning(double normsize);
	void analyze_gsd_histcnt( MPI_VOLSTATS* meta, double* res );

	std::vector<MemRegion*> mrg;
	std::vector<IDRangeBucket> LookupTable;

	SnapshotMetaData meta;
	topoHdl* mytopohdl;
	/*
	* grains are organized in grain ID range
	* so if a thread requests to find a grain by ID, it compares only a few grains 
	* in the corresponding ID range bucket rather than the entire list of grains in the snapshot
	* on success this lookup table provides the MemRegionID and the Position in this list
	* in effect, the total amount of ID comparisons to identify random grain IDs reduces 
	* from O(Settings::LargestGrainID/2) to O(Settings::MaxIDRange/2),
	* i.e. for a 10^7 grains to 10^2 i.e. 5 orders of magnitude on average
	*/
};
typedef class Snapshot * SnapshotP;


class MemRegion
{
public:
	MemRegion( void );
	~MemRegion( void );

	bool writeLookupTable( unsigned int extGID, unsigned int mr, unsigned int lid, unsigned int location);
	unsigned int queryMemRegionViaLookupTable( unsigned int extGID, unsigned int lid );

	void addFaceData( double size, unsigned int refid, unsigned int nborid, unsigned int reflid, unsigned int nborlid );
	unsigned int InWhichMemRegionIsTheGrain( unsigned int id, unsigned int lid );
	bool omp_storeFaces( unsigned int this_mreg );
	bool omp_storeGrainData( void );
	bool omp_allocGrainMemory( unsigned int ngr_per_mr );

	unsigned int* nborid_malloc( unsigned int size );
	void nborid_delete( unsigned int* nb, ksmetaP smeta );
	unsigned int* nborid_expand( unsigned int* nb_old, unsigned int size_old, unsigned int size_new );
	unsigned int* nborid_resize( unsigned int* nb_old, unsigned int size_old, unsigned int copylimit );
	bool nborid_add_noduplicates( unsigned int* nbid, unsigned last_to_compare, unsigned int gid );
	void find_gid( unsigned int ggid, unsigned int* wmr, unsigned int* wid );
	bool omp_knn( unsigned int kkmax );
	bool omp_knn_elimbnd( unsigned int kkmax, std::ofstream &outputfile );
	inline double disori2mobility( double theta );
	void omp_knn_kshell_perimeter( std::vector<unsigned int>& kids, ksmetaP meta, kshell_props* res, unsigned int nk );
	double omp_getArea( unsigned int gid );
	double omp_getBaryCenterX( unsigned int gid );
	double omp_getBaryCenterY( unsigned int gid );

	void memhello( void );
	void dummyCheckGrainBucketConsistency( void );
	bool checkFaceReassignmentConsistency( void );
		
	Snapshot* mys;						//back-reference snapshot the MemRegion belongs to
	unsigned int intID;					//the internal id which -th MemRegion of the Snapshot
	//unsigned int NumberOfGrains;		//size of the data container in this region only
	box Geometry;						//the space of the simulation domain which the MemRegion covers
	std::vector<Grain> GrainBucket;		//actual data of all grains in the container
};
typedef class MemRegion * MemRegionP;



class topoHdl : public mathMethods
{
	//friend class Snapshot;

public:
	topoHdl(void);
	~topoHdl(void);

	//setter
	void set_MPICommunication( int r, int nr );

	//getter
	int get_Rank( void );
	int get_nRanks( void );
	
	//administrative functions
	void hellow( void );
	void init_MPIDatatypes( void );

	bool determineWorkload( void );
	void commitWorkload( void );

	void modelspecific_rawdata_modification( std::vector<SnapshotMetaData>::iterator siii );
	bool read_mpiio_faces_binary_2d3d( std::vector<SnapshotMetaData>::iterator sii );
	bool read_mpiio_grains_binary_2d( std::vector<SnapshotMetaData>::iterator sii );
	bool read_mpiio_grains_binary_3d( std::vector<SnapshotMetaData>::iterator sii );
	bool read_targets_fromlist( void );
	bool readSnapshot( std::vector<SnapshotMetaData>::iterator si, bool ClearDatabase );
	bool initMemoryForOneSnapshot(  std::vector<SnapshotMetaData>::iterator sii, bool cleardb );
	bool nfaces_fcbuf2grbuf( std::vector<SnapshotMetaData>::iterator sii );
	bool readDatasets( void );
	bool readDatasets(unsigned int fd);
	inline unsigned int snapshot2rank( unsigned int extID );

	bool checkIfDatasetIsComplete( void );
	bool queryFileSizes( void );
	bool partitionWorkload( void );
	void writeWorkloadPartitioning( void );

	void dummycontrol( void ); //##MK::developer routine for checking global consistency across processes
	void dummy_ascii_faces_binary( unsigned int extid );
	void dummy_ascii_grains_binary_3d( unsigned int extid );

	unsigned int WhichOfMyDatasets( unsigned int sid );
	struct descr_stats getDescrStats( std::vector<unsigned int>* cand, unsigned int lsid, bool getarea );
	struct target_props getTargetProperties( unsigned int gid, unsigned int lsid, bool getarea, bool getnfaces, bool gethagbfrac, bool getmobdsee, bool getdsee, bool getnfdiff );
	struct rxstats calcRXFraction( unsigned int lsid, std::vector<unsigned int>* dgr, std::vector<unsigned int>* rxg );
	struct nucmodel_bh analyzeNuc_BaileyHirsch( unsigned int gid, unsigned int lsid );
	unsigned int modf_binning(double disori);
	unsigned int see_binning(double see);
	void calc_modf( unsigned int lsid, double* binnedmodf );
	void calc_seedf( unsigned int lsid, double* binnedseedf );
	void calc_localmaxsize(unsigned int lsid, double* localbuffer);
	void calc_localavdrvfrc(unsigned int lsid, seeav* localbuffer);

	std::vector<unsigned int>* analyze_elimbnd_one_mpiioself(std::string logfn_method, std::vector<unsigned int>* OnlyTheseTargets, unsigned int extID);
	std::vector<unsigned int>* analyze_elimbnd_all_self(std::string logfn_method, std::vector<unsigned int>* BlackList, bool UseGrainIDWhiteList);
	std::vector<unsigned int>* analyze_elimbnd_all_db( std::string logfn_method, std::vector<unsigned int>* BlackList, bool UseGrainIDWhiteList );
	std::vector<unsigned int>* analyze_elimbnd_one_db( std::string logfn_method, std::vector<unsigned int>* OnlyTheseTargets, unsigned int extID );
	
	void analyze_grainsize_quantiles_db( void );
	void analyze_see_db( void );
	void analyze_modf_db( void );
	void analyze_gsd_db( void );
	void analyze_drvforce_see_db( void );
	void analyze_maxsizegain_fw_db( void );
	void analyze_rxfraction_db(unsigned int extid);

	void analyze_trajectories_fw_db( void );
	void analyze_trajectories_bk_db( void );
	void analyze_abnormalgraingrowth_db( void );
	void analyze_sizegain_vs_matrix_bk_db( void );
	
	void analyze_knn_naive_db( unsigned int kmax );
	void analyze_classical_nucmodels_db( void );
	

	unsigned int IsIncluded(unsigned int ggid, std::vector<unsigned int>* cand);
	struct aggstats inspectForAbnormalGrains(unsigned int lsid, std::vector<unsigned int>* surv, std::vector<unsigned int>* matr, std::vector<unsigned int>* targ);
	unsigned int ngr_elimbnd( unsigned int extID );
	void gid_elimbnd( unsigned int extID, unsigned int* res );
	void spit_profiling(string mode);
	
	void analyze_grainsize_quantiles_adhoc(unsigned int fd);
	void analyze_see_adhoc(unsigned int fd);
	void analyze_modf_adhoc(unsigned int fd);
	void analyze_gsd_adhoc(unsigned int fd);
	void analyze_drvforce_see_adhoc(unsigned int fd);
	void analyze_maxsizegain_fw_adhoc(unsigned int fd);
	void analyze_meandrvforcesee_adhoc(unsigned int fd);
	void analyze_rxfraction_adhoc(unsigned int fd, std::vector<unsigned int>* matrix, std::vector<unsigned int>* rx );
	void analyze_trajectories_fw_adhoc(unsigned int fd, std::vector<unsigned int>* thetargets);
	void analyze_trajectories_bk_adhoc(unsigned int fd, std::vector<unsigned int>* survivors );
	void analyze_abnormalgraingrowth_adhoc(unsigned int fd, std::vector<unsigned int>* survivors, std::vector<unsigned int>* matrix, std::vector<unsigned int>* matrix0surv );
	void analyze_sizegain_vs_matrix_bk_adhoc(unsigned int fd, std::vector<unsigned int>* survivors, std::vector<unsigned int>* matrix );
	void analyze_topodiff_fw_adhoc(unsigned int fd, std::vector<unsigned int>* thetargets);
	//void analyze_knn_naive_adhoc(unsigned int fd);
	//void analyze_classical_nucmodels_adhoc(unsigned int fd);

	bool readSingleDataset( unsigned int fd );
	bool readSingleDataset( std::vector<SnapshotMetaData>::iterator si, bool loadgrains, bool loadfaces);
	bool analyze_rvebnd_contacttime( void );
	bool load_rvebnd_contacttime(void);
	bool master_mpiio_init_grainsize_quantiles(void);
	bool master_mpiio_init_see(void);
	bool master_mpiio_init_modf(void);
	bool master_mpiio_init_gsd(void);
	bool master_mpiio_init_drvfrc_see(void);
	bool master_mpiio_init_maxszgain_fw(void);
	bool master_mpiio_init_meandrvfrcsee_fw(void);
	bool master_mpiio_init_rxfrac(void);
	bool master_mpiio_init_fw(void);
	bool master_mpiio_init_bk(void);
	bool master_mpiio_init_agg(void);
	bool master_mpiio_init_szgain_bk(void);
	bool master_mpiio_init_topodiff_fw(void);
	//bool master_mpiio_init_knn(void);
	//bool master_mpiio_init_classicalnuc(void);
	bool init_results_container(void);
	bool workpartitioning_fid(unsigned int fid);

	void init_helper(void);
	void consolidate_helper(void);
	
	unsigned int grbufsize;									//how many elements currently in the buffer?
	MPI_3DGrainIO* grbuf;									//a process-local buffer for MPI_IO of Texture_* spread the data into the database such that numa localities can be respected

	unsigned int fcbufsize;									//how many faces currently held in fcbuf?
	MPI_3DFaceIO* fcbuf;									//a process-local buffer for MPI_IO similar as above but for Faces_*

	unsigned int bndbufsize;
	MPI_GrainBndContactIO* bndbuf;

	timeLogger myprofiler;
	timeLogger ioprofiler;

	//MK::helper have always size = 1 + Settings::SnapshotLast if allocated!
	double* helper_maxszgain_fw;
	seeav* helper_meandrvfsee_fw;


private:
	//Process-local database
	std::vector<SnapshotMetaData> DatasetInfo;				//the entirety of snapshots into which the dataset partitions
	std::vector<Snapshot*> LocalDB;
	std::vector<unsigned int> GrainIDWhiteList;				//holding IDs from TargetGrainIDs to implement E_FORWARD_SELECTED during forward tracking

	MPI_Datatype MPI_GrainBndContactIO_Type;
	MPI_Datatype MPI_2DGrainIO_Type;
	MPI_Datatype MPI_3DGrainIO_Type;						//physics data
	MPI_Datatype MPI_3DFaceIO_Type;
	MPI_Datatype MPI_SMDIO_Type;								//management data
	MPI_Datatype MPI_SEEAV_Type;
	MPI_Datatype MPI_VOLSTATS_Type;
	MPI_Datatype MPI_RXSTATS_Type;
	MPI_Datatype MPI_QUANTILES_Type;
	//MPI_Datatype MPI_NUCBAIHIR_Type;
	MPI_Datatype MPI_AGGSTATS_Type;



//MPI-related intracomm identification
	int nRanks;
	int myRank;
};
typedef class topoHdl * topoHdlP;


class spatQueryHdl
{
public:
	spatQueryHdl(void);
	~spatQueryHdl(void);

	//setter and getter
	void set_MPICommunication(int r, int nr) {
		myRank = r; nRanks = nr;
	}
	int get_Rank(void) { return myRank; }
	int get_nRanks(void) { return nRanks; }

	inline unsigned int pos2bucket(double x, double y, double z); //map normalized position on unit cube to container

	bool init_spatquerybox(void);
	bool read_microstructure_uds(void);
	bool read_longrange( unsigned int fd );
	bool read_maxgrainsizes(string fn);
	bool init_fillspatquerybox(void);
	std::vector<sqkey>* rangesearch( unsigned int gID, point3d p, double r);

	void init_neighbordistances(void);
	double quantile_asc_sorted_dbl( std::vector<double>* dat, double qnt );
	void arrivaltime_quantiles( unsigned int fd );


	//query fast all points up to R distance from given position in rve
	std::vector<std::vector<mpoint3d>*> pp3;
	std::vector<std::vector<sqkey>*> nbors;
	sqb rve;

	//array of size 1+Settings::LargestGrainID
	point3d* ppp; //gid2barycenter
	double* chi; //gid2chi
	double* mxsz; //gid2maxsz
	bool* process;	 //gid2bnd

private:
	int myRank;
	int nRanks;
};
typedef class spatQueryHdl * spatQueryHdlP;

class utilityHdl
{
public:
	utilityHdl(void);
	~utilityHdl(void);

	//setter
	void set_MPICommunication( int r, int nr );

	//getter
	int get_Rank( void );
	int get_nRanks( void );

	bool init_KNNDatabase( void );
	bool read_targetfile( void );
	bool read_knnfile( unsigned int kkmax );

	void analyze_apriori_prediction( void );
	void analyze_find_knn_for_targets( void );
	void analyze_Gest_for_targets( void );
	
	std::vector<unsigned int> TargetGIDs;
	std::vector<std::vector<knn_meta>*> KNNDatabase; //thereby generating a hashtable of IDs in the range [0, Settings::LargestGrainID]
	unsigned int kkmax;

private:

//MPI-related intracomm identification
	int nRanks;
	int myRank;
};
typedef class utilityHdl * utilityHdlP;


//#ifdef USE_POLYTRI
class arealanaHdl : public mathMethods
{
public:
	arealanaHdl(void);
	~arealanaHdl(void);

	//setter
	void set_MPICommunication( int r, int nr );
	void init_MPIDatatypes( void );
	
	//getter
	int get_Rank( void );
	int get_nRanks( void );

	void init( void );
	bool read_microstructure_uds_2d( void );
	void read_initialmicrostructure_metadata( void );
	bool read_snapshot_2d( unsigned int fd );


	struct aabb2d calc_polygon_aabb( std::vector<point2d>* thecontour, unsigned int gID );
	double calc_polygon_area( std::vector<point2d>* thecontour );
	void calc_contourarea( void );
	void reverse_contours( void );
	void remove_core_duplicates( void );
	void remove_degenerated_grains( void );
	void remove_lastpoints( void );
	void remove_insufficiently_supported_contours( void );
	void modify_rawdata_gragles( unsigned int fd );
	std::vector<unsigned int>* triangularize_contours( unsigned int fd );
	void calc_barycenter( void );
	void calc_aabb( void );

	struct overlap2d calc_overlap_polygon_circle( unsigned int gidx, double xcir, double ycir, double r, bool render );
	void eliminate_non_triangularized_grains( unsigned int fd, std::vector<unsigned int>* gid_to_kick );
	void envcharacterize_omp_naive( unsigned int fd ); //off the shelf OMP parallelization without considering NUMA issues...
	void report( unsigned int fd );
	void reset_contour_memory( unsigned int fd );
	void metareport( void );
	void spit_profiling(); //MK::first dev solution

	std::vector<isectionmeta> gmetadata;
	//each arealanaHdl has only one storage structure to keep in memory the
	//currently worked on snapshot, otherwise the amount of contour data very likely flood memory
	std::vector<double> gcontour2darea;					//storing the area of a specific contour
	std::vector<std::vector<point2d>*> gcontour2d;		//non-NULL pointer to contour points of remaining grains with grainIDs on interval [1,Settings::LargestGrainID]
	std::vector<std::vector<triangle>*> gsegments2d;	//individual triangles building the 2d contour
	std::vector<barycenter2d> contours_bary_quick;		//shorter array eluding non-existent IDs
	std::vector<aabb2d> contours_aabb_quick;
	std::vector<std::vector<chi>*> ChiTable;

	struct analysismetrics2d fdmeta;

private:

//MPI-related intracomm identification
	int nRanks;
	int myRank;

//MPI-related types
	MPI_Datatype MPI_GBContourPoint_Type;
	//MPI_Datatype MPI_GNUInfo_Type;

//Internal profiling
	timeLogger myprofiler;
};
typedef class arealanaHdl * arealanaHdlP;
//#endif



//currently discrete approximation therefore tetgen not defined
class volanaHdl : public mathMethods
{
public:
	volanaHdl(void);
	~volanaHdl(void);

	//setter
	void set_MPICommunication( int r, int nr );
	void init_MPIDatatypes( void );
	
	//getter
	int get_Rank( void );
	int get_nRanks( void );

	bool read_snapshot_volume_mpi_gragles( std::string fn, unsigned int extx, unsigned int exty, unsigned int extz );
	bool read_snapshot_volume( unsigned int fd, unsigned int size );
	bool read_snapshot_meta_txt_gragles( std::string fn );
	bool read_snapshot_meta( unsigned int fd );
	bool read_snapshot_3d( unsigned int fd );
	void calc_grainvolume( unsigned int fd );
	unsigned int read_snapshot_rve_extent( unsigned int fd );
	void envcharacterize_omp_naive( snapshotmeta thisone );
	void write_envsphere( unsigned int refgid, unsigned int fid, int bx, int by, int bz, int rd, int xmin, int xmax, int ymin, int ymax, int zmin, int zmax, int size );
	void init_discrete_version( void );

	bool read_microstructure_uds_3d( void );
	void read_initialmicrostructure_metadata( void );
	void report( snapshotmeta thisone );
	void metareport( void );
	void reset_xitable( void );
	void spit_profiling();
	
	std::vector<unsigned int> ghullvolume;
	std::vector<isectionmeta> gmetadata;					//##MK::consider in the future not to set up as a hash but rather a map or tree as many hash places remain empty for timesteps with few grains remaining
	std::vector<std::vector<unsigned int>*> AssgnSnapshots;	//collection of pointer to local data/fid's a rank processes pointer point to implicit 3d array of grain IDs telling which grain at which position
	std::vector<snapshotmeta> AssgnSnapshotsMeta;			//AssgnSnapshot[i] dataset with ID AssgnSnapshotsFID[i]
	std::vector<std::vector<UDS_Grain3DInfo>*> MetaSnapshots; //the corresponding updated! metadata
	std::vector<std::vector<xi>*> XiTable;

private:

//MPI-related intracomm identification
	int nRanks;
	int myRank;

//MPI-related types

//Internal profiling
	timeLogger myprofiler;
};
typedef class volanaHdl * volanaHdlP;


class curvapprxHdl : public mathMethods
{
public:
	curvapprxHdl(void);
	~curvapprxHdl(void);

	//setter
	void set_MPICommunication( int r, int nr );
	void init_MPIDatatypes( void );

	//getter
	int get_Rank( void );
	int get_nRanks( void );

	bool init( void );
	bool read_snapshot_contour_mpi_gragles( std::string fn );
	bool read_snapshot_contour( std::string prefix, unsigned int fid );
	void write_snapshot_contour_gnu_gragles( unsigned int fd );
	void modelspecific_rawdata_modification( void );
	struct wghtd_imcurv_res ompcrit_imean_curv_core( unsigned int g );
	bool integral_mean_curvature_approximation( unsigned int fd );
	void report_wghtd_imcurv( std::string prefix, unsigned int fd, bool clear_rbucket );
	void reset_memory( void );

	std::vector<std::vector<GBContourPoint>*> gcontour2d; //hash of piecewise linearized contour segments for threads to operate on
	std::vector<std::vector<wghtd_imcurv_res>*> gresults_data2d; //bucket pointing to individual result snippets
	std::vector<imcurv_meta> gresults_meta2d;

private:

//MPI-related intracomm identification
	int nRanks;
	int myRank;

//MPI-related types
	MPI_Datatype MPI_GBContourPoint_Type;
	MPI_Datatype MPI_GBJunctionPoint_Type;

//Internal profiling
	timeLogger myprofiler;
};
typedef class curvapprxHdl * curvapprxHdlP;



class curvfuseHdl : public mathMethods
{
public:
	curvfuseHdl(void);
	~curvfuseHdl(void);

	//setter
	void set_MPICommunication( int r, int nr );
	void init_MPIDatatypes( void );

	//getter
	int get_Rank( void );
	int get_nRanks( void );

	bool init_targets( void );
	bool read_targets_fromlist( void );
	//unsigned int check_dumpdata_existence();
	void reset_dumpdata_todefaults( void );
	unsigned int init_localbuffer();
	unsigned int init_results_file();
	bool readDatasets(unsigned int fid );
	bool workpartitioning_fid( unsigned int fid);
	void fuse_capillary_activity(unsigned int fid);
	void localaveraging_capillary_activity(unsigned int fid);
	void spit_average_results(void);
	void spit_profiling( void );

	std::vector<unsigned int> GrainIDWhiteList;				//holding IDs from TargetGrainIDs to work on
	std::vector<capdump> dumpbuffer;

	std::vector<capav> localacc;				//local accumulation of averaging


	timeLogger myprofiler;

private:

//MPI-related intracomm identification
	int nRanks;
	int myRank;

//MPI-related types
	MPI_Datatype MPI_CapillaryDump_Type;
};
typedef class curvfuseHdl * curvfuseHdlP;


#endif

#endif