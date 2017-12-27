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

#ifndef __TOPOLOGYTRACER2D3D_DATATYPES_H__
#define __TOPOLOGYTRACER2D3D_DATATYPES_H__

#include "TopologyTracer2D3D_Settings.h"
//#include "TopologyTracer2D3D_Defs.h"

/*
	MK::The key idea of the database is as follows: a simulation with individual timesteps fid, the snapshots,	
	are partitioned onto the MPI processes (ranks). Each of which holds a number of snapshots for a contiguous interval 
	in time. Each snapshot is partitioned further into a number of so-called local memory regions that host a 
	certain contiguous volume of the simulation domain. The purpose of this design is to enable ccNUMA-aware spatial 
	locality for the processing of each snapshot with via for instance OpenMP threading.
*/

struct descr_stats					//generic output structure for descriptive statistics
{
	double sum;
	double mean;
	unsigned int n;
	//4B padding left
	descr_stats(double _s, double _m, unsigned int _n) : sum(_s), mean(_m), n(_n) {}
	descr_stats() : sum(0.0), mean(0.0), n(0) {}
};
typedef descr_stats * descr_statsP;


struct target_props
{
	double area;					//area of a grain in 2d, volume in 3d
	double hagbfraction;			//boundary length (2d), face area fraction (3d) with HAGB
	double mob_dsee;				//mobility times stored elastic energy (SEE)
	double dsee;					//stored elastic energy difference over perimeter
	unsigned int nfacesk1;			//first order number neighbors of the grain
	int nfacesk1diff;				//sum of difference in faces to nearest neighbors
	//4B padding left
	target_props(double _sz, double _hagbf, double _mse, double _dsee, unsigned int _nf1, int _nfd1) : area(_sz), 
		hagbfraction(_hagbf), mob_dsee(_mse), dsee(_dsee), nfacesk1(_nf1), nfacesk1diff(_nfd1) {}
	target_props() : area(0.0), hagbfraction(0.0), mob_dsee(0.0), dsee(0.0), nfacesk1(0), nfacesk1diff(0) {}
};
typedef target_props* target_propsP;


struct rxstats
{
	double totalsize_elimbnd;					//domain area covered by grains without boundary contact
	double totalsize_targets;					//domain area covered by the target grains
	double X;									//X := totalsize_targets/totalsize_elimbnd;
	double ndefg;								//how many matrix grains
	double nrxg;								//how many considered as to be rxed grains
	rxstats(double _tsze, double _tszt, double _X, double _ndg, double _nrg) : totalsize_elimbnd(_tsze),
		totalsize_targets(_tszt), X(_X), ndefg(_ndg), nrxg(_nrg) {}
	rxstats() : totalsize_elimbnd(0.0), totalsize_targets(0.0), X(0.0), ndefg(0), nrxg(0) {}
};
typedef rxstats* rxstatsP;


struct nucmodel_bh
{
	//instantaneous analysis of a target on 1-st order neighbors
	double hagbfraction;						//HAGB fraction of this grain

	//actual determined EHR, if grain does not have 50% HAGB for instance the r_hagb50_maxdrho = BH_RADIUS_INVALID;
	double r;									//EHR size (area, vol) of the target, set irrespective whether grain has HAGB
	double r_hagb25;							//EHR size (area, vol) of the target, only different from 0.0  when hagbfraction >= 25%
	double r_hagb50;							//... only different from 0.0 when hagbfraction >= 50%
	double r_hagb75;							//... only different from 0.0 when hagbfraction >= 75%

	double hv_hagb_drho;						//evaluating SE only about perimeter with HAGB character
	double hv_all_drho;							//evaluating stored energy difference about the perimeter with Heaviside function, i.e. only positive driving forces count (eliminate possibility for negative values in the rc terms)
	//hence it may often appear that for a grain with only negative net rho differences over the perimeter the result is zero

	//critical values according to the criterion
	//SE model Bailey Hirsch
	double rc_baihir_hagb;						//critical size for this grain considering as the driving force difference the evaluate of hv_hagb_drho
	double rc_baihir_all;						//critical size for this grain considering as the driving force difference the evaluate of hv_all_drho
	//SE model Bate Hutchinson
	double rc_bathut_hagb;						// ... see above
	double rc_bathut_all;						// ...

	bool grainfound;							//allows to avoid criterion comparisons for nuclei that were not found
	//7B padding left
	nucmodel_bh() : hagbfraction(0.0), r(0.0), r_hagb25(0.0), r_hagb50(0.0), r_hagb75(0.0), hv_hagb_drho(0.0), hv_all_drho(0.0), rc_baihir_hagb(INFINITE_RADIUS), rc_baihir_all(INFINITE_RADIUS), rc_bathut_hagb(INFINITE_RADIUS), rc_bathut_all(INFINITE_RADIUS), grainfound(false) {}
};
typedef nucmodel_bh* nucmodel_bhP;


struct aggstats
{
	double totalsize_survivors;		//area/volume covered by survivors
	double totalsize_matrix;		//..the matrix considered, i.e. excluding survivors
	double totalsize_targets;		//..the targets

	double rmean_survivors;			// EHR of survivors
	double rmean_matrix;			// EHR of matrix grains
	double rmean_targets;			// EHR of target grains, there will be a difference to the rmean_matrix because the targets account also for the abnormal grains and hence bias the mean value

	//counts also in double to simplify reading in with MATLAB
	unsigned int n_survivors;		// how many grains still there consistency
	unsigned int n_matrix;			// will decrease strong than from global coarsening because population narrows down to survivors only
	unsigned int n_targets;			// bound to n_targets > n_matrix
	unsigned int n_agg_m;			// with respect to mean of matrix
	unsigned int n_agg_t;			// with respect to mean of targets
	unsigned int n_agg_m_insurv;	// ... but also included in survivors
	unsigned int n_agg_t_insurv;	// ... but also included
	//4B padding left
	aggstats( double _tsz_s, double _tsz_m, double _tsz_t, double _ehr_s, double _ehr_m, double _ehr_t, 
		unsigned int _n_s, unsigned int _n_m, unsigned int _n_t, unsigned int _na_m, unsigned int _na_t, 
		unsigned int _na_mi, unsigned int _na_ti) : totalsize_survivors(_tsz_s), totalsize_matrix(_tsz_m),
		totalsize_targets(_tsz_t), rmean_survivors(_ehr_s), rmean_matrix(_ehr_m), rmean_targets(_ehr_t),
		n_survivors(_n_s), n_matrix(_n_m), n_targets(_n_t), n_agg_m(_na_m), n_agg_t(_na_t), n_agg_m_insurv(_na_mi),
		n_agg_t_insurv(_na_ti) {}
	aggstats() : totalsize_survivors(0.0), totalsize_matrix(0.0), totalsize_targets(0.0), rmean_survivors(0.0), rmean_targets(0.0), n_survivors(0), n_matrix(0), n_targets(0), n_agg_m(0), n_agg_t(0), n_agg_m_insurv(0), n_agg_t_insurv(0) {}
};
typedef aggstats* aggstatsP;


struct kshell_props
{
	double ksbnd_peri;			//region boundary perimeter or area
	double ref_phi1;			//orientation of the reference
	double ref_PHI;
	double ref_phi2;
	double ref_see;				//stored elastic energy of the reference

	double area_hagbthresholded;//accounts for the perimeter length with only a HAGB
	double area_velocity_w;		//accounts for assuming pure stored elastic energy driven only velocity perimeter how beneficial is the environment drive expansion of the grain
	double area_ksgrains;		//the total area of all grains in this shell
	//i.e. if a four-sided grain has two opposite boundaries same mobility (lets assume 1.0) and same dSEE positive, these faces will cause the grain to expand
	//however if the two others are of same mobility but negative dSEE but abs(dSEE) and only if the geometry is the same the grain will maintain its shape
	//double area_disori2ref_w;	//weighted by disorientation of the neighbors to a reference (the target grain in the 0-th shell) over the segments
	//double area_dsee2ref_w;		//weighted by difference in stored elastic energy  of the neighbors to a reference (the target grain in the 0-th shell) over the segments
	//kshell_props() : area(INVALID_LENGTH), ref_phi1(0.0), ref_PHI(0.0), ref_phi2(0.0), ref_see(0.0), area_hagbthresholded(0.0), area_disori2ref_w(0.0), area_dsee2ref_w(0.0) {}
	kshell_props(double _ks, double _e1, double _e2, double _e3, double _see, double _ah, double _avw, double _akg ) :
		ksbnd_peri(_ks), ref_phi1(_e1), ref_PHI(_e2), ref_phi2(_e3), ref_see(_see), area_hagbthresholded(_ah),
			area_velocity_w(_avw), area_ksgrains(_akg) {}
	kshell_props() : ksbnd_peri(INVALID_LENGTH), ref_phi1(0.0), ref_PHI(0.0), ref_phi2(0.0), ref_see(0.0), area_hagbthresholded(0.0), area_velocity_w(0.0), area_ksgrains(0.0) {}
};
typedef kshell_props* kshell_propsP;


struct knn_meta
{
	unsigned int gid;
	unsigned int pad;
	double x;
	double y;

	unsigned int kshellbased;	//data arrays with elements [0, kkmax]
	unsigned int peribased;		//data arrays with elements [0, kkmax) MK::mind exclusive!

	unsigned int* GrainsInShell;
	double* AreaOfKShell;
	double* BndLen;
	double* HAGB;
	double* MdSEE;
	knn_meta(unsigned int _gid, double _x, double _y, unsigned int _ksbsd, unsigned int _pebsd ) :
		gid(_gid), pad(0), x(_x), y(_y), kshellbased(_ksbsd), peribased(_pebsd), GrainsInShell(NULL), AreaOfKShell(NULL),
			BndLen(NULL), HAGB(NULL), MdSEE(NULL) {}
	knn_meta() : gid(THE_DOMAIN_ITSELF), pad(0), x(0.0), y(0.0), kshellbased(0), peribased(0), GrainsInShell(NULL), 
		AreaOfKShell(NULL), BndLen(NULL), HAGB(NULL), MdSEE(NULL) {}
};
typedef knn_meta * knn_metaP;


struct lukey
{
	unsigned int wmr;
	unsigned int wid;
	lukey(unsigned int _mr, unsigned _id) : wmr(_mr), wid(_id) {}
	lukey() : wmr(NOT_ASSIGNED_YET), wid(NOT_ASSIGNED_YET) {}
};
typedef lukey* lukeyP;


//metadata for k-nearest neighbor construction
struct ksmeta
{
	unsigned int first;							//enables to mark a range of consecutive elements in an array of grain IDs which are members of a certain k-shell
	unsigned int last;							//inclusive! [first, last]
	ksmeta(unsigned int _f, unsigned int _l ) : first(_f), last(_l) {}
	ksmeta() : first(0), last(0) {}
};
typedef ksmeta * ksmetaP;


struct box										//3d axis-aligned bounding box
{
	double xmi;
	double xmx;
	double ymi;
	double ymx;
	double zmi;
	double zmx;
	box(double _xmi, double _xmx, double _ymi, double _ymx, double _zmi, double _zmx ) : 
		xmi(_xmi), xmx(_xmx), ymi(_ymi), ymx(_ymx), zmi(_zmi), zmx(_zmx) {}
	box(): xmi(0.0), xmx(0.0), ymi(0.0), ymx(0.0), zmi(0.0), zmx(0.0) {}
};
typedef box * boxP;


//database entries
struct Face
{
	double size;								//length in 2D, area in 3D
	unsigned int extNBID;						//which neighbor is it?
												//indices aiding the quick accessing of the neighbor
	unsigned int mrg_idx;						//identifying the logical memory region of the snapshot
	unsigned int idx;							//identifying the position in the sub-list inside the memory region
	//4B padding
	Face( double _sz, unsigned int _extnbid, unsigned int _mrgidx, unsigned int _idx) : size(_sz), extNBID(_extnbid),
		mrg_idx(_mrgidx), idx(_idx) {}
	Face() : size(0.0), extNBID(THE_DOMAIN_ITSELF), mrg_idx(INVALID_MEMORYREGION), idx(INVALID_INDEX) {}
};
typedef Face * FaceP;


struct Grain
{
	unsigned int extGID;						//unique id \in \mathbb{N} and > 0 which the simulation assigned to the grain
	unsigned int nfaces;						//how many faces does the grain have
	unsigned int nfaces_identified;				//capture potential deviations in counting of faces
	//MK::depending on at which point the GraGLeS solver outputs the number of faces there may be faces bounded only by a single contour point, i.e. internal solver length scale h...
	unsigned int boundary;						//1 yes (grain has contact with the boundary), 0 no (inside bulk)

	double size;								//area in 2D, volume in 3D
	double surface;								//perimeter in 2D, surface area in 3D

	double phi1;								//grain orientation Bunge, "ZXZ" Euler angle in entire Euler space
	double Phi;
	double phi2;

	double x;									//barycenter of the grain
	double y;
	double z;									//is 0.0 for a 2D simulation

	//further properties
	double see;									//stored elastic energy
		
	Face* neighbors;							//a list of all neighboring grains (as such the faces)

	//##MK::add functionality for k-shell and junctions
	//nbortopoP kshell;			//shells of local neighbors
	Grain(unsigned int _extgid, unsigned int _nf, unsigned int _nfv, unsigned int _boundary, double _sz, double _surf,
		double _e1, double _e2, double _e3, double _x, double _y, double _z, double _see ) :
	extGID(_extgid), nfaces(_nf), nfaces_identified(_nfv), boundary(_boundary), size(_sz), surface(_surf), 
		phi1(_e1), Phi(_e2), phi2(_e3), x(_x), y(_y), z(_z), see(_see), neighbors(NULL) {}
	Grain() : extGID(THE_DOMAIN_ITSELF), nfaces(0), nfaces_identified(0), boundary(BOUNDARY_CONTACT), size(0.0), surface(0.0),
		phi1(0.0), Phi(0.0), phi2(0.0), x(0.0), y(0.0), z(0.0), see(0.0), neighbors(NULL) {}
};
typedef Grain * GrainP;

//##tripleline
//##quadruple point


struct SnapshotMetaData
{
	unsigned int extID;							//an ID the simulation assigned to the file of which this snapshot carries the data, it is external in the sense that someone else decided how the file should be named

	unsigned int nMemRegionsXYZ;
	unsigned int nMemRegionsX;
	unsigned int nMemRegionsY;
	unsigned int nMemRegionsZ;					//into how many logical sub-sets are the data for the snapshot split to assure spatial locality on ccNUMA nodes

	unsigned int ng;							//how many grains in this snapshot?
	unsigned int nf;							//how many faces in this snapshot?

	//internal organization of the snapshots on the processes
	unsigned int RID;							//the rank in which the snapshot is stored
	double size_ng;								//file size in byte for Texture_* bin contain grains and
	double size_nf;								//Faces_*.bin containing faces
	SnapshotMetaData( unsigned int _extid, unsigned int _nxyz, unsigned int _nx, unsigned int _ny, unsigned int _nz,
		unsigned int _ng, unsigned int _nf, unsigned int _rid, double _szng, double _sznf ) :
		extID(_extid), nMemRegionsXYZ(_nxyz), nMemRegionsX(_nx), nMemRegionsY(_ny), nMemRegionsZ(_nz),
			ng(_ng), nf(_nf), RID(_rid), size_ng(_szng), size_nf(_sznf) {}
	SnapshotMetaData() : extID(UNKNOWN_ID), nMemRegionsXYZ(1), nMemRegionsX(1), nMemRegionsY(1), nMemRegionsZ(1), 
			ng(0), nf(0), RID(UNKNOWN_ID), size_ng(0.0), size_nf(0.0) {}
};
typedef SnapshotMetaData * SnapshotMetaDataP;


struct GrainHash
{
	unsigned int gid;							//a grain ID, was externally defined during dataset creation
	unsigned int mrg_idx;						//the MemRegion to which the grain was assigned
	unsigned int idx;							//the position in the IDRange bucket in which the grain is located
	GrainHash( unsigned int _gid, unsigned int _mrgidx, unsigned int _idx) : 
	gid(_gid), mrg_idx(_mrgidx), idx(_idx) {}
	GrainHash() : gid(THE_DOMAIN_ITSELF), mrg_idx(INVALID_MEMORYREGION), idx(INVALID_INDEX) {}
};
typedef GrainHash * GrainHashP;


struct IDRangeBucket
{
	unsigned int next;							//next location that is filled in that will be filled in
	unsigned int size;							//identifies how many hash entries were found for these ID values

	GrainHash* bucket;							//points to a container which the threads can scan to identify a grain via if ( gid == bucket[cand].gid) found return
	IDRangeBucket( unsigned int _nx, unsigned int _sz, GrainHash* _bu ) : next(_nx), size(_sz), bucket(_bu) {}
	IDRangeBucket() : next(0), size(0), bucket(NULL) {}
};
typedef IDRangeBucket * IDRangeBucketP;



struct isectionmeta
{
	//MK::better structures will be necessary!
	unsigned int gid;
	unsigned int x0; //initial barycenter from synthetic microstructure generator, MK::necessary to avoid numerical drift of circular inspection area
	unsigned int y0;
	unsigned int z0;
	double q0;
	double q1;
	double q2;
	double q3;
	double see;
	isectionmeta(double _gid, unsigned int _x, unsigned int _y, unsigned int _z, double _q0, double _q1, double _q2, double _q3, double _see) : 
	gid(_gid), x0(_x), y0(_y), z0(_z), q0(_q0), q1(_q1), q2(_q2), q3(_q3), see(_see) {}
	isectionmeta() : gid(THE_DOMAIN_ITSELF), x0(0), y0(0), z0(0), q0(1.0), q1(0.0), q2(0.0), q3(0.0), see(0.0) {}
};
typedef isectionmeta * isectionmetaP;


struct snapshotmeta
{
	unsigned int fid;
	unsigned int size;
	snapshotmeta( unsigned int _f, unsigned int _s) : fid(_f), size(_s) {}
	snapshotmeta() : fid(0), size(0) {}
};
typedef snapshotmeta * snapshotmetaP;

struct point2d
{
	double x;
	double y;
	point2d(double _x, double _y) : x(_x), y(_y) {}
	point2d() : x(-1.0), y(-1.0) {} //initial assignment to check invalid
	//as we work only on positive quadrant
};
typedef point2d * point2dP;


struct point3d
{
	double x;
	double y;
	double z;
	point3d(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}
	point3d() : x(-1.0), y(-1.0), z(-1.0) {} //initial assignment to check invalid
	//as we work only on positive quadrant
};
typedef point3d * point3dP;


struct mpoint3d
{
	double x;
	double y;
	double z;
	unsigned int gid;
	bool exclude; //(==true) because of boundary contact?
	mpoint3d(double _x, double _y, double _z, unsigned int _gid, bool _exc) :
		x(_x), y(_y), z(_z), gid(_gid), exclude(_exc) {}
	mpoint3d() : 
		x(-1.0), y(-1.0), z(-1.0), gid(THE_DOMAIN_ITSELF), exclude(false) {} //initial assignment to check invalid
	//as we work only on positive quadrant
};
typedef mpoint3d * mpoint3dP;


struct sqkey
{
	double dist;
	unsigned int pid;
	bool exclude;
	bool pad1;
	bool pad2;
	bool pad3;

	sqkey(double _d, unsigned int _pid, bool _exc) : 
		dist(_d), pid(_pid), exclude(_exc), pad1(false), pad2(false), pad3(false) {}
	sqkey() : dist(0.0), pid(std::numeric_limits<unsigned int>::max()), exclude(false),
		pad1(false), pad2(false), pad3(false) {}
};
typedef sqkey * sqkeyP;


struct arrivaltime_results
{
	unsigned int gid;
	unsigned int nt;
	double intval;
	double mxsize;
	arrivaltime_results( unsigned int _gid, unsigned int _nt, double _i, double _m) :
		gid(_gid), nt(_nt), intval(_i), mxsize(_m) {}
	arrivaltime_results() : gid(THE_DOMAIN_ITSELF), nt(0), 
		intval(std::numeric_limits<double>::lowest()), 
		mxsize(std::numeric_limits<double>::lowest()) {}
};



struct sqb
{
	unsigned int nx;
	unsigned int ny;
	
	unsigned int nz;
	unsigned int nxy;
	
	unsigned int nxyz;
	unsigned int np;
	
	sqb(unsigned int _nx, unsigned int _ny, unsigned int _nz, unsigned int _np) :
		nx(_nx), ny(_ny), nz(_nz), nxy(_nx*_ny), nxyz(_nx*_ny*_nz), np(_np) {}
	sqb() : nx(1), ny(1), nz(1), nxy(1), nxyz(1), np(0) {}
};



struct idtriplet
{
	unsigned int u;
	unsigned int v;
	unsigned int w;
	idtriplet(unsigned int _u, unsigned int _v, unsigned int _w) : u(_u), v(_v), w(_w) {}
	idtriplet() : u(0), v(0), w(0) {}
};
typedef idtriplet * idtripletP;


struct triangle
{
	//###MK::check implementation later
	double x1;
	double y1;
	double x2;
	double y2;
	double x3;
	double y3;
	triangle(double _x1, double _y1, double _x2, double _y2, double _x3, double _y3) :
		x1(_x1), y1(_y1), x2(_x2), y2(_y2), x3(_x3), y3(_y3) {}
	triangle() : x1(-1.0), y1(-1.0), x2(-1.0), y2(-1.0), x3(-1.0), y3(-1.0) {}
};
typedef triangle * triangleP;


struct barycenter2d
{
	double x;
	double y;
	unsigned int gid;
	unsigned int pad;
	barycenter2d(double _x, double _y, unsigned int _gid) : x(_x), y(_y), gid(_gid), pad(0) {}
	barycenter2d() : x(0.0), y(0.0), gid(THE_DOMAIN_ITSELF), pad(0) {}
};
typedef barycenter2d * barycenter2dP;


struct aabb2d
{
	double xmi; //lower left corner
	double ymi;
	double xmx; //diagonal upper right corner
	double ymx;
	unsigned int gid;
	unsigned int pad;
	aabb2d(double _xmi, double _ymi, double _xmx, double _ymx, unsigned int _gid) : 
		xmi(_xmi), ymi(_ymi), xmx(_xmx), ymx(_ymx), gid(_gid), pad(0) {}
	aabb2d() : xmi(std::numeric_limits<double>:: max()), ymi(std::numeric_limits<double>:: max()), 
		xmx(std::numeric_limits<double>:: lowest()), ymx(std::numeric_limits<double>:: lowest()), 
		gid(THE_DOMAIN_ITSELF), pad(0) {}
};
typedef aabb2d * aabb2dP;


struct overlap2d
{
	double area;
	double leakage;		//if area computation was unsuccessful, i.e. returnINTERSECTION_ERROR
						//we return the area of the triangle to account in a threadsafe manner
						//for the accumulation of potential leaks in the reference area
	overlap2d( double _res, double _leak ) : 
		area(_res), leakage(_leak) {}
	overlap2d() : 
		area(0.0), leakage(0.0) {}
};
typedef overlap2d * overlap2dP;


struct analysismetrics2d
{
	double TotalAreaEnclosedByContour;			//total area enclosed by all contour lines from 
	double TotalAreaTriangularizable;			//total area triangularizable
	double TotalAreaNonTriangularizable;		//...
	double OverlapLeakageTotal;					//total area of triangles for which leakage occurred
	unsigned long OverlapLeakageCnts;			//total number of incipiences of leakages
	unsigned long TotalTriangleCnts;
	unsigned long TotalPolygonsTriangularized;	//...successfully
	unsigned long TotalPolygonsNonTriangularizable;
	unsigned long TotalTriangleVisitsOverlap;
	analysismetrics2d(double _Aen, double _A1tri, double _A0tri, double _Oleak, unsigned long _Oleakcnts, 
		unsigned long _Tcnts, unsigned long _P1cnts, unsigned long _P0cnts, unsigned long _Visits) :
		TotalAreaEnclosedByContour(_Aen), TotalAreaTriangularizable(_A1tri), TotalAreaNonTriangularizable(_A0tri), 
			OverlapLeakageTotal(_Oleak), OverlapLeakageCnts(_Oleakcnts), 
			TotalTriangleCnts(_Tcnts), TotalPolygonsTriangularized(_P1cnts), 
			TotalPolygonsNonTriangularizable(_P0cnts), TotalTriangleVisitsOverlap(_Visits) {}
	analysismetrics2d() : 
		TotalAreaEnclosedByContour(0.0), TotalAreaTriangularizable(0.0), TotalAreaNonTriangularizable(0.0),
			OverlapLeakageTotal(0.0), OverlapLeakageCnts(0), 
			TotalTriangleCnts(0), TotalPolygonsTriangularized(0), TotalPolygonsNonTriangularizable(0),
			TotalTriangleVisitsOverlap(0) {}
};	


struct barycenter3d
{
	double x;
	double y;
	double z;
	unsigned int gid;
	//4B padding
	barycenter3d(double _x, double _y, double _z, unsigned int _gid) : 
	x(_x), y(_y), z(_z), gid(_gid) {}
	barycenter3d() : x(0.0), y(0.0), z(0.0), gid(THE_DOMAIN_ITSELF) {}
};
typedef barycenter3d * barycenter3dP;


struct tetrahedron
{
	//###MK::check implementation later
	double x1;
	double y1;
	double x2;
	double y2;
	double x3;
	double y3;
	double x4;
	double y4;
	tetrahedron(double _x1, double _y1, double _x2, double _y2, double _x3, double _y3, double _x4, double _y4) :
		x1(_x1), y1(_y1), x2(_x2), y2(_y2), x3(_x3), y3(_y3), x4(_x4), y4(_y4) {}
	tetrahedron() : x1(-1.0), y1(-1.0), x2(-1.0), y2(-1.0), x3(-1.0), y3(-1.0), x4(-1.0), y4(-1.0) {}
};
typedef tetrahedron * tetrahedronP;


struct aabb3d
{
	double xmi; //lower left corner
	double ymi;
	double zmi;
	double xmx; //diagonal upper right corner
	double ymx;
	double zmx;
	unsigned int gid;
	//4B padding
	aabb3d(double _xmi, double _ymi, double _zmi, double _xmx, double _ymx, double _zmx, unsigned int _gid ) : 
		xmi(_xmi), ymi(_ymi), zmi(_zmi), xmx(_xmx), ymx(_ymx), zmx(_zmx), gid(_gid) {}
	aabb3d() : xmi(std::numeric_limits<double>:: max()), ymi(std::numeric_limits<double>:: max()), 
		zmi(std::numeric_limits<double>:: max() ), xmx(std::numeric_limits<double>:: lowest()), 
		ymx(std::numeric_limits<double>:: lowest()), zmx(std::numeric_limits<double>:: lowest()), gid(THE_DOMAIN_ITSELF) {}
};
typedef aabb3d * aabb3dP;


struct rectangle {
	double ox;	//origin, i.e. xmin,ymin
	double oy;
	double dx;	//diagonal upper edge xmax, ymax, i.e. clock-wise ordering
	double dy;
	rectangle(double _ox, double _oy, double _dx, double _dy) : ox(_ox), oy(_oy), dx(_dx), dy(_dy) {}
	rectangle() : ox(std::numeric_limits<double>:: max()), oy(std::numeric_limits<double>:: max()), 
		dx(std::numeric_limits<double>:: lowest()), dy(std::numeric_limits<double>:: lowest()) {}
};
typedef rectangle * rectangleP;


struct affine {
	double a;	//image scaling factor to arrive at unit circle
	double tx;	//translation vector
	double ty;
	affine(double _a, double _tx, double _ty) : a(_a), tx(_tx), ty(_ty) {}
	affine() : a(1.0), tx(0.0), ty(0.0) {} //neutral elements of an affine transformation
};
typedef affine * affineP;

struct chi
{
	double x;			//Center of the analysis circle, defined by initial microstructure!
	double y;
	double r;			//Radius of analysis circle, area is implicit as PI*SQR(r)
	double Aself;		//area which the target grain's polygon covers of this circle
	double Anbors;		//total area which the neighbors protruding into the circle cover
	double Adiff;		//numerical PI*SQR(r)-Aself-Anbors, keeps track of numerics, i.e. deviation of computed total and expectation
	double Pdrv;		//\sum_{i=1}^{N} 0.5Gb^2*(\rho_{nbor} - \rho_{me}) * A_{nbor}^{circle}/(Acircle-Aself) effective area weighted driving force distribution, positive for SEE-driving growth of target, negative for SEE-driving shrinkage
						//unit: N/m^2*m^2*1/m^2 = J/m^3
	double Pspd;		//\sum_{i=1}^{N} m_{me,nbors} 0.5Gb^2*(\rho_{nbor} - \rho_{me}) * A_{nbor}^{circle}/(Acircle-Aself) effective area weighted migration speed distribution 
						//unit: J/m^3*m^4/Js = m/s
	double Pdis;		//\sum_{i=1}^{N} m_{me,nbors} * A_{nbor}^{circle}/(Acircle-Aself) effective area weighted mobility distribution (m^4/Js)
	unsigned int gid;
	unsigned int pad;	//4B padding
	chi( double _x, double _y, double _r, double _as, double _an, double _ad, double _Pdrv, double _Pspd, double _Pdis, unsigned int _gid) :
		x(_x), y(_y), r(_r), Aself(_as), Anbors(_an), Adiff(_ad), Pdrv(_Pdrv), 
			Pspd(_Pspd), Pdis(_Pdis), gid(_gid), pad(0) {}
	chi() : x(-1.0), y(-1.0), r(-1.0), Aself(-1.0), Anbors(-1.0), Adiff(-1.0), 
		Pdrv( std::numeric_limits<double>::lowest() ), 
		Pspd( std::numeric_limits<double>::lowest() ), Pdis( std::numeric_limits<double>::lowest() ),
		gid(THE_DOMAIN_ITSELF), pad(0) {}
};
typedef chi * chiP;


struct mind //helper, bookkeepping structure for MC solver
{
	unsigned int gid;
	unsigned int cnt;
	mind( unsigned int _gid, unsigned int _cnt) : gid(_gid), cnt(_cnt) {}
	mind() : gid(THE_DOMAIN_ITSELF), cnt(0) {}
};


struct xi
{
	unsigned int x;			//Center of the analysis circle
	unsigned int y;
	unsigned int z;
	unsigned int Vtotal;	//discrete cnts of sphere
	unsigned int Vself;		//discrete cnts grain self
	unsigned int Vnbors;	//discrete cnts neighbors
	unsigned int gid;
	
	unsigned int pad;		//4B padding

	double r;				//Radius of analysis sphere, volume is discrete converging to 4.0/3.0*PI*CUBE(r)*CUBE(1.0/Settings::CurrentDomainEdgeLength)
	double Pdrv;			//\sum_{i=1}^N 0.5Gb^2*(\rho_{nbor} - \rho_{me}) * V_{nbor}^{sphere}/(Vsphere-Vself)
							//this is an unbiased volume-averaged driving force prognosis
							//unit is N/m^2*m^2*1/m^2*1 = N/m^2 = Nm/m^3=J/m^3
	double Pspd;			//\sum_{i=1}^N m_{me,nbor}*0.5Gb^2*(\rho_{nbor} - \rho_{me}) * V_{nbor}^{sphere}/(Vsphere-Vself)
							//m mobility between the analysis grain me and the n(eigh)bor and the enclosed volume of nbor in the sphere
							//Vsphere-Vself is the "volume" count covered discretely by the neighbors
							//this is an unbiased volume-averaged integral migration speed prognosis, i.e. in addition to Pdrv
							//taking into account that the central grain me may not be able to feed from a driving force reservoir
							//because it can sweep it only with a low-mobile boundary...
							//unit is : m^4/Js*N/m^2*m^2*1/m^2*1 = m/s
	double Pdis;			//\sum_{i=1}^N m_{me,nbor} * V_{nbor}^{sphere}/(Vsphere-Vself) effective mobility environment
	xi(unsigned int _x, unsigned int _y, unsigned int _z, unsigned int _vt, unsigned int _vs, unsigned int _vn, 
		unsigned int _gid, double _r, double _Pdrv, double _Pspd, double _Pdis) : x(_x), y(_y), z(_z), 
		Vtotal(_vt), Vself(_vs), Vnbors(_vn), gid(_gid), pad(0), r(_r), Pdrv(_Pdrv), Pspd(_Pspd), Pdis(_Pdis) {}
	xi() : x(0), y(0), z(0), Vtotal(0), Vself(0), Vnbors(0), gid(THE_DOMAIN_ITSELF), pad(0), r(0.0), 
		Pdrv( std::numeric_limits<double>::lowest() ), Pspd( std::numeric_limits<double>::lowest() ), Pdis( std::numeric_limits<double>::lowest() ) {}
};
typedef xi * xiP;


struct fastprops
{
	double q0;
	double q1;
	double q2;
	double q3;
	double see;
	fastprops(double _q0, double _q1, double _q2, double _q3, double _see) : 
		q0(_q0), q1(_q1), q2(_q2), q3(_q3), see(_see) {}
	fastprops() : q0(1.0), q1(0.0), q2(0.0), q3(0.0), see(0.0) {}
};
typedef fastprops * fastpropsP;


struct UDS_Grain3DInfo
{
	double q0;
	double q1;
	double q2;
	double q3;
	double see;
/*	unsigned int xmi;
	unsigned int xmx;
	unsigned int ymi;
	unsigned int ymx;
	unsigned int zmi;
	unsigned int zmx;*/
	unsigned int gid;
	//4B padding
	UDS_Grain3DInfo(double _q0, double _q1, double _q2, double _q3, double _see, unsigned int _gid) : 
		q0(_q0), q1(_q1), q2(_q2), q3(_q3), see(_see), gid(_gid) {}
	UDS_Grain3DInfo() : q0(1.0), q1(0.0), q2(0.0), q3(0.0), see(0.0), gid(THE_DOMAIN_ITSELF) {}
};
typedef UDS_Grain3DInfo * UDS_Grain3DInfoP;


struct GBContourPoint {
	double x;
	double y;
	double energy;
	double mobility;
	unsigned int myid;
	unsigned int nborid;
	GBContourPoint() : x(0.0),y(0.0),energy(0.0),mobility(0.0),myid(THE_DOMAIN_ITSELF),nborid(THE_DOMAIN_ITSELF) {}
	GBContourPoint( double _x, double _y, double _en, double _mob,
		unsigned int _myid, unsigned int _nborid ) : 
		x(_x), y(_y), energy(_en), mobility(_mob), myid(_myid), nborid(_nborid) {};
	~GBContourPoint() {}; //
};
typedef GBContourPoint * GBContourPointP;


struct GBJunctionPoint {
	double x;
	double y;
	unsigned int myid;
	unsigned int jtype;
	GBJunctionPoint() : x(0.0),y(0.0),myid(THE_DOMAIN_ITSELF),jtype(UNKNOWN_JUNCTION_TYPE) {}
	GBJunctionPoint( double _x, double _y, unsigned int _myid, unsigned int _jtype ) : 
		x(_x), y(_y), myid(_myid), jtype(_jtype) {};
	~GBJunctionPoint() {};
};
typedef GBJunctionPoint * GBJunctionPointP;


struct wghtd_imcurv_res {
	double imcurv_with_tjp;					//integral mean curvature including vertices at triple junctions
	double imcurv_without_tjp;				//integral mean curvature excluding vertices at triple junctions
	double effcapindspeed_with_tjp;			//weighted integral mean curvature with tjp
	double effcapindspeed_without_tjp;		//and without
	double effcapdrvforce_with_tjp;			//weighted integral mean curvature with tjp
	double effcapdrvforce_without_tjp;		//and without
	unsigned int cp_with_tjp;				//total number of vertices
	unsigned int cp_without_tjp;			//only virtual vertices not immediately adjacent to triple junctions
	unsigned int gID;
	bool status;
	bool pad1;
	bool pad2;
	bool pad3;
	//three padding bytes remain
	wghtd_imcurv_res(double _i1t, double _i0t, double _wi1t, double _wi0t, double _xi1t, double _xi0t, unsigned int _n1t, unsigned int _n0t, unsigned int _gid, bool _stat) : 
		imcurv_with_tjp(_i1t), imcurv_without_tjp(_i0t), effcapindspeed_with_tjp(_wi1t), effcapindspeed_without_tjp(_wi0t), 
			effcapdrvforce_with_tjp(_xi1t), effcapdrvforce_without_tjp(_xi0t), cp_with_tjp(_n1t), cp_without_tjp(_n0t), gID(_gid), 
			status(_stat), pad1(false), pad2(false), pad3(false) {}
	wghtd_imcurv_res() : imcurv_with_tjp(0.0), imcurv_without_tjp(0.0), effcapindspeed_with_tjp(0.0), effcapindspeed_without_tjp(0.0), 
		effcapdrvforce_with_tjp(0.0), effcapdrvforce_without_tjp(0.0), cp_with_tjp(0), cp_without_tjp(0), gID(0), 
		status(false), pad1(false), pad2(false), pad3(false) {}
};
typedef struct wghtd_imcurv_res* wghtd_imcurv_resP;


struct imcurv_meta {
	unsigned int fid;
	unsigned int ngr_processed;
	imcurv_meta( unsigned int _fd, unsigned int _ngr) : fid(_fd), ngr_processed(_ngr) {}
	imcurv_meta() : fid(std::numeric_limits<unsigned int>::max()), ngr_processed(0) {}
	~imcurv_meta() {}
};
typedef struct imcurv_meta * imcurv_metaP;


struct capdump {
	double imcurv;
	double pcurv;
	double vcurv;
	unsigned int nvt;
	unsigned int nvv;
	capdump( double _ic, double _pc, double _vc, unsigned int _nvt, unsigned int _nvv ) :
		imcurv(_ic), pcurv(_pc), vcurv(_vc), nvt(_nvt), nvv(_nvv) {}
	capdump() : imcurv(std::numeric_limits<double>::max()), pcurv(std::numeric_limits<double>::max()), vcurv(std::numeric_limits<double>::max()),
		nvt(0), nvv(0) {}
	~capdump(){}
};
typedef struct capdump * capdumpP;


struct capav { //averaging temporary
	double imcurv;
	double pcurv;
	double vcurv;
	unsigned int nvt;
	unsigned int nvv;
	unsigned int gid;
	unsigned int n;
	capav(const unsigned int _gid) : imcurv(0.0), pcurv(0.0), vcurv(0.0), nvt(0), nvv(0), gid(_gid), n(0) {}
	~capav(){}

	void add( const capdump _dmp) {
		this->imcurv += _dmp.imcurv;
		this->pcurv += _dmp.pcurv;
		this->vcurv += _dmp.vcurv;
		this->nvt += _dmp.nvt;
		this->nvv += _dmp.nvv;
		this->n++;
	}
};
typedef struct capav * capavP;


struct seeav {
	double psee;
	double vsee;
	unsigned int gid;
	unsigned int n;
	seeav() : psee(0.0), vsee(0.0), gid(THE_DOMAIN_ITSELF), n(0) {}
	seeav(const unsigned int _gid) : psee(0.0), vsee(0.0), gid(_gid), n(0) {}
	~seeav(){}

	void add( const double _pincr, const double _vincr ) {
		this->psee += _pincr;
		this->vsee += _vincr;
		this->n++;
	}
};
typedef struct seeav * seeavP;

#endif