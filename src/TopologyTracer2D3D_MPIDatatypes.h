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


#ifndef __TOPOLOGYTRACER2D3D_MPIDATATYPES_H__
#define __TOPOLOGYTRACER2D3D_MPIDATATYPES_H__

#include "TopologyTracer2D3D_Datatypes.h"

//MPI IO container for file interaction
typedef struct
{
	unsigned int extGID;
	unsigned int whenhitboundary;	//snapshot at which the grain hits the boundary, default value std::numeric_limits<unsigned int>::max() means never
} MPI_GrainBndContactIO;


typedef struct 
{
	double volume;
	double perimeter;
	double GBEnergy;
	double BulkEnergy;

	double phi1;
	double Phi;
	double phi2;

	double x;
	double y;

	unsigned int id;
	unsigned int NeighbourCount;
	unsigned int intersectsBoundaryGrain;
	unsigned int pad;
} MPI_2DGrainIO;


typedef struct
{
	double size;								//size is a grain area for 2D sim, and the grain volume in the 3D sim
	double surfaceArea;
	double GBEnergy;
	double BulkEnergy;

	double phi1;								//Bunge, "ZXZ" Euler angle in entire Euler space
	double Phi;
	double phi2;

	double x;									//barycenter of the grain
	double y;
	double z;									//is 0.0 for a 2D simulation

	unsigned int id;
	unsigned int NeighbourCount;
	unsigned int intersectsBoundaryGrain;		//1 yes (grain has contact with the boundary), 0 no (inside bulk)
	unsigned int pad;							//MK::currently not occupied
} MPI_3DGrainIO;


typedef struct
{
	double size;								//size is a grain boundary length for 2D sim, and an area for 3D sim
	unsigned int gA;							//identifier of the grain gA > gB and gA, gB \in \mathbb{N} but no one 0
	unsigned int gB;
} MPI_3DFaceIO;


typedef struct
{
	double size_ng;								//file size in byte for Texture_* bin contain grains and
	double size_nf;								//Faces_*.bin containing faces
	unsigned int extSID;						//a once externally defined ID for the snapshot
	unsigned int ng;							//how many grains in this snapshot?
	unsigned int nf;							//how many faces in this snapshot?
	unsigned int RID;

	unsigned int nmrx;							//region partitioning
	unsigned int nmry;
	unsigned int nmrz;
	unsigned int nmrxyz;
	//MPI_SMDIO( double _szng, double _sznf, unsigned int _esid, unsigned int _ng, unsigned int _nf, unsigned int _rid,
	//	unsigned int _nmrx, unsigned int _nmry, unsigned int _nmrz, unsigned int _nmrxyz) : size_ng(_szng), size_nf(_sznf), 
	//	extSID(_esid), ng(_ng), nf(_nf), RID(_rid), nmrx(_nmrx), nmry(_nmry), nmrz(_nmrz), nmrxyz(_nmrxyz) {}
	//MPI_SMDIO() : size_ng(0.0), size_nf(0.0), extSID(UNKNOWN_ID), ng(0), nf(0), RID(UNKNOWN_ID), nmrx(0), nmry(0), nmrz(0), nmrxyz(0) {}
	//region partitioning not transferred because in all practical cases the same for each Snapshot
	//the rank in which the snapshot is stored
} MPI_SMDIO;


typedef struct
{
	double totalsize_ng_elimbnd;
	double meansize_ng_elimbnd;
	double varsize_ng_elimbnd;

	double ng; //unsigned int ng;
	double ng_elimbnd; //unsigned int ng_elimbnd;
} MPI_VOLSTATS;


typedef struct
{
	double totalsize_elimbnd;					//domain area covered by grains without boundary contact
	double totalsize_targets;					//domain area covered by the target grains
	double X;									//X := totalsize_targets/totalsize_elimbnd;
	double ndefg;								//how many matrix grains, double values to aid reading in matlab
	double nrxg;								//how many considered as to be rxed grains
} MPI_RXSTATS;


typedef struct
{
	double totalsize_survivors;
	double totalsize_matrix;
	double totalsize_targets;

	double rmean_survivors;			// EHR of survivors
	double rmean_matrix;			// EHR of matrix grains
	double rmean_targets;			// EHR of target grains, there will be a difference to the rmean_matrix because the targets account also for the abnormal grains and hence bias the mean value

	//counts also in double to simplify reading in with MATLAB
	double n_survivors;				// how many grains still there consistency
	double n_matrix;				// will decrease strong than from global coarsening because population narrows down to survivors only
	double n_targets;				// bound to n_targets > n_matrix

	double n_agg_m;					// with respect to mean of matrix
	double n_agg_t;					// with respect to mean of targets
	double n_agg_m_insurv;			// ... also in the list of survivors
	double n_agg_t_insurv;			// ... also in ...
} MPI_AGGSTATS;


typedef struct
{
	double bx;
	double by;
	double gid;
} MPI_GNUInfo;


typedef struct
{
	double x;
	double y;
	double en;
	double mob;
	unsigned int gid;
	unsigned int nborid;
} MPI_GBContourPoint;


typedef struct
{
	double x;
	double y;
	unsigned int gid;
	unsigned int jtype;
} MPI_GBJunctionPoint;


typedef struct
{
	double imcurv;
	double pcurv;
	double vcurv;
	unsigned int nvt;
	unsigned int nvv;
	unsigned int gid;
	unsigned int padding;
} MPI_CapillaryDump;

#endif