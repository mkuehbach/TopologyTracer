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

//*********************************************************************************************************************************************************************************************
#ifndef __TOPOLOGYTRACER2D3D_SETTINGS_H__
#define __TOPOLOGYTRACER2D3D_SETTINGS_H__

//mount entire STL stuff
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstring>
#include <stdexcept>
#include <iomanip>
#include <string>

#include <cmath>
#include <vector>
#include <bitset>
#include <limits>
#include <algorithm>

#include "TopologyTracer2D3D_Defs.h"

//parallelization and thirdparties
#include <mpi.h>
#include <omp.h>

#include "thirdparty/RapidXML/rapidxml.hpp"

using namespace std;
using namespace rapidxml;



//global user input validity check limits
#define PI								(3.141592653589793238462643383279502884197169399375105820974)
#define LOCALDBMAXIMUMSIZE				((90.0) * (1024.0) * (1024.0) * (1024.0) )		//Byte

#define MPIIO_DEFAULT_BLOCKLENGTH		((100.0) * (1024.0) * (1024.0))					//Byte
#define EPSILON							(1.0e-12)										//meter
#define EMPIRICAL_LOOKUP_MAXIDRANGE		100
#define EMPIRICAL_LOOKUP_MINIDRANGE		10
#define EMPIRICAL_LARGEST_GRAINID		1000000000

//MODF related
#define DEFAULT_MINIMUM_DISORI_FCC		(0.0)											//deg
#define DEFAULT_MAXIMUM_DISORI_FCC		(62.8)
#define DEFAULT_DISORI_BINWIDTH			(0.1)

//SEE related
#define	DEFAULT_MINIMUM_SEE				(1.0e11)
#define DEFAULT_MAXIMUM_SEE				(1.0e17)
#define DEFAULT_SEE_BINWIDTH			(1.0e11)

//GSD normalized to instantaneous mean related
#define DEFAULT_MINIMUM_GSD				(0.0)
#define DEFAULT_MAXIMUM_GSD				(27.0)
#define DEFAULT_GSD_BINWIDTH			(0.01)

//k-nearest neighbors related
#define	DEFAULT_KSHELL_MAX				16
#define DEFAULT_HAGBTHRESHOLD_MIN		(0.0)
#define DEFAULT_HAGBTHRESHOLD_MAX		(62.8)
#define DEFAULT_HAGBTHRESHOLD			((15.0) / (180.0) * (PI))

//dimensionality
#define TWO_DIMENSIONS					2
#define THREE_DIMENSIONS				3


//LongRange related
#define E_LEVELSET_RVE_EXTENT_MIN		10
#define E_LEVELSET_RVE_EXTENT_MAX		10000			//##MK::very no 3d Levelset larger than this^3 likely not in the next years....

enum E_COARSENING_MODEL {
	E_LEVELSET				//C Mie{\ss}en et. al. 2d and 3d level-set code GraGLeS
	//add further models here
};


enum E_ANALYSIS_MODE {
	E_UNKNOWN_MODE,
	E_TRACKING_PARALLEL,		//track grain evolution with multiple MPI processes
	E_TRACKING_SEQUENTIAL,		//track grains with only one MPI process
	E_KNNCORRELATION,			//
	E_NUCSITES_LONGRANGE2D,		//implements exemplarily Kuehbachs approach to account for long-range effects in 2D
	E_NUCSITES_LONGRANGE3D,		//implements the above in 3D
	E_APPROX_CURVATURE_CONTOUR,
	E_CAPILLARY_ACTIVITY		//beta functionality to fuse single-grain-resolved temporal trajectory of capillary activity data
};


enum E_ANALYZE_FORWARD_MODE {
	E_FORWARD_ALL,				//take all grains
	E_FORWARD_SELECTED			//only certain IDs
};


enum E_CAPTRACKING_MODE {
	E_CAP_AVERAGE,
	E_CAP_HISTORY
};


enum E_CAPTRACKING_WHICH {
	E_CAP_ALL,
	E_CAP_SEL
};


class Settings {
public:

	static E_ANALYSIS_MODE AnalysisMode;				//the work to do
	static E_COARSENING_MODEL SimModel;					//the model used to generate the raw data
	static std::string TargetGIDFromFilename;
	static std::string KNNDataFromFilename;
	static std::string GNUDataFromFilename;
	static std::string UDSDataFromFilename;
	static std::string GrainBndContactFilename;
	static unsigned int SimID;							//consistent identifier of the data analysis run
	static double LocalDatabaseMaximumSize;				//what is the maximum memory for storing data in each process

	static unsigned int SnapshotFirst;					//the data range on which the TopologyTracer should operate
	static unsigned int SnapshotOffset;
	static unsigned int SnapshotLast;

	static unsigned int Dimensionality;					//analysis scope

	static double MPIReadBlockLength;					//how many byte are read at most with a single IO call?
	static double HAGBMobility;							//maximum mobility (m^4/Js) used in simulation
	static double HAGBEnergy;							//maximum specific boundary energy (J/m^2) used in simulation
	static double DislocEnPerM;							//0.5Gbb (J/m) fixed temperature line length used in simulation
	static double PhysicalDomainSize;					//(m) how large in real space is the edge of the domain

	static unsigned int LargestGrainID;					//the largest ID in the dataset
	static unsigned int LookupMaxGrainIDRange;			//how many grain id are at most gathered in a lookup table bucket? the smaller the quantity the less checks to find grains but the more cache misses
	static unsigned int MemRegionsX;
	static unsigned int MemRegionsY;
	static unsigned int MemRegionsZ;
	static unsigned int MaxNumberOfKShells;
	static double HAGBDetectionThreshold;
	static double DisoriAngleBinMin;
	static double DisoriAngleBinMax;
	static double DisoriAngleBinWidth;					//degrees
	static double SeeBinMin;
	static double SeeBinMax;
	static double SeeBinWidth;
	static double GSDBinMin;
	static double GSDBinMax;
	static double GSDBinWidth;

	static double LongRangeRadiusMin;
	static double LongRangeRadiusIncr;
	static double LongRangeRadiusMax;
	static double InitialDomainEdgeLength;

	static bool DeveloperMode;
	static bool ProbeWorkPartMode;
	static bool ProbeBoundaryContact;

	static bool AnalyzeTrajectoriesForward;
	static bool AnalyzeTrajectoriesBackward;
	static bool AnalyzeSizeGainVsMatrixBackward;
	static bool AnalyzeMaxSizeGainForward;
	static bool AnalyzeMeanDrvForceSEEForward;
	static bool AnalyzeKNN;
	static bool AnalyzeGrainSizeQuantiles;
	static bool AnalyzeMODF;
	static bool AnalyzeSEE;
	static bool AnalyzeGSD;
	static bool AnalyzeDrvForceSEE;
	static bool AnalyzeClassicalNucModels;
	static bool AnalyzeApproxRXFraction;
	static bool AnalyzeAbnormalGrains;
	static bool AnalyzeTopologyDifferenceForward;

	static unsigned int AnalyzeTrajectoriesForwardMode;

	static bool ComputeCurvatureTJP;
	static bool TranslateBinary2GNU;

	static E_CAPTRACKING_WHICH CapillaryTrackingTargets;
	static E_CAPTRACKING_MODE CapillaryTrackingMode;

	//functions
	static void readXML(string filename = "");
	static bool checkUserInput( void );

	static unsigned int get_GSDNQuantiles(void);
	static unsigned int get_SEEBinCount(void);
	static unsigned int get_MODFBinCount(void);
	static unsigned int get_GSDBinCount(void);
	static double get_SEEBinEnd(void);
	static double get_MODFBinEnd(void);
	static double get_GSDBinEnd(void);
};

#endif
