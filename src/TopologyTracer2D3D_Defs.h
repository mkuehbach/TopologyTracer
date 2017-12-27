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


#ifndef __TOPOLOGYTRACER2D3D_DEFS_H__
#define __TOPOLOGYTRACER2D3D_DEFS_H__


#define THREEDIMENSIONAL			3
#define TWODIMENSIONAL				2
#define FLATTEN_THIRD_DIMENSION		(0.0)

//several other macros
#define SQR(a)						((a)*(a))
#define CUBE(a)						((a)*(a)*(a))


//MPI-related
#define SINGLE_PROCESS				1
#define MASTER						0

//MPI status handling
#define	DOING_WELL					1			//leave as they are for error management!
#define I_CRASHED					0


//MPI file IO offset increments
#define MPIIO_OFFSET_INCR_UINT		4
#define MPIIO_OFFSET_INCR_DOUBLE	8


//identifying names for clarity of the code
#define PHI1						0
#define PHI							1
#define PHI2						2

//filesize
#define TEXTUREFILE					0
#define FACEFILE					1


//flags
#define NOT_ASSIGNED_YET			2000000000
#define NO_GRAINS					0
#define NO_FACES					0
#define NOT_YET_KNOWN				0
#define NOT_ACCESSIBLE_YET			2111111111

#define FIRST_REGION				0
#define UNKNOWN_ID					0
#define NOT_YET_FOUND				2000000000


#define NOT_ONEOFMINE				3000000000

#define MAXIMUM_NUMBER_OF_RANKS		1000000000
#define UNABLE_TO_IDENTIFY_RANK		1000000001	//MK:: > MAXIMUM_NUMBER_OF_RANKS required!


#define THE_DOMAIN_ITSELF			0
#define UNKNOWN_JUNCTION_TYPE		0
#define ONE_FOR_THE_DOMAIN			1
#define BOUNDARY_CONTACT			1

//empirical factors
#define	EMPIRICAL_OVERHEAD_FACTOR	(0.95)		//accounts for the fact that a snapshot-based partitioning of the data will not comply exactly with maximum capacity

//MK::definition of the relevant binary files
//Faces_<no>.bin, no header but simply block of MPI_3DFaceIO structs
//Texture_<no>.bin, no header but simply block of MPI_3DTextureIO structs

//index handling and memory regions
#define INVALID_MEMORYREGION		2000000000	//compatible with both uint32 and and int32 and very likely neither so many ranks nor threads nor grains in the near future...
#define INVALID_INDEX				2000000000
#define UNKNOWN_CANDIDATE			2000000000	//very likely not so many MPI ranks


#define DO_NOT_ANALYZE		0x00
#define DO_ANALYZE			0x01

#define OMP_GET_NUM_THREADS					1
#define OMP_GET_NUM_THREAD					1

//return values for descriptive statics when an internal error occurred


//return values for grain properties when an internal error occurred
#define GRAINVOLUME_ON_ERROR		(0.0)
#define NFACES_ON_ERROR				0
#define NFACESDIFF_ON_ERROR			0
//#define SEE_ON_ERROR				(0.0)
//#define MEANDISORI_ON_ERROR			(0.0)
//#define MEANSEE_ON_ERROR			(0.0)
#define HAGBFRACTION_ON_ERROR		(0.0)
#define MOBDSEE_ON_ERROR			(0.0)
#define DSEE_ON_ERROR				(0.0)

//return values for RXSTATS when internal error occurred
#define VOL_ON_ERROR				(0.0)
#define CNT_ON_ERROR				0
#define X_ON_ERROR					(0.0)
#define RXEVO_EPSILON				(1.0e-14)

//trajectory tracing
#define INVALID						0

//quantile value descriptive statistics
#define MINIMUM_SAMPLE_SIZE_FOR_VARIANCE	100		//assuming that beyond this population size the Student t distribution converged already significant to the normal distribution MK::must not be below 2!
#define QUANTILES_EPSILONLIMIT		(1.0e-14)


//evaluation results
#define CRITERION_METHODS		8

#define R_RC_BAIHIR_HAGB		0
#define R_RC_BAIHIR_ALL			1
#define R_RC_BATHUT_HAGB		2
#define R_RC_BATHUT_ALL			3

#define R50_RC_BAIHIR_HAGB		4
#define R50_RC_BAIHIR_ALL		5
#define R50_RC_BATHUT_HAGB		6
#define R50_RC_BATHUT_ALL		7


//k-nearest neighbor search
#define DEFAULT_KNN_BUFFERSIZE		1024		//choose larger than an unlikely number of neighbors likely split 32KB in 1/8 (8-way associativity!)
#define DEFAULT_KNN_RECACHE			1024
#define REASONABLE_CLEARING_VALUE	0

#define KNN_STATS_SUSPICIOUS		2
#define KNN_STATS_EXPELLED			1
#define KNN_STATS_TOTAL				0

//k-nearest neighbor region properties
#define INVALID_LENGTH				(-666.666)


//Bailey Hirsch nuclei analysis
#define BH_RADIUS_INVALID			(0.0)
#define INFINITE_RADIUS				(1.0e24)
#define MAXIMUM_NUMBER_OF_TARGETS	(15000000)	//to avoid that MAXIMUM_NUMBER_OF_TARGETS * CRITICAL RADIUS CONCEPTS exceeds 2^32-1
#define BH_MINIMUM_RHO				(1.0e0)

//abnormal grains
#define CRITICAL_RADIUS_RATIO		(3.0)		//when grains are categorized as abnormal grains	
//#define METER2MICRON				(1.0e-6)

//statistics
#define STATISTICS_VOL_NQUANTILES	100			//0.01 binning, i.e. 1 percentile increments


//MODF computations
#define MODF_EPSILON				(1.0e-9)
#define MODF_INFTY					((2.1)*(PI))

//SEE computations
#define SEE_EPSILON					(1.0e-9)
#define SEE_INFTY					(1.0e20)

//GSD computations
#define GSD_EPSILON					(1.0e-9)
#define GSD_INFTY					(1.0e3)

//LongRange Persistence analysis
//#define PERSISTENCE_VERBOSE		//define when interested in prompting triangle count per grain
#define MINIMUM_REQUIRED_TOTALAREA					(0.99)	//absolute values below which to much inaccuracy for unbiased statistics above which too much geometrical inaccuracy
#define MAXIMUM_ADMISSIBLE_NONTRIANGULARIZED		(0.01)




#define INVALID_HASHKEY				4096000666

//SI unit scaling
#define METER2MICRON(x)				((x)*(1.0e6))
#define MICRON2METER(x)				((x)*(1.0e-6))

//conversions
#define DEG2RAD(deg)					((deg)/(180.0)*(PI))
#define RAD2DEG(rad)					((rad)*(180.0)/(PI))

#endif