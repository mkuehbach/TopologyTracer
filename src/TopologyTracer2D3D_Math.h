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

#ifndef __TOPOLOGYTRACER2D3D_MATH_H__
#define __TOPOLOGYTRACER2D3D_MATH_H__

#include "TopologyTracer2D3D_MPIDatatypes.h"

#define PI						(3.141592653589793238462643383279502884197169399375105820974)
#define NPRIMES					32
#define DOUBLE_EPSILON			(1.0e-15)

//inheriting from <float.h>
#define DBL_EPSILON				2.2204460492503131e-016 /* smallest such that 1.0+DBL_EPSILON != 1.0 */ 
#define MIN(a,b)				(((a)<(b))?(a):(b))
#define MAX(a,b)				(((a)>(b))?(a):(b))

#define NO_SOLUTION				0
#define ONE_SOLUTION			1
#define TWO_SOLUTIONS			2
#define INVALID_SOL				(-1.0) //MK::do not set to on [0,0, 1,0]

//specific defines for point in triangle test
#define TRIANGLE_EPSILON				1.0e-12
#define PQ_EPSILON						1.0e-32
#define COLLINEARITY_ACCURACY			(1.0e-12)

//triangle detection error code
#define CRITICAL_WARNING_LEVEL			(1.0e-8)
#define INTERSECTION_DETECTION_ERROR	(-666.0)		//MK::required to be negative!
//LongRangePersistence
#define	MINIMUM_CRITICAL_INTERSECTION_AREA	( ((1.0)/(100000.0)) * ((1.0)/(100000.0)) ) //below which triangle area do we no longer attempt to compute the intersection area of the circle and triangle in particular the cap areas very accurately?

//handle numerics of small numbers heuristically?
#define	HEURISTIC_HANDLING_OF_SMALL_NUMBERS


#define DISTANCE_TOLERANCE				(1.0e-14)

struct pqsolve
{
	double sol1;
	double sol2;
	unsigned int nsol;
	pqsolve() : sol1(INVALID_SOL), sol2(INVALID_SOL), nsol(NO_SOLUTION) {}
};
typedef pqsolve * pqsolveP;


//Math library version defining conventions of Euler angles and rotations
#define SCORE_LIBRARY		//chose either or
//#define DEGRAEF_LIBRARY

class mathMethods
{
public:
	mathMethods();
	~mathMethods();

	void multiplyQuaternions( double *q, double* p, double* r );
	void euler2quaternion( double * euler, double * q );
	double misorientationCubic( double pa1, double Pa, double pa2, double pb1, double Pb, double pb2 );
	double misorientationCubic( double* p, double* q );

	double Heaviside( double x );
	double Heaviside1Minusf( double nbrho, double merho );

	double dis2mob( double theta ) {
		return ( 1.0 - (0.99 * exp( -5.0 * (pow(theta / DEG2RAD(15.0), 9.0)))) );
	}

	bool WorkPartitioning( unsigned int grainID, unsigned int nthr, unsigned int mythr ) {
		if ( grainID % nthr == mythr )	return true;
		else							return false;
	}

	//functions for calculating exact overlap between shapes is under a 3-clause BSD style license
	//original cython version by Thomas Robitaille. Converted to C by Kyle Barbary. */

	struct pqsolve stable_pqsolver( double a, double b, double c );

	//geometry
	void outer_triedge_normal( double x1, double y1, double x2, double y2, double x3, double y3, double bx, double by, double tcandx, double tcandy, double* n );
	void outer_quadedge_normal( double centerx, double centery, double bx, double by, double tcandx, double tcandy, double* n );
	bool point_in_circle( double x, double y, double xc, double yc, double r );
	double distance_point_to_edge(double x, double y, double x1, double y1, double x2, double y2);
	bool point_in_triangle_fast( double x, double y, double x1, double y1, double x2, double y2, double x3, double y3 );
	bool point_in_triangle( double x, double y, double x1, double y1, double x2, double y2, double x3, double y3 );
	double circular_cap_area( double h, double r );
	double triangle_area( double ux, double uy, double vx, double vy );

	double helper_bisec_a_t( double nx, double ny );
	double helper_bisec_b_t( double bx, double by, double nx, double ny, double xc, double yc );
	double helper_bisec_c_t( double bx, double by, double xc, double yc, double R );
	double helper_intersect_a_t( double tx, double ty );
	double helper_intersect_b_t( double xx1, double yy1, double tx, double ty, double xc, double yc );
	double helper_intersect_c_t( double xx1, double yy1, double xc, double yc, double R);

	double tri_cir_v0e1( double xa, double ya, double xb, double yb, double x, double y, double r, double xc, double yc );
	double tri_cir_v0e2( double xa, double ya, double xb, double yb, double xc, double yc, double x, double y, double r, bool ab, bool bc, bool ca );
	double tri_cir_v0e3( double xa, double ya, double xb, double yb, double xc, double yc, double x, double y, double r );
	double tri_cir_v1e2( double xa, double ya, double xb, double yb, double xc, double yc, double x, double y, double r, bool a, bool b, bool c, bool ab, bool bc, bool ca );
	double tri_cir_v1e3( double xa, double ya, double xb, double yb, double xc, double yc, double x, double y, double r, bool a, bool b, bool c, bool ab, bool bc, bool ca );
	double tri_cir_v2e3(  double xa, double ya, double xb, double yb, double xc, double yc, double x, double y, double r, bool a, bool b, bool c, bool ab, bool bc, bool ca );
	double intersection_area_triangle_circle( double x1, double y1, double x2, double y2, double x3, double y3, double x, double y, double r );
	double intersection_area_triangle_circle_mc( double x1, double y1, double x2, double y2, double x3, double y3, double x, double y, double r, double ld );
	bool intersection_area_triangle_check( double A, double x1, double y1, double x2, double y2, double x3, double y3, double x, double y, double r );
};

#endif


