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

#include "TopologyTracer2D3D_Math.h"

using namespace std;


mathMethods::mathMethods( void ){}


mathMethods::~mathMethods( void ){}


#ifdef SCORE_LIBRARY
void mathMethods::multiplyQuaternions( double *q, double* p, double* r )
{
	//MK::ok, mathematically multiplies quaternions q and p, active or passive rotation convention does not matter
	//verified via http://mathworld.wolfram.com/Quaternion.html as well as http://www.mathworks.de/de/help/aeroblks/quaternionmultiplication.html
	//equivalent to the multiplication given in Grimmer, H, 1974 Acta Cryst A30, 685-688
	//mind vector cross product notation qp = q0p0 - qbar*pbar + q0 pbar + p0 qbar + qbar cross pbar with bar vector quantities parts of the quaternion
	r[0] = + q[0] *	p[0]	- q[1] * p[1]	- q[2] *	p[2]	- q[3] *	p[3];
	r[1] = + q[1] *	p[0]	+ q[0] * p[1]	- q[3] *	p[2]	+ q[2] *	p[3];
	r[2] = + q[2] *	p[0]	+ q[3] * p[1]	+ q[0] *	p[2]	- q[1] *	p[3];
	r[3] = + q[3] *	p[0]	- q[2] * p[1]	+ q[1] *	p[2]	+ q[0] *	p[3];
}


void mathMethods::euler2quaternion( double * euler, double * q )
{
	//OK, 20130326MK convention: Bunge ZXZ, which represents the (3,1,3) case analyzed in: Diebel 2006 
	//Representing Attitude: Euler Angles, Unit Quaternions, and Rotation Vectors, James Diebel, Stanford University, Stanford, California 94301-9010, Email: diebel@stanford.edu
	//mind utilization of addition theorems as validated with www.wolframalpha.com 
	//cos(a+b) = c(a+b) = cacb-sasb
	//cos(a-b) = c(a-b) = cacb+sasb
	//sin(a+b) = s(a+b) = sacb+casb
	//sin(a-b) = s(a-b) = sacb-casb

	double p1 = euler[0]; //Diebel PSI
	double t  = euler[1]; //Diebel theta
	double p2 = euler[2]; //Diebel PHI

	double co1 = cos(t/2);
	double s1 = sin(t/2);

	double p[4] = {co1*cos((p1+p2)/2), s1*cos((p1-p2)/2), s1*sin((p1-p2)/2), co1*sin((p1+p2)/2)}; //applying sin, cos addition theorems

	q[0] = p[0];
	q[1] = p[1];
	q[2] = p[2];
	q[3] = p[3];
}


double mathMethods::misorientationCubic( double pa1, double Pa, double pa2, double pb1, double Pb, double pb2 )
{
	//MK::ok
	double oria[3] = { pa1, Pa, pa2 };
	double orib[3] = { pb1, Pb, pb2 };

	double p[4];
	double q[4];

	euler2quaternion( oria, p );
	euler2quaternion( orib, q );

	double qm1[4];

	//Inverse of quaternion q utilizing unity quaternions
	qm1[0] = q[0];
	qm1[1] = -1.0 * q[1];
	qm1[2] = -1.0 * q[2];
	qm1[3] = -1.0 * q[3];

	double r[4]; //Resulting quaternion, rotation of the two previous quaternions pq-1

	multiplyQuaternions( qm1, p, r );
	//MK:: was p post applied to qm1 but that was not consistent with MTex, in particular not with -(q^-1) = -q0,q1,q2,q3
	///#######MK:::::ONLY THIS BRING AGREEMENT WITH MTEX!!!!but not p, qm1 !!!!!!!


	//Now, we have to determine the smallest angle, following Grimmer H, Acta Cryst 1974, A30 pp685-688
	double r0[6][4];    //There are 12 possible angles

	//Note: The notation r0 is due to the definition of the quaternion which lie
	//in the fundamental zone, this vector possesses the smallest angle, in such a way
	//that r0 is actually the scalar part of this quaternion

	double a = r[0];
	double b = r[1];
	double c = r[2];
	double d = r[3];

	double fac = 1.0 / sqrt( 2.0 ); //0.70710678;

	//six fundamental quaternions																					//Grimmer
	r0[0][0] =(a+b)*fac; 		r0[0][1]=(a-b)*fac; 		r0[0][2] = (c+d)*fac; 		r0[0][3] = (c-d)*fac;		//7b
	r0[1][0] =(a+c)*fac; 		r0[1][1]=(a-c)*fac;	 		r0[1][2] = (b+d)*fac; 		r0[1][3] = (b-d)*fac;		//7c
	r0[2][0] =(a+d)*fac; 		r0[2][1]=(a-d)*fac; 		r0[2][2] = (b+c)*fac; 		r0[2][3] = (b-c)*fac;		//7d
	r0[3][0] =(a+b+c+d)*0.5; 	r0[3][1]=(a+b-c-d)*0.5; 	r0[3][2] = (a-b+c-d)*0.5; 	r0[3][3] = (a-b-c+d)*0.5;	//7e
	r0[4][0] =(a+b+c-d)*0.5; 	r0[4][1]=(a+b-c+d)*0.5; 	r0[4][2] = (a-b+c+d)*0.5; 	r0[4][3] = (a-b-c-d)*0.5;	//7fGrimmer
	r0[5][0] = a;				r0[5][1] = b;				r0[5][2] = c;				r0[5][3] = d;				//Grimmer7a

	double omega = 0.0;

	for(int i=0;i<6;i++)
		for( int j=0;j<4;j++ )
			if( fabs(r0[i][j]) > omega ) //as the sign changes are arbitrary
				omega=fabs(r0[i][j]);

	if( omega > 1 ) //avoid singularity of acos function
		omega = (double) (int) omega;

	omega = 2 * acos(omega);
	//QUICKASSERT( omega <= MAXIMUM_MISORI_FCC );
	return omega;
}


double mathMethods::misorientationCubic( double* p, double* q )
{
	double qm1[4];

	//Inverse of quaternion q utilizing unity quaternions
	qm1[0] = q[0];
	qm1[1] = -1.0 * q[1];
	qm1[2] = -1.0 * q[2];
	qm1[3] = -1.0 * q[3];

	double r[4]; //Resulting quaternion, rotation of the two previous quaternions pq-1

	multiplyQuaternions( qm1, p, r );
	//MK:: was p post applied to qm1 but that was not consistent with MTex, in particular not with -(q^-1) = -q0,q1,q2,q3
	///#######MK:::::ONLY THIS BRING AGREEMENT WITH MTEX!!!!but not p, qm1 !!!!!!!


	//Now, we have to determine the smallest angle, following Grimmer H, Acta Cryst 1974, A30 pp685-688
	double r0[6][4];    //There are 12 possible angles

	//Note: The notation r0 is due to the definition of the quaternion which lie
	//in the fundamental zone, this vector possesses the smallest angle, in such a way
	//that r0 is actually the scalar part of this quaternion

	double a = r[0];
	double b = r[1];
	double c = r[2];
	double d = r[3];

	double fac = 1.0 / sqrt( 2.0 ); //0.70710678;

	//six fundamental quaternions
	r0[0][0] =(a+b)*fac; 		r0[0][1]=(a-b)*fac; 		r0[0][2] = (c+d)*fac; 		r0[0][3] = (c-d)*fac;
	r0[1][0] =(a+c)*fac; 		r0[1][1]=(a-c)*fac;	 		r0[1][2] = (b+d)*fac; 		r0[1][3] = (b-d)*fac;
	r0[2][0] =(a+d)*fac; 		r0[2][1]=(a-d)*fac; 		r0[2][2] = (b+c)*fac; 		r0[2][3] = (b-c)*fac;
	r0[3][0] =(a+b+c+d)*0.5; 	r0[3][1]=(a+b-c-d)*0.5; 	r0[3][2] = (a-b+c-d)*0.5; 	r0[3][3] = (a-b-c+d)*0.5;
	r0[4][0] =(a+b+c-d)*0.5; 	r0[4][1]=(a+b-c+d)*0.5; 	r0[4][2] = (a-b+c+d)*0.5; 	r0[4][3] = (a-b-c-d)*0.5;
	r0[5][0] = a;				r0[5][1] = b;				r0[5][2] = c;				r0[5][3] = d;

	double omega = 0.0;

	for(int i=0;i<6;i++)
		for( int j=0;j<4;j++ )
			if( fabs(r0[i][j]) > omega ) //as the sign changes are arbitrary
				omega=fabs(r0[i][j]);

	if( omega > 1 ) //avoid singularity of acos function
		omega = (double) (int) omega;

	omega = 2 * acos(omega);
	//QUICKASSERT( omega <= MAXIMUM_MISORI_FCC );
	return omega;
}

#else //DEGRAEF_LIBRARY adhering to the definitions in ###
//consistency with DAMASK v2.0.1 ? ###
void mathMethods::multiplyQuaternions( double *q, double* p, double* r ) {
}


void mathMethods::euler2quaternion( double * euler, double * q ) {
}


double mathMethods::misorientationCubic( double pa1, double Pa, double pa2, double pb1, double Pb, double pb2 ) {
}


double mathMethods::misorientationCubic( double* p, double* q ) {
}

#endif

double mathMethods::Heaviside( double x ) {
	if ( x < 0 )
		return 0.0;
	//else
		return x; //implicit 1*x
}


double mathMethods::Heaviside1Minusf( double nbrho, double merho ) {
	if ( nbrho > merho && nbrho > BH_MINIMUM_RHO ) 
		return (1.0 - merho / nbrho);
	//else
		return 0.0;
}


struct pqsolve mathMethods::stable_pqsolver( double a, double b, double c ) {
	//solves ad^2+bd+c=0 for d where d is for instance a scaling coefficient to obtain an intersection point on a precomputed line
	struct pqsolve res;
	if ( fabs(a) < PQ_EPSILON ) { //a close to zero
		if ( fabs(b) < PQ_EPSILON ) //b close to zero
			return res;
		else {
			res.nsol = ONE_SOLUTION;
			res.sol1 = -1.0*c / b;
			res.sol1 = -1.0*c / b;
			return res;
		}
	}

	//classical quadratic solution equation approach
	double D = SQR(b) - 4.0*a*c;
	if ( D < 0.0 ) //res.nsol = NO_SOLUTION is set the default via struct constructor
		return res;

	if ( D < PQ_EPSILON ) { //##MK::stability?
		res.nsol = ONE_SOLUTION;
		res.sol1 = -0.5*b/a; //if only one return solution the same value
		res.sol2 = -0.5*b/a;
		return res;
	}

	res.nsol = TWO_SOLUTIONS;
	if ( b < 0.0 ) {
		double s1 = 2.0*c / (-1.0*b + sqrt(D));
		double s2 = (-1.0*b + sqrt(D)) / (2.0*a);
		res.sol1 = MIN(s1,s2);
		res.sol2 = MAX(s1,s2);
	}
	else {
		double s3 = (-1.0*b - sqrt(D)) / (2.0*a);
		double s4 = 2.0*c / (-1.0*b - sqrt(D));
		res.sol1 = MIN(s3,s4);
		res.sol2 = MAX(s3,s4);
	}
	return res;
	//solutions sorted ascendingly
}


void mathMethods::outer_triedge_normal( double x1, double y1, double x2, double y2, double x3, double y3, double bx, double by, double tcandx, double tcandy, double* n )
{
	//determine outer triedge normal to a triangle, works for convex polygons
	double center[2] = { (x1+x2+x3)/3.0, (y1+y2+y3)/3.0 };
	//dot product is negative when outer normal, numericvally zero impossible because the candidate is not a zero vector
	if ( ((center[0]-bx)*tcandy + (center[1]-by)*(-1.0)*tcandx) < 0.0 ) {
		n[0] = tcandy;
		n[1] = -1.0*tcandx;
	}
	else {
		n[0] = -1.0*tcandy;
		n[1] = tcandx;
	}
}


void mathMethods::outer_quadedge_normal( double centerx, double centery, double bx, double by, double tcandx, double tcandy, double* n )
{
	//determine outer triedge normal to a quad polygon, works for convex polygons
	//dot product is negative when outer normal, numericvally zero impossible because the candidate is not a zero vector
	if ( ((centerx-bx)*tcandy + (centery-by)*(-1.0)*tcandx) < 0.0 ) {
		n[0] = tcandy;
		n[1] = -1.0*tcandx;
	}
	else {
		n[0] = -1.0*tcandy;
		n[1] = tcandx;
	}
}


bool mathMethods::point_in_circle( double x, double y, double xc, double yc, double r )
{
	return ( ((SQR(x-xc)+SQR(y-yc)) > SQR(r)) ? false : true );
}


double mathMethods::distance_point_to_edge(double x, double y, double x1, double y1, double x2, double y2)
{
	double dist = SQR(x2 - x1) + SQR(y2 - y1);
	double dot = ((x - x1)*(x2 - x1) + (y - y1)*(y2 - y1)) / dist;

	if ( dot < 0.0 )
		return (SQR(x - x1) + SQR(y - y1));
	else if ( dot <= 1.0 )
		return ( SQR(x1-x) + (SQR(y1-y) - (SQR(dot)*dist)) );
	else
		return(SQR(x-x2) + SQR(y-y2));
}


bool mathMethods::point_in_triangle_fast( double x, double y, double x1, double y1, double x2, double y2, double x3, double y3 )
{
	//numerical safe procedure to test whether point x,y is in or exactly on the boundary of a triangle
	//based on the approach of http://totologic.blogspot.de/2014/01/accurate-point-in-triangle-test.html
	//faster because assumes that pretest against triangle AABB was already performed

	//catch most inclusions with naive point in triangle test
	double A = ((y2 - y3)*(x1 - x3) + (x3 - x2)*(y1 - y3));
	double a = ((y2 - y3)*(x - x3) + (x3 - x2)*(y - y3)) / A;
	double b = ((y3 - y1)*(x - x3) + (x1 - x3)*(y - y3)) / A;
	double c = 1.0 - a - b;
	if ( 0.0 <= a && a <= 1.0 && 0.0 <= b && b <= 1.0 && 0.0 <= c && c <= 1.0 ) 
		return true;

	//check for numerical ill defined close positioning on the edge with a finite edge thickness
	if ( distance_point_to_edge( x, y, x1, y1, x2, y2 ) <= SQR(TRIANGLE_EPSILON) ) return true;
	if ( distance_point_to_edge( x, y, x2, y2, x3, y3 ) <= SQR(TRIANGLE_EPSILON) ) return true;
	if ( distance_point_to_edge( x, y, x3, y3, x1, y1 ) <= SQR(TRIANGLE_EPSILON) ) return true;

	return false;
}


bool mathMethods::point_in_triangle( double x, double y, double x1, double y1, double x2, double y2, double x3, double y3 )
{
	//numerical safe procedure to test whether point x,y is in or exactly on the boundary of a triangle
	//based on the approach of http://totologic.blogspot.de/2014/01/accurate-point-in-triangle-test.html
	//pretest via finite AABB
	//std::numeric_limits<double>:: lowest(), max()
	if ( x < (MIN(x1,MIN(x2,x3)) - TRIANGLE_EPSILON)) return false;
	if ( x > (MAX(x1,MAX(x2,x3)) + TRIANGLE_EPSILON)) return false;
	if ( y < (MIN(y1,MIN(y2,y3)) - TRIANGLE_EPSILON)) return false;
	if ( y > (MAX(y1,MAX(y2,y3)) + TRIANGLE_EPSILON)) return false;

	//catch most inclusions with naive point in triangle test
	double A = ((y2 - y3)*(x1 - x3) + (x3 - x2)*(y1 - y3));
	double a = ((y2 - y3)*(x - x3) + (x3 - x2)*(y - y3)) / A;
	double b = ((y3 - y1)*(x - x3) + (x1 - x3)*(y - y3)) / A;
	double c = 1.0 - a - b;
	if ( 0.0 <= a && a <= 1.0 && 0.0 <= b && b <= 1.0 && 0.0 <= c && c <= 1.0 ) 
		return true;

	//check for numerical ill defined close positioning on the edge with a finite edge thickness
	if ( distance_point_to_edge( x, y, x1, y1, x2, y2 ) <= SQR(TRIANGLE_EPSILON) ) return true;
	if ( distance_point_to_edge( x, y, x2, y2, x3, y3 ) <= SQR(TRIANGLE_EPSILON) ) return true;
	if ( distance_point_to_edge( x, y, x3, y3, x1, y1 ) <= SQR(TRIANGLE_EPSILON) ) return true;

	return false;
}


double mathMethods::circular_cap_area( double h, double r ) {
	//computes area of circular cap with height h of circle with radius R
	//##MK::numerical stability of the square root?
	//argument h \in [0,2r]
	return ( SQR(r) * acos((r-h)/r) - (r-h)*sqrt(2*r*h-SQR(h)) );
	//##MK::use identity acos(x) = PI/2.0 - asin(x)
}


double mathMethods::triangle_area( double ux, double uy, double vx, double vy ) 
{
	return ( fabs(0.5*(ux*vy - uy*vx)) );
}

double mathMethods::helper_bisec_a_t( double nx, double ny ) {
	return ( SQR(nx) + SQR(ny) );
}

double mathMethods::helper_bisec_b_t( double bx, double by, double nx, double ny, double xc, double yc ) {
	return (2.0*(nx*bx - nx*xc + ny*by - ny*yc));
}

double mathMethods::helper_bisec_c_t( double bx, double by, double xc, double yc, double R ) {
	return ( SQR(bx) - 2.0*xc*bx + SQR(xc) + SQR(by) - 2.0*yc*by + SQR(yc) - SQR(R) );
}

double mathMethods::helper_intersect_a_t( double tx, double ty ) {
	return ( SQR(tx) + SQR(ty) );
}

double mathMethods::helper_intersect_b_t( double xx1, double yy1, double tx, double ty, double xc, double yc ) {
	return ( 2.0*(xx1*tx - xc*tx + yy1*ty - yc*ty) );
}

double mathMethods::helper_intersect_c_t( double xx1, double yy1, double xc, double yc, double R) {
	return( SQR(xx1) + SQR(xc) - 2.0*xc*xx1 + SQR(yy1) + SQR(yc) - 2.0*yy1*yc - SQR(R) );
}


double mathMethods::tri_cir_v0e1( double xa, double ya, double xb, double yb, double x, double y, double r, double xc, double yc ) {
#ifdef DEBUG_TRIANGLE
	std::cout << "v0e1" << std::endl;
	#ifdef DEBUG_EARLY_OUT
		return 0.0;
	#endif
#endif
	//computes the intersection area of a circle at point x, y with radius r with a triangle
	//whose edge (xa,ya) (xb, yb) intersects the circle while all vertices a,b,c lay outside thus intersecting via cap or the negative of a cap
	double aa, bb, cc, tx, ty;
	tx = xb-xa;
	ty = yb-ya;
	aa = this->helper_intersect_a_t( tx, ty );
	bb = this->helper_intersect_b_t( xa, ya, tx, ty, x, y );
	cc = this->helper_intersect_c_t( xa, ya, x, y , r );
	struct pqsolve ar;
	ar = stable_pqsolver( aa, bb, cc );
	if ( ar.nsol == TWO_SOLUTIONS ) { //two applicable solutions expected but only one gives height of cap
		if ( 0.0 <= ar.sol1 && ar.sol1 <= 1.0 && 0.0 <= ar.sol2 && ar.sol2 <= 1.0 ) {
			double bx, by, nx, ny, dist;
			bx = 0.5*((xa+ar.sol1*tx) + (xa+ar.sol2*tx));
			by = 0.5*((ya+ar.sol1*ty) + (ya+ar.sol2*ty));
			//find where connecting line of perpendicular bisectors cut the circle investigate bisector that points through the center
			//we know that no other edge intersects with the circle so either one of them is out the other inside, i.e.
			//we need to find the correct side of circular cap
			nx = ty; //possible because we are interested just a perpendicular to the line tx, ty, MK::which is not normalized and not necessarily the inner or outer normal to the segment!
			ny = -1.0*tx;
			aa = helper_bisec_a_t( nx, ny );
			bb = helper_bisec_b_t( bx, by, nx, ny, x, y );
			cc = helper_bisec_c_t( bx, by, x, y, r );
			struct pqsolve cp;
			cp = stable_pqsolver( aa, bb, cc);
			if ( cp.nsol != TWO_SOLUTIONS ) { std::cerr << "ERROR::v0e1 fewer solutions than expected!" << std::endl; return INTERSECTION_DETECTION_ERROR; }
			//only one can be in the triangle (because only one edge cuts the circle)
			double A = 0.0;
			short consistency = 0;
			if ( point_in_triangle( bx+cp.sol1*nx, by+cp.sol1*ny, xa, ya, xb, yb, xc, yc ) == true ) {
				consistency++;
				dist = sqrt( SQR(cp.sol1*nx) + SQR(cp.sol1*ny) ); //how strongly does the circle stretch into the triangle?
				A = (dist>r) ? (PI*SQR(r)-circular_cap_area(2.0*r-dist, r)) : circular_cap_area( dist, r );
			}
			if ( point_in_triangle( bx+cp.sol2*nx, by+cp.sol2*ny, xa, ya, xb, yb, xc, yc ) == true ) {
				consistency++;
				dist = sqrt( SQR(cp.sol2*nx) + SQR(cp.sol2*ny) );
				A = (dist>r) ? (PI*SQR(r)-circular_cap_area(2.0*r-dist, r)) : circular_cap_area( dist, r );
			}
			if ( consistency == 0 || consistency == 2 ) { std::cerr << "ERROR::Inconsistency v0e1 TWO_SOLUTIONS!" << std::endl; return INTERSECTION_DETECTION_ERROR; }
			return A;
		}
		else { std::cerr << "ERROR::v0e1 no applicable solution for cap computation!" << std::endl; return INTERSECTION_DETECTION_ERROR; }
	}
	else if ( ar.nsol == ONE_SOLUTION ) { //numerical case of circle touching edge of the triangle
		//potentially x,y inside the triangle
		if ( point_in_triangle( x, y, xa, ya, xb, yb, xc, yc ) == true ) //circle touching triangle edge from the inside
			return PI*SQR(r);
		else //circle laying outside but touching the edge
			return 0.0;
	}
	else { std::cerr << "ERROR::v0e1 no solutions found cap!" << std::endl; return INTERSECTION_DETECTION_ERROR; }
	
}


double mathMethods::tri_cir_v0e2(  double xa, double ya, double xb, double yb, double xc, double yc, double x, double y, double r, bool ab, bool bc, bool ca ) {
#ifdef DEBUG_TRIANGLE
	std::cout << "v0e2" << std::endl;
	#ifdef DEBUG_EARLY_OUT
		return 0.0;
	#endif
#endif
	double res = PI*SQR(r);

	//which caps to substract, only those which are cut, function assumes that exactly two of the arguments ab, bc, ca are true the other false
	double x1,y1,x2,y2, tx, ty, a, b, c, bx, by, nx, ny, dist;
	double outernormal[2];
	//short consistency = 0;
	struct pqsolve ar, cp;
	if ( ab == true ) { //all three ifs have the same structure
		x1 = xa; //we only change the point association
		y1 = ya;
		x2 = xb;
		y2 = yb;

		tx = x2-x1;
		ty = y2-y1;
		a = helper_intersect_a_t( tx, ty );
		b = helper_intersect_b_t( x1, y1, tx, ty, x, y );
		c = helper_intersect_c_t( x1, y1, x, y, r );
		ar.nsol = NO_SOLUTION; ar.sol1 = INVALID_SOL; ar.sol2 = INVALID_SOL;
		ar = stable_pqsolver( a, b, c );
		if ( ar.nsol == TWO_SOLUTIONS ) {
			if ( 0.0 <= ar.sol1 && ar.sol1 <= 1.0 && 0.0 <= ar.sol2 && ar.sol2 <= 1.0 ) { //applicable solution?
				bx = 0.5*((x1 + ar.sol1*tx) + (x1 + ar.sol2*tx));
				by = 0.5*((y1 + ar.sol1*ty) + (y1 + ar.sol2*ty));
				outer_triedge_normal( xa,ya,xb,yb,xc,yc, bx, by, tx, ty, outernormal );
				nx = outernormal[0];
				ny = outernormal[1];
				a = helper_bisec_a_t( nx, ny );
				b = helper_bisec_b_t( bx, by, nx, ny, x, y );
				c = helper_bisec_c_t( bx, by, x, y, r );
				cp.nsol = NO_SOLUTION; cp.sol1 = INVALID_SOL; cp.sol2 = INVALID_SOL;
				cp = stable_pqsolver( a, b, c);
				//the cap height to cutoff is the distance from bx, by to the circle intersection point b+ cp.sol * n but 
				//only the positive solution because n is an outernormal!
				if ( cp.nsol != TWO_SOLUTIONS ) { std::cerr << "ERROR::v0e2 ab fewer solutions than expected!" << std::endl; return INTERSECTION_DETECTION_ERROR; }
				//so there are two solutions, find the positive one
				if ( cp.sol1 >= 0.0 )
					dist = sqrt(SQR(cp.sol1*nx)+SQR(cp.sol1*ny));
				else if ( cp.sol2 >= 0.0 )
					dist = sqrt(SQR(cp.sol2*nx)+SQR(cp.sol2*ny));
				else { std::cerr << "ERROR::inconsistency v0e2 1.detected in computation of cap area!" << std::endl; return INTERSECTION_DETECTION_ERROR; }
				//MK::the calling of the circular cap area computation is stable and possible for height arguments in [0,2r]!
				res = res - circular_cap_area( dist, r );
				/*consistency = 0;
				dist = std::numeric_limits<double>:: max();
				if ( point_in_triangle( bx+cp.sol1*nx, by+cp.sol1*ny, xa, ya, xb, yb, xc, yc ) == false ) { //find the point outside the triangle to compute the cap height
					consistency++;
					if ( sqrt(SQR(cp.sol1*nx) + SQR(cp.sol1*ny)) < dist )
						dist = sqrt( SQR(cp.sol1*nx) + SQR(cp.sol1*ny) );
				}
				if ( point_in_triangle( bx+cp.sol2*nx, by+cp.sol2*ny, xa, ya, xb, yb, xc, yc ) == false ) { //MK::need to distinguish cases because nx,ny is neither an outer nor a normal vector to the edge!
					consistency++;
					if ( sqrt(SQR(cp.sol2*nx) + SQR(cp.sol2*ny)) < dist )
						dist = sqrt( SQR(cp.sol2*nx) + SQR(cp.sol2*ny) );
				}
				if ( consistency == 0 ) { std::cerr << "ERROR::inconsistency v0e2 1.detected in computation of cap area!" << std::endl; return 0.0; }
				res = res - circular_cap_area( dist, r );*/
			}
			else { std::cerr << "ERROR::tri_cir_v0e2 ab no applicable solution bisector point found!" << std::endl; }
		} 
		else if ( ar.nsol == NO_SOLUTION ) { std::cerr << "ERROR::tri_cir_v0e2 ab does not find at all a solution for computing the outside cap!" << std::endl; return INTERSECTION_DETECTION_ERROR; } //##MK::add output
		else {} //nothing to do ar.nsol == ONE_SOLUTION ) { //numerical touching of circle from the inside at the edge --> degeneracy of the cap, do nothing because we assume area negligible
	}
	if ( bc == true ) {
		x1 = xb;
		y1 = yb;
		x2 = xc;
		y2 = yc;

		tx = x2-x1;
		ty = y2-y1;
		a = helper_intersect_a_t( tx, ty );
		b = helper_intersect_b_t( x1, y1, tx, ty, x, y );
		c = helper_intersect_c_t( x1, y1, x, y, r );
		ar.nsol = NO_SOLUTION; ar.sol1 = INVALID_SOL; ar.sol2 = INVALID_SOL;
		ar = stable_pqsolver( a, b, c );
		if ( ar.nsol == TWO_SOLUTIONS ) {
			if ( 0.0 <= ar.sol1 && ar.sol1 <= 1.0 && 0.0 <= ar.sol2 && ar.sol2 <= 1.0 ) { //applicable solution?
				bx = 0.5*((x1 + ar.sol1*tx) + (x1 + ar.sol2*tx));
				by = 0.5*((y1 + ar.sol1*ty) + (y1 + ar.sol2*ty));
				outer_triedge_normal( xa,ya,xb,yb,xc,yc, bx, by, tx, ty, outernormal );
				nx = outernormal[0];
				ny = outernormal[1];
				a = helper_bisec_a_t( nx, ny );
				b = helper_bisec_b_t( bx, by, nx, ny, x, y );
				c = helper_bisec_c_t( bx, by, x, y, r );
				cp.nsol = NO_SOLUTION; cp.sol1 = INVALID_SOL; cp.sol2 = INVALID_SOL;
				cp = stable_pqsolver( a, b, c);
				if ( cp.nsol != TWO_SOLUTIONS ) { std::cerr << "ERROR::v0e2 ab fewer solutions than expected!" << std::endl; return INTERSECTION_DETECTION_ERROR; }
				/*consistency = 0;
				dist = std::numeric_limits<double>:: max();
				if ( point_in_triangle( bx+cp.sol1*nx, by+cp.sol1*ny, xa, ya, xb, yb, xc, yc ) == false ) { //find the point outside the triangle to compute the cap height
					consistency++;
					if ( sqrt( SQR(cp.sol1*nx) + SQR(cp.sol1*ny) ) < dist )
						dist = sqrt( SQR(cp.sol1*nx) + SQR(cp.sol1*ny) );
				}
				if ( point_in_triangle( bx+cp.sol2*nx, by+cp.sol2*ny, xa, ya, xb, yb, xc, yc ) == false ) { //MK::need to distinguish cases because nx,ny is neither an outer nor a normal vector to the edge!
					consistency++;
					if ( sqrt( SQR(cp.sol2*nx) + SQR(cp.sol2*ny) ) < dist )
						dist = sqrt( SQR(cp.sol2*nx) + SQR(cp.sol2*ny) );
				}
				if ( consistency == 0 ) { std::cerr << "ERROR::inconsistency v0e2 2. detected in computation of cap area!" << std::endl; return 0.0; }
				res = res - circular_cap_area( dist, r );*/
				if ( cp.sol1 >= 0.0 )
					dist = sqrt(SQR(cp.sol1*nx)+SQR(cp.sol1*ny));
				else if ( cp.sol2 >= 0.0 )
					dist = sqrt(SQR(cp.sol2*nx)+SQR(cp.sol2*ny));
				else { std::cerr << "ERROR::inconsistency v0e2 2. detected in computation of cap area!" << std::endl; return INTERSECTION_DETECTION_ERROR; }
				res = res - circular_cap_area( dist, r );
			}
			else { std::cerr << "ERROR::tri_cir_v0e2 ab no applicable solution bisector point found!" << std::endl; }
		} 
		else if ( ar.nsol == NO_SOLUTION ) { std::cerr << "ERROR::tri_cir_v0e2 ab does not find at all a solution for computing the outside cap!" << std::endl; return INTERSECTION_DETECTION_ERROR; } //##MK::add output
		else {} //nothing to do ar.nsol == ONE_SOLUTION ) { //numerical touching of circle from the inside at the edge --> degeneracy of the cap, do nothing because we assume area negligible
	}
	if ( ca == true ) {
		x1 = xc;
		y1 = yc;
		x2 = xa;
		y2 = ya;

		tx = x2-x1;
		ty = y2-y1;
		a = helper_intersect_a_t( tx, ty );
		b = helper_intersect_b_t( x1, y1, tx, ty, x, y );
		c = helper_intersect_c_t( x1, y1, x, y, r );
		ar.nsol = NO_SOLUTION; ar.sol1 = INVALID_SOL; ar.sol2 = INVALID_SOL;
		ar = stable_pqsolver( a, b, c );
		if ( ar.nsol == TWO_SOLUTIONS ) {
			if ( 0.0 <= ar.sol1 && ar.sol1 <= 1.0 && 0.0 <= ar.sol2 && ar.sol2 <= 1.0 ) { //applicable solution?
				bx = 0.5*((x1 + ar.sol1*tx) + (x1 + ar.sol2*tx));
				by = 0.5*((y1 + ar.sol1*ty) + (y1 + ar.sol2*ty));
				outer_triedge_normal( xa,ya,xb,yb,xc,yc, bx, by, tx, ty, outernormal );
				nx = outernormal[0];
				ny = outernormal[1];
				a = helper_bisec_a_t( nx, ny );
				b = helper_bisec_b_t( bx, by, nx, ny, x, y );
				c = helper_bisec_c_t( bx, by, x, y, r );
				cp.nsol = NO_SOLUTION; cp.sol1 = INVALID_SOL; cp.sol2 = INVALID_SOL;
				cp = stable_pqsolver( a, b, c);
				if ( cp.nsol != TWO_SOLUTIONS ) { std::cerr << "ERROR::v0e2 ab fewer solutions than expected!" << std::endl; return INTERSECTION_DETECTION_ERROR; }
				/*consistency = 0;
				dist = std::numeric_limits<double>:: max();
				if ( point_in_triangle( bx+cp.sol1*nx, by+cp.sol1*ny, xa, ya, xb, yb, xc, yc ) == false ) { //find the point outside the triangle to compute the cap height
					consistency++;
					if ( sqrt( SQR(cp.sol1*nx) + SQR(cp.sol1*ny) ) < dist )
						dist = sqrt( SQR(cp.sol1*nx) + SQR(cp.sol1*ny) );
				}
				if ( point_in_triangle( bx+cp.sol2*nx, by+cp.sol2*ny, xa, ya, xb, yb, xc, yc ) == false ) { //MK::need to distinguish cases because nx,ny is neither an outer nor a normal vector to the edge!
					consistency++;
					if ( sqrt( SQR(cp.sol2*nx) + SQR(cp.sol2*ny) ) < dist )
						dist = sqrt( SQR(cp.sol2*nx) + SQR(cp.sol2*ny) );
				}
				if ( consistency == 0 ) { std::cerr << "ERROR::inconsistency v0e2 3. detected in computation of cap area!" << std::endl; return 0.0; }
				res = res - circular_cap_area( dist, r );*/
				if ( cp.sol1 >= 0.0 )
					dist = sqrt(SQR(cp.sol1*nx)+SQR(cp.sol1*ny));
				else if ( cp.sol2 >= 0.0 )
					dist = sqrt(SQR(cp.sol2*nx)+SQR(cp.sol2*ny));
				else { std::cerr << "ERROR::inconsistency v0e2 2. detected in computation of cap area!" << std::endl; return INTERSECTION_DETECTION_ERROR; }
				res = res - circular_cap_area( dist, r );
			}
			else { std::cerr << "ERROR::tri_cir_v0e2 ab no applicable solution bisector point found!" << std::endl; }
		} 
		else if ( ar.nsol == NO_SOLUTION ) { std::cerr << "ERROR::tri_cir_v0e2 ab does not find at all a solution for computing the outside cap!" << std::endl; return INTERSECTION_DETECTION_ERROR; } //##MK::add output
		else {} //nothing to do ar.nsol == ONE_SOLUTION ) { //numerical touching of circle from the inside at the edge --> degeneracy of the cap, do nothing because we assume area negligible
	}

	return res;
}


double mathMethods::tri_cir_v0e3( double xa, double ya, double xb, double yb, double xc, double yc, double x, double y, double r ) {
#ifdef DEBUG_TRIANGLE
	std::cout << "v0e3" << std::endl;
	#ifdef DEBUG_EARLY_OUT
		return 0.0;
	#endif
#endif
	//computes the intersection area of a circle at point x, y with radius r with a triangle
	//whose edges (xa,ya), (xb, yb), (xc, yc) are outside the triangle, the edges of cut off three circular caps from the circle
	//strategy: inclusion area is circle area minus area of these three caps
	double res = PI*SQR(r);

	double a, b, c;
	double x1,y1, x2, y2, tx, ty;
	double bx, by, nx, ny, dist;
	short consistency;
	struct pqsolve ar, cp;
	double outernormal[2];
	x1 = xa;
	y1 = ya;
	x2 = xb;
	y2 = yb;
//--->substract first cap from now on reutilizable code snippets

	tx = x2-x1;
	ty = y2-y1;
	a = helper_intersect_a_t( tx, ty );
	b = helper_intersect_b_t( x1, y1, tx, ty, x, y );
	c = helper_intersect_c_t( x1, y1, x, y, r );
	ar.nsol = NO_SOLUTION; ar.sol1 = INVALID_SOL; ar.sol2 = INVALID_SOL;
	ar = stable_pqsolver( a, b, c );
	if ( ar.nsol == TWO_SOLUTIONS ) {
		if ( 0.0 <= ar.sol1 && ar.sol1 <= 1.0 && 0.0 <= ar.sol2 && ar.sol2 <= 1.0 ) { //applicable bisector point found
			bx = 0.5*((x1 + ar.sol1*tx) + (x1 + ar.sol2*tx));
			by = 0.5*((y1 + ar.sol1*ty) + (y1 + ar.sol2*ty));
			outer_triedge_normal(xa,ya,xb,yb,xc,yc,bx,by,tx,ty, outernormal );
			nx = outernormal[0]; //ty;
			ny = outernormal[1]; //-1.0*tx;
			a = helper_bisec_a_t( nx, ny );
			b = helper_bisec_b_t( bx, by, nx, ny, x, y );
			c = helper_bisec_c_t( bx, by, x, y, r );
			cp.nsol = NO_SOLUTION; cp.sol1 = INVALID_SOL; cp.sol2 = INVALID_SOL;
			cp = stable_pqsolver( a, b, c);
			if ( cp.nsol != TWO_SOLUTIONS ) { std::cerr << "ERROR::v0e3 1. fewer solutions than expected in a  bisector search!" << std::endl; return INTERSECTION_DETECTION_ERROR; }
			//compute overlapping cap area knowing there are two solutions for cp
			consistency = 0;
			if ( cp.sol1 > 0.0 ) {
				consistency++;
				dist = sqrt(SQR(cp.sol1*nx)+SQR(cp.sol1*ny));
			}
			if ( cp.sol2 > 0.0 ) {
				consistency++;
				dist = sqrt(SQR(cp.sol2*nx)+SQR(cp.sol2*ny));
			}
			if ( consistency == 0 || consistency == 2 ) { std::cerr << "ERROR::v0e3 1. fewer solutions than expected in a bisector search, no positive solutions!" << std::endl; return INTERSECTION_DETECTION_ERROR; }
			res = res - circular_cap_area( dist, r );
			/*dist = std::numeric_limits<double>:: max();
			if ( point_in_triangle( bx+cp.sol1*nx, by+cp.sol1*ny, xa, ya, xb, yb, xc, yc ) == false ) { //find the point outside the triangle to compute the cap height
				consistency++;
				if ( sqrt( SQR(cp.sol1*nx) + SQR(cp.sol1*ny) ) < dist ) 
					dist = sqrt( SQR(cp.sol1*nx) + SQR(cp.sol1*ny) );
			}
			if ( point_in_triangle( bx+cp.sol2*nx, by+cp.sol2*ny, xa, ya, xb, yb, xc, yc ) == false ) { //MK::need to distinguish cases because nx,ny is neither an outer nor a normal vector to the edge!
				consistency++;
				if ( sqrt( SQR(cp.sol2*nx) + SQR(cp.sol2*ny) ) < dist ) 
					dist = sqrt( SQR(cp.sol2*nx) + SQR(cp.sol2*ny) );
			}
			if ( consistency == 0 ) { std::cerr << "ERROR::inconsistency v0e3 1. detected in computation of cap area!" << std::endl; return 0.0; }
			res = res - circular_cap_area( dist, r );*/
		}
		else { std::cerr << "ERROR::tri_cir_v0e3 1. no applicable bisector point found!" << std::endl; return INTERSECTION_DETECTION_ERROR; }
	} 
	else if ( ar.nsol == NO_SOLUTION ) { std::cerr << "ERROR::tri_cir_v0e3 1. does not find two solutions for computing the outside cap!" << std::endl; return INTERSECTION_DETECTION_ERROR; } //##MK::add output
	else {} //nothing to do ar.nsol == ONE_SOLUTION ) { //numerical touching of circle from the inside at the edge --> degeneracy of the cap, do nothing

//--->substract second cap
	x1 = xb;
	y1 = yb;
	x2 = xc;
	y2 = yc;

	tx = x2-x1;
	ty = y2-y1;
	a = helper_intersect_a_t( tx, ty );
	b = helper_intersect_b_t( x1, y1, tx, ty, x, y );
	c = helper_intersect_c_t( x1, y1, x, y, r );
	ar.nsol = NO_SOLUTION; ar.sol1 = INVALID_SOL; ar.sol2 = INVALID_SOL;
	ar = stable_pqsolver( a, b, c );
	if ( ar.nsol == TWO_SOLUTIONS ) {
		if ( 0.0 <= ar.sol1 && ar.sol1 <= 1.0 && 0.0 <= ar.sol2 && ar.sol2 <= 1.0 ) { //applicable bisector point found
			bx = 0.5*((x1 + ar.sol1*tx) + (x1 + ar.sol2*tx));
			by = 0.5*((y1 + ar.sol1*ty) + (y1 + ar.sol2*ty));
			outer_triedge_normal(xa,ya,xb,yb,xc,yc,bx,by,tx,ty, outernormal );
			nx = outernormal[0]; //ty;
			ny = outernormal[1]; //-1.0*tx;
			a = helper_bisec_a_t( nx, ny );
			b = helper_bisec_b_t( bx, by, nx, ny, x, y );
			c = helper_bisec_c_t( bx, by, x, y, r );
			cp.nsol = NO_SOLUTION; cp.sol1 = INVALID_SOL; cp.sol2 = INVALID_SOL;
			cp = stable_pqsolver( a, b, c);
			if ( cp.nsol != TWO_SOLUTIONS ) { std::cerr << "ERROR::v0e3 2. fewer solutions than expected in a  bisector search!" << std::endl; return INTERSECTION_DETECTION_ERROR; }
			//compute overlapping cap area knowing there are two solutions for cp
			consistency = 0;
			if ( cp.sol1 > 0.0 ) {
				consistency++;
				dist = sqrt(SQR(cp.sol1*nx)+SQR(cp.sol1*ny));
			}
			if ( cp.sol2 > 0.0 ) {
				consistency++;
				dist = sqrt(SQR(cp.sol2*nx)+SQR(cp.sol2*ny));
			}
			if ( consistency == 0 || consistency == 2 ) { std::cerr << "ERROR::v0e3 2. fewer solutions than expected in a bisector search, no positive solutions!" << std::endl; return INTERSECTION_DETECTION_ERROR; }
			res = res - circular_cap_area( dist, r );
		}
		else { std::cerr << "ERROR::tri_cir_v0e3 2. no applicable bisector point found!" << std::endl; return INTERSECTION_DETECTION_ERROR; }
	} 
	else if ( ar.nsol == NO_SOLUTION ) { std::cerr << "ERROR::tri_cir_v0e3 2. does not find two solutions for computing the outside cap!" << std::endl; return INTERSECTION_DETECTION_ERROR; } //##MK::add output
	else {} //nothing to do ar.nsol == ONE_SOLUTION ) { //numerical touching of circle from the inside at the edge --> degeneracy of the cap, do nothing

//--->substract last cap
	x1 = xc;
	y1 = yc;
	x2 = xa;
	y2 = ya;

	tx = x2-x1;
	ty = y2-y1;
	a = helper_intersect_a_t( tx, ty );
	b = helper_intersect_b_t( x1, y1, tx, ty, x, y );
	c = helper_intersect_c_t( x1, y1, x, y, r );
	ar.nsol = NO_SOLUTION; ar.sol1 = INVALID_SOL; ar.sol2 = INVALID_SOL;
	ar = stable_pqsolver( a, b, c );
	if ( ar.nsol == TWO_SOLUTIONS ) {
		if ( 0.0 <= ar.sol1 && ar.sol1 <= 1.0 && 0.0 <= ar.sol2 && ar.sol2 <= 1.0 ) { //applicable bisector point found
			bx = 0.5*((x1 + ar.sol1*tx) + (x1 + ar.sol2*tx));
			by = 0.5*((y1 + ar.sol1*ty) + (y1 + ar.sol2*ty));
			outer_triedge_normal(xa,ya,xb,yb,xc,yc,bx,by,tx,ty, outernormal );
			nx = outernormal[0]; //ty;
			ny = outernormal[1]; //-1.0*tx;
			a = helper_bisec_a_t( nx, ny );
			b = helper_bisec_b_t( bx, by, nx, ny, x, y );
			c = helper_bisec_c_t( bx, by, x, y, r );
			cp.nsol = NO_SOLUTION; cp.sol1 = INVALID_SOL; cp.sol2 = INVALID_SOL;
			cp = stable_pqsolver( a, b, c);
			if ( cp.nsol != TWO_SOLUTIONS ) { std::cerr << "ERROR::v0e3 3. fewer solutions than expected in a  bisector search!" << std::endl; return INTERSECTION_DETECTION_ERROR; }
			//compute overlapping cap area knowing there are two solutions for cp
			consistency = 0;
			if ( cp.sol1 > 0.0 ) {
				consistency++;
				dist = sqrt(SQR(cp.sol1*nx)+SQR(cp.sol1*ny));
			}
			if ( cp.sol2 > 0.0 ) {
				consistency++;
				dist = sqrt(SQR(cp.sol2*nx)+SQR(cp.sol2*ny));
			}
			if ( consistency == 0 || consistency == 2 ) { std::cerr << "ERROR::v0e3 3. fewer solutions than expected in a bisector search, no positive solutions!" << std::endl; return INTERSECTION_DETECTION_ERROR; }
			res = res - circular_cap_area( dist, r );
		}
		else { std::cerr << "ERROR::tri_cir_v0e3 3. no applicable bisector point found!" << std::endl; return INTERSECTION_DETECTION_ERROR; }
	} 
	else if ( ar.nsol == NO_SOLUTION ) { std::cerr << "ERROR::tri_cir_v0e3 3. does not find two solutions for computing the outside cap!" << std::endl; return INTERSECTION_DETECTION_ERROR; } //##MK::add output
	else {} //nothing to do ar.nsol == ONE_SOLUTION ) { //numerical touching of circle from the inside at the edge --> degeneracy of the cap, do nothing

	return res;
}


double mathMethods::tri_cir_v1e2( double xa, double ya, double xb, double yb, double xc, double yc, double x, double y, double r, bool a, bool b, bool c, bool ab, bool bc, bool ca ) {
#ifdef DEBUG_TRIANGLE
	std::cout << "v1e2" << std::endl;
	#ifdef DEBUG_EARLY_OUT
		return 0.0;
	#endif
#endif
	//computes intersection area of arbitrary nondegenerated triangle with a circle with one of its vertex in[2] in the circle
	//(the one for which a,b,c == true) and two adjoning edges cutting the circle while the opposite laying remaining edge is outside the circlee
	//partition triangle into a) old triangle truncated by new edge (connecting the edge intersection points) + b) circular cap area additional intersection with the triangle

	//which vertex is inside?
	double in[2], e1[6], e2[6]; //e1 inclosed point other point is outside
	if ( a == true ) {
		in[0] = xa;		in[1] = ya;
		e1[0] = xa;		e1[1] = ya;		e1[2] = xb;		e1[3] = yb;		e1[4] = e1[2]-e1[0];		e1[5] = e1[3]-e1[1];
		e2[0] = xa;		e2[1] = ya;		e2[2] = xc;		e2[3] = yc;		e2[4] = e2[2]-e1[0];		e2[5] = e2[3]-e2[1];
	} else if ( b == true ) {
		in[0] = xb;		in[1] = yb;
		e1[0] = xb;		e1[1] = yb;		e1[2] = xa;		e1[3] = ya;		e1[4] = e1[2]-e1[0];		e1[5] = e1[3]-e1[1];
		e2[0] = xb;		e2[1] = yb;		e2[2] = xc;		e2[3] = yc;		e2[4] = e2[2]-e1[0];		e2[5] = e2[3]-e2[1];
	} else if ( c == true ) {
		in[0] = xc;		in[1] = yc;
		e1[0] = xc;		e1[1] = yc;		e1[2] = xa;		e1[3] = ya;		e1[4] = e1[2]-e1[0];		e1[5] = e1[3]-e1[1];
		e2[0] = xc;		e2[1] = yc;		e2[2] = xb;		e2[3] = yb;		e2[4] = e2[2]-e1[0];		e2[5] = e2[3]-e2[1];
	} else { std::cerr << "ERROR:: in tri_cir_v1e2 without any vertex inside!" << std::endl; return INTERSECTION_DETECTION_ERROR; }

	//which edges cut the triangle?
	//find intersection points
	double aa, bb, cc;
	double ip1[2];
	aa = helper_intersect_a_t( e1[4], e1[5] );
	bb = helper_intersect_b_t( e1[0], e1[1], e1[4], e1[5], x, y );
	cc = helper_intersect_c_t( e1[0], e1[1], x, y, r );
	struct pqsolve i1;
	i1 = stable_pqsolver( aa, bb, cc ); //only solutions "on the edge segment" in [0,1] accepted
	if ( i1.nsol == ONE_SOLUTION ) {
		if ( 0.0 <= i1.sol1 && i1.sol1 <= 1.0 ) { //first solution on segment?
			ip1[0] = e1[0] + i1.sol1*e1[4];
			ip1[1] = e1[1] + i1.sol1*e1[5];
		}
		else { std::cerr << "ERROR::Inconsistency v1e2 i1 ONESOLUTION!" << std::endl; return INTERSECTION_DETECTION_ERROR; }
	}
	else if ( i1.nsol == TWO_SOLUTIONS ) { //##MK::should not be possible as one of the two triangle vertices is inside the other outside
		if ( 0.0 <= i1.sol1 && i1.sol1 <= 1.0 ) {
			ip1[0] = e1[0] + i1.sol1*e1[4];
			ip1[1] = e1[1] + i1.sol1*e1[5];
		}
		else if ( 0.0 <= i1.sol2 && i1.sol2 <= 1.0 ) {
			ip1[0] = e1[0] + i1.sol2*e1[4];
			ip1[1] = e1[1] + i1.sol2*e1[5];
		}
		else { std::cerr << "ERROR::Inconsistency v1e2 i1 TWOSOLUTIONS!" << std::endl; return INTERSECTION_DETECTION_ERROR; }
	}
	else { std::cerr << "ERROR::No solution v1e2 i1 NOSOLUTION!" << std::endl; return INTERSECTION_DETECTION_ERROR; }

	double ip2[2];
	aa = helper_intersect_a_t( e2[4], e2[5] );
	bb = helper_intersect_b_t( e2[0], e2[1], e2[4], e2[5], x, y );
	cc = helper_intersect_c_t( e2[0], e2[1], x, y, r );
	struct pqsolve i2;
	i2 = stable_pqsolver( aa, bb, cc ); //only solution "on the edge segment" in [0,1] accepted
	if ( i2.nsol == ONE_SOLUTION ) {
		if ( 0.0 <= i2.sol1 && i2.sol1 <= 1.0 ) { //the only solution on the segment?
			ip2[0] = e2[0] + i2.sol1*e2[4];
			ip2[1] = e2[1] + i2.sol1*e2[5];
		}
		else { std::cerr << "ERROR::Inconsistency v1e2 i2 ONESOLUTION!" << std::endl; return INTERSECTION_DETECTION_ERROR; }
	}
	else if ( i2.nsol == TWO_SOLUTIONS ) {
		if ( 0.0 <= i2.sol1 && i2.sol1 <= 1.0 ) {
			ip2[0] = e2[0] + i2.sol1*e2[4];
			ip2[1] = e2[1] + i2.sol1*e2[5];
		}
		else if ( 0.0 <= i2.sol2 && i2.sol2 <= 1.0 ) { //no, okay but the second?
			ip2[0] = e2[0] + i2.sol2*e2[4];
			ip2[1] = e2[1] + i2.sol2*e2[5];
		}
		else { std::cerr << "ERROR::Inconsistency v1e2 i2 TWOSOLUTIONS!" << std::endl; return INTERSECTION_DETECTION_ERROR; }
	}
	else { std::cerr << "ERROR::No solution v1e2 i2 NOSOLUTION!" << std::endl; return INTERSECTION_DETECTION_ERROR; }

	//knowing now the intersection points ix1,iy1 and ix2,iy2 we can compute the triangle area spanned by
	//in[0],in[1] , ix1,iy1 , ix2,iy2
	double A = triangle_area( ip1[0]-in[0], ip1[1]-in[1], ip2[0]-in[0], ip2[1]-in[1] );

	//for very flat almost numerically degenerated triangles a heuristics is necessary to avoid to strong numerical inaccuracies
	//when attempting to compute the circular cap area with the cap height (dist) approaching double precision
#ifdef HEURISTIC_HANDLING_OF_SMALL_NUMBERS
	if ( A < 0.0 ) {
		std::cerr << "ERROR::Inscribed triangle in v1e2 has negative area!" << std::endl; return INTERSECTION_DETECTION_ERROR;
	}
	if ( A <= MINIMUM_CRITICAL_INTERSECTION_AREA ) {
//		std::cerr << "WARNING::v1e2 expects very small cap height therefore returning less accurate value!" << std::endl;
		return A;
	}
#endif

	//it remains to get the area of the circular cap which is cut of from the circle by the segment line i1i2
	double bx = 0.5*(ip1[0]+ip2[0]);
	double by = 0.5*(ip1[1]+ip2[1]);
	double tx = ip2[0]-ip1[0];
	double ty = ip2[1]-ip1[1];
	double outernormal[2];
	outer_triedge_normal(in[0],in[1],ip1[0],ip1[1],ip2[0],ip2[1],bx,by,tx,ty,outernormal);
	double nx = outernormal[0]; //ty; //ip2[1]-ip1[1];
	double ny = outernormal[1]; //tx; //-1.0*(ip2[0]-ip1[0]);

	//MK::the included vertex in[0], in[1] is not necessarily located at the center of the circle!
	struct pqsolve bi;
	aa = helper_bisec_a_t( nx, ny );
	bb = helper_bisec_b_t( bx, by, nx, ny, x, y );
	cc = helper_bisec_c_t( bx, by, x, y, r );
	bi = stable_pqsolver( aa, bb, cc);
	if ( bi.nsol != TWO_SOLUTIONS ) { std::cerr << "ERROR::Inconsistency in v1e2 at cap calculation!" << std::endl; return INTERSECTION_DETECTION_ERROR; }
	short consistency = 0;
	double dist;
	if ( bi.sol1 > 0.0 ) {
		consistency++;
		dist = sqrt(SQR(bi.sol1*nx)+SQR(bi.sol1*ny));
	}
	if ( bi.sol2 > 0.0 ) {
		consistency++;
		dist = sqrt(SQR(bi.sol2*nx)+SQR(bi.sol2*ny));
	}
	if ( consistency == 0 || consistency == 2 ) { std::cerr << "ERROR::Inconsistency in v1e2 finding positive cap solutions!" << std::endl; return INTERSECTION_DETECTION_ERROR; }
	A = A + circular_cap_area( dist, r );
	return A;
	/*double dist;
	if ( bi.nsol == TWO_SOLUTIONS ) {
		short consistency = 0;
		if ( point_in_triangle( bx+bi.sol1*nx, by+bi.sol1*ny, xa,ya,xb,yb,xc,yc ) == true ) {
			consistency++;
			dist = sqrt( SQR(bi.sol1*nx) + SQR(bi.sol1*ny) );
			if ( dist>r )	A = A + (PI*SQR(r) - circular_cap_area( 2.0*r-dist, r ));
			else			A = A + circular_cap_area( dist, r );
		}
		if ( point_in_triangle( bx+bi.sol2*nx, by+bi.sol2*ny, xa,ya,xb,yb,xc,yc ) == true ) {
			consistency++;
			dist = sqrt( SQR(bi.sol2*nx) + SQR(bi.sol2*ny) );
			if ( dist>r)	A = A + (PI*SQR(r) - circular_cap_area( 2.0*r-dist, r ));
			else			A = A + circular_cap_area( dist,r );
		}
		if ( consistency == 0 || consistency == 2 ) { std::cerr << "ERROR::Inconsistency v1e2 bi two solutions!" << std::endl; return 0.0; }
		return A;
	}
	else { std::cerr << "WARNING::Hardly circle/triangle overlap v1e2 bisector calculation!" << std::endl; return 0.0; }
	//not return yet, so any of the impossible configration, complain!
	std::cerr << "ERROR::No return value for v1e2 bisector calculation!" << std::endl;
	return 0.0;*/
}


double mathMethods::tri_cir_v1e3( double xa, double ya, double xb, double yb, double xc, double yc, double x, double y, double r, bool a, bool b, bool c, bool ab, bool bc, bool ca ) {
#ifdef DEBUG_TRIANGLE
	std::cout << "v1e3" << std::endl;
	#ifdef DEBUG_EARLY_OUT
		return 0.0;
	#endif
#endif
	//computes intersection area of arbitrary nondegenerated triangle with a circle with one of its vertices in the circle
	//(the one for which a,b,c == true) and three edges cutting the circle, we can obtain the area by partitioning into
	//a) triangle formed by the included vertex and the intersection points on the two adjoining edges +
	//b) a quad polygon formed by these intersection points and the two other insections of the remaining opposite edge with the circle +
	//c,d) the two circular caps above or below the connection lines of these vertices connecting the opposite edge with the partitioning edge of triangle a)

	double in[2], adj[4]; //structure edge1 opposite point x,y | edge 2 opposite point x,y
	if ( a == true ) {
		in[0] = xa;		in[1] = ya;
		adj[0] = xb;	adj[1] = yb; //order of triangle vertex storage is not relevant as there are only two vertices that can be opposite to another in a triangle
		adj[2] = xc;	adj[3] = yc;
	} else if ( b == true ) {
		in[0] = xb;		in[1] = yb;
		adj[0] = xa;	adj[1] = ya;
		adj[2] = xc;	adj[3] = yc;
	} else if ( c == true ) {
		in[0] = xc;		in[1] = yc;
		adj[0] = xa;	adj[1] = ya;
		adj[2] = xb;	adj[3] = yb;
	} else { std::cerr << "ERROR::v1e3 No vertex inside!" << endl; return INTERSECTION_DETECTION_ERROR; }

	//get intersection points of adjointing edges as i1 and i2
	bool found = false;
	double aa, bb, cc;
	double t1[2], ip1[2];
	t1[0] = adj[0]-in[0];
	t1[1] = adj[1]-in[1];
	aa = helper_intersect_a_t( t1[0], t1[1] );
	bb = helper_intersect_b_t( in[0], in[1], t1[0], t1[1], x, y );
	cc = helper_intersect_c_t( in[0], in[1], x, y, r );
	struct pqsolve i1;
	i1 = stable_pqsolver( aa, bb, cc ); //only solutions "on the edge segment" in [0,1] accepted
	if ( i1.nsol == ONE_SOLUTION ) {
		if ( 0.0 <= i1.sol1 && i1.sol1 <= 1.0 ) { //first solution on segment?
			ip1[0] = in[0] + i1.sol1*t1[0];
			ip1[1] = in[1] + i1.sol1*t1[1];
		}
		else { std::cerr << "ERROR::Inconsistency v1e3 i1 ONESOLUTION!" << std::endl; return INTERSECTION_DETECTION_ERROR; }
	}
	else if ( i1.nsol == TWO_SOLUTIONS ) { //##MK::should not be possible as one vertex is inside the endpoint outside
		found = false;
		if ( 0.0 <= i1.sol1 && i1.sol1 <= 1.0 ) { //first solution on segment?
			found = true;
			ip1[0] = in[0] + i1.sol1*t1[0];
			ip1[1] = in[1] + i1.sol1*t1[1];
		}
		//was else if
		if ( 0.0 <= i1.sol2 && i1.sol2 <= 1.0 ) { //MK::because the roots are ordered ascendingly if both roots evaluate to intersection points this will be the farther apart from the point in
			found = true;
			ip1[0] = in[0] + i1.sol2*t1[0];
			ip1[1] = in[1] + i1.sol2*t1[1];
		}
		if ( found == false ) { std::cerr << "ERROR::Inconsistency v1e3 i1 TWOSOLUTIONS!" << std::endl; return INTERSECTION_DETECTION_ERROR; }
	}
	else { std::cerr << "ERROR::No solution v1e3 i1 NOSOLUTION!" << std::endl; return INTERSECTION_DETECTION_ERROR; }

	double t2[2], ip2[2];
	t2[0] = adj[2]-in[0];
	t2[1] = adj[3]-in[1];
	aa = this->helper_intersect_a_t( t2[0], t2[1] );
	bb = this->helper_intersect_b_t( in[0], in[1], t2[0], t2[1], x, y );
	cc = this->helper_intersect_c_t( in[0], in[1], x, y, r );
	struct pqsolve i2;
	i2 = stable_pqsolver( aa, bb, cc ); //only solution "on the edge segment" in [0,1] accepted
	if ( i2.nsol == ONE_SOLUTION ) {
		if ( 0.0 <= i2.sol1 && i2.sol1 <= 1.0 ) { //the only solution on the segment?
			ip2[0] = in[0] + i2.sol1*t2[0];
			ip2[1] = in[1] + i2.sol1*t2[1];
		}
		else { std::cerr << "ERROR::Inconsistency v1e3 i2 ONESOLUTION!" << std::endl; return INTERSECTION_DETECTION_ERROR; }
	}
	if ( i2.nsol == TWO_SOLUTIONS ) { //##MK::shold not be possible
		found = false;
		if ( 0.0 <= i2.sol1 && i2.sol1 <= 1.0 ) {
			found = true;
			ip2[0] = in[0] + i2.sol1*t2[0];
			ip2[1] = in[1] + i2.sol1*t2[1];
		}
		if ( 0.0 <= i2.sol2 && i2.sol2 <= 1.0 ) {
			found = true;
			ip2[0] = in[0] + i2.sol2*t2[0];
			ip2[1] = in[1] + i2.sol2*t2[1];
		}
		if ( found == false ) { std::cerr << "ERROR::Inconsistency v1e3 i2 TWOSOLUTIONS!" << std::endl; return INTERSECTION_DETECTION_ERROR; }
	}
	else { std::cerr << "ERROR::No solution v1e3 i2 NOSOLUTION!" << std::endl; return INTERSECTION_DETECTION_ERROR; }

	//furthermore we need the two intersection points with the opposite edge, this edge however connects adj[0],adj[1] with adj[2], adj[3]
	double t3[2], ip3[2], ip4[2];
	double nx, ny, bx, by, tx, ty;
	t3[0] = adj[2] - adj[0];
	t3[1] = adj[3] - adj[1];
	aa = helper_intersect_a_t( t3[0], t3[1] );
	bb = helper_intersect_b_t( adj[0], adj[1], t3[0], t3[1], x, y );
	cc = helper_intersect_c_t( adj[0], adj[1], x, y, r );
	struct pqsolve i3;
	i3 = stable_pqsolver( aa, bb, cc );
	if ( i3.nsol == TWO_SOLUTIONS ) { //vertices should span an edge connecting two opposites sites of the circle
		if ( 0.0 <= i3.sol1 && i3.sol1 <= 1.0 ) {
			ip3[0] = adj[0] + i3.sol1*t3[0];
			ip3[1] = adj[1] + i3.sol1*t3[1];
		}
		else { std::cerr << "ERROR::Inconsistency v1e3 ip3!" << std::endl; return INTERSECTION_DETECTION_ERROR; }
		if ( 0.0 <= i3.sol2 && i3.sol2 <= 1.0 ) {
			ip4[0] = adj[0] + i3.sol2*t3[0];
			ip4[1] = adj[1] + i3.sol2*t3[1];
		}
		else { std::cerr << "ERROR::Inconsistency v1e3 ip3!" << std::endl; return INTERSECTION_DETECTION_ERROR; }
	}
	else if ( i3.nsol == ONE_SOLUTION ) {
		//MK::as the vertices in in and adj are the vertices referred to by a, b, c the function intersection_area_triangle_circle already evaluated for intersection or touching
		//of an edge with the circle, however there we probe != NO_SOLUTION, i.e. ONE or TWO SOLUTIONS because of this one can have here the
		//special case in which the opposite edge touches the circle at a point adj[i]+i3.sol1*t3[i] which is not covered by the v1e2 case because the aforementioned 
		//detections in intersection_area evaluated in such case of != NO_SOLUTION to edge intersects/touches
		//so it remains to compute the triangle in, ip2, ip3 and the cap atop
		double A = triangle_area( ip1[0]-in[0], ip1[1]-in[1], ip2[0]-in[0], ip2[1]-in[1] );

		//now the circular cap which the connecting line ip1-ip2 cuts of from the circle...
		bx = 0.5*(ip1[0]+ip2[0]);
		by = 0.5*(ip1[1]+ip2[1]);
		tx = ip2[0]-ip1[0];
		ty = ip2[1]-ip1[1];
		double outernormal[2];
		outer_triedge_normal( in[0],in[1], ip1[0], ip1[1], ip2[0], ip2[1], bx, by, tx, ty, outernormal );
		nx = outernormal[0]; //ty;
		ny = outernormal[1]; //-1.0*tx;
		aa = helper_bisec_a_t( nx, ny );
		bb = helper_bisec_b_t( bx, by, nx, ny, x, y );
		cc = helper_bisec_c_t( bx, by, x, y, r );
		struct pqsolve cap;
		cap = this->stable_pqsolver( aa, bb, cc );
		if ( cap.nsol != TWO_SOLUTIONS ) { std::cerr << "ERROR::Inconsistency in v1e3 special case tri/cap I dont find TWO SOLUTIONS!" << std::endl; return INTERSECTION_DETECTION_ERROR; }
		short consistency = 0;
		double dist;
		if ( cap.sol1 > 0.0 ) {
			consistency++;
			dist = sqrt(SQR(cap.sol1*nx)+SQR(cap.sol1*ny));
		}
		if ( cap.sol2 > 0.0 ) {
			consistency++;
			dist = sqrt(SQR(cap.sol2*nx)+SQR(cap.sol2*ny));
		}
		if ( consistency == 0 || consistency == 2 ) { std::cerr << "ERROR::Inconsistency in v1e3 special case tri/cap TWOSOLUTION!" << std::endl; return INTERSECTION_DETECTION_ERROR; }
		A = A + circular_cap_area(dist, r );
		return A;
		/*if ( cap.nsol == TWO_SOLUTIONS ) {
			short consistency = 0;
			double dist = std::numeric_limits<double>:: max(); //shortest
			if ( point_in_triangle( bx+cap.sol1*nx, by+cap.sol1*ny, xa, ya, xb, yb, xc, yc ) == true ) {
				consistency++;
				if ( sqrt( SQR(cap.sol1*nx) + SQR(cap.sol1*ny) ) < dist ) 
					dist = sqrt( SQR(cap.sol1*nx) + SQR(cap.sol1*ny) );
			}
			if ( point_in_triangle( bx+cap.sol2*nx, by+cap.sol2*ny, xa, ya, xb, yb, xc, yc ) == true ) {
				consistency++;
				if ( sqrt( SQR(cap.sol2*nx) + SQR(cap.sol2*ny) ) < dist ) 
					dist = sqrt( SQR(cap.sol2*nx) + SQR(cap.sol2*ny) );
			}
			if ( consistency == 0 ) { std::cerr << "ERROR::Inconsistency in v1e3 special case tri/cap TWOSOLUTION!" << std::endl; return 0.0; }
			A = A + circular_cap_area(dist, r );
			return A;
		}
		else if ( cap.nsol == ONE_SOLUTION ) { std::cerr << "ERROR::Inconsistency in v1e3 special case tri/cap ONESOLUTION!" << std::endl; return 0.0; }
		else { std::cerr << "ERROR::Inconsistency in v1e3 special case tri/cap NOSOLUTION!" << std::endl; return 0.0; }*/
	}
	else { std::cerr << "ERROR::Inconsistency v1e3 ip3 and ip4 less ONE_SOLUTION, so it would be a v1e2 case but wasnt detected as such!" << std::endl; return INTERSECTION_DETECTION_ERROR; }

	//now we have a triangle in, ip1, ip2 surplus a quad polygon comprising the points ip1, ip2, ip3, ip4
	double A = triangle_area( ip1[0]-in[0], ip1[1]-in[1], ip2[0]-in[0], ip2[1]-in[1] );
	A = A + triangle_area( ip4[0]-ip1[0], ip4[1]-ip1[1], ip2[0]-ip1[0], ip2[1]-ip1[1] ) + triangle_area( ip4[0]-ip1[0], ip4[1]-ip1[1], ip3[0]-ip1[0], ip3[1]-ip1[1] );

#ifdef HEURISTIC_HANDLING_OF_SMALL_NUMBERS
	if ( A < 0.0 ) {
		std::cerr << "ERROR::Inscribed polygon of triangle surplus quad in v1e3 has negative area!" << std::endl; return INTERSECTION_DETECTION_ERROR;
	}
	if ( A <= MINIMUM_CRITICAL_INTERSECTION_AREA ) {
//		std::cerr << "WARNING::v1e3 expects very small cap height therefore returning less accurate value!" << std::endl;
		return A;
	}
#endif

	//it remains to add the two circular caps on the connecting lines ip1ip3 and ip2ip4
	bx = 0.5*(ip1[0]+ip3[0]);
	by = 0.5*(ip1[1]+ip3[1]);
	tx = ip3[0]-ip1[0];
	ty = ip3[1]-ip1[1];
	double outernormal[2];
	outer_quadedge_normal( (ip1[0]+ip2[0]+ip3[0]+ip4[0])/4.0, (ip1[1]+ip2[1]+ip3[1]+ip4[1])/4.0, bx, by, tx, ty, outernormal );
	nx = outernormal[0]; //ty; //ip3[1] - ip1[1];
	ny = outernormal[1]; //-1.0*tx; //(ip3[0]-ip1[0]); //perpendicular to the vector point from ip1 to ip3 not normalized, neither guaranteed to be inner or outer normal!

	aa = helper_bisec_a_t( nx, ny );
	bb = helper_bisec_b_t( bx, by, nx, ny, x, y );
	cc = helper_bisec_c_t( bx, by ,x , y, r );
	struct pqsolve cap1;
	cap1 = stable_pqsolver( aa, bb, cc );
	/*if ( cap1.nsol == TWO_SOLUTIONS ) {
		short consistency = 0;
		double dist = std::numeric_limits<double>:: max(); //shortest
		if ( point_in_triangle( bx+cap1.sol1*nx, by+cap1.sol1*ny, xa,ya,xb,yb,xc,yc ) == true ) {
			consistency++;
			if ( sqrt( SQR(cap1.sol1*nx) + SQR(cap1.sol1*ny) ) < dist ) 
				dist = sqrt( SQR(cap1.sol1*nx) + SQR(cap1.sol1*ny) );
		}
		if ( point_in_triangle( bx+cap1.sol2*nx, by+cap1.sol2*ny, xa,ya,xb,yb,xc,yc ) == true ) {
			consistency++;
			if ( sqrt( SQR(cap1.sol2*nx) + SQR(cap1.sol2*ny) ) < dist )
				dist = sqrt( SQR(cap1.sol2*nx) + SQR(cap1.sol2*ny) );
		}
		if ( consistency == 0 ) { std::cerr << "ERROR::Unexpected inconsistency in cap determination v1e3!" << std::endl; return 0.0; }
		A = A + circular_cap_area( dist, r );
	}
	else { std::cerr << "ERROR::Inconsistency cap1 v1e3?" << std::endl; return 0.0; }*/
	//if only one solution then circle touches at intersection point ie ip1 and ip3 overlap however this is inconsistent with the fact that the 1. adjacent point was detected outside!
	if ( cap1.nsol != TWO_SOLUTIONS ) { std::cerr << "ERROR::Inconsistency cap1 v1e3?" << std::endl; return INTERSECTION_DETECTION_ERROR; }
	short consistency = 0;
	double dist;
	if ( cap1.sol1 > 0.0 ) {
		consistency++;
		dist = sqrt(SQR(cap1.sol1*nx)+SQR(cap1.sol1*ny));
	}
	if ( cap1.sol2 > 0.0 ) {
		consistency++;
		dist = sqrt(SQR(cap1.sol2*nx)+SQR(cap1.sol2*ny));
	}
	//TWO SOLUTIONS BUT NONE WITH POSITIVE VALUE, there is evidence that this is a numerical issue, for instance see:
	//this configuration where one cap's height approaches again double precision limit, therefore
	//192065;14;5.97488499999976439e-08;-666;0.753942000000000001;0.0186133999999999987;0.75586399999999998;0.0187019000000000005;0.754125000000000045;0.0186839999999999992;0.754787449151937628;0.0180793110830934668;0.0010000000000000000';
#ifdef HEURISTIC_HANDLING_OF_SMALL_NUMBERS
	if ( consistency == 0 ) { std::cerr << "WARNING::Very small cap1 height in v1e3 therefore less accuracy area!" << std::endl; return A; }
	if ( consistency == 2 ) { std::cerr << "ERROR::Inconsistency cap1 v1e3!" << std::endl; return INTERSECTION_DETECTION_ERROR; }
#else
	if ( consistency == 0 || consistency == 2 ) { std::cerr << "ERROR::Inconsistency cap1 v1e3!" << std::endl; return INTERSECTION_DETECTION_ERROR; }
#endif
	A = A + circular_cap_area( dist, r );


	//cap2 analogously
	bx = 0.5*(ip4[0]+ip2[0]);
	by = 0.5*(ip4[1]+ip2[1]);
	tx = ip4[0]-ip2[0];
	ty = ip4[1]-ip2[1];
	outer_quadedge_normal( (ip1[0]+ip2[0]+ip3[0]+ip4[0])/4.0, (ip1[1]+ip2[1]+ip3[1]+ip4[1])/4.0, bx, by, tx, ty, outernormal );
	nx = outernormal[0]; //ty; //ip4[1]-ip2[1]; //##MK::be careful, this can connect the wrong edges...
	ny = outernormal[1]; //-1.0*tx; //(ip4[0]-ip2[0]);

	aa = helper_bisec_a_t( nx, ny );
	bb = helper_bisec_b_t( bx, by, nx, ny, x, y );
	cc = helper_bisec_c_t( bx, by ,x , y, r );
	struct pqsolve cap2;
	cap2 = stable_pqsolver( aa, bb, cc );
	/*if ( cap2.nsol == TWO_SOLUTIONS ) {
		short consistency = 0;
		double dist = std::numeric_limits<double>:: max(); //shortest
		if ( point_in_triangle( bx+cap2.sol1*nx, by+cap2.sol1*ny, xa,ya,xb,yb,xc,yc ) == true ) {
			consistency++;
			if ( sqrt( SQR(cap2.sol1*nx) + SQR(cap2.sol1*ny) ) < dist ) 
				dist = sqrt( SQR(cap2.sol1*nx) + SQR(cap2.sol1*ny) );
		}
		if ( point_in_triangle( bx+cap2.sol2*nx, by+cap2.sol2*ny, xa,ya,xb,yb,xc,yc ) == true ) {
			consistency++;
			if ( sqrt( SQR(cap2.sol2*nx) + SQR(cap2.sol2*ny) ) < dist )
				dist = sqrt( SQR(cap2.sol2*nx) + SQR(cap2.sol2*ny) );
		}
		if ( consistency == 0 || consistency == 2 ) { std::cerr << "ERROR::Unexpected inconsistency in cap determination v1e3!" << std::endl; return 0.0; }
		A = A + circular_cap_area( dist, r );
	}
	else { std::cerr << "ERROR::Inconsistency cap1 v2e3?" << std::endl; return 0.0; }*/
	if ( cap2.nsol != TWO_SOLUTIONS ) { std::cerr << "ERROR::Inconsistency cap2 v1e3?" << std::endl; return INTERSECTION_DETECTION_ERROR; }
	consistency = 0;
	if ( cap2.sol1 > 0.0 ) {
		consistency++;
		dist = sqrt(SQR(cap2.sol1*nx)+SQR(cap2.sol1*ny));
	}
	if ( cap2.sol2 > 0.0 ) {
		consistency++;
		dist = sqrt(SQR(cap2.sol2*nx)+SQR(cap2.sol2*ny));
	}
	//same story as for cap1
#ifdef HEURISTIC_HANDLING_OF_SMALL_NUMBERS
	if ( consistency == 0 ) { std::cerr << "WARNING::Very small cap2 height in v1e3 therefore less accuracy area!" << std::endl; return A; }
	if ( consistency == 2 ) { std::cerr << "ERROR::Inconsistency cap2 v1e3!" << std::endl; return INTERSECTION_DETECTION_ERROR; }
#else
	if ( consistency == 0 || consistency == 2 ) { std::cerr << "ERROR::Inconsistency cap2 v1e3!" << std::endl; return INTERSECTION_DETECTION_ERROR; }
#endif
	A = A + circular_cap_area( dist, r );

	return A;
}


double mathMethods::tri_cir_v2e3( double xa, double ya, double xb, double yb, double xc, double yc, double x, double y, double r, bool a, bool b, bool c, bool ab, bool bc, bool ca ) {
#ifdef DEBUG_TRIANGLE
	std::cout << "v2e3" << std::endl;
	#ifdef DEBUG_EARLY_OUT
		return 0.0;
	#endif
#endif
	//computes intersection area of arbitrary nondegenerated triangle with a circle with two of its vertices in the circle the remaining outside
	//(the ones for which a,b,c == true) and two edges cutting the circle while the remaining edge inside the circle
	//partition triangle into a) a polygon completely inside by cutting the old triangle with a new edge (connecting the two edges' intersection points)
	//+ b) circular cap area additional intersection with the triangle

	//which vertex is outside?
	double out[2], in[4];
	unsigned int find = 0;
	if ( a == true ) {
		in[find+0] = xa;		in[find+1] = ya;
		find = find + 2; //my first inner vertex
	}
	else { out[0] = xa;		out[1] = ya; }
	if ( b == true ) {
		in[find+0] = xb;		in[find+1] = yb;
		find = find + 2;
	}
	else { out[0] = xb;	out[1] = yb; }
	if ( c == true  ) {
		in[find+0] = xc;		in[find+1] = yc;
		find = find + 2;
	}
	else { out[0] = xc;	out[1] = yc; }
	if ( find != 4 ) { std::cerr << "ERROR::v2e3 Logic error!" << std::endl; return INTERSECTION_DETECTION_ERROR; }

	//we know now already which edges cut the triangle

	//which edges cut the triangle, must be those connecting the inner points in with the outer out
	double t1[2], t2[2], ip1[2], ip2[2];
	double aa, bb, cc;
	t1[0] = out[0] - in[0];	t1[1] = out[1] - in[1];
	t2[0] = out[0] - in[2];	t2[1] = out[1] - in[3];

	aa = helper_intersect_a_t( t1[0], t1[1] );
	bb = helper_intersect_b_t( in[0], in[1], t1[0], t1[1], x, y );
	cc = helper_intersect_c_t( in[0], in[1], x, y, r );
	struct pqsolve i1;
	i1 = stable_pqsolver( aa, bb, cc ); //only solutions "on the edge segment" in [0,1] accepted
	if ( i1.nsol == ONE_SOLUTION ) {
		if ( 0.0 <= i1.sol1 && i1.sol1 <= 1.0 ) { //first solution on segment?
			ip1[0] = in[0] + i1.sol1*t1[0];
			ip1[1] = in[1] + i1.sol1*t1[1];
		}
		else { std::cerr << "ERROR::Inconsistency v2e3 i1 ONESOLUTION!" << std::endl; return INTERSECTION_DETECTION_ERROR; }
	}
	else if ( i1.nsol == TWO_SOLUTIONS ) { //##MK::should not be possible as one vertex is inside the other outside
		if ( 0.0 <= i1.sol1 && i1.sol1 <= 1.0 ) { //first solution on segment?
			ip1[0] = in[0] + i1.sol1*t1[0];
			ip1[1] = in[1] + i1.sol1*t1[1];
		}
		else if ( 0.0 <= i1.sol2 && i1.sol2 <= 1.0 ) { //no, okay but the second?
			ip1[0] = in[0] + i1.sol2*t1[0];
			ip1[1] = in[1] + i1.sol2*t1[1];
		}
		else { std::cerr << "ERROR::Inconsistency v2e3 i1 TWOSOLUTIONS!" << std::endl; return INTERSECTION_DETECTION_ERROR; }
	}
	else { std::cerr << "ERROR::No solution v2e3 i1 NOSOLUTION!" << std::endl; return INTERSECTION_DETECTION_ERROR; }

	aa = helper_intersect_a_t( t2[0], t2[1] );
	bb = this->helper_intersect_b_t( in[2], in[3], t2[0], t2[1], x, y );
	cc = this->helper_intersect_c_t( in[2], in[3], x, y, r );
	struct pqsolve i2;
	i2 = stable_pqsolver( aa, bb, cc ); //only solution "on the edge segment" in [0,1] accepted
	if ( i2.nsol == ONE_SOLUTION ) {
		if ( 0.0 <= i2.sol1 && i2.sol1 <= 1.0 ) { //the only solution on the segment?
			ip2[0] = in[2] + i2.sol1*t2[0];
			ip2[1] = in[3] + i2.sol1*t2[1];
		}
		else { std::cerr << "ERROR::Inconsistency v2e3 i2 ONESOLUTION!" << std::endl; return INTERSECTION_DETECTION_ERROR; }
	}
	else if ( i2.nsol == TWO_SOLUTIONS ) {
		if ( 0.0 <= i2.sol1 && i2.sol1 <= 1.0 ) { //first solution on segment?
			ip2[0] = in[2] + i2.sol1*t2[0];
			ip2[1] = in[3] + i2.sol1*t2[1];
		}
		else if ( 0.0 <= i2.sol2 && i2.sol2 <= 1.0 ) { //no, okay but the second?
			ip2[0] = in[2] + i2.sol2*t2[0];
			ip2[1] = in[3] + i2.sol2*t2[1];
		}
		else { std::cerr << "ERROR::Inconsistency v2e3 i2 TWOSOLUTIONS!" << std::endl; return INTERSECTION_DETECTION_ERROR; }
	}
	else { std::cerr << "ERROR::No solution v2e3 i2 NOSOLUTION!" << std::endl; return INTERSECTION_DETECTION_ERROR; }

	//knowing now the included points in[0], in[1] and in[2], in[3] and the new intersection points i1[2] and i2[2]
	//we can compute the polygon area spanned by these vertices as the sum of two triangles
	double A = triangle_area( in[0]-ip1[0], in[1]-ip1[1], in[2]-ip1[0], in[3]-ip1[1] );
	A = A + triangle_area( ip2[0]-ip1[0], ip2[1]-ip1[1], in[2]-ip1[0], in[3]-ip1[1] );

#ifdef HEURISTIC_HANDLING_OF_SMALL_NUMBERS
	if ( A < 0.0 ) {
		std::cerr << "ERROR::Inscribed quad in v2e3 has negative area!" << std::endl; return INTERSECTION_DETECTION_ERROR;
	}
	if ( A <= MINIMUM_CRITICAL_INTERSECTION_AREA ) {
//		std::cerr << "WARNING::v2e3 expects very small cap height therefore returning less accurate value!" << std::endl;
		return A;
	}
#endif

	//it remains to get the area of the circular cap which is cut off from the circle by the segment line i1i2
	double dist;
	double tx = ip2[0]-ip1[0];
	double ty = ip2[1]-ip1[1];
	double bx = 0.5*(ip1[0]+ip2[0]);
	double by = 0.5*(ip1[1]+ip2[1]);
	//double outernormal[2];
	//outer_triedge_normal( xa,ya,xb,yb,xc,yc,bx,by,tx,ty, outernormal );
	double nx = ty; //ip2[1]-ip1[1];
	double ny = -1.0*tx; //(ip2[0]-ip1[0]);

	//MK::the included vertex inx, iny is not necessarily located at the center of the circle!
	struct pqsolve bi;
	aa = helper_bisec_a_t( nx, ny );
	bb = helper_bisec_b_t( bx, by, nx, ny, x, y );
	cc = helper_bisec_c_t( bx, by, x, y, r );
	bi = stable_pqsolver( aa, bb, cc);
	if ( bi.nsol == TWO_SOLUTIONS ) { //which point is in triangle
		if ( point_in_triangle( bx+bi.sol1*nx, by+bi.sol1*ny, xa,ya,xb,yb,xc,yc ) == true ) {
			dist = sqrt( SQR(bi.sol1*nx) + SQR(bi.sol1*ny) );
			A = A + circular_cap_area( dist, r );
			return A;
		}
		else if ( point_in_triangle( bx+bi.sol2*nx, by+bi.sol2*ny, xa,ya,xb,yb,xc,yc ) == true ) {
			dist = sqrt( SQR(bi.sol2*nx) + SQR(bi.sol2*ny) );
			A = A + circular_cap_area( dist, r );
			return A;
		}
		//##MK::what about the spherical cap degenerated to its negative?
		else { std::cerr << "ERROR::Inconsistency v2e3 bi two solutions!" << std::endl; return INTERSECTION_DETECTION_ERROR; }
	}
	else { std::cerr << "WARNING::Inconsistency because ONE_SOLUTION, i.e. circle edge touching for bisector!" << std::endl; return INTERSECTION_DETECTION_ERROR; }
}


double mathMethods::intersection_area_triangle_circle( double x1, double y1, double x2, double y2, double x3, double y3, double x, double y, double r ) {
#ifdef DEBUG_TRIANGLE
	//std::cout << "Analysing triangle << x1y1x2y2x3y3 are all = " << x1 << ";" << y1 << ";" << x2 << ";" << y2 << ";" << x3 << ";" << y3 << std::endl;
#endif
	//determine how many points are in the circle
	bool vert[3] = { false, false, false };
	vert[0] = point_in_circle( x1, y1, x, y, r );
	vert[1] = point_in_circle( x2, y2, x, y, r );
	vert[2] = point_in_circle( x3, y3, x, y, r );
	unsigned int nvert = 0;
	if ( vert[0] == true ) nvert++;
	if ( vert[1] == true ) nvert++;
	if ( vert[2] == true ) nvert++;

	//three vertices inside? then triangle completely inside the circle
	if ( nvert == 3 ) {
#ifdef DEBUG_TRIANGLE
	std::cout << "v3e3" << std::endl;
#endif
		return triangle_area( x2-x1, y2-y1, x3-x1, y3-y1 );
	}

	//all other cases, check whether triangle edges cut or are completely inside the circle
	bool edge[3] = { false, false, false }; //check for complete enclosure first
	if ( vert[0] == true && vert[1] == true ) edge[0] = true; //edge ab
	if ( vert[1] == true && vert[2] == true ) edge[1] = true; //edge bc
	if ( vert[2] == true && vert[0] == true ) edge[2] = true; //edge ca
	
	double tx, ty, a, b, c; //partial enclosure, compute intersection segment with circle points
	tx = x2-x1; //probe edge ab
	ty = y2-y1;
	a = helper_intersect_a_t( tx, ty );
	b = helper_intersect_b_t( x1, y1, tx, ty, x, y );
	c = helper_intersect_c_t( x1, y1, x, y, r );
	struct pqsolve ar;
	ar = stable_pqsolver( a, b, c); //only accept solutions \in [0,1]
	if ( ar.nsol != NO_SOLUTION ) { //return -1.0 for wrong or nonexistent solutions...
		if ( 0.0 <= ar.sol1 && ar.sol1 <= 1.0 ) edge[0] = true;
		if ( 0.0 <= ar.sol2 && ar.sol2 <= 1.0 ) edge[0] = true;
	}
	ar.nsol = NO_SOLUTION; ar.sol1 = INVALID_SOL; ar.sol2 = INVALID_SOL; 
	tx = x3-x2; //probe edge bc
	ty = y3-y2;
	a = helper_intersect_a_t( tx, ty );
	b = helper_intersect_b_t( x2, y2, tx, ty, x, y );
	c = helper_intersect_c_t( x2, y2, x, y, r );
	ar = stable_pqsolver( a, b, c);
	if ( ar.nsol != NO_SOLUTION ) {
		if ( 0.0 <= ar.sol1 && ar.sol1 <= 1.0 ) edge[1] = true;
		if ( 0.0 <= ar.sol2 && ar.sol2 <= 1.0 ) edge[1] = true;
	}
	ar.nsol = NO_SOLUTION; ar.sol1 = INVALID_SOL; ar.sol2 = INVALID_SOL;
	tx = x1-x3; //probe edge ca
	ty = y1-y3;
	a = helper_intersect_a_t( tx, ty);
	b = helper_intersect_b_t( x3, y3, tx, ty, x, y );
	c = helper_intersect_c_t( x3, y3, x, y, r );
	ar = stable_pqsolver( a, b, c);
	if ( ar.nsol != NO_SOLUTION ) {
		if ( 0.0 <= ar.sol1 && ar.sol1 <= 1.0 ) edge[2] = true;
		if ( 0.0 <= ar.sol2 && ar.sol2 <= 1.0 ) edge[2] = true;
	}

	//summarize how many intersecting or enclosed edges the triangle has with the circle
	unsigned int nedge = 0;
	if ( edge[0] == true ) nedge++;
	if ( edge[1] == true ) nedge++;
	if ( edge[2] == true ) nedge++;

#ifdef DEBUG_TRIANGLE
	if ( nvert > 0 || nedge > 0 )
		std::cout << "Configuration vert__" << vert[0] << vert[1] << vert[2] << "__edge__" << edge[0] << edge[1] << edge[2] << std::endl;
#endif

	//the only possible case of nvert == 3 (complete triangle enclosure) was already handled, now cater the remaining ones
	if ( nvert == 2 ) {
		if ( nedge == 3 ) {
			return tri_cir_v2e3( x1, y1, x2, y2, x3, y3, x, y, r, vert[0], vert[1], vert[2], edge[0], edge[1], edge[2] );
		}
		else {
			if ( nedge == 0 ) { std::cerr << "ERROR::v2e0 should be impossible!" << std::endl; return INTERSECTION_DETECTION_ERROR; }
			else if ( nedge == 1 ) { std::cerr << "ERROR::v2e1 should be impossible!" << std::endl; return INTERSECTION_DETECTION_ERROR; }
			else { std::cerr << "ERROR::v2e2 should be impossible!" << std::endl; return INTERSECTION_DETECTION_ERROR; }
		}
	}
	else if ( nvert == 1 ) {
		if ( nedge == 3 ) {
			return tri_cir_v1e3( x1, y1, x2, y2, x3, y3, x, y, r, vert[0], vert[1], vert[2], edge[0], edge[1], edge[2] );
		}
		else if ( nedge == 2 ) {
			return tri_cir_v1e2( x1, y1, x2, y2, x3, y3, x, y, r, vert[0], vert[1], vert[2], edge[0], edge[1], edge[2] );
		}
		else {
			if ( nedge == 1 ) { std::cerr << "ERROR::v1e1 should be impossible!" << std::endl; return INTERSECTION_DETECTION_ERROR; }
			else { std::cerr << "ERROR::v1e0 should be impossible!" << std::endl; return INTERSECTION_DETECTION_ERROR; }
		}
	}
	else if ( nvert == 0 ) {
		if ( nedge == 0 ) {
			if ( point_in_triangle( x, y, x1, y1, x2, y2, x3, y3 ) == true ) //circle fully in triangle
				return PI * SQR(r);
			else //triangle completely outside
				return 0.0;
		}
		else if ( nedge == 1 ) {
			if ( edge[0] == true ) { //cap or cap negative protruding into triangle
				return tri_cir_v0e1( x1, y1, x2, y2, x, y, r , x3, y3 );
			} 
			else if ( edge[1] == true ) {
				return tri_cir_v0e1( x2, y2, x3, y3, x, y, r , x1, y1 );
			}
			else if ( edge[2] == true ) {
				return tri_cir_v0e1( x3, y3, x1, y1, x, y, r , x2, y2 );
			}
			else { std::cerr << "ERROR::v0e1 without any edge cutting!" << std::endl; return INTERSECTION_DETECTION_ERROR; }
		}
		else if ( nedge == 2 ) { //only two caps to substract
			return tri_cir_v0e2( x1, y1, x2, y2, x3, y3, x, y, r, edge[0], edge[1], edge[2] );
		}
		else if ( nedge == 3 ) { //three caps to substract
			return tri_cir_v0e3( x1, y1, x2, y2, x3, y3, x, y, r );
		}
		else { std::cerr << "ERROR::Unknown case for v0e..." << std::endl; return INTERSECTION_DETECTION_ERROR; } 
	}
	else { std::cerr << "ERROR::Inconsistency to identify an analysis case for triangle << x1y1x2y2x3y3 are all = " << x1 << ";" << y1 << ";" << x2 << ";" << y2 << ";" << x3 << ";" << y3 << endl; return INTERSECTION_DETECTION_ERROR; }
}

#define RASTERING_ACCURACY		(2.0)
double mathMethods::intersection_area_triangle_circle_mc( double x1, double y1, double x2, double y2, double x3, double y3, double x, double y, double r, double ld ) {
	//rasterization Monte Carlo approach to approximate circle triangle overlap area
	//double xmin = std::numeric_limits<double>:: max();
	//double ymin = std::numeric_limits<double>:: max();
	//double xmax = std::numeric_limits<double>:: lowest();
	//double ymax = std::numeric_limits<double>:: lowest();
	double xmin = MIN(x1,MIN(x2,x3)); //[0,1]
	double xmax = MAX(x1,MAX(x2,x3));
	double ymin = MIN(y1,MIN(y2,y3));
	double ymax = MAX(y1,MAX(y2,y3));

	//raster window is *max - *min
	double refld = ld * RASTERING_ACCURACY;
	double h = 1.0 / refld;
	double promote, owin[4]; //xmi,xmx,ymi,ymx
	int truncate;
	//1.) floating point pixel on which the value falls, cast to int, substract/add safety, make to double again
	promote = xmin * refld;		truncate = promote;		truncate--;		owin[0] = h*((double) truncate) + 0.5*h;
	promote = xmax * refld;		truncate = promote;		truncate++;		owin[1] = h*((double) truncate) + 0.5*h;
	promote = ymin * refld;		truncate = promote;		truncate--;		owin[2] = h*((double) truncate) + 0.5*h;
	promote = ymax * refld;		truncate = promote;		truncate++;		owin[3] = h*((double) truncate) + 0.5*h;

	unsigned int total = 0;
	unsigned int in = 0;
	for ( double yy = owin[2]; yy < owin[3]; yy += h ) {
		for ( double xx = owin[0]; xx < owin[1]; xx += h ) {
			total++;

			//quick pretest for fast rejects
			if ( xx < (xmin - TRIANGLE_EPSILON)) continue;
			if ( xx > (xmax + TRIANGLE_EPSILON)) continue;
			if ( yy < (ymin - TRIANGLE_EPSILON)) continue;
			if ( yy > (ymax + TRIANGLE_EPSILON)) continue;

			//test for inclusion in circle is quicker than in triangle so do not bother testing both which will more likely become rejected though having put already effort into it
			if ( point_in_circle( xx, yy, x, y, r ) == false ) continue;

			//so point is definately in the circle and possibly in the triangle
			if ( point_in_triangle_fast( xx, yy, x1, y1, x2, y2, x3, y3 ) == true )
				in++;
		}
	}
	return ( ((double) in) * SQR(h) );
	//return ( ((double) in) / ((double) total) );
}


bool mathMethods::intersection_area_triangle_check( double A,  double x1, double y1, double x2, double y2, double x3, double y3, double x, double y, double r ) {
	//A cannot be larger than the triangle area
	if ( A > triangle_area( x2-x1, y2-y1, x3-x1, y3-y1 ) )
		return false;
	if ( A < 0.0 )
		return false;
	if ( A > PI*SQR(r) )
		return false;

	return true;
}


//TESTING CODE
/*
//calling example
struct pqsolve r;
r = worker->CircleLineSegmentIntersection( 1.5, 1.5, 3.0, 1.5, 1.0, 1.0, 2.000001 );
cout << "nsol/sol1/sol2 = " << r.nsol << "\t\t" << setprecision(18) << r.sol1 << "\t\t" << setprecision(18) << r.sol2 << endl;
r = worker->CircleLineSegmentIntersection( 3.0, 1.5, 3.0, 4.5, 1.0, 1.0, 2.000001 );
cout << "nsol/sol1/sol2 = " << r.nsol << "\t\t" << setprecision(18) << r.sol1 << "\t\t" << setprecision(18) << r.sol2 << endl;
r = worker->CircleLineSegmentIntersection( 1.5, 1.5, 3.0, 4.5, 1.0, 1.0, 2.000001 );
cout << "nsol/sol1/sol2 = " << r.nsol << "\t\t" << setprecision(18) << r.sol1 << "\t\t" << setprecision(18) << r.sol2 << endl;
*/