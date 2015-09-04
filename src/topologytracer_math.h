//MK::TopologyTracer is a software developed by Markus K{\u}hbach in 2015
//for questions and comments please contact markus.kuehbach at rwth-aachen.de, last major revision 20150720
//##MK::marks comments that where the code is in rather quick and dirty,i.e. parts of the code which have been validated operative as is but offer certainly room for further improvements


//###DEBUGpragma once
#ifndef _TOPOLOGYTRACER_MATH_H_
#define _TOPOLOGYTRACER_MATH_H_

#include "topologytracer_defs.h"
#include "topologytracer_random.h"


#define DOUBLE_ACCURACY 		5e-16				//2^-51 as determined with MATLAB eps(2)
#define RGBRANGE				255
#define MYMATH_STANDARD_SEED	-4256
#define MINIMUM_ANGULAR_SPREAD	(0.017453292)		//1.0/180.0_PI_
#define LAGBtoHAGBTransition	(0.261799388)		//15.0/180.0_PI_
#define MAXIMUM_MISORI_FCC		(1.09606677)		//FCC 62.8/180.0*_PI_
#define SYMMETRIES_IN_FCC		24
#define _PI_					3.1415926535897932384626433832795
#define SQR(a)					((a)*(a))


//IPFZ Coloring
#define FAIL_MYMATH_NUMERIC		(-1.0)
#define EPS_PROXIMITY_IPFZ		(0.01)
#define IPF_COLOR_STRETCH_R		(0.5)
#define IPF_COLOR_STRETCH_G		(0.5)
#define IPF_COLOR_STRETCH_B		(0.5)
#define EPS_ENVIRONMENT			(1e-7)


class randomClass;


class mathMethods
{
public:
	mathMethods( void );
	~mathMethods( void );


	//helper functions utilized in CORe
	void bubbleSort ( double arr [ ], int size );
	void sortInt( long arr [], long size );
	void swap ( double& x, double& y );
	void swapInt ( long& x, long& y );

	void factorizeIn3( long n, long * factors );
	char isPrime( long n );
	void sort(int n, double *ra);


	//statistical tests
	void K_S_Test( double * data1, unsigned long n1, double * data2, unsigned long n2, double * d, double * prob );
	double probks( double alam );


	//quaternion algebra
	void multiplyQuaternions( double *q, double* p, double* r );
	
	//convert orientations among parametrizations, such as Bunge (3,1,3) rotation convention and quaternion
	void euler2quaternion( double * euler, double * q );
	void quaternion2Euler( double * quat, double * euler );

	//calculate disorientation among two orientations in various parametrizations
	double misorientationCubic( double pa1, double Pa, double pa2, double pb1, double Pb, double pb2 );
	void misorientationQuaternionCubic( double* p, double* q, double* quat  );
	double misorientationCubicQxQ( double q01, double q11, double q21, double q31, double q02, double q12, double q22, double q32 );

	//generate new orientations with some scatter
	void randomOrientationShoemake( double * result );
	void randomMisorientationShoemake( double theta, double* qr );

	void rotateOrientation( double *oriOri, double angle, double u, double v, double w, double *newOri );
	void newOrientationFromReference( double *oriOri, double deviation, double *newOri );
	void newOrientationFromReferenceQuat( double *qbunge, double deviation, double *newquat );

	//##better not safe to use functions
	void randomMisorientationAxisConsideredLuis(  double * qref, double * qr, double maxDev );
	void newOrientationFromReferenceFixedAngularCone(double * oriOri, double maxDev,double angle, double u, double v, double w,double * newOri);


	//identify orientation by RGB scheme
	void devtorefEuler2RGB ( double *bunge, double *ideal, double maxDev, unsigned char *rgb); //blue channel stretch from 0.0 to maxDev in radian, all other orientations white

	//an own PRNG
	randomClass r;
};

#endif