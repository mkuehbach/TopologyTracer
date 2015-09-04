//MK::TopologyTracer is a software developed by Markus K{\u}hbach in 2015
//for questions and comments please contact markus.kuehbach at rwth-aachen.de, last major revision 20150720
//##MK::marks comments that where the code is in rather quick and dirty,i.e. parts of the code which have been validated operative as is but offer certainly room for further improvements

//MKLB::quaternion based math library written by Luis Barrales-Mora and Markus K{\u}hbach
//MK: Unless indicated, all functions in radians, quaternion algebra: q = -q define an equivalent rotation

#include "topologytracer_math.h"


mathMethods::mathMethods( void )
{
	r.init( DEFAULT_SEED );
}


mathMethods::~mathMethods( void )
{
}


char mathMethods::isPrime( long n )
{
	long n_root = (long) ( sqrt((double)n) + 0.5 );

	for( long i = 2; i <= n_root; i++ )
		if( n % i == 0 ) return 0;

	return 1;
}


void mathMethods::factorizeIn3( long n, long * factors )
{
	long f[3]={0};

	if( isPrime( n ) )
	{
		factors[0]=n;
		factors[1]=1;
		factors[2]=1;
		return;
	}

	long nhalf = (long) (0.5*n);
	long *Fac = (long *) calloc( nhalf, sizeof( long ) );

	//if( !Fac ) exitus("ERROR: Cannot allocate more memory :(");

	int count=0;

	for( long i = nhalf; i > 0; i-- )
	{
		if( n % i==0 )
		{
			Fac[count]=i;
			count++;
		}
	}

	long ft = count;
	count=0;
	long sum = 3*n;

	for( long i=0;i<ft;i++ )
		for( long j=0;j<ft;j++ )
		{
			long p = Fac[i]*Fac[i]*Fac[j];
			long s = Fac[i]+Fac[i]+Fac[j];

			if( p == n && s < sum )
			{
				f[0]=Fac[i];
				f[1]=Fac[i];
				f[2]=Fac[j];
				sum = s;
			}
		}

		for( long i=0;i<ft;i++ )
			for( long j=i+1;j<ft;j++ )
				for( long k=j+1;k<ft;k++ )
				{
					long p = Fac[i]*Fac[j]*Fac[k];
					long s = Fac[i]+Fac[j]+Fac[k];
					if( p == n && s < sum )
					{
						f[0]=Fac[i];
						f[1]=Fac[j];
						f[2]=Fac[k];
						sum = s;
					}
				}

		sortInt( f,3 );
		factors[0]=f[2];
		factors[1]=f[1];
		factors[2]=f[0];
		free(Fac);
}


void mathMethods::bubbleSort ( double arr [ ], int size ) // Sort components of a quaternion in ASCENDING order
{
	int last = size - 2; 
	int isChanged = 1; 

	while ( last >= 0 && isChanged ) {
		isChanged = 0; 
		for ( int k = 0; k <= last; k++ )
			if ( arr[k] > arr[k+1] ) {
				swap ( arr[k], arr[k+1] );
				isChanged = 1; 
			}

		last--;
	}
}


void mathMethods::sortInt( long arr [ ], long size ) // Sort integers
{
	long last = size - 2; 
	long isChanged = 1; 

	while ( last >= 0 && isChanged ) {
		isChanged = 0; 
		for ( long k = 0; k <= last; k++ )
			if ( arr[k] > arr[k+1] ) {
				swapInt ( arr[k], arr[k+1] );
				isChanged = 1; 
			}

		last--;
	}
}


void mathMethods::swapInt( long& x, long& y ) //Required for the sorting
{
	long temp;
	temp = x;
	x = y;
	y = temp;
}


void mathMethods::swap( double& x, double& y ) //Required for the sorting
{
	double temp;
	temp = x;
	x = y;
	y = temp;
}


void mathMethods::K_S_Test( double * data1, unsigned long n1, double * data2, unsigned long n2, double * d, double * prob )
{
	unsigned long j1 = 1;
	unsigned long j2 = 1;
	double d1, d2, dt, en1, en2, en;
	double fn1 = 0.0;
	double fn2 = 0.0;

	sort( n1, data1 );
	sort( n2, data2 );

	en1 = n1;
	en2 = n2;
	*d = 0.0;

	while( j1 <= n1 && j2 <= n2 ) {
		if( (d1=data1[j1]) <= (d2=data2[j2]) ) fn1 = j1++ / en1;
		if( d2 <= d1 ) fn2 = j2++/en2;
		if( (dt=fabs(fn2-fn1)) > *d ) *d = dt;
	}
	en = sqrt(en1*en2/(en1+en2));
	*prob = probks( (en+0.12+0.11/en)*(*d) );
}


double mathMethods::probks( double alam )
{
	int j;
	double a2, term;
	double fac = 2.0;
	double sum = 0.0;
	double termbf = 0.0;

	a2 = -2.0 * SQR(alam);
	for( j=1; j<=100; j++ )
	{
		term = fac * exp( a2 * SQR( j ) );
		sum += term;
		if( fabs(term) <= EPS1 * termbf || fabs(term) <= EPS2 * sum ) return sum;
		fac = -fac;
		termbf = fabs( fac );
	}
	return 1.0;
}


void mathMethods::sort(int n, double *ra)
{
	int l, j, ir, i;
	float rra;

	l = (n >> 1) + 1;
	ir = n;

	for ( ; ; )
	{
		if (l > 1)					/* still in hiring phase */
			rra = ra[--l];
		else						/* in retirement-and-promotion phase */
		{
			rra = ra[ir];			/* clear a space at end of array */
			ra[ir]=ra[1];			/* retire the top of the heap into it */
			if (--ir == 1) 			/* done with last promotion */
			{
				ra[1] = rra;
				return;
			}						/* end if */
		}							/* end else */
		i = l;						/* whether we are in the hiring phase */
		j = l << 1;					/* or promotion phase, we here set up */
		while ( j <= ir )
		{
			if ( (j < ir) && (ra[j] < ra[j + 1]) )
				++j;				/* compare to the better underling */
				if ( rra < ra[j] )	/* demote rra */
				{
					ra[i] = ra[j];
					j += (i = j);
				}
				else
					j = ir + 1;		/* this is rra's level; set j to */
		}							/* terminate the sift-down */
		ra[i] = rra;				/* put rra into its slot */
	}
}



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


void mathMethods::quaternion2Euler( double * quat, double * euler )
{
	//convention: Bunge, ZXZ, equal to case (3,1,3) as analyzed in Diebel, James, 2006:
	//Representing Attitude: Euler Angles, Unit Quaternions and Rotation Vectors
	//Gimbal lock situation analyzed following the line of Melcher et. al., Conversion of EBSD data by a quaternion based algorithm....
	//TECHNISCHE MECHANIK, 30, 4, (2010), 401 - 413

	double q0 = quat[0];
	double q1 = quat[1];
	double q2 = quat[2];
	double q3 = quat[3];
	double PHI, sP, phi1, phi2;

	double cosPHI = SQR(q3) - SQR(q2) - SQR(q1) + SQR(q0);
	double y0 =	2*q1*q3	-	2*q0*q2; //as  following Diebel equation 434
	double x0 =	2*q2*q3	+	2*q0*q1;
	double y1 =	2*q1*q3	+	2*q0*q2;
	double x1 = -	2*q2*q3	+	2*q0*q1; //but atan2(y,x) yields a domain error when both SQR(x) and SQR(y) <= DOUBLE_ACCURACY zero!

	//this approach works only for properly normalized unit quaternions...

	//acos(x) has goes to +pi for the argument x approaching from right to -1.0 that is when q3 and q0 are numerically zero
	if( cosPHI > (1.0 - DOUBLE_ACCURACY) ) cosPHI = 1.0; 
	if( cosPHI < (-1.0 + DOUBLE_ACCURACY) ) cosPHI = -1.0;

	//now application of acos function is safe to use... 
	PHI = acos( cosPHI );


	//special case: PHI=0.0, q1=q2=0 -->cosPHI=1
	//special case:and PHI=_PI_ q0=q3=0 -->cosPHI=-1
	//in both of which equation 434 would cause atan2(0,0) nevertheless then sin(PHI) is a singularity free indicator of gimbal lock occurs
	sP = sin(PHI);


	if( SQR(sP) > DOUBLE_ACCURACY ) {
		phi2 = atan2( y0 / sP, x0 / sP ); //##atan2( y0, x0 )
		phi1 = atan2( y1 / sP, x1 / sP ); //##atan2( y1, x1 )
	}
	else {
		//gimbal lock in the case PHI=0.0 first rotation lets ND' || ND0 and next as well so as if at all additive rotation about phi1+phi2 about ND0, choice for either phi1 or phi2 is arbitrary
		phi1 = atan2( 2*(q1*q2 + q0*q3), (SQR(q0) + SQR(q1) - SQR(q2) - SQR(q3)) );
		phi2 = 0.0; //arbitrary choice, Rollett for instance partitions equally phi1 and phi2 as in the case PHI = 0 but:: atan2(a12,a11) for phi1 is for instance inconsistent with Diebel

		//more safe but equally heuristic is to make use of Melchers acos(q0^2 - q3^2) = phi1+phi2 and partition randomly the angular contribution in the case PHI=0 
		//or acos(q1^2 - q2^2) = phi1 - phi2
	}


	//it always holds the under m-3m that the Euler space is symmetric to 0 <= phi1 <= 2*_PI_ , 0 <= PHI <= _PI_, 0 <= phi2 <= 2*_PI_
	if (phi1 < (0.0 - DOUBLE_ACCURACY) )
		phi1 += 2 * _PI_;
	if (phi2 < (0.0 - DOUBLE_ACCURACY) )
		phi2 += 2 * _PI_;


	euler[2] = phi2; //following the notation order used in Diebel, James 2006
	euler[1] = PHI;
	euler[0] = phi1;
}



void mathMethods::misorientationQuaternionCubic( double* p, double* q, double* quat  )
{
	//MK::ok
	double qm1[4];    //Inverse of quaternion q

	//Inverse of quaternion q is the same like the conjugate for unit quaternions
	qm1[0] = q[0];
	qm1[1] = -1.0 * q[1];
	qm1[2] = -1.0 * q[2];
	qm1[3] = -1.0 * q[3];

	double r[4]; //Resulting misorientation quaternion, m = pq-1

	multiplyQuaternions( qm1, p, r ); //MK:: was p post applied to qm1 but that was not consistent with MTex, in particular not with -(q^-1) = -q0,q1,q2,q3


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

	double rq[4];

	int mi;
	double max=0.0;

	for( int i=0;i<6;i++ )						//Determing the quaternion with the maximal component and the component itself
		for( int j=0;j<4;j++ )
		{
			if( fabs(r0[i][j]) > max )
			{
				max=fabs(r0[i][j]);
				mi=i;
			}
		}

	rq[0] = fabs( r0[mi][0] );					//Disorientation requires all components positive
	rq[1] = fabs( r0[mi][1] );
	rq[2] = fabs( r0[mi][2] );
	rq[3] = fabs( r0[mi][3] );

	bubbleSort( rq,4 );						//Sorting into ascending order, because a desorientation in the SST
									//requires a quaternion with q0>=q1>=q2>=q3 which represents a minimal 
	quat[0] = rq[3];						//rotation angle and an axis fulfilling h>k>l
	quat[1] = rq[2];
	quat[2] = rq[1];
	quat[3] = rq[0];
	//additionally it is required that rq3 >= sum of all others and rq3 * (sqrt(2)-1) >= rq2
}


double mathMethods::misorientationCubicQxQ( double q01, double q11, double q21, double q31, double q02, double q12, double q22, double q32 )
{
	//MK::ok
	double p[4] = {q01,q11,q21,q31};
	double q[4] = {q02,q12,q22,q32};

	double qm1[4];    //Inverse of quaternion q
	qm1[0] = q[0];
	qm1[1] = -1.0 * q[1];
	qm1[2] = -1.0 * q[2];
	qm1[3] = -1.0 * q[3];

	double r[4]; //Resulting quaternion, rotation of the two previous quaternions pq-1

	multiplyQuaternions( qm1, p, r ); //MK:: was p post applied to qm1 but that was not consistent with MTex, in particular not with -(q^-1) = -q0,q1,q2,q3

	//Now, we have to determine the smallest angle.

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
			if( fabs(r0[i][j]) > omega )
				omega=fabs(r0[i][j]);


	if( omega > 1.0 ) //avoid singularity of acos function
		omega = (double) (int) omega;

	omega=2*acos(omega);
	//QUICKASSERT( omega <= 1.099 );
	return omega;
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

	multiplyQuaternions( qm1, p, r ); //MK:: was p post applied to qm1 but that was not consistent with MTex, in particular not with -(q^-1) = -q0,q1,q2,q3
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


void mathMethods::randomOrientationShoemake( double * result )
{
	//K. Shoemake, Graphic Gems III (editor D. Kirk) CalTech pp124-134
	//##MK::mind quaternion order of Shoemake w  + i*v + j*x + k*z <=> 0 +  1 2 3 with 
	double q[4]={0,0,0,0};
	double qnorm;

	double X0 = r.leEcuyer();
	double X1 = r.leEcuyer();
	double X2 = r.leEcuyer();

	double r1 = sqrt(1-X0);
	double r2 = sqrt(X0);
	double theta1 = 2 * _PI_ * X1;
	double theta2 = 2 * _PI_ * X2;

	q[0] = r1 * sin(theta1); //w
	q[1] = r1 * cos(theta1); //v
	q[2] = r2 * sin(theta2); //x
	q[3] = r2 * cos(theta2); //z

	qnorm = sqrt( SQR(q[0]) + SQR(q[1]) + SQR(q[2]) + SQR(q[3]) );

	//normalize
	q[0] = q[0] / qnorm;
	q[1] = q[1] / qnorm;
	q[2] = q[2] / qnorm;
	q[3] = q[3] / qnorm;

	double angles[3] = {0.0};

	quaternion2Euler( q, angles );

	result[0] = angles[0];
	result[1] = angles[1];
	result[2] = angles[2];
}


void mathMethods::randomMisorientationShoemake( double theta, double* qr  )
{
	//##MK::practically possible, but maybe not mathematically sound
	double q[4] = {0.0, 0.0, 0.0, 0.0};

	double qcrit = cos( 0.5 * theta);
	if ( theta < MINIMUM_ANGULAR_SPREAD ) { //limit to assure the while loop to finish
		qcrit = cos( 0.5 * MINIMUM_ANGULAR_SPREAD );
	}

	while( q[0] < qcrit ) {

		double X0 = r.leEcuyer();
		double X1 = r.leEcuyer();
		double X2 = r.leEcuyer();

		double r1 = sqrt(1-X0);
		double r2 = sqrt(X0);
		double theta1 = 2 * _PI_ * X1;
		double theta2 = 2 * _PI_ * X2;

		q[0] = r1 * sin(theta1); //w
		q[1] = r1 * cos(theta1); //v
		q[2] = r2 * sin(theta2); //x
		q[3] = r2 * cos(theta2); //z

		double qnorm = sqrt( SQR(q[0]) + SQR(q[1]) + SQR(q[2]) + SQR(q[3]) );

		//normalize
		q[0] = q[0] / qnorm;
		q[1] = q[1] / qnorm;
		q[2] = q[2] / qnorm;
		q[3] = q[3] / qnorm;
		//definately this algorithm samples random on the SO(3) the resulting quaternion however is not necessarily a disorientation!
	}

	qr[0] = q[0];
	qr[1] = q[1];
	qr[2] = q[2];
	qr[3] = q[3];
}


void mathMethods::rotateOrientation( double *oriOri, double angle, double u, double v, double w, double *newOri ) //Angle in radians
{
	//MK::should not be used as not properly validated but with the additional comment is correct?
	double qori[4], qrot[4];
	double _norm = 1.0 / sqrt(SQR(u)+SQR(v)+SQR(w));

	euler2quaternion( oriOri, qori );

	double qref[4] = { cos( 0.5* angle ), u * _norm * sin( 0.5 * angle ), v * _norm * sin( 0.5 * angle ), w * _norm * sin( 0.5 * angle ) };

	//MK::ok, an orientation qref in triclinic space can be used to post rotated another rotation which is parameterized as a quaternion as here qori
	multiplyQuaternions( qref, qori, qrot );

	double euler[3] = {0};

	quaternion2Euler( qrot, euler );

	newOri[0] = euler[0];
	newOri[1] = euler[1];
	newOri[2] = euler[2];
}


void mathMethods::newOrientationFromReference( double *bunge, double deviation, double *newOri )
{
	//MK::ok
	double qrndmisori[4];
	double qbunge[4], rotated[4];

	euler2quaternion( bunge, qbunge );

	//deviation is internally limited to MINIMUM_ANGULAR_SPREAD if called
	randomMisorientationShoemake( deviation, qrndmisori );

	//applying a rotation on another quaternion parameterized orientation
	multiplyQuaternions( qrndmisori , qbunge, rotated );


	double newEuler[3];
	quaternion2Euler( rotated, newEuler );

	newOri[0] = newEuler[0];
	newOri[1] = newEuler[1];
	newOri[2] = newEuler[2];
}


void mathMethods::newOrientationFromReferenceQuat( double *qbunge, double deviation, double *newquat )
{
	//MK::ok, does the same as newOrientationFromReference but overloaded with quaternions
	double qrndmisori[4];

	//this is really a problem when deviation approaches zero
	randomMisorientationShoemake( deviation, qrndmisori );

	//applying a rotation on another quaternion parameterized orientation
	multiplyQuaternions( qrndmisori , qbunge, newquat );
}


void mathMethods::randomMisorientationAxisConsideredLuis(  double * qref, double * qr, double maxDev  )
{
	double theta = cos( 0.5 * _PI_ );

	double q[4] = {0.0};

	double dev = cos(0.5 * maxDev);
	if ( maxDev < MINIMUM_ANGULAR_SPREAD ) {
		dev = cos( 0.5 * MINIMUM_ANGULAR_SPREAD );
	}

	double refNorm = sqrt( SQR(qref[0]) + SQR(qref[1]) + SQR(qref[2]) + SQR(qref[3]) );

	qref[0] /= refNorm;
	qref[1] /= refNorm;
	qref[2] /= refNorm;
	qref[3] /= refNorm;

	//idea is to provide a random scatter about an angle as well that has not more than a disorientation in misorientation space against the 
	//already prescribed axis
	//###MK::this sampling method is wrong??
	while( theta < dev  ) {
		double s = r.leEcuyer();
		double sigma1 = sqrt(1-s);
		double sigma2 = sqrt(s);
		double theta1 = 2 * _PI_ * r.leEcuyer();
		double theta2 = 2 * _PI_ * r.leEcuyer();

		q[0] = fabs(sigma2*cos(theta2));
		q[1] = fabs(sigma1*sin(theta1));
		q[2] = fabs(sigma1*cos(theta1));
		q[3] = fabs(sigma2*sin(theta2));

		double norm = sqrt( SQR(q[0])+SQR(q[1])+SQR(q[2])+SQR(q[3]) );

		q[0] /= norm;
		q[1] /= norm;
		q[2] /= norm;
		q[3] /= norm;

		bubbleSort( q, 4 );

		//##MK::because it tries to establish a lower than dev disorientation but among qref, but the multiplicative terms here do not test for a disorientation
		//##MK::anyway it is here tacitly assumed that the disorientation among misorientations is comparable to that of pure orientations, i.e rotations mapping sample frame to crystal frame??
		theta = q[3]*qref[0] + q[2]*qref[1] + q[1]*qref[2] + q[0]*qref[3];

		//##MK::needs potential correction by misorientationCubicQxQ( qref, q ) and prior the Shoemake algorithm without fabs and bubbleSort...
	}

	qr[0] = q[3];
	qr[1] = q[2];
	qr[2] = q[1];
	qr[3] = q[0];
}


void mathMethods::newOrientationFromReferenceFixedAngularCone(double * oriOri, double maxDev,double angle, double u, double v, double w,double * newOri)
{
	//##MK::??heuristic approach but not expected completely right as what the sign convention of the rotation about the particular directed axis
	double qr[4];
	double ori[4], qideal[4];

	double _norm = 1.0 / sqrt( SQR(u)+SQR(v)+SQR(w) );

	euler2quaternion( oriOri, qideal );

	double quvw[4] = { cos( 0.5* angle ), u * _norm * sin( 0.5 * angle ), v * _norm * sin( 0.5 * angle ), w * _norm * sin( 0.5 * angle ) };

	//#####MK::seems not correct - sample a misorientation quaternion in the cone about the quaternion quvw
	randomMisorientationAxisConsideredLuis(  quvw, qr, maxDev  );

	//MK::??, if one can be sure that at least here qr is the particular angular constrainted rotation quaternion qr then indeed post multiplication on qideal seems correct
	multiplyQuaternions( qr, qideal, ori );

	double euler[3];
	quaternion2Euler( ori, euler );

	newOri[0] = euler[0];
	newOri[1] = euler[1];
	newOri[2] = euler[2];
}

void mathMethods::devtorefEuler2RGB ( double *bunge, double *ideal, double maxDev, unsigned char *rgb)
{
	//blue channel stretch from 0.0 to maxDev in radian, all other orientations white
	rgb[0] = RGBRANGE;
	rgb[1] = RGBRANGE;
	rgb[2] = RGBRANGE;
}
