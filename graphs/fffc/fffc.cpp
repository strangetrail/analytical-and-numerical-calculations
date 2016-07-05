#include <stdlib.h>
//
// Remove this in final release
//
//include stdio_h
//
#include <math.h>
#include "./libUtils.h"
#include "./libFFTcstm.h"
#include "./libCZTcstm.h"
#include "./fffc.h"
// Moved to header for test purposes:
/*
const char qRotation[][2][2] =
{
	{
		{  1,  1 },
		{  1,  1 },
	},
	{
		{  1, -1 },
		{ -1,  1 },
	},
	{
		{ -1, -1 },
		{ -1, -1 },
	},
	{
		{ -1,  1 },
		{  1, -1 },
	},
};
const double cos225 = cos( M_PI / 8 );
const double sin225 = sin( M_PI / 8 );
const double tri45 = cos( M_PI / 4 );
*/
// Moved to header for test purposes.
Complex **jQWP;

double Lmn( int alpha, int n, double x ) // alpha, m - azimuthal, n - radial
{
	// int    i, j;
	double LmnRes = 0;

	if ( n > 1 )
		LmnRes = (( alpha + 1 - x ) / n ) * Lmn( alpha + 1, n - 1, x ) - ( x / n ) * Lmn( alpha + 2, n - 2, x );
	else
		if ( n > 0 )
			LmnRes = alpha + 1 - x;
		else
			LmnRes = 1;
	return LmnRes;
}
void getHEvector( TETM vectype, LGBeam params, double x, double y, double z, Complex **EHlocal )
{
// DEBUG:
/*
int *iptr_tmp = new int;
*iptr_tmp = 0;
double dtmp = 0;
*/
// DEBUG.
	//U_CHAR     tmp_uc;
	int        ipow, ipowodd, ipow2odd,
	           itmp1, itmp2;
	double     imreCommon,
	           ro2, phi, r2, r,
	           tmp1, tmp2, tmp3, yztmp;
	iComplex   mult;
	Complex    xyCommon;
//
// define A001                          
// 	__BGNINIT                             
// 		__PPT(Complex) PPD(2) PPD(3) PPA(EH)
// 	__ENDINIT
//
// printf(STRINGIZE_TOKEN(A001));
//
////// tmp_uc = 2 * 3 * sizeof ( Complex );
////// EHlocal = (Complex **)malloc( (size_t)tmp_uc );
	// getEmptyVectors<Complex *>( 2, &EH, NULL );
	// getEmptyVectors<Complex>( 3, &( EH[0] ), NULL );
	// getEmptyVectors<Complex>( 3, &( EH[1] ), NULL );
	if (vectype)
	{
		EHlocal[0][2].Re = 0;
		EHlocal[0][2].Im = 0;
		ipow = 2 * params.n - params.m + 1;
		ipowodd = abs( ipow ) % 2;
		ipow2odd = (( abs( ipow ) / 2 ) % 2 ); // `abs' was ommited here for some reason, I cann't remember which. It should be `floor()' instead.
		phi = atan( y / x );
		ro2 = x * x + y * y;
		r2 = ro2 + z * z;
		r = sqrt( r2 );
		tmp1 = params.k0 * pow( params.w0, 2 ) / ( 2 * r2 );
		tmp2 = ro2 * params.k0 * tmp1; // CHECK THIS !!!!!! (see compiler warnings on -pedantic -Wextra -Wall)
		yztmp = y * z;
		imreCommon = -yztmp * tmp1 * pow( sqrt( tmp2 ), params.m ) * Lmn( params.m, params.n, tmp2) / ro2 * exp( -tmp2 / 2 ) * cos( params.m * phi );
//	// Optimized:
		if ( ipow > -1 )
		{
			itmp1 = -ipow2odd + 1 - ipow2odd; // itmp1 = ( ipow2odd * (-1) + 1 - ipow2odd ); // V
			mult.Re = ( 1 - ipowodd ) * itmp1;
			mult.Im = ipowodd * itmp1;
		}
		else
		{
			itmp2 = ( ipow2odd + ipow2odd - 1 );
			mult.Re = ( 1 - ipowodd ) * itmp2;
			mult.Im = ipowodd * itmp2;
		}
//	// Optimized.
		tmp3 = params.k0 * r;
		xyCommon.Re = imreCommon * cos( tmp3 );
		xyCommon.Im = imreCommon * sin( tmp3 );
		xyCommon = mult * xyCommon;
		EHlocal[0][0] = y * xyCommon;
		EHlocal[0][1] = -x * xyCommon;
		xyCommon = xyCommon / r * sqrt( params.epsilon / params.mu );
		EHlocal[1][0] = xyCommon * x * z;
		EHlocal[1][1] = xyCommon * yztmp;
		EHlocal[1][2] = xyCommon * (-ro2);
//	// DEBUG:
/*
		dtmp = 0;
*/
//	// DEBUG.
	}
	else
	{
//	// Same as above, see article source.
	}
	//return EH;
}
//
template <typename T> T *rotateXY( const double sinR, const double cosR, T x, T y, short ind, const char rMatrix [4][2][2] )
{
	static T *xy;
//
// #define A002             
// BGNINIT                  
// 	PPT(T), PPD(2), PPA(xy)
// ENDINIT
//
// printf(STRINGIZE_TOKEN(A002));
//
//
xy = (T *)malloc( 2 * sizeof ( T ) );
//
//
	// getEmptyVectors<T>( 2, &xy, NULL );
	xy[0] = x * cosR * rMatrix[ind][0][0] + y * sinR * rMatrix[ind][0][1];
	xy[1] = y * cosR * rMatrix[ind][1][1] - x * sinR * rMatrix[ind][1][0];
	return xy;
}
void doQWP( Complex *Xte,/* Complex *Xte45,*/ Complex *XteQWP, Complex **jQWP )
{
	int      i, j;
//Complex *tmp;
//
//tmp = rotateXY<Complex>( tri45, tri45, Xte[0], Xte[1], qRotation[0] ); // ????????
//Xte45[0] = tmp[0];
//Xte45[1] = tmp[1];
//XteQWP[2] = Xte[2];
//te45[0] = Xte[0] * tri225 + Xte[1] * tri225;
//te45[1] = Xte[1] * tri225 - Xte[0] * tri225;
	for ( i = 0; i < 2; i++ )
		for ( j = 0; j < 2; j++ )
		{
			//XteQWP[i] += (Xte[j]) * jQWP[i][j];// `+=' not working
			XteQWP[i] = XteQWP[i] + ( Xte[j] * jQWP[i][j] );
		}
	XteQWP[2] = Xte[2];
}
//
double *getLensSurfPointAbstract( double x, double y )
{
	double *xy;
//
// printf( 
// STRINGIZE_TOKEN( 
// BGNINIT                      
// 	PPT(double) PPD(2) PPA(xy) 
// ENDINIT 
// ) 
// );
//
//
xy = ( double * )malloc( 2 * sizeof ( double ) );
//
//
	// getEmptyVectors<double>( 2, &xy, NULL );
	xy[0] = x;
	xy[1] = y;
	return xy;
}
//
double getWarg( double lensF, double z, double NA )
{
	return -2 * M_PI / ( ( lensF + z ) * NA );
}
//
Complex *****ffft( double nt, double NA, double f, double z, int a, int b, LGBeam **iBeams, double * (*getLensSurfPoint)( double, double ) )
{
	int       i, j,/* k,*/ ii, jj,
	          M;
	double    QWParg,
	          phi, eta,
	          x, y, xLocal, yLocal,
	          tmp,
	         *xy,
	         *tmpRef;
	Complex **EHte, *****f_cplx,
	         *Ete,/* *Ete45,*/ *EteQWP/*,
	         *Hte, *Hte45, *HteQWP*/;
//
// #define PPTLA PPT(Complex)
// #define PPLGet() PPTLA PPD(3)
//
// #define A003()                           
// BGNINIT                                  
// 	PPTLA PPDL(a,2,3) PPDR(2,M) PPA(f_cplx)
// 	DDL(PPLGet,Ete,EteQWP)                 
// 	PPTLA PPD(2) PPA(jQWP)                 
// ENDINIT
//
//
f_cplx = ( Complex ***** )malloc( a * sizeof ( Complex **** ) );
Ete = ( Complex * )malloc( 3 * sizeof ( Complex ) );
EteQWP = ( Complex * )malloc( 3 * sizeof ( Complex ) );
//
//
// printf(STRINGIZE_TOKEN(A003()));
//
	// getEmptyVectors<Complex ****>( a, &f_cplx, NULL );
	// getEmptyVectors<Complex>( 3, &Ete,/* &Ete45,*/ &EteQWP, NULL );
	// getEmptyArrays<Complex>( 2, &jQWP, NULL );
	QWParg = M_PI / 4;
// This here below should be initializaed in the first plase:
	jQWP[0][0].Re = cos( QWParg );
	jQWP[0][0].Im = sin( QWParg );
	jQWP[1][1].Re = sin( QWParg );
	jQWP[1][1].Im = -cos( QWParg );
// This here below should be initializaed in the first plase.
	for ( ii = 0; ii < a; ii++ )
	{
		// getEmptyVectors<Complex ***>( 2, &( f_cplx[ii] ), NULL );
		// getEmptyVectors<Complex **>( 3, &( f_cplx[ii][0] ), &( f_cplx[ii][1] ), NULL );
		M = pow( 2, 1 + ceil( log( 2 * NA * NA * abs( z ) * iBeams[ii][0].k0 / ( 2 * M_PI * sqrt( nt * nt - NA * NA ))) / log( 2 )));

	f_cplx[i] = ( Complex **** )malloc( 2 * sizeof ( Complex *** ) );
	for ( j = 0; j < 2; j++  )
	{
		f_cplx[i][j] = ( Complex *** )malloc( 3 * sizeof ( Complex ** ) );
		for ( ii = 0; ii < 3; ii++ )
		{
			f_cplx[i][j][ii] = ( Complex ** )malloc( M * sizeof ( Complex * ) );
			for ( jj = 0; jj < M; jj++ )
			{
				f_cplx[i][j][ii][jj] = (Complex *)malloc( M * sizeof ( Complex ) );
			}
		}
	}


			// getEmptyArrays<Complex>( M, &( f_cplx[ii][0][i] ), &( f_cplx[ii][1][i] ), NULL );
		for ( i = 0; i < M; i++ )
			for ( j = 0; j < M; j++ )
			{
				phi = atan( i / j );//??????
				eta = asin( NA * sqrt( i * i + j * j ) / ( M * nt ) );
				tmp = f * tan( eta );
				x = cos( phi ) * tmp;
				y = sin( phi ) * tmp;
				xy = ( *getLensSurfPoint )( x, y );
				for ( jj = 0; jj < b; jj++ )
				{
					xLocal = xy[0] - iBeams[ii][jj].x0;
					yLocal = xy[1] - iBeams[ii][jj].y0;
					tmpRef = rotateXY<double>( sin225, cos225, xLocal, yLocal, (short)( iBeams[ii][jj].q ), qRotation );
					xLocal = tmpRef[0];
					yLocal = tmpRef[1];
					EHte = iBeams[ii][jj].getHE( TE, iBeams[ii][jj], xLocal, yLocal, z ); // HERE
					Ete = EHte[0];
					doQWP( Ete,/* Ete45,*/ EteQWP, jQWP );
//				doQWP( Hte,/* Hte45,*/ HteQWP );
					f_cplx[ii][0][0][i][j] = EteQWP[0];//actually this is a sum of Ete[jj], but in the case of thin Gaussian beam set it's negligeable
					f_cplx[ii][0][1][i][j] = EteQWP[1];//likewise
					f_cplx[ii][0][2][i][j] = EteQWP[2];//likewise
				}
			}
		doCZT( M, f_cplx[ii][0][0], getWarg( f, z, NA ), iBeams[ii][0].k0 );// OPTIMIZE all this and similar lines!!!
		doCZT( M, f_cplx[ii][0][1], getWarg( f, z, NA ), iBeams[ii][0].k0 );
		doCZT( M, f_cplx[ii][0][2], getWarg( f, z, NA ), iBeams[ii][0].k0 );
	}
	return f_cplx;
}
