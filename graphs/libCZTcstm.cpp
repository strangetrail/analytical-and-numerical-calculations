/****************************************************************************/
/*                                                                          */
/*                                                                          */
/*                                                                          */
/*                            libCZTcstm.cpp                                */
/*                                                                          */
/*                                                                          */
/*                                                                          */
/****************************************************************************/
//
#include <stdlib.h>
//
// Remove this in final release
//
//include stdio_h
//
#include <math.h>
#include "./libFFTcstm.h"
#include "./libUtils.h"

//#define PERF_DOMAIN_CHECK

Complex **doCZT( int n, Complex **z_cplx, double w_arg, double k0)
{
	int i,j/*k,a,b*/;
	Complex   a_const,
	        **w_cplx,
	        **a_arg,
	        **a_fft,
	        **arg_ifft,
	        **cplx_res,
	        **w_res,
	        **final_res,
	        **w_domain,
	        **w_domainfft;
	double /* *cplx_arg, */ // Why this variable unused? Is this a bug?
	       tmp_1, tmp_2, tmp_3,
	       central_axis_delta =  0.5 - n / 2,
	       ddi, ddj;
//
/*
	Remove this after successfull testing:
	// #define PPL PPT(Complex),PPD(n),PPD(n),
	// getEmptyVectors<double>( 2, &cplx_arg, NULL );
	// getEmptyArrays<Complex>( n, &w_cplx, &a_arg, &arg_ifft, &w_res, &final_res, &w_domain, NULL );
*/
	// #define PPL PPT(Complex) PPDR(2,n)
// printf( 
// STRINGIZE_TOKEN( 
// 	BGNINIT                                                   
// 		PPT(double) PPD(2) PPA(cplx_arg)                        
// 		DDL(PPL,w_cplx,a_arg,arg_ifft,w_res,final_res,w_domain) 
// 	ENDINIT 
// ) 
// );
// 	#undef PPL
//
//
//cplx_arg = ( double * )malloc( 2 * sizeof ( double ) );// Why this variable unused? Is this a bug?
w_cplx = ( Complex ** )malloc( n * sizeof ( Complex * ) );
a_arg = ( Complex ** )malloc( n * sizeof ( Complex * ) );
arg_ifft = ( Complex ** )malloc( n * sizeof ( Complex * ) );
w_res = ( Complex ** )malloc( n * sizeof ( Complex * ) );
final_res = ( Complex ** )malloc( n * sizeof ( Complex * ) );
w_domain = ( Complex ** )malloc( n * sizeof ( Complex * ) );
for ( i = 0; i < n; i++ )
{
	w_cplx[i] = ( Complex * )malloc( n * sizeof ( Complex ) );
	a_arg[i] = ( Complex * )malloc( n * sizeof ( Complex ) );
	arg_ifft[i] = ( Complex * )malloc( n * sizeof ( Complex ) );
	w_res[i] = ( Complex * )malloc( n * sizeof ( Complex ) );
	final_res[i] = ( Complex * )malloc( n * sizeof ( Complex ) );
	w_domain[i] = ( Complex * )malloc( n * sizeof ( Complex ) );
}
//
//
	for ( i = 0; i < n; ddi = i + central_axis_delta, i++ )
		for ( j = 0; j < n; ddj = j + central_axis_delta, j++ )
		{
#ifndef PERF_DOMAIN_CHECK
			tmp_1 = -( ddi + ddj ) * k0;
			tmp_2 = ( ddi * ddi + ddj * ddj ) * w_arg / 2;
			w_res[i][j].Re = cos( tmp_2 );
			w_res[i][j].Im = sin( tmp_2 );
			a_const.Re = cos( tmp_1 );
			a_const.Im = sin( tmp_1 );
			a_arg[i][j] = z_cplx[i][j] * w_res[i][j] * a_const;
#endif
// Domain check:
//
//#ifdef PERF_DOMAIN_CHECK
			tmp_3 = w_arg * ( -( ddi * ddi + ddj * ddj ) / 2 );
			w_domain[i][j].Re = cos( tmp_3 );
			w_domain[i][j].Im = sin( tmp_3 );
//#endif
// Domain check^
//		w_cplx[i][j].Re = M_PI / ( w_arg / 2 ) * cos( M_PI * M_PI / ( 0.5 * w_arg ) * ( -( ddi * ddi + ddj * ddj )));
//		w_cplx[i][j].Im = M_PI / ( w_arg / 2 ) * sin( M_PI * M_PI / ( 0.5 * w_arg ) * ( -( ddi * ddi + ddj * ddj )));
		}
// Domain check:
//
		w_domainfft = doFFT( 0, n, w_domain, (char *)"output_fft1", NULL );
#ifdef PERF_DOMAIN_CHECK
		doFFT( 1, n, w_domainfft, (char *)"output_domain_ifft", NULL );
		printResults( (char *)"output_domain_direct", (char *)const_outFormatPrecision, 6, n, w_cplx );
#endif
// Domain check^
#ifndef PERF_DOMAIN_CHECK
		a_fft = doFFT( 0, n, a_arg, (char *)"output_fft2", NULL );
		for ( i = 0; i < n; i++ )
			for (j = 0; j < n; j++ )
				arg_ifft[i][j] = a_fft[i][j] * w_domainfft[i][j];
		cplx_res = doFFT( 1, n, arg_ifft, (char *)"output_ifft", NULL );
		for ( i = 0; i < n; i++ )
			for (j = 0; j < n; j++ )
				final_res[i][j] = cplx_res[i][j] * w_res[i][j];
/*
To apply the phase shift on infinitesimal propagation distance solve below equation for z_cplx and then recalculate distribution with same or larger z
final_res = doIFFT( doFFT( z_cplx * exp( sqrt( -1 ) * ( w_arg / 2 - k0 ) ) ) * doFFT( exp( ( - sqrt( -1 ) ) * w_arg / 2 ) ) ) * w_res;
*/
#endif
		return final_res;
}
