//
/****************************************************************************/
/*                                                                          */
/*                                                                          */
/*                                                                          */
/*                              libUtils.cpp                                */
/*                                                                          */
/*                                                                          */
/*                                                                          */
/****************************************************************************/
//
#include <stdlib.h>
/*
#include <stdarg.h>
*/
#include "./libUtils.h"
//
#define POLY_MAXTERMS 256
//
/*
struct polyCoordFrame
{
	short powConst,
	      powFrameRate,
	      powX,
	      powY,
	      powSin,
	      powCos;
}
polyCoordFrame operator* (const polyCoordFrame& x, const polyCoordFrame& y)
{
	polyCoordFrame res_pcf;
//
	res_pcf.powConst = x.powConst + y.powConst;
	res_pcf.powFrameRate = x.powFrameRate + y.powFrameRate;
	res_pcf.powU = x.powU + y.powU;
	res_pcf.powV = x.powV + y.powV;
	res_pcf.powV = x.powSin + y.powSin;
	res_pcf.powV = x.powCos + y.powCos;
//
	return res_pcf;
}
polyCoordFrame* operator+ (const polyCoordFrame& x, const polyCoordFrame& y)
{
	polyCoordFrame* res_pcfRef;
//
	res_pcfRef = new polyCoordFrame[POLY_MAXTERMS];
	if ( ( x.powConst == y.powConst ) &&  ( x.powFrameRate == y.powFrameRate ) && ( x.powX == y.powX ) && ( x.powY == y.powY ) && ( x.pow == y.pow ) )
	re
	c_res.Re = x.Re + y.Re;
	c_res.Im = x.Im + y.Im;
//
	return c_res;
}
*/
//
/*void getEmptyArrays( unsigned int typeSize, ... ) // [unsigned int typeSize; unsigned int [d1, d2, .., dN], 0; long address; ...]; unsigned int 0
{
	if ( atoi(getenv( "USER_GDB_DEBUG_THROW_EXPLICIT" )) > 0 )
		throw;

	unsigned char  dims;
	U_LONG          num,
	                  l[25];
	int               i,          j;
	U_LONG         size,     length,  ptrlen,
	              ioffs,      poffs, pstep,
	              vaarg,
	              varef;
	va_list          vl;
//
	va_start( vl, typeSize );
	size = typeSize;
	ptrlen = sizeof ( U_LONG * );
	do
	{
		i = 0;
		num = 1;
		l[i] = (U_INT)va_arg( vl, int );
		// if ( l[0] == 0 ) throw eid; // TODO: implement `try{}catch(eid){}' for this exception.
		while ( ( l[++i] = (U_INT)va_arg( vl, int ) ) != 0 )
		// {
			num *= l[i - 1];
			// i++;
		// }
		vaarg = (long)va_arg( vl, long ); // !!! `vaarg' is a pointer to pointer.
		dims = i - 1;
		( *( (U_LONG *)vaarg ) ) = (U_LONG)calloc( (( dims > 0 ) ? ptrlen * num : 0 ) + size * l[dims], 1 );
		varef = *( (U_LONG *)vaarg );
		if ( dims > 0 )
		{
			ioffs = 0;
			poffs = ptrlen * l[0];
			// istep = 1;
			pstep = (( dims > 1 ) ? ptrlen : size ) * l[1]; // HERE
			num = 1;
			for ( i = 0; i < dims; i++ )
			{
				length = num * l[i];
				for ( j = 0; j < length; j++ )
					( (U_LONG *)varef )[ioffs + j] = varef + poffs + j * pstep;
				num *= l[i];
				ioffs = poffs;
				poffs = poffs + ptrlen * length;
				pstep = (( i < dims - 1 ) ? ptrlen : size) * l[i + 1];
			}
		}
	}
	while ( ( size = (unsigned char)va_arg( vl, int )) != 0 );
	va_end( vl );

}
*/
//
