#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "./fffc/libUtils.h"
/*
#define TEST_VALUE_TO_CHANGE 123
#define MACRO_PREPEND(a) ((TEST_VALUE_TO_CHANGE)+(a))
#define NAME1() 123455
#define CONCATMACROS(X,...) NAME ## X ()
#define SOMECHARS (U_LONG)&wewewewewe
*/
//
typedef U_SHORT GLbyte;
//
int main(int argc, char *argv[])
{
/*
	int value = 0;
//
	value = MACRO_PREPEND(value);
	printf( "\n%i\n", value );
	#undef TEST_VALUE_TO_CHANGE
	#define TEST_VALUE_TO_CHANGE 77
//
	printf(STRINGIZE_TOKEN(SOMECHARS));
	printf(STRINGIZE_TOKEN(CONCATMACROS(1,bla,bla,bla)));
*/
//
#define N 256
#define GRAPHNUM 13
GLbyte **graphContainer;
U_LONG NN;
L_DOUBLE *zmin, *zmax;
//
#define GetPPL() PPT(L_DOUBLE) PPD(GRAPHNUM)
//
/*
//
#define A001(PPL,...) PPL () PP_NARG(__VA_ARGS__)
#define A002 A001(GetPPL,zmin,zmax)
//
printf("\n\n%i\n\n", PP_NARG(GetPPL,zmin,zmax));
//
printf("\n\n");
printf(STRINGIZE_TOKEN(( A001(GetPPL,zmin,zmax) )));
printf("\n\n");
*/
//
#define A000    \
			BGNINIT                                                  \
				PPT(GLbyte) PPD(GRAPHNUM) PPD(N*N) PPA(graphContainer) \
				DDL(GetPPL,zmin,zmax)                                  \
			ENDINIT
//
printf("\n\n");
printf(STRINGIZE_TOKEN(A000));
printf("\n\n");
//
BGNINIT                                                  \
	PPT(GLbyte) PPD(GRAPHNUM) PPD(N*N) PPA(graphContainer) \
	DDL(GetPPL,zmin,zmax)                                  \
ENDINIT
NN = N*N;
for ( U_SHORT i = 0; i < GRAPHNUM; i++ )
	zmin[i] = i;
for ( U_SHORT i = 0; i < GRAPHNUM; i++ )
	zmax[i] = i;
int fff;
for ( U_SHORT i = 0; i < GRAPHNUM; i++ )
{
	fff = 4;
	for ( U_INT j = 0; j < NN; j++ )
		graphContainer[i][j] = 1234;
}
printf( "\n\n" );
//
	return 0;
}
