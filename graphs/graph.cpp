

/* graph.cpp */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// #define GLEW_STATIC
// #define GL_GLEXT_PROTOTYPES
#include <GL/glew.h>
#include <GL/glut.h>
// #include <GL/freeglut.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <ft2build.h>
#include FT_FREETYPE_H
#include "./fffc/libUtils.h"
#include "./common/shader_utils.h"
#include "./fffc/libFFTcstm.h"
#include "./fffc/libCZTcstm.h"
#include "./fffc/fffc.h"
//
#define EPSILON_0 0.88541878176e-11L // F / m
/*
   Default values (1/3 of longest side of 1024x768 0.9e-6 m pitch SLM LC cell, gren DPSS, air)
*/
#define BEAMWAIST 0.306e-2L // m
#define WAVELENGTH 0.532e-6L // m
#define EPSILON_R 1.00058986L
#define DELTA_EPSILON_R 0.5e-6L
#define MU 1.25663753e-6 // H / m
//
struct point_text {
	GLfloat x;
	GLfloat y;
	GLfloat s;
	GLfloat t;
};
//aux:
GLint uniform_color_test;
//aux.
GLuint program;
GLint attribute_coord2d;
//axis:
// GLint attribute_coord2d_test;
// GLint uniform_transform;
GLint uniform_switch_transform;
//axis.
GLint uniform_vertex_transform;
GLint uniform_vertex_transform2;
GLint uniform_vertex_transform90;
GLint uniform_vertex_transform180_90;
GLint uniform_vertex_transform180_90_vert;
GLint uniform_vertex_transform_textX;
GLint uniform_texture_transform;
//axis:
GLint uniform_color;
// GLint uniform_color_red;
GLint uniform_switch_color;
//axis.
GLuint texture_id;
GLint uniform_mytexture;
//text:
GLint attribute_coord_text;
GLint uniform_tex;
GLint uniform_color_text;
//text.
FT_Library ft;
FT_Face face;
const char *fontfilename;
float offset_x = 0.0;
float offset_y = 0.0;
float offset_xx = 0;
float offset_yy = 0;
float scale = 1.0;
float scale_x = 1;
float scale_y = 1;
L_DOUBLE scale_math = 1.0;
L_DOUBLE default_beam_scale;
bool interpolate = false;
bool clamp = true;
bool rotate = false;
/* volatile */bool switch_scale = false;
bool polygonoffset = true;
/* volatile */bool switch_scale_do_once = false;
unsigned char datatype = 0;
		const char graphnames[13][30] =
		{
			"Magnetic field, X axis, Re",
			"Magnetic field, X axis, Im",
			"Magnetic field, Y axis, Re",
			"Magnetic field, Y axis, Im",
			"Magnetic field, Z axis, Re",
			"Magnetic field, Z axis, Im",
			"Electric field, X axis, Re",
			"Electric field, X axis, Im",
			"Electric field, Y axis, Re",
			"Electric field, Y axis, Im",
			"Electric field, Z axis, Re",
			"Electric field, Z axis, Im",
			"Energy Flux"
		};
		GLuint vbo[9];
		static GLbyte **graphContainer;
		static L_DOUBLE **mathGraphContainer;
		         /***znormRe, **znormIm,
		          *zmin, *zmax, *znorm*/
	static L_DOUBLE *zmin;
	static L_DOUBLE *zmax;
	static L_DOUBLE *znorm;
	static L_DOUBLE **znormRe;
	static L_DOUBLE **znormIm;
//
static Complex **EH;
struct point
{
	GLfloat x;
	GLfloat y;
};
const int border = 10;
const int ticksize = 10;
const short nticksmax = 40;
//
const L_DOUBLE ZMAXLIM = 0.1e+1L;
const L_DOUBLE ZMINLIM = -ZMAXLIM;
#define N 256
#define XYSCALE 20.0
#define GRAPHNUM 12
//
void normalize( L_DOUBLE &z, L_DOUBLE &min, L_DOUBLE &max )
{
	if ( ( z < ZMINLIM ) || ( z > ZMAXLIM ) )
	{
		if ( z < 0 )
			z = ZMINLIM;
		else
			z = ZMAXLIM;
	}
	if ( z < min )
		min = z;
	if ( z > max )
		max = z;
}
//
//
//
//
//
void limit( L_DOUBLE zmin_local, L_DOUBLE zmax_local, L_DOUBLE &z3 )
{
	if ( ( zmin_local < 0 ) && ( zmax_local > 0 ) )
	{
		if ( abs( zmin_local ) < abs( zmax_local ) )
		{
			z3 = zmax_local;
//printf( "z3 1max\n\n" );
		}
		else
		{
			z3 = -zmin_local;
//printf( "z3 1min\n\n" );
		}
	}
	else
	{
		if ( ( zmin_local < 0 ) && ( zmax_local <= 0 ) )
		{
			z3 = -zmin_local;
//printf( "z3 2min\n\n" );
		}
		else
		{
			if ( ( zmin_local >= 0 ) && ( zmax_local > 0 ) )
			{
				z3 = zmax_local;
//printf( "z3 3max\n\n" );
			}
			else
			{
				z3 = abs( zmax_local );
//printf( "z3 4other\n\n" );
			}
		}
	}
//printf( "z3 = %e\n\n", z3 );
}


// Reevaluating:


void reevaluate_graph ( L_DOUBLE **graph, short ind_graph )
{
	int indtmp1;
	// Removed due to preceding normalization and limitation of scaled values:
	/*
	double z1/_*, z2*_/;
	*/
	/// HERE
// Removed due to preceding normalization and limitation of scaled values:
////
/*
	for ( short i = 0; i < GRAPHNUM; i++ )
	{
		zmin[ind_graph] = ZMAXLIM;
		zmax[ind_graph] = ZMINLIM;
	}
*/
////
// Removed due to preceding normalization and limitation of scaled values.
// Removed due to preceding normalization and limitation of scaled values:
/*
	for ( int i = 0; i < N; i++ )
	{
		for ( int j = 0; j < N; j++ )
		{
			indtmp1 = N * i + j;
			z1 = graph[ind_graph    ][indtmp1] * scale_math;
			//if ( ind_graph < GRAPHNUM )
				//z2 = graph[ind_graph + 1][indtmp1] * scale_math;
			normalize( z1, zmin[ind_graph    ], zmax[ind_graph    ] );
			//if ( ind_graph < GRAPHNUM )
				//normalize( z2, zmin[ind_graph + 1], zmax[ind_graph + 1] );
			znormRe[ind_graph    ][indtmp1] = z1;
			//if ( ind_graph < GRAPHNUM )
				//znormIm[ind_graph + 1][indtmp1] = z2;
		}
	}
	//limit( zmin, zmax, ind_graph, z1 );
	// if ( ind_graph < GRAPHNUM )
		// limit( zmin, zmax, ind_graph + 1, z2 );
*/
// Removed due to preceding normalization and limitation of scaled values.
	for ( int i = 0; i < N; i++ )
		for ( int j = 0; j < N; j++ )
		{
			indtmp1 = N * i + j;
			/*
			Changed due to preceding normalization \
			 and limitation of scaled values:      \
			*/
			graphContainer[ind_graph][indtmp1] =                \
			 roundf ( 127 * /*( */( scale_math/*V*/             \
			                        * graph[ind_graph][indtmp1] \
			                      )/* / z )*/                   \
			          + 128                                     \
			        );                                           
			/*
			HERE: ambiguity in `[index]' & `[index + 1]' and in                 \
			 `znormRe' & `znormIm' notations - check the same issue below.      \
			 That's serious error. Also: save in graph[][] only znormX[] values \
			 and remove ZXXXLIM verification in `normalize' for better scaling  \
			 along Z axis.                                                       
			*/
			/*
			Changed due to preceding normalization and limitation of scaled values:
			*/
			// if ( ind_graph < GRAPHNUM )
				// graphContainer[ind_graph + 1][indtmp1] = roundf( ( znormIm[ind_graph + 1][indtmp1] / z2 ) * 127 + 128 );
		}
}
//
//
// Reevaluating.
//
//
int init_resources()
{
//data:
// U_SHORT tmp_uc;
	int        i,      j,       k,
					 ind, indtmp, indtmp1;
	L_DOUBLE   /*z,      d,*/          // 8 bytes, 2.2e-308L < x < 1.8e+308L
						x1,     y1,
// TEST:
						x2,     y2,
// TEST.
						z1,     z2,      z3;
// TEST:
	Complex ****f_cplx;
// TEST.
//
//// Moved to globals:
/*
	static L_DOUBLE *zmin;
	static L_DOUBLE *zmax;
	static L_DOUBLE *znorm;
	static L_DOUBLE *znormRe;
	static L_DOUBLE *znormIm;
*/
//// Moved to globals.
	// GLbyte graph[N][N];
//data.
//text:
//
/* Initialize the FreeType2 library */
if (FT_Init_FreeType(&ft)) {
fprintf(stderr, "Could not init freetype library\n");
return 0;
}
/* Load a font */
fontfilename = "FreeSans.ttf";
if (FT_New_Face(ft, fontfilename, 0, &face)) {
fprintf(stderr, "Could not open font %s\n", fontfilename);
return 0;
}
//
//text.
	program = create_program("graph.v.glsl", "graph.f.glsl");
	if (program == 0)
		return 0;
	attribute_coord2d = get_attrib(program, "coord2d");
//axis:
//
	uniform_switch_transform = get_uniform( program, "switch_transform" );
//
//axis.
	uniform_vertex_transform = get_uniform(program, "vertex_transform");
	uniform_vertex_transform2 = get_uniform(program, "vertex_transform2");
	uniform_vertex_transform90 = get_uniform(program, "vertex_transform90");
	uniform_vertex_transform180_90 = get_uniform(program, "vertex_transform180_90");
	uniform_vertex_transform180_90_vert = get_uniform(program, "vertex_transform180_90_vert");
	uniform_vertex_transform_textX = get_uniform(program, "vertex_transform_textX");
	uniform_texture_transform = get_uniform(program, "texture_transform");
	uniform_mytexture = get_uniform(program, "mytexture");
//axis:
//
	uniform_color = get_uniform(program, "color");
	uniform_switch_color = get_uniform( program, "switch_color" );
//
//axis.
//aux:
	uniform_color_test = get_uniform( program, "colortest" );
//aux.
//text:
//
	attribute_coord_text = get_attrib(program, "coordtext");
	uniform_tex = get_uniform(program, "tex");
	uniform_color_text = get_uniform(program, "colortext");
//
//text.
	if ( attribute_coord2d == -1                                              \
       || uniform_vertex_transform == -1 || uniform_vertex_transform2 == -1 \
       || uniform_vertex_transform90 == -1                                  \
       || uniform_vertex_transform180_90 == -1                              \
       || uniform_vertex_transform180_90_vert == -1                         \
       || uniform_vertex_transform_textX == -1                              \
       || uniform_texture_transform == -1 || uniform_mytexture == -1        \
       || uniform_switch_color == -1 || uniform_switch_transform == -1      \
       || uniform_color == -1 || uniform_color_test == -1                   \
       || attribute_coord_text == -1 || uniform_tex == -1                   \
       || uniform_color_text == -1                                          \
     )                                                                       
		return 0;
// Create our datapoints, store it as bytes
//
// 	 ttt = sizeof ( GLbyte * );
// 	 graphContainer = (GLbyte **)malloc( GRAPHNUM * ttt );
// 
// 	 #define GetPPL() PPT(L_DOUBLE) PPD(GRAPHNUM)
// printf( 
// STRINGIZE_TOKEN( 
// 	 BGNINIT                                                  
// 		 PPT(GLbyte) PPD(GRAPHNUM) PPD(N*N) PPA(graphContainer) 
// 		 DDL(GetPPL,zmin,zmax)                                  
// 	 ENDINIT 
// ) 
// );
//
// `alloc' replaced with c++ `new' operators:
//graphContainer = (GLbyte **)malloc( ( GRAPHNUM + 1 /* additional graph for testing */ ) * sizeof ( GLbyte * ) );
graphContainer = new GLbyte * [ GRAPHNUM + 1 /* additional graph for testing */ ];
//mathGraphContainer = (L_DOUBLE **)malloc( ( GRAPHNUM + 1 /* additional graph for testing */ ) * sizeof ( L_DOUBLE * ) );
mathGraphContainer = new L_DOUBLE * [ GRAPHNUM + 1 /* additional graph for testing */ ];
// Fixed for znormRe & znormIm and graph [i] & graph[i + 1] ambiguity:
//znormRe = (L_DOUBLE **)malloc( ( GRAPHNUM / 2 )/* V */ * sizeof ( L_DOUBLE * ) );
znormRe = new L_DOUBLE * [ GRAPHNUM / 2/* V */ ];
//znormIm = (L_DOUBLE **)malloc( ( GRAPHNUM / 2 )/* V */ * sizeof ( L_DOUBLE * ) );
znormIm = new L_DOUBLE * [ GRAPHNUM / 2/* V */ ];
// Fixed for znormRe & znormIm and graph [i] & graph[i + 1] ambiguity.
//zmin = (L_DOUBLE *)malloc( GRAPHNUM + 1 /* additional min-max pair for testing */ * sizeof (L_DOUBLE) );
zmin = new L_DOUBLE [ GRAPHNUM + 1 /* additional min-max pair for testing */ ];
//zmax = (L_DOUBLE *)malloc( GRAPHNUM + 1 /* additional min-max pair for testing */ * sizeof (L_DOUBLE) );
zmax = new L_DOUBLE [ GRAPHNUM + 1 /* additional min-max pair for testing */ ];
//znorm = (L_DOUBLE *)malloc( N*N * sizeof ( L_DOUBLE ) );
znorm = new L_DOUBLE [ N*N ];
// Fixed for znormRe & znormIm and graph [i] & graph[i + 1] ambiguity:
for ( i = 0; i < GRAPHNUM / 2; i++ )
{
	//graphContainer[i] = (GLbyte *)malloc( N*N * sizeof ( GLbyte ) );
	graphContainer[i] = new GLbyte [ N*N ];
	//mathGraphContainer[i] = (L_DOUBLE *)malloc( N*N * sizeof ( L_DOUBLE ) );
	mathGraphContainer[i] = new L_DOUBLE [ N*N ];
	//znormRe[i] = (L_DOUBLE *)malloc( N*N * sizeof ( L_DOUBLE ) );/* V */
	znormRe[i] = new L_DOUBLE [ N*N ];
	//znormIm[i] = (L_DOUBLE *)malloc( N*N * sizeof ( L_DOUBLE ) );/* V */
	znormIm[i] = new L_DOUBLE [ N*N ];
}
for ( i = GRAPHNUM / 2/*Remove this later*/; i < GRAPHNUM; i++ )
{
	//graphContainer[i] = (GLbyte *)malloc( N*N * sizeof ( GLbyte ) );
	graphContainer[i] = new GLbyte [ N*N ];
	//mathGraphContainer[i] = (L_DOUBLE *)malloc( N*N * sizeof ( L_DOUBLE ) );
	mathGraphContainer[i] = new L_DOUBLE [ N*N ];
}
// Fixed for znormRe & znormIm and graph [i] & graph[i + 1] ambiguity.
//graphContainer[GRAPHNUM] = (GLbyte *)malloc( N*N * sizeof ( GLbyte ) );/* additional graph for testing */
graphContainer[GRAPHNUM] = new GLbyte [ N*N ];
//mathGraphContainer[GRAPHNUM] = (L_DOUBLE *)malloc( N*N * sizeof ( L_DOUBLE ) );/* additional graph for testing */
mathGraphContainer[GRAPHNUM] = new L_DOUBLE [ N*N ];
//
//printf( "And now errors are going to happen ...\n\n" );
//
//EH = (Complex **)malloc( 2 * sizeof ( Complex * ) );
EH = new Complex * [ 2 ];
for (i = 0; i < 2; i++)
{
	//EH[i] = (Complex *)malloc( 3 * sizeof ( Complex ) );
	EH[i] = new Complex [ 3 ];
}
// TEST:
f_cplx = new Complex *** [2];// Number of fields.
for ( short i = 0; i < 2; i++ )
{
	f_cplx[i] = new Complex ** [3];// Number of components.
	for ( short j = 0; j < 3; j++ )
	{
		f_cplx[i][j] = new Complex * [N];
		for ( short k = 0; k < N; k++ )
		{
			f_cplx[i][j][k] = new Complex [N];
		}
	}
}
// TEST.
// `alloc' replaced with c++ `new' operators.
//printf( "\n\nAnd it's awesome - it didn't turn ugly.\n\n" );
	for ( i = 0; i < GRAPHNUM; i++ )
	{
		zmin[i] = ZMAXLIM;
		zmax[i] = ZMINLIM;
	}
//printf( "zmin=%e; zmax=%e.\n\n", zmin[0], zmax[0] );
	default_beam_scale = 1 / BEAMWAIST;
	LGBeam probe =
	{
		BEAMWAIST,
		0, 0, 0,
		2 * M_PI /  WAVELENGTH,
		MU, EPSILON_0 * ( EPSILON_R + 0.5 * DELTA_EPSILON_R ),
		4, 1,
		0,
		NULL
	};
//
	for ( i = 0; i < N; i++ )
	{
		for ( j = 0; j < N; j++ )
		{
			indtmp1 = N * i + j;
			float x = (i - N / 2) / (N / 2.0); // -1.0 <= x <= 1.0
			float y = (j - N / 2) / (N / 2.0); // -1.0 <= y <= 1.0
			x1 = x * 2 * BEAMWAIST;
			y1 = y * 2 * BEAMWAIST;
			// d = hypotf(x, y) * 4.0; // test
/*
			if ( y > 0.95 )
				z = -1;
			else
				if ( y < -0.95 )
					z = 1;
				else
					z = y; // (1 - d * d) * expf(d * d / -2.0); // test
*/
			z3 = 0.25 + Lmn( 0, 3, 2 * y + 1 ) / 2;
			normalize( z3, zmin[GRAPHNUM], zmax[GRAPHNUM] );
			//
			//// EH = probe.getHE( TE, probe, x1, y1, 30.0 );
			//
//		// TEST:
			x2 = x1;
			y2 = y1;
//		// TEST.
			double *tmpRef = rotateXY<double>( sin225, cos225, x1, y1, 1, qRotation );
			x1 = tmpRef[0];
			y1 = tmpRef[1];
			getHEvector( TE, probe, x1, y1, 30.0, EH );
			for ( k = 0; k < 2; k++ )
				for ( ind = 0; ind < 3; ind++ )
				{
					indtmp  = 6 * k + 2 * ind;
					z1 = EH[1/* k */][ind].Re;
					z2 = EH[1/* k */][ind].Im;
					normalize( z1, zmin[indtmp    ], zmax[indtmp    ] );
					normalize( z2, zmin[indtmp + 1], zmax[indtmp + 1] );
//				// Fixed for znormRe & znormIm and graph [i] & graph[i + 1] ambiguity:
					znormRe[3 * k + ind/* V */][indtmp1] = z1;/* V */
					znormIm[3 * k + ind/* V */][indtmp1] = z2;/* V */
//				// Fixed for znormRe & znormIm and graph [i] & graph[i + 1] ambiguity.
				}
//
//
//
//		// TEST:
			Complex **jQWP = new Complex * [2];
			for ( int iii = 0; iii < 2; iii++ )
				jQWP[iii] = new Complex [2];
			for ( short iii = 0; iii < 2; iii++ )
				for ( short jjj = 0; jjj < 2; jjj++ )
				{
					jQWP[iii][jjj].Re = 0;
					jQWP[iii][jjj].Im = 0;
				}
			double QWParg = M_PI / 4;
			jQWP[0][0].Re = cos( QWParg );
			jQWP[0][0].Im = sin( QWParg );
			jQWP[1][1].Re = sin( QWParg );
			jQWP[1][1].Im = -cos( QWParg );
			getHEvector( TE, probe, x2, y2, 30.0, EH );
			Complex *EteQWP = ( Complex * )malloc( 3 * sizeof ( Complex ) );
			for ( short iii = 0; iii < 3; iii++ )
			{
				EteQWP[iii].Re = 0;
				EteQWP[iii].Im = 0;
			}
			doQWP( EH[1], EteQWP, jQWP );
//
//
			z2 = EteQWP[0].Re;
			f_cplx[1][0][i][j].Re = z2;
			//normalize( z2, zmin[6], zmax[6] );// Not required here for CZT testing.
			//znormRe[3][indtmp1] = z2;
//
			z2 = EteQWP[1].Re;
			f_cplx[1][1][i][j].Re = z2;
			//normalize( z2, zmin[8], zmax[8] );// Not required here for CZT testing.
			//znormRe[4][indtmp1] = z2;
//
			z2 = EteQWP[2].Re;
			f_cplx[1][2][i][j].Re = z2;
			//normalize( z2, zmin[10], zmax[10] );// Not required here for CZT testing.
			//znormRe[5][indtmp1] = z2;
//
//
			z1 = EteQWP[0].Im;
			f_cplx[1][0][i][j].Im = z2;
			//normalize( z1, zmin[7], zmax[7] );// Not required here for CZT testing.
			//znormIm[3][indtmp1] = z1;
//
			z1 = EteQWP[1].Im;
			f_cplx[1][1][i][j].Im = z2;
			//normalize( z1, zmin[9], zmax[9] );// Not required here for CZT testing.
			//znormIm[4][indtmp1] = z1;
//
			z1 = EteQWP[2].Im;
			f_cplx[1][2][i][j].Im = z2;
			//normalize( z1, zmin[11], zmax[11] );// Not required here for CZT testing.
			//znormIm[5][indtmp1] = z1;
//
//
			free( EteQWP );
			for ( int iii = 0; iii < 2; iii++ )
				free( jQWP[iii] );
			free( jQWP );
//		// TEST.
//
//
//
			// graphContainer[1][indtmp1] = z;
			znorm[i * N + j] = z3;
			// graph[i][j] = roundf(z * 127 + 128);
			// graphContainer[0][i * N + j] = roundf( z3 * 127 + 128 );
		}
	}
//
//
// TEST:
// DEBUG:
//
//
const char *debugOutputFileName = "debugTestOutput_CZT",
           *debugFormatPrecision = "  % 1.4e %c %1.4ei";
      char *fileName,
           *fileFormat;
//
// DEBUG.
fileName = ( char * )malloc( 5 + strlen( debugOutputFileName ) );
fileFormat = (char *) debugFormatPrecision;
//
sprintf( fileName, "%s%s%s", debugOutputFileName, "0", ".txt" );
printResults( fileName, fileFormat, 6, N, f_cplx[1][0] );
doCZT( N, f_cplx[1][0], getWarg( 0.05, 0.035, 0.54041950027058415544 ), probe.k0 );
//
sprintf( fileName, "%s%s%s", debugOutputFileName, "1", ".txt" );
printResults( fileName, fileFormat, 6, N, f_cplx[1][1] );
doCZT( N, f_cplx[1][1], getWarg( 0.05, 0.035, 0.54041950027058415544 ), probe.k0 );
//
sprintf( fileName, "%s%s%s", debugOutputFileName, "2", ".txt" );
printResults( fileName, fileFormat, 6, N, f_cplx[1][2] );
doCZT( N, f_cplx[1][2], getWarg( 0.05, 0.035, 0.54041950027058415544 ), probe.k0 );
//
for ( i = 0; i < N; i++ )
{
	for ( j = 0; j < N; j++ )
	{
		indtmp1 = i * N + j;
//
//
		normalize( f_cplx[1][0][i][j].Re, zmin[6], zmax[6] );
		znormRe[3][indtmp1] = f_cplx[1][0][i][j].Re;
//
		normalize( f_cplx[1][1][i][j].Re, zmin[8], zmax[8] );
		znormRe[4][indtmp1] = f_cplx[1][1][i][j].Re;
//
		normalize( f_cplx[1][2][i][j].Re, zmin[10], zmax[10] );
		znormRe[5][indtmp1] = f_cplx[1][2][i][j].Re;
//
//
		normalize( f_cplx[1][0][i][j].Im, zmin[7], zmax[7] );
		znormIm[3][indtmp1] = f_cplx[1][0][i][j].Im;
//
		normalize( f_cplx[1][1][i][j].Im, zmin[9], zmax[9] );
		znormIm[4][indtmp1] = f_cplx[1][1][i][j].Im;
//
		normalize( f_cplx[1][2][i][j].Im, zmin[11], zmax[11] );
		znormIm[5][indtmp1] = f_cplx[1][2][i][j].Im;
//
//
	}
}
// TEST.
//
//
for ( k = 0; k < 2; k++ )
	for ( ind = 0; ind < 3; ind++ )
	{
		indtmp  = 6 * k + 2 * ind;
		limit( zmin[indtmp], zmax[indtmp], z1 );
		limit( zmin[indtmp + 1], zmax[indtmp + 1], z2 );
	}
//
/*limit_inline:*/
//
	if ( ( zmin[GRAPHNUM] < 0 ) && ( zmax[GRAPHNUM] > 0 ) )
	{
		if ( abs( zmin[GRAPHNUM] ) < abs( zmax[GRAPHNUM] ) )
		{
			z3 = zmax[GRAPHNUM];
//printf( "z3 1max\n\n" );
		}
		else
		{
			z3 = -zmin[GRAPHNUM];
//printf( "z3 1min\n\n" );
		}
	}
	else
	{
		if ( ( zmin[GRAPHNUM] < 0 ) && ( zmax[GRAPHNUM] <= 0 ) )
		{
			z3 = -zmin[GRAPHNUM];
//printf( "z3 2min\n\n" );
		}
		else
		{
			if ( ( zmin[GRAPHNUM] >= 0 ) && ( zmax[GRAPHNUM] > 0 ) )
			{
				z3 = zmax[GRAPHNUM];
//printf( "z3 3max\n\n" );
			}
			else
			{
				z3 = abs( zmax[GRAPHNUM] );
//printf( "z3 4other\n\n" );
			}
		}
	}
//printf( "z3 = %e\n\n", z3 );
//
/*limit_inline.*/
//
	for ( i = 0; i < N; i++ )
		for ( j = 0; j < N; j++ )
		{
			indtmp1 = N * i + j;
			graphContainer[GRAPHNUM][i * N + j] = roundf( ( znorm[i * N + j] / z3 ) * 127 + 128 ) ;
			for ( k = 0; k < 2; k++ )
				for ( ind = 0; ind < 3; ind++ )
				{
					indtmp  = 6 * k + 2 * ind;
//				// Fixed for znormRe & znormIm and graph [i] & graph[i + 1] ambiguity:
//				// Saving for further reevaluation:
					mathGraphContainer[indtmp    ][indtmp1] = znormRe[3 * k + ind/* V */][indtmp1] / z1;
					mathGraphContainer[indtmp + 1][indtmp1] = znormIm[3 * k + ind/* V */][indtmp1] / z2;
//				// Saving for further reevaluation.
					graphContainer[indtmp    ][indtmp1] = roundf( ( mathGraphContainer[indtmp    ][indtmp1]/* V *//*znormRe[3 * k + ind/_* V *_/][indtmp1] / z1*/ ) * 127 + 128 );
					graphContainer[indtmp + 1][indtmp1] = roundf( ( mathGraphContainer[indtmp + 1][indtmp1]/* V *//*znormIm[3 * k + ind/_* V *_/][indtmp1] / z2*/ ) * 127 + 128 );
//				// Fixed for znormRe & znormIm and graph [i] & graph[i + 1] ambiguity.
				}
		}
/**/
// `alloc' replaced with c++ `new' operators:
	for ( i = 0; i < 2; i++ )
	{
		//free( EH[i] );
		delete [] EH[i];
	}
	//free( EH );
	delete [] EH;
	for ( i = 0; i < GRAPHNUM / 2; i++ )
	{
		//free( znormIm[i] );
		delete [] znormIm[i];
	}
	//free( znormIm );
	delete [] znormIm;
	for ( i = 0; i < GRAPHNUM / 2; i++ )
	{
		//free( znormRe[i] );
		delete [] znormRe[i];
	}
	//free( znormRe );
	delete [] znormRe;
	//free( znorm );
	delete [] znorm;
	//free( zmax );
	delete [] zmax;
	//free( zmin );
	delete [] zmin;
for ( short i = 0; i < 2; i++ )
{
	for ( short j = 0; j < 3; j++ )
	{
		for ( short k = 0; k < N; k++ )
			delete [] f_cplx[i][j][k];
		delete [] f_cplx[i][j];
	}
	delete [] f_cplx[i];// Number of components.
}
delete [] f_cplx;// Number of fields.
// `alloc' replaced with c++ `new' operators.
	//printf( "\n\nAnd it's double awesome - it didn't turn ugly this time.\n\n" );
	/* Upload the texture with our datapoints */
	glActiveTexture(GL_TEXTURE0);
	glGenTextures(1, &texture_id);
	glBindTexture(GL_TEXTURE_2D, texture_id);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_LUMINANCE, N, N, 0, GL_LUMINANCE, GL_UNSIGNED_BYTE, graphContainer[0] );
//axis:
	// Create two vertex buffer objects // ... four vertex bufffer ...
	// 0 - active grid
	// 1 - active graph
	// 2 - text
	// 3 - border
	// 4 - axis
	// ...
	// 8 - meshes for graph
	glGenBuffers(9, vbo);
	// Create a VBO for the border
	static const point border[4] = { {-1, -1}, {1, -1}, {1, 1}, {-1, 1} };
	glBindBuffer(GL_ARRAY_BUFFER, vbo[3]);
	glBufferData(GL_ARRAY_BUFFER, sizeof border, border, GL_STATIC_DRAW);
//axis.
	// Create an array for 101 * 101 vertices
	glm::vec2 vertices[101][101];
	for ( i = 0; i < 101; i++ )
	{
		for ( j = 0; j < 101; j++ )
		{
			vertices[i][j].x = (j - 50) / 50.0;
			vertices[i][j].y = (i - 50) / 50.0;
		}
	}
	// Tell OpenGL to copy our array to the buffer objects
	glBindBuffer(GL_ARRAY_BUFFER, vbo[0]);
	glBufferData(GL_ARRAY_BUFFER, sizeof vertices, vertices, GL_STATIC_DRAW);
	// Create an array of indices into the vertex array that traces both horizontal and vertical lines
	GLushort indices[100 * 101 * 6];
	i = 0;
	for (int y = 0; y < 101; y++)
	{
		for (int x = 0; x < 100; x++)
		{
			indices[i++] = y * 101 + x;
			indices[i++] = y * 101 + x + 1;
		}
	}
	for (int x = 0; x < 101; x++)
	{
		for (int y = 0; y < 100; y++)
		{
			indices[i++] = y * 101 + x;
			indices[i++] = (y + 1) * 101 + x;
		}
	}
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo[1]);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, 100 * 101 * 4 * sizeof *indices, indices, GL_STATIC_DRAW);

	// Create another array of indices that describes all the triangles needed to create a completely filled surface
	i = 0;

	for (int y = 0; y < 101; y++) {
		for (int x = 0; x < 100; x++) {
			indices[i++] = y * 101 + x;
			indices[i++] = y * 101 + x + 1;
			indices[i++] = (y + 1) * 101 + x + 1;

			indices[i++] = y * 101 + x;
			indices[i++] = (y + 1) * 101 + x + 1;
			indices[i++] = (y + 1) * 101 + x;
		}
	}

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo[8]);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof indices, indices, GL_STATIC_DRAW);

	return 1;
}
//
// Create a projection matrix that has the same effect as glViewport().
// Optionally return scaling factors to easily convert normalized device coordinates to pixels.
//
glm::mat4 viewport_transform(float x, float y, float width, float height, float *pixel_x = 0, float *pixel_y = 0)
{
	// Map OpenGL coordinates (-1,-1) to window coordinates (x,y),
	// (1,1) to (x + width, y + height).
	// First, we need to know the real window size:
	//float window_width = glutGet(GLUT_WINDOW_WIDTH);// Unused.
	//float window_height = glutGet(GLUT_WINDOW_HEIGHT);// Unused.
	// Calculate how to translate the x and y coordinates:
	//float offset_xx = (2.0 * x + (width - window_width)) / window_width;// Unused.
	//float offset_yy = (2.0 * y + (height - window_height)) / window_height;// Unused.
	// Calculate how to rescale the x and y coordinates:
	//float scale_x = /*scale **/ width / window_width;// Unused.
	//float scale_y = /*scale **/ height / window_height;// Unused.
	// Calculate size of pixels in OpenGL coordinates
	if (pixel_x)
		*pixel_x = 2.0 / width;
	if (pixel_y)
		*pixel_y = 2.0 / height;
	// return glm::scale(glm::translate(glm::mat4(1), glm::vec3(offset_xx, offset_yy, 0)), glm::vec3(scale_x/*temp:*//*scale*//*temp.*/, scale_y/*temp:*//*scale*//*temp.*/, 1));
	return glm::scale(glm::translate(glm::mat4(1), glm::vec3(x, y, 0)), glm::vec3(width/*scale_x*//*temp:*//*scale*//*temp.*/, height/*scale_y*//*temp:*//*scale*//*temp.*/, 1));
}
//text:
/**
 * Render text using the currently loaded font and currently set font size.
 * Rendering starts at coordinates (x, y), z is always 0.
 * The pixel coordinates that the FreeType2 library uses are scaled by (sx, sy).
 */
void render_text(const char *text, float x, float y, float sx, float sy) {
	const char *p;
	FT_GlyphSlot g = face->glyph;

	/* Create a texture that will be used to hold one "glyph" */
	GLuint tex;

	glActiveTexture(GL_TEXTURE0);
	glGenTextures(1, &tex);
	glBindTexture(GL_TEXTURE_2D, tex);
	glUniform1i(uniform_tex, 0);

	/* We require 1 byte alignment when uploading texture data */
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

	/* Clamping to edges is important to prevent artifacts when scaling */
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

	/* Linear filtering usually looks best for text */
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

	/* Set up the VBO for our vertex data */
	glEnableVertexAttribArray(attribute_coord_text);
	glBindBuffer(GL_ARRAY_BUFFER, vbo[2]);
	glVertexAttribPointer(attribute_coord_text, 4, GL_FLOAT, GL_FALSE, 0, 0);

	/* Loop through all characters */
	for (p = text; *p; p++) {
		/* Try to load and render the character */
		if (FT_Load_Char(face, *p, FT_LOAD_RENDER))
			continue;

		/* Upload the "bitmap", which contains an 8-bit grayscale image, as an alpha texture */
		glTexImage2D(GL_TEXTURE_2D, 0, GL_ALPHA, g->bitmap.width, g->bitmap.rows, 0, GL_ALPHA, GL_UNSIGNED_BYTE, g->bitmap.buffer);

		/* Calculate the vertex and texture coordinates */
		float x2 = x + g->bitmap_left * sx;
		float y2 = -y - g->bitmap_top * sy;
		float w = g->bitmap.width * sx;
		float h = g->bitmap.rows * sy;

		point_text box[4] = {
			{x2, -y2, 0, 0},
			{x2 + w, -y2, 1, 0},
			{x2, -y2 - h, 0, 1},
			{x2 + w, -y2 - h, 1, 1},
		};

		/* Draw the character on the screen */
		glBufferData(GL_ARRAY_BUFFER, sizeof box, box, GL_DYNAMIC_DRAW);
		glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);

		/* Advance the cursor to the start of the next character */
		x += (g->advance.x >> 6) * sx;
		y += (g->advance.y >> 6) * sy;
	}

	glDisableVertexAttribArray(attribute_coord_text);
	glDeleteTextures(1, &tex);
}
//text.
void display()
{
//text:
	float sx = 1.0 / glutGet(GLUT_WINDOW_WIDTH);
	float sy = 1.0 / glutGet(GLUT_WINDOW_HEIGHT);
//text.
//axis:
	float window_width = glutGet(GLUT_WINDOW_WIDTH);// Was unused for some reason lately.
	float window_height = glutGet(GLUT_WINDOW_HEIGHT);// Was unused for some reason lately.
//axis.
//		//graph reevaluation:
			if ( switch_scale )
			{
				reevaluate_graph ( mathGraphContainer, datatype );

				// DEBUG:
				/*
				printf ( "Reevaluating with scale factor." );
				*/
				// DEBUG.

			}
//		//graph reevaluation.
//		// reset math scale:
			if ( switch_scale_do_once )
			{
				short datatype_prev = datatype > 0 ? datatype - 1 : GRAPHNUM;
				scale_math = 1.0;
				switch_scale = false; // Ensure that mathematical scaling is disabled.
				reevaluate_graph( mathGraphContainer, datatype_prev );
				switch_scale_do_once = false;
				printf( "Scaling switched to use 3D mesh sizes\n" );
			}
//		// reset math scale.
//		//graph selection:
			glTexImage2D( GL_TEXTURE_2D, 0, GL_LUMINANCE, N, N, 0, GL_LUMINANCE, GL_UNSIGNED_BYTE, graphContainer[datatype] );
//		//graph selection.
	glUseProgram(program);
//axis:
	// glClearColor(1, 1, 1, 1);
	// glClear(GL_COLOR_BUFFER_BIT);
//axis.
			// glClearColor(0.0, 0.0, 0.0, 0.0);
			// glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
			glUniform1i(uniform_mytexture, 0);
			glm::mat4 model;
			if (rotate)
				model = glm::rotate(glm::mat4(1.0f), float (glutGet(GLUT_ELAPSED_TIME) / 1000.0), glm::vec3(0.0f, 0.0f, 1.0f));
			else
				model = glm::mat4(1.0f);
			glm::mat4 view = glm::lookAt(glm::vec3(0.0, 2.1, 2.1), glm::vec3(0.0, 0.0, 0.0), glm::vec3(0.0, 0.0, 1.0));
			glm::mat4 projection = glm::perspective(45.0f, 1.0f * 640 / 480, 0.1f, 10.0f);
			glm::mat4 vertex_transform = projection * view * model;
			glm::mat4 texture_transform = glm::translate(glm::scale(glm::mat4(1.0f), glm::vec3(scale, scale, 1)), glm::vec3(offset_x, offset_y, 0));
			glUniformMatrix4fv(uniform_vertex_transform, 1, GL_FALSE, glm::value_ptr(vertex_transform));
//
//additional_transform_180_deg:
glm::mat4 model2;
model2 = glm::rotate( glm::mat4( 1.0f ), float ( 180.0/180.0 * M_PI ), glm::vec3( 0.0f, 0.0f, 1.0f ) );
glm::mat4 vertex_transform2 = vertex_transform * model2;
glUniformMatrix4fv( uniform_vertex_transform2, 1, GL_FALSE, glm::value_ptr( vertex_transform2 ) );
//additional_transform_180_deg.
//
//additional_transform_90_deg:
glm::mat4 model90;
model90 = glm::rotate( glm::mat4( 1.0f ), float ( 90.0/180.0 * M_PI ), glm::vec3( 1.0f, 0.0f, 0.0f ) );
glm::mat4 vertex_transform90 = vertex_transform * model90;
glUniformMatrix4fv( uniform_vertex_transform90, 1, GL_FALSE, glm::value_ptr( vertex_transform90 ) );
//additional_transform_90_deg.
//
//additional_transform_180_90_deg:
glm::mat4 vertex_transform180_90 = vertex_transform2 * model90;
glUniformMatrix4fv( uniform_vertex_transform180_90, 1, GL_FALSE, glm::value_ptr( vertex_transform180_90 ) );
//additional_transform_180_90_deg.
//
//additional_transform_180_90_deg_vert:
glm::mat4 model_vert;
model_vert = glm::rotate( glm::mat4( 1.0f ), float ( 180.0/180.0 * M_PI ), glm::vec3( 1.0f, 0.0f, 0.0f ) );
glm::mat4 vertex_transform180_90_vert = vertex_transform180_90 * model_vert;
glUniformMatrix4fv( uniform_vertex_transform180_90_vert, 1, GL_FALSE, glm::value_ptr( vertex_transform180_90_vert ) );
//additional_transform_180_90_deg_vert.
//
//additional_transform_textX:
glm::mat4 model_textX;
model_textX = glm::rotate( glm::mat4( 1.0f ), float ( 90.0/180.0 * M_PI ), glm::vec3( 0.0f, 0.0f, 1.0f ) );
glm::mat4 vertex_transform_textX = vertex_transform * model_textX;
glUniformMatrix4fv( uniform_vertex_transform_textX, 1, GL_FALSE, glm::value_ptr( vertex_transform_textX ) );
//additional_transform_textX.
//
			glUniformMatrix4fv(uniform_texture_transform, 1, GL_FALSE, glm::value_ptr(texture_transform));

			glClearColor( 0.0, 0.0, 0.0, 0.0 );
			glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

			/* Set texture wrapping mode */
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, clamp ? GL_CLAMP_TO_EDGE : GL_REPEAT);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, clamp ? GL_CLAMP_TO_EDGE : GL_REPEAT);
			/* Set texture interpolation mode */
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, interpolate ? GL_LINEAR : GL_NEAREST);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, interpolate ? GL_LINEAR : GL_NEAREST);

			GLfloat grey[4] = { 0.5, 0.5, 0.5, 1 };
			glUniform4fv( uniform_color, 1, grey );

			glEnable(GL_DEPTH_TEST);
			if ( polygonoffset )
			{
				glPolygonOffset( 1, 0 );
				glEnable( GL_POLYGON_OFFSET_FILL );
			}
			/* Draw the grid using the indices to our vertices using our vertex buffer objects */
			glEnableVertexAttribArray(attribute_coord2d);
			glBindBuffer(GL_ARRAY_BUFFER, vbo[0]);
			glVertexAttribPointer(attribute_coord2d, 2, GL_FLOAT, GL_FALSE, 0, 0);

			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo[8]);
			glDrawElements(GL_TRIANGLES, 100 * 100 * 6, GL_UNSIGNED_SHORT, 0);

			glPolygonOffset(0, 0);
			glDisable(GL_POLYGON_OFFSET_FILL);

			/* Draw the grid, very bright */
			GLfloat bright[4] = { 2, 2, 2, 1 };
			glUniform4fv(uniform_color, 1, bright);

			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo[1]);
			glDrawElements(GL_LINES, 100 * 101 * 4, GL_UNSIGNED_SHORT, 0);
			/* Stop using the vertex buffer object */
			glDisableVertexAttribArray(attribute_coord2d);
			glBindBuffer(GL_ARRAY_BUFFER, 0);
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
			// glClearColor(0.0, 0.0, 0.0, 0.0);
			// glClear(GL_COLOR_BUFFER_BIT);

//text:
//switch_transform_and_color:
//aux:
	//glUniform1i( uniform_switch_color, (GLint)3 );
//aux.
	glUniform1i( uniform_switch_transform, (GLint)2 );
	glUniform1i( uniform_switch_color, (GLint)2 );
//switch_transform_and_color:
	/* White background */
	// glClearColor(1, 1, 1, 1);
	// glClear(GL_COLOR_BUFFER_BIT);
	/* Enable blending, necessary for our alpha texture */
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	GLfloat white[4] = { 1, 1, 1, 1 };
	glUniform4fv(uniform_color_test, 1, white);
	//GLfloat red[4] = { 1, 0, 0, 1 };// Unused.
	//GLfloat transparent_green[4] = { 0, 1, 0, 0.5 };// Unused.
	/* Set font size to 48 pixels, color to white */
	FT_Set_Pixel_Sizes(face, 0, 48);
	glUniform4fv(uniform_color_text, 1, white);
	/* Effects of alignment */
	render_text("-1", -1/* + 2 * sx*/, -1 - face->glyph->bitmap.rows * sy - 1.0 / ( border + ticksize )/* + 2 * sy*/, sx, sy);
	render_text("0", 0/* + 2 * sx*/, -1 - face->glyph->bitmap.rows * sy - 1.0 / ( border + ticksize )/* + 2 * sy*/, sx, sy);
	render_text("Y", 0/* + 2 * sx*/, -1 - 2.1 * face->glyph->bitmap.rows * sy - 1.0 / ( border + ticksize )/* + 2 * sy*/, sx, sy);
	render_text("1", 1/* + 2 * sx*/, -1 - face->glyph->bitmap.rows * sy - 1.0 / ( border + ticksize )/* + 2 * sy*/, sx, sy);
//
// additional transform:
//
	glUniform1i( uniform_switch_transform, (GLint)7 );
	render_text("0", 0/* + 2 * sx*/, -1 - face->glyph->bitmap.rows * sy - 1.0 / ( border + ticksize )/* + 2 * sy*/, sx, sy);
	render_text("X", 0/* + 2 * sx*/, -1 - 2.1 * face->glyph->bitmap.rows * sy - 1.0 / ( border + ticksize )/* + 2 * sy*/, sx, sy);
	glUniform1i( uniform_switch_transform, (GLint)2 );
//
// additional transform.
//
/*
	render_text("The Misaligned Fox Jumps Over The Lazy Dog", -1 + 8.5 * sx, 1 - 100.5 * sy, sx, sy);
*/
	/* Scaling the texture versus changing the font size */
/*
	render_text("The Small Texture Scaled Fox Jumps Over The Lazy Dog", -1 + 8 * sx, 1 - 175 * sy, sx * 0.5, sy * 0.5);
	FT_Set_Pixel_Sizes(face, 0, 24);
	render_text("The Small Font Sized Fox Jumps Over The Lazy Dog", -1 + 8 * sx, 1 - 200 * sy, sx, sy);
	FT_Set_Pixel_Sizes(face, 0, 48);
	render_text("The Tiny Texture Scaled Fox Jumps Over The Lazy Dog", -1 + 8 * sx, 1 - 235 * sy, sx * 0.25, sy * 0.25);
	FT_Set_Pixel_Sizes(face, 0, 12);
	render_text("The Tiny Font Sized Fox Jumps Over The Lazy Dog", -1 + 8 * sx, 1 - 250 * sy, sx, sy);
	FT_Set_Pixel_Sizes(face, 0, 48);
*/
	/* Colors and transparency */
/*
	render_text("The Solid Black Fox Jumps Over The Lazy Dog", -1 + 8 * sx, 1 - 430 * sy, sx, sy);
	glUniform4fv(uniform_color_text, 1, red);
	render_text("The Solid Red Fox Jumps Over The Lazy Dog", -1 + 8 * sx, 1 - 330 * sy, sx, sy);
	render_text("The Solid Red Fox Jumps Over The Lazy Dog", -1 + 28 * sx, 1 - 450 * sy, sx, sy);
	glUniform4fv(uniform_color_text, 1, transparent_green);
	render_text("The Transparent Green Fox Jumps Over The Lazy Dog", -1 + 8 * sx, 1 - 380 * sy, sx, sy);
	render_text("The Transparent Green Fox Jumps Over The Lazy Dog", -1 + 18 * sx, 1 - 440 * sy, sx, sy);
*/
//switch_transform_and_color:
//aux:
	//glUniform1i( uniform_switch_color, (GLint)3 );
//aux.
	glUniform1i( uniform_switch_transform, (GLint)0 );
	glUniform1i( uniform_switch_color, (GLint)0 );
//switch_transform_and_color:
//text.
//axis:
//transformation:
	// glUniformMatrix4fv(uniform_vertex_transform, 1, GL_FALSE, glm::value_ptr(vertex_transform));
//transformation.
	glEnableVertexAttribArray(attribute_coord2d);
	// glEnableVertexAttribArray(attribute_coord2d_test);
	// glVertexAttribPointer(attribute_coord2d, 2, GL_FLOAT, GL_FALSE, 0, 0);
//switch_transform:
	glUniform1i( uniform_switch_transform, (GLint)1 );
//switch_transform.
//switch_color:
//aux:
	//glUniform1i( uniform_switch_color, (GLint)3 );
//aux.
	glUniform1i( uniform_switch_color, (GLint)1 );
//switch_color.
	float pixel_x, pixel_y;
	/* ---------------------------------------------------------------- */
	/* Draw the borders */
	// Calculate a transformation matrix that gives us the same normalized device coordinates as above
	//glm::mat4 transform =// Unused (get its value from the function below).
	viewport_transform(border + ticksize, border + ticksize, window_width - border * 2 - ticksize, window_height - border * 2 - ticksize, &pixel_x, &pixel_y);
	// Tell our vertex shader about it
	// glUniformMatrix4fv(uniform_vertex_transform, 1, GL_FALSE, glm::value_ptr(transform));
	glUniformMatrix4fv(uniform_vertex_transform, 1, GL_FALSE, glm::value_ptr(vertex_transform));
	// Set the color to black
	// GLfloat white[4] = { 1, 1, 1, 1 };
	glUniform4fv(uniform_color, 1, white);
	// GLfloat red[4] = { 1, 0, 0, 1 };
	// glUniform4fv(uniform_color_red, 1, red);
	// Draw a border around our graph
	glBindBuffer(GL_ARRAY_BUFFER, vbo[3]);
	// glVertexAttribPointer(attribute_coord2d_test, 2, GL_FLOAT, GL_FALSE, 0, 0);
	glVertexAttribPointer(attribute_coord2d, 2, GL_FLOAT, GL_FALSE, 0, 0);
	glDrawArrays(GL_LINE_LOOP, 0, 4);
	/* ---------------------------------------------------------------- */
	point ticks[82];
	// point *ticks;
	/* ---------------------------------------------------------------- */
	/* Draw the y tick marks */

	float tickspacing = 0.05 * powf(10, -floor(log10(scale_y)));	// desired space between ticks, in graph coordinates
	float bottom = -1.0 / scale_y - offset_yy;	// bottom edge, in graph coordinates
	float top = 1.0 / scale_y - offset_yy;	// top edge, in graph coordinates
	int bottom_i = ceil(bottom / tickspacing);	// index of bottomleft tick, counted from the origin
	int top_i = floor(top / tickspacing);	// index of topright tick, counted from the origin
	float rem = bottom_i * tickspacing - bottom;	// space between bottom edge of graph and the first tick

	float firsttick = -1.0 + rem * scale_y;	// first tick in device coordinates

	int nticks = top_i - bottom_i + 1;	// number of ticks to show

	if (nticks > 41)
		nticks = 41;	// should not happen

	for (int i = 0; i < nticks; i++) {
		float y = firsttick + i * tickspacing * scale_y;
		float tickscale = ((i + bottom_i) % 10) ? 0.5 : 1;
		ticks[i * 2].y = y;
		ticks[i * 2].x = -1;
		ticks[i * 2 + 1].y = y;
		ticks[i * 2 + 1].x = -1 - ticksize * tickscale * pixel_x;
	}
	glBindBuffer(GL_ARRAY_BUFFER, vbo[4]);
	glBufferData(GL_ARRAY_BUFFER, sizeof ticks, ticks, GL_DYNAMIC_DRAW);
	// glVertexAttribPointer(attribute_coord2d_test, 2, GL_FLOAT, GL_FALSE, 0, 0);
	glVertexAttribPointer(attribute_coord2d, 2, GL_FLOAT, GL_FALSE, 0, 0);
	glDrawArrays(GL_LINES, 0, nticks * 2);

//switch_transform:
	glUniform1i( uniform_switch_transform, (GLint)3 );
//switch_transform.

//additional_transform:
	glBindBuffer(GL_ARRAY_BUFFER, vbo[5]);
	glBufferData(GL_ARRAY_BUFFER, sizeof ticks, ticks, GL_DYNAMIC_DRAW);
	glVertexAttribPointer(attribute_coord2d, 2, GL_FLOAT, GL_FALSE, 0, 0);
	glDrawArrays(GL_LINES, 0, nticks * 2);
//additional_transform.

//switch_transform:
	// glUniform1i( uniform_switch_transform, (GLint)1 );
//switch_transform.

//switch_transform:
	glUniform1i( uniform_switch_transform, (GLint)4 );
//switch_transform.

//additional_transform:
	glBindBuffer(GL_ARRAY_BUFFER, vbo[6]);
	glBufferData(GL_ARRAY_BUFFER, sizeof ticks, ticks, GL_DYNAMIC_DRAW);
	glVertexAttribPointer(attribute_coord2d, 2, GL_FLOAT, GL_FALSE, 0, 0);
	glDrawArrays(GL_LINES, 0, nticks * 2);
//additional_transform.

//switch_transform:
	// glUniform1i( uniform_switch_transform, (GLint)1 );
//switch_transform.

//switch_transform:
	glUniform1i( uniform_switch_transform, (GLint)5 );
//switch_transform.

//additional_transform:
	glBindBuffer(GL_ARRAY_BUFFER, vbo[7]);
	glBufferData(GL_ARRAY_BUFFER, sizeof ticks, ticks, GL_DYNAMIC_DRAW);
	glVertexAttribPointer(attribute_coord2d, 2, GL_FLOAT, GL_FALSE, 0, 0);
	glDrawArrays(GL_LINES, 0, nticks * 2);
//additional_transform.

//switch_transform:
	glUniform1i( uniform_switch_transform, (GLint)1 );
//switch_transform.

	/* ---------------------------------------------------------------- */
	/* Draw the x tick marks */
	tickspacing = 0.05 * powf(10, -floor(log10(scale_x)));	// desired space between ticks, in graph coordinates
	float left = -1.0 / scale_x - offset_xx;	// left edge, in graph coordinates
	float right = 1.0 / scale_x - offset_xx;	// right edge, in graph coordinates
	int left_i = ceil(left / tickspacing);	// index of left tick, counted from the origin
	int right_i = floor(right / tickspacing);	// index of right tick, counted from the origin
	rem = left_i * tickspacing - left;	// space between left edge of graph and the first tick

	firsttick = -1.0 + rem * scale_x;	// first tick in device coordinates

	nticks = right_i - left_i + 1;	// number of ticks to show

	if (nticks > 41)
		nticks = 41;	// should not happen

	for (int i = 0; i < nticks; i++) {
		float x = firsttick + i * tickspacing * scale_x;
		float tickscale = ((i + left_i) % 10) ? 0.5 : 1;

		ticks[i * 2].x = x;
		ticks[i * 2].y = -1;
		ticks[i * 2 + 1].x = x;
		ticks[i * 2 + 1].y = -1 - ticksize * tickscale * pixel_y;
	}
	glBindBuffer(GL_ARRAY_BUFFER, vbo[4]);
	glBufferData(GL_ARRAY_BUFFER, sizeof ticks, ticks, GL_DYNAMIC_DRAW);
	// glVertexAttribPointer(attribute_coord2d_test, 2, GL_FLOAT, GL_FALSE, 0, 0);
	glVertexAttribPointer(attribute_coord2d, 2, GL_FLOAT, GL_FALSE, 0, 0);
	glDrawArrays(GL_LINES, 0, nticks * 2);

//switch_transform:
	glUniform1i( uniform_switch_transform, (GLint)3 );
//switch_transform.
//additional_transform:
	glBindBuffer(GL_ARRAY_BUFFER, vbo[5]);
	glBufferData(GL_ARRAY_BUFFER, sizeof ticks, ticks, GL_DYNAMIC_DRAW);
	// glVertexAttribPointer(attribute_coord2d_test, 2, GL_FLOAT, GL_FALSE, 0, 0);
	glVertexAttribPointer(attribute_coord2d, 2, GL_FLOAT, GL_FALSE, 0, 0);
	glDrawArrays(GL_LINES, 0, nticks * 2);
//additional_transform.

//switch_transform:
	// glUniform1i( uniform_switch_transform, (GLint)1 );
//switch_transform.

//switch_transform:
	glUniform1i( uniform_switch_transform, (GLint)4 );
//switch_transform.

//additional_transform:
	glBindBuffer(GL_ARRAY_BUFFER, vbo[6]);
	glBufferData(GL_ARRAY_BUFFER, sizeof ticks, ticks, GL_DYNAMIC_DRAW);
	// glVertexAttribPointer(attribute_coord2d_test, 2, GL_FLOAT, GL_FALSE, 0, 0);
	glVertexAttribPointer(attribute_coord2d, 2, GL_FLOAT, GL_FALSE, 0, 0);
	glDrawArrays(GL_LINES, 0, nticks * 2);
//additional_transform.

//switch_transform:
	// glUniform1i( uniform_switch_transform, (GLint)1 );
//switch_transform.

//switch_transform:
	glUniform1i( uniform_switch_transform, (GLint)6 );
//switch_transform.

//additional_transform:
	glBindBuffer(GL_ARRAY_BUFFER, vbo[7]);
	glBufferData(GL_ARRAY_BUFFER, sizeof ticks, ticks, GL_DYNAMIC_DRAW);
	// glVertexAttribPointer(attribute_coord2d_test, 2, GL_FLOAT, GL_FALSE, 0, 0);
	glVertexAttribPointer(attribute_coord2d, 2, GL_FLOAT, GL_FALSE, 0, 0);
	glDrawArrays(GL_LINES, 0, nticks * 2);
//additional_transform.

//switch_transform:
	glUniform1i( uniform_switch_transform, (GLint)1 );
//switch_transform.

	// And we are done.

//switch_color:

//aux:
	//glUniform1i( uniform_switch_color, (GLint)3 );
//aux.

glUniform1i( uniform_switch_color, (GLint)0 );
//switch_color.

// glDisableVertexAttribArray(attribute_coord2d_test);
	glDisableVertexAttribArray(attribute_coord2d);

//switch_transform:
	glUniform1i( uniform_switch_transform, (GLint)0 );
//switch_transform.

//axis.

			glutSwapBuffers();
		}

		void special(int key, int x, int y)
		{
//	// TEMPORARY:
		if ( x == 0 )
			printf( "\n" );
		if ( y == 0 )
			printf( "\n" );
//	// TEMPORARY.
		switch ( key )
		{
			case GLUT_KEY_F1:
				interpolate = !interpolate;
				printf("Interpolation is now %s\n", interpolate ? "on" : "off");
			break;
		case GLUT_KEY_F2:
		clamp = !clamp;
		printf("Clamping is now %s\n", clamp ? "on" : "off");
		break;
		case GLUT_KEY_F3:
		rotate = !rotate;
		printf("Rotation is now %s\n", rotate ? "on" : "off");
		break;
		case GLUT_KEY_F4:
			if ( switch_scale )
			{
				switch_scale = false;
				switch_scale_do_once = true;
			}
			datatype = ( datatype + 1 ) % ( GRAPHNUM + 1 /* additional item for test purposes */ );
			printf( "Displaying graph for %s\n", graphnames[datatype] );
		break;
		case GLUT_KEY_F5:
			switch_scale = !switch_scale;
			printf( "Scaling switched to use %s\n", switch_scale ? "mathematical sizes" : "3D mesh sizes" );
		break;
		case GLUT_KEY_F6:
			polygonoffset = !polygonoffset;
			printf("Polygon offset is now %s\n", polygonoffset ? "on" : "off");
		break;
		case GLUT_KEY_LEFT:
		offset_x -= 0.03;
		break;
		case GLUT_KEY_RIGHT:
		offset_x += 0.03;
		break;
		case GLUT_KEY_UP:
		offset_y += 0.03;
		break;
		case GLUT_KEY_DOWN:
		offset_y -= 0.03;
		break;
		case GLUT_KEY_PAGE_UP:
			if ( switch_scale )
			{
				scale_math *= 1.5;
				printf( "Z axis scale now equal to %e\n", scale_math );
			}
			else
				scale *= 1.5;
			// scale_x *= scale;
			// scale_y *= scale;
		break;
		case GLUT_KEY_PAGE_DOWN:
			if ( switch_scale )
			{
				scale_math /= 1.5;  /// HERE
				printf( "Z axis scale now equal to %e\n", scale_math );
			}
			else
				scale /= 1.5;
			// scale_x *= scale;
			// scale_y *= scale;
		break;
		case GLUT_KEY_HOME:
		offset_x = 0.0;
		offset_y = 0.0;
		scale = 1.0;
		scale_math = 1.0;
		break;
		}
	glutPostRedisplay();
    }
//
void free_resources()
{
	glDeleteProgram( program );
	//free( zmax );
	//delete [] zmax;
	//free( zmin );
	//delete [] zmin;
	//free( mathGraphContainer );
	delete [] mathGraphContainer;
	//free( graphContainer );
	delete [] graphContainer;
}
//
    int main(int argc, char *argv[]) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
    glutInitWindowSize(640, 480);
    glutCreateWindow("");
    GLenum glew_status = glewInit();
    if (GLEW_OK != glew_status) {
    fprintf(stderr, "Error: %s\n", glewGetErrorString(glew_status));
    return 1;
    }
    if (!GLEW_VERSION_2_0) {
    fprintf(stderr, "No support for OpenGL 2.0 found\n");
    return 1;
    }
    GLint max_units;
    glGetIntegerv(GL_MAX_VERTEX_TEXTURE_IMAGE_UNITS, &max_units);
    if (max_units < 1) {
    fprintf(stderr, "Your GPU does not have any vertex texture image units\n");
    return 1;
    }
    printf("Use left/right/up/down to move.\n");
    printf("Use pageup/pagedown to change the horizontal scale.\n");
    printf("Press home to reset the position and scale.\n");
    printf("Press F1 to toggle interpolation.\n");
    printf("Press F2 to toggle clamping.\n");
    printf("Press F3 to toggle rotation.\n");
    printf("Press F4 to switch graph.\n");
    if (init_resources()) {
    glutDisplayFunc(display);
    glutIdleFunc(display);
    glutSpecialFunc(special);
    glutMainLoop();
    }
    free_resources();
    return 0;
    }
