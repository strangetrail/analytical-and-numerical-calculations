typedef unsigned char UCHAR;
struct LGBeam;
enum TETM
{
	TM = 0,
	TE = 1
};
struct LGBeam
{
	double w0,
	       x0,
	       y0,
	       z0,
	       k0,
	       mu,
	       epsilon;
	int    m,
	       n;
	UCHAR  q;
	Complex ** (*getHE)( TETM, LGBeam, double, double, double );
};
Complex *****ffft( double nt, double NA, double f, double z, int a, int b, LGBeam **iBeams, double * (*getLensSurfPoint)( double, double ) );
/*
  Export these functions nexported for debug purposes only
*/
const char qRotation[4][2][2] =
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
double Lmn( int alpha, int n, double x );
void doQWP( Complex *Xte,/* Complex *Xte45,*/ Complex *XteQWP, Complex **jQWP );
void getHEvector( TETM vectype, LGBeam params, double x, double y, double z, Complex **EHlocal );
template <typename T> T *rotateXY( const double sinR, const double cosR, T x, T y, short ind, const char rMatrix [4][2][2] );
double getWarg( double lensF, double z, double NA );
//rotateXY<double>( double sinR, double cosR, double x, double y, char **rMatrix );
