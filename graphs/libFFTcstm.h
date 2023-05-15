#include <stdarg.h>
template <typename T> struct xComplex
{
	T Re,
	  Im;
};
typedef xComplex<double> Complex;
typedef xComplex<int> iComplex;
void getEmptyArrays( unsigned char typeSize, ... );
Complex operator* (const iComplex& x, const Complex& y);
Complex operator* (const Complex& x, const iComplex& y);
Complex operator+ (const Complex& x, const Complex& y);
Complex operator* (const Complex& x, const Complex& y);
Complex operator* (const Complex& x, const double& y);
Complex operator* (const double& x, const Complex& y);
Complex operator/ (const Complex& x, const double& y);
Complex operator+= (const Complex& x, const Complex& y);
Complex **doFFT( short inv, int n, Complex **f_cplx, char *chrs_filename, char *expFormat );
// DEBUG:
void printResults( char *chrs_filename, char *outFormatPrecision, short ncols, int n, Complex **f_cplx );
// DEBUG.
