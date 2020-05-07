#ifndef quadqags2d_internals_hpp__
#define quadqags2d_internals_hpp__


#include <math.h>
#include <gsl/gsl_integration.h>
#include "quadqags2d.hpp"


// The volume of GSL integration workspace:
#define GIWSIZE 20000
#define NUMBER_OF_SAMPLES 1000
#define GSL_INTEG_TYPE GSL_INTEG_GAUSS61
// Desired relative integration error:
#define EPSREL 1.0e-5
// Desired absolute integration error:
#define EPSABS 0.0


enum IntegrandParameters
{
	THEOMEGA,
	ZZ,
	TT,
	THEPHI0,
	THEALPHAEULER
};


class VariablesAndFunction
{
	public:
		double  y,
		        lower_boundary,
		        upper_boundary,
		       *params;

		dd2d_fun_t function;

		VariablesAndFunction ();
		VariablesAndFunction ( const VariablesAndFunction &vafReference );
};


typedef class VariablesAndFunction VarsAndFun_t;


#endif  // quadqags2d_internals_hpp__
