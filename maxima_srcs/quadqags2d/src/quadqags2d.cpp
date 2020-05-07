#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include "quadqags2d.hpp"
#include "quadqags2d_internals.hpp"
#include "integrand_export.hpp"


gsl_integration_workspace *giwQAGS2D;


double f_1d (double y, void *params);
double f_2d ( double x, void *params );


void allocGIW ()
{
	giwQAGS2D = gsl_integration_workspace_alloc (GIWSIZE);
}

void freeGIW ()
{
	gsl_integration_workspace_free ( giwQAGS2D );
}


double quad_qags_2d ( double l_lower_x, double l_lower_y, \
                      double l_upper_x, double l_upper_y, \
                      double *params,                     \
                      dd2d_fun_t integrand                \
                    )                                      
{
	double       result,
	             error;
	gsl_function F;
	VarsAndFun_t vafParams;


	vafParams.y = 0.0;
	vafParams.lower_boundary = l_lower_x;
	vafParams.upper_boundary = l_upper_x;
	vafParams.function = integrand;
	vafParams.params = params;

	F.function = &f_1d;
	F.params = &vafParams;

	gsl_integration_qag ( &F,                \
	                      l_lower_y,         \
	                      l_upper_y,         \
	                      EPSABS,            \
	                      EPSREL,            \
	                      NUMBER_OF_SAMPLES, \
	                      GSL_INTEG_TYPE,    \
	                      giwQAGS2D,         \
	                      &result,           \
	                      &error             \
	                    );                    


	return (result);
}

// Integration in 2D frame - additional free variables (e.g. `y')
//  and limits are passed through `params':
double f_2d ( double x, void *params )
{
	double y,
	       f_xy;
	VarsAndFun_t vafIntegrand ( *(VarsAndFun_t *)params );
	dd2d_fun_t integrand = vafIntegrand.function;


	y = vafIntegrand.y;
	/*
	f_xy = (*integrand) (x, y, vafIntegrand.params);
	*/
	f_xy = (double) integrand2d ( x, y,                   \
	                     vafIntegrand.params[0], \
	                     vafIntegrand.params[1], \
	                     vafIntegrand.params[2], \
	                     vafIntegrand.params[3], \
	                     vafIntegrand.params[4]  \
	                   );                         


	return f_xy;
}

// Integration in 1D frame, The variable `params' contains
//  only integrand function reference and limits:
double f_1d (double y, void *params)
{
	double       result,
	             error;
	gsl_function F_xy;
	// TODO : Does it really needs to implement "move constructor" concept
	//         for this?
	VarsAndFun_t vafIntegrand ( *(VarsAndFun_t *)params );


	vafIntegrand.y = y;

	F_xy.function = &f_2d;
	F_xy.params = &vafIntegrand;

	gsl_integration_qag ( &F_xy,                       \
	                      vafIntegrand.lower_boundary, \
	                      vafIntegrand.upper_boundary, \
	                      EPSABS,                      \
	                      EPSREL,                      \
	                      NUMBER_OF_SAMPLES,           \
	                      GSL_INTEG_TYPE,              \
	                      giwQAGS2D,                   \
	                      &result,                     \
	                      &error                       \
	                    );                              

	return result;
}
