#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include "quadpack_test.h"
#include "quadpack_test_internals.h"


gsl_integration_workspace *w;

double f_2d (double x, void *params)
{
	VariablesAndFunction_t vafIntegrand = *(VariablesAndFunction_t *)params;
	xy2d_function_t integrand = vafIntegrand.function;
	double y = vafIntegrand.y;
	/*double f_xy = 1/sqrt(1+pow(x,2)+pow(y,2));*/
	double f_xy = (*integrand) (x, y);
	return (f_xy);
}

double f (double y, void *params)
{
	VariablesAndFunction_t vafParameters = { y, (xy2d_function_t)params };
	gsl_function F_xy;
	F_xy.function = &f_2d;
	F_xy.params = &vafParameters;
	double result, error;
	gsl_integration_qags(&F_xy, -1, 1, 0, 1e-7, 1000, w, &result, &error);
	return result;
}

int launch_quadpack_test ( const int argc, const unsigned char *argv, double * points, xy2d_function_t integrand)
{
	(void)argc; (void)argv;
	double   result, error,
	         expected = -4.0;

	xy2d_function_t params = integrand;

	gsl_function F;
	w = gsl_integration_workspace_alloc (1000);

	F.function = &f;
	F.params = &params;

	printf ("The length of C float = %ld\n", sizeof(float));
	printf ("The length of C double float = %ld\n", sizeof(double));
	for (int i = 0; i < 500; i++)
	{
		gsl_integration_qags (&F, -1, 1, 0, 1e-7, 1000, w, &result, &error);
		if ( i < 4 )
		{
			printf ( "point[%d] = % .18f\n", i, points[i] );
			points[i] = result;
		}
	}

	printf ("result = % .18f\n", result);
	printf ("exact result = % .18f\n", expected);
	printf ("estimated error = % .18f\n", error);
	printf ("actual error = % .18f\n", result - expected);
	printf ("intervals = %zu\n", w->size);

	gsl_integration_workspace_free (w);

	return (0);
}
