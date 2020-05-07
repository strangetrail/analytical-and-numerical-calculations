#ifndef quadqags2d_hpp__
#define quadqags2d_hpp__


#include <gsl/gsl_integration.h>


typedef double ( *dd2d_fun_t ) ( double , double, double * );


extern "C" gsl_integration_workspace *giwQAGS2D;


extern "C" void allocGIW ();
extern "C" void freeGIW ();

extern "C" double quad_qags_2d ( double l_lower_x, double l_lower_y, \
                                 double l_upper_x, double l_upper_y, \
                                 double *params,                     \
                                 dd2d_fun_t integrand                \
                               );                                     


#endif  // quadqags2d_hpp__
