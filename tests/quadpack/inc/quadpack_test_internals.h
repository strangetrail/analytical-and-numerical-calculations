#ifndef launch_quadpack_test_internals_h__
#define launch_quadpack_test_internals_h__

typedef double ( *xy2d_function_t ) ( double , double );

struct VariablesAndFunction
{
	double y;
	xy2d_function_t function;
};

typedef struct VariablesAndFunction VariablesAndFunction_t;

#endif  // launch_quadpack_test_internals_h__
