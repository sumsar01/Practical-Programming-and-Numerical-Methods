#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_errno.h>

double rosenbrock(const gsl_vector* x, void* params){
	const double x0 = gsl_vector_get(x,0);
	const double y0 = gsl_vector_get(x,1);
	return (1-x0)*(1-x0) + 100*(y0-x0*x0)*(y0-x0*x0);
}

int main(){

	const gsl_multimin_fminimizer_type* type =
		gsl_multimin_fminimizer_nmsimplex2;
	gsl_multimin_fminimizer* S = gsl_multimin_fminimizer_alloc(type,2);

	gsl_multimin_function F;
	F.f = &rosenbrock;
	F.n = 2;
	F.params = NULL;

	gsl_vector *x = gsl_vector_alloc(2);
	gsl_vector_set(x, 0, 0.5);
	gsl_vector_set(x, 1, 0.5);

	gsl_vector *step = gsl_vector_alloc(2);
	gsl_vector_set(step, 0, 0.01);
	gsl_vector_set(step, 1, 0.01);

	printf("Starting guess at (%g,%g)\n",gsl_vector_get(x,0), gsl_vector_get(x,1));

	gsl_multimin_fminimizer_set(S, &F, x, step);

	double eps = 1e-6;
	int flag;
	int iterations = 0;

	do{
		gsl_multimin_fminimizer_iterate(S);
		flag = gsl_multimin_test_size(S->size,eps);
		iterations++;
		printf("(%g,%g)\n", gsl_vector_get(S->x,0), gsl_vector_get(S->x,1));

	}while(flag == GSL_CONTINUE && iterations < 1000);

	if(iterations == 1000){printf("Root did not converge after 1000 iterations.\n");}
	else{printf("Root found after %i iterations.\n",iterations);
	printf("Extremum found at (%g,%g).\n", gsl_vector_get(S->x,0),gsl_vector_get(S->x,1));}

	gsl_multimin_fminimizer_free(S);
	gsl_vector_free(x);

return 0;
}
