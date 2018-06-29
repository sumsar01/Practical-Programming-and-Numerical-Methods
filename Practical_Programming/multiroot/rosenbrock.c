#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int rosenbrock
(const gsl_vector*  x, void *params, gsl_vector* f){
	const double x0 = gsl_vector_get(x,0);
	const double x1 = gsl_vector_get(x,1);

	const double dfdx = -2*(1-x0) - 400*(x1-x0*x0);
	const double dfdy = 200*(x1 - x0*x0);
	gsl_vector_set(f,0,dfdx);
	gsl_vector_set(f,1,dfdy);
return GSL_SUCCESS;
}

int main(){

const gsl_multiroot_fsolver_type* T = gsl_multiroot_fsolver_hybrids;
gsl_multiroot_fsolver* S = gsl_multiroot_fsolver_alloc(T,2);

gsl_multiroot_function F;
F.f = rosenbrock;
F.n = 2;
F.params = NULL;

gsl_vector* start = gsl_vector_alloc(2);
gsl_vector_set(start,0,0);
gsl_vector_set(start,1,0);
gsl_multiroot_fsolver_set(S,&F,start);

int flag;
double eps = 1e-12;
int iterations = 0;

do{
	gsl_multiroot_fsolver_iterate(S);
	flag = gsl_multiroot_test_residual(S->f,eps);
	iterations++;
	printf("%g,%g\n", gsl_vector_get(S->x,0),gsl_vector_get(S->x,1));
}
while(flag == GSL_CONTINUE && iterations < 1000);

if(iterations == 1000){printf("Root did not converge after 1000 iterations.\n");}
else{printf("Root found after %i iterations.\n",iterations);
printf("Extremum found at (%g,%g).\n", gsl_vector_get(S->x,0),gsl_vector_get(S->x,1));}


gsl_vector_free(start);
gsl_multiroot_fsolver_free(S);

return 0;
}
