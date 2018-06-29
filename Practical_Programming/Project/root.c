
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>

int nat_log(const gsl_vector* t, void* params, gsl_vector* f){
	double x = *(double*)params;
	gsl_vector_set(f,0,exp(gsl_vector_get(t,0))-x);
return GSL_SUCCESS;
}

double mylog(double x){
assert(x>0);

/* reduce argument! */

const gsl_multiroot_fsolver_type* T = gsl_multiroot_fsolver_hybrids;
gsl_multiroot_fsolver* S = gsl_multiroot_fsolver_alloc(T,1);

gsl_multiroot_function F;
F.f = nat_log;
F.n = 1;
F.params = (void*) &x;

gsl_vector* start = gsl_vector_alloc(1);
gsl_vector_set(start,0,1);
gsl_multiroot_fsolver_set(S,&F,start);

double eps = 1e-3;
int status;

do{
	gsl_multiroot_fsolver_iterate(S);
	status = gsl_multiroot_test_residual(S->f,eps);
}while(status==GSL_CONTINUE);


double result = gsl_vector_get(S->x,0);
gsl_vector_free(start);
gsl_multiroot_fsolver_free(S);
return result;
}
