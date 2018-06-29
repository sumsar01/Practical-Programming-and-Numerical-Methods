#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

double Fe(double e, double r);

int shooter(const gsl_vector* x, void* params, gsl_vector* f){
	gsl_vector_set(f, 0, Fe(gsl_vector_get(x, 0), 8.));
	return GSL_SUCCESS;
}

int main(){

FILE* file;
file = fopen("data.dat","w");


const gsl_multiroot_fsolver_type* T = gsl_multiroot_fsolver_hybrids;
gsl_multiroot_fsolver* S = gsl_multiroot_fsolver_alloc(T,1);

gsl_multiroot_function F;
	F.f = shooter;
	F.n = 1;
	F.params = NULL;

gsl_vector *x = gsl_vector_alloc(1);
gsl_vector_set(x, 0, -1);

gsl_multiroot_fsolver_set(S, &F, x);
	double eps = 1e-6;
	int iterations = 0;
	int flag;

do{
	gsl_multiroot_fsolver_iterate(S);
	flag = gsl_multiroot_test_residual(S->f,eps);
	iterations++;
}while(flag == GSL_CONTINUE && iterations < 1000);

if(iterations == 1000){printf("Root did not converge after 1000 iterations.\n");}
else{printf("Root finder success.\n");}

double E = gsl_vector_get(S->x,0);


for(double r = 0; r<8+1e-5; r+=8.0/1000.){
	fprintf(file,"%g %g\n", r, Fe(E, r));
}
	fprintf(file,"\n\n");

for(double r = 0; r<8+1e-5; r+=8./20.){
	fprintf(file,"%g %g\n",r, r*exp(-r));
}

printf("The energy was found to be E = %g\n\n", E);

gsl_multiroot_fsolver_free(S);
gsl_vector_free(x);

}
