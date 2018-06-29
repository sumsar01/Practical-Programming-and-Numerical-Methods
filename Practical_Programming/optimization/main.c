#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_errno.h>
#include <assert.h>

struct experimental_data {int n; double* t, *y, *e;};

double deviation_from_model(const gsl_vector *x, void *params){
	double A = gsl_vector_get(x,0);
	double T = gsl_vector_get(x,1);
	double B = gsl_vector_get(x,2);
	struct experimental_data *p = (struct experimental_data*) params;
	int 	n = p->n;
	double *t = p->t;
	double *y = p->y;
	double *e = p->e;
	double sum = 0;
	#define MODEL(t) A*exp(-(t)/T) + B
	for(int i = 0; i<n;i++) sum+=pow( (MODEL(t[i])-y[i]) /e[i],2);
	return sum;
}

int main(){

double t[] = {0.47,1.41,2.36,3.30,4.24,5.18,6.13,7.07,8.01,8.95};
double y[] = {5.49,4.08,3.54,2.61,2.09,1.91,1.55,1.47,1.45,1.25};
double e[] = {0.26,0.12,0.27,0.10,0.15,0.11,0.13,0.07,0.15,0.09};

	int n = sizeof(t)/sizeof(t[0]);

	printf("#  t[i]  y[i]  e[i]\n");
	for(int i=0;i<n;i++) printf("%g %g %g\n",t[i],y[i],e[i]);
	printf("\n\n");

	struct experimental_data params;
	params.n = n;
	params.t = t;
	params.y = y;
	params.e = e;

	gsl_multimin_function F;
	F.f = deviation_from_model;
	F.n = 3;
	F.params = (void*) &params;

	gsl_vector *start = gsl_vector_alloc(F.n);
	gsl_vector_set(start,0,3);
	gsl_vector_set(start,1,2);
	gsl_vector_set(start,2,1);

	gsl_vector *step = gsl_vector_alloc(F.n);
	gsl_vector_set(step,0,2);
	gsl_vector_set(step,1,2);
	gsl_vector_set(step,2,2);

	gsl_multimin_fminimizer* state =
		gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2,F.n);
	gsl_multimin_fminimizer_set(state, &F, start, step);

	int iteration = 0, status;
	double acc = 0.01;

	do{
		iteration++;
		int iteration_status = gsl_multimin_fminimizer_iterate(state);

		if(iteration_status != 0){
			fprintf(stderr,"unable to improve\n");
			break;
			}

		status = gsl_multimin_test_size(state->size,acc);
		if(status == GSL_SUCCESS) fprintf(stderr,"converged\n");
		fprintf(stderr,"iteration = %3i   ",iteration);
		fprintf(stderr,"A = %.3f   ",gsl_vector_get(state->x,0));
		fprintf(stderr,"T = %.3f   ",gsl_vector_get(state->x,1));
		fprintf(stderr,"B = %.3f   ",gsl_vector_get(state->x,2));
		fprintf(stderr,"size = %g   ",state->size);
		fprintf(stderr,"fval/n = %.3f   ",state->fval/n);
		fprintf(stderr,"\n");

	}while( status == GSL_CONTINUE && iteration < 100);

	double A = gsl_vector_get(state->x,0);
	double T = gsl_vector_get(state->x,1);
	double B = gsl_vector_get(state->x,2);

	printf("# t, A*exp(-t/T)+B\n");
	double dt = (t[n-1]-t[0])/50;
	for(double ti = t[0]; ti < t[n-1]+dt; ti+=dt)
		printf("%g %g\n",ti,MODEL(ti));

gsl_multimin_fminimizer_free(state);
gsl_vector_free(start);
gsl_vector_free(step);
return 0;
}
