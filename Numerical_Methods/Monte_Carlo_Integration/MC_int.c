#include"MC_int.h"

void randomx(gsl_vector* a, gsl_vector* b, gsl_vector* x){
	for(int i = 0; i < a->size; ++i){
		double a_i = gsl_vector_get(a, i);
		double b_i = gsl_vector_get(b, i);
		double RND = ((double)rand()/RAND_MAX);
		double x_i = a_i + RND*(b_i - a_i);
		gsl_vector_set(x, i, x_i);
	}
}


int mc_plain(double fun(gsl_vector*), gsl_vector* a, gsl_vector* b, int N, double* result, double* error){

	if(N <= 0){
		fprintf(stderr, "error in mc_plain: N not positive.");
		return -1;
	}

	if(a->size != b->size){
		fprintf(stderr, "error in mc_plain: vectors must be same size.");
		return -1;
	}

	double V = 1;
	for(int i = 0; i < a->size; ++i){
		double a_i = gsl_vector_get(a, i);
		double b_i = gsl_vector_get(b, i);
		V *= b_i - a_i;
	}

	gsl_vector* x = gsl_vector_alloc(a->size);
	double sum = 0, sum2 = 0;
	for(int i = 0; i < N; ++i){
		randomx(a, b, x);
		double fx = fun(x);
		sum += fx;
		sum2 += fx*fx;
	}

	double avr = sum/N;
	double var = sum2/N - avr*avr;
	*result = avr*V;
	*error = sqrt(var/N)*V;

	return 0;
}
