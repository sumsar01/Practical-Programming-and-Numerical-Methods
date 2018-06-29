#include"MC_int.h"

//Testing functions
double fun1(gsl_vector* x){
	double x1 = gsl_vector_get(x, 0);
	double x2 = gsl_vector_get(x, 1);
	double x3 = gsl_vector_get(x, 2);
	return 1/(1- cos(x1)*cos(x2)*cos(x3))/pow(M_PI,3);
}

double fun2(gsl_vector* x){
	return atan(exp(gsl_vector_get(x, 0)));
}

double circle(gsl_vector* x){
	return gsl_vector_get(x, 0);
}



int test_integrals(void){

	int N = 1e6;
	printf("\nf(x) = arctan(exp(x)) integrated from -2 to 2 is found to be\nI = ");
	gsl_vector* a = gsl_vector_alloc(1);
	gsl_vector_set(a, 0, -2);
	gsl_vector* b = gsl_vector_alloc(1);
	gsl_vector_set(b, 0, 2);
	double result, error;
	int error_message = mc_plain(fun2, a, b, N, &result, &error);
	if(error_message == -1) return -1;
	printf("%g +- %g\nWith Monte-carlo using N = %i points\n", result, error, N);
	gsl_vector_free(a);
	gsl_vector_free(b);


	printf("\nf(x) = 1 integrated for r = [0,1] and phi = [0, 2*Pi] found to be\nI = ");
	a = gsl_vector_alloc(2);
	gsl_vector_set(a, 0, 0.);
	gsl_vector_set(a, 1, 0.);
	b = gsl_vector_alloc(2);
	gsl_vector_set(b, 0, 1);
	gsl_vector_set(b, 1, 2*M_PI);
	error_message = mc_plain(circle, a, b, N, &result, &error);
	if(error_message == -1) return -1;
	printf("%g +- %g\nWith Monte-carlo using N = %i pooints\n", result, error, N);
	gsl_vector_free(a);
	gsl_vector_free(b);


	printf("\nf(x) = 1 / (1 - cos(x)*cos(y)*cos(z)) / Pi^3 integrated from -Pi to Pi in x, y, and z is found to be\nI = ");
	a = gsl_vector_alloc(3);
	gsl_vector_set(a, 0, 0.);
	gsl_vector_set(a, 1, 0.);
	gsl_vector_set(a, 2, 0.);
	b = gsl_vector_alloc(3);
	gsl_vector_set(b, 0, M_PI);
	gsl_vector_set(b, 1, M_PI);
	gsl_vector_set(b, 2, M_PI);
	error_message = mc_plain(fun1, a, b, N, &result, &error);
	if(error_message == -1) return -1;
	printf("%g +- %g\nWith Monte-carlo using N = %i points\n\n", result, error, N);
	gsl_vector_free(a);
	gsl_vector_free(b);


	return 0;

}



int test_error(void){

	FILE* file = fopen("error_mc.txt", "w");
	double result, error;
	gsl_vector* a = gsl_vector_alloc(1);
	gsl_vector_set(a, 0, -2);
	gsl_vector* b = gsl_vector_alloc(1);
	gsl_vector_set(b, 0, 2);


	for(double N = 1; N < 5; N += 0.05){
		int error_message = mc_plain(fun2, a, b, (int)round(pow(N,10)), &result, &error);
		if(error_message == -1) return -1;
		fprintf(file, "%i %g\n", (int)round(pow(N,10)), fabs(result-M_PI));
	}

	return 0;
}


int main(void){

	int error_message = test_integrals();
	if(error_message == -1) return -1;

	error_message = test_error();
	if(error_message == -1) return -1;
}



















