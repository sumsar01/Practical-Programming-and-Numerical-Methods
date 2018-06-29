#include"adapt_int.h"
#include<gsl/gsl_integration.h>

int fun_count;
//test functions

double fun1(double x){
	return sqrt(x);
}

double fun2(double x){
        return 1./sqrt(x);
}

double fun3(double x){
        return log(x)/sqrt(x);
}

double fun4(double x){
	fun_count++;
        return 4*sqrt(1-pow(1-x,2));
}

double fun5(double x){
        fun_count++;
	return exp(-x);
}

double fun5_gsl(double x, void* params){
	fun_count++;
        return exp(-x);
}

double fun6(double x){
	fun_count++;
        return exp(-pow(x,2))/sqrt(M_PI);
}

double fun6_gsl(double x, void* params){
	fun_count++;
        return exp(-pow(x,2))/sqrt(M_PI);
}

int main(void){

	printf("\nf(x) = sqrt(x) integrated from 0 to 1 is found to be\nI = ");
	double Q = adapt(fun1, 0., 1., 1e-6, 1e-6);
	printf("%g\n",Q);

	 printf("\nf(x) = 1/sqrt(x) integrated from 0 to 1 is found to be\nI = ");
        Q = adapt(fun2, 0., 1., 1e-6, 1e-6);
        printf("%g\n",Q);

	printf("\nf(x) = log(x)/sqrt(x) integrated from 0 to 1 is found to be\nI = ");
        Q = adapt(fun3, 0., 1., 1e-6, 1e-6);
        printf("%g\n",Q);

	fun_count = 0;
	printf("\nf(x) = 4*sqrt(1 - (1-x)^2) integrated from 0 to 1 is found to be\nI = ");
        Q = adapt(fun4, 0., 1., 1e-6, 1e-6);
        printf("%.20g\n",Q);
	printf("I-Pi = %.20g\n", Q-M_PI);
	printf("In %i function calls, with an accuracy of eps = 1e-6 and acc = 1e-6.\n\n", fun_count);

	fun_count = 0;
        printf("\nf(x) = exp(-x) integrated from 0 to inf is found to be\nI = ");
        Q = adapt(fun5, 0., INFINITY, 1e-6, 1e-6);
        printf("%g\n",Q);
        printf("In %i function calls, with an accuracy of eps = 1e-6 and acc = 1e-6.\n\n", fun_count);

	double result, abserr;
	fun_count = 0;
	gsl_function F;
	F.function = fun5_gsl;
	F.params = NULL;
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(100);
	gsl_integration_qagiu(&F, 0, 1e-6, 1e-6, 100, w, &result, &abserr);
	printf("And with gsl integrator.\nI = %g\n", result);
	printf("In %i function calls\n\n", fun_count);

	fun_count = 0;
	printf("\nf(x) = exp(-x^2)/sqrt(Pi) integrated from -inf to inf is found to be\nI = ");
	Q = adapt(fun6, -INFINITY, INFINITY, 1e-6, 1e-6);
	printf("%g\n", Q);
	printf("In %i function calls, with an accuracy of eps = 1e-6 and acc = 1e-6.\n", fun_count);

	fun_count = 0;
	F.function = fun6_gsl;
	gsl_integration_qagi(&F, 1e-6, 1e-6, 100, w, &result, &abserr);
	printf("And with gsl integrator.\nI = %g\n", result);
	printf("In %i function calls\n\n", fun_count);


	gsl_integration_workspace_free(w);

	return 0;
}










