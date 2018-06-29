#include "adapt_int.h"

double transform(double fun(double), double x, int type, double limit){
	if(type == 0) return fun(x);
	if(type == 1) return fun(limit - (1-x)/x)/pow(x,2);
	if(type == 2) return fun(limit + (1-x)/x)/pow(x,2);
	if(type == 3) return (fun((1-x)/x) + fun(-(1-x)/x))/pow(x,2);
	else{
		fprintf(stderr,"error in tranform.");
	}assert(0);
}


double adapt24(double fun(double), double a, double b, double acc, double eps, double f_2,
	 double f_3, double iterations, int type, double limit){

	if(iterations >= 1e6){
		fprintf(stderr,"error in adapt24.");
	}assert(iterations < 1e6);

	double f_1 = transform(fun, a + (b-a)/6, type, limit);
	double f_4 = transform(fun, a + 5*(b-a)/6, type, limit);
	double Q = (2*f_1 + f_2 + f_3 + 2*f_4)/6*(b-a);
	double q = (f_1 + f_2 + f_3 + f_4)/4*(b-a);
	double tol = acc + eps*fabs(Q);
	double error = fabs(Q-q);
	if(error < tol) return Q;
	else{
		double Q1 = adapt24(fun, a, (a+b)/2, acc/sqrt(2), eps, f_1, f_2, iterations+1, type, limit);
		double Q2 = adapt24(fun, (a+b)/2, b, acc/sqrt(2), eps, f_3, f_4, iterations+1, type, limit);
		return Q1 + Q2;
	}
}


double adapt(double fun(double), double a, double b, double acc, double eps){
	int iterations = 0, type;
	double limit = 0;

	if(a > b) return -adapt(fun, b, a, acc, eps);
	if(isinf(a) == 0 && isinf(b) == 0) type = 0;
	if(isinf(-a) == 1 && isinf(b) == 0){type = 1; limit = b; a = 0; b = 1;}
	if(isinf(a) == 0 && isinf(b) == 1){type = 2; limit = a; a = 0; b = 1;}
	if(isinf(-a) == 1 && isinf(b) == 1){type = 3; a = 0; b = 1;}

	double f_2 = transform(fun, a + 2*(b-a)/6, type, limit);
	double f_3 = transform(fun, a + 4*(b-a)/6, type, limit);

	return adapt24(fun, a, b, acc, eps, f_2, f_3, iterations, type, limit);
}



























