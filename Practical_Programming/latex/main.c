#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gsl/gsl_odeiv2.h>
#include<gsl/gsl_errno.h>

int error_ode(double x, const double y[], double dydx[], void* params){
	dydx[0] = 2/sqrt(M_PI) * exp(-pow(x,2));
	return GSL_SUCCESS;
}

double my_error_fun(double x){
	gsl_odeiv2_system sys;
	sys.function = error_ode;
	sys.jacobian = NULL;
	sys.dimension = 1;
	sys.params = NULL;

	double hstart = copysign(0.1,x);
	double acc = 1e-6;
	double eps = 1e-6;
	gsl_odeiv2_driver* driver = gsl_odeiv2_driver_alloc_y_new
		(&sys, gsl_odeiv2_step_rkf45,hstart,acc,eps);

	double x_0 = 0;
	double y[1] = {0};
	gsl_odeiv2_driver_apply(driver,&x_0,x,y);

	gsl_odeiv2_driver_free(driver);
	return y[0];
}

int main(int argc, char* argv[]){
	int nr_of_inputs = argc;

	double start = atof(argv[1]);
	double stop = atof(argv[2]);
	double interval = atof(argv[3]);

	for(double x = start; x < stop+1e-5; x += interval)
		printf("%g %g\n", x, my_error_fun(x));

	return 0;
}
