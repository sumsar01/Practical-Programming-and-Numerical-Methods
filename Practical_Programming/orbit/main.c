#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

int orbit_ode(double x, const double y[], double yprime[], void *params){
//	double epsilon = *(double *)params;
	yprime[0] = y[1];
	yprime[1] = 1 - y[0] + (*(double*)params)*y[0]*y[0];
return GSL_SUCCESS;
}


double orbit_equation(double x,double start_val[],double epsilon){
	gsl_odeiv2_system sys;
	sys.function = orbit_ode;
	sys.jacobian = NULL;
	sys.dimension = 2;
	sys.params = (void*) &epsilon;
	double y[2]; y[0] = start_val[0]; y[1] = start_val[1];

	double hstart = 1e-3;
	double acc = 1e-8;
	double eps = 1e-8;
	gsl_odeiv2_driver* driver =
		gsl_odeiv2_driver_alloc_y_new
				(&sys, gsl_odeiv2_step_rkf45, hstart, acc, eps);

	double x_0 = 0;
	gsl_odeiv2_driver_apply(driver, &x_0, x, y);

	gsl_odeiv2_driver_free(driver);
	return y[0];

}


int main(){

//problem 2

	FILE* file = fopen("orbit.dat","w");
	double start_val[2];
	start_val[0] = 1.0;
	start_val[1] = 0.0;
	double epsilon = 0.0;

	for(double x=0; x<= 100; x++)
		fprintf(file, "%g %g\n",2*M_PI/100*x,
			orbit_equation(2*M_PI/100*x, start_val, epsilon));



	fprintf(file,"\n\n");
	start_val[1] = -0.5;

	for(double x=0; x<= 100; x++)
		fprintf(file, "%g %g\n",2*M_PI/100*x,
			orbit_equation(2*M_PI/100*x, start_val, epsilon));

	fprintf(file,"\n\n");
	epsilon = 0.01;

	for(double x=0; x<= 1500; x++)
		fprintf(file, "%g %g\n",2*M_PI/100*x,
			orbit_equation(2*M_PI/100*x, start_val, epsilon));

return 0;
}






