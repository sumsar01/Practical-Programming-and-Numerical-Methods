#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multiroots.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#define STEPPER gsl_odeiv2_step_rkf45
#define ABSERR 1e-6
#define RELERR 1e-6
#define STARTSTEP 1e-3

int ode_wave(double r, const double y[], double dy[], void* params){
	double e = *(double*)params;
	dy[0] = y[1];
	dy[1] = 2*(-1/r-e) * y[0];
return GSL_SUCCESS;
}

double Fe(double e, double r){
	assert(r >= 0);
	const double rmin = 1e-3;
	if(r<rmin) return r-r*r;

	gsl_odeiv2_system system;
	system.function = ode_wave;
	system.jacobian = NULL;
	system.dimension = 2;
	system.params = (void*)&e;

	gsl_odeiv2_driver* driver =
		gsl_odeiv2_driver_alloc_y_new
			(&system, STEPPER, STARTSTEP, ABSERR, RELERR);

	double t = rmin, y[] = {t-t*t, 1-2*t};
	int status = gsl_odeiv2_driver_apply( driver, &t, r, y);
	if(status != GSL_SUCCESS) fprintf(stderr,"Fe: odeiv2 error: %d\n",status);

	gsl_odeiv2_driver_free(driver);

return y[0];
}


