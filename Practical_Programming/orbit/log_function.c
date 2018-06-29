#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

int log_ode(double x, const double y[], double yprime[], void *params){
        yprime[0] = y[0]*(1-y[0]);
return GSL_SUCCESS;
}

double mylog(double x){
        gsl_odeiv2_system sys;
        sys.function = log_ode;
        sys.jacobian = NULL;
        sys.dimension = 1;
        sys.params = NULL;

        double hstart = 0.1;
        double acc = 1e-6;
        double eps = 1e-6;
        gsl_odeiv2_driver* driver =
                gsl_odeiv2_driver_alloc_y_new
                        (&sys, gsl_odeiv2_step_rkf45, hstart, acc, eps);

        double t = 0;
        double y[2] = {0.5};
        gsl_odeiv2_driver_apply(driver, &t,x,y);

        gsl_odeiv2_driver_free(driver);
        return y[0];
}

int main(){
        for(double x=0;x<3;x+=0.1)
                printf("%g %g\n",x,mylog(x));
        printf("\n");
        printf("\n");

        for(double x=0; x<3;x+=0.3)
                printf("%g %g\n",x,1/(1+exp(-x)));
}
