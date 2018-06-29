#include <stdio.h>
#include <gsl/gsl_sf_airy.h>
#include <math.h>

int main(){
	for(double x = -2*M_PI; x < 2*M_PI; x += 0.013)
//		printf("%g\n", x);


//		int mode =1 //
//  GSL_VEGAS_MODE_IMPORTANCE_ONLY

		printf("%g %g %g \n", x,
		gsl_sf_airy_Ai(x,GSL_PREC_DOUBLE),
		gsl_sf_airy_Bi(x,GSL_PREC_DOUBLE));

return 0;

}
