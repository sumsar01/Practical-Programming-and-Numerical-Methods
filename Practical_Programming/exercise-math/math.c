#include "tgmath.h"
#include "complex.h"
#include "stdio.h"
#include "complex.h"
#include "math.h"

int main(){
	printf("\nProblem 1\n");

	double x=5;
	double g=tgamma(x);
	printf("gamma(5) = %g\n",g);

	double y=0.5;
	printf("j1(0.5) = %g\n",j1(y));

	double complex  z=csqrt(-2.0);
	printf("sqrt(-2) = %fi\n",cimag(z));

	double complex w=I;
	double complex c=cexp(w);
	printf("e^i= %f + i %f\n",creal(c),cimag(c));

	double complex c_pi=M_PI*I;
	double complex e_pi=cexp(c_pi);
	printf("e^i*pi = %g + i %g\n",creal(e_pi),cimag(e_pi));

	double complex p=cpow(I,M_E);
	printf("i^e = %g + i %g\n",creal(p),cimag(p));

	float b=1.0/9;
	double n=1.0/9;
	long double m=1.0/9;
	printf("\nProblem 2\n");
	printf("b = %.25g \nn = %.25lg \nm = %.25Lg\n",b,n,m);
}
