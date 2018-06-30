#include "ODE.h"

//test functions
void trig_via_ODE(double t, gsl_vector* y, gsl_vector* dydt){
	double y_1 = gsl_vector_get(y, 0);
	double y_2 = gsl_vector_get(y, 1);
	double dydt_1 = y_2;
	double dydt_2 = -y_1;
	gsl_vector_set(dydt, 0, dydt_1);
	gsl_vector_set(dydt, 1, dydt_2);
}


void airy_fun(double t, gsl_vector* y, gsl_vector* dydt){
	double y_1 = gsl_vector_get(y, 0);
	double y_2 = gsl_vector_get(y, 1);
	double dydt_1 = y_2;
	double dydt_2 = t*y_1;
	gsl_vector_set(dydt, 0, dydt_1);
	gsl_vector_set(dydt, 1, dydt_2);
}



int main(void){

	printf("\nSolving d2ydt2 = t*y to get airy function with y(0) = 0.355028 and dydt(0) = -0.258819\n");
	printf("to make sure that y->0 for t->inf\n");
	printf("by solving ODE for every point needed.\n\n");
	FILE* file = fopen("airy.txt", "w");

	double t, h, b;
	gsl_vector* y = gsl_vector_alloc(2);

	for(double x = -15; x < 6+1e-6; x+=0.1){
		gsl_vector_set(y, 0, 0.355028);
		gsl_vector_set(y, 1, -0.258819);
		h = 0.01*x/fabs(x);
		t = 0;

		int error = driver(&t, x, &h, y, 1., 1e-6, &rkstep4, &airy_fun);
		if(error == -1) return -1;

		fprintf(file, "%g %g %g\n", t, gsl_vector_get(y, 0), gsl_vector_get(y, 1));
	}


	file = fopen("sin.txt", "w");


	for(double x = 0; x < 2*M_PI; x+=0.1){

	t = 0.;
	b = 2*M_PI;
	h = 0.01*x/fabs(x);
	gsl_vector_set(y, 0, 0.);
	gsl_vector_set(y, 1, 1.);

	int error = driver(&t, x, &h, y, 1., 1e-6, &rkstep4, &trig_via_ODE);
	if(error == -1) return -1;

		fprintf(file, "%g %g %g\n", t, gsl_vector_get(y, 0), gsl_vector_get(y, 1));


	}

	printf("Integrating d2ydt2 = -y for y(0) = 0 and dydt(0) = 1 showing every point calculated on the way\n");
	printf("with an accuracy of 1e-2\n\n");

	return 0;
}








