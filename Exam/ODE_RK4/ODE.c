#include"ODE.h"


int rkstep4(double t, double h, gsl_vector* y, void f(double t, gsl_vector* y, gsl_vector* dydt), gsl_vector* yh, gsl_vector* err, gsl_vector* y2){

	int error = vector_copy(y, y2);
	if(error == -1) return -1;

//doing full-step

	f(t, y, err); // err being k_0 here
	error = vector_sum(y, 1., err, 1./6*h);//y = y + 1/6*k_0
	if(error == -1) return -1;


	error = vector_sum(err, 1/2.*h, y, 1.);
	if(error == -1) return -1;
	f(t + h/2., err, yh);//err being k_1 here
	error = vector_sum(y, 1., yh, 1./3*h);//y = y + 1/6*k_0 + 1/3*k_1
	if(error == -1) return -1;

	error = vector_sum(yh, 1/2.*h, y, 1.);
	if(error == -1) return -1;
	f(t + h/2., yh, err); // err being k_2 here
	error = vector_sum(y, 1., err, 1./3*h);//y = y + 1/6*k_0 + 1/3*k_1 + 1/3*k_2
	if(error == -1) return -1;


	error = vector_sum(err, h, y, 1.);
	if(error == -1) return -1;
	f(t + h, err, yh); // yh being k_3 here
	error = vector_sum(y, 1., yh, 1./6*h);//y = y + 1/6*k_0 + 1/3*k_1 + 1/3*k_2 + 1/6*k_3
	if(error == -1) return -1;

//doing two half-step

	f(t, y2, err); // err being k_0 here
	error = vector_sum(y2, 1., err, 1./6*h);//y = y + 1/6*k_0
	if(error == -1) return -1;

	error = vector_sum(err, 1/4.*h, y2, 1.);
	if(error == -1) return -1;
	f(t + h/4., err, yh);//err being k_1 here
	error = vector_sum(y2, 1., yh, 1./3*h);//y = y + 1/6*k_0 + 1/3*k_1
	if(error == -1) return -1;

	error = vector_sum(yh, 1/4.*h, y2, 1.);
	if(error == -1) return -1;
	f(t + h/4., yh, err); // err being k_2 here
	error = vector_sum(y2, 1., err, 1./3*h);//y = y + 1/6*k_0 + 1/3*k_1 + 1/3*k_2
	if(error == -1) return -1;

	error = vector_sum(err, 1./2*h, y2, 1.);
	if(error == -1) return -1;
	f(t + h/2, err, yh); // yh being k_3 here
	error = vector_sum(y2, 1., yh, 1./6*h);//y = y + 1/6*k_0 + 1/3*k_1 + 1/3*k_2 + 1/6*k_3
	if(error == -1) return -1;

	//estimating error with by doing delta_y = (y_full_step - y_two_half_steps)/(2^p-1) where the order p = 4.


	error = vector_sum(y2, -1./15, y, 1./15);
	if(error == -1) return -1;

	error = vector_copy(y2, err);
	if(error == -1) return -1;


	return 0;
}


int driver(double* t, double b, double* h, gsl_vector* y, double acc, double eps, int stepper(double t, double h, gsl_vector* y, void f(double t, gsl_vector* y, gsl_vector* dydt), gsl_vector* yh, gsl_vector* err, gsl_vector* y2), void f(double t, gsl_vector* y, gsl_vector* dydt)){

	if(*t == b) return 0;
	if(fabs(*t + *h) > fabs(b)){*h = b - *t;}

	gsl_vector* yh = gsl_vector_alloc(y->size);
	gsl_vector* err = gsl_vector_alloc(y->size);
	gsl_vector* y2 = gsl_vector_alloc(y->size);
	int iterations = 0;
	double err_norm, y_norm, tau, a = *t;

	do{

		iterations++;
		int error = stepper(*t, *h, y, f, yh, err, y2);
		if(error == -1) return -1;

		err_norm = sqrt(dot_product(err, err));
		y_norm = sqrt(dot_product(y, y));


		tau = (eps*y_norm + acc)*sqrt(*h/(b-a));

		if(tau > err_norm){
			*t = *t + *h;
		}

	//	*h *= pow(tau/err_norm, 0.25)* 0.95;
		if(fabs(*t + *h) > fabs(b)){*h = b - *t;}

	}while(iterations < 1e6 && fabs(*t - b) > 1e-6);

	if(iterations == 1e6){
		fprintf(stderr, "error in driver: did not converge after %i cycles\n", iterations);
		return -1;
	}

	return 0;
}








