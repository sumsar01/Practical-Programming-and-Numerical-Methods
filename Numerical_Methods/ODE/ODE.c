#include"ODE.h"


int rkstep12(double t, double h, gsl_vector* y, void f(double t, gsl_vector* y, gsl_vector* dydt), gsl_vector* yh, gsl_vector* err){

	f(t, y, err); // err being k_0 here
	int error = vector_sum(y, 1., err, h/2);//y = 1.*y+h/2*err
	if(error == -1) return -1;

	f(t + h/2, y, yh);//yh being k_1/2 here
	error = vector_sum(y, 1., err, -h/2);
	if(error == -1) return -1;

	error = vector_sum(err, h/2, yh, -h/2);
	if(error == -1) return -1;

	error = vector_sum(yh, h, y, 1.);
	if(error == -1) return -1;

	return 0;
}


int driver(double* t, double b, double* h, gsl_vector* y, double acc, double eps, int stepper(double t, double h, gsl_vector* y, void f(double t, gsl_vector* y, gsl_vector* dydt), gsl_vector* yh, gsl_vector* err), void f(double t, gsl_vector* y, gsl_vector* dydt)){

	if(*t == b) return 0;
	if(fabs(*t + *h) > fabs(b)){*h = b - *t;}

	gsl_vector* yh = gsl_vector_alloc(y->size);
	gsl_vector* err = gsl_vector_alloc(y->size);
	int iterations = 0;
	double err_norm, y_norm, tau, a = *t;

	do{

		iterations++;
		int error = stepper(*t, *h, y, f, yh, err);
		if(error == -1) return -1;

		err_norm = sqrt(dot_product(err, err));
		y_norm = sqrt(dot_product(y, y));


		tau = (eps*y_norm + acc)*sqrt(*h/(b-a));
		if(tau > err_norm){
			vector_copy(yh, y);
			*t = *t + *h;
		}

		*h *= pow(tau/err_norm, 0.25)* 0.95;
		if(fabs(*t + *h) > fabs(b)){*h = b - *t;}


	}while(iterations < 1e6 && fabs(*t - b) > 1e-12);

	if(iterations == 1e6){
		fprintf(stderr, "error in driver: did not converge after %i cycles\n", iterations);
		return -1;
	}

	return 0;
}



int driver_with_path(double* t, double b, double* h, gsl_vector* y, double acc, double eps, int stepper(double t, double h, gsl_vector* y, void f(double t, gsl_vector* y, gsl_vector* dydt), gsl_vector* yh, gsl_vector* err), void f(double t, gsl_vector* y, gsl_vector* dydt), gsl_matrix* Y, int* num_saved){

	if(y->size != (Y->size2)-1){
		fprintf(stderr, "error in driver_with_path: wrong dimensions.");
		return -1;
	}


	gsl_vector* yh = gsl_vector_alloc(y->size);
	gsl_vector* err = gsl_vector_alloc(y->size);
	int iterations = 0;
	*num_saved = 0;
	double err_norm, y_norm, tau, a = *t;

	do{
		iterations++;
		int error = stepper(*t, *h, y, f, yh, err);
		if(error == -1) return -1;

		err_norm = sqrt(dot_product(err, err));
		y_norm = sqrt(dot_product(y, y));
		tau = (eps*y_norm + acc)*sqrt(*h/(b-a));

		if(tau > err_norm){
			vector_copy(yh, y);
			*t = *t + *h;

			if(*num_saved < Y->size1){
				gsl_matrix_set(Y, *num_saved, 0, *t);
				for(int i = 0; i < (y->size); ++i){
					gsl_matrix_set(Y, *num_saved, i+1, gsl_vector_get(y, i));
				}
			}
			(*num_saved)++;
	 	}

		*h *= pow(tau/err_norm, 0.25)* 0.95;

		if(*t + *h > b){*h = b - *t;}

	}while(iterations < 1e6 && fabs(*t -b) > 1e-12);


	if(iterations == 1e6){
		fprintf(stderr, "error in driver: ODE did not converge after %i cycles", iterations);
		return -1;
	}

	return 0;
}























