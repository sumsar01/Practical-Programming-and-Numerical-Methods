#include "root_finding.h"

double vector_norm(gsl_vector* x){
	double norm = 0.;
	for(int i=0; i< x->size; ++i){
		norm += pow(gsl_vector_get(x,i),2);
	}

	return sqrt(norm);
}

int vector_sum(gsl_vector* x, gsl_vector* y, double a){
	assert(x->size = y->size);

	for(int i = 0; i< x->size; ++i){
		double x_i = gsl_vector_get(x,i);
		double y_i = gsl_vector_get(y,i);

		x_i = x_i + a*y_i;
		gsl_vector_set(x, i, x_i);
	}

	return 0;
}

//calculates jacobian of fun with step dx

int jacobian_num(void (*fun)(gsl_vector* x, gsl_vector* fx), gsl_matrix* J,
	 gsl_vector* x, gsl_vector* fx, gsl_vector* fx_dx, double dx){

	assert(J->size1 = J->size2);
	assert(J->size2 = x->size);
	assert(x->size = fx->size);
	assert(fx->size = fx_dx->size);

	for(int i = 0; i< x->size; ++i){
		double x_i = gsl_vector_get(x,i);
		gsl_vector_set(x,i,x_i+dx);
		fun(x, fx_dx);
		gsl_vector_set(x, i, x_i);

		for(int j = 0; j< x->size; ++j){
			double df_jdx_i =
			(gsl_vector_get(fx_dx,j) - gsl_vector_get(fx,j))/dx;
			gsl_matrix_set(J,j,i,df_jdx_i);
		}
	}

	return 0;
}

int newton_num(void (*fun)(gsl_vector* x, gsl_vector* fx), gsl_vector* x,
	 double dx, double eps){

	assert(dx > 0);
	assert(eps > 0);

	int n = x->size;
	int num_iter = 0;

	gsl_matrix* J = gsl_matrix_alloc(n,n);
	gsl_matrix* R = gsl_matrix_alloc(n,n);
	gsl_vector* delta_x = gsl_vector_alloc(n);
	gsl_vector* fx = gsl_vector_alloc(n);
	//fx_dx used both f_i(x_1,...,x_n) and f_i(x+delta_x)
	gsl_vector* fx_dx = gsl_vector_alloc(n);

	fun(x, fx);
	double norm_fx = vector_norm(fx);

	do{

		int error = jacobian_num(fun, J, x, fx, fx_dx, dx);
		if(error == -1) return -1;

		QR_gs_decomp(J, R);
		QR_gs_solver(J, R, fx, delta_x);
	//delta_x is -delta_x, to save computations
	//for this reason the sign of lambda is flipped
	//such that x = x -lambda*delta_x is implemented as vector_sum(x, delta_x, lambda)
		double lambda = 1.;

		vector_sum(x, delta_x, -lambda);//x = x+lambda*delta_x
		fun(x, fx_dx);
		double norm_fx_dx = vector_norm(fx_dx);

		while(norm_fx_dx > (1-lambda/2.)*norm_fx && lambda > dx/10.){
			lambda /= 2.;
			vector_sum(x,delta_x,lambda);//x = x - lambda*delta_x
			fun(x,fx_dx);
			norm_fx_dx = vector_norm(fx_dx);
		}

		fun(x,fx);
		norm_fx = vector_norm(fx);

		num_iter++;

	}while(norm_fx > eps && num_iter < 1e6);

	if(num_iter >= 1e6){
		fprintf(stderr,"Error newton_jacobian did not converge after  %i iterations.\n", num_iter);
		return -1;
	}

	gsl_matrix_free(J);
	gsl_matrix_free(R);
	gsl_vector_free(delta_x);
	gsl_vector_free(fx);
	gsl_vector_free(fx_dx);

	return num_iter;
}


int newton_jacobian(void (*fun)(gsl_vector* x, gsl_vector* fx),
	void (*jacobian)(gsl_vector* x, gsl_matrix* J), gsl_vector* x, double eps){

	int n = x->size;
	int num_iter = 0;

        gsl_matrix* J = gsl_matrix_alloc(n,n);
        gsl_matrix* R = gsl_matrix_alloc(n,n);
        gsl_vector* delta_x = gsl_vector_alloc(n);
        gsl_vector* fx = gsl_vector_alloc(n);
        //fx_dx used both f_i(x_1,...,x_n) and f_i(x+delta_x)
        gsl_vector* fx_dx = gsl_vector_alloc(n);

        fun(x, fx);
        double norm_fx = vector_norm(fx);

	do{

		jacobian(x,J);

		QR_gs_decomp(J, R);
                QR_gs_solver(J, R, fx, delta_x);
        //delta_x is -delta_x, to save computations
        //for this reason the sign of lambda is flipped
        //such that x = x -lambda*delta_x is implemented as vector_sum(x, delta_x,$
                double lambda = 1.;

                vector_sum(x, delta_x, -lambda);//x = x+lambda*delta_x
                fun(x, fx_dx);
                double norm_fx_dx = vector_norm(fx_dx);

                while(norm_fx_dx > (1-lambda/2.)*norm_fx && lambda > 0.01){
                        lambda /= 2.;
                        vector_sum(x,delta_x,lambda);//x = x - lambda*delta_x
                        fun(x,fx_dx);
                        norm_fx_dx = vector_norm(fx_dx);
                }

                fun(x,fx);
                norm_fx = vector_norm(fx);

                num_iter++;

        }while(norm_fx > eps && num_iter < 1e6);

        if(num_iter >= 1e6){
                fprintf(stderr,"Error newton_jacobian did not converge after  %i iterations.\n",num_iter);
                return -1;
        }

        gsl_matrix_free(J);
        gsl_matrix_free(R);
        gsl_vector_free(delta_x);
        gsl_vector_free(fx);
        gsl_vector_free(fx_dx);

        return num_iter;
}




















