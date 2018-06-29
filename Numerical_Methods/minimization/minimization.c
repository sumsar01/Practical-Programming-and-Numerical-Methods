#include"minimization.h"

int min_newton(double (*f)(gsl_vector* x), void gradient(gsl_vector* x, gsl_vector* df),void hessian(gsl_vector* x, gsl_matrix* H), gsl_vector* xstart, double eps){

	if(eps <= 0){
		fprintf(stderr, "error in min_newton\n");
		return -1;
	}

	int n = xstart->size;
	int iterations = 0;

	gsl_matrix* H = gsl_matrix_alloc(n,n);
	gsl_matrix* R = gsl_matrix_alloc(n,n);
	gsl_vector* grad_f = gsl_vector_alloc(n);
	gsl_vector* Delta_x = gsl_vector_alloc(n);


	gradient(xstart, grad_f);
	double norm_grad_f = sqrt(dot_product(grad_f, grad_f));


	do{
		hessian(xstart, H);


		QR_gs_decomp(H,R);
		QR_gs_solver(H, R, grad_f, Delta_x);
		//Delta_x is actually -Delta_x, for faster computation
		//sing of lambda is flipped, to correct for this.
		//thus x = x - lambda*Delta_x is implemented as vector_sum(x, Delta_x, lambda)
		double lambda = 1.;

		double fx = f(xstart);
		vector_sum(xstart, Delta_x, -lambda); //x = x + lambda*Delta_x


		while(f(xstart) > (fx + 1e-2*lambda*dot_product(Delta_x, grad_f))){
			lambda /= 2.;
			vector_sum(xstart, Delta_x, lambda); //x = x + lambda*Delta_x
		}

		gradient(xstart, grad_f);
		norm_grad_f = sqrt(dot_product(grad_f, grad_f));
		iterations++;

	}while(norm_grad_f > eps && iterations < 1e6);

	if(iterations >= 1e6){
		fprintf(stderr, "error in min_newton.\n");
		return -1;
	}

	gsl_matrix_free(H);
	gsl_matrix_free(R);
	gsl_vector_free(Delta_x);
	gsl_vector_free(grad_f);

	return iterations;

}


int gradient_num(double (*f)(gsl_vector* x), gsl_vector* grad_f, gsl_vector* x, double dx){

	if(grad_f->size != x->size){
		fprintf(stderr,"error in gradient_num: wrong dimensions");
		return -1;
	}


	double fx = f(x);
	for(int i = 0; i < grad_f->size; ++i){
		double x_i = gsl_vector_get(x,i);
		gsl_vector_set(x, i, x_i + dx);
		double fx_dx = f(x);
		gsl_vector_set(x, i, x_i);
		double grad_f_i = (fx_dx - fx) / dx;
		gsl_vector_set(grad_f, i, grad_f_i);
	}

	return 0;
}



int hessian_update(gsl_matrix* H, gsl_matrix* H_new, gsl_vector* y, gsl_vector* s, double lambda){

	assert(H->size1 = H->size2);
	assert(H_new->size1 = H_new->size2);
	assert(H->size1 = y->size);
	assert(y->size = s->size);


	int error = vector_sum(s, s, -lambda-1); //s = -lambda*s = -lambda*Delta_x
	if(error == -1) return -1;

	double y_trans_H_inv_s = 0.; //y_trans*H_inv*s
	for(int i = 0; i < y->size; ++i){
		double y_i = gsl_vector_get(y,i);
		for(int j = 0; j < y->size; ++j){
			double s_j = gsl_vector_get(y,j);
			double H_ij = gsl_matrix_get(H, i, j);
			y_trans_H_inv_s += y_i*s_j*H_ij;
		}
	}

	if(fabs(y_trans_H_inv_s) < 1e-3) {set_identity(H); return 0;}

	//result_ij = s_i * sum_k s_k*H_kj - sum_q H_iq*y_q sum_k s_k*H_kj
	for(int i = 0; i < H->size1; ++i){
		for(int j = 0; j < H->size1; ++j){
			double s_i = gsl_vector_get(s, i);
			double result = 0.;
			for(int k = 0; k < H->size1; ++k){
				double H_y_i = 0.;
				for(int q = 0; q < H->size1; ++q){
					double H_iq = gsl_matrix_get(H, i, q);
					double y_q = gsl_vector_get(y, q);
					H_y_i += H_iq*y_q;
				}
				double s_k = gsl_vector_get(s, k);
				double H_kj = gsl_matrix_get(H, k, j);
				result += s_i * s_k * H_kj - H_y_i * s_k * H_kj;
			}
			gsl_matrix_set(H_new, i, j, result/y_trans_H_inv_s);
		}
	}

	matrix_copy(H_new, H);

	return 0;

}


int min_newton_num(double (*f)(gsl_vector* x), gsl_vector* xstart, double dx, double eps){

	if(eps <= 0){
                fprintf(stderr, "error in min_newton\n");
                return -1;
        }


	int n = xstart->size;
	int iterations = 0;

	gsl_matrix* H = gsl_matrix_alloc(n,n); //H is H^-1
	gsl_matrix* extra_mem = gsl_matrix_alloc(n,n);
	gsl_vector* grad_f = gsl_vector_alloc(n);
	gsl_vector* grad_f_dx = gsl_vector_alloc(n);
	gsl_vector* Delta_x = gsl_vector_alloc(n);

	int error = set_identity(H);
	if(error == -1) return -1;

	error = gradient_num(f, grad_f, xstart, dx);
	if(error == -1) return -1;


	do{

		error = mult_matrix_vector(H, grad_f, Delta_x); //finds H^-1*grad_f = Delta_x
		if(error == -1) return -1;

                //Delta_x is actually -Delta_x, for faster computation
                //sing of lambda is flipped, to correct for this.
                //thus x = x - lambda*Delta_x is implemented as vector_sum(x, Delta_x, lambda)
		double lambda = 1.;

		double fx = f(xstart);
		error = vector_sum(xstart, Delta_x, -lambda);  //x = x + lambda*Delta_x
		if(error == -1) return -1;


		while(f(xstart) > (fx + 1e-2*lambda*dot_product(Delta_x, grad_f)) &&
				 lambda > dx/10){
			lambda /= 2.;
			error = vector_sum(xstart, Delta_x, lambda); // x = x - lambda*Delta_x
			if(error == -1) return -1;
		}

		error = gradient_num(f, grad_f_dx, xstart, dx);
		if(error == -1) return -1;

		double norm_grad_f_dx = sqrt(dot_product(grad_f_dx, grad_f_dx));
		iterations++;


		if(norm_grad_f_dx < eps){break;}
		else{
			error = vector_sum(grad_f, grad_f_dx, -1.);//grad_f = grad_f - grad_f_dx
			if(error == -1) return -1;
			error = vector_sum(grad_f, grad_f, -2.); //grad_f_dx = -grad_f_dx
			if(error == -1) return -1;
			error = hessian_update(H, extra_mem, grad_f, Delta_x, lambda);
			if(error == -1) return -1;
			error = vector_copy(grad_f_dx, grad_f);
			if(error == -1) return -1;

		}


	}while(iterations < 1e6);

	if(iterations >= 1e6){
		return -1;
	}


	gsl_matrix_free(H);
	gsl_matrix_free(extra_mem);
	gsl_vector_free(Delta_x);
	gsl_vector_free(grad_f);
	gsl_vector_free(grad_f_dx);

	return iterations;
}


























