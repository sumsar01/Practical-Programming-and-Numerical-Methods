#include "fit.h"

int lin_ls_sq_fit_ini(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b,
	 const fit_function* fun, const gsl_vector* x, const gsl_vector* y,
		 const gsl_vector* dy){
// setting b_i = y_i/dy_i and A_ik = f_k(x_i)/dy_i and then decomp into A = QR

	int nr_data = x->size;
	int nr_fun  = fun->nr_fun;

	for(int i = 0; i < nr_data; ++i){
		double y_i = gsl_vector_get(y, i);
		double dy_i = gsl_vector_get(dy, i);
		double b_i = y_i/dy_i;
		gsl_vector_set(b, i, b_i);
	}

	for(int k = 0; k< nr_fun; ++k){
		for(int i = 0; i < nr_data; ++i){
			double x_i = gsl_vector_get(x,i);
			double dy_i = gsl_vector_get(dy,i);
			double f_ki = (fun->function)(x_i, k);
			f_ki = f_ki/dy_i;
			gsl_matrix_set(Q, i, k, f_ki);
		}
	}

	int error = QR_gs_decomp(Q,R);
	if(error == -1) return -1;

	return 0;
}


int lin_ls_sq_fit(const fit_function* fun, const gsl_vector* x,
	const gsl_vector* y, const gsl_vector* dy, gsl_vector* c, gsl_matrix* S){

	if(x->size != y->size || dy->size != y->size){
		fprintf(stderr,"Error in lin_ls_sq_fit: x, y, and dy must be same size\n");
		return -1;
	}

	if(c->size != fun->nr_fun){
        	fprintf(stderr,"Error in lin_ls_sq_fit: c must be same length as fun in (*fit_function)\n");
                return -1;
        }


	if(S->size1 != fun->nr_fun || S->size1 != S->size2){
		fprintf(stderr,"Error in lin_ls_sq_fit: c must be same length as fun in (*fit_function)\n");
                return -1;
        }


        int nr_data = x->size;
        int nr_fun  = fun->nr_fun;

	gsl_matrix* Q = gsl_matrix_alloc(nr_data, nr_fun);
	gsl_matrix* R = gsl_matrix_alloc(nr_fun, nr_fun);
	gsl_vector* b = gsl_vector_alloc(nr_data);


	int error = lin_ls_sq_fit_ini(Q, S, b, fun, x, y, dy);
	if(error == -1) return -1;

	error = QR_gs_solver(Q, S, b, c);
	if(error == -1) return -1;


	//finding covarience S

	error = QR_gs_inverse(S, R);
	if(error == -1) return -1;


	for(int i = 0; i< S->size1; ++i){
		for(int j = 0; j < S->size2; ++j){
			double S_ij = 0;
			for(int k = 0; k < S->size1; ++k){
				double R_ik = gsl_matrix_get(R, i, k);
				double R_trans_kj = gsl_matrix_get(R, j, k);
				S_ij += R_ik*R_trans_kj;
			}
			gsl_matrix_set(S, i, j, S_ij);
		}
	}

	gsl_matrix_free(Q);
	gsl_matrix_free(R);
	gsl_vector_free(b);

	return 0;
}














