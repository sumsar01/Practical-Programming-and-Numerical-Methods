#include "root_finding.h"


//performs QR decomposistion using Gram_Schmidt orthogonalization

int QR_gs_decomp(gsl_matrix *A, gsl_matrix *R){
	int n = A->size1;
	int m = A->size2;

	if(m != R->size1){
		fprintf(stderr,"error in QR_gs_decomp: wrong dim of triangular matrix R\n");
		return -1;
	}
	if(m != R->size2){
		fprintf(stderr,"error in QR_gs_decomp: wrong dim of triangular matrix R\n");
		return -1;
	}


	for(int i = 0; i<m; i++){

		double dot_prod = 0;
		for(int k = 0; k<n; ++k){
			double A_ki = gsl_matrix_get(A, k, i);
			dot_prod += A_ki*A_ki;
		}

		double R_ii = sqrt(dot_prod);
		gsl_matrix_set(R,i,i,R_ii);



		for(int k = 0; k<n; ++k){
			double Q_ki = gsl_matrix_get(A, k, i)/R_ii;
			gsl_matrix_set(A, k, i, Q_ki);
		}



		for(int j = i+1; j<m; ++j){

			dot_prod = 0;
			for(int k = 0; k<n; ++k){
				double Q_ki = gsl_matrix_get(A, k, i);
				double A_kj = gsl_matrix_get(A, k, j);
				dot_prod += Q_ki*A_kj;
			}

			double R_ij = dot_prod;
			gsl_matrix_set(R, i, j, R_ij);
			gsl_matrix_set(R, j, i, 0);


			for(int k = 0; k<n; ++k){
				double Q_ki = gsl_matrix_get(A, k, i);
				double A_kj = gsl_matrix_get(A, k, j);
				A_kj -= Q_ki*R_ij;
				gsl_matrix_set(A, k, j, A_kj);
			}

		}
	}



// returns 1 if A is singular and successful and 0 if just successful.

	if(gsl_matrix_get(R, m-1, m-1) < 1e-12) return 1;

	return 0;

}


int QR_gs_solver(const gsl_matrix *Q, const gsl_matrix *R, const gsl_vector *b, gsl_vector *x){
	int n = Q->size1;
	int m = Q->size2;

	if(m != R->size1){
		fprintf(stderr,"error in QR_gs_solver: wrong dim of right triangular matrix R\n");
		return -1;
	}
	if(m != R->size2){
		fprintf(stderr,"error in QR_gs_solver: wrong dim of right triangular matrix R\n");
		return -1;
	}
	if(m != x->size){
		fprintf(stderr,"error in QR_gs_solver: wrong dim of result vector x\n");
		return -1;
	}
	if(n != b->size){
		fprintf(stderr,"error in QR_gs_solver: wrong dim of vector b in QRx = b\n");
		return -1;
	}


	for(int i = 0; i<m; i++){
		double x_i = 0;
		for(int j=0; j<n; ++j){
			double Q_trans_ij = gsl_matrix_get(Q, j, i);
			double b_i = gsl_vector_get(b, j);
			x_i += Q_trans_ij*b_i;
		}
		gsl_vector_set(x, i, x_i);
	}

	for(int i = m-1; i >= 0; i--){
		double x_i = gsl_vector_get(x, i);
		for(int j = i+1; j<m; ++j){
			double R_ij = gsl_matrix_get(R, i, j);
			double x_j = gsl_vector_get(x, j);
			x_i -= R_ij*x_j;
		}
		x_i /=gsl_matrix_get(R,i,i);
		gsl_vector_set(x, i, x_i);
	}

	return 0;

}


int QR_gs_inverse(const gsl_matrix *R, gsl_matrix *B){
	int n = R->size1;

	if(n != B->size1){
		fprintf(stderr,"error in QR_gs_inverse: wrong dim of result matrix B\n");
		return -1;
	}
	if(B->size1 != B->size2){
		fprintf(stderr,"error in QR_gs_inverse: Result matrix must be square \n");
		return -1;
	}
	if(n != R->size1){
		fprintf(stderr,"error in QR_gs_inverse: wrong dim of right triangular matrix R\n");
		return -1;
	}
	if(R->size1 != R->size2){
		fprintf(stderr,"error in QR_gs_inverse: Upper right triangular ~/prog/numeriske_metoder/interpolationmatrix R must be square\n");
		return -1;
	}

	gsl_matrix *I = gsl_matrix_alloc(n,n);
	for(int i = 0; i<n; ++i){
		gsl_matrix_set(I, i, i, 1);
		for(int j = i+1; j<n; ++j){
			gsl_matrix_set(I, i, j, 0);
			gsl_matrix_set(I, j, i, 0);
		}
	}


	for(int i = 0; i<n; i++){
		gsl_vector_view B_i = gsl_matrix_column(B,i);
		gsl_vector_const_view e_i = gsl_matrix_const_column(I, i);
		QR_gs_solver(I, R, &(e_i.vector), &(B_i.vector));
	}

	gsl_matrix_free(I);
	return 0;
}



























