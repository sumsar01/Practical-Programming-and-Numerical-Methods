#include "MC_int.h"


void print_matrix(gsl_matrix* A){
	for (int i = 0; i < A->size1; ++i) {
		for (int j = 0; j < A->size2; ++j) {
		double A_ij = gsl_matrix_get(A, i, j);
		A_ij = round(A_ij * 1e12) * 1e-12;
		printf("%.2g\t", A_ij);
		}
	printf("\n");
	}
}



void print_vector(gsl_vector* x){
	for(int i = 0; i < x->size; ++i) {
		double x_i = gsl_vector_get(x, i);
		x_i = round(x_i * 1e12) * 1e-12;
		printf("%.2g\n", x_i);
	}
}





int matrix_copy(gsl_matrix* A, gsl_matrix* B){
	if(A->size1 != B->size1 || A->size2 != B->size2){
		fprintf(stderr, "Error in matrix_copy: dimention mismatch\n");
		return -1;
	}


	for (int i = 0; i < A->size1; ++i) {
		for (int j = 0; j < A->size2; ++j) {
			gsl_matrix_set(B, i, j, gsl_matrix_get(A, i, j));
		}
	}

	return -1;
}





int vector_copy(gsl_vector* v, gsl_vector* u) {
	if (v->size != u->size) {
		fprintf(stderr, "Error in vector_copy: dimention mismatch\n");
		return -1;
	}

	for(int i = 0; i < v->size; ++i) {
		gsl_vector_set(u, i, gsl_vector_get(v, i));
	}

	return 0;
}



double dot_product(gsl_vector* x, gsl_vector* y){
	double dot_prod = 0.;

	for (int i = 0; i < x->size; ++i) {
		dot_prod += gsl_vector_get(x, i) * gsl_vector_get(y, i);
	}

	return dot_prod;
}



int vector_sum(gsl_vector* x, gsl_vector* y, double a){
	if (x->size != y->size){
		fprintf(stderr, "Error in vector_sum: vector x and y must be same length.\n");
		return -1;
	}

	for (int i = 0; i < x->size; ++i) {
		double x_i = gsl_vector_get(x, i);
		double y_i = gsl_vector_get(y, i);
		x_i = x_i + a*y_i;
		gsl_vector_set(x, i, x_i);
	}

	return 0;
}


int set_identity(gsl_matrix* A){
	if (A->size1 != A->size2){
		fprintf(stderr, "Error in set_identity: A must be square.\n");
		return -1;
	}

	for (int i = 0; i < A->size1; ++i) {
		gsl_matrix_set(A, i, i, 1.);
		for (int j = i+1; j < A->size2; ++j) {
			gsl_matrix_set(A, i, j, 0.);
			gsl_matrix_set(A, j, i, 0.);
		}
	}

return 0;
}





int mult_matrix_vector(gsl_matrix* A, gsl_vector* x, gsl_vector* b){
	if(A->size1 != A->size2 || A->size1 != x->size || x->size != b->size){
		fprintf(stderr, "Error in multiplication_matrix_vector: Dimentions mismatch.\n");
		return -1;
	}

	for (int i = 0; i < x->size; ++i) {
		double b_i = 0;
		for (int j = 0; j < x->size; ++j) {
			double A_ij = gsl_matrix_get(A, i, j);
			double x_j = gsl_vector_get(x, j);
			b_i += A_ij*x_j;
		}

		gsl_vector_set(b, i, b_i);
	}

	return 0;
}



int mult_matrix_matrix(gsl_matrix* A, gsl_matrix* B, gsl_matrix* C){
	if(A->size1 != A->size2 || B->size1 != B->size2 || C->size1 != C->size2 || A->size1 != B->size1 || B->size1 != C->size1){
		fprintf(stderr, "Error in multiplication_matrix_matrix: Dimentions mismatch.\n");
		return -1;
	}

	for (int i = 0; i < A->size1; ++i) {
		for (int j = 0; j < A->size1; ++j) {
			double C_ij = 0;
			for (int k = 0; k < A->size1; ++k) {
				double A_ik = gsl_matrix_get(A, i, k);
				double B_kj = gsl_matrix_get(A, k, j);
				C_ij += A_ik*B_kj;
			}
			gsl_matrix_set(C, i, j, C_ij);
		}
	}

	return 0;
}
