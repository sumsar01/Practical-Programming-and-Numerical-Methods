#include "lin_eq.h"

void matrix_print(gsl_matrix *A){
	for(int i = 0; i < A->size1; ++i){
		for(int j = 0; j < A->size2; ++j){
			double A_ij = gsl_matrix_get(A, i, j);
			A_ij = round(A_ij * 1e12)*1e-12;
			printf("%.2g\t",A_ij);
		}
		printf("\n");
	}
}

void vector_print(gsl_vector *x){
	for(int i = 0; i < x->size; ++i){
		double x_i = gsl_vector_get(x, i);
		x_i = round(x_i*1e12)*1e-12;
		printf("%.2g\n", x_i);
	}
}

int tall_matrix_test(){
	gsl_matrix* A = gsl_matrix_alloc(4,3);
	gsl_matrix* R = gsl_matrix_alloc(3,3);

	for(int i = 0; i < A->size1; ++i){
		for(int j = 0; j < A->size2; ++j){
			gsl_matrix_set(A, i, j, ((double) rand())/((double)RAND_MAX)*10-5);
		}
	}

	printf("\nTall matrix test \n");
	printf("\nRandom matrix A found to be \n A = \n");
	matrix_print(A);

	int status = QR_gs_decomp(A, R);
	if(status == -1) return -1;

	printf("\n A is QR facturiced into \n Q = \n");
	matrix_print(A);
	printf("R = \n");
	matrix_print(R);

	gsl_matrix * prod = gsl_matrix_alloc(3,3);
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1., A, A, 0, prod);

	printf("\nQ^T Q = \n");
	matrix_print(prod);
	gsl_matrix_free(prod);

	printf("\nQR = \n");
	prod = gsl_matrix_alloc(4,3);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., A, R, 0, prod);
	matrix_print(prod);

	gsl_matrix_free(prod);
	gsl_matrix_free(A);
	gsl_matrix_free(R);

	return 0;
}

int square_matrix_test(void){

	gsl_matrix *A = gsl_matrix_alloc(3,3);
	gsl_matrix *R = gsl_matrix_alloc(3,3);
	gsl_vector *b = gsl_vector_alloc(3);
	gsl_vector *x = gsl_vector_alloc(3);

	for(int i = 0; i < b->size; ++i){
		gsl_vector_set(b, i, ((double)rand())/((double)RAND_MAX)*10-5);
	}
	for(int i = 0; i < A->size1; ++i){
		for( int j = 0; j < A->size2; ++j){
			gsl_matrix_set(A, i, j, ((double) rand())/((double)RAND_MAX)*10-5);
		}
	}

	printf("\nsquare matrix test \n");

	printf("\nA =\n");
	matrix_print(A);
	printf("b = \n");
	vector_print(b);

	int status = QR_gs_decomp(A, R);
	if(status == -1) return -1;
	if(status == 1) return -1;

	status = QR_gs_solver(A, R, b, x);
	if(status == -1) return -1;

	printf("\nAx = b solution is\nx=\n");
	vector_print(x);

	gsl_matrix* prod = gsl_matrix_alloc(3,3);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., A, R, 0, prod);
	gsl_blas_dgemv(CblasNoTrans, 1., prod, x, 0., b);
	printf("\nAx found to be \nAx =\n");
	vector_print(b);

	gsl_matrix* B = gsl_matrix_alloc(3, 3);
	status = QR_gs_inverse(A, R, B);
	if(status == -1) return -1;



	printf("\n Inverse of A found to be\n A⁻¹ = \n");
	matrix_print(B);

	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., A, B, 0, prod);
	printf("\n A A⁻¹ found to be \n A A⁻¹ = \n");
	matrix_print(prod);

	gsl_matrix_free(A);
	gsl_matrix_free(R);
	gsl_matrix_free(B);
	gsl_matrix_free(prod);
	gsl_vector_free(x);
	gsl_vector_free(b);

	return 0;
}




















