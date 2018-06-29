#include "jacobi.h"
#include <gsl/gsl_blas.h>

int main_sweep(void){

	gsl_matrix* A = gsl_matrix_alloc(3,3);
	gsl_vector* e = gsl_vector_alloc(3);
	gsl_matrix* V = gsl_matrix_alloc(3,3);
	gsl_matrix* D = gsl_matrix_alloc(3,3);


	double rand_num = ((double) rand())/((double)RAND_MAX)*10-5;


	for(int i=0; i < A->size1; ++i){
		gsl_matrix_set(A, i, i, rand_num);
		for(int j = i+1; j < A->size2; ++j){
			double rand_num = ((double) rand())/((double)RAND_MAX)*10-5;
			gsl_matrix_set(A, i, j, rand_num);
			gsl_matrix_set(A, j, i, rand_num);
		}
	}


	printf("\nTesting Jacobi diagonalization with sweeps.\n");
	printf("Random 3x3 matrix A found to be: \nA = \n");
	print_matrix(A);


	int error = jacobi_diag_sweep(A, V, e);
	if(error == -1) return -1;

	for(int i = 0; i< D->size1; ++i){
		double ev_i = gsl_vector_get(e, i);
		gsl_matrix_set(D, i, i, ev_i);
		for(int j=i+1; j< D->size1; ++j){
			gsl_matrix_set(D, i, j, 0);
			gsl_matrix_set(D, j, i, 0);
		}
	}

	printf("\nA decomposed into V D V^T via Jacobi diagonalization using sweeps.\n V = \n");
	print_matrix(V);
	printf("D = \n");
	print_matrix(D);

	printf("\nV^T A V = \n");

	gsl_matrix* B = gsl_matrix_alloc(3,3);
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1., V, A, 0., B);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., B, V, 0., A);

	print_matrix(A);

	gsl_matrix_free(A);
	gsl_matrix_free(B);
	gsl_matrix_free(V);
	gsl_matrix_free(D);
	gsl_vector_free(e);

	return 0;
}

int main_ev_by_ev_test(void){
        gsl_matrix* A = gsl_matrix_alloc(3,3);
        gsl_vector* e = gsl_vector_alloc(3);
        gsl_matrix* V = gsl_matrix_alloc(3,3);
        gsl_matrix* D = gsl_matrix_alloc(3,3);

       double rand_num = ((double) rand())/((double)RAND_MAX)*10-5;


        for(int i=0; i < A->size1; ++i){
                gsl_matrix_set(A, i, i, rand_num);
                for(int j = i+1; j < A->size2; ++j){
                        double rand_num = ((double) rand())/((double)RAND_MAX)*10-5;
                        gsl_matrix_set(A, i, j, rand_num);
                        gsl_matrix_set(A, j, i, rand_num);
                }
        }

	printf("\nTesting Jacobi diagonalization with ev by ev.\n");
        printf("Random 3x3 matrix A found to be: \nA = \n");
        print_matrix(A);

        int error = jacobi_diag_ev_by_ev(A, V, e, 3);
        if(error == -1) return -1;

        for(int i = 0; i< D->size1; ++i){
                double ev_i = gsl_vector_get(e, i);
                gsl_matrix_set(D, i, i, ev_i);
                for(int j=i+1; j< D->size1; ++j){
                        gsl_matrix_set(D, i, j, 0);
                        gsl_matrix_set(D, j, i, 0);
                }
        }

        printf("\nA decomposed into V D V^T via Jacobi diagonalization using ev by ev.\n V = \n");
        print_matrix(V);
        printf("D = \n");
        print_matrix(D);

        printf("\nV^T A V = \n");

        gsl_matrix* B = gsl_matrix_alloc(3,3);
        gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1., V, A, 0., B);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., B, V, 0., A);

        print_matrix(A);

        gsl_matrix_free(A);
        gsl_matrix_free(B);
        gsl_matrix_free(V);
        gsl_matrix_free(D);
        gsl_vector_free(e);

        return 0;
}

int main(void){
	int error = main_sweep();
	if (error == -1) return -1;

	error = main_ev_by_ev_test();
	if(error == -1) return -1;

	return 0;
}
