#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <stdio.h>
#include <stdlib.h>

int main(void){

double matrix[] = 	{ 6.13, -2.90, 5.86,
			 8.08, -6.31, -3.89,
			 -4.36, 1.00, 0.19 };

double vector[] = 	{6.23, 5.37, 2.29};

gsl_matrix_view m = gsl_matrix_view_array(matrix,3,3);

gsl_vector_view b = gsl_vector_view_array(vector,3);

gsl_vector *x = gsl_vector_alloc(3);

gsl_matrix* A = gsl_matrix_alloc(3,3);

gsl_matrix_memcpy(A,&m.matrix);

int s;

gsl_permutation * p = gsl_permutation_alloc(3);

gsl_linalg_LU_decomp(&m.matrix, p ,&s);

gsl_linalg_LU_solve(&m.matrix, p ,&b.vector, x);

// test

gsl_vector* y = gsl_vector_alloc(3);
gsl_blas_dgemv(CblasNoTrans, 1, A,x,0,y);


printf("solution x to Mx=b\n");
gsl_vector_fprintf(stdout,x,"%g");

printf("result of b with solution\n");
gsl_vector_fprintf(stdout,y,"%g");

printf("b should be\n 6.23\n 5.37 \n 2.29\n");


gsl_permutation_free(p);
gsl_vector_free(x);
gsl_vector_free(y);
gsl_matrix_free(A);



return 0;
}
