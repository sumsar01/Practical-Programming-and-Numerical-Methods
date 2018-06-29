#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

int QR_gs_decomp(gsl_matrix *A, gsl_matrix *R);

int QR_gs_solver(const gsl_matrix *Q, const gsl_matrix *R, const gsl_vector *b, gsl_vector *x);

int QR_gs_inverse(const gsl_matrix *Q, const gsl_matrix *R, gsl_matrix *B);

void matrix_print(gsl_matrix *A);

void vector_print(gsl_vector *x);

int tall_matrix_test();

int square_matrix_test(void);
