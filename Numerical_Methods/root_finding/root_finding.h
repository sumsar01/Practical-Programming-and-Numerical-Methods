#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

int QR_gs_decomp(gsl_matrix *A, gsl_matrix *R);

int QR_gs_solver(const gsl_matrix *Q, const gsl_matrix *R, const gsl_vector *b, gsl_vector *x);

int QR_gs_inverse(const gsl_matrix *R, gsl_matrix *B);

void print_matrix(gsl_matrix *A);

void print_vector(gsl_vector *x);

int newton_num(void fun(gsl_vector* x, gsl_vector* fx), gsl_vector* x, double dx, double eps);

int newton_jacobian(void (*fun)(gsl_vector* x, gsl_vector* fx), void (*jacobian)(gsl_vector* x, gsl_matrix* J), gsl_vector* x, double eps);
