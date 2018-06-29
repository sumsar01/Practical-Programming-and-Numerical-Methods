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

double dot_product(gsl_vector* x, gsl_vector* y);

int min_newton(double (*f)(gsl_vector* x), void grad(gsl_vector* x, gsl_vector* df),void hessian(gsl_vector* x, gsl_matrix* H),gsl_vector* xstart, double eps);

int min_newton_num(double (*f)(gsl_vector* x), gsl_vector* xstart, double dx, double eps);

double dot_product(gsl_vector* x, gsl_vector* y);

int vector_sum(gsl_vector* x, gsl_vector* y, double a);

int set_identity(gsl_matrix* A);

int mult_matrix_vector(gsl_matrix* A, gsl_vector* x, gsl_vector* b);

int mult_matrix_matrix(gsl_matrix* A, gsl_matrix* B, gsl_matrix* C);

int matrix_copy(gsl_matrix* A, gsl_matrix* B);

int vector_copy(gsl_vector* v, gsl_vector* u);
