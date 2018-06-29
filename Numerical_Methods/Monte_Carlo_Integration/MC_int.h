#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<assert.h>

int mc_plain(double fun(gsl_vector* A), gsl_vector* a, gsl_vector* b, int N, double* result, double* error);

void print_matrix(gsl_matrix *A);

void print_vector(gsl_vector *x);

double dot_product(gsl_vector* x, gsl_vector* y);

double dot_product(gsl_vector* x, gsl_vector* y);

int vector_sum(gsl_vector* x, gsl_vector* y, double a);

int set_identity(gsl_matrix* A);

int mult_matrix_vector(gsl_matrix* A, gsl_vector* x, gsl_vector* b);

int mult_matrix_matrix(gsl_matrix* A, gsl_matrix* B, gsl_matrix* C);

int matrix_copy(gsl_matrix* A, gsl_matrix* B);

int vector_copy(gsl_vector* v, gsl_vector* u);


