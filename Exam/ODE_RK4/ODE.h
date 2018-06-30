#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<assert.h>

int rkstep4(double t, double h, gsl_vector* y, void f(double t, gsl_vector* y, gsl_vector* dydt), gsl_vector* yh, gsl_vector* err, gsl_vector* y2);

int driver(double* t, double b, double* h, gsl_vector* y, double acc, double eps, int stepper(double t, double h, gsl_vector* y, void f(double t, gsl_vector* y, gsl_vector* dydt), gsl_vector* yh, gsl_vector* err, gsl_vector* y2), void f(double t, gsl_vector* y, gsl_vector* dydt));

int driver_with_path(double* t, double b, double* h, gsl_vector* y, double acc, double eps, int stepper(double t, double h, gsl_vector* y, void f(double t, gsl_vector* y, gsl_vector* dydt), gsl_vector* yh, gsl_vector* err, gsl_vector* y2), void f(double t, gsl_vector* y, gsl_vector* dydt), gsl_matrix* Y, int* num_saved);

int QR_gs_decomp(gsl_matrix *A, gsl_matrix *R);

int QR_gs_solver(const gsl_matrix *Q, const gsl_matrix *R, const gsl_vector *b, gsl_vector *x);

int QR_gs_inverse(const gsl_matrix *R, gsl_matrix *B);

void print_matrix(gsl_matrix *A);

void print_vector(gsl_vector *x);

double dot_product(gsl_vector* x, gsl_vector* y);

double dot_product(gsl_vector* x, gsl_vector* y);

int vector_sum(gsl_vector* x, double a, gsl_vector* y, double b);

int vector_add(gsl_vector* x, double a);

int set_identity(gsl_matrix* A);

int mult_matrix_vector(gsl_matrix* A, gsl_vector* x, gsl_vector* b);

int mult_matrix_matrix(gsl_matrix* A, gsl_matrix* B, gsl_matrix* C);

int matrix_copy(gsl_matrix* A, gsl_matrix* B);

int vector_copy(gsl_vector* v, gsl_vector* u);


