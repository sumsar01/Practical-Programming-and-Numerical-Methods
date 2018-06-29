#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <assert.h>

int jacobi_diag(gsl_matrix* A, gsl_matrix* V, gsl_vector* e);

int jacobi_diag_sweep(gsl_matrix* A, gsl_matrix* V, gsl_vector* e);

int jacobi_diag_ev_by_ev(gsl_matrix* A, gsl_matrix* V, gsl_vector* e, int num_of_ev);

void print_matrix(gsl_matrix* A);

void print_vector(gsl_vector* x);


