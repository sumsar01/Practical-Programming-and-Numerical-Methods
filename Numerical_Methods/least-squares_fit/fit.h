#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

typedef struct {double (*function) (double, int); int nr_fun;} fit_function;

int QR_gs_decomp(gsl_matrix *A, gsl_matrix *R);

int QR_gs_solver(const gsl_matrix *Q, const gsl_matrix *R, const gsl_vector *b, gsl_vector *x);

int QR_gs_inverse(const gsl_matrix *R, gsl_matrix *B);

int lin_ls_sq_fit(const fit_function* fun, const gsl_vector* x, const gsl_vector* y, const 		gsl_vector* dy, gsl_vector* c, gsl_matrix* S);
