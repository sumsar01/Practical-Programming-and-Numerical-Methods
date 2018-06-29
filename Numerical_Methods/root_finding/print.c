#include "root_finding.h"

void print_matrix(gsl_matrix* A){
	for(int i=0; i < A->size1; ++i){
		for(int j=0; j < A->size2; ++j){
			double A_ij = gsl_matrix_get(A, i, j);
			A_ij = round(A_ij*1e12)*1e-12;
			printf("%.2g\t",A_ij);
		}
		printf("\n");
	}

}

void print_vector(gsl_vector* x){
	for(int i=0; i< x->size; ++i){
		double x_i = gsl_vector_get(x,i);
		x_i = round(x_i*1e12)*1e-12;
		printf("%.5g\n",x_i);
	}
}
