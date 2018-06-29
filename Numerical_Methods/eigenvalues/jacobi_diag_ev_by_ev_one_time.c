#include "jacobi.h"

int main(int argc, char* argv[]){

        int n = atof(argv[1]);

        gsl_matrix* A = gsl_matrix_alloc(n,n);
        gsl_matrix* V = gsl_matrix_alloc(n,n);
        gsl_vector* e = gsl_vector_alloc(n);

        double rand_num = ((double) rand())/((double)RAND_MAX)*10-5;


        for(int i=0; i< A->size1; ++i){
                gsl_matrix_set(A, i, i, rand_num);
                for(int j = i+1; j < A->size1; ++j){
                        gsl_matrix_set(A, i, j, rand_num);
                        gsl_matrix_set(A, j, i, rand_num);
                }
        }

        int error = jacobi_diag_ev_by_ev(A, V, e, 1);
        if(error == -1) return -1;

        return 0;
}

