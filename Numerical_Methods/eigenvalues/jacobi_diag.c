#include "jacobi.h"


//jacobi diagonalization: upper triangle of A is destroyed;
//e and V accumulate eigenvalues and eigenvector
/*
int jacobi(gsl_matrix* A, gsl_vector* e, gsl_matrix* V){
assert(A->size1==A->size2);
assert(V->size1==V->size2);
assert(A->size1==V->size1);
assert(A->size2==V->size2);


	int changed, sweeps = 0, n = A->size1;

	for(int i=0;i<n;i++) gsl_vector_set(e,i,gsl_matrix_get(A,i,i));
		gsl_matrix_set_identity(V);

	do{changed=0; sweeps++; int p, q;



		for(p=0;p<n;p++)for(q=p+1;q<n;q++){
			double app=gsl_vector_get(e,p);
			double aqq=gsl_vector_get(e,q);
			double apq=gsl_matrix_get(A,p,q);
			double phi=0.5*atan2(2*apq,aqq-app);
			double c=cos(phi),s=sin(phi);
				double app1=c*c*app-2*s*c*apq+s*s*aqq;
				double aqq1=s*s*app+2*s*c*apq+c*c*aqq;

			if(app1!=app || aqq1!=aqq){changed=1;
				gsl_vector_set(e,p,app1);
				gsl_vector_set(e,q,aqq1);
				gsl_matrix_set(A,p,q,0.0);

				for(int i=0;i<p;i++){
					double aip=gsl_matrix_get(A,i,p);
					double aiq=gsl_matrix_get(A,i,q);
					gsl_matrix_set(A,i,p,c*aip-s*aiq);
					gsl_matrix_set(A,i,q,c*aiq+s*aip);
				}

				for(int i=p+1;i<q;i++){
                                        double api=gsl_matrix_get(A,p,i);
                                        double aiq=gsl_matrix_get(A,q,i);
                                        gsl_matrix_set(A,p,i,c*api-s*aiq);
                                        gsl_matrix_set(A,i,q,c*aiq+s*api);
				}

				for(int i=p+1;i<n;i++){
                                        double api=gsl_matrix_get(A,p,i);
                                        double aqi=gsl_matrix_get(A,q,i);
                                        gsl_matrix_set(A,p,i,c*api-s*aqi);
                                        gsl_matrix_set(A,q,i,c*aqi+s*api);
                                }

				for(int i=0;i<n;i++){
                                        double vip=gsl_matrix_get(V,i,p);
                                        double viq=gsl_matrix_get(V,i,q);
                                        gsl_matrix_set(V,i,p,c*vip-s*viq);
                                        gsl_matrix_set(V,i,q,c*viq+s*vip);
                                }
			}
		}
	}

	while(changed!=0);

	return sweeps;
}

*/


	//borrowed from Casper Poulsen.
int jacobi(gsl_matrix* A, gsl_matrix* V, gsl_vector* diag, int p, int q){
    //Since A is symmetric, we only sweep the upper triangular part,
    //and keep the original lower triangular part. So we let A come
    //from the lower part and A'=Am be the upper part.

 //Tests
    if (A->size1 != A->size2){
        fprintf(stderr, "Error in jacobi_diag_rotation: A must be square.\n");
        return -1;
    }
    if (q <= p){
        fprintf(stderr, "Error in jacobi_diag_rotation: Must have q>p\n");
        return -1;
    }

    //Running Algorithm
    int n = A->size1;

    double A_pq = gsl_matrix_get(A, p, q);
    double A_pp = gsl_vector_get(diag, p);
    double A_qq = gsl_vector_get(diag, q);
    double phi = atan2(2*A_pq, (A_qq - A_pp)) / 2;
    double c = cos(phi);
    double s = sin(phi);

    double Am_pp = pow(c,2)*A_pp - 2*s*c*A_pq + pow(s,2)*A_qq;
    double Am_qq = pow(s,2)*A_pp + 2*s*c*A_pq + pow(c,2)*A_qq;
    double Am_pq = s*c*(A_pp - A_qq) + (pow(c,2) - pow(s,2))*A_pq;

    if(abs(Am_pq) > 1e-16){
        fprintf(stderr, "Error in jacobi_diag_rotation: Jacobi rotation failed to set A_pq = 0");
        return -1;
    }

    gsl_vector_set(diag, p, Am_pp);
    gsl_vector_set(diag, q, Am_qq);
    gsl_matrix_set(A, p, q, Am_pq);


    for (int i = 0; i < p; ++i) {
        double A_pi = gsl_matrix_get(A, i, p);
        double A_qi = gsl_matrix_get(A, i, q);
        double Am_pi = c*A_pi - s*A_qi;
        double Am_qi = s*A_pi + c*A_qi;
        gsl_matrix_set(A, i, p, Am_pi);
        gsl_matrix_set(A, i, q, Am_qi);

    }

    for (int i = p+1; i < q; ++i) {
        double A_pi = gsl_matrix_get(A, p, i);
        double A_qi = gsl_matrix_get(A, i, q);
        double Am_pi = c*A_pi - s*A_qi;
        double Am_qi = s*A_pi + c*A_qi;
        gsl_matrix_set(A, p, i, Am_pi);
        gsl_matrix_set(A, i, q, Am_qi);

    }

    for (int i = q+1; i < n; ++i) {
        double A_pi = gsl_matrix_get(A, p, i);
        double A_qi = gsl_matrix_get(A, q, i);
        double Am_pi = c*A_pi - s*A_qi;
        double Am_qi = s*A_pi + c*A_qi;
        gsl_matrix_set(A, p, i, Am_pi);
        gsl_matrix_set(A, q, i, Am_qi);
    }

    //Refreces V.
    for (int i = 0; i < n; ++i) {
        double V_ip = gsl_matrix_get(V, i, p);
        double V_iq = gsl_matrix_get(V, i, q);
        double Vm_ip = c*V_ip - s*V_iq;
        double Vm_iq = s*V_ip + c*V_iq;
        gsl_matrix_set(V, i, p, Vm_ip);
        gsl_matrix_set(V, i, q, Vm_iq);
    }

    return 0;
}





int jacobi_conv_test(gsl_matrix* A){
	double maximum = 0.;
	int n = A->size1;

	for(int i = 0; i<n; ++i){
		for(int j = i+1; j<n; ++j){
			double A_ij = abs(gsl_matrix_get(A, i, j));
			if(A_ij > maximum) maximum = A_ij;
		}
	}

	if(maximum < 1e-12) return 0;
	else return 1;
}

int jacobi_diag_sweep(gsl_matrix* A, gsl_matrix* V, gsl_vector* e){
	int n = A->size1;

	//setting up eigenvalues as diagonal of A.
	for(int i=0; i<n; ++i){
		gsl_vector_set(e, i, gsl_matrix_get(A, i, i));
	}


	//sets V=I
/*	for(int i=0; i < n; ++i){
		gsl_matrix_set(V, i, i, 1.);
		for(int j = i+1; j < n; ++j){
			gsl_matrix_set(V, i, j, 0);
			gsl_matrix_set(V, j, i, 0);
		}
	}

*/
	gsl_matrix_set_identity(V);



	//algorithm
	int result, converged, count = 0;
	do{
		for(int i = 0; i<n-1; ++i){
			for(int j = i+1; j<n; ++j){
				result = jacobi(A, V, e, i, j);
				if(result == -1) return -1;
			}
		}

		converged = jacobi_conv_test(A);
		if(converged == -1) return -1;
		count++;
	}

	while(converged == 1 && count < 1000*n);



	if(count == 1000*n){
		fprintf(stderr,"Error jacobi did not converge within 1000*n sweeps.");
		return 1;
	}


	//restor A
	for(int i = 0; i<n; ++i){
		for(int j = i+1; j<n; ++j){
			double A_ij = gsl_matrix_get(A, j, i);
			gsl_matrix_set(A, i, j, A_ij);
		}
	}

	return 0;
}


int jacobi_conv_test_single_ev(gsl_matrix* A, int num_of_ev){
	assert(A->size1 = A->size2);

	double maximum = 0.;
	int n = A->size1;

	for(int j = num_of_ev+1; j<n; ++j){
		double A_ij = abs(gsl_matrix_get(A, num_of_ev, j));
		if(A_ij > maximum) maximum = A_ij;
	}

	if(maximum < 1e-12) return 0;
	else return 1;
}


int jacobi_diag_ev_by_ev
	(gsl_matrix* A, gsl_matrix* V, gsl_vector* e, int num_of_ev){

	//tests
	assert(A->size1 = A->size2);
	assert(V->size1 = V->size2);
	assert(V->size2 = A->size1);
	assert(e->size = V->size1);

	int n = A->size1;

	//set eigenvalues as diagonal of A
	for(int i = 0; i<n; ++i){
		gsl_vector_set(e, i, gsl_matrix_get(A, i, i));
	}

	//set V = I
	gsl_matrix_set_identity(V);

	//Algorithm
	int result, converged, count = 0;
	for(int i=0; i<num_of_ev;++i){
		do{
			for(int j=i+1; j<n; ++j){
				result = jacobi(A, V, e, i, j);
				if(result == -1) return -1;
			}
			converged = jacobi_conv_test_single_ev(A, i);
			if(converged == -1) return -1;
			count++;
		}while(converged == 1 && count < 1000*n);
	}

	if(count == 1000*n){
		fprintf(stderr,"jacobi_diag_sweep did not converge within 1000*n sweeps\n");
		return -1;
	}

	//set unusable values of eigenvalues and V to 0
	for(int i = num_of_ev; i<n; ++i){
		gsl_vector_set(e, i, 0);
		for(int j=0; j<n; ++j){
			gsl_matrix_set(V, j, i, 0);
		}
	}

	//restores A
	for(int i = 0; i<n; ++i){
		for(int j = i+1; j<n; ++j){
			double A_ji = gsl_matrix_get(A, j, i);
			gsl_matrix_set(A, i, j, A_ji);
		}
	}

	return 0;
}














