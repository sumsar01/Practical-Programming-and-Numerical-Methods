#include "fit.h"


void data_1(gsl_vector* x, gsl_vector* y, gsl_vector* dy){

	gsl_vector_set(x, 0, 0.1);
	gsl_vector_set(x, 1, 1.33);
	gsl_vector_set(x, 2, 2.55);
	gsl_vector_set(x, 3, 3.78);
	gsl_vector_set(x, 4, 5.);
	gsl_vector_set(x, 5, 6.22);
	gsl_vector_set(x, 6, 7.45);
	gsl_vector_set(x, 7, 8.68);
	gsl_vector_set(x, 8, 9.9);

	gsl_vector_set(y, 0, -15.3);
        gsl_vector_set(y, 1, 0.32);
        gsl_vector_set(y, 2, 2.45);
        gsl_vector_set(y, 3, 2.75);
        gsl_vector_set(y, 4, 2.27);
        gsl_vector_set(y, 5, 1.35);
        gsl_vector_set(y, 6, 0.157);
        gsl_vector_set(y, 7, -1.23);
        gsl_vector_set(y, 8, -2.75);

        gsl_vector_set(dy, 0, 1.04);
        gsl_vector_set(dy, 1, 0.594);
        gsl_vector_set(dy, 2, 0.983);
        gsl_vector_set(dy, 3, 0.998);
        gsl_vector_set(dy, 4, 1.11);
        gsl_vector_set(dy, 5, 0.398);
        gsl_vector_set(dy, 6, 0.535);
        gsl_vector_set(dy, 7, 0.968);
        gsl_vector_set(dy, 8, 0.478);

}

double fitting_fun_1(double x, int i){
	switch(i){
		case 0: return log(x); break;
		case 1: return 1.0; break;
		case 2: return x; break;
		default: {fprintf(stderr,"Error in fit_fun: wrong i:%d",i); return NAN;}
	}
}

int test_fit_1(void){

	gsl_vector* x = gsl_vector_alloc(9);
	gsl_vector* y = gsl_vector_alloc(9);
	gsl_vector* dy = gsl_vector_alloc(9);
	data_1(x, y, dy);

	fit_function fun;
	fun.function = fitting_fun_1;
	fun.nr_fun = 3;

	gsl_vector* c = gsl_vector_alloc(3);
	gsl_matrix* S = gsl_matrix_alloc(3,3);

	int error = lin_ls_sq_fit(&fun, x, y, dy, c, S);
	if(error == -1) return -1;


	FILE* file = fopen("fit1.txt","w");

	for(int i=0; i < x->size; ++i){
		fprintf(file,"%g %g %g\n", gsl_vector_get(x,i), gsl_vector_get(y,i),
			gsl_vector_get(dy, i));
	}

	fprintf(file, "\n\n");
	double z_min = gsl_vector_get(x, 0);
	double z_max = gsl_vector_get(x, (x->size)-1)+1e-5;
	double dz = (z_max - z_min)/1000;

	for(double z = z_min; z < z_max; z += dz){
		double f = 0.;
		double df = 0.;
		for( int i = 0; i < c->size; ++i){
			double c_i = gsl_vector_get(c, i);
			for(int j=0; j< c->size; ++j){
				double S_ij = gsl_matrix_get(S, i, j);
				df += S_ij*fitting_fun_1(z, i)*fitting_fun_1(z,j);
			}
			f += c_i*fitting_fun_1(z, i);
		}
		fprintf(file, "%g %g %g\n", z, f, sqrt(df));
	}

	return 0;
}


void data_2(gsl_vector* x, gsl_vector* y, gsl_vector* dy){

        gsl_vector_set(x, 0, 0.1);
        gsl_vector_set(x, 1, 0.145);
        gsl_vector_set(x, 2, 0.211);
        gsl_vector_set(x, 3, 0.307);
        gsl_vector_set(x, 4, 0.447);
        gsl_vector_set(x, 5, 0.649);
        gsl_vector_set(x, 6, 0.944);
        gsl_vector_set(x, 7, 1.372);
        gsl_vector_set(x, 8, 1.995);
        gsl_vector_set(x, 9, 2.900);


        gsl_vector_set(y, 0, 12.644);
        gsl_vector_set(y, 1, 9.235);
        gsl_vector_set(y, 2, 7.377);
        gsl_vector_set(y, 3, 6.460);
        gsl_vector_set(y, 4, 5.555);
        gsl_vector_set(y, 5, 5.896);
        gsl_vector_set(y, 6, 5.673);
        gsl_vector_set(y, 7, 6.964);
        gsl_vector_set(y, 8, 8.896);
        gsl_vector_set(y, 9, 11.355);


        gsl_vector_set(dy, 0, 0.858);
        gsl_vector_set(dy, 1, 0.359);
        gsl_vector_set(dy, 2, 0.505);
        gsl_vector_set(dy, 3, 0.403);
        gsl_vector_set(dy, 4, 0.683);
        gsl_vector_set(dy, 5, 0.605);
        gsl_vector_set(dy, 6, 0.856);
        gsl_vector_set(dy, 7, 0.351);
        gsl_vector_set(dy, 8, 1.083);
        gsl_vector_set(dy, 9, 1.002);

}


double fitting_fun_2(double x, int i){
        switch(i){
                case 0: return 1/x; break;
                case 1: return 1.0; break;
                case 2: return x; break;
                default: {fprintf(stderr,"Error in fit_fun: wrong i:%d",i); return NAN;}
        }
}


int test_fit_2(void){

	gsl_vector* x = gsl_vector_alloc(10);
        gsl_vector* y = gsl_vector_alloc(10);
        gsl_vector* dy = gsl_vector_alloc(10);
        data_2(x, y, dy);


        fit_function fun;
        fun.function = fitting_fun_2;
        fun.nr_fun = 3;

        gsl_vector* c = gsl_vector_alloc(3);
        gsl_matrix* S = gsl_matrix_alloc(3,3);

        int error = lin_ls_sq_fit(&fun, x, y, dy, c, S);
        if(error == -1) return -1;

        FILE* file = fopen("fit2.txt","w");

        for(int i=0; i < x->size; ++i){
                fprintf(file,"%g %g %g\n", gsl_vector_get(x,i),
			gsl_vector_get(y,i),gsl_vector_get(dy, i));
        }

        fprintf(file, "\n\n");
        double z_min = gsl_vector_get(x, 0);
	double z_max = gsl_vector_get(x, (x->size)-1)+1e-5;
        double dz = (z_max - z_min)/1000;

       for(double z = z_min; z < z_max; z += dz){
                double f = 0.;
                double df = 0.;
                for( int i = 0; i < c->size; ++i){
                        double c_i = gsl_vector_get(c, i);
                        for(int j=0; j< c->size; ++j){
                                double S_ij = gsl_matrix_get(S, i, j);
                                df += S_ij*fitting_fun_2(z, i)*fitting_fun_2(z,j);
                        }
                        f += c_i*fitting_fun_2(z, i);
                }
                fprintf(file, "%g %g %g\n", z, f, sqrt(df));
        }

        return 0;
}


int main(void){

	int error = test_fit_1();
	if(error == -1) return -1;

	error = test_fit_2();
	if(error == -1) return -1;

	return 0;
}












