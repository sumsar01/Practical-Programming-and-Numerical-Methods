#include"minimization.h"

double rosenbrock(gsl_vector* x){
	double x_1 = gsl_vector_get(x, 0);
        double x_2 = gsl_vector_get(x, 1);
	return pow(1 - x_1, 2) + 100*pow(x_2 - x_1*x_1, 2);
}

void rosenbrock_grad(gsl_vector* x, gsl_vector* fx){
        double x_1 = gsl_vector_get(x, 0);
        double x_2 = gsl_vector_get(x, 1);
        double f_1 = -2*(1-x_1) - 400*(x_2-x_1*x_1)*x_1;
        double f_2 = 200*(x_2-x_1*x_1);
        gsl_vector_set(fx, 0, f_1);
        gsl_vector_set(fx, 1, f_2);
}

void rosenbrock_hessian(gsl_vector* x, gsl_matrix* H){
        double x_1 = gsl_vector_get(x,0);
        double x_2 = gsl_vector_get(x,1);
        double H_11 = 2-400*(x_2-x_1*x_1)+800*x_1*x_1;
        double H_12 = -400*x_1;
        double H_21 = -400*x_1;
        double H_22 = 200;
        gsl_matrix_set(H,0,0,H_11);
        gsl_matrix_set(H,0,1,H_12);
        gsl_matrix_set(H,1,0,H_21);
        gsl_matrix_set(H,1,1,H_22);
}


double himmelblau(gsl_vector* x){
	double x_1 = gsl_vector_get(x, 0);
	double x_2 = gsl_vector_get(x, 1);
	return pow(x_1*x_1 + x_2 -11, 2) + pow(x_1 + x_2*x_2 -7, 2);
}

void himmelblau_grad(gsl_vector* x, gsl_vector* fx){
	double x_1 = gsl_vector_get(x, 0);
	double x_2 = gsl_vector_get(x, 1);
	double f_1 = 4*(x_1*x_1 + x_2 - 11)*x_1 + 2*(x_1 + x_2*x_2 -7);
	double f_2 = 2*(x_1*x_1 + x_2 - 11) + 4*(x_1 + x_2*x_2 - 7)*x_2;
	gsl_vector_set(fx, 0, f_1);
	gsl_vector_set(fx, 1, f_2);
}

void himmelblau_hessian(gsl_vector* x, gsl_matrix* H){
	double x_1 = gsl_vector_get(x, 0);
	double x_2 = gsl_vector_get(x, 1);
	double H_11 = 4 * (x_1 * x_1 + x_2 - 11) + 8*x_1*x_1+2;
	double H_12 = 4 * x_1 + 4 * x_2;
	double H_21 = 4 * x_1 + 4 * x_2;
	double H_22 = 4 * (x_1 + x_2 * x_2 - 7) + 8 * x_2 * x_2 + 2;
	gsl_matrix_set(H, 0, 0, H_11);
	gsl_matrix_set(H, 0, 1, H_12);
	gsl_matrix_set(H, 1, 0, H_21);
	gsl_matrix_set(H, 1, 1, H_22);
}

double fit_fun(gsl_vector* x){
	double x_1 = gsl_vector_get(x, 0);
	double x_2 = gsl_vector_get(x, 1);
	double x_3 = gsl_vector_get(x, 2);

	double t[] = {0.23,1.29,2.35,3.41,4.47,5.53,6.59,7.65,8.71,9.77};
	double y[] = {4.64,3.38,3.01,2.55,2.29,1.67,1.59,1.69,1.38,1.46};
	double e[] = {0.42,0.37,0.34,0.31,0.29,0.27,0.26,0.25,0.24,0.24};
	int N = sizeof(t)/sizeof(t[0]);

	double F = 0.;

	for(int i = 0; i < N; ++i){
		F += pow((x_1*exp(-t[i]/x_2) + x_3 - y[i])/e[i], 2);
	}

	return F;
}


//testing Newton for minimum with gradient and hessian given by user.
int tesing_min_newton(void){

	gsl_vector* xstart = gsl_vector_alloc(2);

	gsl_vector_set(xstart, 0, 0);
	gsl_vector_set(xstart, 1, 0);

	int iterations = min_newton(&rosenbrock, &rosenbrock_grad, &rosenbrock_hessian, xstart, 1e-4);
	if(iterations == -1) return -1;

	printf("\nFinding minimum of Rosenbrock function f(x,y) = (1-x)^2+100*(y-x^2)^2 \n");
	printf("The solution is found to be (x,y) = \n");
	print_vector(xstart);
	printf("converging in %i steps.\n", iterations);

	gsl_vector_set(xstart, 0, -3);
	gsl_vector_set(xstart, 1, 3);

	iterations = min_newton(&himmelblau, &himmelblau_grad, &himmelblau_hessian, xstart, 1e-4);
	if(iterations == -1) return -1;

        printf("\nFinding minimum of Himmelblau function f(x,y) = (x^2+y-11)^2 + (x+y^2-7)^2 \n");
        printf("The solution is found to be (x,y) = \n");
        print_vector(xstart);
        printf("converging in %i steps.\n", iterations);

	gsl_vector_free(xstart);

	return 0;
}



//testing Newton for minimum with gradient calculated numerically and hessian updated.
int testing_min_newton_num(void){

	gsl_vector* xstart = gsl_vector_alloc(2);

	gsl_vector_set(xstart, 0, 0);
	gsl_vector_set(xstart, 1, 0);

	int iterations = min_newton_num(&rosenbrock, xstart, 1e-6, 1e-4);
	if(iterations == -1) return -1;

        printf("\nFinding minimum of Rosenbrock function f(x,y) = (1-x)^2+100*(y-x^2)^2 \n");
        printf("The solution is found to be (x,y) = \n");
        print_vector(xstart);
        printf("converging in %i steps.\n", iterations);

        gsl_vector_set(xstart, 0, -3);
        gsl_vector_set(xstart, 1, 3);

	iterations = min_newton_num(&himmelblau, xstart, 1e-6, 1e-4);
        if(iterations == -1) return -1;

        printf("\nFinding minimum of Himmelblau function f(x,y) = (x^2+y-11)^2 + (x+y^2-7)^2 \n");
        printf("The solution is found to be (x,y) = \n");
        print_vector(xstart);
        printf("converging in %i steps.\n", iterations);

        gsl_vector_free(xstart);

        return 0;
}


int testing_min_newton_fit(void){

	gsl_vector* xstart = gsl_vector_alloc(3);
	gsl_vector_set(xstart, 0, 1);
	gsl_vector_set(xstart, 1, 1);
	gsl_vector_set(xstart, 2, 1);

	int iterations = min_newton_num(&fit_fun, xstart, 1e-6, 1e-4);
	if(iterations == -1) return -1;

	printf("\nFitting function to data \n");
	printf("The solution is found to be (A,T,B) = \n");
	print_vector(xstart);
	printf("converging in %i steps.\n\n",iterations);

	FILE* file = fopen("data.txt","w");

        double t[] = {0.23,1.29,2.35,3.41,4.47,5.53,6.59,7.65,8.71,9.77};
        double y[] = {4.64,3.38,3.01,2.55,2.29,1.67,1.59,1.69,1.38,1.46};
        double e[] = {0.42,0.37,0.34,0.31,0.29,0.27,0.26,0.25,0.24,0.24};
        int N = sizeof(t)/sizeof(t[0]);

	for (int i = 0; i  < N; ++i){
		fprintf(file, "%g %g %g\n", t[i], y[i], e[i]);
	}
	fprintf(file, "\n\n");
	double x_1 = gsl_vector_get(xstart, 0);
	double x_2 = gsl_vector_get(xstart, 1);
	double x_3 = gsl_vector_get(xstart, 2);
	for( double x = t[0]; x < t[N-1]+1e-5; x+=0.1){
		fprintf(file, "%g %g\n", x, x_1*exp(-x/x_2) + x_3);
	}

	gsl_vector_free(xstart);

	return 0;
}


int main(void){

	int error = tesing_min_newton();
	if(error == -1) return -1;

	error = testing_min_newton_num();
	if(error == -1) return -1;

	error = testing_min_newton_fit();
	if(error == -1) return -1;

	return 0;
}
























