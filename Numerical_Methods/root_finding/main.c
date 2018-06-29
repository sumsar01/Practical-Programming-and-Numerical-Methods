#include "root_finding.h"

//test functions
void system_of_equations(gsl_vector* x, gsl_vector* fx){
	double A = 10000;
	double x_1 = gsl_vector_get(x,0);
	double x_2 = gsl_vector_get(x,1);
	double f_1 = A*x_1*x_2-1.;
	double f_2 = exp(-x_1)+exp(-x_2)-1.-1./A;
	gsl_vector_set(fx, 0, f_1);
	gsl_vector_set(fx, 1, f_2);
}


void rosenbrock_grad(gsl_vector* x, gsl_vector* fx){
	double x_1 = gsl_vector_get(x,0);
        double x_2 = gsl_vector_get(x,1);
        double f_1 = -2*(1-x_1)-400*(x_2-x_1*x_1)*x_1;
        double f_2 = 200*(x_2-x_1);
        gsl_vector_set(fx, 0, f_1);
        gsl_vector_set(fx, 1, f_2);
}

void rosenbrock_grad_jacobian(gsl_vector* x, gsl_matrix* J){
	double x_1 = gsl_vector_get(x,0);
        double x_2 = gsl_vector_get(x,1);
	double J_11 = 2-400*(x_2-x_1*x_1)+800*x_1*x_1;
	double J_12 = -400*x_1;
	double J_21 = -400*x_1;
	double J_22 = 200;
	gsl_matrix_set(J,0,0,J_11);
	gsl_matrix_set(J,0,1,J_12);
	gsl_matrix_set(J,1,0,J_21);
	gsl_matrix_set(J,1,1,J_22);
}


void himmelblau_grad(gsl_vector* x, gsl_vector* fx){
	double x_1 = gsl_vector_get(x,0);
        double x_2 = gsl_vector_get(x,1);
        double f_1 = 4*(x_1*x_1+x_2-11)*x_1+2*(x_1+x_2*x_2-7);
        double f_2 = 2*(x_1*x_1+x_2-11)+4*(x_1+x_2*x_2-7)*x_2;
        gsl_vector_set(fx, 0, f_1);
        gsl_vector_set(fx, 1, f_2);
}

void my_eq(gsl_vector* x, gsl_vector* fx){
	double x_1 = gsl_vector_get(x,0);
	double f_1 = cos(x_1) - x_1;
	gsl_vector_set(fx, 0, f_1);
}

void my_eq_jacobian(gsl_vector* x, gsl_matrix* J){
	double x_1 = gsl_vector_get(x,0);
	double J_1 = -sin(x_1) -1;
	gsl_matrix_set(J,0,0,J_1);
}

//Newton method with numerical Jacobian testing
int numerical_test(void){
	gsl_vector* x =gsl_vector_alloc(2);
	gsl_vector_set(x,0,0);
	gsl_vector_set(x,1,1);

	int num_iter = newton_num(&system_of_equations, x, 0.001, 1e-6);
	if(num_iter == -1) return -1;

	printf("\nTesting Newton method with numerical Jacobian.\n");
	printf("Solving A*x*y = 1, exp(-x) + exp(-y) = 1+ 1/A where A = 10000\n");
	printf("The solution is found to be [x,y] = \n");
	print_vector(x);
	printf("Converging in %i steps.\n", num_iter);


	gsl_vector_set(x,0,0);
	gsl_vector_set(x,1,0);

        num_iter = newton_num(&rosenbrock_grad, x, 0.001, 1e-6);
        if(num_iter == -1) return -1;

        printf("\nSolving for the minimum of Rosenbrock fucntion.\n");
        printf("f(x,y) = (1-x)²+100*(y-x²)²\n");
        printf("The solution is found to be [x,y] = \n");
        print_vector(x);
        printf("Converging in %i steps.\n", num_iter);


        gsl_vector_set(x,0,-3);
        gsl_vector_set(x,1,3);

        num_iter = newton_num(&himmelblau_grad, x, 0.001, 1e-6);
        if(num_iter == -1) return -1;

        printf("\nSolving for the minimum of Himmelblau fucntion.\n");
        printf("f(x,y) = (x²+y-11)²+(x+y²-7)²\n");
        printf("The solution is found to be [x,y] = \n");
        print_vector(x);
        printf("Converging in %i steps.\n", num_iter);


	gsl_vector_free(x);
	x = gsl_vector_alloc(1);
	gsl_vector_set(x,0,1);

        num_iter = newton_num(&my_eq, x, 0.001, 1e-6);
        if(num_iter == -1) return -1;

        printf("\nSolving cos(x) = x.\n");
        printf("The solution is found to be x = \n");
        print_vector(x);
        printf("Converging in %i steps.\n", num_iter);

        gsl_vector_free(x);

	return 0;
}



//Newton method with analytical Jacobian testing

int analytical_test(void){

	gsl_vector* x = gsl_vector_alloc(2);
	gsl_vector_set(x,0,0);
	gsl_vector_set(x,1,0);


        int num_iter = newton_jacobian(&rosenbrock_grad,&rosenbrock_grad_jacobian,
		 x, 1e-6);
        if(num_iter == -1) return -1;

        printf("\nTesting Newton method with analytical Jacobian.\n");
        printf("\nSolving for the minimum of Rosenbrock fucntion.\n");
        printf("f(x,y) = (1-x)²+100*(y-x²)²\n");
        printf("The solution is found to be [x,y] = \n");
        print_vector(x);
        printf("Converging in %i steps.\n", num_iter);


        gsl_vector_free(x);
        x = gsl_vector_alloc(1);
        gsl_vector_set(x,0,1);

	num_iter = newton_jacobian(&my_eq,&my_eq_jacobian,
                 x, 1e-6);
        if(num_iter == -1) return -1;

        printf("\nSolving cos(x) = x.\n");
        printf("The solution is found to be x = \n");
        print_vector(x);
        printf("Converging in %i steps.\n", num_iter);

        gsl_vector_free(x);

	return 0;
}

int main(void){

	int error = numerical_test();
	if(error == -1) return -1;

	error = analytical_test();
	if(error == -1) return -1;

	return 0;
}



















