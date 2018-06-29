#include "lin_eq.h"

int main(){

	int status = tall_matrix_test();
	if(status == -1) return -1;

	status = square_matrix_test();
	if(status == -1) return -1;

	return 0;
}
