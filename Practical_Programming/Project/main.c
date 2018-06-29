#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double mylog(double);

int main(){
	for(double x = 0.1; x<10; x+=0.01)
		printf("%g %g\n",x,mylog(x));

	printf("\n\n");

	for(double x=0.1; x<10; x+=0.3)
		printf("%g %g\n",x,log(x));


return 0;
}
