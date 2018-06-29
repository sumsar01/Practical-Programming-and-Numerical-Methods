#include "tgmath.h"
#include "stdio.h"
#include "limits.h"

int equal(double a, double b, double tau, double epsilon){

if(abs(a-b) < tau)
{
	printf("1\n");
}
else if(fabs(a-b)/(fabs(a)+fabs(b)) < epsilon/2)
{
	printf("1\n");
}
else
{
printf("0\n");
}
}
