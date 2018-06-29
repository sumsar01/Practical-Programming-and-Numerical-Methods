#include "limits.h"
#include "float.h"
#include "stdio.h"
#include "tgmath.h"

int equal(double a, double b, double tau, double epsilon);

int main(){



printf("1.1)\n");


int i = 1;
while(i+1>i){i++;}
printf("my max int using while = %i\n",i);

int a=1;
for(a = 1; a+1>a;a++){}
printf("my max int using for = %i\n",a);

int b = 1;
do{b++;}
while(b+1>b);
printf("my max int using do = %i\n",b);

 printf("INT_MAX using limits.h is %i \n\n",INT_MAX);


printf("1.2)\n");

int c = 1;
while(c-1<c){c--;}
printf("My min int using while = %i\n",c);

int d = 1;
for(d=1; d-1<d; d--){}
printf("My min int using for = %i\n",d);

i = 1;
do{i--;}
while(i-1<i);
printf("My min int using do = %i\n\n",i);


printf("INT_MIN using limits.h is %i\n",INT_MIN);


printf("1.3)\n");

float x = 1;
while(1+x!=1){x/=2;}
x*=2;
printf("my float epsilon using while is = %g\n",x);

double y = 1;
while(1+y!=1){y/=2;}
y*=2;
printf("my double epsilon using while is = %g\n",y);

long double z = 1;
while(1+z!=1){z/=2;} z*=2;
printf("my long double epsilong using while is = %Lg\n",z);


 printf("FLT_EPSILON using flout.h is %g\n", FLT_EPSILON);

double e;
for(e=1; 1+e!=1; e/=2){} e*=2;
printf("my double epsilon using for is = %g\n",e);


double h = 1;
do
{h/=2;}
while(1+h!=1);
printf("my double epsilon using do is = %g\n",h);

double f = 1;
for(f = 1;1+f!=1;){f/=2;}
printf("my double epsilong using for = %g\n\n",f);


printf("2.1)\n");

int max = INT_MAX/3;
//int max = 100;

float sum_up_float = 0;
for(int i = 1; i <= max; i++){sum_up_float+=1.0f/i;}
printf("sum_up_float = %g\n",sum_up_float);


float sum_down_float = 0;
for(int i = max; i > 0; i--){sum_down_float += 1.0/i;}
printf("sum_down_float = %g\n", sum_down_float);


printf("\n2.2)\n");
printf("we get below the precicion of the float.");



printf("\n2.3)\n");
printf("The sum 1/x does not converge. The same goes for s1 and s2\n");

printf("\n2.4)\n");

double sum_up_double = 0;
for(int i = 1; i <= max; i++){sum_up_double += 1.0/i;}
printf("sum_up_double = %g\n", sum_up_double);

double sum_down_double = 0;
for(int i = max; i > 0; i--){sum_down_double += 1.0/i;}
printf("sum_down_double = %g\n", sum_down_double);



printf("\nProblem 3\n");
	equal(1,2,3,4);
	equal(100,1,1,1);

}

