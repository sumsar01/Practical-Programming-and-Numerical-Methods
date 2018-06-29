#include <stdio.h>
#include "nvector.h"
#include <stdlib.h>
#include <math.h>
#include <assert.h>

nvector* nvector_alloc(int n){
	nvector* v = malloc(sizeof(nvector));
	(*v).size = n;
	(*v).data = malloc(n*sizeof(double));
	if( v==NULL ) fprintf(stderr,"error in nvector_alloc\n");
	return v;
}

void nvector_free(nvector* v){
	free(v->data);
	free(v);
}

void nvector_set(nvector* v, int i, double value){
	(*v).data[i] = value;
}

double nvector_get(nvector* v, int i){
	return v->data[i];
}

double nvector_dot_product(nvector* u, nvector*v)
{
	assert(u->size == v->size);
	for(int i = 0; i < u->size; i++)
	{
		double s = nvector_get(u, i) * nvector_get(v, i);
		nvector_set(u, i, s);
	}
}

void nvector_print(char *s, nvector* v){
	printf("%s",s);
	for(int i = 0; i < v->size; i++)
		printf("%g",v->data[i]);
	printf("\n");
}


