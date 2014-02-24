#ifndef POLY_LIB
#define POLY_LIB
#include <math.h>

typedef struct{
	int deg;
	double* coef;
}poly;

poly* psum(poly* A, poly* B);
poly* pmult(poly* A, poly* B);
poly* ppow(poly* X, int Y);
poly* coef2poly(double coef);
void displayPoly(poly* pnm);
double evalMultiVariable(poly* polyn, double* var_values);

#endif
