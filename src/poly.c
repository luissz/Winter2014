#include "poly.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

poly* psum(poly* A, poly* B){
	poly* sres=(poly*)malloc(sizeof(poly));

	if(A->deg > B->deg){
		int i;
		sres->deg = A->deg;
		sres->coef = (double*)malloc((sres->deg+1)*sizeof(double));
		for(i=0;i < sres->deg+1;i++){
			if(i < (sres->deg - B->deg))
				sres->coef[i]=A->coef[i];
			else
				sres->coef[i]=A->coef[i]+B->coef[i-(sres->deg - B->deg)];
		}
	}
	else if(A->deg < B->deg){
		int i;
		sres->deg = B->deg;
		sres->coef = (double*)malloc((sres->deg+1)*sizeof(double));
		for(i=0;i < sres->deg+1;i++){
			if(i < (sres->deg - A->deg))
				sres->coef[i]=B->coef[i];
			else
				sres->coef[i]=B->coef[i]+A->coef[i-(sres->deg - A->deg)];
		}
	}
	else{
		int i;
		sres->deg = A->deg;
		sres->coef = (double*)malloc((sres->deg+1)*sizeof(double));
		for(i=0;i < sres->deg+1;i++){
			sres->coef[sres->deg-i]=A->coef[A->deg-i]+B->coef[B->deg-i];
		}
	}

	return sres;
}

poly* pmult(poly* A, poly* B){
	int i;
	int j;
	int pos;
	poly* mres=(poly*)malloc(sizeof(poly));
	mres->deg = A->deg + B->deg;
	mres->coef = (double*)malloc((mres->deg+1)*sizeof(double));

	for(i=0;i < mres->deg+1;i++)
		mres->coef[i]=0.0;

	for(i=0;i < A->deg+1;i++){
		for(j=0;j < B->deg+1;j++){
			pos = (A->deg - i) + (B->deg - j);
			mres->coef[mres->deg - pos] += A->coef[i] * B->coef[j];
		}
	}

	return mres;
}

poly* ppow(poly* X, int Y){
	int i;
	poly* pres=(poly*)malloc(sizeof(poly));

	if(Y==0){
		pres->deg = 0;
		pres->coef = (double*)malloc(sizeof(double));
		pres->coef[0] = 1.0;
	}
	else{
		pres->deg = X->deg;
		pres->coef = (double*)malloc((pres->deg+1)*sizeof(double));
		for(i=0;i < pres->deg+1;i++)
			pres->coef[i] = X->coef[i];
	}

	for(i=1;i < Y;i++){
		pres = pmult(pres,X);
	}

	return pres;
}

poly* coef2poly(double coef){
	poly* aux = (poly*)malloc(sizeof(poly));

	aux->deg = 0;
	aux->coef = (double*)malloc(sizeof(double));
	aux->coef[0] = coef;

	return aux;
}

void displayPoly(poly* pnm){
	int i;
	for (i=0; i<pnm->deg+1; i++) {
		printf("%e s^%d ",pnm->coef[i],pnm->deg-i);
	}
	printf("\n");
}

double evalMultiVariable(poly* polyn, double* var_values){
	int i;
	double total;

	for(i = 0, total = 0.0;i < polyn->deg;i++){
		//printf("%lf x %lf\n",polyn->coef[i],(i==polyn->deg-1) ? 1.0 : var_values[i]);
		total += (i==polyn->deg-1) ? polyn->coef[i] : polyn->coef[i]*var_values[i];
	}

	return total;
}
