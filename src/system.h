#ifndef SYSTEM_LIB
#define SYSTEM_LIB
#include "poly.h"
#include <stdlib.h>
#include <stdio.h>

typedef struct{
	poly num;
	poly den;
}ctf;

typedef struct{
	poly num;
	poly den;
	double ts;
}dtf;

typedef struct{
	dtf tfunc;
	int dstate;
	double* in;
	double* out;
}dsystem;

void readTfData(FILE* fh, poly* N, poly* D);
void readTs(FILE* fh, dtf* discrete);
void initBilinear(poly* num, poly* den, double ts);
void subztitution(ctf* con, dtf* dis);
void toNegPows(dtf* dis);
void contSystemMult(ctf* a, ctf* b);
void discSystemMult(dtf* a, dtf* b);
void unityNegFeedback(ctf* sys);
void contZOH(poly* num, poly* den);
void discZOH(poly* num, poly* den);
void cleanCSystem(ctf* sys);

void evalSystemStep(dsystem* sys);
void evalSystemResponse(dsystem* sys, unsigned int nlim);

void displayCSystem(ctf* sys);
void displayDSystem(dtf* sys);

void initializeProcess(dsystem* proc, FILE* fh);
#endif
