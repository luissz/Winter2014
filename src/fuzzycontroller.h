//Fuzzy Control library
#ifndef FUZZY_LIB
#define FUZZY_LIB
#include <stdio.h>
#include "system.h"

#define ZERO_TH 0.000000001

//Redefinition of the MF type: Uses four
//parameters, defining a trapezoid.
//Down means u(x) = 0.0
//Up means u(x) = 1.0
typedef struct{
	double downleft;
	double upleft;
	double upright;
	double downright;
	double center;
	double core;
	double uplim;
	double lolim;
}memfunc;

/*typedef struct{
	//functype type;
	double center;
	double core;
	double uplim;
	double lolim;
}memfunc;*/

typedef struct{
	int ndims;
	//int* udivs;
	int* utype;//union type (for each rule)
	int* out_funcs;
}famtable;

typedef struct{
	double x;
	double y;
}vertex;

typedef struct{
	int nvert;
	vertex* varr;
	double area;
}polygon;

typedef struct{
	int ninps;
	int nouts;
	double* in_vars;
	double* out_vars;
	int* udivs;
	memfunc** func_mat;
	double** fuzzy_in;
	double** fuzzy_out;
	famtable* fam;
	polygon* conseq;
}fuzzycontroller;

typedef struct{
	int ninps;
	int nouts;
	int n_rules;
	double* in_vars;
	double* out_vars;
	int* udivs;
	memfunc** func_mat;
	double** fuzzy_in;
	double** fuzzy_out;
	famtable* fam;
	poly* rules;
}fuzzyTScontroller;

void fuzzify(int ninps, int* ndivs, double* crisp, memfunc** func_mat, double** fuzzy);
void fuzzifyTS(fuzzyTScontroller* ts_fuzzy);
void infere2_1(int* udivs, double** ifuzzy, famtable* fam, double** ofuzzy);
void infereTS(fuzzyTScontroller* ts_fuzzy);
double max(double a, double b);
double min(double a, double b);
polygon* createUnionPolygon(int nouts, int* udivs, memfunc** func_mat, double** fuzzy, char* inference);
polygon* mamdaniPolygon(memfunc* out_func, double fuzzy);
polygon* initializePolygon(int verts);
polygon* mamdaniUnion(polygon* a, polygon* b);
vertex* linesIntersect(vertex* a1, vertex* a2, vertex* b1, vertex* b2);
void polygonArea(polygon* poly);
void defuzzify(char* method, polygon* poly, int nouts, double* out_vars);
void defuzzifyTS(fuzzyTScontroller* ts_fuzzy);
void clearController(fuzzycontroller* fc);
void clearTSController(fuzzyTScontroller* fc);
void initializeController(fuzzycontroller* fuzzyc, FILE* fh);
void initializeTSController(fuzzyTScontroller* fuzzyc, FILE* fh);
//void initializeProcess(dsystem* proc, FILE* fh);
#endif
