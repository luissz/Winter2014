//Fuzzy TS controller
/*
Calculates the output of a TS fuzzy system
Usage:
tscontroller ts_system_filename inputs_filename
*/
#include "fuzzycontroller.h"
#include <stdlib.h> 
#include <string.h>

#define MAX 50

int main(int argc, char** argv){
	FILE* fh;
	fuzzyTScontroller ts_example1;
	//double y;
	double** x;
	int i;
	int j;
	int c;
	int d;
	int ns;
	char dummy[MAX];

	if(argc != 3){
		printf("%s\n",(argc < 3)?"Not enough input arguments":"Too many input arguments");
		exit(1);
	}

//------TS system Initialization-----
	fh = fopen(argv[1],"r");//Controller filename

	if(fh == NULL){
		printf("Fuzzy System file not found!\n");
		exit(1);
	}

	initializeTSController(&ts_example1,fh);

	fclose(fh);
//------------------------------------

//------Simulation setup--------------
	fh = fopen(argv[2],"r");//Inputs filename

	ns = -1;
	while(!feof(fh)){
		fgets(dummy,MAX,fh);
		ns++;
	}
	fclose(fh);

	//x = (double**)malloc(ts_example1.nouts*sizeof(double*));
	x = (double**)malloc(ts_example1.ninps*sizeof(double*));
	for(i = 0;i < ts_example1.ninps;i++)
		x[i] = (double*)malloc(ns*sizeof(double));

	fh = fopen(argv[2],"r");//Inputs filename
	if(fh == NULL){
		printf("Input file not found!\n");
		exit(1);
	}

	for(i = 0;i < ns;i++){//Reading inputs
		for(j = 0;j < ts_example1.ninps;j++){
			if(j == ts_example1.ninps-1)
				fscanf(fh,"%lf\n",&x[j][i]);
			else
				fscanf(fh,"%lf\t",&x[j][i]);
		}
	}
	fclose(fh);
//------------------------------------

	//printf("Simul cycle, pd.in_vars=%lf\n",pd.in_vars[0]);
//------Simulation cycle--------------
	i=0;
	/*iae = 0.0;
	os = 0.0;
	mp = 0.0;*/
	c=0;
	d=0;

	//printf(" x1\t\t x2\t\t y\n-------------------------------------------\n");
	while(i < ns){
		clearTSController(&ts_example1);

		ts_example1.in_vars[0] = x[0][i];
		ts_example1.in_vars[1] = x[1][i];
		//printf("Input variables--------\nx1 %lf\tx2 %lf\n",ts_example1.in_vars[0],ts_example1.in_vars[1]);

		fuzzifyTS(&ts_example1);
		/*printf("Fuzzification-------\n");
		for(c=0;c<ts_example1.ninps;c++)
			for(d=0;d<ts_example1.udivs[c];d++)
				printf((d==ts_example1.udivs[c]-1)?"%lf\n":"%lf ",ts_example1.fuzzy_in[c][d]);*/

		infereTS(&ts_example1);
		/*printf("Inference-----------\n");
		for(c=0;c<ts_example1.nouts;c++)
			for(d=0;d<ts_example1.n_rules;d++)
				printf((d==ts_example1.n_rules-1)?"%lf\n":"%lf ",ts_example1.fuzzy_out[c][d]);*/

		defuzzifyTS(&ts_example1);

		//printf("%lf\t%lf\t",ts_example1.in_vars[0],ts_example1.in_vars[1]);
		//printf("%lf\n",ts_example1.out_vars[0]);
		printf(((i+1)%11==0)?"%lf\n":"%lf\t",ts_example1.out_vars[0]);

		//iae += (i==0) ? 0.0 : (fabs(e[i]) + fabs(e[i-1])) * proc.tfunc.ts / 2.0;
		//os += (proc.out[i] < ysp) ? 0.0 : (fabs(e[i]) + fabs(e[i-1])) * proc.tfunc.ts * 1.0 / 2.0;
		//mp += ((e[i] - e[i-1]) < 0.0) ? 0.0 : (fabs(e[i]) + fabs(e[i-1])) * proc.tfunc.ts * 10.0 / 2.0;
		i++;
	}
	//printf("%e\n",iae+os+mp);

	return 0;
}
