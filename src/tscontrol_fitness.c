//Optimization of a TS controller
/*
Returns the fitness value of a fuzzy controller given an ideal output
Usage:
tscontrol_fitness ts_system_filename inputs_filename idealout_filename
*/
//#include "fuzzycontroller.h"
#include <stdio.h>
#include <stdlib.h> 
#include <string.h>
#include <math.h>

#define MAX 500
#define OUT_FILE ("/Users/sanchez/Documents/MIE/MIE-IAR/TSinference/data/tsoutput.dat")
#define CTRL_COM ("/Users/sanchez/Documents/MIE/MIE-IAR/TSinference/bin/tscontroller")
#define RES_FILE ("/Users/sanchez/Documents/MIE/MIE-IAR/TSinference/data/fitness_out.tmp")

int main(int argc, char** argv){
	FILE* fh;
	double accum;
	double* out;
	double* z;
	int i;
	//int j;
	int ns = 0;
	char command[MAX];

	if(argc < 4){
		printf("Not enough input arguments\n");
		exit(1);
	}

//------Ideal output Initialization---
	fh = fopen(argv[3],"r");

	if(fh == NULL){
		printf("Ideal output file not found\n");
		exit(1);
	}

	z = (double*)malloc(ns*sizeof(double));
	out = (double*)malloc(ns*sizeof(double));

	/*for(i = 0;i < ns;i++){//Reading inputs
		fscanf(fh,"%lf\n",&z[i]);
	}*/
	while(!feof(fh)){
		fscanf(fh,"%lf\n",&z[i]);
		ns++;
	}

	printf("ns is %d\n",ns);

	fclose(fh);
//------------------------------------

//------Error calculation cycle-------
	//Use the controller and input files to call the tscontroller program and produce the outputs
	sprintf(command,"%s %s %s > %s",CTRL_COM,argv[1],argv[2],OUT_FILE);
	system(command);

	//printf("Running: %s\n",command);

	//Calculate the error of the output against the ideal output
	fh = fopen(OUT_FILE,"r");
	if(fh == NULL){
		printf("Testing output file not found\n");
		exit(1);
	}

	i=0;
	accum=0.0;
	while(i < ns){
		fscanf(fh,"%lf\n",&out[i]);
		//printf("%lf\n",out[i]);
		accum += pow(z[i] - out[i],2.0f);
		i++;
	}

	fclose(fh);

	fh = fopen(RES_FILE,"w");
	if(fh == NULL){
		printf("Fitness output file not created\n");
		exit(1);
	}

	fprintf(fh,"%e\n",sqrt(accum/ns));
	printf("Wrote: %e\n",sqrt(accum/ns));

	fclose(fh);

	return 0;
}
