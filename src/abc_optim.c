//ABC optimization for fuzzy TS system
/*
Finds the controller best suited for a given ideal output
Usage:
abc_optim optim_cfg_filename input_system_filename reference_filename
*/
/* Incluir las bibliotecas ABC */

#include "ccABC.h"
#include <stdlib.h>
#include <string.h>

/* Parámetros del algoritmo de colonia de abejas artificiales */

//#define SN  (10)					// Número de posiciones de alimento
//#define MCN (100)				// Número de iteraciones máximo
//#define DIM	(4)					// Dimensionalidad del problema

/* Definión del problema de optimización */

#define PROBLEM (0)				// Tipo de problema (minimización 0, maximización 1)
#define CONSTRAINT (1)			// Variables a optimizar con 1, o sin 0 restricción

//String constants
#define FIT_APP ("~/Documents/MIE/MIE-IAR/TSinference/bin/tscontrol_fitness")
#define TEMP_SYS ("/Users/sanchez/Documents/MIE/MIE-IAR/TSinference/data/TSconUT.tmp")
#define RES_FILE ("/Users/sanchez/Documents/MIE/MIE-IAR/TSinference/data/fitness_out.tmp")
#define MAX 500//Max string length

//Global variables
int SN;
int MCN;
int DIM;
char input_file[MAX];
char ideal_file[MAX];

double TakagiSugenoProblem(double *x, int n){
	FILE* ts_sys;
	FILE* res_fh;
	char com[500];
	double fit=1000.0;

	if(n!=DIM){
		printf("Not the needed input parameters.\n");
		return fit;
	}

	ts_sys=fopen(TEMP_SYS,"w");
	if(ts_sys == NULL){
		printf("Could not create temporary file.\n");
		return fit;
	}

	fprintf(ts_sys,"2 1\n2 2\n0.0 0.0 0.0 %f\n%f 20.0 20.0 20.0\n0.0 0.0 0.0 %f\n%f 10.0 10.0 10.0\n0 1 2 3\n%f %f %f\n%f %f %f\n%f %f %f\n%f %f %f",x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],x[12],x[13],x[14],x[15]);

	fclose(ts_sys);

	sprintf(com,"%s %s %s %s",FIT_APP,TEMP_SYS,input_file,ideal_file);
	//printf("%s\n",com);
	system(com);

	res_fh=fopen(RES_FILE,"r");
	if(res_fh == NULL){
		printf("Could not create temporary file.\n");
		return fit;
	}

	fscanf(res_fh,"%lf",&fit);
	printf("Read: %lf\n",fit);

	fclose(res_fh);

	return fit;
}

int main(int argc, char** argv) {
	/*Declaración de varibles */
	FILE* conf;
	abc A;
	int i;
	float lim;

	//Arguments:
	//0 - program name
	//1 - optimization configuration filename
	//2 - input system filename
	//3 - reference output filename
	if(argc < 4){
		printf("Not enough input arguments.\n");
		exit(1);
	}

	conf=fopen(argv[1],"r");
	if(conf==NULL){
		printf("Could not open configuration file.\n");
		exit(1);
	}

	//Reading configuration parameters
	fscanf(conf,"%s\n%s\n",input_file,ideal_file);
	fscanf(conf,"%d %d\n%d\n",&SN,&MCN,&DIM);

	srand(time(NULL));		/* Inicializar la semilla de números aleatorios */
	srand48(time(NULL));

	A = abcNew();			/* Inicializar los datos del sitema ABC */
	abcSet(&A,SN,MCN);		/* Enviar los parámetros al sistema ABC */
	abcSetProblem(&A,TakagiSugenoProblem,DIM,PROBLEM,CONSTRAINT);	/* Enviar un problema de optimización */
	abcComputeLimit(&A);		/* Calcular el parámetro limit del sistem ABC */
	abcAllocateMemory(&A);		/* Asignar la memoria al sistema ABC */
	//abcDisplaySettings(&A);	/* Imprimir los parámetros del sistema ABC*/

	//Dimensions:
	//0 - x1, left
	//1 - x1, right
	//2 - x2, left
	//3 - x2, right
	//4 - R1, x1
	//5 - R1, x2
	//6 - R1, c
	//7 - R2, x1
	//8 - R2, x2
	//9 - R2, c
	//10 - R3, x1
	//11 - R3, x2
	//12 - R3, c
	//13 - R4, x1
	//14 - R4, x2
	//15 - R4, c
	for(i = 0; i < A.D; i++){
		fscanf(conf,"%f\n",&lim);
		if(i < 4)//Antecedent parameters
			abcAddLimitAtDimension(&A,i,0.0,lim);	/* Agregar los límites a las variables a optimizar */ 
		else//Consequent parameters
			abcAddLimitAtDimension(&A,i,0.0,lim);					/* Agregar los límites a las variables a optimizar */ 
	}

	abcInit(&A);				/* Inicializar la posiciones de alimento */

	for(i = 0;i < 1;i++){//Solutions
		while((A.x[i][0] < A.x[i][1]) || (A.x[i][2] < A.x[i][3])){//Parameters/input
			abcInit(&A);
		}
		//A.[i][0] remains the same
		A.xmax[1] = A.x[i][0];
		//A.[i][2] remains the same
		A.xmax[3] = A.x[i][2];
	}

	abcMemoryzeBestSolution(&A);		/* Memorizar la mejor solución alcanzada hasta el momento */
	//abcDisplay(&A);			/* Imprimir los datos de las posiciones de alimento */

	//printf("The program reaches here...\n");
	while(!abcMaxCycle(&A)){		/* Hasta que el número máximo de ciclos sea alcanzado */
		abcEmployedBeesPhase(&A);	/* Fase de las abejas obreras */
		abcOnlookerBeesPhase(&A);	/* Fase de las abejas observadoras */
		abcScoutBeesPhase(&A);		/* Fase de las abejas exploradoras */
		abcMemoryzeBestSolution(&A);	/* Memorizar la mejor solución alcanzada hasta el momento */
    	}
	printf("Best solution found:\n");
	abcDisplay(&A);			/* Imprimir los datos de las posiciones de alimento */

	//Showing the best fitness found
	printf("Best fitness found:\n");
	printf("%.10e\t%.10e\n",A.bestF,A.bestFit);

	/*sprintf(str,"bestp_%s.txt",argv[2]);
	best=fopen(str,"w");// Writing the best parameters aka Original
	if(best==NULL){
		printf("Could not open output file.\nBest: ");
		for(i=0;i<A.D;i++)
			printf("x[%d]=%lf%s",i,A.bestX[i],(i==A.D-1)?"\n":", ");
	}
	else{
		for(i=0;i<A.D;i++)
			fprintf(best,"%.10e%s",A.bestX[i],(i==A.D-1)?"":"\n");
	}*/

	abcDelete(&A);				/* Librerar memoria del sistema ABC */

	return 0;
}
