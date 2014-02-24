/*!
    \file ccABC.h
    \brief Archivo de documentación.
    Librería con las herramientas necesarias para la implementación de un sistema de colonia de abejas artificiales.
    \author Marco Antonio Contreras Cruz
    \date 22-12-2012
    \version 1.0
    \copyright Universidad de Guanajuato, DICIS, LaViRIA.
*/

#ifndef CC_ABC_H
#define CC_ABC_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/*!
    \struct abc
    \brief Estructura para la implementación de un algoritmo de colonia de abejas artificiales.
*/
typedef struct abcTag{
    int SN;          					//!< Número de soluciones.
    int MCN;         					//!< Número Máximo de ciclos.
    int D;           					//!< Número de variables a optimizar.
    int limit;       					//!< Número de ciclos antes de que una solución se considere como abandona.
    int cycle;       					//!< Ciclo actual del algoritmo.
    int *C;           					//!< Contador de prueba para cada una de las posiciones de alimento.
	int problem;	 					//!< Tipo de problema 0/1 minimización/maximización.
    int isSet;       					//!< El sistema fue configurado 1/0 si/no.
	int isSetTheProblem; 				//!< Se configuro el problema de optimización.
	int isAllocateTheMemory;			//!< Se ha reservado la memoria 1/0 si/no.
    double **x;      					//!< Posiciones de las fuentes de alimento.
    double *xmin;    					//!< Límites inferiores de las variables a optimizar.
    double *xmax;    					//!< Límites superiores de las variables a optimizar.
    double *f;     						//!< Valor de la posiciones en la función a optimizar.
    double *fit;    					//!< Calidad de las posiciones de alimento.
    double *p;     						//!< Probabilidad de las posiciones de alimento de ser seleccionadas por las abejas observadoras.
	double *bestX;       				//!< Posición de la mejor solucón.
	double bestF;						//!< Evaluación de la mejor fuente.
	double bestFit;						//!< Mejor calidad de la posición de alimento.
    double min;      					//!< Mímimo valor en la función a optimizar.
    double max;      					//!< Máximo valor la función a optimizar.
    double avg;      					//!< Promedio de los valores de la función a optimizar.
    double std;      					//!< Desviación estándar de los valores de la función a optimizar.
    double sumF;   						//!< Suma de lo valores de la función a optimizar.
	double (*function)(double*,int);	//!< Puntero a un función a ser optimizada.
	int constraint;						//!< Con o sin restricción en las variables
}abc;

/*! \fn double abcFitMinimize(double f)
    \brief Calcular la calidad de la posición de alimento en un problema de minimización.
    \param f Evaluación de la posición de alimento en la función a optimizar.
*/
double abcFitMinimize(double f);

/*! \fn double abcFitMaximize(double f)
    \brief Calcular la calidad de la posición de alimento en un problema de maximización.
    \param f Evaluación de la posición de alimento en la función a optimizar.
*/
double abcFitMaximize(double f);

/*! \fn double abcAvg(abc *ABC)
    \brief Calcular el promedio de valores de las posiciones de alimento evaluadas en la función a ser optimizada.
    \param ABC Una colonia de abejas artificiales enviada por referencia.
*/
double abcAvg(abc *ABC);

/*! \fn double abcMin(abc *ABC)
    \brief Calcular el valor mínimo de las posiciones de alimento evaluadas en la función a ser optimizada.
    \param ABC Una colonia de abejas artificiales enviada por referencia.
*/
double abcMin(abc *ABC);

/*! \fn double abcMax(abc *ABC)
    \brief Calcular el valor máximo de las posiciones de alimento evaluadas en la función a ser optimizada.
    \param ABC Una colonia de abejas artificiales enviada por referencia.
*/
double abcMax(abc *ABC);

/*! \fn double abcSumfunction(abc *ABC)
    \brief Calcular la suma de valores de las posiciones de alimento evaluadas en la función a ser optimizada.
    \param ABC Una colonia de abejas artificiales enviada por referencia.
*/
double abcSumfunction(abc *ABC);

/*! \fn double abcStd(abc *ABC)
    \brief Calcular la desviación estándar de los valores de las posiciones de alimento evaluadas en la función a ser optimizada.
    \param ABC Una colonia de abejas artificiales enviada por referencia.
*/
double abcStd(abc *ABC);

/*! \fn int abcSetProblem(abc *ABC,double (*function)(double*,int),int D,int problem)
    \brief Enviar un problema de optimización a la colonia de abejas artificiales.
    \param ABC Una colonia de abejas artificiales enviada por referencia.
    \param function Puntero a la función a optimizar.
    \param D Número de variables a optimizar.
    \param problem Tipo de problema: uno para problemas de maximización y cero para problemas de minimización.
	\param constraint Variables a optimizar con o sin restricción.
*/
int abcSetProblem(abc *ABC,double (*function)(double*,int),int D,int problem,int constraint);

/*! \fn int abcAddLimitAtDimension(abc *ABC,int d,double min,double max)
    \brief Agregar los límites a un variable a optimizar.
    \param ABC Una colonia de abejas artificiales enviada por referencia.
    \param d Número de variable.
    \param min Límite inferior.
    \param max Límite superior.
*/
int abcAddLimitAtDimension(abc *ABC,int d,double min,double max);

/*! \fn int abcSet(abc *ABC,int SN,int MCN)
    \brief Configurar el sistema de colonia de abejas artificiales.
    \param ABC Una colonia de abejas artificiales enviada por referencia.
    \param SN Número de posiciones de alimento.
    \param min Número máximo de ciclos dentro del algoritmo.
*/
int abcSet(abc *ABC,int SN,int MCN);

/*! \fn int abcAddDimension(abc *ABC,int D)
    \brief Agregar el número de variables a optmizar.
    \param ABC Una colonia de abejas artificiales enviada por referencia.
    \param D Número de variables.
*/
int abcAddDimension(abc *ABC,int D);

/*! \fn int abcComputeLimit(abc *ABC)
    \brief Calcular el número de ciclos antes de que una fuente se considere como abandonada.
    \param ABC Una colonia de abejas artificiales enviada por referencia.
    \param D Número de variables.
*/
int abcComputeLimit(abc *ABC);

/*! \fn int abcSetLimit(abc *ABC,int limit)
    \brief Enviar el número de ciclos antes de que una fuente se considere como abandonada.
    \param ABC Una colonia de abejas artificiales enviada por referencia.
    \param limit Número de ciclos.
*/
int abcSetLimit(abc *ABC,int limit);

/*! \fn int abcRoulette(abc *ABC)
    \brief Regresa el índice de una posición de alimento seleccianada por el método de la ruleta.
    \param ABC Una colonia de abejas artificiales enviada por referencia.
    \return Índice de una posición de alimento.
*/
int abcRoulette(abc *ABC);

/*! \fn int abcMaxCycle(abc *ABC)
    \brief El algoritmo ha llegado al número máximo de ciclos.
    \param ABC Una colonia de abejas artificiales enviada por referencia.
    \return uno si se ha llegado al número máximo de ciclos y cero en caso contrario.
*/
int abcMaxCycle(abc *ABC);

/*! \fn int abcAllocateMemory(abc *ABC)
    \brief Asignar la memoria al sistema de abejas artificiales.
    \param ABC Una colonia de abejas artificiales enviada por referencia.
*/
int abcAllocateMemory(abc *ABC);

/*! \fn int abcDelete(abc *ABC)
    \brief Liberar la memoria reservada en el algoritmo.
    \param ABC Una colonia de abejas artificiales enviada por referencia.
*/
int abcDelete(abc *ABC);

/*! \fn int abcDisplay(abc *ABC)
    \brief Imprimir la información de las posiciones de alimento.
    \param ABC Una colonia de abejas artificiales enviada por referencia.
*/
int abcDisplay(abc *ABC);

/*! \fn int abcDisplay(abc *ABC)
    \brief Imprimir la información de la posición de alimento index.
    \param ABC Una colonia de abejas artificiales enviada por referencia.
*/
int abcDisplayFoodPosition(abc *ABC,int index);

/*! \fn int abcLoadSettings(abc *ABC,const char *name)
    \brief Cargar los datos de configuración del sistema.
    \param ABC Una colonia de abejas artificiales enviada por referencia.
    \param name Nombre del archivo de configuración.
*/
int abcLoadSettings(abc *ABC,const char *name);

/*! \fn int abcDisplaySettings(abc *ABC)
    \brief Imprimir los datos de configuración del sistema.
    \param ABC Una colonia de abejas artificiales enviada por referencia.
*/
int abcDisplaySettings(abc *ABC);

/*! \fn abc abcNew(void)
    \brief Crear una nueva variable del tipo abc.
    \return Variable del tipo abc.
*/
abc abcNew(void);

/*! \fn abc int abcInit(abc *ABC)
    \brief Inicializar las posiciones de alimento de forma aleatoria.
    \param ABC Una colonia de abejas artificiales enviada por referencia.
*/
int abcInit(abc *ABC);

/*! \fn int abcEmployedBeesPhase(abc *ABC)
    \brief Implementar la fase de las abejas obreras.
    \param ABC Una colonia de abejas artificiales enviada por referencia.
*/
int abcEmployedBeesPhase(abc *ABC);

/*! \fn int abcOnlookerBeesPhase(abc *ABC)
    \brief Implementar la fase de las abejas observadoras.
    \param ABC Una colonia de abejas artificiales enviada por referencia.
*/
int abcOnlookerBeesPhase(abc *ABC);

/*! \fn int abcScoutBeesPhase(abc *ABC)
    \brief Implementar la fase de las abejas exploradoras.
    \param ABC Una colonia de abejas artificiales enviada por referencia.
*/
int abcScoutBeesPhase(abc *ABC);

/*! \fn int abcMemoryzeBestSolution(abc *ABC)
    \brief Memorizar la mejor posición de alimento encontrada hasta el momento.
    \param ABC Una colonia de abejas artificiales enviada por referencia.
*/
int abcMemoryzeBestSolution(abc *ABC);

/*! \fn void abcError(const char *error)
    \brief Enviar un mensaje de error al usuario.
    \param error Cadena con el mensaje de error.
*/
void abcError(const char *error);

#endif
