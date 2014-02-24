#include "system.h"
#include <math.h>

void cleanCSystem(ctf* sys){
	int i;

	//while(fabs(sys->num.coef[sys->num.deg])<1.0e-12 && fabs(sys->den.coef[sys->den.deg])<1.0e-12){
	while(sys->num.coef[sys->num.deg]==0.0 && sys->den.coef[sys->den.deg]==0.0){
		sys->num.deg--;
		sys->den.deg--;
	}

	//while(fabs(sys->num.coef[0])<1.0e-12){
	while(sys->num.coef[0]==0.0){
		for(i=0;i<sys->num.deg-1;i++){
			sys->num.coef[i]=sys->num.coef[i+1];
		}
		sys->num.deg--;
	}

	//while(fabs(sys->den.coef[0])<1.0e-12){
	while(sys->den.coef[0]==0.0){
		for(i=0;i<sys->den.deg-1;i++){
			sys->den.coef[i]=sys->den.coef[i+1];
		}
		sys->den.deg--;
	}
}

void readTfData(FILE* fh, poly* N, poly* D){
	int i;

	fscanf(fh,"%d",&(N->deg));
	N->coef=(double*)malloc((N->deg+1)*sizeof(double));
	for(i=0;i < N->deg+1;i++){
		fscanf(fh,"%lf",&(N->coef[i]));
	}

	fscanf(fh,"%d",&(D->deg));
	D->coef=(double*)malloc((D->deg+1)*sizeof(double));
	for(i=0;i < D->deg+1;i++){
		fscanf(fh,"%lf",&(D->coef[i]));
	}
}

void readTs(FILE* fh, dtf* discrete){
	fscanf(fh,"%lf\n",&(discrete->ts));
}

void initBilinear(poly* num, poly* den, double ts){
	num->deg=1;
	num->coef=(double*)malloc(2*sizeof(double));
	num->coef[0]=2.0/ts;
	num->coef[1]=-2.0/ts;
	den->deg=1;
	den->coef=(double*)malloc(2*sizeof(double));
	den->coef[0]=1.0;
	den->coef[1]=1.0;
}

void subztitution(ctf* con, dtf* dis){
	int i;
	int p;
	poly* bi_n;
	poly* bi_d;
	poly* vecn[con->num.deg+1];
	poly* vecd[con->den.deg+1];
	poly* aux;

	bi_n =  (poly*)malloc(sizeof(poly));
	bi_d =  (poly*)malloc(sizeof(poly));
	initBilinear(bi_n,bi_d,dis->ts);

	if(con->num.deg < con->den.deg)
		p = con->den.deg;
	else
		p = con->num.deg;

	for(i=0;i < con->num.deg+1;i++){//Multiplication by the numerator of the bilinear function (numerator)
		if(i == con->num.deg)
			vecn[i] = coef2poly(con->num.coef[i]);
		else{
			aux = coef2poly(con->num.coef[i]);
			vecn[i] = pmult(aux, ppow(bi_n, con->num.deg - i));
		}
	}

	for(i=0;i < con->den.deg+1;i++){//(denominator)
		if(i == con->den.deg)
			vecd[i] = coef2poly(con->den.coef[i]);
		else{
			aux = coef2poly(con->den.coef[i]);
			vecd[i] = pmult(aux, ppow(bi_n, con->den.deg - i));
		}
	}

	for(i=0;i < con->num.deg+1;i++){//Multiplication by the denominator of the bilinear function
		vecn[i] = pmult(vecn[i], ppow(bi_d, p - (con->num.deg - i)));
	}

	for(i=0;i < con->den.deg+1;i++)
		vecd[i] = pmult(vecd[i], ppow(bi_d, p - (con->den.deg - i)));

	aux = vecn[0];//Simplification
	for(i=1;i < con->num.deg+1;i++)
		aux = psum(aux, vecn[i]);
	dis->num = *aux;

	aux = vecd[0];
	for(i=1;i < con->den.deg+1;i++)
		aux = psum(aux, vecd[i]);
	dis->den = *aux;
}

void toNegPows(dtf* dis){
	int i;
	int p = dis->num.deg;
	int lim = (p%2) ? (p+1)/2 : p/2;
	double aux;

	for(i=0;i < lim;i++){
		aux = dis->num.coef[i];
		dis->num.coef[i] = dis->num.coef[p - i];
		dis->num.coef[p - i] = aux;

		aux = dis->den.coef[i];
		dis->den.coef[i] = dis->den.coef[p - i];
		dis->den.coef[p - i] = aux;
	}
}

void contSystemMult(ctf* a, ctf* b){
	a->num = *(pmult(&a->num,&b->num));
	a->den = *(pmult(&a->den,&b->den));
}

void discSystemMult(dtf* a, dtf* b){
	a->num = *(pmult(&a->num,&b->num));
	a->den = *(pmult(&a->den,&b->den));
	if(a->ts != b->ts)
		printf("Warning: DSystem multiplicands have different ts\n");
}

void unityNegFeedback(ctf* sys){
	sys->den = *(psum(&sys->num,&sys->den));
}

void contZOH(poly* num, poly* den){
	num->deg=0;
	num->coef=(double*)malloc(sizeof(double));
	num->coef[0]=1.0;
	den->deg=1;
	den->coef=(double*)malloc(2*sizeof(double));
	den->coef[0]=1.0;
	den->coef[1]=0.0;
}

void discZOH(poly* num, poly* den){
	num->deg=1;
	num->coef=(double*)malloc(2*sizeof(double));
	num->coef[0]=1.0;
	num->coef[1]=-1.0;
	den->deg=1;
	den->coef=(double*)malloc(2*sizeof(double));
	den->coef[0]=1.0;
	den->coef[1]=0.0;
}

void evalSystemStep(dsystem* sys){
	int j;
	int step = sys->dstate;
	sys->out[step] = (sys->tfunc.num.coef[sys->tfunc.num.deg] * sys->in[step]) / sys->tfunc.den.coef[sys->tfunc.den.deg];

	for(j=1;j < sys->tfunc.num.deg+1 && step > j-1;j++){
		sys->out[step] += ((sys->tfunc.num.coef[sys->tfunc.num.deg-j] * sys->in[step-j]) - (sys->tfunc.den.coef[sys->tfunc.den.deg-j] * sys->out[step-j])) / sys->tfunc.den.coef[sys->tfunc.den.deg];
	}
	sys->dstate++;
}

void evalSystemResponse(dsystem* sys, unsigned int nlim){
	sys->dstate = 0;
	while(sys->dstate < nlim){
		evalSystemStep(sys);
	}
}

void displayCSystem(ctf* sys){
	int i;
	for (i=0; i<sys->num.deg+1; i++) {
		printf("%e s^%d ",sys->num.coef[i],sys->num.deg-i);
	}
	printf("\n");
	for (i=0; i<sys->den.deg+1; i++) {
		printf("%e s^%d ",sys->den.coef[i],sys->den.deg-i);
	}
	printf("\n");
}

void displayDSystem(dtf* sys){
	int i;
	for (i=0; i<sys->num.deg+1; i++) {
		printf("%e z^%d ",sys->num.coef[i],sys->num.deg-i);
	}
	printf("\n");
	for (i=0; i<sys->den.deg+1; i++) {
		printf("%e z^%d ",sys->den.coef[i],sys->den.deg-i);
	}
	printf("\n");
}

void initializeProcess(dsystem* proc, FILE* fh){
	int j;

	//ts
	fscanf(fh,"%lf\n",&(proc->tfunc.ts));

	//numerator coefs
	fscanf(fh,"%d\n",&(proc->tfunc.num.deg));
	proc->tfunc.num.coef = (double*)malloc((proc->tfunc.num.deg+1)*sizeof(double));
	for(j=0;j < proc->tfunc.num.deg+1;j++)
		fscanf(fh,"%lf\n",&(proc->tfunc.num.coef[j]));

	//denominator coefs
	fscanf(fh,"%d\n",&(proc->tfunc.den.deg));
	proc->tfunc.den.coef = (double*)malloc((proc->tfunc.den.deg+1)*sizeof(double));
	for(j=0;j < proc->tfunc.den.deg+1;j++)
		fscanf(fh,"%lf\n",&(proc->tfunc.den.coef[j]));

	fclose(fh);

	if(proc->tfunc.num.deg!=proc->tfunc.den.deg){
		printf("Error: ndeg and ddeg should match\n");
		exit(1);
	}
}
