//Fuzzy Controller library
#include "fuzzycontroller.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

void fuzzify(int ninps, int* ndivs, double* crisp, memfunc** func_mat, double** fuzzy){
	int i;
	int j;
	double m;
	double b;

	i=0;
	while(i < ninps){
		j=0;
		while(j < ndivs[i]){
			if(crisp[i] > func_mat[i][j].lolim && crisp[i] < func_mat[i][j].uplim){
				if(crisp[i] > (func_mat[i][j].center - (func_mat[i][j].core/2.0)) && crisp[i] < (func_mat[i][j].center + (func_mat[i][j].core/2.0))){
					fuzzy[i][j] = 1.0;
				}
				else{
					if(crisp[i] < func_mat[i][j].center){
						m = 1.0/((func_mat[i][j].center - func_mat[i][j].core / 2.0) - func_mat[i][j].lolim);
						b = -1.0*m*func_mat[i][j].lolim;
					}
					else{
						m = -1.0/(func_mat[i][j].uplim - (func_mat[i][j].center + func_mat[i][j].core / 2.0));
						b = -1.0*m*func_mat[i][j].uplim;
					}
					fuzzy[i][j] = m * crisp[i] + b;
				}
			}
			else{//From here
				if(crisp[i] <= func_mat[i][j].lolim && func_mat[i][j].lolim == (func_mat[i][j].center - func_mat[i][j].core / 2.0) && j==0){
					fuzzy[i][j] = 1.0;
				}
				else if(crisp[i] >= func_mat[i][j].uplim && func_mat[i][j].uplim == (func_mat[i][j].center + func_mat[i][j].core / 2.0) && j==ndivs[i]-1){
					fuzzy[i][j] = 1.0;
				}
				else{//to here
					fuzzy[i][j] = 0.0;
				}
			}
			j++;
		}
		i++;
	}
}

void fuzzifyTS(fuzzyTScontroller* ts_fuzzy){
	int i;
	int j;
	double m;
	double b;

	i=0;
	while(i < ts_fuzzy->ninps){
		j=0;
		while(j < ts_fuzzy->udivs[i]){
			if(ts_fuzzy->in_vars[i] > ts_fuzzy->func_mat[i][j].downleft && ts_fuzzy->in_vars[i] < ts_fuzzy->func_mat[i][j].downright){
				if(ts_fuzzy->in_vars[i] >= ts_fuzzy->func_mat[i][j].upleft && ts_fuzzy->in_vars[i] <= ts_fuzzy->func_mat[i][j].upright){
					ts_fuzzy->fuzzy_in[i][j] = 1.0;
				}
				else{
					if(ts_fuzzy->in_vars[i] < ts_fuzzy->func_mat[i][j].upleft){
						m = 1.0/(ts_fuzzy->func_mat[i][j].upleft - ts_fuzzy->func_mat[i][j].downleft);
						b = -1.0 * m * ts_fuzzy->func_mat[i][j].downleft;
					}
					else{
						m = -1.0/(ts_fuzzy->func_mat[i][j].downright - ts_fuzzy->func_mat[i][j].upright);
						b = -1.0 * m * ts_fuzzy->func_mat[i][j].downright;
					}
					ts_fuzzy->fuzzy_in[i][j] = m * ts_fuzzy->in_vars[i] + b;
				}
			}
			else{//From here
				if(ts_fuzzy->in_vars[i] <= ts_fuzzy->func_mat[i][j].downleft && j==0){
					ts_fuzzy->fuzzy_in[i][j] = 1.0;
				}
				else if(ts_fuzzy->in_vars[i] >= ts_fuzzy->func_mat[i][j].downright && j==ts_fuzzy->udivs[i]-1){
					ts_fuzzy->fuzzy_in[i][j] = 1.0;
				}
				else{//to here
					ts_fuzzy->fuzzy_in[i][j] = 0.0;
				}
			}
			j++;
		}
		i++;
	}
}

// 0 - max,min union
// 1 - max,max union
void infere2_1(int* udivs, double** ifuzzy, famtable* fam, double** ofuzzy){
	int i;
	int j;
	int offset=0;

	i=0;
	while(i < udivs[0]){
		j=0;
		while(j < udivs[1]){
			if(fam->out_funcs[offset+j] > -1){
				switch(fam->utype[offset+j]){
					case 0:
						if(ofuzzy[0][fam->out_funcs[offset+j]] < min(ifuzzy[0][i],ifuzzy[1][j]))
							ofuzzy[0][fam->out_funcs[offset+j]] = min(ifuzzy[0][i],ifuzzy[1][j]);
						break;
					case 1:
						if(ofuzzy[0][fam->out_funcs[offset+j]] < max(ifuzzy[0][i],ifuzzy[1][j]))
							ofuzzy[0][fam->out_funcs[offset+j]] = max(ifuzzy[0][i],ifuzzy[1][j]);
						break;
					default:
						printf("Error: Not a valid union type entry %d %d\n",i,j);
						break;
				}
			}
			j++;
		}
		offset += udivs[1];
		i++;
	}
}

void infereTS(fuzzyTScontroller* ts_fuzzy){
	int i;
	int j;
	int k;
	int l;

	for(i = 0;i < ts_fuzzy->ninps-1;i++){
		for(j = 0;j < ts_fuzzy->udivs[i];j++){
			for(k = i+1;k < ts_fuzzy->ninps;k++){
				for(l = 0;l < ts_fuzzy->udivs[k];l++){//0 means it just works for 1 output
					ts_fuzzy->fuzzy_out[0][ts_fuzzy->fam->out_funcs[(i*(ts_fuzzy->n_rules/ts_fuzzy->udivs[i]))+(j*ts_fuzzy->udivs[k])+l]] = evalMultiVariable(&ts_fuzzy->rules[ts_fuzzy->fam->out_funcs[(i*(ts_fuzzy->n_rules/ts_fuzzy->udivs[i]))+(j*ts_fuzzy->udivs[k])+l]],ts_fuzzy->in_vars);
//					printf("Evaluating rule %d: %lf\n",ts_fuzzy->fam->out_funcs[(i*(ts_fuzzy->n_rules/ts_fuzzy->udivs[i]))+(j*ts_fuzzy->udivs[k])+l],ts_fuzzy->fuzzy_out[0][ts_fuzzy->fam->out_funcs[(i*(ts_fuzzy->n_rules/ts_fuzzy->udivs[i]))+(j*ts_fuzzy->udivs[k])+l]]);
				}
			}
		}
	}
}

double max(double a, double b){
	return (a > b) ? a : b;
}

double min(double a, double b){
	return (a > b) ? b : a;
}

polygon* createUnionPolygon(int nouts, int* udivs, memfunc** func_mat, double** fuzzy, char* inference){
	int first;
	int i;
	int j;
	int ninps=2;
	polygon* total;
	polygon* inc;

	if(!strcmp(inference,"mamdani")){
		i=0;
		while(i < nouts){
			j=0;
			first=1;
			while(j < udivs[ninps+i]){
				if(fuzzy[i][j] > 0.0){
					if(fuzzy[i][j+1] > 0.0){
						inc = mamdaniPolygon(&func_mat[i+ninps][j+1],fuzzy[i][j+1]);
						if(first == 1){
							total = mamdaniPolygon(&func_mat[i+ninps][j],fuzzy[i][j]);
							first = 0;
						}
						total = mamdaniUnion(total,inc);
					}
					else if(first == 1){
						total = mamdaniPolygon(&func_mat[i+ninps][j],fuzzy[i][j]);
						first = 0;
					}
				}
				j++;
			}
			i++;
		}
	}
	else if(!strcmp(inference,"smamdani")){
		//printf("Warning: Creating a polygon from singletons\n");
		i=0;
		while(i < nouts){
			j=0;
			first=1;
			while(j < udivs[ninps+i]){
				if(fuzzy[i][j] > 0.0){
					if(first == 1){
						total = initializePolygon(3);
						total->varr[0].x = func_mat[ninps+i][j].center;
						total->varr[0].y = 0.0;
						total->nvert++;
						total->varr[1].x = func_mat[ninps+i][j].center;
						total->varr[1].y = fuzzy[i][j];
						total->nvert++;
						total->varr[2].x = func_mat[ninps+i][j].center;
						total->varr[2].y = 0.0;
						total->nvert++;
						first = 0;
					}
					else{
						int k=0;
						inc = initializePolygon(total->nvert+1);
						while(k < total->nvert-1){
							inc->varr[k].x = total->varr[k].x;
							inc->varr[k].y = total->varr[k].y;
							inc->nvert++;
							k++;
						}
						inc->varr[k].x = func_mat[ninps+i][j].center;
						inc->varr[k].y = fuzzy[i][j];
						inc->nvert++;
						inc->varr[k+1].x = func_mat[ninps+i][j].center;
						inc->varr[k+1].y = 0.0;
						inc->nvert++;

						free(total);
						total = inc;
						inc = NULL;
					}
				}
				j++;
			}
			i++;
		}
	}
	else{
		printf("Error: Inference method not recognized\n");
	}

	return total;
}

polygon* mamdaniPolygon(memfunc* out_func, double fuzzy){
	polygon* poly;
	double m;
	double b;

	poly = initializePolygon(4);

	poly->varr[0].x = out_func->lolim;
	poly->varr[0].y = 0.0;
	poly->nvert++;

	m = 1.0 / ((out_func->center - out_func->core / 2.0) - out_func->lolim);
	b = -1.0 * m * out_func->lolim;
	poly->varr[1].x = (fuzzy - b) / m;
	poly->varr[1].y = fuzzy;
	poly->nvert++;

	m = -1.0 / (out_func->uplim - (out_func->center + out_func->core / 2.0));
	b = -1.0 * m * out_func->uplim;
	poly->varr[2].x = (fuzzy - b) / m;
	poly->varr[2].y = fuzzy;
	poly->nvert++;

	poly->varr[3].x = out_func->uplim;
	poly->varr[3].y = 0.0;
	poly->nvert++;

	polygonArea(poly);

	return poly;
}

polygon* initializePolygon(int verts){
	polygon* poly =  (polygon*)malloc(sizeof(polygon));

	poly->nvert = 0;
	poly->varr = (vertex*)malloc(verts*sizeof(vertex));
	poly->area = 0.0;

	return poly;
}

polygon* mamdaniUnion(polygon* a, polygon* b){
	polygon* total;
	vertex* isec;
	double ua;
	double ub;
	double intval=0.5;
	int i;
	int j;
	int v;
	int ulim;
	int llim;

	if(a->nvert < 4 || b->nvert < 4){
		printf("Error: Not Mamdani polygons to merge\n");
		return total;
	}

	ua = a->varr[a->nvert - 2].y;
	ub = b->varr[b->nvert - 2].y;

	if(ua > intval && ub > intval){//Takes mem-func intersection at 0.5
		isec = linesIntersect(&a->varr[a->nvert-2],&a->varr[a->nvert-1],&b->varr[0],&b->varr[1]);

		if(a->varr[a->nvert - 2].y == a->varr[a->nvert - 3].y){
			v = a->nvert + b->nvert - 2;
			ulim = 2;
			llim = 1;
		}
		else if(b->varr[1].y == b->varr[2].y){
			v = a->nvert + b->nvert - 2;
			ulim = 1;
			llim = 2;
		}
		else{
			//isec = linesIntersect(&a->varr[a->nvert-2],&a->varr[a->nvert-1],&b->varr[0],&b->varr[1]);
			v = a->nvert + b->nvert - 1;
			ulim = 1;
			llim = 1;
		}
	}
	else{
		if(ua < ub){
			isec = linesIntersect(&a->varr[a->nvert-3],&a->varr[a->nvert-2],&b->varr[0],&b->varr[1]);
			v = a->nvert + b->nvert - 2;
			ulim = 2;
			llim = 1;
		}
		else if(ua > ub){
			isec = linesIntersect(&a->varr[a->nvert-2],&a->varr[a->nvert-1],&b->varr[1],&b->varr[2]);
			v = a->nvert + b->nvert - 2;
			ulim = 1;
			llim = 2;
		}
		else{
			v = a->nvert + b->nvert - 4;
			ulim = 2;
			llim = 2;
		}
	}

	total = initializePolygon(v);
	i=0;
	while(i < a->nvert-ulim){
		total->varr[i] = a->varr[i];
		total->nvert++;
		i++;
	}

	if(ua != ub || (ua > intval && ub > intval)){
		total->varr[i] = *isec;
		total->nvert++;
		i++;
	}

	j = llim;
	while(j < b->nvert){
		total->varr[i] = b->varr[j];
		total->nvert++;
		i++;
		j++;
	}

	polygonArea(total);

	return total;
}

vertex* linesIntersect(vertex* a1, vertex* a2, vertex* b1, vertex* b2){
	double ma;
	double ba;
	double mb;
	double bb;
	vertex* isec = (vertex*)malloc(sizeof(vertex));

	ma = (a1->y - a2->y) / (a1->x - a2->x);
	ba = a1->y - ma * a1->x;
	mb = (b1->y - b2->y) / (b1->x - b2->x);
	bb = b1->y - mb * b1->x;

	isec->x = (bb - ba) / (ma - mb);
	isec->y = ma * isec->x + ba;

	return isec;
}

void polygonArea(polygon* poly){//The polygon is defined CW, change last sign for CCW
	int i;
	int next;

	poly->area = 0.0;
	i=0;
	while(i < poly->nvert){
		next = (i == poly->nvert-1) ? 0 : i + 1;
		poly->area += (poly->varr[i].x * poly->varr[next].y - poly->varr[next].x * poly->varr[i].y) / 2.0;
		i++;
	}

	poly->area *= -1.0;
}

void defuzzify(char* method, polygon* poly, int nouts, double* out_vars){
	int i;
	int j;
	double r;
	double m;
	double b;
	//double aux;

	if(!strcmp(method,"centroid")){
		i=0;
		while(i < nouts){
			j=0;
			r=0.0;
			while(j < poly[i].nvert-1){
				m = (poly[i].varr[j+1].y - poly[i].varr[j].y) / (poly[i].varr[j+1].x - poly[i].varr[j].x);
				b = poly[i].varr[j].y - m * poly[i].varr[j].x;
				r += (m / 3.0) * (pow(poly[i].varr[j+1].x,3.0) - pow(poly[i].varr[j].x,3.0)) + (b / 2.0) * (pow(poly[i].varr[j+1].x,2.0) - pow(poly[i].varr[j].x,2.0));
				/*aux = (m / 3.0) * (pow(poly[i].varr[j+1].x,3.0) - pow(poly[i].varr[j].x,3.0)) + (b / 2.0) * (pow(poly[i].varr[j+1].x,2.0) - pow(poly[i].varr[j].x,2.0));
				r += aux;
				printf("Defuzzy vtx %d, m: %lf b: %lf aux: %lf\n",j,m,b,aux);*/
				j++;
			}
			out_vars[i] = r / poly[i].area;
			i++;
		}
	}
	else if(!strcmp(method,"sentroid")){
		i=0;
		while(i < nouts){
			j=0;
			r=0.0;
			m=0.0;
			//printf("Nvert %d\n",poly[0].nvert);
			while(j < poly[i].nvert-2){//CHANGE: Include all the membership values in m
				r += (poly[i].varr[j+1].x) * (poly[i].varr[j+1].y);
				m += poly[i].varr[j+1].y;
				j++;
			}
			out_vars[i] = r / m;
			i++;
		}
	}
	else{
		printf("Error: Not a valid defuzzification method\n");
	}
}

void defuzzifyTS(fuzzyTScontroller* ts_fuzzy){
	int i;
	int j;
	int k;
	int l;
	double num_sum;
	double den_sum;
	//double valid_in1;
	//double valid_in2;
	double minwt;
	double* w = (double*)malloc(ts_fuzzy->n_rules*sizeof(double));

	for(i = 0;i < ts_fuzzy->ninps-1;i++){
		for(j = 0;j < ts_fuzzy->udivs[i];j++){
			for(k = i+1;k < ts_fuzzy->ninps;k++){
				for(l = 0;l < ts_fuzzy->udivs[k];l++){
					//printf("%lf > 0.0\t%lf > 0.0\n",(ts_fuzzy->rules[(i*(ts_fuzzy->n_rules/ts_fuzzy->udivs[i]))+(j*ts_fuzzy->udivs[k])+l].coef[i]),(ts_fuzzy->rules[(i*(ts_fuzzy->n_rules/ts_fuzzy->udivs[i]))+(j*ts_fuzzy->udivs[k])+l].coef[k]));
					if((ts_fuzzy->rules[(i*(ts_fuzzy->n_rules/ts_fuzzy->udivs[i]))+(j*ts_fuzzy->udivs[k])+l].coef[i] != 0.0) && (ts_fuzzy->rules[(i*(ts_fuzzy->n_rules/ts_fuzzy->udivs[i]))+(j*ts_fuzzy->udivs[k])+l].coef[k] != 0.0)){
						minwt = min(ts_fuzzy->fuzzy_in[i][j],ts_fuzzy->fuzzy_in[k][l]);
					}
					else{
						if(ts_fuzzy->rules[(i*(ts_fuzzy->n_rules/ts_fuzzy->udivs[i]))+(j*ts_fuzzy->udivs[k])+l].coef[i] > 0.0)
							minwt = ts_fuzzy->fuzzy_in[i][j];
						else if(ts_fuzzy->rules[(i*(ts_fuzzy->n_rules/ts_fuzzy->udivs[i]))+(j*ts_fuzzy->udivs[k])+l].coef[k] > 0.0)
							minwt = ts_fuzzy->fuzzy_in[k][l];
						else
							minwt = 0.0;
					}
					w[(i*(ts_fuzzy->n_rules/ts_fuzzy->udivs[i]))+(j*ts_fuzzy->udivs[k])+l] = minwt;
					//printf("idx: %d\twght: %lf\n",(i*(ts_fuzzy->n_rules/ts_fuzzy->udivs[i]))+(j*ts_fuzzy->udivs[k])+l,w[(i*(ts_fuzzy->n_rules/ts_fuzzy->udivs[i]))+(j*ts_fuzzy->udivs[k])+l]);
				}
			}
		}
	}

	for(i = 0;i < ts_fuzzy->nouts;i++){
		for(j = 0,num_sum = 0.0,den_sum = 0.0;j < ts_fuzzy->n_rules;j++){
			num_sum += w[j]*ts_fuzzy->fuzzy_out[i][j];
			den_sum += w[j];
		}
		ts_fuzzy->out_vars[i] = (den_sum < ZERO_TH) ? 0.0 : num_sum/den_sum;
	}
}

void clearController(fuzzycontroller* fc){
	int i;
	int j;

	i=0;
	while(i < (fc->ninps + fc->nouts)){
		j=0;
		while(j < fc->udivs[i]){
			if(i < fc->ninps)
				fc->fuzzy_in[i][j] = 0.0;
			else
				fc->fuzzy_out[i - fc->ninps][j] = 0.0;
			j++;
		}
		i++;
	}
}

void clearTSController(fuzzyTScontroller* fc){
	int i;
	int j;

	i=0;
	while(i < fc->ninps){
		j=0;
		while(j < fc->udivs[i]){
			fc->fuzzy_in[i][j] = 0.0;
			j++;
		}
		i++;
	}

	i=0;
	while(i < fc->nouts){
		j=0;
		while(j < fc->n_rules){
			fc->fuzzy_out[i][j] = 0.0;
			j++;
		}
		i++;
	}
}

void initializeController(fuzzycontroller* fuzzyc, FILE* fh){
	int i;
	int j;

	fscanf(fh,"%d ",&fuzzyc->ninps);
	fscanf(fh,"%d\n",&fuzzyc->nouts);
	fuzzyc->in_vars = (double*)malloc(fuzzyc->ninps*sizeof(double));
	fuzzyc->out_vars = (double*)malloc(fuzzyc->nouts*sizeof(double));
	fuzzyc->udivs = (int*)malloc((fuzzyc->ninps+fuzzyc->nouts)*sizeof(int));

	//Assign udivs
	for(i=0;i < (fuzzyc->ninps+fuzzyc->nouts);i++){
		if(i < fuzzyc->ninps)
			fscanf(fh,"%d ",&fuzzyc->udivs[i]);//inputs
		else
			fscanf(fh,"%d\n",&fuzzyc->udivs[i]);//outputs (expecting ONE output '\n')
	}

	fuzzyc->func_mat = (memfunc**)malloc((fuzzyc->ninps+fuzzyc->nouts)*sizeof(memfunc));
	for(i=0;i<(fuzzyc->ninps+fuzzyc->nouts);i++)
		fuzzyc->func_mat[i] = (memfunc*)malloc(fuzzyc->udivs[i]*sizeof(memfunc));

	fuzzyc->fuzzy_in = (double**)malloc(fuzzyc->ninps*sizeof(double));
	for(i=0;i<fuzzyc->ninps;i++)
		fuzzyc->fuzzy_in[i] = (double*)malloc(fuzzyc->udivs[i]*sizeof(double));

	fuzzyc->fuzzy_out = (double**)malloc(fuzzyc->nouts*sizeof(double));
	for(i=0;i<fuzzyc->nouts;i++)
		fuzzyc->fuzzy_out[i] = (double*)malloc(fuzzyc->udivs[i+fuzzyc->ninps]*sizeof(double));

	//Assign membership functions
	for(i=0;i < (fuzzyc->ninps + fuzzyc->nouts);i++){
		if(i < fuzzyc->ninps)
			for(j=0;j < fuzzyc->udivs[i];j++){
				fscanf(fh,"%lf ",&fuzzyc->func_mat[i][j].center);
				fscanf(fh,"%lf ",&fuzzyc->func_mat[i][j].core);
				fscanf(fh,"%lf ",&fuzzyc->func_mat[i][j].uplim);
				fscanf(fh,"%lf\n",&fuzzyc->func_mat[i][j].lolim);
			}
		else
			for(j=0;j < fuzzyc->udivs[i];j++){
				fscanf(fh,"%lf ",&fuzzyc->func_mat[i][j].center);
				fscanf(fh,"%lf ",&fuzzyc->func_mat[i][j].core);
				fscanf(fh,"%lf ",&fuzzyc->func_mat[i][j].uplim);
				fscanf(fh,"%lf\n",&fuzzyc->func_mat[i][j].lolim);
			}
	}

	//Assign FAM (ninps sensitive; 2 inps)
	fuzzyc->fam = (famtable*)malloc(sizeof(famtable));

	fuzzyc->fam->ndims = fuzzyc->ninps;
	fuzzyc->fam->utype = (int*)malloc((fuzzyc->udivs[0]*fuzzyc->udivs[1])*sizeof(int));//dim(in1) x dim(in2)
	fuzzyc->fam->out_funcs = (int*)malloc((fuzzyc->udivs[0]*fuzzyc->udivs[1])*sizeof(int));

	for(i=0;i<(fuzzyc->udivs[0]*fuzzyc->udivs[1]);i++)
		fscanf(fh,"%d",&fuzzyc->fam->utype[i]);

	for(i=0;i<(fuzzyc->udivs[0]*fuzzyc->udivs[1]);i++)
		fscanf(fh,"%d",&fuzzyc->fam->out_funcs[i]);

	fclose(fh);
}

void initializeTSController(fuzzyTScontroller* fuzzyc, FILE* fh){
	int i;
	int j;

	//Assign ninps, nouts
	fscanf(fh,"%d ",&fuzzyc->ninps);
	fscanf(fh,"%d\n",&fuzzyc->nouts);

	//Allocate in_vars, out_vars, udivs
	fuzzyc->in_vars = (double*)malloc(fuzzyc->ninps*sizeof(double));
	fuzzyc->out_vars = (double*)malloc(fuzzyc->nouts*sizeof(double));
	fuzzyc->udivs = (int*)malloc(fuzzyc->ninps*sizeof(int));

	//Assign udivs
	for(i = 0;i < fuzzyc->ninps;i++){
		if(i == fuzzyc->ninps-1)
			fscanf(fh,"%d\n",&fuzzyc->udivs[i]);//inputs
		else
			fscanf(fh,"%d ",&fuzzyc->udivs[i]);//inputs
	}

	//Allocate MF matrix
	fuzzyc->func_mat = (memfunc**)malloc(fuzzyc->ninps*sizeof(memfunc*));
	for(i = 0;i < fuzzyc->ninps;i++)
		fuzzyc->func_mat[i] = (memfunc*)malloc(fuzzyc->udivs[i]*sizeof(memfunc));

	//Number of rules
	fuzzyc->n_rules = 1;
	for(i = 0;i < fuzzyc->ninps;i++)
		fuzzyc->n_rules *= fuzzyc->udivs[i];

	//Allocate rules
	fuzzyc->rules = (poly*)malloc(fuzzyc->n_rules*sizeof(poly));

	//Allocate fuzzy variables
	fuzzyc->fuzzy_in = (double**)malloc(fuzzyc->ninps*sizeof(double));
	for(i = 0;i < fuzzyc->ninps;i++)
		fuzzyc->fuzzy_in[i] = (double*)malloc(fuzzyc->udivs[i]*sizeof(double));

	fuzzyc->fuzzy_out = (double**)malloc(fuzzyc->nouts*sizeof(double));
	for(i = 0;i < fuzzyc->nouts;i++)
		fuzzyc->fuzzy_out[i] = (double*)malloc(fuzzyc->n_rules*sizeof(double));

	//Assign membership functions
	for(i = 0;i < fuzzyc->ninps;i++){
		for(j=0;j < fuzzyc->udivs[i];j++){
			fscanf(fh,"%lf ",&fuzzyc->func_mat[i][j].downleft);
			fscanf(fh,"%lf ",&fuzzyc->func_mat[i][j].upleft);
			fscanf(fh,"%lf ",&fuzzyc->func_mat[i][j].upright);
			fscanf(fh,"%lf\n",&fuzzyc->func_mat[i][j].downright);
		}
	}

	//Allocate & assign FAM (ninps (in)sensitive)
	fuzzyc->fam = (famtable*)malloc(sizeof(famtable));

	fuzzyc->fam->ndims = fuzzyc->ninps;
	fuzzyc->fam->utype = (int*)malloc(fuzzyc->n_rules*sizeof(int));
	fuzzyc->fam->out_funcs = (int*)malloc(fuzzyc->n_rules*sizeof(int));

	for(i = 0;i < fuzzyc->n_rules;i++){
		fuzzyc->fam->utype[i] = 0;
		//fscanf(fh,"%d",&fuzzyc->fam->utype[i]);
		fscanf(fh,(i == fuzzyc->n_rules-1)?"%d\n":"%d ",&fuzzyc->fam->out_funcs[i]);
	}

	//Allocate & assign rule polynomials
	for(i = 0;i < fuzzyc->n_rules;i++){
		fuzzyc->rules[i].deg = fuzzyc->ninps+1;//Really not the degree but the numer of terms
		fuzzyc->rules[i].coef = (double*)malloc(fuzzyc->rules[i].deg*sizeof(double));

		for(j = 0;j < fuzzyc->rules[i].deg;j++)
			fscanf(fh,(i == fuzzyc->rules[i].deg-1)?"%lf\n":"%lf ",&fuzzyc->rules[i].coef[j]);
	}

	fclose(fh);
}
