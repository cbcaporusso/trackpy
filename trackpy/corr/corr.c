#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "main.h"

int hex_correlations(int n, double *rx, double *ry, double *hex_re, double *hex_im, double lx, double ly, double dr, char *hex_name) {

    int i,j;
	double xij, yij;
    
	double dd;
    double psiRe_i,psiRe_j;
	double psiIm_i,psiIm_j;
    
	int shell;

    double max_lenght;
    max_lenght = lx;
	if (ly < max_lenght) max_lenght = ly;

    int xdiv, ydiv, ndiv;
    
	xdiv = (int)((.5*lx/dr)+1);
    ydiv = (int)((.5*ly/dr)+1);
    
	ndiv = xdiv;
	if (ydiv < ndiv) ndiv = ydiv;

    double **g;
    g = malloc(n*sizeof(double));
    for (i=0;i<n;i++) g[i] = malloc(n*sizeof(double));

	// print status
	fprintf(stderr,"Correlations for %s with %d particles in a box of size %f x %f with dr = %f and %d shells (max lenght = %f) ... " , hex_name, n, lx, ly, dr, ndiv, max_lenght);

	for(i=0;i<n;i+=1) {
		if(i%1000==0) fprintf(stderr,"%.2f\n",100.*i/n);
		psiRe_i = hex_re[i];
		psiIm_i = hex_im[i];

		for(j=0;j<n;j++) {
			//distanza con PBC
			//
			xij = rx[i]-rx[j];
			if (xij < -0.5*lx) xij+=lx;
			if (xij >  0.5*lx) xij-=lx;
			yij = ry[i]-ry[j];
			if (yij < -0.5*ly) yij+=ly;
			if (yij >  0.5*ly) yij-=ly;

			dd=sqrt((xij*xij)+(yij*yij));
			
			psiRe_j = hex_re[j];
			psiIm_j = hex_im[j];

			if(dd<.5*max_lenght) {
				shell=(int)(dd/dr);
				g[shell][1] += 1.;
				g[shell][2] += (psiRe_i*psiRe_j + psiIm_i*psiIm_j)/( psiRe_i*psiRe_i + psiIm_i*psiIm_i );
				g[shell][3] += (psiRe_j*psiIm_i - psiRe_i*psiIm_j)/( psiRe_i*psiRe_i + psiIm_i*psiIm_i);
			}
		}
		
	}

	FILE *fout;
	char corr_file[255];
	sprintf(corr_file,"%s.corr", hex_name);

	fout=fopen(corr_file,"w");

	for(i=0;i<ndiv;i++) {
		g[i][0]=i*dr+dr/2.;
		if(g[i][1]>0){
			fprintf(fout,"%lf %lf %lf\n",g[i][0],g[i][2]/g[i][1],g[i][3]/g[i][1]);
			//printf("%lf %lf %lf\n",g[i][0],g[i][2]/g[i][1],g[i][3]/g[i][1]);
		}
	}
    
    fclose(fout);
    free(g);
    
	return ndiv;

}
