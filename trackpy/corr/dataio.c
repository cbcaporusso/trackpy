#include <stdio.h>

#include "main.h"


int read_hexatic_conf(char *hex_name,double *rx,double *ry,double *hex_re, double *hex_im) {
	
	int bead_id;
	double xx,yy;
	double trash1,trash2;
	double re,img;
	int n;


	FILE *fin;
	fin = fopen(hex_name,"r");


	while(fscanf(fin,"%d %lf %lf %lf %lf %lf %lf",&bead_id,&xx,&yy,&trash1,&re,&img,&trash2)!=EOF) {
		printf("reading bead %d\n", bead_id);
		rx[n] = xx;
		ry[n] = yy;
		hex_re[n] = re;
		hex_im[n] = img;
		n++;
		
	}

	fclose(fin);
	
	return n;
}


