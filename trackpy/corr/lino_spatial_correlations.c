#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.14159265358979323846
#define maxbeads 25500000
#define box_l 0
#define box_r 320.848419153
#define box_b 0
#define box_t 320.848419153

/* SYNOPSIS
 * 
 * ./executable.x <file>.dat.hexatic
 *
 * <file>.dat.hexatic :	1 x1 z1 frac1 Re_hex1 Im_hex1
			2 x2 z2 frac2 Re_hex2 Im_hex2
			...
 *
 */

double **allodouble(int nx,int ny){

	int i;
	double **m;

	m=(double **)malloc(nx*sizeof(double *));
	if(!m){
		printf("allocation failure 1 in allodouble\n");
		exit(1);
	}

	m[0]=(double *)malloc(nx*ny*sizeof(double));
	if(!m[0]){
		printf("allocation failure 2 in allodouble\n");
		exit(1);
	}

	for(i=0;i<nx;i++)m[i]=m[0]+i*ny;

	return m;

}


int main(int argc, char *argv[]) {

	char namein[255], nameout[255];
	FILE *fpin, *fpout;
	double xx, zz, xi, zi, xj, zj;
	double trash1, trash2, dist, re, img, check;
	int bead_id;

	int shell;
        int nbeads=0;
        int i,j,k,n;

	double x_box_size = (double)(box_r-box_l);
	double y_box_size = (double)(box_t-box_b);
	double max_length = x_box_size;
	if (y_box_size>max_length) max_length = y_box_size;
//	if(box_r-box_l != box_t-box_b) {
//		fprintf(stderr,"no square box ?\n");
//		exit(0);
//	}

	double dr = 1.;						//spacing for the histogram: dr*sigma
//	fprintf(stderr,"    > dr = %lf\n",dr);
	int x_ndiv = (int)((0.5*x_box_size)/dr) + 1;			//number of division per-side for the histogram
	int y_ndiv = (int)((0.5*y_box_size)/dr) + 1;			//number of division per-side for the histogram
	int ndiv = x_ndiv;
	if (y_ndiv>ndiv) ndiv = y_ndiv;
//	fprintf(stderr,"    > ndiv = %d\n",ndiv);

	double **order;							//allocation for the input
	order = allodouble(maxbeads,4);
	
//	double phi = XXXphiXXX;
//	double a = sqrt(PI/(2.*phi*sqrt(3)));
//	double kx_perf = 2*PI/a;
//	double kz_perf = -2*PI/(a*sqrt(3));

	int n_peaks = (int)((argc-2)/2.);
//	int n_peaks = 1;

	double **g;							//allocation for the output
	g = allodouble(ndiv,4+2*n_peaks);
	for(i=0;i<ndiv;i++) {
                for(k=0;k<4+2*n_peaks;k++) g[i][k]=0.;
        }
//	for (k=0;k<3;k++) {
//		g[k] = allodouble(ndiv[k],4+2*n_peaks);
//		for(i=0;i<ndiv[k];i++) {
//			for(j=0;j<4+2*n_peaks;j++) g[k][i][j]=0.;
//		}
//	}

        double kx[20], kz[20];

	fprintf(stderr,"    > n peaks = %d:",n_peaks);
	for (i=0;i<n_peaks;i++) {
		sscanf(argv[2+i],"%lf",&kx[i]);
		sscanf(argv[2+i+1],"%lf",&kz[i]);
		fprintf(stderr," (%lf,%lf)",kx[i],kz[i]);
	}
	fprintf(stderr,"\n");

//	kx[0] = kx_perf;
//	kz[0] = kz_perf;

//	fprintf(stderr,"    > n peaks = %d: (%lf,%lf)",n_peaks,kx[0],kz[0]);
//	for (k=0;k<n_peaks-1;k++) {
//		sscanf(argv[2+2*k], "%lf", &kx[k+1]);
//		sscanf(argv[2+2*k+1], "%lf", &kz[k+1]);
//		fprintf(stderr," (%lf,%lf)",kx[k+1],kz[k+1]);
//	}
//	kx[1]=kx[0]-1.; kz[1]=kz[0]; fprintf(stderr," (%lf,%lf)",kx[1],kz[1]);
//	kx[2]=kx[0]/(2.*PI); kz[2]=kz[0]/(2.*PI); fprintf(stderr," (%lf,%lf)",kx[2],kz[2]);
//	fprintf(stderr,"\n");

	if (argc < 2) {
		fprintf(stderr,"usage: %s <file>.dat.hexatic \n",argv[0]);
		exit(1);
	}

	sprintf(namein,"%s",argv[1]);
	fprintf(stderr,"    > Opening file: %s\n",argv[1]);
	fpin = fopen(namein,"r");

	while(fscanf(fpin,"%d %lf %lf %lf %lf %lf %lf",&bead_id,&xx,&zz,&trash1,&re,&img,&trash2)!=EOF) {

		order[nbeads][0]=xx;
		order[nbeads][1]=zz;
		order[nbeads][2]=re;
		order[nbeads][3]=img;

		nbeads++;
	}

	fclose(fpin);

	fprintf(stderr,"    > Number of beads: %4d. \n",nbeads);

	for(i=0;i<nbeads;i++) {

		xi=order[i][0];
		zi=order[i][1];

		for(j=0;j<nbeads;j++) {

			xj=order[j][0];
			if(xj-xi<-0.5*x_box_size) xj = xj+x_box_size;
			if(xj-xi>0.5*x_box_size) xj = xj-x_box_size;
			zj=order[j][1];
			if(zj-zi<-0.5*y_box_size) zj = zj+y_box_size;
			if(zj-zi>0.5*y_box_size) zj = zj-y_box_size;

			dist=sqrt((xi-xj)*(xi-xj) + (zi-zj)*(zi-zj));

			if(dist<0.5*max_length){
				shell=(int)(dist/dr);
				g[shell][1]+=1.;
				g[shell][2]+= (order[i][2]*order[j][2] + order[i][3]*order[j][3])/(order[i][2]*order[i][2] + order[i][3]*order[i][3]);
                                g[shell][3]+= (order[j][2]*order[i][3] - order[i][2]*order[j][3])/(order[i][2]*order[i][2] + order[i][3]*order[i][3]);
				for(k=0;k<n_peaks;k++) {
					g[shell][4+2*k]+= cos(kx[k]*(xj-xi)+kz[k]*(zj-zi));
					g[shell][4+2*k+1]+= sin(kx[k]*(xj-xi)+kz[k]*(zj-zi));
				}
			}

		}

		if (i%10000==0) fprintf(stderr,"    > %d of %d\n",i,nbeads);

	}

	sprintf(nameout,"%s.correlations_dr%3.1lf",namein,dr);
	fpout = fopen(nameout,"w");
	for(i=0;i<ndiv;i++){
		g[i][0]=i*dr+dr/2.;
		fprintf(fpout,"%8.4lf %6.0lf",g[i][0],g[i][1]);
		if(g[i][1]!=0){
			for(k=2;k<4+2*n_peaks;k++) fprintf(fpout," %7.4lf",g[i][k]/g[i][1]);
		}
		else for(k=2;k<4+2*n_peaks;k++) fprintf(fpout," %7.4lf",0.);
		fprintf(fpout,"\n");
	}
	fclose(fpout);

	free(order);
//	for (n=0;n<3;n++) free(g[n]);
	free(g);

	return 0;

}
