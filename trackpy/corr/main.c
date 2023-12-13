#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "main.h"
#include "disks.h"
#include "dataio.h"

int main(int argc, char* argv[]) {

    const int N = atoi(argv[1]); //da impostare
    const double phi = atof(argv[2]);
    double const dr = 1.; //spacing
    char *hex_name;
    hex_name = argv[3];
    int i;
    int Neff;
    int ndiv;
   
    double lx = sqrt(N*PI/4/phi);
    double ly = lx;

    printf("%s conf with %d beads, density %lf and size %lf x %lf\n",hex_name,N,phi,lx,ly);

    double *rx;
    double *ry;
    rx = malloc(N*sizeof(double));
    ry = malloc(N*sizeof(double));

    double *hex_re;
    double *hex_im;

    hex_re = malloc(N*sizeof(double));
    hex_im = malloc(N*sizeof(double));

    printf("alloccated memory for hexatic order parameter\n");
 
    Neff = read_hexatic_conf(hex_name, rx, ry, hex_re, hex_im);
    //for (i=0;i<Neff;i++) printf("%lf %lf %lf %lf\n",rx[i],ry[i],hex[i][1],hex[i][2]);
    printf("Neff = %d",Neff);
    ndiv = hex_correlations(Neff, rx, ry, hex_re, hex_im, lx, ly, dr, hex_name); 

    free(rx);
    free(ry);
    free(hex_re);
    free(hex_im);

}
