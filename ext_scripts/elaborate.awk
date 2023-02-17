#!/usr/bin/awk -f

BEGIN{PI = atan2(0, -1)}

{

	x[$1]=$2;
	y[$1]=$3;
	ro[$1]=PI/4./$4;
	nfield[$1]=NF;
	for(i=5;i<=NF;i++)neigh[$1,i]=$i;
	npart++;

}

END{

	for(i=1;i<=npart;i++){

		k=0;
		psi_re=0;
		psi_im=0;
		for(j=5;j<=nfield[i];j++){

			if(neigh[i,j]>=0){

				dx=x[neigh[i,j]]-x[i]
				dy=y[neigh[i,j]]-y[i]	
#				bc
				if(dx>L/2.) dx=dx-L
				if(dx<-L/2.) dx=dx+L
				if(dy>L/2.) dy=dy-L
				if(dy<-L/2.) dy=dy+L
				### calcoliamo angolo tra -180 e 180 usando atan2 (vanno prese direzioni normalizzate!!) ###		
				dist=sqrt(dx**2+dy**2)
				alpha=atan2(dy/dist,dx/dist)
				if(alpha<0) alpha+=2*PI
				psi_re+=cos(alpha*6.)
				psi_im+=sin(alpha*6.)
				k++
#				print i,x[i],y[i],neigh[i],dist,cosalpha,alpha,k,psi_re,psi_im

			}	

		}

		psi_re=psi_re/k
		psi_im=psi_im/k
		print i,x[i],y[i],ro[i],psi_re,psi_im,sqrt(psi_re*psi_re+psi_im*psi_im)

	}

}
