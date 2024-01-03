#!/usr/bin/env awk -f

function PBC_1d(x2,x1,L) {			
            dx12 = x1-x2;
			if (dx12 < -0.5*L) dx12+=L;
			if (dx12 >  0.5*L) dx12-=L;
            return dx12
            }

BEGIN {
    i=0;
    j=0;
    l=0;
    fnr=0;
    L=254.364110000000012;
    dt
}

{

if (FNR==1) fnr++

if (fnr == 1) {
    if (FNR<=9) {
        header[l]=$0
        l++
    }
    if (FNR>9) {
        id1[i]    =   $1
        x1[$1]    =   $2
        y1[$1]    =   $3
        vx1[$1]   =   $4
        vy1[$1]   =   $5
        i++
    }
}

if (fnr == 2) {
    if (FNR>9) {

        id2[j]    =   $1
        x2[$1]    =   $2
        y2[$1]    =   $3
        vx2[$1]   =   $4
        vy2[$1]   =   $5
        j++
    }
}

}


END{
    id1[i]=id1[i-1]+1
    for (k=0;k<l;k++){
        print header[k]
    }
    for (k=1;k<=i;k++) {

        dx=PBC_1d(x2[k],x1[k],L)
        dy=PBC_1d(y2[k],y1[k],L)

        print id1[k]-1,x1[k],y1[k],dx/dt,dy/dt
    }
}