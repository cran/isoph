#include <R.h>
#include <math.h>

/*1. RPA_ti*/
void RPA_ti(int *n, int *nt, int *m, double *z, double *zobs, double *Y, int *dN, double *Y2, int *dN2)
{
    int i, j, h;

    for(i=0; i<*n; i++){
        h=0;
        if (z[i]<zobs[1]){
            h=0;
        }else if (z[i]>=zobs[*m-1]){
            h=*m-1;
        }else{
            do{
                h=h+1;
            }while( !(zobs[h]<=z[i] && z[i]<zobs[h+1]) ) ;
        }
        for(j=0; j<*nt; j++){
            Y2[h+j*(*m)] = Y2[h+j*(*m)]+ Y[i+j*(*n)];
            dN2[h+j*(*m)]=dN2[h+j*(*m)]+dN[i+j*(*n)];
        }
   }
}


/*2. RPA_tid*/
void RPA_td(int *n, int *nt, int *m, double *zt, double *zobs, double *Y, int *dN, double *Y2, int *dN2)
{
    int i, j, h;

    for(i=0; i<*n; i++){
        for(j=0; j<*nt; j++){
            if( (Y[i+j*(*n)]==1) ){ /*since NA for zt is not allowed*/
                h=0;
                if (zt[i+j*(*n)]<zobs[1]){
                    h=0;
                }else if (zt[i+j*(*n)]>=zobs[*m-1]){
                    h=*m-1;
                }else{
                    do{
                        h=h+1;
                    }while( !(zobs[h]<=zt[i+j*(*n)] && zt[i+j*(*n)]<zobs[h+1]) ) ;
                }
                Y2[h+j*(*m)] = Y2[h+j*(*m)]+ Y[i+j*(*n)];
                dN2[h+j*(*m)]=dN2[h+j*(*m)]+dN[i+j*(*n)];
            }
        }
    }
}

