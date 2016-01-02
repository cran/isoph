#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector BTFft(int m, int n, List interval, NumericVector z, NumericVector psi){

	NumericVector psifull(n);

	int i, j, h;

    psifull(0)=psi(0);
    psifull(n-1)=psi(m-1);

    for(i=1; i<n-1; i++){
        for(h=0; h<m; h++){
            SEXP interval2= interval[h];
            Rcpp::NumericVector interval3(interval2);
            if(interval3(0)<=z(i) && z(i)<interval3(1)){
                psifull(i)=psi(h);
            }
        }
    }

    return Rcpp::wrap(psifull);
}

