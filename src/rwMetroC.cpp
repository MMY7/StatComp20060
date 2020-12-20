#include <Rcpp.h>
using namespace Rcpp;

//' @title  The density function of the standard Laplace distribution 
//' @description Calculate density function of the standard Laplace distribution using Rcpp
//' @param x vector of quantiles.
//' @return the density
//' @examples
//' \dontrun{
//' f <- pdfLaplaceC(1)
//' }
//' @export
// [[Rcpp::export]]

double pdfLaplaceC(double x) 
{
  return exp(-1*abs(x))/2;
 }

#include <Rcpp.h>
using namespace Rcpp;

//' @title A random walk Metropolis sampler for the standard Laplace distribution using Rcpp
//' @description A random walk Metropolis sampler for generating the standard Laplace distribution using Rcpp
//' @param sigma the variance
//' @param N the length of the chain
//' @param xinitial the initial value
//' @return a list with two components:x and k give the chain and the number of candidate points rejected
//' @examples
//' \dontrun{
//' rw<- Metro.rw(2,3000,25)
//' }
//' @export
// [[Rcpp::export]]
List rwMetroC(double sigma, int N, double xinitial){
  NumericVector u = runif(N);
  NumericVector x(N); 
  x[0] = xinitial;
  int k = 0;
  NumericVector y;
  for(int i = 1;i <= (N-1);i++){
    y = rnorm(1, x[i-1], sigma);
    if (u[i] <= (pdfLaplaceC(y[0])/pdfLaplaceC(x[i-1]))) {
      x[i] = y[0];
    }
    else {
      x[i] = x[i-1];
	  k++;
	}
  }
  List res = List::create(x,k);
  return (res);
} 