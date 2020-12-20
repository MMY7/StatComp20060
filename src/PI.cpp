#include <Rcpp.h>
using namespace Rcpp;
//' @title PI
//' @description Use Monte Carlo to calculate PI using Rcpp
//' @param N the number of samples
//' @return PI
//' @examples
//' \dontrun{
//' pi <- piSugar(1000)
//' }
//' @export
// [[Rcpp::export]]
double piSugar(const int N) {
  RNGScope scope; // ensure RNG gets set/reset
  NumericVector x = runif(N);
  NumericVector y = runif(N);
  NumericVector d = sqrt(x*x+y*y);
  return 4.0 * sum(d < 1.0) / N;
}