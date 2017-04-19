#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
double predprobCpp(int y, int n, int nmax, double alpha_e, double beta_e, double p_s, double theta_t) {

  double prob = 0.0;
  double eps = std::numeric_limits<double>::epsilon();
  double pxy;

  for (int x = 0; x < nmax - n + 1; x++) {
    pxy = (1.0 - R::pbeta(p_s, alpha_e + y + x, beta_e + nmax - y - x, 1, 0));
    if (pxy > theta_t || std::abs(pxy - theta_t) < eps) {
      prob += exp(
        R::lchoose(nmax - n, x) +
          R::lbeta(alpha_e + y + x, beta_e + nmax - y - x) -
          R::lbeta(alpha_e + y, beta_e + n - y)
      );
    }
  }

  return prob;

}
