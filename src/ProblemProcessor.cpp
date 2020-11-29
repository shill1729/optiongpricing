#include <Rcpp.h>
#include "ProblemProcessor.h"

Rcpp::NumericVector ProblemProcessor::discretize(double a, double b, unsigned int n)
{
  Rcpp::NumericVector x(n+1);
  double h = (b-a)/n;
  for(unsigned int i = 0; i<n+1; i++)
  {
    x[i] = a+i*h;
  }
  return x;
}

Rcpp::NumericMatrix ProblemProcessor::initial_cond(Rcpp::NumericMatrix u, Rcpp::NumericVector p, std::vector<int> res)
{
  for (int j = 0; j < res[1] + 1; j++)
  {
    u(0, j) = p[j];
  }
  return u;
}

Rcpp::NumericMatrix ProblemProcessor::boundary_cond(Rcpp::NumericMatrix u, Rcpp::NumericVector p, std::vector<int> res)
{
  int M = res[1];
  for (int i = 0; i < res[0] + 1; i++)
  {
    u(i, 0) = p[0];
    u(i, M) = p[M];
  }
  return u;
}

Rcpp::NumericMatrix ProblemProcessor::post_step(unsigned int i, unsigned int L, Rcpp::NumericMatrix u, Rcpp::NumericVector sol, Rcpp::NumericVector ic, bool stopping)
{
  unsigned int M = u.cols() - 1 - 2 * L;
  for (unsigned int j = 1 + L; j <= M - 1 + L; j++)
  {
    // For variational inequalities resulting from optional stopping problems
    if (stopping)
    {
      u(i, j) = std::max<double>(sol[j - 1 - L], ic[j]);

    }
    else if (!stopping)
    {
      u(i, j) = sol[j - 1 - L];
    }

  }
  return u;
}
