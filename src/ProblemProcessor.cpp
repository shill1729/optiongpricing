#include <RcppEigen.h>
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


// Enter M or m-1
Rcpp::NumericMatrix ProblemProcessor::initial_cond(Rcpp::NumericMatrix u, Rcpp::NumericVector p)
{

  for (int j = 0; j < u.ncol(); j++)
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

// Assuming the region is already extended
Rcpp::NumericMatrix ProblemProcessor::jump_boundaries(Rcpp::NumericMatrix u, Rcpp::NumericVector ic, unsigned int N, unsigned int M, unsigned int jump_points)
{


  // Set boundary conditions
  for (unsigned int i = 0; i < N + 1; ++i)
  {
    // Lower half
    for (unsigned int j = 0; j <= jump_points; ++j)
    {
      u(i, j) = ic[j];

    }
    // Upper half
    for (unsigned int j = M + jump_points; j <= M + 2 * jump_points; ++j)
    {
      u(i, j) = ic[j];
    }

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

Eigen::MatrixXd ProblemProcessor::populate_jump_matrix(int i, Rcpp::NumericMatrix u, unsigned int N, unsigned int M, unsigned int jump_points)
{
  int L = jump_points;
  Eigen::MatrixXd jumpMatrix(M - 1, 2 * L + 1);
  // The jump coefficient matrix is a wrap around of the solution matrix
  for (unsigned int l = 0; l < M - 1; ++l)
  {
    for (int r = -L; r <= L; ++r)
    {
      jumpMatrix(l, L + r) = u(i - 1, l + L + 1 + r);
    }
  }
  return jumpMatrix;
}

Eigen::VectorXd ProblemProcessor::populate_jump_pdf(int L, double h, std::string distr, Rcpp::NumericVector param)
{
  Eigen::VectorXd jump_size_pdf(2*L+1);


  for (int i = -L; i < L; i++)
  {
    if(distr == "norm")
    {
      jump_size_pdf[i + L] = 2 * JumpSizePdf::pdf_norm(i * h, param);
    } else if(distr == "kou")
    {
      Rcpp::NumericVector jp = {param[0], param[1], param[2], 0, 0};
      jump_size_pdf[i + L] = 2 * JumpSizePdf::pdf_dkou(i * h, jp);
    } else if(distr == "dkou")
    {
      jump_size_pdf[i + L] = 2 * JumpSizePdf::pdf_dkou(i * h, param);
    }

  }
  jump_size_pdf[0] = 0.5 * jump_size_pdf[0];
  jump_size_pdf[2 * L] = 0.5 * jump_size_pdf[2 * L];
  return jump_size_pdf;
}

