#ifndef PROBLEMPROCESSOR_H
#define PROBLEMPROCESSOR_H
#include <RcppEigen.h>
#include "JumpSizePdf.h"


class ProblemProcessor{
public:

  ProblemProcessor();
  virtual ~ProblemProcessor();
  static Rcpp::NumericVector discretize(double a, double b, unsigned int n);
  static Rcpp::NumericMatrix initial_cond(Rcpp::NumericMatrix u, Rcpp::NumericVector p);
  static Rcpp::NumericMatrix boundary_cond(Rcpp::NumericMatrix u, Rcpp::NumericVector p, std::vector<int> res);
  static Rcpp::NumericMatrix jump_boundaries(Rcpp::NumericMatrix u, Rcpp::NumericVector ic, unsigned int N, unsigned int M, unsigned int jump_points);
  static Rcpp::NumericMatrix post_step(unsigned int i, unsigned int L, Rcpp::NumericMatrix u, Rcpp::NumericVector sol, Rcpp::NumericVector ic, bool stopping);
  static Eigen::MatrixXd populate_jump_matrix(int i, Rcpp::NumericMatrix u, unsigned int N, unsigned int M, unsigned int jump_points);
  static Eigen::VectorXd populate_jump_pdf(int L, double h, std::string distr, Rcpp::NumericVector param);

};

#endif
