#ifndef PROBLEMPROCESSOR_H
#define PROBLEMPROCESSOR_H


class ProblemProcessor{
public:

  ProblemProcessor();
  virtual ~ProblemProcessor();
  static Rcpp::NumericVector discretize(double a, double b, unsigned int n);
  static Rcpp::NumericMatrix initial_cond(Rcpp::NumericMatrix u, Rcpp::NumericVector p, std::vector<int> res);
  static Rcpp::NumericMatrix boundary_cond(Rcpp::NumericMatrix u, Rcpp::NumericVector p, std::vector<int> res);
  static Rcpp::NumericMatrix post_step(unsigned int i, unsigned int L, Rcpp::NumericMatrix u, Rcpp::NumericVector sol, Rcpp::NumericVector ic, bool stopping);

};

#endif
