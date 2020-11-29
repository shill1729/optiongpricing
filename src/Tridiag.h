#ifndef TRIDIAG_H
#define TRIDIAG_H


class Tridiag{
public:

  Tridiag();
  virtual ~Tridiag();

  // Solve tridiagonal linear systems
  static Rcpp::NumericVector thomasAlgorithm(Rcpp::NumericVector a, Rcpp::NumericVector b, Rcpp::NumericVector c, Rcpp::NumericVector d);

  // For constant tridiagonals
  static Rcpp::NumericVector thomasAlgorithm(double a, double b, double c, Rcpp::NumericVector d);
};

#endif
