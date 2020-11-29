#ifndef PAYOFFS_H
#define PAYOFFS_H

class Payoffs{
public:

  Payoffs();
  virtual ~Payoffs();

  static double put(double K, double S);
  static double call(double K, double S);
  static Rcpp::NumericVector computePayoff(std::string type, double strike, double spot, Rcpp::NumericVector x, std::vector<int> res);
};


#endif
