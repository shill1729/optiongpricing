#include <Rcpp.h>
#include "Payoffs.h"

double Payoffs::put(double K, double S)
{
  return std::max<double>(K - S, 0.0);
}

double Payoffs::call(double K, double S)
{
  return std::max<double>(S - K, 0.0);
}

Rcpp::NumericVector Payoffs::computePayoff(std::string type, double strike, double spot, Rcpp::NumericVector x, std::vector<int> res)
{
  Rcpp::NumericVector payoff(res[1] + 1);
  for (int i = 0; i < res[0]+1; i++)
  {
    if (type == "put")
    {
      payoff[i] = put(strike, spot*std::exp(x[i]));
    }
    else if (type == "call")
    {
      payoff[i] = call(strike, spot * std::exp(x[i]));
    }
  }
  return payoff;
}
