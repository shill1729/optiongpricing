#include <Rcpp.h>
#include <cmath>
#include "JumpSizePdf.h"


double JumpSizePdf::pdf_norm(double x, Rcpp::NumericVector param)
{
  double mean = param[0];
  double sd = param[1];
  double f;
  f = (1 / (std::sqrt(2 * std::atan(1) * 4) * sd)) * std::exp(-std::pow(x - mean, 2.0) / (2 * std::pow(sd, 2.0)));
  return f;
}

double JumpSizePdf::pdf_exp(double x, double lambda)
{
  double f = 0.0;
  if (x < 0.0)
  {
    f = 0.0;
  }
  else if (x >= 0.0)
  {
    f = lambda * std::exp(-lambda * x);
  }
  return f;
}

double JumpSizePdf::mgf_exp(double t, double lambda)
{
  double m = 0.0;
  if (t < lambda)
  {
    m = lambda / (lambda - t);
  }
  else
  {
    abort();
    //Rcpp::print("MGF for EXP RVS argument must be less than lambda");
    //debug::print("MGF for Exp RVs' argument must be less than lambda");
  }
  return m;
}

double JumpSizePdf::pdf_dkou(double x, Rcpp::NumericVector param)
{
  double prob = param[0];
  double alpha = param[1];
  double beta = param[2];
  double ku = param[3];
  double kd = param[4];

  double a = 1.0 / alpha;
  double b = 1.0 / beta;
  double f = 0.0;
  double y;
  if (x >= ku)
  {
    y = x - ku;
    f = prob * pdf_exp(y, a);
  }
  else if (x < kd)
  {
    y = x - kd;
    f = (1.0 - prob) * pdf_exp(-y, b);

  }
  return f;
}

double JumpSizePdf::mgf_dkou(double t, Rcpp::NumericVector param)
{
  double prob = param[0];
  double alpha = param[1];
  double beta = param[2];
  double ku = param[3];
  double kd = param[4];

  double a = 1.0 / alpha;
  double b = 1.0 / beta;
  double q = 1.0 - prob;
  double m;
  m = q * std::exp(kd) * mgf_exp(-t, b) + prob * std::exp(ku) * mgf_exp(t, a);
  return m;
}
