#ifndef JUMPSIZEPDF_H
#define JUMPSIZEPDF_H

class JumpSizePdf
{
public:
  // PDF for normally distributed jump sizes
  static double pdf_norm(double x, Rcpp::NumericVector param);


  // PDF and MGF for exponential RVs
  static double pdf_exp(double x, double lambda);
  static double mgf_exp(double t, double lambda);

  // PDF for displaced double-Kou distributed jump-sizes
  static double pdf_dkou(double x, Rcpp::NumericVector param);
  static double mgf_dkou(double t, Rcpp::NumericVector param);

};

#endif
