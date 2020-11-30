// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// implicitScheme
Rcpp::NumericMatrix implicitScheme(double strike, double maturity, std::string type, double spot, double r, double q, double v, std::vector<int> res, bool american);
RcppExport SEXP _optionpricing_implicitScheme(SEXP strikeSEXP, SEXP maturitySEXP, SEXP typeSEXP, SEXP spotSEXP, SEXP rSEXP, SEXP qSEXP, SEXP vSEXP, SEXP resSEXP, SEXP americanSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type strike(strikeSEXP);
    Rcpp::traits::input_parameter< double >::type maturity(maturitySEXP);
    Rcpp::traits::input_parameter< std::string >::type type(typeSEXP);
    Rcpp::traits::input_parameter< double >::type spot(spotSEXP);
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    Rcpp::traits::input_parameter< double >::type v(vSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type res(resSEXP);
    Rcpp::traits::input_parameter< bool >::type american(americanSEXP);
    rcpp_result_gen = Rcpp::wrap(implicitScheme(strike, maturity, type, spot, r, q, v, res, american));
    return rcpp_result_gen;
END_RCPP
}
// blackScholesPDE
double blackScholesPDE(double strike, double maturity, std::string type, double spot, double r, double q, double v, std::vector<int> res, bool american);
RcppExport SEXP _optionpricing_blackScholesPDE(SEXP strikeSEXP, SEXP maturitySEXP, SEXP typeSEXP, SEXP spotSEXP, SEXP rSEXP, SEXP qSEXP, SEXP vSEXP, SEXP resSEXP, SEXP americanSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type strike(strikeSEXP);
    Rcpp::traits::input_parameter< double >::type maturity(maturitySEXP);
    Rcpp::traits::input_parameter< std::string >::type type(typeSEXP);
    Rcpp::traits::input_parameter< double >::type spot(spotSEXP);
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    Rcpp::traits::input_parameter< double >::type v(vSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type res(resSEXP);
    Rcpp::traits::input_parameter< bool >::type american(americanSEXP);
    rcpp_result_gen = Rcpp::wrap(blackScholesPDE(strike, maturity, type, spot, r, q, v, res, american));
    return rcpp_result_gen;
END_RCPP
}
// cpoisScheme
Rcpp::NumericMatrix cpoisScheme(double strike, double maturity, std::string type, double spot, double r, double a, double b, std::vector<int> res, bool american);
RcppExport SEXP _optionpricing_cpoisScheme(SEXP strikeSEXP, SEXP maturitySEXP, SEXP typeSEXP, SEXP spotSEXP, SEXP rSEXP, SEXP aSEXP, SEXP bSEXP, SEXP resSEXP, SEXP americanSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type strike(strikeSEXP);
    Rcpp::traits::input_parameter< double >::type maturity(maturitySEXP);
    Rcpp::traits::input_parameter< std::string >::type type(typeSEXP);
    Rcpp::traits::input_parameter< double >::type spot(spotSEXP);
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type res(resSEXP);
    Rcpp::traits::input_parameter< bool >::type american(americanSEXP);
    rcpp_result_gen = Rcpp::wrap(cpoisScheme(strike, maturity, type, spot, r, a, b, res, american));
    return rcpp_result_gen;
END_RCPP
}
// cpoisPDE
double cpoisPDE(double strike, double maturity, std::string type, double spot, double r, double a, double b, std::vector<int> res, bool american);
RcppExport SEXP _optionpricing_cpoisPDE(SEXP strikeSEXP, SEXP maturitySEXP, SEXP typeSEXP, SEXP spotSEXP, SEXP rSEXP, SEXP aSEXP, SEXP bSEXP, SEXP resSEXP, SEXP americanSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type strike(strikeSEXP);
    Rcpp::traits::input_parameter< double >::type maturity(maturitySEXP);
    Rcpp::traits::input_parameter< std::string >::type type(typeSEXP);
    Rcpp::traits::input_parameter< double >::type spot(spotSEXP);
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type res(resSEXP);
    Rcpp::traits::input_parameter< bool >::type american(americanSEXP);
    rcpp_result_gen = Rcpp::wrap(cpoisPDE(strike, maturity, type, spot, r, a, b, res, american));
    return rcpp_result_gen;
END_RCPP
}
// imexScheme
Rcpp::NumericMatrix imexScheme(double strike, double maturity, std::string type, double spot, double r, double q, double v, double lambda, std::string jumpDistr, Rcpp::NumericVector jumpParam, std::vector<int> res, bool american);
RcppExport SEXP _optionpricing_imexScheme(SEXP strikeSEXP, SEXP maturitySEXP, SEXP typeSEXP, SEXP spotSEXP, SEXP rSEXP, SEXP qSEXP, SEXP vSEXP, SEXP lambdaSEXP, SEXP jumpDistrSEXP, SEXP jumpParamSEXP, SEXP resSEXP, SEXP americanSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type strike(strikeSEXP);
    Rcpp::traits::input_parameter< double >::type maturity(maturitySEXP);
    Rcpp::traits::input_parameter< std::string >::type type(typeSEXP);
    Rcpp::traits::input_parameter< double >::type spot(spotSEXP);
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    Rcpp::traits::input_parameter< double >::type v(vSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< std::string >::type jumpDistr(jumpDistrSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type jumpParam(jumpParamSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type res(resSEXP);
    Rcpp::traits::input_parameter< bool >::type american(americanSEXP);
    rcpp_result_gen = Rcpp::wrap(imexScheme(strike, maturity, type, spot, r, q, v, lambda, jumpDistr, jumpParam, res, american));
    return rcpp_result_gen;
END_RCPP
}
// pricerPIDE
double pricerPIDE(double strike, double maturity, std::string type, double spot, double r, double q, double v, double lambda, std::string jumpDistr, Rcpp::NumericVector jumpParam, std::vector<int> res, bool american);
RcppExport SEXP _optionpricing_pricerPIDE(SEXP strikeSEXP, SEXP maturitySEXP, SEXP typeSEXP, SEXP spotSEXP, SEXP rSEXP, SEXP qSEXP, SEXP vSEXP, SEXP lambdaSEXP, SEXP jumpDistrSEXP, SEXP jumpParamSEXP, SEXP resSEXP, SEXP americanSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type strike(strikeSEXP);
    Rcpp::traits::input_parameter< double >::type maturity(maturitySEXP);
    Rcpp::traits::input_parameter< std::string >::type type(typeSEXP);
    Rcpp::traits::input_parameter< double >::type spot(spotSEXP);
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    Rcpp::traits::input_parameter< double >::type v(vSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< std::string >::type jumpDistr(jumpDistrSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type jumpParam(jumpParamSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type res(resSEXP);
    Rcpp::traits::input_parameter< bool >::type american(americanSEXP);
    rcpp_result_gen = Rcpp::wrap(pricerPIDE(strike, maturity, type, spot, r, q, v, lambda, jumpDistr, jumpParam, res, american));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_optionpricing_implicitScheme", (DL_FUNC) &_optionpricing_implicitScheme, 9},
    {"_optionpricing_blackScholesPDE", (DL_FUNC) &_optionpricing_blackScholesPDE, 9},
    {"_optionpricing_cpoisScheme", (DL_FUNC) &_optionpricing_cpoisScheme, 9},
    {"_optionpricing_cpoisPDE", (DL_FUNC) &_optionpricing_cpoisPDE, 9},
    {"_optionpricing_imexScheme", (DL_FUNC) &_optionpricing_imexScheme, 12},
    {"_optionpricing_pricerPIDE", (DL_FUNC) &_optionpricing_pricerPIDE, 12},
    {NULL, NULL, 0}
};

RcppExport void R_init_optionpricing(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
