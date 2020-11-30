#include <RcppEigen.h>
#include "Tridiag.h"
#include "Payoffs.h"
#include "ProblemProcessor.h"

Rcpp::NumericMatrix imexSolver(bool american, Rcpp::NumericVector payoff, std::vector<int> res, double r, double q, double v, double lambda, std::string jumpDistr, Rcpp::NumericVector jumpParam, double h, double k)
{
  int N = res[0];
  int M = res[1];
  int L = M/2+1;
  int m = M+1+2*L;
  double eta = 0;

  if(jumpDistr == "norm")
  {
    eta = std::exp(jumpParam[0]+0.5*jumpParam[1]*jumpParam[1])-1;

  } else if(jumpDistr == "kou")
  {
    Rcpp::NumericVector jp = {jumpParam[0], jumpParam[1], jumpParam[2], 0, 0};
    eta = JumpSizePdf::mgf_dkou(1, jp) - 1;
  } else if(jumpDistr == "dkou")
  {
    eta = JumpSizePdf::mgf_dkou(1, jumpParam) - 1;
  }
  Eigen::VectorXd jump_size_pdf = ProblemProcessor::populate_jump_pdf(L, h, jumpDistr, jumpParam);

  Rcpp::NumericMatrix u(N + 1, m);
  //u = u.Zero(N + 1, M + 1);
  u = ProblemProcessor::initial_cond(u, payoff);
  u = ProblemProcessor::jump_boundaries(u, payoff, N, M, L);
  double a = (r - q - 0.5 * v * v-lambda*eta) / (2 * h);
  double b = (v * v) / (2 * h * h);
  double alpha = b - a;
  double beta = -r - 2 * b-lambda;
  double delta = a + b;
  // Eigen::MatrixXd jumpMatrix(M - 1, 2 * L + 1);
  for (int i = 1; i < N + 1; i++)
  {
    // Populate jump-matrix
    Eigen::MatrixXd jumpMatrix = ProblemProcessor::populate_jump_matrix(i, u, N, M, L);
    Eigen::VectorXd ju(M - 1);
    ju = jumpMatrix * jump_size_pdf;
    // RHS of tridiagonal system
    Rcpp::NumericVector d(M - 1);
    for (int j = 1; j <= M - 1; j++)
    {
      d[j-1] = u(i-1, j + L)+k*lambda*h*0.5*ju[j-1];
    }
    d[0] += k*alpha * u(0, L);
    d[M - 2] += k*delta * u(0, M+L);
    // Solving the system and repopulating
    Rcpp::NumericVector sol = Tridiag::thomasAlgorithm(-k * alpha, 1 - k * beta, -k * delta, d);
    u = ProblemProcessor::post_step(i, L, u, sol, payoff, american);
  }
  return u;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix imexScheme(double strike, double maturity, std::string type, double spot, double r, double q, double v, double lambda, std::string jumpDistr, Rcpp::NumericVector jumpParam, std::vector<int> res, bool american = true)
{
  int N = res[0];
  int M = res[1];
  int L = M/2+1;
  int m = M+1+2*L;
  double B = 0.5*v*M*std::sqrt(3*maturity/N);
  double h = (2*B)/M;
  double k = (maturity/N);
  B = B + L*h;
  Rcpp::NumericVector x = ProblemProcessor::discretize(-B, B, m-1);
  Rcpp::NumericVector payoff = Payoffs::computePayoff(type, strike, spot, x, res);
  Rcpp::NumericMatrix u = imexSolver(american, payoff, res, r, q, v, lambda, jumpDistr, jumpParam, h, k);
  return u;
}

//' Price options via Black-Scholes PDE
//' @param strike the strike price of the option contract
//' @param maturity the time until maturity
//' @param type either "put" or "call"
//' @param spot the spot price
//' @param r the risk free rate
//' @param q the continuous dividend yield
//' @param v the volatiltiy
//' @param lambda the mean number of jumps annually
//' @param jumpDistr name of jump-size distribution: "norm", "kou" or "dkou"
//' @param jumpParam vector of jump parameters, see details
//' @param res the time-space resolution
//' @param american bool for american options
//'
//' @description {Compute European/American options via a finite difference solver for the
//' jump diffusion PIDE.}
//' @return numeric
//' @export pricerPIDE
// [[Rcpp::export]]
double pricerPIDE(double strike, double maturity, std::string type, double spot, double r, double q, double v, double lambda, std::string jumpDistr, Rcpp::NumericVector jumpParam, std::vector<int> res, bool american = true)
{

  int M = res[1];
  int L = M/2+1;
  int m = M+1+2*L;
  int sj = round(m/2);

  Rcpp::NumericMatrix u = imexScheme(strike, maturity, type, spot, r, q, v, lambda, jumpDistr, jumpParam, res, american);
  return u(res[0], sj);
}
