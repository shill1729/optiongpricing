#include <Rcpp.h>
#include "Tridiag.h"
#include "Payoffs.h"
#include "ProblemProcessor.h"

Rcpp::NumericMatrix cpoisSolver(bool american, Rcpp::NumericVector payoff, std::vector<int> res, double r, double a, double b, double h, double k)
{
  int N = res[0];
  int M = res[1];
  double lambda = (r+b)/(std::exp(a)-1.0);
  Rcpp::NumericMatrix u(N + 1, M + 1);
  //u = u.Zero(N + 1, M + 1);
  u = ProblemProcessor::initial_cond(u, payoff);
  u = ProblemProcessor::boundary_cond(u, payoff, res);
  double alpha = -(a*lambda-b)/(2*h);
  double beta = -r;
  double delta = -alpha;
  for (int i = 1; i < N + 1; i++)
  {
    // RHS of tridiagonal system
    Rcpp::NumericVector d(M - 1);
    for (int j = 0; j < M - 1; j++)
    {
      d[j] = u(i-1, j + 1);
    }
    d[0] += k*alpha * u(0, 0);
    d[M - 2] += k*delta * u(0, M);
    // Solving the system and repopulating
    Rcpp::NumericVector sol = Tridiag::thomasAlgorithm(-k * alpha, 1 - k * beta, -k * delta, d);
    u = ProblemProcessor::post_step(i, 0, u, sol, payoff, american);
  }
  return u;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix cpoisScheme(double strike, double maturity, std::string type, double spot, double r, double a, double b, std::vector<int> res, bool american = true)
{
  int N = res[0];
  int M = res[1];
  double lambda = (r+b)/(std::exp(a)-1.0);
  double B = maturity*lambda*a-b*maturity;
  double h = (2.0*B)/M;
  double k = (maturity/N);

  Rcpp::NumericVector x = ProblemProcessor::discretize(-B, B, M);
  Rcpp::NumericVector payoff = Payoffs::computePayoff(type, strike, spot, x, res);
  Rcpp::NumericMatrix u = cpoisSolver(american, payoff, res, r, a, b, h, k);
  return u;
}

//' Price options via Compensated Poisson PDE
//' @param strike the strike price of the option contract
//' @param maturity the time until maturity
//' @param type either "put" or "call
//' @param spot the spot price
//' @param r the risk free rate
//' @param a coefficient of the Poisson process
//' @param b coefficient of the compensating drift
//' @param res the grid resolution
//' @param american bool for american options
//'
//' @description {Compute European/American options via a finite difference solver for the
//' PDE for compensated Poisson log-returns.}
//' @return numeric
//' @export cpoisPDE
// [[Rcpp::export]]
double cpoisPDE(double strike, double maturity, std::string type, double spot, double r, double a, double b, std::vector<int> res, bool american = true)
{
  Rcpp::NumericMatrix u = cpoisScheme(strike, maturity, type, spot, r, a, b, res, american);
  return u(res[0], res[1]/2);
}
