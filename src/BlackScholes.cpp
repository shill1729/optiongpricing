#include <Rcpp.h>
#include "Tridiag.h"

double put(double K, double S)
{
  return std::max<double>(K - S, 0.0);
}

double call(double K, double S)
{
  return std::max<double>(S - K, 0.0);
}

Rcpp::NumericVector computePayoff(std::string type, double strike, double spot, Rcpp::NumericVector x, std::vector<int> res)
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

Rcpp::NumericVector discretize(double a, double b, unsigned int n)
{
  Rcpp::NumericVector x(n+1);
  double h = (b-a)/n;
  for(unsigned int i = 0; i<n+1; i++)
  {
    x[i] = a+i*h;
  }
  return x;
}

Rcpp::NumericMatrix initial_cond(Rcpp::NumericMatrix u, Rcpp::NumericVector p, std::vector<int> res)
{
  for (int j = 0; j < res[1] + 1; j++)
  {
    u(0, j) = p[j];
  }
  return u;
}

Rcpp::NumericMatrix boundary_cond(Rcpp::NumericMatrix u, Rcpp::NumericVector p, std::vector<int> res)
{
  int M = res[1];
  for (int i = 0; i < res[0] + 1; i++)
  {
    u(i, 0) = p[0];
    u(i, M) = p[M];
  }
  return u;
}

Rcpp::NumericMatrix post_step(unsigned int i, unsigned int L, Rcpp::NumericMatrix u, Rcpp::NumericVector sol, Rcpp::NumericVector ic, bool stopping)
{
  unsigned int M = u.cols() - 1 - 2 * L;
  for (unsigned int j = 1 + L; j <= M - 1 + L; j++)
  {
    // For variational inequalities resulting from optional stopping problems
    if (stopping)
    {
      u(i, j) = std::max<double>(sol[j - 1 - L], ic[j]);

    }
    else if (!stopping)
    {
      u(i, j) = sol[j - 1 - L];
    }

  }
  return u;
}

Rcpp::NumericMatrix implicitSolver(bool american, Rcpp::NumericVector payoff, std::vector<int> res, double r, double q, double v, double h, double k)
{
  int N = res[0];
  int M = res[1];
  Rcpp::NumericMatrix u(N + 1, M + 1);
  //u = u.Zero(N + 1, M + 1);
  u = initial_cond(u, payoff, res);
  u = boundary_cond(u, payoff, res);
  double a = (r - q - 0.5 * v * v) / (2 * h);
  double b = (v * v) / (2 * h * h);
  double alpha = b - a;
  double beta = -r - 2 * b;
  double delta = a + b;
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
    u = post_step(i, 0, u, sol, payoff, american);
  }
  return u;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix implicitScheme(double strike, double maturity, std::string type, double spot, double r, double q, double v, std::vector<int> res, bool american = true)
{
  int N = res[0];
  int M = res[1];
  double B = 0.5*v*M*std::sqrt(3*maturity/N);
  double h = (2*B)/M;
  double k = (maturity/N);
  Rcpp::NumericVector x = discretize(-B, B, M);
  Rcpp::NumericVector payoff = computePayoff(type, strike, spot, x, res);
  Rcpp::NumericMatrix u = implicitSolver(american, payoff, res, r, q, v, h, k);
  return u;
}

//' Price options via Black-Scholes PDE
//' @param strike the strike price of the option contract
//' @param maturity the time until maturity
//' @param type either "put" or "call
//' @param spot the spot price
//' @param r the risk free rate
//' @param q the continuous dividend yield
//' @param v the volatiltiy
//' @param res the grid resolution
//' @param american bool for american options
//'
//' @description {Compute European/American options via a finite difference solver for the
//' Black-Scholes PDE.}
//' @return numeric
//' @export blackScholesPDE
// [[Rcpp::export]]
double blackScholesPDE(double strike, double maturity, std::string type, double spot, double r, double q, double v, std::vector<int> res, bool american = true)
{
  Rcpp::NumericMatrix u = implicitScheme(strike, maturity, type, spot, r, q, v, res, american);
  return u(res[0], res[1]/2);
}


Rcpp::NumericMatrix cpoisSolver(bool american, Rcpp::NumericVector payoff, std::vector<int> res, double r, double a, double b, double h, double k)
{
  int N = res[0];
  int M = res[1];
  double lambda = (r+b)/(std::exp(a)-1);
  Rcpp::NumericMatrix u(N + 1, M + 1);
  //u = u.Zero(N + 1, M + 1);
  u = initial_cond(u, payoff, res);
  u = boundary_cond(u, payoff, res);
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
    u = post_step(i, 0, u, sol, payoff, american);
  }
  return u;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix cpoisScheme(double strike, double maturity, std::string type, double spot, double r, double a, double b, std::vector<int> res, bool american = true)
{
  int N = res[0];
  int M = res[1];
  double lambda = (r+b)/(std::exp(a)-1);
  double B = maturity*lambda*a-b*maturity;
  double h = (2*B)/M;
  double k = (maturity/N);

  Rcpp::NumericVector x = discretize(-B, B, M);
  Rcpp::NumericVector payoff = computePayoff(type, strike, spot, x, res);
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
