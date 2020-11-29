# Geometric Compensated Poisson Process price dynamics:
# S_t=S_0 e^X_t where X_t=aN_t-bt for a Poisson process N_t

#' Compute European option prices via analytic formulas for Poisson driven prices.
#'
#' @param strike the strike price of the option
#' @param maturity the maturity in trading years until expiration
#' @param type type of option, "call" or "put"
#' @param spot the current spot price
#' @param rate the risk-neutral rate
#' @param a the coefficient of the Poisson process
#' @param b the drift-compensation coefficient
#'
#' @description {European option pricing via exact formulas under Levy-Poisson model dynamics for
#' vanilla options like calls and puts.}
#' @return numeric
#' @export analyticPoissonPrice
analyticPoissonPrice <- function(strike, maturity, type, spot, rate, a, b)
{
  if(a == 0)
  {
    stop("'a' must be non-zero")
  }
  lambdaStar <- (rate+b)/(exp(a)-1)
  rn <- findistr::pgcpp(strike, maturity, spot, a, b, lambdaStar)
  fr <- findistr::pgcpp(strike, maturity, spot, a, b, exp(a)*lambdaStar)
  if(type == "put")
  {

    z <- strike*exp(-rate*maturity)*rn-spot*fr
    return(z)
  } else if(type == "call")
  {
    z <- spot*(1-fr)-strike*exp(-rate*maturity)*(1-rn)
    return(z)
  } else
  {
    stop("Invalid contract type")
  }
}


#' Compute European option prices via Monte-Carlo means under Levy-Poisson dynamics
#'
#' @param strike the strike price of the option
#' @param maturity the maturity in trading years until expiration
#' @param type type of option, "call" or "put"
#' @param spot the current spot price
#' @param rate the risk-neutral rate
#' @param a the coefficient of the Poisson process
#' @param b the drift-compensation coefficient
#' @param n number of variates to use in Monte-Carlo mean
#'
#' @description {European option pricing via Monte-Carlo means under Levy-Poisson model dynamics for
#' vanilla options like calls and puts.}
#' @return numeric
#' @export monteCarloPoisson
monteCarloPoisson <- function(strike, maturity, type, spot, rate, a, b, n = 10^6)
{
  if(a == 0)
  {
    stop("'a' must be non-zero")
  }
  lambda <- (rate+b)/(exp(a)-1)
  s <- findistr::rgcpp(n, maturity, spot, a, b, lambda)
  if(type == "put")
  {
    payoff <- pmax(strike-s, 0)

  } else if(type == "call")
  {
    payoff <- pmax(s-strike, 0)
  }
  return(exp(-rate*maturity)*mean(payoff))
}

