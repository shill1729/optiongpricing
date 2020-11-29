# Geometric Compensated Poisson Process price dynamics:
# S_t=S_0 e^X_t where X_t=aN_t-bt for a Poisson process N_t

#' Compute European option prices via analytic formulas for Poisson driven prices.
#'
#' @param strike the strike price of the option
#' @param maturity the maturity in trading years until expiration
#' @param type type of option, "call" or "put"
#' @param spot the current spot price
#' @param rate the risk-neutral rate
#' @param div the continuous dividend yield rate
#' @param volat the volatility level
#'
#' @description {European option pricing via exact formulas under Levy-Poisson model dynamics for
#' vanilla options like calls and puts.}
#' @return numeric
#' @export analyticBlackScholesPrice
analyticBlackScholesPrice <- function(strike, maturity, type, spot, rate, div, volat)
{
  rn <- findistr::pblackscholes(strike, maturity, spot, rate, volat)
  fr <- findistr::pblackscholes(strike, maturity, spot, rate+volat^2, volat)
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


#' Compute European option prices via Monte-Carlo means under Black-Scholes dynamics
#'
#' @param strike the strike price of the option
#' @param maturity the maturity in trading years until expiration
#' @param type type of option, "call" or "put"
#' @param spot the current spot price
#' @param rate the risk-neutral rate
#' @param div the (continuous) dividend yield rate
#' @param volat the volatility level
#' @param n number of variates to use in Monte-Carlo mean
#'
#' @description {European option pricing via Monte-Carlo means under Black-Scholes dynamics for
#' vanilla options like calls and puts.}
#' @return numeric
#' @export monteCarloBlackScholes
monteCarloBlackScholes <- function(strike, maturity, type, spot, rate, div, volat, n = 10^6)
{
  s <- findistr::rblackscholes(n, maturity, spot, rate-div, volat)
  if(type == "put")
  {
    payoff <- pmax(strike-s, 0)

  } else if(type == "call")
  {
    payoff <- pmax(s-strike, 0)
  }
  return(exp(-rate*maturity)*mean(payoff))
}
