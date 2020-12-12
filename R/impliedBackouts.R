
#' Back out implied rate via Black-Scholes PDE
#'
#' @param fee the given market price of the option
#' @param strike the given strike price of the option
#' @param maturity the maturity of the option in trading years until expiration
#' @param type the type of option: "call" or "put"
#' @param spot the spot price of the underlying
#' @param volat the volatility level
#'
#' @description {Call to base r \code{uniroot} to backout the implied rate of an option given
#'  all other parameters of the Black-Scholes model.}
#' @return numeric
#' @export impliedRate
impliedRate <- function(fee, strike, maturity, type, spot, volat)
{
  f <- function(r)
  {
    optionpricing::blackScholesPDE(strike, maturity, type, spot, r, 0, volat, res = c(400, 400), american = TRUE)-fee
  }

  z <- stats::uniroot(f, c(-100, 100))
  return(z$root)
}


#' Back out implied volatility via Black-Scholes PDE
#'
#' @param fee the given market price of the option
#' @param strike the given strike price of the option
#' @param maturity the maturity of the option in trading years until expiration
#' @param type the type of option: "call" or "put"
#' @param spot the spot price of the underlying
#' @param rate the risk-free rate
#'
#' @description {Call to base r \code{uniroot} to backout the implied volatility of an option given
#'  all other parameters of the Black-Scholes model.}
#' @return numeric
#' @export impliedVolatility
impliedVolatility <- function(fee, strike, maturity, type, spot, rate)
{
  f <- function(v)
  {
    optionpricing::blackScholesPDE(strike, maturity, type, spot, rate, 0, v, res = c(400, 400), american = TRUE)-fee
  }

  z <- stats::uniroot(f, c(0.01, 5))
  return(z$root)
}
