#' Merton's jump-diffusion model analytic price formula for European options
#'
#' @param strike The strike price
#' @param maturity the year until expirations
#' @param type What type of option, put or call
#' @param spot the spot price
#' @param rate the drift rate
#' @param div The dividend yield rate
#' @param volat the diffusion coefficients of each component (a vector)
#' @param lambda the yearly jump-rate
#' @param jm the jump mean
#' @param jv the jump volatility
#' @param num_terms the number of terms to sum in the Poisson expectation
#'
#' @description {Semi-nalytic series formula for European options under Merton's jump diffusion dynamics}
#' @return numeric
#' @export analyticMertonPrice
analyticMertonPrice <- function(strike, maturity, type, spot, rate, div, volat, lambda, jm, jv, num_terms = 200)
{
  eta <- exp(jm+jv^2/2)-1
  lam <- lambda*(1+eta)
  n  <- 0:num_terms
  sigma_n <- sqrt(volat^2+(n*jv^2)/maturity)
  rtilde <- rate-div-lambda*eta-0.5*volat^2+n*jm/maturity
  r_n <- rtilde+sigma_n^2/2
  bs_n <- mapply(function(X, Y){
    analyticBlackScholesPrice(strike, maturity, type, spot, Y, div = 0, X)
  }, X = sigma_n, Y = r_n)
  pp <- stats::dpois(n, lambda = lam*maturity)
  price <- sum(bs_n*pp)
  return(price)
}

#' Merton's jump-diffusion model analytic price formula for European options
#'

#' @param strike The strike price
#' @param maturity the year until expirations
#' @param type What type of option, put or call
#' @param spot the spot price
#' @param rate the drift rate
#' @param div The dividend yield rate
#' @param volat the diffusion coefficients of each component (a vector)
#' @param lambda the yearly jump-rate
#' @param jm the jump mean
#' @param jv the jump volatility
#' @param n number of variates
#'
#' @description {Monte-Carlo approximation to European price under Merton's jump diffusion.}
#' @return numeric
#' @export monteCarloMertonPrice
monteCarloMertonPrice <- function(strike, maturity, type, spot, rate, div, volat, lambda, jm, jv, n = 10^4)
{
  eta <- exp(jm+0.5*jv^2)-1
  variates <- matrix(0, nrow = n)
  for(i in 1:n)
  {
    numjumps <- stats::rpois(1, lambda = lambda*maturity)
    jumpsizes <- stats::rnorm(numjumps, mean = jm, sd = jv)
    jumps <- sum(jumpsizes)
    variates[i] <- spot*exp((rate-div-lambda*eta-0.5*volat^2)*maturity+volat*sqrt(maturity)*stats::rnorm(1)+jumps)
  }
  if(type == "call")
  {
    variates <- pmax(variates-strike, 0.0)
  } else if(type == "put")
  {
    variates <- pmax(strike - variates, 0.0)
  } else if(type == "straddle")
  {
    variates <- pmax(strike-variates, 0.0)+pmax(variates-strike, 0.0)
  } else if(type == "cdf")
  {
    variates <- ifelse(variates <= strike, 1, 0)
  }
  if(type != "cdf")
  {
    variates <- exp(-(rate-div)*maturity)*variates
  }
  estimate <- mean(variates)
  return(estimate)
}
