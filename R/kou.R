#' Comptue a conditional expectation via Monte-Carlo simulation/averaging
#'

#' @param strike the strike price
#' @param maturity the maturity of the option contract
#' @param type the payoff to use
#' @param spot the spot price
#' @param rate the risk-neutral rate of return
#' @param div the continuous dividend yield
#' @param volat the volatility level
#' @param lambda the mean number of jumps per annum
#' @param parameters vector of model parameters
#' @param n number of variates to sample
#'
#' @description {Compute a conditional expectation via discounting and averaging}
#' @export monteCarloKouPrice
monteCarloKouPrice <- function(strike, maturity, type = "call", spot, rate, div, volat, lambda, parameters, n = 20000)
{
  # rate <- parameters[1]
  # div <- parameters[2]
  # volat <- parameters[3]
  # lambda <- parameters[4]
  prob <- parameters[1]
  alpha <- parameters[2]
  beta <- parameters[3]
  ku <- parameters[4]
  kd <- parameters[5]

  eta <- kellyfractions:::mgfdkou(1, prob, alpha, beta, ku, kd)-1
  variates <- matrix(0, nrow = n)
  for(i in 1:n)
  {
    numjumps <- stats::rpois(1, lambda = lambda*maturity)
    jumpsizes <- findistr::rdkou(numjumps, prob, alpha, beta, ku, kd)
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
