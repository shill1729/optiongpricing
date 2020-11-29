#' Compare accuracies of Monte-Carlo and analytic formulas for option prices
#'
#' @param strike the strike price of the option
#' @param maturity the maturity in trading years until expiration
#' @param type type of option, "call" or "put"
#' @param spot the current spot price
#' @param rate the risk-neutral rate
#' @param div the continuous dividend yield rate
#' @param volat the volatility level
#' @param a the coefficient of the Poisson process
#' @param b the drift-compensation coefficient
#' @param lambda the mean number of jumps annually
#' @param jm the mean jump-size
#' @param jv the standard deviation of jump-size
#' @param n number of variates to use in Monte-Carlo mean
#'
#' @description {Compare convergence of Monte-Carlo method to analytic formulas for all models.}
#' @return numeric
#' @export testAnalyticVsMonteCarlo
testAnalyticVsMonteCarlo <- function(strike, maturity, type, spot, rate, div, volat, a, b, lambda, jm, jv, n = 10^6)
{
  # Analytic prices
  analyticPrices <- data.frame(bs = 0, pois = 0, merton = 0)
  analyticPrices$bs <- analyticBlackScholesPrice(strike, maturity, type, spot, rate, div, volat)
  analyticPrices$pois <- analyticPoissonPrice(strike, maturity, type, spot, rate, a, b)
  analyticPrices$merton <- analyticMertonPrice(strike, maturity, type, spot, rate, div, volat, lambda, jm, jv)
  # Monte-Carlo prices
  monteCarloPrices <- data.frame(bs = 0, pois = 0)
  monteCarloPrices$bs <- monteCarloBlackScholes(strike, maturity, type, spot, rate, div, volat, n = n)
  monteCarloPrices$pois <- monteCarloPoisson(strike, maturity, type, spot, rate, a, b, n = n)
  monteCarloPrices$merton <- monteCarloMertonPrice(strike, maturity, type, spot, rate, div, volat, lambda, jm, jv, n = 10^4)

  # Combine into one data set and check errors
  prices <- rbind(analyticPrices, monteCarloPrices)
  prices[3, ] <- abs(prices[1, ]-prices[2, ])
  prices[4, ] <- prices[3, ]/prices[1, ]
  prices[5, ] <- prices[4, ]*100
  rownames(prices) <- c("analytic", "monte-carlo", "abs error", "rel. error", "% error")
  return(prices)
}



#' Test Black-Scholes European option pricers
#'
#' @param strike the strike price of the option
#' @param maturity the maturity in trading years until expiration
#' @param type type of option, "call" or "put"
#' @param spot the current spot price
#' @param rate the risk-neutral rate
#' @param div the continuous dividend yield rate
#' @param volat the volatility level
#' @param N time resolution
#' @param M space resolution
#'
#' @description {Verify convergence between the three methods: analytic, Monte-Carlo and PDE}
#' @return matrix
#' @export testBlackScholesEuro
testBlackScholesEuro <- function(strike, maturity, type, spot, rate, div, volat, N = 300, M = 300)
{
  prices <- matrix(0, 3)
  prices[1] <- analyticBlackScholesPrice(strike, maturity, type, spot, rate, div, volat)
  prices[2] <- monteCarloBlackScholes(strike, maturity, type, spot, rate, div, volat)
  prices[3] <- blackScholesPDE(strike, maturity, type, spot, rate, div, volat, c(N, M), FALSE)
  absErrors <- abs(prices-prices[1])
  relErrors <- absErrors/prices[1]
  perErrors <- relErrors*100
  output <- cbind(prices, absErrors, relErrors, perErrors)
  rownames(output) <- c("exact", "monte-carlo", "pde")
  colnames(output) <- c("price", "abs-error", "rel-error", "% error")
  return(output)

}

#' Test Poisson-pricer European option pricers
#'
#' @param strike the strike price of the option
#' @param maturity the maturity in trading years until expiration
#' @param type type of option, "call" or "put"
#' @param spot the current spot price
#' @param rate the risk-neutral rate
#' @param a the Poisson coefficient
#' @param b the drift coefficient
#' @param N time resolution
#' @param M space resolution
#'
#' @description {Verify convergence between the three methods: analytic, Monte-Carlo and PDE}
#' @return matrix
#' @export testPoissonEuro
testPoissonEuro <- function(strike, maturity, type, spot, rate, a, b, N = 300, M = 300)
{
  prices <- matrix(0, 3)
  prices[1] <- analyticPoissonPrice(strike, maturity, type, spot, rate, a, b)
  prices[2] <- monteCarloPoisson(strike, maturity, type, spot, rate, a, b)
  prices[3] <- cpoisPDE(strike, maturity, type, spot, rate, a, b, res = c(N, M), FALSE)
  absErrors <- abs(prices-prices[1])
  relErrors <- absErrors/prices[1]
  perErrors <- relErrors*100
  output <- cbind(prices, absErrors, relErrors, perErrors)
  rownames(output) <- c("exact", "monte-carlo", "pde")
  colnames(output) <- c("price", "abs-error", "rel-error", "% error")
  return(output)

}
