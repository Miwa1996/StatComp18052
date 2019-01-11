#' @title A Metropolis-Hastings sampler using R
#' @description A Metropolis-Hastings sampler using R
#' @param a Mean Value of Normal Distribution
#' @param b Variance of Normal Distribution
#' @param N Test number
#' @return Histogram of Sampling data
#' @examples
#' \dontrun{
#' Norm(3,5,10000)
#' }
#' @export
Norm = function(a,b,N){

  x <- c()
  x0 <-100
  x[1] <- x0
  k <- 0
  u <- runif(N)
  for (i in 2:N) {
    y <- rnorm(1, x[i - 1],sqrt(b))
    if (u[i] <(dnorm(y,a,b)/dnorm(x[i - 1],a,b)))
      x[i] <- y else {
        x[i] <- x[i - 1]
        k <- k + 1
      }
  }

  hist(x[500:N],freq = F)
  curve(dnorm(x,a,b),add=TRUE)

}
