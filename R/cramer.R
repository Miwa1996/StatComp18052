#' @title function:two-sample Cramer-von Mises test for equal distributions
#' @description function:two-sample Cramer-von Mises test for equal distributions
#' @param x a vector of data
#' @param y a vector of data
#' @param data data sources
#' @return Histogram of  data
#' @examples
#' \dontrun{
#' attach(chickwts)
#' install.packages('latticeExtra')
#' library(latticeExtra)
#' x <- sort(as.vector(weight[feed == "soybean"]))
#' y <- sort(as.vector(weight[feed == "linseed"]))
#' cvm1 <- cvm(x,y,"Example 8.1")
#' }
#' @export



cvm <- function(x,y,data){

  r <- 1000 #permutation samples
  reps <- numeric(r)
  n <- length(x)
  m <- length(y)
  v.n <- numeric(n)
  v1.n <- numeric(n)
  v.m <- numeric(m)
  v1.m <- numeric(m)
  z <- c(x,y)
  N <- length(z)
  Ix <- seq(1:n)
  Iy <- seq(1:m)
  v.n <- (x-Ix)**2
  v.m <- (y-Iy)**2
  #test statistic
  reps_0 <- ((n * sum(v.n)+m * sum(v.m))/(m * n * N))-(4 * m * n - 1)/(6 * N)
  for (k in 1:r){#permutation samples
    w <- sample(N,size=n,replace=FALSE)
    x1 <- sort(z[w])
    y1 <- sort(z[-w])
    v1.n <- (x1-Ix)**2
    v1.m <- (y1-Iy)**2
    reps[k] <- ((n * sum(v1.n)+m * sum(v1.m))/(m * n * N))-(4 * m * n - 1)/(6 * N)
  }
  p <- mean(c(reps_0,reps) >= reps_0)
  return(
    histogram(c(reps_0,reps),
              type="density",
              col="#0080ff",
              xlab="Replicates of Cramer-Von Mises statistic",
              ylab=list(rot=0),
              main=paste0("Data:",data),
              sub=list(substitute(paste(hat(p),"=",pvalue),list(pvalue=p)),col=2),
              panel=function(...){
                panel.histogram(...)
                panel.abline(v=reps_0,col=2,lwd=2)
              })
  )
}
