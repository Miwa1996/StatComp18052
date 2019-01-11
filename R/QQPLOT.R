#' @title A test of normality using R
#' @description A test of normality using R
#' @param alpha Threshold of test
#' @param data A vector of data
#' @return QQ-plot, histogram and test result of shapiro.test
#' @examples
#' \dontrun{
#' data = rnorm(3000)
#' normtest(data)
#' }
#' @export
normtest = function(data,alpha=0.05) {

  if (0 == (n <- length(data)))
    stop("data is empty or has only NAs")

  par(mfrow=c(2,1))
  drop = which(data==0)
  if(length(drop)>0)
    x = sort(data[-drop],decreasing = T) else
      x = sort(data,decreasing = T)
    n = length(x)
    unif.p = -log10(ppoints(n))
    plot(unif.p,x,pch=16,main = 'QQ-plot',xlab = 'Expected',ylab = 'observed',col = 'red')
    lines(c(0,max(unif.p)),c(0,max(unif.p)),lty=2,col = 'blue')

    hist(data,freq=F,main="Histogram and Density Estimation Curve")
    lines(density(data),col="blue")
    x = c(round(min(data)):round(max(data)))
    lines(x,dnorm(x,mean(data),sd(data)),col="red")

    test = shapiro.test(data)
    if(test$p.value>alpha){
      print(paste("success:Subject to normal distribution,p.value=",test$p.value,">",alpha))
    }else{
      print(paste("error:Not subject to normal distribution,p.value=",test$p.value,"<=",alpha))
    }
    test
}










