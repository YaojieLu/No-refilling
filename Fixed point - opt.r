
options(digits=20)
source("Functions.r")

# Results
ca <- 400
k <- 0.05
MAP <- 1825
pkx <- 0.5

f1 <- Vectorize(function(wLr){
  averAif1 <- Vectorize(function(wLi)averAif(wLi, wLr))
  res1 <- averAif1(wLr)
  res2 <- optimize(averAif1, c(0.1, 0.3), tol=.Machine$double.eps, maximum=T)
  res <- (res2$objective-res1)^2
  message(wLr, res2$maximum)
  return(res)
})

optwL <- optimize(f1, c(0.14, 0.2), tol=.Machine$double.eps^0.25)
optwL