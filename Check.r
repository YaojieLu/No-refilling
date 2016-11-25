
source("Functions.r")
data <- read.csv("Results/ESS - no refilling 1.csv")

ca <- 400
k <- 0.025
MAP <- 3650

wLr <- as.numeric(subset(data, k==0.025 & MAP==3650, select=optwL))

f1 <- Vectorize(function(wLi)averAirelf(wLi, wLr=wLr))

curve(f1, wLr-1e-4, wLr+1e-4)
abline(v=wLr)
