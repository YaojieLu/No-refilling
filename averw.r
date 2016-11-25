
source("Functions.r")
data <- read.csv("Results/ESS - no refilling.csv")

ca <- 400
b <- 182.5
MAP <- b

# lower k
a <- 0.025
k <- a
wL <- as.numeric(subset(data, k==a & MAP==b, select=optwL))
wL
averwf(wL)

# lower k
a <- 0.1
k <- a
wL <- as.numeric(subset(data, k==a & MAP==b, select=optwL))
wL
averwf(wL)