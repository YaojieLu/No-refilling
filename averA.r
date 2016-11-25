
source("Functions.r")
data <- read.csv("Results/ESS - no refilling.csv")

averA <- numeric(length=nrow(data))

for(i in 1:nrow(data)){
  ca <- data[i, 1]
  k <- data[i, 2]
  MAP <- data[i, 3]
  averA[i] <- averAf(data[i, 4])
  message(k, " ", MAP)
}

# Collect results
res <- cbind(data, averA)
colnames(res) <- c("ca", "k", "MAP", "optwL", "averA")

write.csv(res, "Results/Results - no refilling.csv", row.names = FALSE)
