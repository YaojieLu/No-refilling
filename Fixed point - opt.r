
options(digits=20)
source("Functions.r")

#environmental conditions
ca <- 400
k <- c(0.025, 0.05, 0.1)
MAP <- seq(0.5, 10, by=0.5)*365
env <- as.vector(expand.grid(ca, k, MAP))

# Initialize
optwL <- matrix(nrow=nrow(env), ncol=1)

# Sensitivity Analysis
for(i in 1:nrow(env)){
  ca <- env[i, 1]
  k <- env[i, 2]
  MAP <- env[i, 3]
  optwL[i, ] <- optimize(optwLf, c(0.14, 0.2), tol=.Machine$double.eps)$minimum
}

# Collect results
res <- cbind(env, optwL)
colnames(res) <- c("ca", "k", "MAP", "optwL")

write.csv(res, "Results/ESS - no refilling.csv", row.names = FALSE)
