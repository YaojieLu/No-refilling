
source("Functions.r")
source("Functions (S1).r")
data <- read.csv("Results/ESS - no refilling 1.csv")

# Initializing
ca <- 400
ESSf1 <- Vectorize(function(w)ESSf(w))
ESSBf1 <- Vectorize(function(w)ESSBf(w))
wL <- uniroot(ESSBf1, c(0.12, 1), tol=.Machine$double.eps)$root

# Figures
Cols <- c("purple", "forestgreen", "orange")
Ltys <- rep(seq(1, length(Cols), by=1), length(Cols))
windows(8, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, mar=c(3.5, 4, 2, 1), mfrow=c(1, 1))
# S1
curve(ESSf1, wL, 1,
      xaxt="n", yaxt="n", xlab=NA, ylab=NA,
      xlim=c(0, 1), ylim=c(0, 0.3),
      cex.lab=1.5, lwd=4, main="ESS in Scenario 2")

segments(0, 0, wL, 0, lty=2, lwd=4)
curve(gsmaxf, 0, 1, add=T, lty=4, lwd=4)

axis(1, xlim=c(0, 1), pos=0, lwd=2)
mtext(expression(italic(w)),side=1, line=2, cex=1.3)
axis(2, ylim=c(0, 0.3), pos=0, lwd=2)
mtext(expression(italic(g[s])~(mol~m^-2~s^-1)), side=2, line=1.9, cex=1.3)

# ESS
for(i in 1:nrow(data)){
  with(list(wL=data$optwL[i], Col=Cols[ceiling(i/(length(Cols)))], Lty=Ltys[i]), {
    gswLf1 <- Vectorize(function(w)gswLf(w, wL))
    curve(gswLf1, wL, 1, add=T, col=Col, lty=Lty)
  })
}

legend("topleft", title="Scenario 1", c("ESS", "Water use envelope"), lty=c(1, 4), lwd=4)
legend("bottomleft", title="MAP (mm)", c("365", "1825", "3650"), pch=15, col=Cols, bg="white")
legend("bottomright", title="Rainfall frequency (per day)", c("0.025", "0.05", "0.1"), lty=Ltys, bg="white")
box()

dev.copy2pdf(file = "Figures/ESS.pdf")
