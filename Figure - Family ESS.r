
source("Functions.r")

SA <- c(0.2, 0.19, 0.18, 0.15, 0.14)

# Figure
Cols <- rainbow(length(SA))
windows(8, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, mar=c(3.5, 4, 1, 1), mfrow=c(1, 1))
plot(0, 0, xlim=c(0, 1), ylim=c(0, 0.3), cex.lab=1.5,
     type="n", xaxt="n", yaxt="n", xlab=NA, ylab=NA)

axis(1, xlim=c(0, 1), pos=0, lwd=2)
mtext(expression(italic(w)),side=1, line=2, cex=1.3)
axis(2, ylim=c(0, 0.3), pos=0, lwd=2)
mtext(expression(italic(g[s])~(mol~m^-2~s^-1)), side=2, line=1.9, cex=1.3)

# Family ESS
for(i in 1:length(SA)){
  with(list(wL=SA[i], Col=Cols[i]), {
    gswLf1 <- Vectorize(function(w)gswLf(w, wL))
    curve(gswLf1, wL, 1, add=T, col=Col)
  })
}

legend("bottomright", title=expression(italic(w[L])), as.character(SA), lty=1, col=Cols, bg="white")
box()
