
# family ESS
gswLf <- function(w, wL,
                  a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02,
                  pe=-1.58*10^-3, b=4.38, h2=l*LAI/nZ*p/1000, kxmax=5, c=2.64, d=3.54){
  
  # ps(w)
  psf <- function(w)pe*w^(-b)
  # xylem conductance function
  kxf <- function(x)kxmax*exp(-(-x/d)^c)
  
  ps <- psf(w)
  pxmin <- psf(wL)
  kxmin <- kxf(pxmin)
  res <- (ps-pxmin)*h2*kxmin/(h*VPD)
  return(res)
}

# Af(gs)
Af <- function(gs, ca=400, Vcmax=50, cp=30, Km=703, Rd=1, LAI=1)LAI*1/2*(Vcmax+(Km+ca)*gs-Rd-((Vcmax)^2+2*Vcmax*(Km-ca+2*cp)*gs+((ca+Km)*gs+Rd)^2-2*Rd*Vcmax)^(1/2))

# averA for invader
averAif <- function(wLi, wLr,
                    a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02,
                    pe=-1.58*10^-3, b=4.38, h2=h/1000, kxmax=5, c=2.64, d=3.54,
                    gamma=1/((MAP/365/k)/1000)*nZ){
  
  gswLfr <- Vectorize(function(w)ifelse(w<wLr, 0, gswLf(w, wLr)))
  gswLfi <- Vectorize(function(w)ifelse(w<wLi, 0, gswLf(w, wLi)))
  
  Evf <- function(w)h*VPD*gswLfr(w)
  Lf <- function(w)Evf(w)+w/200
  rLf <- function(w)1/Lf(w)
  integralrLf <- Vectorize(function(w)integrate(rLf, w, 1, rel.tol=.Machine$double.eps^0.3)$value)
  fnoc <- function(w)1/Lf(w)*exp(-gamma*w-k*integralrLf(w))
  
  f1 <- Vectorize(function(w)Af(gswLfi(w))*fnoc(w))
  res <- integrate(f1, wLi, 1, rel.tol=.Machine$double.eps^0.3)$value
  return(res)
}

optwLf <- Vectorize(function(wLr){
  averAif1 <- Vectorize(function(wLi)averAif(wLi, wLr))
  res <- optimize(averAif1, c(0.1, 0.3), tol=.Machine$double.eps, maximum=T)$maximum
  message(wLr, res)
  return(res)
})
