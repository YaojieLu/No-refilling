
source("Functions.r")
data <- read.csv("Results/ESS - no refilling 1.csv")

ca <- 400
a <- 365
MAP <- a

# lower k
b <- 0.025
k <- b
wL <- as.numeric(subset(data, k==b & MAP==a, select=optwL))
cPDF <- cPDFf(wL)
PDFf1 <- function(w)PDFf(w, wL, cPDF)
curve(PDFf1, 1e-3, 1)
abline(v=wL)

# higher k
b <- 0.1
k <- b
wL <- as.numeric(subset(data, k==b & MAP==a, select=optwL))
cPDF <- cPDFf(wL)
PDFf1 <- function(w)PDFf(w, wL, cPDF)
curve(PDFf1, 1e-3, 1, add=T, col="red")
abline(v=wL, col="red")
