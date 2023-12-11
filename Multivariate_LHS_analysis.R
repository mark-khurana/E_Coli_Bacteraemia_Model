#Multivariate LHS analysis

#Load libraries
library(deSolve)
library(tidyverse)
library('lhs')
library('ggplot2')
library('sensitivity')
library(patchwork)
library(cowplot)

#Model code
new.basic.model <- function(t, state, parms) { 
  with(as.list(c(state, parms)),{ 
    #Parameters
    t <- parms["t"]
    a <- parms["a"]
    a1 <- parms["a1"]
    a2 <- parms["a2"]
    a3 <- parms["a3"]
    a4 <- parms["a4"]
    a5 <- parms["a5"]
    a6 <- parms["a6"]
    a7 <- parms["a7"]
    d <- parms["d"]
    D <- parms["D"]
    L <-parms["L"]
    D1 <- parms["D1"]
    D2 <- parms["D2"]
    D3 <- parms["D3"]
    NS <- parms["NS"]
    e <- parms["e"]
    e1 <- parms["e1"]
    e2 <- parms["e2"]
    e3 <- parms["e3"]
    e4 <- parms["e4"]
    e5 <- parms["e5"]
    e6 <- parms["e6"]
    e7 <- parms["e7"]
    e8 <- parms["e8"]
    gamma <- parms["gamma"]
    R1 <- parms["R1"]
    R2 <- parms["R2"]
    R3 <- parms["R3"]
    R4 <- parms["R4"]
    R5 <- parms["R5"]
    R6 <- parms["R6"]
    SA <- parms["SA"]
    SB <- parms["SB"]
    SC <- parms["SC"]
    
    #State Basics
    WA <- state["WA"]
    WB <- state["WB"]
    WC <- state["WC"]
    TA <- state["TA"]
    TB <- state["TB"]
    TC <- state["TC"]
    NE <- state["NE"]
    ND <- state["NE"]
    NR <- state["NE"]
    TA <- state["TA"]
    TB <- state["TB"]
    TC <- state["TC"]
    h <- state["h"]
    
    
    #Tree Diagram states
    AES <- state["AES"]
    AESD <- state["AESD"]
    AESO  <- state["AESO"]
    AESOR <- state["AESOR"]
    AESOD <- state["AESOD"]
    AESOEF <- state["AESOEF"]
    AESOES <- state["AESOES"]
    AESOEFE <- state["AESOEFE"]
    AESOEFO <- state["AESOEFO"]
    AESOEFD <- state["AESOEFD"]
    AESOESD <- state["AESOESD"]
    AESOESO <- state["AESOESO"]
    AESOEFED <- state["AESOEFED"]
    AESOEFEO <- state["AESOEFEO"]
    AESOEFOR <- state["AESOEFOR"]
    AESOEFOD <- state["AESOEFOD"]
    AESOEFOE <- state["AESOEFOE"]
    AESOESOR <- state["AESOESOR"]
    AESOESOD <- state["AESOESOD"]
    AESOEFEOR <- state["AESOEFEOR"]
    AESOEFEOD <- state["AESOEFEOD"]
    AESOEFOED <- state["AESOEFOED"]
    AESOEFOEO <- state["AESOEFOEO"]
    AESOEFOEOR <- state["AESOEFOEOR"]
    AESOEFOEOD <- state["AESOEFOEOD"]
    
    AESEF  <- state["AESEF"]
    AESEFE <- state["AESEFE"]
    AESEFO <- state["AESEFO"]
    AESEFD <- state["AESEFD"]
    AESEFED <- state["AESEFED"]
    AESEFEO <- state["AESEFEO"]
    AESEFOR <- state["AESEFOR"]
    AESEFOD <- state["AESEFOD"]
    AESEFOE <- state["AESEFOE"]
    AESEFOEO <- state["AESEFOEO"]
    AESEFEOR <- state["AESEFEOR"]
    AESEFEOD <- state["AESEFEOD"]
    AESEFOED <- state["AESEFOED"]
    AESEFOEOR <- state["AESEFOEOR"]
    AESEFOEOD <- state["AESEFOEOD"]
    AESES <- state["AESES"]
    AESESD <- state["AESESD"]
    AESESO <- state["AESESO"]
    AESESOR <- state["AESESOR"]
    AESESOD <- state["AESESOD"]
    
    BEF <- state["BEF"]
    BEFE <- state["BEFE"]
    BEFO <- state["BEFO"]
    BEFD <- state["BEFD"]
    BEFED <- state["BEFED"]
    BEFEO <- state["BEFEO"]
    BEFOR <- state["BEFOR"]
    BEFOD <- state["BEFOD"]
    BEFOE <- state["BEFOE"]
    BEFEOR <- state["BEFEOR"]
    BEFEOD <- state["BEFEOD"]
    BEFOED <- state["BEFOED"]
    BEFOEO <- state["BEFOEO"]
    BEFOEOR <- state["BEFOEOR"]
    BEFOEOD <- state["BEFOEOD"]
    
    CES <- state["CES"]
    CESD <- state["CESD"]
    CESO <- state["CESO"]
    CESOR <- state["CESOR"]
    CESOD <- state["CESOD"]
    
    
    #Equations
    dAES = NS*h*WA - AES*(((e^((1/d)/2))^e1)*D*D1 + d + a*a1 + a*a2)
    dAESD = AES*((e^((1/d)/2))^e1)*D*D1
    dAESO = AES*d - AESO*(((e^((1/d)/2))^e1)*D*D1 + (1/(t-(1/d)*R1)) + a*a4 + a*a5)
    dAESOR = AESO*(1/(t-(1/d)*R1))
    dAESOD = AESO*((e^((1/d)/2))^e1)*D*D1
    dAESOEF = AESO*a*a4 - AESOEF*(a*a6 + d + ((e^((1/d)/2))^e2)*D*D2)
    dAESOES = AESO*a*a5 - AESOES*(d + ((e^((1/d)/2))^e2)*D*D3)
    dAESOEFE = AESOEF*a*a6 - AESOEFE*(d + ((e^((1/d)/2))^e3)*D*D3)
    dAESOEFO = AESOEF*d - AESOEFO*((1/(t-(1/d)*R2)) + ((e^((1/d)/2))^e2)*D*D2 + a*a7)
    dAESOEFD = AESOEF*((e^((1/d)/2))^e2)*D*D2
    dAESOESD = AESOES*((e^((1/d)/2))^e2)*D*D3
    dAESOESO = AESOES*d - AESOESO*((1/(t-(1/d)*R2)) + ((e^((1/d)/2))^e2)*D*D3)
    dAESOEFED = AESOEFE*((e^((1/d)/2))^e3)*D*D3
    dAESOEFEO = AESOEFE*d - AESOEFEO*((1/(t-(1/d)*R2)) + ((e^((1/d)/2))^e3)*D*D3)
    dAESOEFOR = AESOEFO*(1/(t-(1/d)*R2))
    dAESOEFOD = AESOEFO*((e^((1/d)/2))^e2)*D*D2
    dAESOEFOE = AESOEFO*a*a7 - AESOEFOE*(d + ((e^((1/d)/2))^e3)*D*D3)
    dAESOESOR = AESOESO*(1/(t-(1/d)*R2))
    dAESOESOD = AESOESO*((e^((1/d)/2))^e2)*D*D3
    dAESOEFEOR = AESOEFEO*(1/(t-(1/d)*R2))
    dAESOEFEOD = AESOEFEO*((e^((1/d)/2))^e3)*D*D3
    dAESOEFOED = AESOEFOE*((e^((1/d)/2))^e3)*D*D3
    dAESOEFOEO = AESOEFOE*d - AESOEFOEO*((1/(t-(1/d)*R2)) + ((e^((1/d)/2))^e3)*D*D3)
    dAESOEFOEOR = AESOEFOEO*(1/(t-(1/d)*R2))
    dAESOEFOEOD = AESOEFOEO*((e^((1/d)/2))^e3)*D*D3
    
    dAESEF = AES*a*a1 - AESEF*(a*a6 + d + ((e^((1/d)/2))^e4)*D*D2)
    dAESEFE = AESEF*a*a6 - AESEFE*(((e^((1/d)/2))^e5)*D*D3 + d)
    dAESEFO = AESEF*d - AESEFO*(a*a7 + (1/(t-(1/d)*R3)) + ((e^((1/d)/2))^e4)*D*D2)
    dAESEFD = AESEF*((e^((1/d)/2))^e4)*D*D2
    dAESEFED = AESEFE*((e^((1/d)/2))^e5)*D*D3
    dAESEFEO = AESEFE*d - AESEFEO*((1/(t-(1/d)*R2)) + ((e^((1/d)/2))^e5)*D*D3)
    dAESEFOR = AESEFO*(1/(t-(1/d)*R3))
    dAESEFOD = AESEFO*((e^((1/d)/2))^e4)*D*D2
    dAESEFOE = AESEFO*a*a7 - AESEFOE*(((e^((1/d)/2))^e5)*D*D3 + d)
    dAESEFOEO = AESEFOE*d - AESEFOEO*((1/(t-(1/d)*R2)) + ((e^((1/d)/2))^e5)*D*D3)
    dAESEFEOR = AESEFEO*(1/(t-(1/d)*R2))
    dAESEFEOD = AESEFEO*((e^((1/d)/2))^e5)*D*D3
    dAESEFOED = AESEFOE*((e^((1/d)/2))^e5)*D*D3
    dAESEFOEOR = AESEFOEO*(1/(t-(1/d)*R2))
    dAESEFOEOD = AESEFOEO*((e^((1/d)/2))^e5)*D*D3
    
    dAESES = AES*a*a2 - AESES*(((e^((1/d)/2))^e4)*D*D3 + d)
    dAESESD = AESES*((e^((1/d)/2))^e4)*D*D3
    dAESESO = AESES*d - AESESO*((1/(t-(1/d)*R4)) + ((e^((1/d)/2))^e4)*D*D3)
    dAESESOR = AESESO*(1/(t-(1/d)*R4))
    dAESESOD = AESESO*((e^((1/d)/2))^e4)*D*D3
    
    dBEF = NS*h*WB - BEF*(a*a3 + d + ((e^((1/d)/2))^e6)*D*D2)
    dBEFE = BEF*a*a3 - BEFE*(d + ((e^((1/d)/2))^e7)*D*D3)
    dBEFO = BEF*d - BEFO*((1/(t-(1/d)*R5)) + a*a7 + ((e^((1/d)/2))^e6)*D*D2)
    dBEFD = BEF*((e^((1/d)/2))^e6)*D*D2
    dBEFED = BEFE*((e^((1/d)/2))^e7)*D*D3
    dBEFEO = BEFE*d - BEFEO*((1/(t-(1/d)*R2)) + ((e^((1/d)/2))^e7)*D*D3)
    dBEFOR = BEFO*(1/(t-(1/d)*R5))
    dBEFOD = BEFO*((e^((1/d)/2))^e6)*D*D2
    dBEFOE = BEFO*a*a7 - BEFOE*(d + ((e^((1/d)/2))^e7)*D*D3)
    dBEFEOR = BEFEO*(1/(t-(1/d)*R2))
    dBEFEOD = BEFEO*((e^((1/d)/2))^e7)*D*D3
    dBEFOED = BEFOE*((e^((1/d)/2))^e7)*D*D3
    dBEFOEO = BEFOE*d - BEFOEO*((1/(t-(1/d)*R2)) + ((e^((1/d)/2))^e7)*D*D3)
    dBEFOEOR = BEFOEO*(1/(t-(1/d)*R2))
    dBEFOEOD = BEFOEO*((e^((1/d)/2))^e7)*D*D3
    
    dCES = NS*h*WC - CES*(d + ((e^((1/d)/2))^e8)*D*D3)
    dCESD = CES*((e^((1/d)/2))^e8)*D*D3
    dCESO = CES*d - CESO*((1/(t-(1/d)*R6)) + ((e^((1/d)/2))^e8)*D*D3)
    dCESOR = CESO*(1/(t-(1/d)*R6))
    dCESOD = CESO*((e^((1/d)/2))^e8)*D*D3
    
    #Outcome equations
    dh = h*L
    
    dND = CESO*((e^((1/d)/2))^e8)*D*D3 +  CES*((e^((1/d)/2))^e8)*D*D3 + BEFOEO*((e^((1/d)/2))^e7)*D*D3 + BEFOE*((e^((1/d)/2))^e7)*D*D3 + BEFEO*((e^((1/d)/2))^e7)*D*D3 + BEF*((e^((1/d)/2))^e6)*D*D2 + BEFE*((e^((1/d)/2))^e7)*D*D3 + BEFO*((e^((1/d)/2))^e6)*D*D2 + AESESO*((e^((1/d)/2))^e4)*D*D3 + AESEFE*((e^((1/d)/2))^e5)*D*D3 + AESEFO*((e^((1/d)/2))^e4)*D*D2 + AESEFEO*((e^((1/d)/2))^e5)*D*D3 + AESEFOE*((e^((1/d)/2))^e5)*D*D3 + AESEFOEO*((e^((1/d)/2))^e5)*D*D3 + AESES*((e^((1/d)/2))^e4)*D*D3 + AES*((e^((1/d)/2))^e1)*D*D1 + AESO*((e^((1/d)/2))^e1)*D*D1 + AESOEF*((e^((1/d)/2))^e2)*D*D2 + AESOES*((e^((1/d)/2))^e2)*D*D3 + AESOEFE*((e^((1/d)/2))^e3)*D*D3 + AESOEFO*((e^((1/d)/2))^e2)*D*D2 + AESOESO*((e^((1/d)/2))^e2)*D*D3 + AESOEFEO*((e^((1/d)/2))^e3)*D*D3 + AESOEFOE*((e^((1/d)/2))^e3)*D*D3 + AESOEFOEO*((e^((1/d)/2))^e3)*D*D3
    
    dNR = AESO*(1/(t-(1/d)*R1)) + AESOEFO*(1/(t-(1/d)*R2)) + AESOESO*(1/(t-(1/d)*R2)) + AESOEFEO*(1/(t-(1/d)*R2)) + AESOEFOEO*(1/(t-(1/d)*R2)) + AESEFO*(1/(t-(1/d)*R3)) + AESEFEO*(1/(t-(1/d)*R2)) + AESEFOEO*(1/(t-(1/d)*R2)) + AESESO*(1/(t-(1/d)*R4)) + BEFO*(1/(t-(1/d)*R5)) + BEFEO*(1/(t-(1/d)*R2)) + BEFOEO*(1/(t-(1/d)*R2)) + CESO*(1/(t-(1/d)*R6))
    
    dTA = AES*d - AESO*(((e^((1/d)/2))^e1)*D*D1 + (1/(t-(1/d)*R1)) + a*a4 + a*a5) + AES*a*a1 - AESEF*(a*a6 + d + ((e^((1/d)/2))^e4)*D*D2) + AES*a*a2 - AESES*(((e^((1/d)/2))^e4)*D*D3 + d) + AESO*a*a4 - AESOEF*(a*a6 + d + ((e^((1/d)/2))^e2)*D*D2) + AESO*a*a5 - AESOES*(d + ((e^((1/d)/2))^e2)*D*D3) + SA*(NS*h*WA - AES*(((e^((1/d)/2))^e1)*D*D1 + d +a*a1 +a*a2) + NS*h*WB - BEF*(a*a3 + d + ((e^((1/d)/2))^e6)*D*D2) + NS*h*WC - CES*(d + ((e^((1/d)/2))^e8)*D*D3))
    
    dTB = AESOEF*a*a6 - AESOEFE*(d + ((e^((1/d)/2))^e3)*D*D3) + AESOEF*d - AESOEFO*((1/(t-(1/d)*R2)) + ((e^((1/d)/2))^e2)*D*D2 + a*a7) + AESOEFO*a*a7 - AESOEFOE*(d + ((e^((1/d)/2))^e3)*D*D3) + AESEF*a*a6 - AESEFE*(((e^((1/d)/2))^e5)*D*D3 + d) + AESEF*d - AESEFO*(a*a7 + (1/(t-(1/d)*R3)) + ((e^((1/d)/2))^e4)*D*D2) + AESEFO*a*a7 - AESEFOE*(((e^((1/d)/2))^e5)*D*D3 + d) + BEF*a*a3 - BEFE*(d + ((e^((1/d)/2))^e7)*D*D3) + BEF*d - BEFO*((1/(t-(1/d)*R5)) +a*a7 + ((e^((1/d)/2))^e6)*D*D2) + BEFO*a*a7 - BEFOE*(d + ((e^((1/d)/2))^e7)*D*D3) + SB*(NS*h*WA - AES*(((e^((1/d)/2))^e1)*D*D1 + d + a*a1 + a*a2) + NS*h*WB - BEF*(a*a3 + d + ((e^((1/d)/2))^e6)*D*D2) + NS*h*WC - CES*(d + ((e^((1/d)/2))^e8)*D*D3))
    
    dTC = AESOEFE*d - AESOEFEO*((1/(t-(1/d)*R2)) + ((e^((1/d)/2))^e3)*D*D3) + AESOEFOE*d - AESOEFOEO*((1/(t-(1/d)*R2)) + ((e^((1/d)/2))^e3)*D*D3) + AESOES*d - AESOESO*((1/(t-(1/d)*R2)) + ((e^((1/d)/2))^e2)*D*D3) + AESEFE*d - AESEFEO*((1/(t-(1/d)*R2)) + ((e^((1/d)/2))^e5)*D*D3) + AESEFOE*d - AESEFOEO*((1/(t-(1/d)*R2)) + ((e^((1/d)/2))^e5)*D*D3) + AESES*d - AESESO*((1/(t-(1/d)*R4)) + ((e^((1/d)/2))^e4)*D*D3) + BEFE*d - BEFEO*((1/(t-(1/d)*R2)) + ((e^((1/d)/2))^e7)*D*D3) + BEFOE*d - BEFOEO*((1/(t-(1/d)*R2)) + ((e^((1/d)/2))^e7)*D*D3) + CES*d - CESO*((1/(t-(1/d)*R6)) + ((e^((1/d)/2))^e8)*D*D3) + SC*(NS*h*WA - AES*(((e^((1/d)/2))^e1)*D*D1 + d +a*a1 +a*a2) + NS*h*WB - BEF*(a*a3 + d + ((e^((1/d)/2))^e6)*D*D2) + NS*h*WC - CES*(d + ((e^((1/d)/2))^e8)*D*D3))
    
    dNE = NS*h*WA - AES*(((e^((1/d)/2))^e1)*D*D1 + d +a*a1 +a*a2) + NS*h*WB - BEF*(a*a3 + d + ((e^((1/d)/2))^e6)*D*D2) + NS*h*WC - CES*(d + ((e^((1/d)/2))^e8)*D*D3)
    
    dWA = -1*gamma*TA - (1/3)*gamma*TB - (1/3)*gamma*TC
    dWB = 1*gamma*TA - (1/3)*gamma*TB - (2/3)*gamma*TC
    dWC = (2/3)*gamma*TB + 1*gamma*TC
    
    dTAA = 1*TA
    dTBB = 1*TB
    dTCC = 1*TC
    
    dxdt <- c(dh, dTAA, dTBB, dTCC, dND, dNR, dTA, dTB, dTC, dNE, dWA, dWB, dWC, dAES, dAESD, dAESO, dAESOR, dAESOD, dAESOEF, dAESOES, dAESOEFE, dAESOEFO, dAESOEFD, dAESOESD, dAESOESO, dAESOEFED, dAESOEFEO, dAESOEFO, dAESOEFOD, dAESOEFOE, dAESOESOR, dAESOESOD, dAESOEFEOR, dAESOEFEOD, dAESOEFOED, dAESOEFOEO, dAESOEFOEOR, dAESOEFOEOD, dAESEF, dAESEFE, dAESEFO, dAESEFD, dAESEFED, dAESEFEO, dAESEFOR, dAESEFOD, dAESEFOE, dAESEFOEO, dAESEFEOR, dAESEFEOD, dAESEFOED, dAESEFOEOR, dAESEFOEOD, dAESES, dAESESD, dAESESO, dAESESOR, dAESESOD, dBEF, dBEFE, dBEFO, dBEFD, dBEFED, dBEFEO, dBEFOR, dBEFOD,dBEFOE, dBEFEOR, dBEFEOD, dBEFOED, dBEFOEO, dBEFOEOR, dBEFOEOD, dCES, dCESD, dCESO, dCESOR, dCESOD)
    ## return result as a list
    list(dxdt)}
  )
}

#LHS Code For Second Line Resistance -----------

param_ranges <- as.data.frame(cbind(c(5,14),c(0.0005,0.01), c(0.25,0.99), c(0.005,0.02), 
                                    c(1.1,2), c(2,4), c(0.000001, 0.0001), c(1.1, 5)))
colnames(param_ranges) <- c("min","max")
get.size <- function(out) tail(1,out$WC)
nsamples <- 1000
nparameters <- 8
X <- randomLHS(nsamples, nparameters)

Y <- matrix(0, nrow=nsamples, ncol=35)


Y[,1] <- qunif(X[,1], min = param_ranges[1,1], max = param_ranges[2,1])
Y[,2] <- qunif(X[,2], min = param_ranges[1,2], max = param_ranges[2,2])
Y[,3] <- qunif(X[,3], min = param_ranges[1,3], max = param_ranges[2,3])
Y[,4] <- qunif(X[,4], min = param_ranges[1,4], max = param_ranges[2,4]) 
Y[,5] <- qunif(X[,5], min = param_ranges[1,5], max = param_ranges[2,5])
Y[,6] <- qunif(X[,6], min = param_ranges[1,6], max = param_ranges[2,6]) 
Y[,7] <- qunif(X[,7], min = param_ranges[1,7], max = param_ranges[2,7])
Y[,8] <- qunif(X[,8], min = param_ranges[1,8], max = param_ranges[2,8])
Y[,9] <- 0.0000466
Y[,10] <- 1.375
Y[,11] <- 1.125
Y[,12] <- 1.5
Y[,13] <- 1
Y[,14] <- 0.2
Y[,15] <- 1
Y[,16] <- 1.5
Y[,17] <- 1
Y[,18] <- 100000
Y[,19] <- 0
Y[,20] <- 0
Y[,21] <- 0
Y[,22] <- 0
Y[,23] <- 0
Y[,24] <- 0
Y[,25] <- 0
Y[,26] <- 0
Y[,27] <- 1
Y[,28] <- 0
Y[,29] <- 1
Y[,30] <- 1
Y[,31] <- 1
Y[,32] <- 1
Y[,33] <- 0
Y[,34] <- 0
Y[,35] <- 1

colnames(Y) <- c("t","a", "d", "D", "D2", "D3", "gamma", "e",
                 "L", "a1", "a2", "a3", "a4", "a5", "a6", "a7", "D1", "NS", "e1",
                 "e2", "e3", "e4", "e5", "e6", "e7", "e8", "R1", 
                 "R2", "R3", "R4", "R5", "R6", "SA", "SB", "SC")

store_out_treat <- matrix(NA,nsamples, 2) # where store output
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)

for(i in 1:nsamples){
  out <- as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=Y[i,])) # could also store Y[i,] in store output to save going back to LHS matrix later
  store_out_treat[i,] <- c(i,out[74,"WC"])
}
store_out_treat <- as.data.frame(store_out_treat)
colnames(store_out_treat) <- c("lhs_sample","output")


#Graph
data_mean_treat <- mean(store_out_treat$output)
data_max_treat <- max(store_out_treat$output)
data_min_treat <- min(store_out_treat$output)

lhs_plot_treatment_duration <- ggplot(store_out_treat, aes(x = "Sample", y = output*100)) +
  geom_violin(fill = "black", color = "black", alpha = 0.5) +
  labs(
    x = "Density",
    y = "LHS Output",
    title = "Second-line resistance"
  )

#Combine parameter values and output values
mergedlhs <- cbind(store_out_treat, Y)

#Partial Rank Correlation
bonferroni.alpha <- 0.05/12
prcc_treatment <- pcc(mergedlhs[,3:10], mergedlhs[,2], nboot = 1000, rank=TRUE, conf=1-bonferroni.alpha)

summary_treatment <- print(prcc_treatment)

# Create a data frame for the summary data
summary_df <- data.frame(
  Parameter = row.names(summary_treatment),
  Coefficient = summary_treatment$original,
  Lower_CI = summary_treatment$`min. c.i.`,
  Upper_CI = summary_treatment$`max. c.i.`
)

#plotting
coefficient_plot_treatment_duration <- ggplot(summary_df, aes(x = Parameter, y = Coefficient)) +
  geom_point(size = 2, color = "black") +
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Add horizontal line
  labs(
    x = "Parameter",
    y = "Coefficient",
    title = "Partial rank correlation coefficients, second line resistance"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Display the plot
print(coefficient_plot_treatment_duration)







# LHS Death -----

#Scenario A, Baseline
parms <- c(L=0.0000466, t=7, a=0.001, a1=1, a2=0.2, a3=1, a4=1, a5=0.2, a6=1, a7=1.5, d=0.5, D=0.01, D1=1, D2=1.5, D3=2, NS=100000, e=1.3, e1=0, e2=0, e3=2, e4=1, e5=2, e6=1, e7=2, e8=1, gamma=0.000018, R1=1, R2=0, R3=0, R4=0, R5=0, R6=0, SA=1, SB=0, SC=0)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputBaselineA<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))


#Add new column
outputBaselineA$cum_mort_rate <- outputBaselineA$ND/(outputBaselineA$NR + outputBaselineA$ND)*100

outputBaselineA$Scenario <- "A"

#LHS Code with Treatment Duration
param_ranges <- as.data.frame(cbind(c(5,14),c(0.0005,0.01), c(0.25,0.99), c(0.005,0.02), 
                                    c(1.1,2), c(2,4), c(0.000001, 0.0001), c(1.1, 5)))
colnames(param_ranges) <- c("min","max")
get.size <- function(out) tail(1,out$ND)
nsamples <- 1000
nparameters <- 8
X <- randomLHS(nsamples, nparameters)

Y <- matrix(0, nrow=nsamples, ncol=35)

Y[,1] <- qunif(X[,1], min = param_ranges[1,1], max = param_ranges[2,1])
Y[,2] <- qunif(X[,2], min = param_ranges[1,2], max = param_ranges[2,2])
Y[,3] <- qunif(X[,3], min = param_ranges[1,3], max = param_ranges[2,3])
Y[,4] <- qunif(X[,4], min = param_ranges[1,4], max = param_ranges[2,4]) 
Y[,5] <- qunif(X[,5], min = param_ranges[1,5], max = param_ranges[2,5])
Y[,6] <- qunif(X[,6], min = param_ranges[1,6], max = param_ranges[2,6]) 
Y[,7] <- qunif(X[,7], min = param_ranges[1,7], max = param_ranges[2,7])
Y[,8] <- qunif(X[,8], min = param_ranges[1,8], max = param_ranges[2,8])
Y[,9] <- 0.0000466
Y[,10] <- 1.375
Y[,11] <- 1.125
Y[,12] <- 1.5
Y[,13] <- 1
Y[,14] <- 0.2
Y[,15] <- 1
Y[,16] <- 1.5
Y[,17] <- 1
Y[,18] <- 100000
Y[,19] <- 0
Y[,20] <- 0
Y[,21] <- 0
Y[,22] <- 0
Y[,23] <- 0
Y[,24] <- 0
Y[,25] <- 0
Y[,26] <- 0
Y[,27] <- 1
Y[,28] <- 0
Y[,29] <- 1
Y[,30] <- 1
Y[,31] <- 1
Y[,32] <- 1
Y[,33] <- 0
Y[,34] <- 0
Y[,35] <- 1

colnames(Y) <- c("t","a", "d", "D", "D2", "D3", "gamma", "e",
                 "L", "a1", "a2", "a3", "a4", "a5", "a6", "a7", "D1", "NS", "e1",
                 "e2", "e3", "e4", "e5", "e6", "e7", "e8", "R1", 
                 "R2", "R3", "R4", "R5", "R6", "SA", "SB", "SC")

store_out_death <- matrix(NA,nsamples, 2) # where store output
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)

for(i in 1:nsamples){
  out <- as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=Y[i,])) # could also store Y[i,] in store output to save going back to LHS matrix later
  store_out_death[i,] <- c(i,out[74,"ND"])
}
store_out_death <- as.data.frame(store_out_death)
colnames(store_out_death) <- c("lhs_sample","output")


#Graph
data_mean_death <- mean(store_out_death$output)
data_max_death <- max(store_out_death$output)
data_min_death <- min(store_out_death$output)

lhs_plot_mortality <- ggplot(store_out_death, aes(x = "Sample", y = output)) +
  geom_violin(fill = "black", color = "black", alpha = 0.5) +
  labs(
    x = "Density",
    y = "LHS Output",
    title = "Number of deaths"
  )

#Combine parameter values and output values
mergedlhs <- cbind(store_out_death, Y)


#Partial Rank Correlation
bonferroni.alpha <- 0.05/12
prcc_mortality <- pcc(mergedlhs[,3:10], mergedlhs[,2], nboot = 1000, rank=TRUE, conf=1-bonferroni.alpha)

summary_mortality <- print(prcc_mortality)

# Create a data frame for the summary data
summary_df <- data.frame(
  Parameter = row.names(summary_mortality),
  Coefficient = summary_mortality$original,
  Lower_CI = summary_mortality$`min. c.i.`,
  Upper_CI = summary_mortality$`max. c.i.`
)

#Plotting
coefficient_plot_mortality <- ggplot(summary_df, aes(x = Parameter, y = Coefficient)) +
  geom_point(size = 2, color = "black") +
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Add horizontal line
  labs(
    x = "Parameter",
    y = "Coefficient",
    title = "Partial rank correlation coefficients, Number of deaths"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Display the plot
print(coefficient_plot_mortality)



# LHS Third-line antibiotic use --------

param_ranges <- as.data.frame(cbind(c(5,14),c(0.0005,0.01), c(0.25,0.99), c(0.005,0.02), 
                                    c(1.1,2), c(2,4), c(0.000001, 0.0001), c(1.1, 5)))
colnames(param_ranges) <- c("min","max")
get.size <- function(out) tail(1,out$WC)
nsamples <- 1000
nparameters <- 8
X <- randomLHS(nsamples, nparameters)

Y <- matrix(0, nrow=nsamples, ncol=35)

Y[,1] <- qunif(X[,1], min = param_ranges[1,1], max = param_ranges[2,1])
Y[,2] <- qunif(X[,2], min = param_ranges[1,2], max = param_ranges[2,2])
Y[,3] <- qunif(X[,3], min = param_ranges[1,3], max = param_ranges[2,3])
Y[,4] <- qunif(X[,4], min = param_ranges[1,4], max = param_ranges[2,4]) 
Y[,5] <- qunif(X[,5], min = param_ranges[1,5], max = param_ranges[2,5])
Y[,6] <- qunif(X[,6], min = param_ranges[1,6], max = param_ranges[2,6]) 
Y[,7] <- qunif(X[,7], min = param_ranges[1,7], max = param_ranges[2,7])
Y[,8] <- qunif(X[,8], min = param_ranges[1,8], max = param_ranges[2,8])
Y[,9] <- 0.0000466
Y[,10] <- 1.375
Y[,11] <- 1.125
Y[,12] <- 1.5
Y[,13] <- 1
Y[,14] <- 0.2
Y[,15] <- 1
Y[,16] <- 1.5
Y[,17] <- 1
Y[,18] <- 100000
Y[,19] <- 0
Y[,20] <- 0
Y[,21] <- 0
Y[,22] <- 0
Y[,23] <- 0
Y[,24] <- 0
Y[,25] <- 0
Y[,26] <- 0
Y[,27] <- 1
Y[,28] <- 0
Y[,29] <- 1
Y[,30] <- 1
Y[,31] <- 1
Y[,32] <- 1
Y[,33] <- 0
Y[,34] <- 0
Y[,35] <- 1

colnames(Y) <- c("t","a", "d", "D", "D2", "D3", "gamma", "e",
                 "L", "a1", "a2", "a3", "a4", "a5", "a6", "a7", "D1", "NS", "e1",
                 "e2", "e3", "e4", "e5", "e6", "e7", "e8", "R1", 
                 "R2", "R3", "R4", "R5", "R6", "SA", "SB", "SC")

store_out_3AB <- matrix(NA,nsamples, 2) # where store output
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)

for(i in 1:nsamples){
  out <- as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=Y[i,])) # could also store Y[i,] in store output to save going back to LHS matrix later
  store_out_3AB[i,] <- c(i,out[74,"TCC"])
}
store_out_3AB <- as.data.frame(store_out_3AB)
colnames(store_out_3AB) <- c("lhs_sample","output")


#Graph
data_mean_3AB <- mean(store_out_3AB$output)
data_max_3AB <- max(store_out_3AB$output)
data_min_3AB <- min(store_out_3AB$output)

lhs_plot_3rdline_AB <- ggplot(store_out_3AB, aes(x = "Sample", y = output)) +
  geom_violin(fill = "black", color = "black", alpha = 0.5) +
  labs(
    x = "Density",
    y = "LHS Output",
    title = "Third-line Antibiotic Use"
  )

#Combine parameter values and output values
mergedlhs <- cbind(store_out_3AB, Y)

#Partial Rank Correlation
bonferroni.alpha <- 0.05/12
prcc_3rdlineAB <- pcc(mergedlhs[,3:10], mergedlhs[,2], nboot = 1000, rank=TRUE, conf=1-bonferroni.alpha)

summary_3rdlineAB <- print(prcc_3rdlineAB)

# Create a data frame for the summary data
summary_df <- data.frame(
  Parameter = row.names(summary_3rdlineAB),
  Coefficient = summary_3rdlineAB$original,
  Lower_CI = summary_3rdlineAB$`min. c.i.`,
  Upper_CI = summary_3rdlineAB$`max. c.i.`
)

#Plotting
coefficient_plot_3rdline_AB <- ggplot(summary_df, aes(x = Parameter, y = Coefficient)) +
  geom_point(size = 2, color = "black") +
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Add horizontal line
  labs(
    x = "Parameter",
    y = "Coefficient",
    title = "Partial rank correlation coefficients, 3rd line antibiotic use"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# Combine figures -----
# Combine the individual plots in the top and bottom rows
top_row_plots <- coefficient_plot_treatment_duration + coefficient_plot_mortality + coefficient_plot_3rdline_AB + plot_layout(nrow = 1, widths = c(1/3, 1/3, 1/3))
bottom_row_plots <- lhs_plot_treatment_duration + lhs_plot_mortality + lhs_plot_3rdline_AB + plot_layout(nrow = 1, widths = c(1/3, 1/3, 1/3))

# Create the final patchwork plot with the specified layout
final_plot <- top_row_plots / bottom_row_plots + 
  plot_layout(nrow = 2, heights = c(50, 50))

print(final_plot)




