# Scenario Analysis

#Load packages:
library("cowplot")
library(deSolve)
library(tidyverse)
library(ggplot2)
library(patchwork)

#Basic model
new.basic.model <- function(t, state, parms) { 
  with(as.list(c(state, parms)),{ 
    #Parameters
    t <- parms["t"]
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
    dAES = NS*h*WA - AES*(e1*D1 + d + a1 + a2)
    dAESD = AES*e1*D1
    dAESO = AES*d - AESO*(e1*D1 + R1 + a4 + a5)
    dAESOR = AESO*R1
    dAESOD = AESO*e1*D1
    dAESOEF = AESO*a4 - AESOEF*(a6 + d + e2*D2)
    dAESOES = AESO*a5 - AESOES*(d + e2*D3)
    dAESOEFE = AESOEF*a6 - AESOEFE*(d + e3*D3)
    dAESOEFO = AESOEF*d - AESOEFO*(R2 + e2*D2 + a7)
    dAESOEFD = AESOEF*e2*D2
    dAESOESD = AESOES*e2*D3
    dAESOESO = AESOES*d - AESOESO*(R2 + e2*D3)
    dAESOEFED = AESOEFE*e3*D3
    dAESOEFEO = AESOEFE*d - AESOEFEO*(R2 + e3*D3)
    dAESOEFOR = AESOEFO*R2
    dAESOEFOD = AESOEFO*e2*D2
    dAESOEFOE = AESOEFO*a7 - AESOEFOE*(d + e3*D3)
    dAESOESOR = AESOESO*R2
    dAESOESOD = AESOESO*e2*D3
    dAESOEFEOR = AESOEFEO*R2
    dAESOEFEOD = AESOEFEO*e3*D3
    dAESOEFOED = AESOEFOE*e3*D3
    dAESOEFOEO = AESOEFOE*d - AESOEFOEO*(R2 + e3*D3)
    dAESOEFOEOR = AESOEFOEO*R2
    dAESOEFOEOD = AESOEFOEO*e3*D3
    
    dAESEF = AES*a1 - AESEF*(a6 + d + e4*D2)
    dAESEFE = AESEF*a6 - AESEFE*(e5*D3 + d)
    dAESEFO = AESEF*d - AESEFO*(a7 + R3 + e4*D2)
    dAESEFD = AESEF*e4*D2
    dAESEFED = AESEFE*e5*D3
    dAESEFEO = AESEFE*d - AESEFEO*(R2 + e5*D3)
    dAESEFOR = AESEFO*R3
    dAESEFOD = AESEFO*e4*D2
    dAESEFOE = AESEFO*a7 - AESEFOE*(e5*D3 + d)
    dAESEFOEO = AESEFOE*d - AESEFOEO*(R2 + e5*D3)
    dAESEFEOR = AESEFEO*R2
    dAESEFEOD = AESEFEO*e5*D3
    dAESEFOED = AESEFOE*e5*D3
    dAESEFOEOR = AESEFOEO*R2
    dAESEFOEOD = AESEFOEO*e5*D3
    
    dAESES = AES*a2 - AESES*(e4*D3 + d)
    dAESESD = AESES*e4*D3
    dAESESO = AESES*d - AESESO*(R4 + e4*D3)
    dAESESOR = AESESO*R4
    dAESESOD = AESESO*e4*D3
    
    dBEF = NS*h*WB - BEF*(a3 + d + e6*D2)
    dBEFE = BEF*a3 - BEFE*(d + e7*D3)
    dBEFO = BEF*d - BEFO*(R5 + a7 + e6*D2)
    dBEFD = BEF*e6*D2
    dBEFED = BEFE*e7*D3
    dBEFEO = BEFE*d - BEFEO*(R2 + e7*D3)
    dBEFOR = BEFO*R5
    dBEFOD = BEFO*e6*D2
    dBEFOE = BEFO*a7 - BEFOE*(d + e7*D3)
    dBEFEOR = BEFEO*R2
    dBEFEOD = BEFEO*e7*D3
    dBEFOED = BEFOE*e7*D3
    dBEFOEO = BEFOE*d - BEFOEO*(R2 + e7*D3)
    dBEFOEOR = BEFOEO*R2
    dBEFOEOD = BEFOEO*e7*D3
    
    dCES = NS*h*WC - CES*(d + e8*D3)
    dCESD = CES*e8*D3
    dCESO = CES*d - CESO*(R6 + e8*D3)
    dCESOR = CESO*R6
    dCESOD = CESO*e8*D3
    
    #Outcome equations
    dh = h*L
    
    dND = CESO*e8*D3 +  CES*e8*D3 + BEFOEO*e7*D3 + BEFOE*e7*D3 + BEFEO*e7*D3 + BEF*e6*D2 + BEFE*e7*D3 + BEFO*e6*D2 + AESESO*e4*D3 + AESEFE*e5*D3 + AESEFO*e4*D2 + AESEFEO*e5*D3 + AESEFOE*e5*D3 + AESEFOEO*e5*D3 + AESES*e4*D3 + AES*e1*D1 + AESO*e1*D1 + AESOEF*e2*D2 + AESOES*e2*D3 + AESOEFE*e3*D3 + AESOEFO*e2*D2 + AESOESO*e2*D3 + AESOEFEO*e3*D3 + AESOEFOE*e3*D3 + AESOEFOEO*e3*D3
    
    dNR = AESO*R1 + AESOEFO*R2 + AESOESO*R2 + AESOEFEO*R2 + AESOEFOEO*R2 + AESEFO*R3 + AESEFEO*R2 + AESEFOEO*R2 + AESESO*R4 + BEFO*R5 + BEFEO*R2 + BEFOEO*R2 + CESO*R6
    
    dTA = AES*d - AESO*(e1*D1 + R1 + a4 + a5) + AES*a1 - AESEF*(a6 + d + e4*D2) + AES*a2 - AESES*(e4*D3 + d) + AESO*a4 - AESOEF*(a6 + d + e2*D2) + AESO*a5 - AESOES*(d + e2*D3) + SA*(NS*h*WA - AES*(e1*D1 + d +a1 +a2) + NS*h*WB - BEF*(a3 + d + e6*D2) + NS*h*WC - CES*(d + e8*D3))
    
    dTB = AESOEF*a6 - AESOEFE*(d + e3*D3) + AESOEF*d - AESOEFO*(R2 + e2*D2 + a7) + AESOEFO*a7 - AESOEFOE*(d + e3*D3) + AESEF*a6 - AESEFE*(e5*D3 + d) + AESEF*d - AESEFO*(a7 + R3 + e4*D2) + AESEFO*a7 - AESEFOE*(e5*D3 + d) + BEF*a3 - BEFE*(d + e7*D3) + BEF*d - BEFO*(R5 +a7 + e6*D2) + BEFO*a7 - BEFOE*(d + e7*D3) + SB*(NS*h*WA - AES*(e1*D1 + d + a1 + a2) + NS*h*WB - BEF*(a3 + d + e6*D2) + NS*h*WC - CES*(d + e8*D3))
    
    dTC = AESOEFE*d - AESOEFEO*(R2 + e3*D3) + AESOEFOE*d - AESOEFOEO*(R2 + e3*D3) + AESOES*d - AESOESO*(R2 + e2*D3) + AESEFE*d - AESEFEO*(R2 + e5*D3) + AESEFOE*d - AESEFOEO*(R2 + e5*D3) + AESES*d - AESESO*(R4 + e4*D3) + BEFE*d - BEFEO*(R2 + e7*D3) + BEFOE*d - BEFOEO*(R2 + e7*D3) + CES*d - CESO*(R6 + e8*D3) + SC*(NS*h*WA - AES*(e1*D1 + d +a1 +a2) + NS*h*WB - BEF*(a3 + d + e6*D2) + NS*h*WC - CES*(d + e8*D3))
    
    dNE = NS*h*WA - AES*(e1*D1 + d +a1 +a2) + NS*h*WB - BEF*(a3 + d + e6*D2) + NS*h*WC - CES*(d + e8*D3)
    
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


# IEAT ---------------

#Different Scenarios
#Low, Scenario A
parms <- c(L=0.0000466, t=7, a1=0.001, a2=0.001/5, a3=0.001, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=0.5, D=0.01, D1=0.01, D2=0.01*1.5, D3=0.01*2, NS=100000, e=1.3^((1/0.5)/2), e1=(1.3^((1/0.5)/2))^0, e2=(1.3^((1/0.5)/2))^0, e3=(1.3^((1/0.5)/2))^2, e4=(1.3^((1/0.5)/2)), e5=(1.3^((1/0.5)/2))^2, e6=(1.3^((1/0.5)/2)), e7=(1.3^((1/0.5)/2))^2, e8=(1.3^((1/0.5)/2)), gamma=0.000018, R1=(1/(7-(1/0.5))), R2=1/7, R3=(1/7), R4=(1/7), R5=(1/7), R6=(1/7), SA=1, SB=0, SC=0)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenALow<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenALow$cum_mort_rate <- outputScenALow$ND/(outputScenALow$NR + outputScenALow$ND)*100

outputScenALow$Scenario <- "A"
outputScenALow$ScenAnalysis <- "Low"


#Graphs
options(scipen=999)  # turn off scientific notation like 1e+06

#Low, Scenario B
parms <- c(L=0.0000466, t=7, a1=0.001*1.375, a2=0.001*1.125, a3=0.001*1.5, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=0.5, D=0.01, D1=0.01, D2=0.01*1.5, D3=0.01*2, NS=100000, e=1.3^((1/0.5)/2), e1=(1.3^((1/0.5)/2))^0, e2=(1.3^((1/0.5)/2))^0, e3=(1.3^((1/0.5)/2))^1, e4=(1.3^((1/0.5)/2))^0, e5=(1.3^((1/0.5)/2))^1, e6=(1.3^((1/0.5)/2))^0, e7=(1.3^((1/0.5)/2))^1, e8=(1.3^((1/0.5)/2)), gamma=0.000018, R1=(1/(7-(1/0.5))), R2=1/7, R3=(1/(7-(1/0.5))), R4=(1/7), R5=(1/(7-(1/0.5))), R6=(1/7), SA=0, SB=1, SC=0)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenBLow<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenBLow$cum_mort_rate <- outputScenBLow$ND/(outputScenBLow$NR + outputScenBLow$ND)*100

outputScenBLow$Scenario <- "B"
outputScenBLow$ScenAnalysis <- "Low"

#Low, Scenario C
parms <- c(L=0.0000466, t=7, a1=0.001*1.375, a2=0.001*1.125, a3=0.001*1.5, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=0.5, D=0.01, D1=0.01, D2=0.01*1.5, D3=0.01*2, NS=100000, e=1.3^((1/0.5)/2), e1=(1.3^((1/0.5)/2))^0, e2=(1.3^((1/0.5)/2))^0, e3=(1.3^((1/0.5)/2))^0, e4=(1.3^((1/0.5)/2))^0, e5=(1.3^((1/0.5)/2))^0, e6=(1.3^((1/0.5)/2))^0, e7=(1.3^((1/0.5)/2))^0, e8=(1.3^((1/0.5)/2))^0, gamma=0.000018, R1=(1/(7-(1/0.5))), R2=1/7, R3=(1/(7-(1/0.5))), R4=(1/(7-(1/0.5))), R5=(1/(7-(1/0.5))), R6=(1/(7-(1/0.5))), SA=0, SB=0, SC=1)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenCLow<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenCLow$cum_mort_rate <- outputScenCLow$ND/(outputScenCLow$NR + outputScenCLow$ND)*100

outputScenCLow$Scenario <- "C"
outputScenCLow$ScenAnalysis <- "Low"

##Medium, 0.001
#Medium, Scenario A
parms <- c(L=0.0000466, t=7, a1=0.001, a2=0.001/5, a3=0.001, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=0.5, D=0.01, D1=0.01, D2=0.01*1.5, D3=0.01*2, NS=100000, e=2^((1/0.5)/2), e1=(2^((1/0.5)/2))^0, e2=(2^((1/0.5)/2))^0, e3=(2^((1/0.5)/2))^2, e4=(2^((1/0.5)/2)), e5=(2^((1/0.5)/2))^2, e6=(2^((1/0.5)/2)), e7=(2^((1/0.5)/2))^2, e8=(2^((1/0.5)/2)), gamma=0.000018, R1=(1/(7-(1/0.5))), R2=1/7, R3=(1/7), R4=(1/7), R5=(1/7), R6=(1/7), SA=1, SB=0, SC=0)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenAMed<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenAMed$cum_mort_rate <- outputScenAMed$ND/(outputScenAMed$NR + outputScenAMed$ND)*100

outputScenAMed$Scenario <- "A"
outputScenAMed$ScenAnalysis <- "Medium"


#Graphs
options(scipen=999)  # turn off scientific notation like 1e+06

#Medium, Scenario B
parms <- c(L=0.0000466, t=7, a1=0.001*1.375, a2=0.001*1.125, a3=0.001*1.5, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=0.5, D=0.01, D1=0.01, D2=0.01*1.5, D3=0.01*2, NS=100000, e=2^((1/0.5)/2), e1=(2^((1/0.5)/2))^0, e2=(2^((1/0.5)/2))^0, e3=(2^((1/0.5)/2))^1, e4=(2^((1/0.5)/2))^0, e5=(2^((1/0.5)/2))^1, e6=(2^((1/0.5)/2))^0, e7=(2^((1/0.5)/2))^1, e8=(2^((1/0.5)/2)), gamma=0.000018, R1=(1/(7-(1/0.5))), R2=1/7, R3=(1/(7-(1/0.5))), R4=(1/7), R5=(1/(7-(1/0.5))), R6=(1/7), SA=0, SB=1, SC=0)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenBMed<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenBMed$cum_mort_rate <- outputScenBMed$ND/(outputScenBMed$NR + outputScenBMed$ND)*100

outputScenBMed$Scenario <- "B"
outputScenBMed$ScenAnalysis <- "Medium"

#Medium, Scenario C
parms <- c(L=0.0000466, t=7, a1=0.001*1.375, a2=0.001*1.125, a3=0.001*1.5, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=0.5, D=0.01, D1=0.01, D2=0.01*1.5, D3=0.01*2, NS=100000, e=2^((1/0.5)/2), e1=(2^((1/0.5)/2))^0, e2=(2^((1/0.5)/2))^0, e3=(2^((1/0.5)/2))^0, e4=(2^((1/0.5)/2))^0, e5=(2^((1/0.5)/2))^0, e6=(2^((1/0.5)/2))^0, e7=(2^((1/0.5)/2))^0, e8=(2^((1/0.5)/2))^0, gamma=0.000018, R1=(1/(7-(1/0.5))), R2=1/7, R3=(1/(7-(1/0.5))), R4=(1/(7-(1/0.5))), R5=(1/(7-(1/0.5))), R6=(1/(7-(1/0.5))), SA=0, SB=0, SC=1)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenCMed<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenCMed$cum_mort_rate <- outputScenCMed$ND/(outputScenCMed$NR + outputScenCMed$ND)*100

outputScenCMed$Scenario <- "C"
outputScenCMed$ScenAnalysis <- "Medium"

##High, 0.005

#High, Scenario A
parms <- c(L=0.0000466, t=7, a1=0.001, a2=0.001/5, a3=0.001, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=0.5, D=0.01, D1=0.01, D2=0.01*1.5, D3=0.01*2, NS=100000, e=5^((1/0.5)/2), e1=(5^((1/0.5)/2))^0, e2=(5^((1/0.5)/2))^0, e3=(5^((1/0.5)/2))^2, e4=(5^((1/0.5)/2)), e5=(5^((1/0.5)/2))^2, e6=(5^((1/0.5)/2)), e7=(5^((1/0.5)/2))^2, e8=(5^((1/0.5)/2)), gamma=0.000018, R1=(1/(7-(1/0.5))), R2=1/7, R3=(1/7), R4=(1/7), R5=(1/7), R6=(1/7), SA=1, SB=0, SC=0)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenAHigh<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenAHigh$cum_mort_rate <- outputScenAHigh$ND/(outputScenAHigh$NR + outputScenAHigh$ND)*100
outputScenAHigh$Scenario <- "A"
outputScenAHigh$ScenAnalysis <- "High"

#High, Scenario B
parms <- c(L=0.0000466, t=7, a1=0.001*1.375, a2=0.001*1.125, a3=0.001*1.5, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=0.5, D=0.01, D1=0.01, D2=0.01*1.5, D3=0.01*2, NS=100000, e=5^((1/0.5)/2), e1=(5^((1/0.5)/2))^0, e2=(5^((1/0.5)/2))^0, e3=(5^((1/0.5)/2))^1, e4=(5^((1/0.5)/2))^0, e5=(5^((1/0.5)/2))^1, e6=(5^((1/0.5)/2))^0, e7=(5^((1/0.5)/2))^1, e8=(5^((1/0.5)/2)), gamma=0.000018, R1=(1/(7-(1/0.5))), R2=1/7, R3=(1/(7-(1/0.5))), R4=(1/7), R5=(1/(7-(1/0.5))), R6=(1/7), SA=0, SB=1, SC=0)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenBHigh<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenBHigh$cum_mort_rate <- outputScenBHigh$ND/(outputScenBHigh$NR + outputScenBHigh$ND)*100
outputScenBHigh$Scenario <- "B"
outputScenBHigh$ScenAnalysis <- "High"

#Medium, Scenario C
parms <- c(L=0.0000466, t=7, a1=0.001*1.375, a2=0.001*1.125, a3=0.001*1.5, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=0.5, D=0.01, D1=0.01, D2=0.01*1.5, D3=0.01*2, NS=100000, e=5^((1/0.5)/2), e1=(5^((1/0.5)/2))^0, e2=(5^((1/0.5)/2))^0, e3=(5^((1/0.5)/2))^0, e4=(5^((1/0.5)/2))^0, e5=(5^((1/0.5)/2))^0, e6=(5^((1/0.5)/2))^0, e7=(5^((1/0.5)/2))^0, e8=(5^((1/0.5)/2))^0, gamma=0.000018, R1=(1/(7-(1/0.5))), R2=1/7, R3=(1/(7-(1/0.5))), R4=(1/(7-(1/0.5))), R5=(1/(7-(1/0.5))), R6=(1/(7-(1/0.5))), SA=0, SB=0, SC=1)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenCHigh<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenCHigh$cum_mort_rate <- outputScenCHigh$ND/(outputScenCHigh$NR + outputScenCHigh$ND)*100

outputScenCHigh$Scenario <- "C"
outputScenCHigh$ScenAnalysis <- "High"



#Merging datasets
DatasetScenIEAT <- rbind(outputScenALow, outputScenBLow, outputScenCLow, outputScenAMed, 
                         outputScenBMed, outputScenCMed, outputScenAHigh, outputScenBHigh, outputScenCHigh)


#Mortality
ND <- subset(DatasetScenIEAT, time==1825)
ND$ScenAnalysis <- factor(ND$ScenAnalysis,levels = c("High", "Medium", "Low"))

NDMortality <- ggplot(ND, aes(ScenAnalysis, ND)) +  geom_point(aes(shape = Scenario, colour = Scenario), size = 1.5, stroke = 1, position = position_jitter(0.05)) +
  scale_color_brewer(palette="Set2") +
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
  labs(x="",
       y="# Deaths") + theme_classic()+
  theme(legend.position = "none", axis.text.x=element_blank())+
  coord_cartesian(ylim=c(0, 50)) 
NDMortality

#Mortality Rate
MR <- ggplot(ND, aes(ScenAnalysis, cum_mort_rate)) + geom_point(aes(shape = Scenario, colour = Scenario), size = 1.5, stroke = 1, position = position_jitter(0.05)) +
  scale_color_brewer(palette="Set2") +
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
  labs(x="",
       y="CFR, %") + theme_classic() +
  theme(legend.position = "none", axis.text.x=element_blank())+
  coord_cartesian(ylim=c(0, 15))
MR

#Resistance to 1st Line Therapy
ResistanceLineOne <- ggplot(ND, aes(ScenAnalysis, WB*100)) + geom_point(aes(shape = Scenario, colour = Scenario), size = 1.5, stroke = 1, position = position_jitter(0.05)) +
  scale_color_brewer(palette="Set2") +
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
  labs(x="",
       y="% Res., 1.L.AB.") + theme_classic()+ 
  theme(legend.position = "none", axis.text.x=element_blank())+
  coord_cartesian(ylim=c(0, 8))
ResistanceLineOne

#Resistance to 2nd Line Therapy
ResistanceLineTwo <- ggplot(ND, aes(ScenAnalysis, WC*100)) + geom_point(aes(shape = Scenario, colour = Scenario), size = 1.5, stroke = 1, position = position_jitter(0.05)) +
  scale_color_brewer(palette="Set2") +
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
  labs(x="",
       y="% Res., 2.L.AB.") + theme_classic()+ 
  theme(legend.position = "none", axis.text.x=element_blank()) +
  coord_cartesian(ylim=c(0, 8))
ResistanceLineTwo


#Number Days Treatment, 1st Line Antibiotics
options(scipen=999)
DaysFirstLine <- ggplot(ND, aes(ScenAnalysis, TAA)) + geom_point(aes(shape = Scenario, colour = Scenario), size = 1.5, stroke = 1, position = position_jitter(0.05)) +
  scale_color_brewer(palette="Set2") +
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
  labs(x="",
       y="Days, F.L.AB.") + theme_classic()+ 
  theme(legend.position = "none", axis.text.x=element_blank()) +
  coord_cartesian(ylim=c(0, 2600))
DaysFirstLine

#Number Days Treatment, 2nd Line Antibiotics
DaysSecondLine <- ggplot(ND, aes(ScenAnalysis, TBB)) + geom_point(aes(shape = Scenario, colour = Scenario), size = 1.5, stroke = 1, position = position_jitter(0.05)) +
  scale_color_brewer(palette="Set2") +
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
  labs(x="",
       y="Days, S.L.AB.") + theme_classic()+ 
  theme(legend.position = "none", axis.text.x=element_blank()) +
  coord_cartesian(ylim=c(0, 1000))
DaysSecondLine

#Number Days Treatment, 3rd Line Antibiotics
DaysThirdLine <- ggplot(ND, aes(ScenAnalysis, TCC)) + 
  geom_point(aes(shape = Scenario, colour = Scenario), size = 1.5, stroke = 1, position = position_jitter(0.05)) +
  scale_color_brewer(palette = "Set2") +
  scale_x_discrete(labels = c("High" = "High (ε=5)", "Medium" = "Medium (ε=2)", "Low" = "Low (ε=1.3)")) +
  theme(axis.text.x = element_text(angle = 65, vjust = 0.6)) +
  labs(x = "IEAT Mortality Group",
       y = "Days, 3.L.AB.") +
  theme_classic() +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0, 1000))


# Plot--
ScenIEATCombined <- ggdraw() +
  draw_plot(NDMortality, x = 0.0, y = 0.81, width = 0.99, height = 0.15) +
  draw_plot(MR, x = 0.0, y = 0.68, width = 0.99, height = 0.15) +
  draw_plot(ResistanceLineOne, x = 0.0, y = 0.55, width = 0.99, height = 0.15) +
  draw_plot(ResistanceLineTwo, x = 0.0, y = 0.41, width = 0.99, height = 0.15) +
  draw_plot(DaysFirstLine, x = 0.0, y = 0.28, width = 0.99, height = 0.15) +
  draw_plot(DaysSecondLine, x = 0.0, y = 0.15, width = 0.99, height = 0.15) +
  draw_plot(DaysThirdLine, x = 0.0, y = 0.01, width = 0.99, height = 0.16)





ScenIEATCombined <- ScenIEATCombined +
  draw_label(
    "A) Inappropriate Empiric Antibiotic Therapy (IEAT), ε", fontface = 'bold', x = 0.01, y = .98, hjust = 0)
ScenIEATCombined




# Resistance ---------

#Different Scenarios
#Low, Scenario A
parms <- c(L=0.0000466, t=7, a1=0.0005, a2=0.0005/5, a3=0.0005, a4=0.0005, a5=0.0005/5, a6=0.0005, a7=0.0005*1.5, d=0.5, D=0.01, D1=0.01, D2=0.01*1.5, D3=0.01*2, NS=100000, e=1.3^((1/0.5)/2), e1=(1.3^((1/0.5)/2))^0, e2=(1.3^((1/0.5)/2))^0, e3=(1.3^((1/0.5)/2))^2, e4=(1.3^((1/0.5)/2)), e5=(1.3^((1/0.5)/2))^2, e6=(1.3^((1/0.5)/2)), e7=(1.3^((1/0.5)/2))^2, e8=(1.3^((1/0.5)/2)), gamma=0.000018, R1=(1/(7-(1/0.5))), R2=1/7, R3=(1/7), R4=(1/7), R5=(1/7), R6=(1/7), SA=1, SB=0, SC=0)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenALow<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenALow$cum_mort_rate <- outputScenALow$ND/(outputScenALow$NR + outputScenALow$ND)*100

outputScenALow$Scenario <- "A"
outputScenALow$ScenAnalysis <- "Low"


#Graphs
options(scipen=999)  # turn off scientific notation like 1e+06

#Low, Scenario B
parms <- c(L=0.0000466, t=7, a1=0.0005*1.375, a2=0.0005*1.125, a3=0.0005*1.5, a4=0.0005, a5=0.0005/5, a6=0.0005, a7=0.0005*1.5, d=0.5, D=0.01, D1=0.01, D2=0.01*1.5, D3=0.01*2, NS=100000, e=1.3^((1/0.5)/2), e1=(1.3^((1/0.5)/2))^0, e2=(1.3^((1/0.5)/2))^0, e3=(1.3^((1/0.5)/2))^1, e4=(1.3^((1/0.5)/2))^0, e5=(1.3^((1/0.5)/2))^1, e6=(1.3^((1/0.5)/2))^0, e7=(1.3^((1/0.5)/2))^1, e8=(1.3^((1/0.5)/2)), gamma=0.000018, R1=(1/(7-(1/0.5))), R2=1/7, R3=(1/(7-(1/0.5))), R4=(1/7), R5=(1/(7-(1/0.5))), R6=(1/7), SA=0, SB=1, SC=0)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenBLow<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenBLow$cum_mort_rate <- outputScenBLow$ND/(outputScenBLow$NR + outputScenBLow$ND)*100

outputScenBLow$Scenario <- "B"
outputScenBLow$ScenAnalysis <- "Low"

#Low, Scenario C
parms <- c(L=0.0000466, t=7, a1=0.0005*1.375, a2=0.0005*1.125, a3=0.0005*1.5, a4=0.0005, a5=0.0005/5, a6=0.0005, a7=0.0005*1.5, d=0.5, D=0.01, D1=0.01, D2=0.01*1.5, D3=0.01*2, NS=100000, e=1.3^((1/0.5)/2), e1=(1.3^((1/0.5)/2))^0, e2=(1.3^((1/0.5)/2))^0, e3=(1.3^((1/0.5)/2))^0, e4=(1.3^((1/0.5)/2))^0, e5=(1.3^((1/0.5)/2))^0, e6=(1.3^((1/0.5)/2))^0, e7=(1.3^((1/0.5)/2))^0, e8=(1.3^((1/0.5)/2))^0, gamma=0.000018, R1=(1/(7-(1/0.5))), R2=1/7, R3=(1/(7-(1/0.5))), R4=(1/(7-(1/0.5))), R5=(1/(7-(1/0.5))), R6=(1/(7-(1/0.5))), SA=0, SB=0, SC=1)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenCLow<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenCLow$cum_mort_rate <- outputScenCLow$ND/(outputScenCLow$NR + outputScenCLow$ND)*100

outputScenCLow$Scenario <- "C"
outputScenCLow$ScenAnalysis <- "Low"

##Medium, 0.001
#Medium, Scenario A
parms <- c(L=0.0000466, t=7, a1=0.001, a2=0.001/5, a3=0.001, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=0.5, D=0.01, D1=0.01, D2=0.01*1.5, D3=0.01*2, NS=100000, e=1.3^((1/0.5)/2), e1=(1.3^((1/0.5)/2))^0, e2=(1.3^((1/0.5)/2))^0, e3=(1.3^((1/0.5)/2))^2, e4=(1.3^((1/0.5)/2)), e5=(1.3^((1/0.5)/2))^2, e6=(1.3^((1/0.5)/2)), e7=(1.3^((1/0.5)/2))^2, e8=(1.3^((1/0.5)/2)), gamma=0.000018, R1=(1/(7-(1/0.5))), R2=1/7, R3=(1/7), R4=(1/7), R5=(1/7), R6=(1/7), SA=1, SB=0, SC=0)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenAMed<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenAMed$cum_mort_rate <- outputScenAMed$ND/(outputScenAMed$NR + outputScenAMed$ND)*100

outputScenAMed$Scenario <- "A"
outputScenAMed$ScenAnalysis <- "Medium"


#Graphs
options(scipen=999)  # turn off scientific notation like 1e+06

#Medium, Scenario B
parms <- c(L=0.0000466, t=7, a1=0.001*1.375, a2=0.001*1.125, a3=0.001*1.5, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=0.5, D=0.01, D1=0.01, D2=0.01*1.5, D3=0.01*2, NS=100000, e=1.3^((1/0.5)/2), e1=(1.3^((1/0.5)/2))^0, e2=(1.3^((1/0.5)/2))^0, e3=(1.3^((1/0.5)/2))^1, e4=(1.3^((1/0.5)/2))^0, e5=(1.3^((1/0.5)/2))^1, e6=(1.3^((1/0.5)/2))^0, e7=(1.3^((1/0.5)/2))^1, e8=(1.3^((1/0.5)/2)), gamma=0.000018, R1=(1/(7-(1/0.5))), R2=1/7, R3=(1/(7-(1/0.5))), R4=(1/7), R5=(1/(7-(1/0.5))), R6=(1/7), SA=0, SB=1, SC=0)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenBMed<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenBMed$cum_mort_rate <- outputScenBMed$ND/(outputScenBMed$NR + outputScenBMed$ND)*100

outputScenBMed$Scenario <- "B"
outputScenBMed$ScenAnalysis <- "Medium"

#Medium, Scenario C
parms <- c(L=0.0000466, t=7, a1=0.001*1.375, a2=0.001*1.125, a3=0.001*1.5, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=0.5, D=0.01, D1=0.01, D2=0.01*1.5, D3=0.01*2, NS=100000, e=1.3^((1/0.5)/2), e1=(1.3^((1/0.5)/2))^0, e2=(1.3^((1/0.5)/2))^0, e3=(1.3^((1/0.5)/2))^0, e4=(1.3^((1/0.5)/2))^0, e5=(1.3^((1/0.5)/2))^0, e6=(1.3^((1/0.5)/2))^0, e7=(1.3^((1/0.5)/2))^0, e8=(1.3^((1/0.5)/2))^0, gamma=0.000018, R1=(1/(7-(1/0.5))), R2=1/7, R3=(1/(7-(1/0.5))), R4=(1/(7-(1/0.5))), R5=(1/(7-(1/0.5))), R6=(1/(7-(1/0.5))), SA=0, SB=0, SC=1)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenCMed<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenCMed$cum_mort_rate <- outputScenCMed$ND/(outputScenCMed$NR + outputScenCMed$ND)*100

outputScenCMed$Scenario <- "C"
outputScenCMed$ScenAnalysis <- "Medium"

##High, 0.01

#High, Scenario A
parms <- c(L=0.0000466, t=7, a1=0.01, a2=0.01/5, a3=0.01, a4=0.01, a5=0.01/5, a6=0.01, a7=0.01*1.5, d=0.5, D=0.01, D1=0.01, D2=0.01*1.5, D3=0.01*2, NS=100000, e=1.3^((1/0.5)/2), e1=(1.3^((1/0.5)/2))^0, e2=(1.3^((1/0.5)/2))^0, e3=(1.3^((1/0.5)/2))^2, e4=(1.3^((1/0.5)/2)), e5=(1.3^((1/0.5)/2))^2, e6=(1.3^((1/0.5)/2)), e7=(1.3^((1/0.5)/2))^2, e8=(1.3^((1/0.5)/2)), gamma=0.000018, R1=(1/(7-(1/0.5))), R2=1/7, R3=(1/7), R4=(1/7), R5=(1/7), R6=(1/7), SA=1, SB=0, SC=0)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenAHigh<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenAHigh$cum_mort_rate <- outputScenAHigh$ND/(outputScenAHigh$NR + outputScenAHigh$ND)*100

outputScenAHigh$Scenario <- "A"
outputScenAHigh$ScenAnalysis <- "High"

#High, Scenario B
parms <- c(L=0.0000466, t=7, a1=0.01*1.375, a2=0.01*1.125, a3=0.01*1.5, a4=0.01, a5=0.01/5, a6=0.01, a7=0.01*1.5, d=0.5, D=0.01, D1=0.01, D2=0.01*1.5, D3=0.01*2, NS=100000, e=1.3^((1/0.5)/2), e1=(1.3^((1/0.5)/2))^0, e2=(1.3^((1/0.5)/2))^0, e3=(1.3^((1/0.5)/2))^1, e4=(1.3^((1/0.5)/2))^0, e5=(1.3^((1/0.5)/2))^1, e6=(1.3^((1/0.5)/2))^0, e7=(1.3^((1/0.5)/2))^1, e8=(1.3^((1/0.5)/2)), gamma=0.000018, R1=(1/(7-(1/0.5))), R2=1/7, R3=(1/(7-(1/0.5))), R4=(1/7), R5=(1/(7-(1/0.5))), R6=(1/7), SA=0, SB=1, SC=0)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenBHigh<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenBHigh$cum_mort_rate <- outputScenBHigh$ND/(outputScenBHigh$NR + outputScenBHigh$ND)*100

outputScenBHigh$Scenario <- "B"
outputScenBHigh$ScenAnalysis <- "High"

#High, Scenario C
parms <- c(L=0.0000466, t=7, a1=0.01*1.375, a2=0.01*1.125, a3=0.01*1.5, a4=0.01, a5=0.01/5, a6=0.01, a7=0.01*1.5, d=0.5, D=0.01, D1=0.01, D2=0.01*1.5, D3=0.01*2, NS=100000, e=1.3^((1/0.5)/2), e1=(1.3^((1/0.5)/2))^0, e2=(1.3^((1/0.5)/2))^0, e3=(1.3^((1/0.5)/2))^0, e4=(1.3^((1/0.5)/2))^0, e5=(1.3^((1/0.5)/2))^0, e6=(1.3^((1/0.5)/2))^0, e7=(1.3^((1/0.5)/2))^0, e8=(1.3^((1/0.5)/2))^0, gamma=0.000018, R1=(1/(7-(1/0.5))), R2=1/7, R3=(1/(7-(1/0.5))), R4=(1/(7-(1/0.5))), R5=(1/(7-(1/0.5))), R6=(1/(7-(1/0.5))), SA=0, SB=0, SC=1)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenCHigh<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenCHigh$cum_mort_rate <- outputScenCHigh$ND/(outputScenCHigh$NR + outputScenCHigh$ND)*100

outputScenCHigh$Scenario <- "C"
outputScenCHigh$ScenAnalysis <- "High"



#Merging datasets
DatasetScenAlpha <- rbind(outputScenALow, outputScenBLow, outputScenCLow, outputScenAMed, 
                          outputScenBMed, outputScenCMed, outputScenAHigh, outputScenBHigh, outputScenCHigh)


#Mortality
ND <- subset(DatasetScenAlpha, time==1825)
ND$ScenAnalysis <- factor(ND$ScenAnalysis,levels = c("High", "Medium", "Low"))

NDMortality <- ggplot(ND, aes(ScenAnalysis, ND)) + geom_point(aes(shape = Scenario, colour = Scenario), size = 1.5, stroke = 1, position = position_jitter(0.05)) +
  scale_color_brewer(palette="Set2") +
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
  labs(x="",
       y="") + theme_classic()+
  theme(legend.position = "none", axis.text.x=element_blank())+
  coord_cartesian(ylim=c(0, 50)) 
NDMortality

#Mortality Rate
MR <- ggplot(ND, aes(ScenAnalysis, cum_mort_rate)) + geom_point(aes(shape = Scenario, colour = Scenario), size = 1.5, stroke = 1, position = position_jitter(0.05)) +
  scale_color_brewer(palette="Set2") +
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
  labs(x="",
       y="") + theme_classic() +
  theme(legend.position = "none", axis.text.x=element_blank())+
  coord_cartesian(ylim=c(0, 15))
MR

#Resistance to 1st Line Therapy
ResistanceLineOne <- ggplot(ND, aes(ScenAnalysis, WB*100)) + geom_point(aes(shape = Scenario, colour = Scenario), size = 1.5, stroke = 1, position = position_jitter(0.05)) +
  scale_color_brewer(palette="Set2") +
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
  labs(x="",
       y="") + theme_classic()+ 
  theme(legend.position = "none", axis.text.x=element_blank())+
  coord_cartesian(ylim=c(0, 8))
ResistanceLineOne

#Resistance to 2nd Line Therapy
ResistanceLineTwo <- ggplot(ND, aes(ScenAnalysis, WC*100)) + geom_point(aes(shape = Scenario, colour = Scenario), size = 1.5, stroke = 1, position = position_jitter(0.05)) +
  scale_color_brewer(palette="Set2") +
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
  labs(x="",
       y="") + theme_classic()+ 
  theme(legend.position = "none", axis.text.x=element_blank()) +
  coord_cartesian(ylim=c(0, 8))
ResistanceLineTwo


#Number Days Treatment, 1st Line Antibiotics
options(scipen=999)
DaysFirstLine <- ggplot(ND, aes(ScenAnalysis, TAA)) + geom_point(aes(shape = Scenario, colour = Scenario), size = 1.5, stroke = 1, position = position_jitter(0.05)) +
  scale_color_brewer(palette="Set2") +
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
  labs(x="",
       y="") + theme_classic()+ 
  theme(legend.position = "none", axis.text.x=element_blank()) +
  coord_cartesian(ylim=c(0, 2600))
DaysFirstLine

#Number Days Treatment, 2nd Line Antibiotics
DaysSecondLine <- ggplot(ND, aes(ScenAnalysis, TBB)) + geom_point(aes(shape = Scenario, colour = Scenario), size = 1.5, stroke = 1, position = position_jitter(0.05)) +
  scale_color_brewer(palette="Set2") +
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
  labs(x="",
       y="") + theme_classic()+ 
  theme(legend.position = "none", axis.text.x=element_blank()) +
  coord_cartesian(ylim=c(0, 1000))
DaysSecondLine

#Number Days Treatment, 3rd Line Antibiotics
DaysThirdLine <- ggplot(ND, aes(ScenAnalysis, TCC)) + 
  geom_point(aes(shape = Scenario, colour = Scenario), size = 1.5, stroke = 1, position = position_jitter(0.05)) +
  scale_color_brewer(palette = "Set2") +
  scale_x_discrete(labels = c("High" = "High (a=0.01)", "Medium" = "Medium (a=0.001)", "Low" = "Low (a=0.0005)")) +
  theme(axis.text.x = element_text(angle = 65, vjust = 0.6)) +
  labs(x = "Baseline Resistance Rate",
       y = "") +
  theme_classic() +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0, 1000))

# Plot
ScenResCombined <- ggdraw() +
  draw_plot(NDMortality, x = 0.0, y = 0.81, width = 0.99, height = 0.15) +
  draw_plot(MR, x = 0.0, y = 0.68, width = 0.99, height = 0.15) +
  draw_plot(ResistanceLineOne, x = 0.0, y = 0.55, width = 0.99, height = 0.15) +
  draw_plot(ResistanceLineTwo, x = 0.0, y = 0.41, width = 0.99, height = 0.15) +
  draw_plot(DaysFirstLine, x = 0.0, y = 0.28, width = 0.99, height = 0.15) +
  draw_plot(DaysSecondLine, x = 0.0, y = 0.15, width = 0.99, height = 0.15) +
  draw_plot(DaysThirdLine, x = 0.0, y = 0.01, width = 0.99, height = 0.16)

ScenResCombined <- ScenResCombined +
  draw_label(
    "B) Baseline Resistance Rate, a", fontface = 'bold', x = 0.01, y = .98, hjust = 0)
ScenResCombined


# Treatment and Empiric Duration -------

#Low, Low, Scenario A
parms <- c(L=0.0000466, t=5, a1=0.001, a2=0.001/5, a3=0.001, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=1, D=0.01, D1=0.01, D2=0.01*1.5, D3=0.01*2, NS=100000, e=(1.3^((1/1)/2)), e1=(1.3^((1/1)/2))^0, e2=(1.3^((1/1)/2))^0, e3=(1.3^((1/1)/2))^2, e4=(1.3^((1/1)/2)), e5=(1.3^((1/1)/2))^2, e6=(1.3^((1/1)/2)), e7=(1.3^((1/1)/2))^2, e8=(1.3^((1/1)/2)), gamma=0.000018, R1=(1/(5-(1/1))), R2=1/5, R3=(1/5), R4=(1/5), R5=(1/5), R6=(1/5), SA=1, SB=0, SC=0)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenAEmpShortShort<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))


#Add new columns
outputScenAEmpShortShort$cum_mort_rate <- outputScenAEmpShortShort$ND/(outputScenAEmpShortShort$NR + outputScenAEmpShortShort$ND)*100
outputScenAEmpShortShort$Scenario <- "A"
outputScenAEmpShortShort$ScenAnalysis <- "Low, Low"


#Low, Low, Scenario B
parms <- c(L=0.0000466, t=5, a1=0.001*1.375, a2=0.001*1.125, a3=0.001*1.5, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=1, D=0.01, D1=0.01, D2=0.01*1.5, D3=0.01*2, NS=100000, e=(1.3^((1/1)/2)), e1=(1.3^((1/1)/2))^0, e2=(1.3^((1/1)/2))^0, e3=(1.3^((1/1)/2))^1, e4=(1.3^((1/1)/2))^0, e5=(1.3^((1/1)/2))^1, e6=(1.3^((1/1)/2))^0, e7=(1.3^((1/1)/2))^1, e8=(1.3^((1/1)/2)), gamma=0.000018, R1=(1/(5-(1/1))), R2=1/5, R3=(1/(5-(1/1))), R4=(1/5), R5=(1/(5-(1/1))), R6=(1/5), SA=0, SB=1, SC=0)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenBEmpShortShort<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))


#Add new column
outputScenBEmpShortShort$cum_mort_rate <- outputScenBEmpShortShort$ND/(outputScenBEmpShortShort$NR + outputScenBEmpShortShort$ND)*100
outputScenBEmpShortShort$Scenario <- "B"
outputScenBEmpShortShort$ScenAnalysis <- "Low, Low"

#Low, Low, Scenario C
parms <- c(L=0.0000466, t=5, a1=0.001*1.375, a2=0.001*1.125, a3=0.001*1.5, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=1, D=0.01, D1=0.01, D2=0.01*1.5, D3=0.01*2, NS=100000, e=1.3^((1/1)/2), e1=(1.3^((1/1)/2))^0, e2=(1.3^((1/1)/2))^0, e3=(1.3^((1/1)/2))^0, e4=(1.3^((1/1)/2))^0, e5=(1.3^((1/1)/2))^0, e6=(1.3^((1/1)/2))^0, e7=(1.3^((1/1)/2))^0, e8=(1.3^((1/1)/2))^0, gamma=0.000018, R1=(1/(5-(1/1))), R2=1/5, R3=(1/(5-(1/1))), R4=(1/(5-(1/1))), R5=(1/(5-(1/1))), R6=(1/(5-(1/1))), SA=0, SB=0, SC=1)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenCEmpShortShort<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenCEmpShortShort$cum_mort_rate <- outputScenCEmpShortShort$ND/(outputScenCEmpShortShort$NR + outputScenCEmpShortShort$ND)*100
outputScenCEmpShortShort$Scenario <- "C"
outputScenCEmpShortShort$ScenAnalysis <- "Low, Low"

#Low, Medium, Scenario A
parms <- c(L=0.0000466, t=5, a1=0.001, a2=0.001/5, a3=0.001, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=0.5, D=0.01, D1=0.01, D2=0.01*1.5, D3=0.01*2, NS=100000, e=1.3^((1/0.5)/2), e1=(1.3^((1/0.5)/2))^0, e2=(1.3^((1/0.5)/2))^0, e3=(1.3^((1/0.5)/2))^2, e4=(1.3^((1/0.5)/2)), e5=(1.3^((1/0.5)/2))^2, e6=(1.3^((1/0.5)/2)), e7=(1.3^((1/0.5)/2))^2, e8=(1.3^((1/0.5)/2)), gamma=0.000018, R1=(1/(5-(1/0.5))), R2=1/5, R3=(1/5), R4=(1/5), R5=(1/5), R6=(1/5), SA=1, SB=0, SC=0)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenAEmpShortMed<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenAEmpShortMed$cum_mort_rate <- outputScenAEmpShortMed$ND/(outputScenAEmpShortMed$NR + outputScenAEmpShortMed$ND)*100
outputScenAEmpShortMed$Scenario <- "A"
outputScenAEmpShortMed$ScenAnalysis <- "Low, Medium"

#Low, Medium, Scenario B
parms <- c(L=0.00000188, t=5, a1=0.001*1.375, a2=0.001*1.125, a3=0.001*1.5, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=0.5, D=0.01, D1=0.01, D2=0.01*1.5, D3=0.01*2, NS=100000, e=1.3^((1/0.5)/2), e1=(1.3^((1/0.5)/2))^0, e2=(1.3^((1/0.5)/2))^0, e3=(1.3^((1/0.5)/2))^1, e4=(1.3^((1/0.5)/2))^0, e5=(1.3^((1/0.5)/2))^1, e6=(1.3^((1/0.5)/2))^0, e7=(1.3^((1/0.5)/2))^1, e8=(1.3^((1/0.5)/2)), gamma=0.000018, R1=(1/(5-(1/0.5))), R2=1/5, R3=(1/(5-(1/0.5))), R4=(1/5), R5=(1/(5-(1/0.5))), R6=(1/5), SA=0, SB=1, SC=0)
state <- c(h=0.00000212, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenBEmpShortMed<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenBEmpShortMed$cum_mort_rate <- outputScenBEmpShortMed$ND/(outputScenBEmpShortMed$NR + outputScenBEmpShortMed$ND)*100
outputScenBEmpShortMed$Scenario <- "B"
outputScenBEmpShortMed$ScenAnalysis <- "Low, Medium"

#Low, Medium, Scenario C
parms <- c(L=0.0000466, t=5, a1=0.001*1.375, a2=0.001*1.125, a3=0.001*1.5, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=0.5, D=0.01, D1=0.01, D2=0.01*1.5, D3=0.01*2, NS=100000, e=1.3^((1/0.5)/2), e1=(1.3^((1/0.5)/2))^0, e2=(1.3^((1/0.5)/2))^0, e3=(1.3^((1/0.5)/2))^0, e4=(1.3^((1/0.5)/2))^0, e5=(1.3^((1/0.5)/2))^0, e6=(1.3^((1/0.5)/2))^0, e7=(1.3^((1/0.5)/2))^0, e8=(1.3^((1/0.5)/2))^0, gamma=0.000018, R1=(1/(5-(1/0.5))), R2=1/5, R3=(1/(5-(1/0.5))), R4=(1/(5-(1/0.5))), R5=(1/(5-(1/0.5))), R6=(1/(5-(1/0.5))), SA=0, SB=0, SC=1)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenCEmpShortMed<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))


#Add new column
outputScenCEmpShortMed$cum_mort_rate <- outputScenCEmpShortMed$ND/(outputScenCEmpShortMed$NR + outputScenCEmpShortMed$ND)*100
outputScenCEmpShortMed$Scenario <- "C"
outputScenCEmpShortMed$ScenAnalysis <- "Low, Medium"

#Low, High, Scenario A
parms <- c(L=0.0000466, t=5, a1=0.001, a2=0.001/5, a3=0.001, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=0.25, D=0.01, D1=0.01, D2=0.01*1.5, D3=0.01*2, NS=100000, e=1.3^((1/0.25)/2), e1=(1.3^((1/0.25)/2))^0, e2=(1.3^((1/0.25)/2))^0, e3=(1.3^((1/0.25)/2))^2, e4=(1.3^((1/0.25)/2)), e5=(1.3^((1/0.25)/2))^2, e6=(1.3^((1/0.25)/2)), e7=(1.3^((1/0.25)/2))^2, e8=(1.3^((1/0.25)/2)), gamma=0.000018, R1=(1/(5-(1/0.25))), R2=1/5, R3=(1/5), R4=(1/5), R5=(1/5), R6=(1/5), SA=1, SB=0, SC=0)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenAEmpShortHigh<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenAEmpShortHigh$cum_mort_rate <- outputScenAEmpShortHigh$ND/(outputScenAEmpShortHigh$NR + outputScenAEmpShortHigh$ND)*100
outputScenAEmpShortHigh$Scenario <- "A"
outputScenAEmpShortHigh$ScenAnalysis <- "Low, High"

#Short, High, Scenario B
parms <- c(L=0.0000466, t=5, a1=0.001*1.375, a2=0.001*1.125, a3=0.001*1.5, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=0.25, D=0.01, D1=0.01, D2=0.01*1.5, D3=0.01*2, NS=100000, e=1.3^((1/0.25)/2), e1=(1.3^((1/0.25)/2))^0, e2=(1.3^((1/0.25)/2))^0, e3=(1.3^((1/0.25)/2))^1, e4=(1.3^((1/0.25)/2))^0, e5=(1.3^((1/0.25)/2))^1, e6=(1.3^((1/0.25)/2))^0, e7=(1.3^((1/0.25)/2))^1, e8=(1.3^((1/0.25)/2)), gamma=0.000018, R1=(1/(5-(1/0.25))), R2=1/5, R3=(1/(5-(1/0.25))), R4=(1/5), R5=(1/(5-(1/0.25))), R6=(1/5), SA=0, SB=1, SC=0)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenBEmpShortHigh<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenBEmpShortHigh$cum_mort_rate <- outputScenBEmpShortHigh$ND/(outputScenBEmpShortHigh$NR + outputScenBEmpShortHigh$ND)*100
outputScenBEmpShortHigh$Scenario <- "B"
outputScenBEmpShortHigh$ScenAnalysis <- "Low, High"

#Short, High, Scenario C
parms <- c(L=0.0000466, t=5, a1=0.001*1.375, a2=0.001*1.125, a3=0.001*1.5, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=0.25, D=0.01, D1=0.01, D2=0.01*1.5, D3=0.01*2, NS=100000, e=1.3^((1/0.25)/2), e1=(1.3^((1/0.25)/2))^0, e2=(1.3^((1/0.25)/2))^0, e3=(1.3^((1/0.25)/2))^0, e4=(1.3^((1/0.25)/2))^0, e5=(1.3^((1/0.25)/2))^0, e6=(1.3^((1/0.25)/2))^0, e7=(1.3^((1/0.25)/2))^0, e8=(1.3^((1/0.25)/2))^0, gamma=0.000018, R1=(1/(5-(1/0.25))), R2=1/5, R3=(1/(5-(1/0.25))), R4=(1/(5-(1/0.25))), R5=(1/(5-(1/0.25))), R6=(1/(5-(1/0.25))), SA=0, SB=0, SC=1)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenCEmpShortHigh<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenCEmpShortHigh$cum_mort_rate <- outputScenCEmpShortHigh$ND/(outputScenCEmpShortHigh$NR + outputScenCEmpShortHigh$ND)*100
outputScenCEmpShortHigh$Scenario <- "C"
outputScenCEmpShortHigh$ScenAnalysis <- "Low, High"

##BREAK

#Medium, Low, Scenario A
parms <- c(L=0.0000466, t=7, a1=0.001, a2=0.001/5, a3=0.001, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=1, D=0.01, D1=0.01, D2=0.01*1.5, D3=0.01*2, NS=100000, e=1.3^((1/1)/2), e1=(1.3^((1/1)/2))^0, e2=(1.3^((1/1)/2))^0, e3=(1.3^((1/1)/2))^2, e4=(1.3^((1/1)/2)), e5=(1.3^((1/1)/2))^2, e6=(1.3^((1/1)/2)), e7=(1.3^((1/1)/2))^2, e8=(1.3^((1/1)/2)), gamma=0.000018, R1=(1/(7-(1/1))), R2=1/7, R3=(1/7), R4=(1/7), R5=(1/7), R6=(1/7), SA=1, SB=0, SC=0)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenAEmpMedShort<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenAEmpMedShort$cum_mort_rate <- outputScenAEmpMedShort$ND/(outputScenAEmpMedShort$NR + outputScenAEmpMedShort$ND)*100
outputScenAEmpMedShort$Scenario <- "A"
outputScenAEmpMedShort$ScenAnalysis <- "Medium, Low"

#Medium, Low, Scenario B
parms <- c(L=0.0000466, t=7, a1=0.001*1.375, a2=0.001*1.125, a3=0.001*1.5, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=1, D=0.01, D1=0.01, D2=0.01*1.5, D3=0.01*2, NS=100000, e=1.3^((1/1)/2), e1=(1.3^((1/1)/2))^0, e2=(1.3^((1/1)/2))^0, e3=(1.3^((1/1)/2))^1, e4=(1.3^((1/1)/2))^0, e5=(1.3^((1/1)/2))^1, e6=(1.3^((1/1)/2))^0, e7=(1.3^((1/1)/2))^1, e8=(1.3^((1/1)/2)), gamma=0.000018, R1=(1/(7-(1/1))), R2=1/7, R3=(1/(7-(1/1))), R4=(1/7), R5=(1/(7-(1/1))), R6=(1/7), SA=0, SB=1, SC=0)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenBEmpMedShort<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenBEmpMedShort$cum_mort_rate <- outputScenBEmpMedShort$ND/(outputScenBEmpMedShort$NR + outputScenBEmpMedShort$ND)*100
outputScenBEmpMedShort$Scenario <- "B"
outputScenBEmpMedShort$ScenAnalysis <- "Medium, Low"

#Medium, Low, Scenario C
parms <- c(L=0.0000466, t=7, a1=0.001*1.375, a2=0.001*1.125, a3=0.001*1.5, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=1, D=0.01, D1=0.01, D2=0.01*1.5, D3=0.01*2, NS=100000, e=1.3^((1/1)/2), e1=(1.3^((1/1)/2))^0, e2=(1.3^((1/1)/2))^0, e3=(1.3^((1/1)/2))^0, e4=(1.3^((1/1)/2))^0, e5=(1.3^((1/1)/2))^0, e6=(1.3^((1/1)/2))^0, e7=(1.3^((1/1)/2))^0, e8=(1.3^((1/1)/2))^0, gamma=0.000018, R1=(1/(7-(1/1))), R2=1/7, R3=(1/(7-(1/1))), R4=(1/(7-(1/1))), R5=(1/(7-(1/1))), R6=(1/(7-(1/1))), SA=0, SB=0, SC=1)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenCEmpMedShort<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenCEmpMedShort$cum_mort_rate <- outputScenCEmpMedShort$ND/(outputScenCEmpMedShort$NR + outputScenCEmpMedShort$ND)*100
outputScenCEmpMedShort$Scenario <- "C"
outputScenCEmpMedShort$ScenAnalysis <- "Medium, Low"

#Medium, Medium, Scenario A
parms <- c(L=0.0000466, t=7, a1=0.001, a2=0.001/5, a3=0.001, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=0.5, D=0.01, D1=0.01, D2=0.01*1.5, D3=0.01*2, NS=100000, e=1.3^((1/0.5)/2), e1=(1.3^((1/0.5)/2))^0, e2=(1.3^((1/0.5)/2))^0, e3=(1.3^((1/0.5)/2))^2, e4=(1.3^((1/0.5)/2)), e5=(1.3^((1/0.5)/2))^2, e6=(1.3^((1/0.5)/2)), e7=(1.3^((1/0.5)/2))^2, e8=(1.3^((1/0.5)/2)), gamma=0.000018, R1=(1/(7-(1/0.5))), R2=1/7, R3=(1/7), R4=(1/7), R5=(1/7), R6=(1/7), SA=1, SB=0, SC=0)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenAEmpMedMed<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenAEmpMedMed$cum_mort_rate <- outputScenAEmpMedMed$ND/(outputScenAEmpMedMed$NR + outputScenAEmpMedMed$ND)*100
outputScenAEmpMedMed$Scenario <- "A"
outputScenAEmpMedMed$ScenAnalysis <- "Medium, Medium"

#Medium, Medium, Scenario B
parms <- c(L=0.0000466, t=7, a1=0.001*1.375, a2=0.001*1.125, a3=0.001*1.5, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=0.5, D=0.01, D1=0.01, D2=0.01*1.5, D3=0.01*2, NS=100000, e=1.3^((1/0.5)/2), e1=(1.3^((1/0.5)/2))^0, e2=(1.3^((1/0.5)/2))^0, e3=(1.3^((1/0.5)/2))^1, e4=(1.3^((1/0.5)/2))^0, e5=(1.3^((1/0.5)/2))^1, e6=(1.3^((1/0.5)/2))^0, e7=(1.3^((1/0.5)/2))^1, e8=(1.3^((1/0.5)/2)), gamma=0.000018, R1=(1/(7-(1/0.5))), R2=1/7, R3=(1/(7-(1/0.5))), R4=(1/7), R5=(1/(7-(1/0.5))), R6=(1/7), SA=0, SB=1, SC=0)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenBEmpMedMed<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenBEmpMedMed$cum_mort_rate <- outputScenBEmpMedMed$ND/(outputScenBEmpMedMed$NR + outputScenBEmpMedMed$ND)*100
outputScenBEmpMedMed$Scenario <- "B"
outputScenBEmpMedMed$ScenAnalysis <- "Medium, Medium"

#Medium, Medium, Scenario C
parms <- c(L=0.0000466, t=7, a1=0.001*1.375, a2=0.001*1.125, a3=0.001*1.5, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=0.5, D=0.01, D1=0.01, D2=0.01*1.5, D3=0.01*2, NS=100000, e=1.3^((1/0.5)/2), e1=(1.3^((1/0.5)/2))^0, e2=(1.3^((1/0.5)/2))^0, e3=(1.3^((1/0.5)/2))^0, e4=(1.3^((1/0.5)/2))^0, e5=(1.3^((1/0.5)/2))^0, e6=(1.3^((1/0.5)/2))^0, e7=(1.3^((1/0.5)/2))^0, e8=(1.3^((1/0.5)/2))^0, gamma=0.000018, R1=(1/(7-(1/0.5))), R2=1/7, R3=(1/(7-(1/0.5))), R4=(1/(7-(1/0.5))), R5=(1/(7-(1/0.5))), R6=(1/(7-(1/0.5))), SA=0, SB=0, SC=1)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenCEmpMedMed<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenCEmpMedMed$cum_mort_rate <- outputScenCEmpMedMed$ND/(outputScenCEmpMedMed$NR + outputScenCEmpMedMed$ND)*100
outputScenCEmpMedMed$Scenario <- "C"
outputScenCEmpMedMed$ScenAnalysis <- "Medium, Medium"

#Medium, High, Scenario A
parms <- c(L=0.0000466, t=7, a1=0.001, a2=0.001/5, a3=0.001, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=0.25, D=0.01, D1=0.01, D2=0.01*1.5, D3=0.01*2, NS=100000, e=1.3^((1/0.25)/2), e1=(1.3^((1/0.25)/2))^0, e2=(1.3^((1/0.25)/2))^0, e3=(1.3^((1/0.25)/2))^2, e4=(1.3^((1/0.25)/2)), e5=(1.3^((1/0.25)/2))^2, e6=(1.3^((1/0.25)/2)), e7=(1.3^((1/0.25)/2))^2, e8=(1.3^((1/0.25)/2)), gamma=0.000018, R1=(1/(7-(1/0.25))), R2=1/7, R3=(1/7), R4=(1/7), R5=(1/7), R6=(1/7), SA=1, SB=0, SC=0)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenAEmpMedHigh<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenAEmpMedHigh$cum_mort_rate <- outputScenAEmpMedHigh$ND/(outputScenAEmpMedHigh$NR + outputScenAEmpMedHigh$ND)*100
outputScenAEmpMedHigh$Scenario <- "A"
outputScenAEmpMedHigh$ScenAnalysis <- "Medium, High"

#Medium, High, Scenario B
parms <- c(L=0.0000466, t=7, a1=0.001*1.375, a2=0.001*1.125, a3=0.001*1.5, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=0.25, D=0.01, D1=0.01, D2=0.01*1.5, D3=0.01*2, NS=100000, e=1.3^((1/0.25)/2), e1=(1.3^((1/0.25)/2))^0, e2=(1.3^((1/0.25)/2))^0, e3=(1.3^((1/0.25)/2))^1, e4=(1.3^((1/0.25)/2))^0, e5=(1.3^((1/0.25)/2))^1, e6=(1.3^((1/0.25)/2))^0, e7=(1.3^((1/0.25)/2))^1, e8=(1.3^((1/0.25)/2)), gamma=0.000018, R1=(1/(7-(1/0.25))), R2=1/7, R3=(1/(7-(1/0.25))), R4=(1/7), R5=(1/(7-(1/0.25))), R6=(1/7), SA=0, SB=1, SC=0)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenBEmpMedHigh<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenBEmpMedHigh$cum_mort_rate <- outputScenBEmpMedHigh$ND/(outputScenBEmpMedHigh$NR + outputScenBEmpMedHigh$ND)*100
outputScenBEmpMedHigh$Scenario <- "B"
outputScenBEmpMedHigh$ScenAnalysis <- "Medium, High"

#Medium, High, Scenario C
parms <- c(L=0.0000466, t=7, a1=0.001*1.375, a2=0.001*1.125, a3=0.001*1.5, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=0.25, D=0.01, D1=0.01, D2=0.01*1.5, D3=0.01*2, NS=100000, e=1.3^((1/0.25)/2), e1=(1.3^((1/0.25)/2))^0, e2=(1.3^((1/0.25)/2))^0, e3=(1.3^((1/0.25)/2))^0, e4=(1.3^((1/0.25)/2))^0, e5=(1.3^((1/0.25)/2))^0, e6=(1.3^((1/0.25)/2))^0, e7=(1.3^((1/0.25)/2))^0, e8=(1.3^((1/0.25)/2))^0, gamma=0.000018, R1=(1/(7-(1/0.25))), R2=1/7, R3=(1/(7-(1/0.25))), R4=(1/(7-(1/0.25))), R5=(1/(7-(1/0.25))), R6=(1/(7-(1/0.25))), SA=0, SB=0, SC=1)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenCEmpMedHigh<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenCEmpMedHigh$cum_mort_rate <- outputScenCEmpMedHigh$ND/(outputScenCEmpMedHigh$NR + outputScenCEmpMedHigh$ND)*100
outputScenCEmpMedHigh$Scenario <- "C"
outputScenCEmpMedHigh$ScenAnalysis <- "Medium, High"

##BREAK

#High, Low, Scenario A
parms <- c(L=0.0000466, t=14, a1=0.001, a2=0.001/5, a3=0.001, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=1, D=0.01, D1=0.01, D2=0.01*1.5, D3=0.01*2, NS=100000, e=1.3^((1/1)/2), e1=(1.3^((1/1)/2))^0, e2=(1.3^((1/1)/2))^0, e3=(1.3^((1/1)/2))^2, e4=(1.3^((1/1)/2)), e5=(1.3^((1/1)/2))^2, e6=(1.3^((1/1)/2)), e7=(1.3^((1/1)/2))^2, e8=(1.3^((1/1)/2)), gamma=0.000018, R1=(1/(14-(1/1))), R2=1/14, R3=(1/14), R4=(1/14), R5=(1/14), R6=(1/14), SA=1, SB=0, SC=0)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenAEmpHighShort<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenAEmpHighShort$cum_mort_rate <- outputScenAEmpHighShort$ND/(outputScenAEmpHighShort$NR + outputScenAEmpHighShort$ND)*100
outputScenAEmpHighShort$Scenario <- "A"
outputScenAEmpHighShort$ScenAnalysis <- "High, Low"

#High, short, Scenario B
parms <- c(L=0.0000466, t=14, a1=0.001*1.375, a2=0.001*1.125, a3=0.001*1.5, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=1, D=0.01, D1=0.01, D2=0.01*1.5, D3=0.01*2, NS=100000, e=1.3^((1/1)/2), e1=(1.3^((1/1)/2))^0, e2=(1.3^((1/1)/2))^0, e3=(1.3^((1/1)/2))^1, e4=(1.3^((1/1)/2))^0, e5=(1.3^((1/1)/2))^1, e6=(1.3^((1/1)/2))^0, e7=(1.3^((1/1)/2))^1, e8=(1.3^((1/1)/2)), gamma=0.000018, R1=(1/(14-(1/1))), R2=1/14, R3=(1/(14-(1/1))), R4=(1/14), R5=(1/(14-(1/1))), R6=(1/14), SA=0, SB=1, SC=0)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenBEmpHighShort<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenBEmpHighShort$cum_mort_rate <- outputScenBEmpHighShort$ND/(outputScenBEmpHighShort$NR + outputScenBEmpHighShort$ND)*100
outputScenBEmpHighShort$Scenario <- "B"
outputScenBEmpHighShort$ScenAnalysis <- "High, Low"

#High, Low, Scenario C
parms <- c(L=0.0000466, t=14, a1=0.001*1.375, a2=0.001*1.125, a3=0.001*1.5, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=1, D=0.01, D1=0.01, D2=0.01*1.5, D3=0.01*2, NS=100000, e=1.3^((1/1)/2), e1=(1.3^((1/1)/2))^0, e2=(1.3^((1/1)/2))^0, e3=(1.3^((1/1)/2))^0, e4=(1.3^((1/1)/2))^0, e5=(1.3^((1/1)/2))^0, e6=(1.3^((1/1)/2))^0, e7=(1.3^((1/1)/2))^0, e8=(1.3^((1/1)/2))^0, gamma=0.000018, R1=(1/(14-(1/1))), R2=1/14, R3=(1/(14-(1/1))), R4=(1/(14-(1/1))), R5=(1/(14-(1/1))), R6=(1/(14-(1/1))), SA=0, SB=0, SC=1)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenCEmpHighShort<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenCEmpHighShort$cum_mort_rate <- outputScenCEmpHighShort$ND/(outputScenCEmpHighShort$NR + outputScenCEmpHighShort$ND)*100
outputScenCEmpHighShort$Scenario <- "C"
outputScenCEmpHighShort$ScenAnalysis <- "High, Low"

#High, Medium, Scenario A
parms <- c(L=0.0000466, t=14, a1=0.001, a2=0.001/5, a3=0.001, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=0.5, D=0.01, D1=0.01, D2=0.01*1.5, D3=0.01*2, NS=100000, e=1.3^((1/0.5)/2), e1=(1.3^((1/0.5)/2))^0, e2=(1.3^((1/0.5)/2))^0, e3=(1.3^((1/0.5)/2))^2, e4=(1.3^((1/0.5)/2)), e5=(1.3^((1/0.5)/2))^2, e6=(1.3^((1/0.5)/2)), e7=(1.3^((1/0.5)/2))^2, e8=(1.3^((1/0.5)/2)), gamma=0.000018, R1=(1/(14-(1/0.5))), R2=1/14, R3=(1/14), R4=(1/14), R5=(1/14), R6=(1/14), SA=1, SB=0, SC=0)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenAEmpHighMed<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenAEmpHighMed$cum_mort_rate <- outputScenAEmpHighMed$ND/(outputScenAEmpHighMed$NR + outputScenAEmpHighMed$ND)*100
outputScenAEmpHighMed$Scenario <- "A"
outputScenAEmpHighMed$ScenAnalysis <- "High, Medium"

#High, Medium, Scenario B
parms <- c(L=0.0000466, t=14, a1=0.001*1.375, a2=0.001*1.125, a3=0.001*1.5, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=0.5, D=0.01, D1=0.01, D2=0.01*1.5, D3=0.01*2, NS=100000, e=1.3^((1/0.5)/2), e1=(1.3^((1/0.5)/2))^0, e2=(1.3^((1/0.5)/2))^0, e3=(1.3^((1/0.5)/2))^1, e4=(1.3^((1/0.5)/2))^0, e5=(1.3^((1/0.5)/2))^1, e6=(1.3^((1/0.5)/2))^0, e7=(1.3^((1/0.5)/2))^1, e8=(1.3^((1/0.5)/2)), gamma=0.000018, R1=(1/(14-(1/0.5))), R2=1/14, R3=(1/(14-(1/0.5))), R4=(1/14), R5=(1/(14-(1/0.5))), R6=(1/14), SA=0, SB=1, SC=0)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenBEmpHighMed<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenBEmpHighMed$cum_mort_rate <- outputScenBEmpHighMed$ND/(outputScenBEmpHighMed$NR + outputScenBEmpHighMed$ND)*100
outputScenBEmpHighMed$Scenario <- "B"
outputScenBEmpHighMed$ScenAnalysis <- "High, Medium"

#High, Medium, Scenario C
parms <- c(L=0.0000466, t=14, a1=0.001*1.375, a2=0.001*1.125, a3=0.001*1.5, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=0.5, D=0.01, D1=0.01, D2=0.01*1.5, D3=0.01*2, NS=100000, e=1.3^((1/0.5)/2), e1=(1.3^((1/0.5)/2))^0, e2=(1.3^((1/0.5)/2))^0, e3=(1.3^((1/0.5)/2))^0, e4=(1.3^((1/0.5)/2))^0, e5=(1.3^((1/0.5)/2))^0, e6=(1.3^((1/0.5)/2))^0, e7=(1.3^((1/0.5)/2))^0, e8=(1.3^((1/0.5)/2))^0, gamma=0.000018, R1=(1/(14-(1/0.5))), R2=1/14, R3=(1/(14-(1/0.5))), R4=(1/(14-(1/0.5))), R5=(1/(14-(1/0.5))), R6=(1/(14-(1/0.5))), SA=0, SB=0, SC=1)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenCEmpHighMed<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenCEmpHighMed$cum_mort_rate <- outputScenCEmpHighMed$ND/(outputScenCEmpHighMed$NR + outputScenCEmpHighMed$ND)*100
outputScenCEmpHighMed$Scenario <- "C"
outputScenCEmpHighMed$ScenAnalysis <- "High, Medium"

#High, High, Scenario A
parms <- c(L=0.0000466, t=14, a1=0.001, a2=0.001/5, a3=0.001, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=0.25, D=0.01, D1=0.01, D2=0.01*1.5, D3=0.01*2, NS=100000, e=1.3^((1/0.25)/2), e1=(1.3^((1/0.25)/2))^0, e2=(1.3^((1/0.25)/2))^0, e3=(1.3^((1/0.25)/2))^2, e4=(1.3^((1/0.25)/2)), e5=(1.3^((1/0.25)/2))^2, e6=(1.3^((1/0.25)/2)), e7=(1.3^((1/0.25)/2))^2, e8=(1.3^((1/0.25)/2)), gamma=0.000018, R1=(1/(14-(1/0.25))), R2=1/14, R3=(1/14), R4=(1/14), R5=(1/14), R6=(1/14), SA=1, SB=0, SC=0)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenAEmpHighHigh<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenAEmpHighHigh$cum_mort_rate <- outputScenAEmpHighHigh$ND/(outputScenAEmpHighHigh$NR + outputScenAEmpHighHigh$ND)*100
outputScenAEmpHighHigh$Scenario <- "A"
outputScenAEmpHighHigh$ScenAnalysis <- "High, High"

#High, High, Scenario B
parms <- c(L=0.0000466, t=14, a1=0.001*1.375, a2=0.001*1.125, a3=0.001*1.5, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=0.25, D=0.01, D1=0.01, D2=0.01*1.5, D3=0.01*2, NS=100000, e=1.3^((1/0.25)/2), e1=(1.3^((1/0.25)/2))^0, e2=(1.3^((1/0.25)/2))^0, e3=(1.3^((1/0.25)/2))^1, e4=(1.3^((1/0.25)/2))^0, e5=(1.3^((1/0.25)/2))^1, e6=(1.3^((1/0.25)/2))^0, e7=(1.3^((1/0.25)/2))^1, e8=(1.3^((1/0.25)/2)), gamma=0.000018, R1=(1/(14-(1/0.25))), R2=1/14, R3=(1/(14-(1/0.25))), R4=(1/14), R5=(1/(14-(1/0.25))), R6=(1/14), SA=0, SB=1, SC=0)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenBEmpHighHigh<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenBEmpHighHigh$cum_mort_rate <- outputScenBEmpHighHigh$ND/(outputScenBEmpHighHigh$NR + outputScenBEmpHighHigh$ND)*100
outputScenBEmpHighHigh$Scenario <- "B"
outputScenBEmpHighHigh$ScenAnalysis <- "High, High"

#High, High, Scenario C
parms <- c(L=0.0000466, t=14, a1=0.001*1.375, a2=0.001*1.125, a3=0.001*1.5, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=0.25, D=0.01, D1=0.01, D2=0.01*1.5, D3=0.01*2, NS=100000, e=1.3^((1/0.25)/2), e1=(1.3^((1/0.25)/2))^0, e2=(1.3^((1/0.25)/2))^0, e3=(1.3^((1/0.25)/2))^0, e4=(1.3^((1/0.25)/2))^0, e5=(1.3^((1/0.25)/2))^0, e6=(1.3^((1/0.25)/2))^0, e7=(1.3^((1/0.25)/2))^0, e8=(1.3^((1/0.25)/2))^0, gamma=0.000018, R1=(1/(14-(1/0.25))), R2=1/14, R3=(1/(14-(1/0.25))), R4=(1/(14-(1/0.25))), R5=(1/(14-(1/0.25))), R6=(1/(14-(1/0.25))), SA=0, SB=0, SC=1)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenCEmpHighHigh<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenCEmpHighHigh$cum_mort_rate <- outputScenCEmpHighHigh$ND/(outputScenCEmpHighHigh$NR + outputScenCEmpHighHigh$ND)*100
outputScenCEmpHighHigh$Scenario <- "C"
outputScenCEmpHighHigh$ScenAnalysis <- "High, High"

#Merging Data
DatasetScenEmpTreat <- rbind(outputScenAEmpShortShort, outputScenBEmpShortShort, 
                             outputScenCEmpShortShort, outputScenAEmpMedShort, outputScenBEmpMedShort,outputScenCEmpMedShort,
                             outputScenAEmpHighShort, outputScenBEmpHighShort, 
                             outputScenCEmpHighShort,outputScenAEmpShortMed, 
                             outputScenBEmpShortMed, outputScenCEmpShortMed, 
                             outputScenAEmpMedMed, outputScenBEmpMedMed,
                             outputScenCEmpMedMed,outputScenAEmpHighMed, 
                             outputScenBEmpHighMed, outputScenCEmpHighMed, 
                             outputScenAEmpShortHigh, outputScenBEmpShortHigh, 
                             outputScenCEmpShortHigh,outputScenAEmpMedHigh, 
                             outputScenBEmpMedHigh, outputScenCEmpMedHigh, 
                             outputScenAEmpHighHigh, outputScenBEmpHighHigh, outputScenCEmpHighHigh)


#Graphs
#Mortality
ND <- subset(DatasetScenEmpTreat, time==1825)
ND$ScenAnalysis <- factor(ND$ScenAnalysis,levels = c("High, High", "High, Medium",
                                                     "High, Low","Medium, High",
                                                     "Medium, Medium","Medium, Low",
                                                     "Low, High","Low, Medium",
                                                     "Low, Low"))

NDMortality <- ggplot(ND, aes(ScenAnalysis, ND)) + geom_point(aes(shape = Scenario, colour = Scenario), size = 1.5, stroke = 1, position = position_jitter(0.05)) +
  scale_color_brewer(palette="Set2") +
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
  labs(x="",
       y="") + theme_classic()+
  theme(legend.position = "none", axis.text.x=element_blank())+
  coord_cartesian(ylim=c(0, 60)) 
NDMortality

#Mortality Rate
MR <- ggplot(ND, aes(ScenAnalysis, cum_mort_rate)) + geom_point(aes(shape = Scenario, colour = Scenario), size = 1.5, stroke = 1, position = position_jitter(0.05)) +
  scale_color_brewer(palette="Set2") +
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
  labs(x="",
       y="") + theme_classic() +
  theme(legend.position = "none", axis.text.x=element_blank())+
  coord_cartesian(ylim=c(0, 18))
MR

#Resistance to 1st Line Therapy
ResistanceLineOne <- ggplot(ND, aes(ScenAnalysis, WB*100)) + geom_point(aes(shape = Scenario, colour = Scenario), size = 1.5, stroke = 1, position = position_jitter(0.05)) +
  scale_color_brewer(palette="Set2") +
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
  labs(x="",
       y="") + theme_classic()+ 
  theme(legend.position = "none", axis.text.x=element_blank())+
  coord_cartesian(ylim=c(0, 12))
ResistanceLineOne

#Resistance to 2nd Line Therapy
ResistanceLineTwo <- ggplot(ND, aes(ScenAnalysis, WC*100)) + geom_point(aes(shape = Scenario, colour = Scenario), size = 1.5, stroke = 1, position = position_jitter(0.05)) +
  scale_color_brewer(palette="Set2") +
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
  labs(x="",
       y="") + theme_classic()+ 
  theme(legend.position = "none", axis.text.x=element_blank()) +
  coord_cartesian(ylim=c(0, 8))
ResistanceLineTwo


#Number Days Treatment, 1st Line Antibiotics
options(scipen=999)
DaysFirstLine <- ggplot(ND, aes(ScenAnalysis, TAA)) + geom_point(aes(shape = Scenario, colour = Scenario), size = 1.5, stroke = 1, position = position_jitter(0.05)) +
  scale_color_brewer(palette="Set2") +
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
  labs(x="",
       y="") + theme_classic()+ 
  theme(legend.position = "none", axis.text.x=element_blank()) +
  coord_cartesian(ylim=c(0, 4200))
DaysFirstLine

#Number Days Treatment, 2nd Line Antibiotics
DaysSecondLine <- ggplot(ND, aes(ScenAnalysis, TBB)) + geom_point(aes(shape = Scenario, colour = Scenario), size = 1.5, stroke = 1, position = position_jitter(0.05)) +
  scale_color_brewer(palette="Set2") +
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
  labs(x="",
       y="") + theme_classic()+ 
  theme(legend.position = "none", axis.text.x=element_blank()) +
  coord_cartesian(ylim=c(0, 1700))
DaysSecondLine

#Number Days Treatment, 3rd Line Antibiotics
DaysThirdLine <- ggplot(ND, aes(ScenAnalysis, TCC)) + 
  geom_point(aes(shape = Scenario, colour = Scenario), size = 1.5, stroke = 1, position = position_jitter(0.05)) +
  scale_color_brewer(palette = "Set2") +
  scale_x_discrete(labels = c("High, High" = "14/4", "High, Medium" = "14/2",
                              "High, Low" = "14/1","Medium, High" = "7/4",
                              "Medium, Medium" = "7/2","Medium, Low" = "7/1",
                              "Low, High" = "5/4","Low, Medium" = "5/2",
                              "Low, Low" = "5/1")) +
  theme(axis.text.x = element_text(angle = 65, vjust = 0.6)) +
  labs(x = "Treatment/Empiric Duration (Days)",
       y = "") +
  theme_classic() +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0, 1600))

# Plot --
ScenTreatCombined <- ggdraw() +
  draw_plot(NDMortality, x = 0.0, y = 0.81, width = 0.99, height = 0.15) +
  draw_plot(MR, x = 0.0, y = 0.68, width = 0.99, height = 0.15) +
  draw_plot(ResistanceLineOne, x = 0.0, y = 0.55, width = 0.99, height = 0.15) +
  draw_plot(ResistanceLineTwo, x = 0.0, y = 0.41, width = 0.99, height = 0.15) +
  draw_plot(DaysFirstLine, x = 0.0, y = 0.28, width = 0.99, height = 0.15) +
  draw_plot(DaysSecondLine, x = 0.0, y = 0.15, width = 0.99, height = 0.15) +
  draw_plot(DaysThirdLine, x = 0.0, y = 0.01, width = 0.99, height = 0.16)

ScenTreatCombined <- ScenTreatCombined +
  draw_label(
    "C) Treatment and Empiric Duration, T and d", fontface = 'bold', x = 0.01, y = .98, hjust = 0)
ScenTreatCombined



# Mortality ---------


#Different Scenarios
#Low (0.005), Marginal (1.5,2) Scenario A
parms <- c(L=0.0000466, t=7, a1=0.001, a2=0.001/5, a3=0.001, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=0.5, D=0.005, D1=0.005, D2=0.005*1.5, D3=0.005*2, NS=100000, e=1.3^((1/0.5)/2), e1=(1.3^((1/0.5)/2))^0, e2=(1.3^((1/0.5)/2))^0, e3=(1.3^((1/0.5)/2))^2, e4=(1.3^((1/0.5)/2)), e5=(1.3^((1/0.5)/2))^2, e6=(1.3^((1/0.5)/2)), e7=(1.3^((1/0.5)/2))^2, e8=(1.3^((1/0.5)/2)), gamma=0.000018, R1=(1/(7-(1/0.5))), R2=1/7, R3=(1/7), R4=(1/7), R5=(1/7), R6=(1/7), SA=1, SB=0, SC=0)
state <- c(h=0.00000212, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenALowMarginal<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenALowMarginal$cum_mort_rate <- outputScenALowMarginal$ND/(outputScenALowMarginal$NR + outputScenALowMarginal$ND)*100
outputScenALowMarginal$Scenario <- "A"
outputScenALowMarginal$ScenAnalysis <- "Low, Marginal"


#Low (0.005), Marginal (1.5,2), Scenario B
parms <- c(L=0.0000466, t=7, a1=0.001*1.375, a2=0.001*1.125, a3=0.001*1.5, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=0.5, D=0.005, D1=0.005, D2=0.005*1.5, D3=0.005*2, NS=100000, e=1.3^((1/0.5)/2), e1=(1.3^((1/0.5)/2))^0, e2=(1.3^((1/0.5)/2))^0, e3=(1.3^((1/0.5)/2))^1, e4=(1.3^((1/0.5)/2))^0, e5=(1.3^((1/0.5)/2))^1, e6=(1.3^((1/0.5)/2))^0, e7=(1.3^((1/0.5)/2))^1, e8=(1.3^((1/0.5)/2)), gamma=0.000018, R1=(1/(7-(1/0.5))), R2=1/7, R3=(1/(7-(1/0.5))), R4=(1/7), R5=(1/(7-(1/0.5))), R6=(1/7), SA=0, SB=1, SC=0)
state <- c(h=0.00000212, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenBLowMarginal<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenBLowMarginal$cum_mort_rate <- outputScenBLowMarginal$ND/(outputScenBLowMarginal$NR + outputScenBLowMarginal$ND)*100
outputScenBLowMarginal$Scenario <- "B"
outputScenBLowMarginal$ScenAnalysis <- "Low, Marginal"

#Low (0.005), Marginal (1.5,2), Scenario C
parms <- c(L=0.0000466, t=7, a1=0.001*1.375, a2=0.001*1.125, a3=0.001*1.5, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=0.5, D=0.005, D1=0.005, D2=0.005*1.5, D3=0.005*2, NS=100000, e=1.3^((1/0.5)/2), e1=(1.3^((1/0.5)/2))^0, e2=(1.3^((1/0.5)/2))^0, e3=(1.3^((1/0.5)/2))^0, e4=(1.3^((1/0.5)/2))^0, e5=(1.3^((1/0.5)/2))^0, e6=(1.3^((1/0.5)/2))^0, e7=(1.3^((1/0.5)/2))^0, e8=(1.3^((1/0.5)/2))^0, gamma=0.000018, R1=(1/(7-(1/0.5))), R2=1/7, R3=(1/(7-(1/0.5))), R4=(1/(7-(1/0.5))), R5=(1/(7-(1/0.5))), R6=(1/(7-(1/0.5))), SA=0, SB=0, SC=1)
state <- c(h=0.00000212, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenCLowMarginal<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenCLowMarginal$cum_mort_rate <- outputScenCLowMarginal$ND/(outputScenCLowMarginal$NR + outputScenCLowMarginal$ND)*100
outputScenCLowMarginal$Scenario <- "C"
outputScenCLowMarginal$ScenAnalysis <- "Low, Marginal"


#Break

#Low (0.005), Wide (2,4) Scenario A
parms <- c(L=0.0000466, t=7, a1=0.001, a2=0.001/5, a3=0.001, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=0.5, D=0.005, D1=0.005, D2=0.005*2, D3=0.005*4, NS=100000, e=1.3^((1/0.5)/2), e1=(1.3^((1/0.5)/2))^0, e2=(1.3^((1/0.5)/2))^0, e3=(1.3^((1/0.5)/2))^2, e4=(1.3^((1/0.5)/2)), e5=(1.3^((1/0.5)/2))^2, e6=(1.3^((1/0.5)/2)), e7=(1.3^((1/0.5)/2))^2, e8=(1.3^((1/0.5)/2)), gamma=0.000018, R1=(1/(7-(1/0.5))), R2=1/7, R3=(1/7), R4=(1/7), R5=(1/7), R6=(1/7), SA=1, SB=0, SC=0)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenALowWide<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenALowWide$cum_mort_rate <- outputScenALowWide$ND/(outputScenALowWide$NR + outputScenALowWide$ND)*100
outputScenALowWide$Scenario <- "A"
outputScenALowWide$ScenAnalysis <- "Low, Wide"

#Low (0.005), Wide (2,4), Scenario B
parms <- c(L=0.0000466, t=7, a1=0.001*1.375, a2=0.001*1.125, a3=0.001*1.5, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=0.5, D=0.005, D1=0.005, D2=0.005*2, D3=0.005*4, NS=100000, e=1.3^((1/0.5)/2), e1=(1.3^((1/0.5)/2))^0, e2=(1.3^((1/0.5)/2))^0, e3=(1.3^((1/0.5)/2))^1, e4=(1.3^((1/0.5)/2))^0, e5=(1.3^((1/0.5)/2))^1, e6=(1.3^((1/0.5)/2))^0, e7=(1.3^((1/0.5)/2))^1, e8=(1.3^((1/0.5)/2)), gamma=0.000018, R1=(1/(7-(1/0.5))), R2=1/7, R3=(1/(7-(1/0.5))), R4=(1/7), R5=(1/(7-(1/0.5))), R6=(1/7), SA=0, SB=1, SC=0)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenBLowWide<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenBLowWide$cum_mort_rate <- outputScenBLowWide$ND/(outputScenBLowWide$NR + outputScenBLowWide$ND)*100
outputScenBLowWide$Scenario <- "B"
outputScenBLowWide$ScenAnalysis <- "Low, Wide"

#Low (0.005), Wide (2,4), Scenario C
parms <- c(L=0.0000466, t=7, a1=0.001*1.375, a2=0.001*1.125, a3=0.001*1.5, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=0.5, D=0.005, D1=0.005, D2=0.005*2, D3=0.005*4, NS=100000, e=1.3^((1/0.5)/2), e1=(1.3^((1/0.5)/2))^0, e2=(1.3^((1/0.5)/2))^0, e3=(1.3^((1/0.5)/2))^0, e4=(1.3^((1/0.5)/2))^0, e5=(1.3^((1/0.5)/2))^0, e6=(1.3^((1/0.5)/2))^0, e7=(1.3^((1/0.5)/2))^0, e8=(1.3^((1/0.5)/2))^0, gamma=0.000018, R1=(1/(7-(1/0.5))), R2=1/7, R3=(1/(7-(1/0.5))), R4=(1/(7-(1/0.5))), R5=(1/(7-(1/0.5))), R6=(1/(7-(1/0.5))), SA=0, SB=0, SC=1)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenCLowWide<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenCLowWide$cum_mort_rate <- outputScenCLowWide$ND/(outputScenCLowWide$NR + outputScenCLowWide$ND)*100
outputScenCLowWide$Scenario <- "C"
outputScenCLowWide$ScenAnalysis <- "Low, Wide"

#Break

#Medium (0.01), Marginal (1.5,2) Scenario A
parms <- c(L=0.0000466, t=7, a1=0.001, a2=0.001/5, a3=0.001, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=0.5, D=0.01, D1=0.01, D2=0.01*1.5, D3=0.01*2, NS=100000, e=1.3^((1/0.5)/2), e1=(1.3^((1/0.5)/2))^0, e2=(1.3^((1/0.5)/2))^0, e3=(1.3^((1/0.5)/2))^2, e4=(1.3^((1/0.5)/2)), e5=(1.3^((1/0.5)/2))^2, e6=(1.3^((1/0.5)/2)), e7=(1.3^((1/0.5)/2))^2, e8=(1.3^((1/0.5)/2)), gamma=0.000018, R1=(1/(7-(1/0.5))), R2=1/7, R3=(1/7), R4=(1/7), R5=(1/7), R6=(1/7), SA=1, SB=0, SC=0)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenAMedMarginal<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenAMedMarginal$cum_mort_rate <- outputScenAMedMarginal$ND/(outputScenAMedMarginal$NR + outputScenAMedMarginal$ND)*100
outputScenAMedMarginal$Scenario <- "A"
outputScenAMedMarginal$ScenAnalysis <- "Medium, Marginal"


#Graphs
options(scipen=999)  # turn off scientific notation like 1e+06

#Medium (0.01), Marginal (1.5,2), Scenario B
parms <- c(L=0.0000466, t=7, a1=0.001*1.375, a2=0.001*1.125, a3=0.001*1.5, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=0.5, D=0.01, D1=0.01, D2=0.01*1.5, D3=0.01*2, NS=100000, e=1.3^((1/0.5)/2), e1=(1.3^((1/0.5)/2))^0, e2=(1.3^((1/0.5)/2))^0, e3=(1.3^((1/0.5)/2))^1, e4=(1.3^((1/0.5)/2))^0, e5=(1.3^((1/0.5)/2))^1, e6=(1.3^((1/0.5)/2))^0, e7=(1.3^((1/0.5)/2))^1, e8=(1.3^((1/0.5)/2)), gamma=0.000018, R1=(1/(7-(1/0.5))), R2=1/7, R3=(1/(7-(1/0.5))), R4=(1/7), R5=(1/(7-(1/0.5))), R6=(1/7), SA=0, SB=1, SC=0)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenBMedMarginal<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenBMedMarginal$cum_mort_rate <- outputScenBMedMarginal$ND/(outputScenBMedMarginal$NR + outputScenBMedMarginal$ND)*100
outputScenBMedMarginal$Scenario <- "B"
outputScenBMedMarginal$ScenAnalysis <- "Medium, Marginal"

#Medium (0.01), Marginal (1.5,2), Scenario C
parms <- c(L=0.0000466, t=7, a1=0.001*1.375, a2=0.001*1.125, a3=0.001*1.5, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=0.5, D=0.01, D1=0.01, D2=0.01*1.5, D3=0.01*2, NS=100000, e=1.3^((1/0.5)/2), e1=(1.3^((1/0.5)/2))^0, e2=(1.3^((1/0.5)/2))^0, e3=(1.3^((1/0.5)/2))^0, e4=(1.3^((1/0.5)/2))^0, e5=(1.3^((1/0.5)/2))^0, e6=(1.3^((1/0.5)/2))^0, e7=(1.3^((1/0.5)/2))^0, e8=(1.3^((1/0.5)/2))^0, gamma=0.000018, R1=(1/(7-(1/0.5))), R2=1/7, R3=(1/(7-(1/0.5))), R4=(1/(7-(1/0.5))), R5=(1/(7-(1/0.5))), R6=(1/(7-(1/0.5))), SA=0, SB=0, SC=1)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenCMedMarginal<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenCMedMarginal$cum_mort_rate <- outputScenCMedMarginal$ND/(outputScenCMedMarginal$NR + outputScenCMedMarginal$ND)*100
outputScenCMedMarginal$Scenario <- "C"
outputScenCMedMarginal$ScenAnalysis <- "Medium, Marginal"

#Break

#Medium (0.01), Wide (2,4) Scenario A
parms <- c(L=0.0000466, t=7, a1=0.001, a2=0.001/5, a3=0.001, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=0.5, D=0.01, D1=0.01, D2=0.01*2, D3=0.01*4, NS=100000, e=1.3^((1/0.5)/2), e1=(1.3^((1/0.5)/2))^0, e2=(1.3^((1/0.5)/2))^0, e3=(1.3^((1/0.5)/2))^2, e4=(1.3^((1/0.5)/2)), e5=(1.3^((1/0.5)/2))^2, e6=(1.3^((1/0.5)/2)), e7=(1.3^((1/0.5)/2))^2, e8=(1.3^((1/0.5)/2)), gamma=0.000018, R1=(1/(7-(1/0.5))), R2=1/7, R3=(1/7), R4=(1/7), R5=(1/7), R6=(1/7), SA=1, SB=0, SC=0)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenAMedWide<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenAMedWide$cum_mort_rate <- outputScenAMedWide$ND/(outputScenAMedWide$NR + outputScenAMedWide$ND)*100
outputScenAMedWide$Scenario <- "A"
outputScenAMedWide$ScenAnalysis <- "Medium, Wide"

#Medium (0.01), Wide (2,4), Scenario B
parms <- c(L=0.0000466, t=7, a1=0.001*1.375, a2=0.001*1.125, a3=0.001*1.5, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=0.5, D=0.01, D1=0.01, D2=0.01*2, D3=0.01*4, NS=100000, e=1.3^((1/0.5)/2), e1=(1.3^((1/0.5)/2))^0, e2=(1.3^((1/0.5)/2))^0, e3=(1.3^((1/0.5)/2))^1, e4=(1.3^((1/0.5)/2))^0, e5=(1.3^((1/0.5)/2))^1, e6=(1.3^((1/0.5)/2))^0, e7=(1.3^((1/0.5)/2))^1, e8=(1.3^((1/0.5)/2)), gamma=0.000018, R1=(1/(7-(1/0.5))), R2=1/7, R3=(1/(7-(1/0.5))), R4=(1/7), R5=(1/(7-(1/0.5))), R6=(1/7), SA=0, SB=1, SC=0)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenBMedWide<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenBMedWide$cum_mort_rate <- outputScenBMedWide$ND/(outputScenBMedWide$NR + outputScenBMedWide$ND)*100
outputScenBMedWide$Scenario <- "B"
outputScenBMedWide$ScenAnalysis <- "Medium, Wide"

#Medium (0.01), Wide (2,4), Scenario C
parms <- c(L=0.0000466, t=7, a1=0.001*1.375, a2=0.001*1.125, a3=0.001*1.5, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=0.5, D=0.01, D1=0.01, D2=0.01*2, D3=0.01*4, NS=100000, e=1.3^((1/0.5)/2), e1=(1.3^((1/0.5)/2))^0, e2=(1.3^((1/0.5)/2))^0, e3=(1.3^((1/0.5)/2))^0, e4=(1.3^((1/0.5)/2))^0, e5=(1.3^((1/0.5)/2))^0, e6=(1.3^((1/0.5)/2))^0, e7=(1.3^((1/0.5)/2))^0, e8=(1.3^((1/0.5)/2))^0, gamma=0.000018, R1=(1/(7-(1/0.5))), R2=1/7, R3=(1/(7-(1/0.5))), R4=(1/(7-(1/0.5))), R5=(1/(7-(1/0.5))), R6=(1/(7-(1/0.5))), SA=0, SB=0, SC=1)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenCMedWide<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenCMedWide$cum_mort_rate <- outputScenCMedWide$ND/(outputScenCMedWide$NR + outputScenCMedWide$ND)*100
outputScenCMedWide$Scenario <- "C"
outputScenCMedWide$ScenAnalysis <- "Medium, Wide"

#Break

#High (0.02), Marginal (1.5,2) Scenario A
parms <- c(L=0.0000466, t=7, a1=0.001, a2=0.001/5, a3=0.001, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=0.5, D=0.02, D1=0.02, D2=0.02*1.5, D3=0.02*2, NS=100000, e=1.3^((1/0.5)/2), e1=(1.3^((1/0.5)/2))^0, e2=(1.3^((1/0.5)/2))^0, e3=(1.3^((1/0.5)/2))^2, e4=(1.3^((1/0.5)/2)), e5=(1.3^((1/0.5)/2))^2, e6=(1.3^((1/0.5)/2)), e7=(1.3^((1/0.5)/2))^2, e8=(1.3^((1/0.5)/2)), gamma=0.000018, R1=(1/(7-(1/0.5))), R2=1/7, R3=(1/7), R4=(1/7), R5=(1/7), R6=(1/7), SA=1, SB=0, SC=0)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenAHighMarginal<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenAHighMarginal$cum_mort_rate <- outputScenAHighMarginal$ND/(outputScenAHighMarginal$NR + outputScenAHighMarginal$ND)*100
outputScenAHighMarginal$Scenario <- "A"
outputScenAHighMarginal$ScenAnalysis <- "High, Marginal"


#High (0.02), Marginal (1.5,2), Scenario B
parms <- c(L=0.0000466, t=7, a1=0.001*1.375, a2=0.001*1.125, a3=0.001*1.5, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=0.5, D=0.02, D1=0.02, D2=0.02*1.5, D3=0.02*2, NS=100000, e=1.3^((1/0.5)/2), e1=(1.3^((1/0.5)/2))^0, e2=(1.3^((1/0.5)/2))^0, e3=(1.3^((1/0.5)/2))^1, e4=(1.3^((1/0.5)/2))^0, e5=(1.3^((1/0.5)/2))^1, e6=(1.3^((1/0.5)/2))^0, e7=(1.3^((1/0.5)/2))^1, e8=(1.3^((1/0.5)/2)), gamma=0.000018, R1=(1/(7-(1/0.5))), R2=1/7, R3=(1/(7-(1/0.5))), R4=(1/7), R5=(1/(7-(1/0.5))), R6=(1/7), SA=0, SB=1, SC=0)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenBHighMarginal<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenBHighMarginal$cum_mort_rate <- outputScenBHighMarginal$ND/(outputScenBHighMarginal$NR + outputScenBHighMarginal$ND)*100
outputScenBHighMarginal$Scenario <- "B"
outputScenBHighMarginal$ScenAnalysis <- "High, Marginal"

#High (0.02), Marginal (1.5,2), Scenario C
parms <- c(L=0.0000466, t=7, a1=0.001*1.375, a2=0.001*1.125, a3=0.001*1.5, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=0.5, D=0.02, D1=0.02, D2=0.02*1.5, D3=0.02*2, NS=100000, e=1.3^((1/0.5)/2), e1=(1.3^((1/0.5)/2))^0, e2=(1.3^((1/0.5)/2))^0, e3=(1.3^((1/0.5)/2))^0, e4=(1.3^((1/0.5)/2))^0, e5=(1.3^((1/0.5)/2))^0, e6=(1.3^((1/0.5)/2))^0, e7=(1.3^((1/0.5)/2))^0, e8=(1.3^((1/0.5)/2))^0, gamma=0.000018, R1=(1/(7-(1/0.5))), R2=1/7, R3=(1/(7-(1/0.5))), R4=(1/(7-(1/0.5))), R5=(1/(7-(1/0.5))), R6=(1/(7-(1/0.5))), SA=0, SB=0, SC=1)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenCHighMarginal<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenCHighMarginal$cum_mort_rate <- outputScenCHighMarginal$ND/(outputScenCHighMarginal$NR + outputScenCHighMarginal$ND)*100
outputScenCHighMarginal$Scenario <- "C"
outputScenCHighMarginal$ScenAnalysis <- "High, Marginal"

#Break

#High (0.02), Wide (2,4) Scenario A
parms <- c(L=0.0000466, t=7, a1=0.001, a2=0.001/5, a3=0.001, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=0.5, D=0.02, D1=0.02, D2=0.02*2, D3=0.02*4, NS=100000, e=1.3^((1/0.5)/2), e1=(1.3^((1/0.5)/2))^0, e2=(1.3^((1/0.5)/2))^0, e3=(1.3^((1/0.5)/2))^2, e4=(1.3^((1/0.5)/2)), e5=(1.3^((1/0.5)/2))^2, e6=(1.3^((1/0.5)/2)), e7=(1.3^((1/0.5)/2))^2, e8=(1.3^((1/0.5)/2)), gamma=0.000018, R1=(1/(7-(1/0.5))), R2=1/7, R3=(1/7), R4=(1/7), R5=(1/7), R6=(1/7), SA=1, SB=0, SC=0)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenAHighWide<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenAHighWide$cum_mort_rate <- outputScenAHighWide$ND/(outputScenAHighWide$NR + outputScenAHighWide$ND)*100
outputScenAHighWide$Scenario <- "A"
outputScenAHighWide$ScenAnalysis <- "High, Wide"

#High (0.02), Wide (2,4), Scenario B
parms <- c(L=0.0000466, t=7, a1=0.001*1.375, a2=0.001*1.125, a3=0.001*1.5, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=0.5, D=0.02, D1=0.02, D2=0.02*2, D3=0.02*4, NS=100000, e=1.3^((1/0.5)/2), e1=(1.3^((1/0.5)/2))^0, e2=(1.3^((1/0.5)/2))^0, e3=(1.3^((1/0.5)/2))^1, e4=(1.3^((1/0.5)/2))^0, e5=(1.3^((1/0.5)/2))^1, e6=(1.3^((1/0.5)/2))^0, e7=(1.3^((1/0.5)/2))^1, e8=(1.3^((1/0.5)/2)), gamma=0.000018, R1=(1/(7-(1/0.5))), R2=1/7, R3=(1/(7-(1/0.5))), R4=(1/7), R5=(1/(7-(1/0.5))), R6=(1/7), SA=0, SB=1, SC=0)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenBHighWide<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenBHighWide$cum_mort_rate <- outputScenBHighWide$ND/(outputScenBHighWide$NR + outputScenBHighWide$ND)*100
outputScenBHighWide$Scenario <- "B"
outputScenBHighWide$ScenAnalysis <- "High, Wide"

#High (0.02), Wide (2,4), Scenario C
parms <- c(L=0.0000466, t=7, a1=0.001*1.375, a2=0.001*1.125, a3=0.001*1.5, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=0.5, D=0.02, D1=0.02, D2=0.02*2, D3=0.02*4, NS=100000, e=1.3^((1/0.5)/2), e1=(1.3^((1/0.5)/2))^0, e2=(1.3^((1/0.5)/2))^0, e3=(1.3^((1/0.5)/2))^0, e4=(1.3^((1/0.5)/2))^0, e5=(1.3^((1/0.5)/2))^0, e6=(1.3^((1/0.5)/2))^0, e7=(1.3^((1/0.5)/2))^0, e8=(1.3^((1/0.5)/2))^0, gamma=0.000018, R1=(1/(7-(1/0.5))), R2=1/7, R3=(1/(7-(1/0.5))), R4=(1/(7-(1/0.5))), R5=(1/(7-(1/0.5))), R6=(1/(7-(1/0.5))), SA=0, SB=0, SC=1)
state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
outputScenCHighWide<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))

#Add new column
outputScenCHighWide$cum_mort_rate <- outputScenCHighWide$ND/(outputScenCHighWide$NR + outputScenCHighWide$ND)*100
outputScenCHighWide$Scenario <- "C"
outputScenCHighWide$ScenAnalysis <- "High, Wide"

#Merging datasets
DatasetScenMort <- rbind(outputScenALowMarginal, outputScenBLowMarginal, outputScenCLowMarginal, 
                         outputScenALowWide, outputScenBLowWide, outputScenCLowWide,
                         outputScenAMedMarginal, outputScenBMedMarginal, outputScenCMedMarginal,
                         outputScenAMedWide, outputScenBMedWide, outputScenCMedWide,
                         outputScenAHighMarginal, outputScenBHighMarginal, outputScenCHighMarginal,
                         outputScenAHighWide, outputScenBHighWide, outputScenCHighWide)


#Graphs

#Mortality
ND <- subset(DatasetScenMort, time==1825)
ND$ScenAnalysis <- factor(ND$ScenAnalysis,levels = c("High, Wide", "High, Marginal", "Medium, Wide", "Medium, Marginal",
                                                     "Low, Wide", "Low, Marginal"))


NDMortality <- ggplot(ND, aes(ScenAnalysis, ND)) + geom_point(aes(shape = Scenario, colour = Scenario), size = 1.5, stroke = 1, position = position_jitter(0.05)) +
  scale_color_brewer(palette="Set2") +
  labs(x="",
       y="") + theme_classic() +
  coord_cartesian(ylim=c(0, 70)) +
  theme(legend.position = "none", axis.text.x=element_blank())
NDMortality

#Mortality Rate
MR <- ggplot(ND, aes(ScenAnalysis, cum_mort_rate)) + geom_point(aes(shape = Scenario, colour = Scenario), size = 1.5, stroke = 1, position = position_jitter(0.05)) +
  scale_color_brewer(palette="Set2") +
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
  labs(x="",
       y="") + theme_classic() +
  theme(legend.position = "none", axis.text.x=element_blank()) +
  coord_cartesian(ylim=c(0, 16.5))
MR

#Resistance to 1st Line Therapy
ResistanceLineOne <- ggplot(ND, aes(ScenAnalysis, WB*100)) + geom_point(aes(shape = Scenario, colour = Scenario), size = 1.5, stroke = 1, position = position_jitter(0.05)) +
  scale_color_brewer(palette="Set2") +
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
  labs(x="",
       y="") + theme_classic()+ 
  theme(legend.position = "none", axis.text.x=element_blank()) +
  coord_cartesian(ylim=c(0, 8))
ResistanceLineOne

#Resistance to 2nd Line Therapy
ResistanceLineTwo <- ggplot(ND, aes(ScenAnalysis, WC*100)) + geom_point(aes(shape = Scenario, colour = Scenario), size = 1.5, stroke = 1, position = position_jitter(0.05)) +
  scale_color_brewer(palette="Set2") +
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
  labs(x="",
       y="") + theme_classic()+ 
  theme(legend.position = "none", axis.text.x=element_blank()) +
  coord_cartesian(ylim=c(0, 7))
ResistanceLineTwo


#Number Days Treatment, 1st Line Antibiotics
options(scipen=999)
DaysFirstLine <- ggplot(ND, aes(ScenAnalysis, TAA)) + geom_point(aes(shape = Scenario, colour = Scenario), size = 1.5, stroke = 1, position = position_jitter(0.05)) +
  scale_color_brewer(palette="Set2") +
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
  labs(x="",
       y="") + theme_classic()+ 
  theme(legend.position = "none", axis.text.x=element_blank()) +
  coord_cartesian(ylim=c(0, 2600))
DaysFirstLine

#Number Days Treatment, 2nd Line Antibiotics
DaysSecondLine <- ggplot(ND, aes(ScenAnalysis, TBB)) + geom_point(aes(shape = Scenario, colour = Scenario), size = 1.5, stroke = 1, position = position_jitter(0.05)) +
  scale_color_brewer(palette="Set2") +
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
  labs(x="",
       y="") + theme_classic() + 
  theme(legend.position = "none", axis.text.x=element_blank()) +
  coord_cartesian(ylim=c(0, 1000))

#Number Days Treatment, 3rd Line Antibiotics
DaysThirdLine <- ggplot(ND, aes(ScenAnalysis, TCC)) + geom_point(aes(shape = Scenario, colour = Scenario), size = 1.5, stroke = 1, position = position_jitter(0.05)) +
  scale_color_brewer(palette="Set2") +
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
  labs(x="Baseline Mortality Rate",
       y="") + theme_classic()+ 
  theme(legend.position = "none") +
  coord_cartesian(ylim=c(0, 1000))


# Combining into one
ScenMortalityCombined <- ggdraw() +
  draw_plot(NDMortality, x = 0.0, y = 0.81, width = 0.99, height = 0.15) +
  draw_plot(MR, x = 0.0, y = 0.68, width = 0.99, height = 0.15) +
  draw_plot(ResistanceLineOne, x = 0.0, y = 0.55, width = 0.99, height = 0.15) +
  draw_plot(ResistanceLineTwo, x = 0.0, y = 0.41, width = 0.99, height = 0.15) +
  draw_plot(DaysFirstLine, x = 0.0, y = 0.28, width = 0.99, height = 0.15) +
  draw_plot(DaysSecondLine, x = 0.0, y = 0.15, width = 0.99, height = 0.15) +
  draw_plot(DaysThirdLine, x = 0.0, y = 0.01, width = 0.99, height = 0.16)

ScenMortalityCombined <- ScenMortalityCombined +
  draw_label(
    "D) Baseline Mortality", fontface = 'bold', x = 0.01, y = .98, hjust = 0)
ScenMortalityCombined








# Combined Figure -----------
# Create the legend plot
DaysThirdLine <- ggplot(ND, aes(ScenAnalysis, TCC)) +
  geom_point(aes(shape = Scenario, colour = Scenario), size = 1.5, stroke = 1, position = position_jitter(0.05)) +
  scale_color_brewer(palette = "Set2") +
  theme(axis.text.x = element_text(angle = 65, vjust = 0.6)) +
  labs(x = "Baseline Mortality Rate", y = "") +
  theme_classic() +
  theme(legend.position = "bottom") +
  coord_cartesian(ylim = c(0, 1000))

# Extract the legend from the original plot
legend_plot <- get_legend(DaysThirdLine)
legend_ggplot <- ggdraw() + draw_plot(legend_plot)

# Combine the top row plots into a single row
top_row_plots <- ScenIEATCombined + ScenResCombined + ScenTreatCombined + ScenMortalityCombined + plot_layout(nrow = 1, widths = c(23, 23, 27,27))

# Create the final patchwork plot with the specified layout
final_plot <- top_row_plots / legend_ggplot + plot_layout(nrow = 2, heights = c(97, 3))
print(final_plot)











