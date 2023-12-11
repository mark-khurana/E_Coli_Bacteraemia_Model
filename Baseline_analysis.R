# Baseline model

#Load relevant libraries
library(deSolve)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(cowplot)

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


# Generate Normal Distribution values for parameters
set.seed(123)
num_samples <- 1000

# Adjustable parameters -----
# Function to generate normal distributions for adjustable parameters
generate_normal_distribution <- function(mean, sd, n, min_val, max_val) {
  result <- rnorm(n, mean, sd)
  result[result < min_val] <- min_val
  result[result > max_val] <- max_val
  return(result)
}

# Adjustable parameters with specified means and tighter standard deviations
adjustable_parameters <- list(
  a = list(mean = 0.001, sd = 0.00075, min_val = 0.0005, max_val = 0.01),        # 95% interval: [0.0005, 0.01]
  D = list(mean = 0.01, sd = 0.0025, min_val = 0.005, max_val = 0.02),           # 95% interval: [0.005, 0.02]
  d = list(mean = 0.5, sd = 0.0625, min_val = 0.375, max_val = 0.625),           # 95% interval: [0.375, 0.625]
  gamma_parameter = list(mean = 0.000018, sd = 0.00002, min_val = 0.000001, max_val = 0.0001),  # 95% interval: [0.000001, 0.0001]
  eps = list(mean = 1.3, sd = 0.0375, min_val = 1.0, max_val = 1.5),              # 95% interval: [1.0, 1.5]
  duration = list(mean = 7, sd = 0.5, min_val = 5, max_val = 10)                  # 95% interval: [5, 9]
)

parameter_samples <- lapply(adjustable_parameters, function(param) {
  generate_normal_distribution(param$mean, param$sd, num_samples, param$min_val, param$max_val)
})

# Function to run the ODE model
run_model <- function(parameters) {
  output <- as.data.frame(ode(func = new.basic.model, y = state, times = times, parms = parameters))
  return(output)
}



# Model Set-Up, Baseline, Scenario A --------

# Initialize an empty list to store results
results_Scenario_A <- list()

# Run ODE Model for Each Sample and Collect Results
for (i in 1:num_samples) {
  # Set parameters based on Latin hypercube sample
  a = parameter_samples$a[i]
  D = parameter_samples$D[i]
  d = parameter_samples$d[i]
  gamma_parameter = parameter_samples$gamma_parameter[i]
  eps = parameter_samples$eps[i]
  e = eps^((1/d)/2) 
  duration = parameter_samples$duration[i]
  parms <- c(L=0.0000466, t=duration, a1=a, a2=a/5, a3=a, a4=a, a5=a/5, a6=a, a7=a*1.5, d=d, D1=D, D2=D*1.5, D3=D*2, NS=100000, e=e, e1=e^0, e2=(e)^0, e3=(e)^2, e4=(e), e5=(e)^2, e6=(e), e7=(e)^2, e8=(e), gamma=gamma_parameter, R1=(1/(duration-(1/d))), R2=(1/duration), R3=(1/duration), R4=(1/duration), R5=(1/duration), R6=(1/duration), SA=1, SB=0, SC=0)
  state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
  times <- seq(from=0, to=365*5, by=25)
  sampled_parameters <- parms
  
  # Run the ODE model
  outputBaselineA <- as.data.frame(ode(func = new.basic.model, y = state, times = times, parms = parms))
  
  # Calculate cumulative mortality rate
  outputBaselineA$cum_mort_rate <- outputBaselineA$ND / (outputBaselineA$NR + outputBaselineA$ND) * 100
  
  # Add Scenario column
  outputBaselineA$Scenario <- "A"
  
  # Add Run column
  outputBaselineA$Run <- i
  
  # Store the results in the list
  results_Scenario_A[[i]] <- outputBaselineA
}

# Combine results into a data table
results_df_Scenario_A <- do.call(rbind, results_Scenario_A)

# Extract values at 1825 days (Five years)
values_at_1825_days_ScenarioA <- results_df_Scenario_A[results_df_Scenario_A$time == 1825, ]
# Summary statistics for each variable
summary_stats_ScenarioA <- summary(values_at_1825_days_ScenarioA[, -c(1, ncol(values_at_1825_days_ScenarioA))])
# Print the summary statistics
print(summary_stats_ScenarioA)


# Baseline, Scenario B ----------

# Initialize an empty list to store results
results_Scenario_B <- list()

# Run ODE Model for Each Sample and Collect Results
for (i in 1:num_samples) {
  # Set parameters based on Latin hypercube sample
  a = parameter_samples$a[i]
  D = parameter_samples$D[i]
  d = parameter_samples$d[i]
  gamma_parameter = parameter_samples$gamma_parameter[i]
  eps = parameter_samples$eps[i]
  e = eps^((1/d)/2) 
  duration = parameter_samples$duration[i]
  parms <- c(L=0.0000466, t=duration, a1=a*1.375, a2=a*1.125, a3=a*1.5, a4=a, a5=a/5, a6=a, a7=a*1.5, d=d, D=D, D1=D, D2=D*1.5, D3=D*2, NS=100000, e=e, e1=e^0, e2=e^0, e3=e^1, e4=e^0, e5=e^1, e6=e^0, e7=e^1, e8=e, gamma=gamma_parameter, R1=(1/(duration-(1/d))), R2=1/duration, R3=(1/(duration-(1/d))), R4=(1/duration), R5=(1/(duration-(1/d))), R6=(1/duration), SA=0, SB=1, SC=0)
  state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
  times <- seq(from=0, to=365*5, by=25)
  sampled_parameters <- parms
  
  # Run the ODE model
  outputBaselineB <- as.data.frame(ode(func = new.basic.model, y = state, times = times, parms = parms))
  
  # Calculate cumulative mortality rate
  outputBaselineB$cum_mort_rate <- outputBaselineB$ND / (outputBaselineB$NR + outputBaselineB$ND) * 100
  
  # Add Scenario column
  outputBaselineB$Scenario <- "B"
  
  # Add Run column
  outputBaselineB$Run <- i
  
  # Store the results in the list
  results_Scenario_B[[i]] <- outputBaselineB
}

# Combine results into a data table
results_df_Scenario_B <- do.call(rbind, results_Scenario_B)

# Extract values at 1825 days (Five years)
values_at_1825_days_ScenarioB <- results_df_Scenario_B[results_df_Scenario_B$time == 1825, ]
# Summary statistics for each variable
summary_stats_ScenarioB <- summary(values_at_1825_days_ScenarioB[, -c(1, ncol(values_at_1825_days_ScenarioB))])
# Print the summary statistics
print(summary_stats_ScenarioB)




# Baseline, Scenario C -------

# Initialize an empty list to store results
results_Scenario_C <- list()

# Run ODE Model for Each Sample and Collect Results
for (i in 1:num_samples) {
  # Set parameters based on Latin hypercube sample
  a = parameter_samples$a[i]
  D = parameter_samples$D[i]
  d = parameter_samples$d[i]
  gamma_parameter = parameter_samples$gamma_parameter[i]
  eps = parameter_samples$eps[i]
  e = eps^((1/d)/2) 
  duration = parameter_samples$duration[i]
  parms <- c(L=0.0000466, t=duration, a1=a*1.375, a2=a*1.125, a3=a*1.5, a4=a, a5=a/5, a6=a, a7=a*1.5, d=d, D=D, D1=D, D2=D*1.5, D3=D*2, NS=100000, e=e, e1=e^0, e2=e^0, e3=e^0, e4=e^0, e5=e^0, e6=e^0, e7=e^0, e8=e^0, gamma=gamma_parameter, R1=(1/(duration-(1/d))), R2=(1/duration), R3=(1/(duration-(1/d))), R4=(1/(duration-(1/d))), R5=(1/(duration-(1/d))), R6=(1/(duration-(1/d))), SA=0, SB=0, SC=1)
  state <- c(h=0.00000188, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
  times <- seq(from=0, to=365*5, by=25)
  sampled_parameters <- parms
  
  # Run the ODE model
  outputBaselineC <- as.data.frame(ode(func = new.basic.model, y = state, times = times, parms = parms))
  
  # Calculate cumulative mortality rate
  outputBaselineC$cum_mort_rate <- outputBaselineC$ND / (outputBaselineC$NR + outputBaselineC$ND) * 100
  
  # Add Scenario column
  outputBaselineC$Scenario <- "C"
  
  # Add Run column
  outputBaselineC$Run <- i
  
  # Store the results in the list
  results_Scenario_C[[i]] <- outputBaselineC
}

# Combine results into a data table
results_df_Scenario_C <- do.call(rbind, results_Scenario_C)






# Combining all datasets ------
combined_baseline_df <- bind_rows(
  mutate(results_df_Scenario_A),
  mutate(results_df_Scenario_B),
  mutate(results_df_Scenario_C)
) %>%
  select(time, Run, Scenario, ND, cum_mort_rate, WB, WC, TAA, TBB, TCC)







# Summary results for table -------
# Filter data at time = 1825
summary_df <- combined_baseline_df %>% filter(time == 1825)

# Summarize values (mean, SD, IQR, 95% CI, etc.) for each variable stratified by Scenario
summary_stats <- summary_df %>%
  group_by(Scenario) %>%
  summarise(across(c(ND, cum_mort_rate, WB, WC, TAA, TBB, TCC),
                   list(mean = mean, sd = sd))) # Note that WB and WC are given as decimals, not percentages

# Find the differences between scenarios A and B, B and C, A and C for each Run
diff_df <- summary_df %>%
  group_by(Run) %>%
  summarise(
    diff_AB_ND = (ND[Scenario == "A"] - ND[Scenario == "B"]),
    diff_BC_ND = (ND[Scenario == "B"] - ND[Scenario == "C"]),
    diff_AC_ND = (ND[Scenario == "A"] - ND[Scenario == "C"]),
    diff_AB_cum_mort_rate = (cum_mort_rate[Scenario == "A"] - cum_mort_rate[Scenario == "B"]),
    diff_BC_cum_mort_rate = (cum_mort_rate[Scenario == "B"] - cum_mort_rate[Scenario == "C"]),
    diff_AC_cum_mort_rate = (cum_mort_rate[Scenario == "A"] - cum_mort_rate[Scenario == "C"]),
    diff_AB_WB = ((WB[Scenario == "A"] - WB[Scenario == "B"]) * 100), #To ensure that they are given as percentages
    diff_BC_WB = ((WB[Scenario == "B"] - WB[Scenario == "C"]) * 100), #To ensure that they are given as percentages
    diff_AC_WB = ((WB[Scenario == "A"] - WB[Scenario == "C"]) * 100), #To ensure that they are given as percentages
    diff_AB_WC = ((WC[Scenario == "A"] - WC[Scenario == "B"]) * 100), #To ensure that they are given as percentages
    diff_BC_WC = ((WC[Scenario == "B"] - WC[Scenario == "C"]) * 100), #To ensure that they are given as percentages
    diff_AC_WC = ((WC[Scenario == "A"] - WC[Scenario == "C"]) * 100), #To ensure that they are given as percentages
    diff_AB_TAA = (TAA[Scenario == "A"] - TAA[Scenario == "B"]),
    diff_BC_TAA = (TAA[Scenario == "B"] - TAA[Scenario == "C"]),
    diff_AC_TAA = (TAA[Scenario == "A"] - TAA[Scenario == "C"]),
    diff_AB_TBB = (TBB[Scenario == "A"] - TBB[Scenario == "B"]),
    diff_BC_TBB = (TBB[Scenario == "B"] - TBB[Scenario == "C"]),
    diff_AC_TBB = (TBB[Scenario == "A"] - TBB[Scenario == "C"]),
    diff_AB_TCC = (TCC[Scenario == "A"] - TCC[Scenario == "B"]),
    diff_BC_TCC = (TCC[Scenario == "B"] - TCC[Scenario == "C"]),
    diff_AC_TCC = (TCC[Scenario == "A"] - TCC[Scenario == "C"])
  )

# Ungroup the dataset before calculating summary statistics
diff_df <- diff_df %>% ungroup()

# Calculate mean, SD, IQR, and 95% CI for each column in diff_df
diff_stats_summary <- diff_df %>%
  reframe(across(starts_with("diff"),
                 list(mean = ~mean(.x), sd = ~sd(.x)),
                 .names = "{.col}_{.fn}")
  )

# Print the results
print(summary_stats)
print(diff_stats_summary)





# Visualizing -------------

# Dark Mint color palette
my_palette <- c("#66c0a5", "#fc8d62", "#8da0cb")

# Assuming your data frame is named combined_baseline_df

# Convert Scenario and Run to factors for better plotting
combined_baseline_df$Scenario <- as.factor(combined_baseline_df$Scenario)
combined_baseline_df$Run <- as.factor(combined_baseline_df$Run)

# Create box plots for each variable with Dark Mint palette
boxplot_ND <- ggplot(combined_baseline_df, aes(x = Scenario, y = ND, fill = Scenario)) +
  geom_boxplot(outlier.size=0.5) +
  labs(title = "A) Number of deaths",
       x = "",
       y = "Deaths") +
  scale_fill_manual(values = my_palette) +
  theme_minimal() +
  theme(legend.position = "none")

boxplot_cum_mort_rate <- ggplot(combined_baseline_df, aes(x = Scenario, y = cum_mort_rate, fill = Scenario)) +
  geom_boxplot(outlier.size=0.5) +
  labs(title = "B) CFR",
       x = "",
       y = "%") +
  scale_fill_manual(values = my_palette) +
  theme_minimal() +
  theme(legend.position = "none")

boxplot_WB <- ggplot(combined_baseline_df, aes(x = Scenario, y = WB*100, fill = Scenario)) +
  geom_boxplot(outlier.size=0.5) +
  labs(title = "C) First-line resistance",
       x = "",
       y = "%") +
  scale_fill_manual(values = my_palette) +
  theme_minimal() +
  theme(legend.position = "none")

boxplot_WC <- ggplot(combined_baseline_df, aes(x = Scenario, y = WC*100, fill = Scenario)) +
  geom_boxplot(outlier.size=0.5) +
  labs(title = "D) Second-line resistance",
       x = "",
       y = "%") +
  scale_fill_manual(values = my_palette) +
  theme_minimal() +
  theme(legend.position = "none")

boxplot_TAA <- ggplot(combined_baseline_df, aes(x = Scenario, y = TAA, fill = Scenario)) +
  geom_boxplot(outlier.size=0.5) +
  labs(title = "E) First-line antibiotics",
       x = "Scenario",
       y = "Days") +
  scale_fill_manual(values = my_palette) +
  theme_minimal() +
  theme(legend.position = "none")

boxplot_TBB <- ggplot(combined_baseline_df, aes(x = Scenario, y = TBB, fill = Scenario)) +
  geom_boxplot(outlier.size=0.5) +
  labs(title = "F) Second-line antibiotics",
       x = "Scenario",
       y = "Days") +
  scale_fill_manual(values = my_palette) +
  theme_minimal() +
  theme(legend.position = "none")

boxplot_TCC <- ggplot(combined_baseline_df, aes(x = Scenario, y = TCC, fill = Scenario)) +
  geom_boxplot(outlier.size=0.5) +
  labs(title = "G) Third-line antibiotics",
       x = "Scenario",
       y = "Days") +
  scale_fill_manual(values = my_palette) +
  theme_minimal() +
  theme(legend.position = "none")

# Add a common legend in the 8th slot
boxplot_legend <- ggplot(combined_baseline_df, aes(x = Scenario, y = TCC, fill = Scenario)) +
  geom_boxplot(outlier.size=0.5) +
  labs(title = "Box Plot for TCC",
       x = "Scenario",
       y = "Values") +
  scale_fill_manual(values = my_palette) +
  theme_minimal()
common_legend <- get_legend(boxplot_legend)  # You can use any of the boxplots for the legend since they all share the same legend


# Combine the box plots into a single figure with a common legend
combined_boxplots <- plot_grid(boxplot_ND, boxplot_cum_mort_rate, boxplot_WB,
                               boxplot_WC, boxplot_TAA, boxplot_TBB,
                               boxplot_TCC, common_legend, nrow = 2, rel_widths = c(1, 1, 1, 1),
                               rel_heights = c(1,1.1))


