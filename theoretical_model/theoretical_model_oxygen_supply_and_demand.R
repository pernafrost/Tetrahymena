## List of variables


## List of functions:
# V <- calculate_cell_volume(rs, rl)
# CS <- calculate_cell_cross_section(rs)
# mu <- calculate_water_viscosity(temp)
# Vml <- convert_um3_to_ml(Vum3)
# StP <- calculate_stokes_power(rs, rl, U, mu, useModel)
# Sc <- calculate_stokes_coefficient(eta_c, rs, rl, mu, useModel)
# U <- calculate_speed_given_power(rs, rl, StP, mu, useModel)
# B <- scale_metabolic_rate_with_temperature(BRef, tempRef, tempH, temp, E, ED)
# Bms <- scale_metabolic_rate_with_body_size(Bms0, M0, M, betaExponent)
# B <- scale_metabolic_rate_with_body_size_and_temperature(B0, M0, M, tempRef, tempH, temp, E, ED, betaExponent) 
# feedingRate <- functional_response(eta, f, U)
# metabolicRateW <- respiration_rate_to_W(oxygenConsumption, dPredator)
# oxygenConsumption <- W_to_respiration_rate(metabolicRateW, dPredator) 
# henryConstant <- calculate_Henry_constant_at_T(temp, gas)
# gasConcentration <- calculate_gas_concentration_from_partial_pressure(temp, gas, partialPressure)
# partialPressure <- calculate_partial_pressure_from_gas_concentration(temp, gas, gasConcentration) 
# D <- scale_diffusion_coefficient_with_temperature(D0, tempRef, temp)
# Uopt <- calculate_optimal_speed(s, eta, fMax)
# calculate_feeding_given_speed <- function(U, eta, fMax)




rm(list=ls()) # clean memory

zeroCelsius <- 273.15 # for the conversion between Kelvin and Celsius
verbose = FALSE # This tells the programme to print some information for each function


#################################################################################################################
# FUNCTION calculate_cell_volume
# calculates the volume of a cell given the two radii
# input:
# rs - the small radius of the cell (in a unit of length, e.g. micrometres)
# rl - the long radius of the cell (in a unit of length, e.g. in micrometres)
# output
# V - the volume of the cell (in the same unit of length to the third power, e.g. micrometres^3)
calculate_cell_volume <- function(rs, rl) {
  V <- 4/3 * pi * rl * rs^2
  return(V)
}






#################################################################################################################
# FUNCTION calculate_cell_cross_section
# calculates the cross_section of a cell given the small radius
# input:
# rs - the small radius of the cell (in a unit of length, e.g. micrometres)
# output
# CS - the cross-section of the cell (in the same unit of length squared, e.g. micrometres^2)
calculate_cell_cross_section <- function(rs) {
  CS <- pi * rs^2
  return(CS)
}




#################################################################################################################
# FUNCTION calculate_water_viscosity
# This function calculates the viscosity of water at a given temperature
# input: temperature in Kelvin degrees
# output: the viscosity of water in Kg/(m s)
# This function is based on parameters for water (not sure about salinity)
# water viscosity increases at low temperature
calculate_water_viscosity <- function(temp) {
  # zeroCelsius <- 273.15
  # verbose = TRUE
  # Parameters of water viscosity and its dependence on temperature
  viscA <- 1.856e-14# parameter A of viscosity calculation in Pa s CONSTANT
  viscB <- 4209# parameter B of viscosity calculation in K CONSTANT
  viscC <- 0.04527 # parameter C of viscosity calculation in K^-1 CONSTANT
  viscD <- -3.376e-5 # parameter D of viscosity calculation in K^-2 CONSTANT
  mu <- viscA * exp(viscB/temp + viscC*temp + viscD*temp^2) # dynamic viscosity of water, approx. 10^-3 Kg/(m s) Andreade equation with extra parameters
  if (verbose)
  {
    print(paste("Viscosity at temperature ", temp - zeroCelsius, "C is ", mu, "Kg/(m s)"))
    print("(Viscosity should be around 10^-3 Kg/(m s))")
  }
  return(mu)
}







#################################################################################################################
# FUNCTION convert_um3_to_ml
# converts cubic micrometre to ml
# input a number in units of um^3
# output: the corresponding volume in ml
convert_um3_to_ml <- function(Vum3){
  Vml <- Vum3 * 1e-12
  return (Vml)
}






#################################################################################################################
# FUNCTION calculate_stokes_power
# this function calculates the stokes power of a spheroid moving in a viscous fluid
# input:
# rs - small radius of the object in micrometres (important, this is the radius, not the cross section
# which means -for instance- that if tetrahymena has a cross-section of 30 micrometres, the radius should
# be 15 micrometres here)
# rl - long radius of the object in micrometres
# U - speed of the object in micrometres / s
# mu - dynamic viscosity of water at the test temperature (in Kg/(m*s) or Pa * s)
# useModel -   # We can either approximate the organism as a sphere with radius rs or as a prolate spheroid with long radius rl and short radius rs
# possible values for useModel are "sphere" and "spheroid"; if no value is specified, the spheroid model is used
# output:
# stokes power in Watts
calculate_stokes_power <- function(rs, rl, U, mu, useModel) { 
  
  # If the model is not specified, use the prolate spheroid model
  if (missing(useModel)){
    useModel = "spheroid" # "spheroid" # decide if we use the sphere or the spheroid model for the drag calculation
  }
  
  # convert the units from micrometres and micrometres per second to international system mks units
  rl <- rl * 10^-6 # now rl is in m
  rs <- rs * 10^-6 # now rs is in m
  U <- U * 10^-6 # now U is in m/s
  # This conversion is useful so that the resulting units for power are Watts
  
  if (verbose)
  {
    # run a quick preliminary test that we are in viscous regime by estimating the Reynolds number
    rho <- 1000 # rho is the density of water in kg / m^3
    nu <- mu / rho # nu is the kinematic viscosity of water (m^2/s)
    Re <- rs*2*U/nu # the Reynolds number is L*U/nu, where L is the characteristic length. Here I use the diametre of the object
    print(paste("With a diametre of ", rs*2, " m and a speed of ", U, "m/s, the Reynolds number is ", Re))
  }
  
  
  
  switch(useModel,
         sphere = {
           CD <- 6 * pi * rs # This is the dimensionless drag coefficient
           StP <- CD * mu * U^2 # Stokes power calculated for one cell
           if (verbose)
           {
             print(paste("using sphere model gives Dc=", CD, "and Stokes Power ", StP, "W / cell"))
           }
         },
         spheroid = {
           tau = rl/(rl^2 - rs^2)^0.5 # This is the parameter tau for the accurate calculation of drag for a spheroid
           CD <- 8 * pi * rs / (sqrt(tau^2 - 1) * ((tau^2 + 1)/2 * log((tau + 1)/(tau - 1), base=exp(1)) - tau)) #This is the dimensionless drag coefficient for an oblate spheroid
           StP <- CD * mu * U^2 # Stokes power calculated for one cell
           if (verbose)
           {
             print(paste("using spheroid model gives CD= ", CD, " and Stokes Power= ", StP,  "W / cell"))
           }
         },
         simplifiedSpheroid = { # reference is Clift's book, table 4.1
           E <- rl/rs # elongation of the spheroid
           CD <- 1.2 * pi * rs * (4+E) #This is the dimensionless drag coefficient for an oblate spheroid
           StP <- CD * mu * U^2 # Stokes coefficient calculated for one cell
           if (verbose)
           {
             print(paste("using simplified spheroid model gives CD= ", CD, " and Stokes Power= ", StP,  "W / cell"))
           }
         })
  return(StP)
}








#################################################################################################################
# FUNCTION calculate_stokes_coefficient
# this function calculates the stokes power of a spheroid moving in a viscous fluid
# input:
# eta_c is the efficiency of conversion of metabolic energy to movement
# rs - small radius of the object in micrometres (important, this is the radius, not the cross section
# which means -for instance- that if tetrahymena has a cross-section of 30 micrometres, the radius should
# be 15 micrometres here)
# rl - long radius of the object in micrometres
# mu - dynamic viscosity of water at the test temperature (in Kg/(m*s) or Pa * s)
# useModel -   # We can either approximate the organism as a sphere with radius rs or as a prolate spheroid with long radius rl and short radius rs
# possible values for useModel are "sphere", "spheroid", and "simplifiedSpheroid"; if no value is specified, the spheroid model is used
# output:
# Sc - stokes coefficient (efficiency factor x drag coefficient x viscosity)


calculate_stokes_coefficient <- function(eta_c, rs, rl, mu, useModel) { 
  Sc <- NA
  
  # If the model is not specified, use the prolate spheroid model
  if (missing(useModel)){
    useModel = "spheroid" # "spheroid" # decide if we use the sphere or the spheroid model for the drag calculation
  }
  
  # convert the units from micrometres and micrometres per second to international system mks units
  rl <- rl * 10^-6 # now rl is in m
  rs <- rs * 10^-6 # now rs is in m
  # This conversion is useful so that the resulting units for power are Watts
  
  # An approximation that always works would be:
  # CD = 1.2*pi*rs*(4 + E)
  
  if (rs == rl)
  {
    useModel <- "sphere"
  }
  
  switch(useModel,
         sphere = {
           CD <- 6 * pi * rs # This is the dimensionless drag coefficient
           Sc <- CD * mu / eta_c # Stokes coefficient calculated for one cell
           if (verbose)
           {
             print(paste("using sphere model gives Dc=", CD, "and Stokes coefficient ", Sc))
           }
         },
         spheroid = { # reference is Katsu-Kimura et al. 2009
           tau = rl/(rl^2 - rs^2)^0.5 # This is the parameter tau for the accurate calculation of drag for a spheroid
           CD <- 8 * pi * rs / (sqrt(tau^2 - 1) * ((tau^2 + 1)/2 * log((tau + 1)/(tau - 1), base=exp(1)) - tau)) #This is the dimensionless drag coefficient for an oblate spheroid
           Sc <- CD * mu / eta_c # Stokes coefficient calculated for one cell
           if (verbose)
           {
             print(paste("using spheroid model gives CD= ", CD, " and Stokes coefficient= ", Sc))
           }
         },
         simplifiedSpheroid = { # reference is Clift's book, table 4.1
           E <- rl/rs # elongation of the spheroid
           CD <- 1.2 * pi * rs * (4+E) #This is the dimensionless drag coefficient for an oblate spheroid
           Sc <- CD * mu / eta_c # Stokes coefficient calculated for one cell
           if (verbose)
           {
             print(paste("using simplified spheroid model gives CD= ", CD, " and Stokes coefficient= ", Sc))
           }
         })
  return(Sc)
}




#################################################################################################################
# FUNCTION calculate_speed_given_power
# this function is similar to the function calculate_stokes_power. The main difference is that here we take
# a value of power in input and we calculate the speed that the object would reach given that value of power
# input:
# rs - small radius of the object in micrometres (important, this is the radius, not the cross section
# which means -for instance- that if tetrahymena has a cross-section of 30 micrometres, the radius should
# be 15 micrometres here)
# rl - long radius of the object in micrometres
# StP - propulsive power produced by the swimmer
# mu - dynamic viscosity of water at the test temperature (in Kg/(m*s) or Pa * s)
# useModel -   # We can either approximate the organism as a sphere with radius rs or as a prolate spheroid with long radius rl and short radius rs
# possible values for useModel are "sphere" and "spheroid"; if no value is specified, the spheroid model is used
# output:
# speed U in micrometres / second
calculate_speed_given_power <- function(rs, rl, StP, mu, useModel) {
  
  # If the model is not specified, use the prolate spheroid model
  if (missing(useModel)){
    useModel = "spheroid" # "spheroid" # decide if we use the sphere or the spheroid model for the drag calculation
  }
  
  # convert the units from micrometres and micrometres per second to international system mks units
  rl <- rl * 10^-6 # now rl is in m
  rs <- rs * 10^-6 # now rs is in m
  # This conversion is useful so that the resulting units for power are Watts
  
  switch(useModel,
         sphere = {
           CD <- 6 * pi * rs # This is the dimensionless drag coefficient
           # StP <- CD * mu * U^2 # Stokes power calculated for one cell
           U <- sqrt(StP/CD / mu) # speed calculated based on stokes power and drag
           if (verbose)
           {
             print(paste("using sphere model gives Dc=", Dc, "and Stokes Power ", StP, "W / cell"))
           }
         },
         spheroid = {
           tau = rl/(rl^2 - rs^2)^0.5 # This is the parameter tau for the accurate calculation of drag for a spheroid
           CD <- 8 * pi * rs / (sqrt(tau^2 - 1) * ((tau^2 + 1)/2 * log((tau + 1)/(tau - 1), base=exp(1)) - tau)) #This is the dimensionless drag coefficient for an oblate spheroid
           # StP <- CD * mu * U^2 # Stokes power calculated for one cell
           U <- sqrt(StP/CD / mu) # speed calculated based on stokes power and drag
           if (verbose)
           {
             print(paste("using spheroid model gives CD= ", CD, " and Stokes Power= ", StP,  "W / cell"))
           }
         })
  
  U <- U * 10^6 # now U is in micrometres / second (it was in metres / second)
  return(U)
}
















#################################################################################################################
# FUNCTION scale_metabolic_rate_with_temperature
# input:
# BRef is the metabolic rate at the reference temperature (e.g. in Watts)
# tempRef is the reference temperature in degrees Kelvin
# tempH temperature of half inactivation of enzimatic units (K)
# temp is the temperature at which we want to predict the metabolic rate, also in degrees kelvin
# E is the activation energy (for instance, expressed in electronvolt eV)
# ED is the deactivation energy (for instance, expressed in electronvolt eV)
# kB is the Boltzmann constant (for instance, expressed in electronvolt / K)
# Notice that the activation and deactivation energy could also be expressed in Joules if more convenient
# but in this case the Boltzmann constant should also be in J/K
scale_metabolic_rate_with_temperature <- function(BRef, tempRef, tempH, temp, E, ED) {
  # kB <- 1.380649e-23 # Boltzmann constant (J/K) CONSTANT
  kB <- 8.617333262145e-5 # Boltzmann constant (eV/K) CONSTANT
  if (verbose)
  {
    print("running function scale_metabolic_rate_with_temperature")
    print(paste("Using Boltzmann constant kB=", kB, "eV/K (check that the input activation and deactivation energy were also in eV)"))
  }
  # Here is the scaling of the standard metabolic rate with temperature (Schoolfield equation)
  schoolfieldActivation <- exp(-E / kB *(1/temp - 1/tempRef))
  schoolfieldDeactivation <- 1 + exp(ED / kB *(1/tempH - 1/temp))
  B <- BRef * schoolfieldActivation / schoolfieldDeactivation # Metabolic rate scaled with temperature
  # the units of B are the same as the units of BRef (could be Watt, could be respiration...)
  return(B)
}









#################################################################################################################
# FUNCTION scale_metabolic_rate_with_body_size
# input:
# Bms0 is the metabolic rate for the reference body size M0
# M0 is the reference body size (e.g. in grams, cubic metres, etc.)
# M is the body size for which we want to estimate the metabolic rate
# betaExponent is the allometric scaling exponent (typically 3/4 or 1)
# output:
# Bms is the scaled metabolic rate
scale_metabolic_rate_with_body_size <- function(Bms0, M0, M, betaExponent) {
  if (missing(betaExponent))
  {  
    betaExponent <- 3/4 # betaExponent is the allometric scaling exponent (typically 3/4 or 1)
  }
  if (verbose)
  {
    print(paste("Using an allometric exponent of beta=", betaExponent))
  }
  Bms <- Bms0 * (M/M0)^betaExponent
  return(Bms)
}










#################################################################################################################
# FUNCTION scale_metabolic_rate_with_body_size_and_temperature
# input:
# B0 is the metabolic rate for the reference body size M0 and the reference temperature tempRef
# M0 is the reference body size (e.g. in grams, cubic metres, etc.)
# M is the body size for which we want to estimate the metabolic rate
# tempRef is the reference temperature in degrees Kelvin
# tempH temperature of half inactivation of enzimatic units (K)
# temp is the temperature at which we want to predict the metabolic rate, also in degrees kelvin
# E is the activation energy (for instance, expressed in electronvolt eV)
# ED is the deactivation energy (for instance, expressed in electronvolt eV)
# betaExponent is the allometric scaling exponent (typically 3/4 or 1)
# Notice that the activation and deactivation energy could also be expressed in Joules if more convenient
# but in this case the Boltzmann constant should also be in J/K
# output:
# B is the metabolic rate for the new body size and new temperature (in the same units given in input)
scale_metabolic_rate_with_body_size_and_temperature <- function(B0, M0, M, tempRef, tempH, temp, E, ED, betaExponent) 
{
  Bms <- scale_metabolic_rate_with_body_size(B0, M0, M, betaExponent)
  B <- scale_metabolic_rate_with_temperature(Bms, tempRef, tempH, temp, E, ED)
  return(B)
}













#################################################################################################################
# FUNCTION functional_response
# This function calculates the energy consumed per unit of time as a function of
# speed of predator and of functional response parameters eta and f
# input:
# f - maximum feeding rate in W (= 1/h)
# U - speed in um/s
# s - stokes coefficient that integrates all the costs of swimming
# except for the swimming speed itself
# eta - half saturation speed such that eta = 1/(CS x h x d) has the units of a speed
# with CS an area (cross-section)
# h an handling time (time to process a unit energy, with units of W^-1)
# d is the energy density of the medium energy/volume
# finally f has the units of a power and is the maximum feeding rate

# Note that here we are assuming density of prey d0 = 1
# the handling time
# output:
# Watts of food consumed per unit time
functional_response <- function(eta, f, U) {
  feedingRate <- f * U / (eta + U)
  return(feedingRate)
}











#################################################################################################################
# FUNCTION respiration_rate_to_W
# Converts the respiration rate of the culture or of the cell to metabolic rate (in Watts)
# input:
# respiration rate of the organism or of the culture in micromol / L / min
# dPredator is the density of the predator (tetrahymena) in cells / ml
# output:
# metabolic power consumption (W) of one single cell
# To make a test let's use some values from Katsu-Kimura, who
# report an oxygen consumption for a single cell of paramecium
# of 0.19-4.4 * 10^-9 L/cell/hour and transform the number to get
# 0.38-8.8 * 10^-5 J/cell/hour. I convert L to micromol by considering
# 1L * 1atm = n * R * T with T = 298 and R = 0.082057366080960
# oxygenConsumption_in_mol_per_min <- 0.19e-9 * 1/(293*0.082057366080960)*60
# [1] 4.741543e-10 mol/L/min
# oxygenConsumption <- oxygenConsumption_in_mol_per_min*10^6
# > respiration_rate_to_W(oxygenConsumption)
# [1] 3.781983e-06 # W/cell
# > respiration_rate_to_W(oxygenConsumption) * 3600
# [1] 0.01361514 # J/cell/h
respiration_rate_to_W <- function(oxygenConsumption, dPredator) {
  if (missing(dPredator))
  {
    dPredator <- 1e-3  # if dPredator is omitted, we pretend that there is only one cell per L
    # so that the function calculates the conversion per cell, rather than per L
    if (verbose)
    {
      print(paste("In function respiration_rate_to_W, calculating conversion of respiration rate per cell to W /cell"))
    }
  }
  # Convert rate of respiration in umol /min / L to rate of energy consumption in W/cell using the following assumptions:
  # Aerobic respiration goes as 1G + 6O2 = 6CO2 + 6H2O (1G is one mol of glucose) and produces 2871458 J (about 686 kCal)
  JoulesPerMol <- 2871458/6 # (J/mol) there are 2871458 J in one mol of glucose and 2871458/6 J in one mol of produced O2
  JoulesPerMicroMol <- JoulesPerMol * 1e-6 # (J/umol)
  metabolicRateJoulesPerMinPerL <- oxygenConsumption * JoulesPerMicroMol # Joules / L / min
  # I divide this metabolic rate by the number of cells of predator per litre of suspension
  # to get the metabolic rate per cell
  metabolicRatePerCellPerMin <- metabolicRateJoulesPerMinPerL / (dPredator * 1e3) #  Joules / min /cell
  metabolicRateW <- metabolicRatePerCellPerMin/60 # (Joules / s / cell) metabolicRatePerCellPerMin was in Joules per min, 
  return(metabolicRateW)
}



#################################################################################################################
# FUNCTION W_to_respiration_rate
# Converts the metabolic rate of one cell (in Watts) to the predicted respiration rate of the culture
# input:
# metabolic power consumption (W) of one single cell
# dPredator is the density of the predator (tetrahymena) in cells / ml
# output:
# respiration rate of the population in micromol / L / min
# optional dPredator is the density of the predator (tetrahymena) in cells / ml
W_to_respiration_rate <- function(metabolicRateW, dPredator) 
{
  if (missing(dPredator))
  {
    dPredator <- 1e-3  # if dPredator is omitted, we pretend that there is only one cell per L
    # so that the function calculates the conversion per cell, rather than per L
    if (verbose)
    {
      print(paste("In function W_to_respiration_rate, calculating conversion of metabolic rate in W /cell to respiration in umol / min /cell"))
    }
  }
  # Convert rate of respiration in umol /min / L to rate of energy consumption in W/cell using the following assumptions:
  # Aerobic respiration goes as 1G + 6O2 = 6CO2 + 6H2O (1G is one mol of glucose) and produces 2871458 J (about 686 kCal)
  JoulesPerMol <- 2871458/6 # (J/mol) there are 2871458 J in one mol of glucose and 2871458/6 J in one mol of produced O2
  JoulesPerMicroMol <- JoulesPerMol * 1e-6 # (J/umol)
  
  
  metabolicRatePerCellPerMin <- metabolicRateW*60 # (Joules / min / cell) metabolicRatePerCellPerMin is in Joules per min,
  # but the metabolic rate was in Joules per second
  
  metabolicRateJoulesPerMinPerL <- metabolicRatePerCellPerMin * (dPredator * 1e3) # (Joules /min / L) I multiply the individual metabolic
  # rate times the number of cells in one L of suspension
  
  oxygenConsumption <- metabolicRateJoulesPerMinPerL / JoulesPerMicroMol # micromol / min / L
  return(oxygenConsumption)
}






#################################################################################################################
# FUNCTION calculate_Henry_constant_at_T
# This function calculates the Henry constant for gas solubility
# depending on temperature
# input:
# temp - temperature in Kelvin degrees
# gas - a string that can be "O2", "CO2", "N2"
# note other parameters could be water salinity
# output:
# the value of the Henry constant
calculate_Henry_constant_at_T <- function(temp, gas) {
  
  tempRefHC = 298.15 # the reference temperature (K) at which
  # Henry constant was calculated
  
  vantHoffParameter=1700 # for the scaling of Henry's constant with temperature
  henryConstantAtTRef = 0.001321528 # Henry's law constant at 298.15 K in mol/L/atm
  switch(gas, O2 ={
    vantHoffParameter=1700
    henryConstantAtTRef = 0.001321528
  },
  CO2={
    vantHoffParameter=2400
    henryConstantAtTRef = 7.8e-4
  },
  N2={
    vantHoffParameter=1300
    henryConstantAtTRef=6.1e-4
  })
  
  henryConstant = henryConstantAtTRef * exp(vantHoffParameter*(1/temp - 1/tempRefHC))
  # https://en.wikipedia.org/wiki/Henry%27s_law#Temperature_dependence
  
  return(henryConstant)
}








#################################################################################################################
# FUNCTION calculate_gas_concentration_from_partial_pressure
# This function calculates the concentration of a gas
# in the aqueous phase from its partial pressure in the gas phase
# assuming that the two phases are at equilibrium
# input:
# temp - temperature in Kelvin degrees
# gas - a string that can be "O2", "CO2", "N2"
# partialPressure (partial pressure of the gas in atmospheres)
# partialPressure is the partial pressure of the gas in the gas phase
# note other parameters could be water salinity
# output:
# gasConcentration - the concentration of the gas in mol/L
calculate_gas_concentration_from_partial_pressure <- function(temp=298.15, gas="O2", partialPressure=0.20948) 
{
  henryConstant = calculate_Henry_constant_at_T(temp, gas) 
  # get the value of Henry constant at the temperature of interest
  gasConcentration = partialPressure * henryConstant
  return(gasConcentration)
}





#################################################################################################################
# FUNCTION calculate_partial_pressure_from_gas_concentration
# This function calculates the partial_pressure of a gas
# in the gas phase from its concentration in the aqueous phase
# assuming that the two phases are at equilibrium
# input:
# gasConcentration - the concentration of the gas in mol/L
# temp - temperature in Kelvin degrees
# gas - a string that can be "O2", "CO2", "N2"
# note other parameters could be water salinity
# output:
# partialPressure (partial pressure of the gas in atmospheres)
# partialPressure is the partial pressure of the gas in the gas phase

calculate_partial_pressure_from_gas_concentration <- function(temp=298.15, gas="O2", gasConcentration=374.3128e-6) 
{
  henryConstant = calculate_Henry_constant_at_T(temp, gas) 
  # get the value of Henry constant at the temperature of interest
  partialPressure = gasConcentration / henryConstant
  return(partialPressure)
  # the same in mol/L in the air would be
  # based on PV=nRT and R=0.082
  #partialPressure*1/0.082/temp molO2/L
  # partialPressure*1/0.082/temp*10^6 micromolO2/L in the air
  # example
  # p <- calculate_partial_pressure_from_gas_concentration(298.15, "O2", 300e-6)
  #> p
  # [1] 0.2270099
  # > p*1/0.082/temp*10^6
  # [1] 9285.306
}








#################################################################################################################
# FUNCTION scale_diffusion_coefficient_with_temperature
# This function scales the diffusion coefficient with temperature
# input:
# temp (K) - temperature in Kelvin degrees for which we want the diffusion coefficient
# tempRef (K) - temperature in Kelvin degrees at which the diffusion coefficient was measured
# D0 (m^2/s) - diffusion coefficient at the temperature tempRef
# output:
# D - the diffusion coefficient at the temperature temp (m^2/s)
# This is obtained from the formula D = D0 * temp/tempRef * muT0 / muT
# where muT0 is the viscosity of water at tempRef and 
# muT is the viscosity at T
scale_diffusion_coefficient_with_temperature <- function(D0 = 2.10e-9, tempRef = 293.15, temp)
{
  if (min(temp) < 273.15 | min(tempRef) < 273.15)
  {
    paste("Are you sure that the temperature is in Kelvin degrees?")
  }
  muT0 <- calculate_water_viscosity(tempRef)
  muT <- calculate_water_viscosity(temp)
  D <- D0 * temp/tempRef * muT0 / muT
  return(D)
}



#################################################################################################################
# FUNCTION calculate_optimal_speed
# This function calculates the optimal swimming speed for maximizing foraging
# input:
# s - stokes coefficient that integrates all the costs of swimming
# except for the swimming speed itself
# eta - half saturation density such that eta = 1/(CS x h x d) has the units of a speed
# with CS an area (cross-section)
# h an handling time (time to process a unit energy, with units of W^-1)
# d is the energy density of the medium energy/volume
# finally f has the units of a power and is the maximum feeding rate
# output:
# Uopt - the optimal speed
calculate_optimal_speed <- function(s, eta, fMax)
{
  eta <- eta*10^-6 # convert eta to m/s from um/s
  f <- fMax * 10^-9 # convert the feeding rate from uW to W
  Uopt <- ((eta*sqrt(f*(8*eta^2*s+27*f)))/(4*3^(3/2)*s)+((3*f*eta)/(2*s)+2*eta^3)/6-(8*eta^3)/27)^(1/3)+eta^2/(9*((eta*sqrt(f*(8*eta^2*s+27*f)))/(4*3^(3/2)*s)+((3*f*eta)/(2*s)+2*eta^3)/6-(8*eta^3)/27)^(1/3))-(2*eta)/3
  
  # if at high Reynolds numbers, the equation is different:
  # Uopt <- sqrt(+(4*sqrt((eta*f)/s))/sqrt(3)+eta^2)/2 - eta/2
  # however, I think we must also be careful when calculating the
  # coefficient s, as it would need to be different (the units are
  # proably different because we now multiply by speed one more time)
  # To be checked
  
  Uopt <- Uopt * 10^6 # converts back speed to um/s
  # all these conversions are because the Stokes coefficient is
  # in international units
  return(Uopt)
}









#################################################################################################################
# FUNCTION calculate_feeding_given_speed
# This function calculates the feeding rate given a swimming speed
# input:
# U - speed
# eta - half saturation density such that eta = 1/(CS x h x d0) has the units of a speed
# with CS an area (cross-section)
# h an handling time (time to process a unit energy, with units of W^-1)
# d is the energy density of the medium energy/volume
# finally f has the units of a power
# eta is already normalised to a given metabolic rate and energy density
# and the same for fMax
# eta <- eta0 * d0/d * B0/B
# fMax <- fMax0 * B0/B
# output:
# Uopt - fRate also a power
calculate_feeding_given_speed <- function(U, eta, fMax)
{
  fRate <- fMax*U/(U + eta)
  return(fRate)
}










#################################################################################################################
# Load libraries
library(ggplot2)
library(Cairo)
library(viridis) # colour blind friendly palette, works in B&W also

setwd("~/Tetrahymena/theoretical_model")


if(!is.null(dev.list())) dev.off()


adaptedTemperature <- 20 # degrees Celsius
Bmeasured <- 0.6936737 # nW
Ea <- 0.8898539 # activation energy
C0 <- 10^-4 #(0.00001 is 10 umol/L) # Oxygen concentration
# according to this reference:
# Wilson, D.F., Rumsey, W.L., Green, T.J. and Vanderkooi, J.M., 1988. The oxygen dependence of mitochondrial oxidative phosphorylation measured by a new optical method for measuring oxygen concentration. Journal of Biological Chemistry, 263(6), pp.2712-2718.
# oxygen concentrations above 0.00001 M or 10 uM would be not limiting anymore at the level of mitochondria


logV_ref <- 4.25
V_ref <- 10^logV_ref # a reference volume in the range of those observed
rs <- (3/8/pi*V_ref)^(1/3) # 25/2 # um
rl <- rs * 2 # 50/2 # um

r_eq <- (rs*rs*rl)^(1/3) # for the purpose of this script I
# approximate the cell to a sphere, because it is easier to
# model diffusion in this way
r_eq_test <- r_eq
# The volume of this cell would be 4/3*pi*r_eq^3

testedTemperature <- seq(10, 40, by=0.5)# c(10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40)


# scale_metabolic_rate_with_temperature(Bmeasured, adaptedTemperature + 273.15, tempH=31+273.15, temp=testedTemperature + 273.15, Ea, 7)
B <- scale_metabolic_rate_with_body_size_and_temperature(B0=Bmeasured, M0=4/3*pi*r_eq^3, M=4/3*pi*r_eq_test^3, tempRef=adaptedTemperature + 273.15, tempH=31+273.15, temp=testedTemperature + 273.15, E=Ea, ED=7, betaExponent=1)
# B is the metabolic rate in nW

C <- calculate_gas_concentration_from_partial_pressure(testedTemperature + 273.15, "O2", 0.20948)
# C is the concentration of O2 in saturated water, in mol/L


D <- scale_diffusion_coefficient_with_temperature(D0 = 2.10e-9, tempRef = 293.15, testedTemperature + 273.15)
# D is the diffusion coefficient of oxygen in water expressed in m^2/s

Oxygen_supply <- 4*pi*D*((C - C0)*10^9)*(r_eq*10^-6)*60 # this is the current (flux x cell surface area) of Oxygen into the cell
# This is the oxygen supply in umol/cell/min
Energy_supply_nW <- respiration_rate_to_W(Oxygen_supply)*10^9

df1 <- data.frame(temp = testedTemperature, B = B, C = C, D = D, Energy_supply_nw = Energy_supply_nW)

labels = c("Metabolic rate", "O2 concentration", "Diffusion coefficient", "Energy supply")
colours = c("#000000", "#0000BB", "#118800", "#CC3300")

g0 <- ggplot() +
  geom_line(data = df1, size=2, linetype = 1, aes(x=temp, y = B, color=labels[1])) +
  annotate("text", size=5, x=10.5, y=1.75, label= paste("E=", round(Ea,2), "eV", sep=""), hjust = 0, parse=F) +
  # annotate("text", size=5, x=10.5, y=3.15, label= paste("E=", Ea, "eV", sep=""), hjust = 0, parse=F) +
  scale_x_continuous(name=expression(paste("Temperature °C "))) +
  scale_y_continuous(name="Metabolic rate (nW)", limits=c(0,2), sec.axis = sec_axis( trans=~.*W_to_respiration_rate(1e-9), name=expression(paste("Metabolic rate (",  mu, "mol[O2]/cell/min", ")")))) +
  theme_classic(base_size = 15) + 
  scale_color_manual(name="Contribution", breaks=labels, values=colours) +
  theme(legend.position = "none") # + # This removes the legend above
g0


g1 <- ggplot() +
  geom_line(data = df1, size=2, linetype = 1, aes(x=temp, y = C*10^6, color=labels[2])) +
  scale_x_continuous(name=expression(paste("Temperature °C "))) +
  scale_y_continuous(name=expression(paste("Oxygen (", mu, "mol/L)")), limits=c(0,400)) +
  theme_classic(base_size = 15) + 
  scale_color_manual(name="Contribution", breaks=labels, values=colours) +
  theme(legend.position = "none") # + # This removes the legend above
g1

g2 <- ggplot() +
  geom_line(data = df1, size=2, linetype = 1, aes(x=temp, y = C*32*10^3, color=labels[2])) +
  scale_x_continuous(name=expression(paste("Temperature °C "))) +
  scale_y_continuous(name=expression(paste("Oxygen (mg/L)")), limits=c(0,15)) +
  theme_classic(base_size = 15) + 
  scale_color_manual(name="Contribution", breaks=labels, values=colours) +
  theme(legend.position = "none") # + # This removes the legend above
g2

g3 <- ggplot() +
  geom_line(data = df1, size=2, linetype = 1, aes(x=temp, y = D, color=labels[3])) +
  scale_x_continuous(name=expression(paste("Temperature °C "))) +
  scale_y_continuous(name=expression(paste("Diffusion coeff. (m"^2, "/s)")), limits=c(0, 4e-9)) +
  theme_classic(base_size = 15) + 
  scale_color_manual(name="Contribution", breaks=labels, values=colours) +
  theme(legend.position = "none") # + # This removes the legend above
g3




g4 <- ggplot() +
  geom_line(data = df1, size=2, linetype = 1, aes(x=temp, y = Energy_supply_nW, color=labels[4])) +
  scale_x_continuous(name=expression(paste("Temperature °C "))) +
  scale_y_continuous(name="Energy supply (nW)", limits=c(0,45), sec.axis = sec_axis( trans=~.*W_to_respiration_rate(1e-9), name=expression(paste("Oxygen supply (",  mu, "mol[O2]/cell/min", ")")))) +
  theme_classic(base_size = 15) + 
  scale_color_manual(name="Contribution", breaks=labels, values=colours) +
  theme(legend.position = "none") # + # This removes the legend above
g4
  

currentTitle="oxygen_supply_demand"

# save plot g0 with metabolic rate
ggsave(file=paste(currentTitle, "_MR_logvol", logV_ref, ".png", sep=""), plot=g0, dpi = 600, width = 15, height = 15, units = "cm")
ggsave(file=paste(currentTitle, "_MR_logvol", logV_ref, ".eps", sep=""), plot=g0, device="eps", dpi = 600, width = 15, height = 15, units = "cm")
ggsave(file=paste(currentTitle, "_MR_logvol", logV_ref, ".pdf", sep=""), plot=g0, device=cairo_pdf, dpi = 600, width = 15, height = 15, units = "cm")

# save plot g1 for oxygen solubility
ggsave(file=paste(currentTitle, "_O2_solubility.png", sep=""), plot=g1, dpi = 600, width = 15, height = 15, units = "cm")
ggsave(file=paste(currentTitle, "_O2_solubility.eps", sep=""), plot=g1, device="eps", dpi = 600, width = 15, height = 15, units = "cm")
ggsave(file=paste(currentTitle, "_O2_solubility.pdf", sep=""), plot=g1, device=cairo_pdf, dpi = 600, width = 15, height = 15, units = "cm")

# save plot g3 for the diffusion coefficient of oxygen in water
ggsave(file=paste(currentTitle, "_D.png", sep=""), plot=g3, dpi = 600, width = 15, height = 15, units = "cm")
ggsave(file=paste(currentTitle, "_D.eps", sep=""), plot=g3, device="eps", dpi = 600, width = 15, height = 15, units = "cm")
ggsave(file=paste(currentTitle, "_D.pdf", sep=""), plot=g3, device=cairo_pdf, dpi = 600, width = 15, height = 15, units = "cm")

# save plot g4 with the full supply of oxygen to the cell
ggsave(file=paste(currentTitle, "_supply_logvol", logV_ref, "_C0", C0, ".png", sep=""), plot=g4, dpi = 600, width = 15, height = 15, units = "cm")
ggsave(file=paste(currentTitle, "_supply_logvol", logV_ref, "_C0", C0, ".eps", sep=""), plot=g4, device="eps", dpi = 600, width = 15, height = 15, units = "cm")
ggsave(file=paste(currentTitle, "_supply_logvol", logV_ref, "_C0", C0, ".pdf", sep=""), plot=g4, device=cairo_pdf, dpi = 600, width = 15, height = 15, units = "cm")





# Now I repeat for a large range of body sizes
logV_ref <- 4.25
V_ref <- 10^logV_ref # a reference volume in the range of those observed
rs <- (3/8/pi*V_ref)^(1/3) # 25/2 # um
rl <- rs * 2 # 50/2 # um

r_eq <- (rs*rs*rl)^(1/3) # for the purpose of this script I
# approximate the cell to a sphere, because it is easier to
# model diffusion in this way

r_eq_test <- r_eq
# The volume of this cell would be V = 4/3*pi*r_eq^3

# I want to test the effects of changing the volume over the
# range below:
allLogV_test <- seq(3, 6, by=0.05)
allData <- data.frame(logV_test=double(), temp = double(), B = double(), C = double(), D = double(), Energy_supply_nw = double())

for (vvv in 1:length(allLogV_test))
{

r_eq_test <- (10^allLogV_test[vvv]*3/4/pi)^(1/3)
testedTemperature <- seq(10, 40, by=0.5)# c(10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40)


# scale_metabolic_rate_with_temperature(Bmeasured, adaptedTemperature + 273.15, tempH=31+273.15, temp=testedTemperature + 273.15, Ea, 7)
B <- scale_metabolic_rate_with_body_size_and_temperature(B0=Bmeasured, M0=4/3*pi*r_eq^3, M=4/3*pi*r_eq_test^3, tempRef=adaptedTemperature + 273.15, tempH=31+273.15, temp=testedTemperature + 273.15, E=Ea, ED=7, betaExponent=1)
# B is the metabolic rate in nW

C <- calculate_gas_concentration_from_partial_pressure(testedTemperature + 273.15, "O2", 0.20948)
# C is the concentration of O2 in saturated water, in mol/L

D <- scale_diffusion_coefficient_with_temperature(D0 = 2.10e-9, tempRef = 293.15, testedTemperature + 273.15)
# D is the diffusion coefficient of oxygen in water expressed in m^2/s

Oxygen_supply <- 4*pi*D*((C - C0)*10^9)*(r_eq*10^-6)*60 # this is the current (flux x cell surface area) of Oxygen into the cell
# This is the oxygen supply in umol/cell/min
Energy_supply_nW <- respiration_rate_to_W(Oxygen_supply)*10^9

df1 <- data.frame(logV_test = allLogV_test[vvv], temp = testedTemperature, B = B, C = C, D = D, Energy_supply_nw = Energy_supply_nW)

allData <- rbind(allData, df1)

}

allData$Energy_balance_nW <- allData$Energy_supply_nw/allData$B


pB <-ggplot(allData,aes(x=temp,y=logV_test,z=B,fill=B)) +
  # geom_tile() + 
  geom_contour(breaks=10^seq(-2,2,by=0.5), aes(colour = after_stat(level)), na.rm=TRUE, size=1.2) +
  scale_fill_viridis(name="Metabolic rate (nW)", option ="C", limits = c(0, 10)) +
  scale_color_viridis(name="Energy supply/demand", option="C", trans = "log10", limits = c(0.01, 100))

pB <-pB + scale_y_continuous(name=expression(paste('log'[10]*'(volume)'," ",  mu, 'm'^3)), limits=c(min(allLogV_test), max(allLogV_test)), expand = c(0, 0))
pB <-pB + scale_x_continuous(name="Temperature (°C)", breaks =seq(10,40,by=5), limits=c(10,35),expand = c(0, 0))
pB <-pB + theme_classic(base_size = 12)
pB <-pB + ggtitle("Metabolic rate (nW)")
pB <-pB + theme(legend.position = "none")

pB #the plot

# save plot pB with the metabolic rate of the cell
ggsave(file=paste(currentTitle, "_MR_C0", C0, ".png", sep=""), plot=pB, dpi = 600, width = 15, height = 15, units = "cm")
ggsave(file=paste(currentTitle, "_MR_C0", C0, ".eps", sep=""), plot=pB, device="eps", dpi = 600, width = 15, height = 15, units = "cm")
ggsave(file=paste(currentTitle, "_MR_C0", C0, ".pdf", sep=""), plot=pB, device=cairo_pdf, dpi = 600, width = 15, height = 15, units = "cm")



p <-ggplot(allData,aes(x=temp,y=logV_test,z=Energy_balance_nW, fill=Energy_balance_nW)) +
  # geom_tile() + 
  geom_contour(breaks=10^seq(0, 11, by=0.5), aes(colour = after_stat(level)), na.rm=TRUE, size=1.2) +
  geom_contour(breaks=1, colour="black", size=2, na.rm=TRUE) +
  scale_fill_viridis(name="Energy supply/demand", option ="C", limits = c(-200, 0)) +
  scale_color_viridis(name="Energy supply/demand", option="C", trans = "log10", limits = c(1, 10000))

p <-p + scale_y_continuous(name=expression(paste('log'[10]*'(volume)'," ",  mu, 'm'^3)), limits=c(min(allLogV_test), max(allLogV_test)), expand = c(0, 0))
p <-p + scale_x_continuous(name="Temperature (°C)", breaks =seq(10,40,by=5), limits=c(10,35),expand = c(0, 0))
p <-p + theme_classic(base_size = 12)
p <-p + ggtitle("Energy supply/demand")
p <-p + theme(legend.position = "none")


p #the plot

# save plot p with the full supply/demand ratio of oxygen to the cell
ggsave(file=paste(currentTitle, "_ratio_C0", C0, ".png", sep=""), plot=p, dpi = 600, width = 15, height = 15, units = "cm")
ggsave(file=paste(currentTitle, "_ratio_C0", C0, ".eps", sep=""), plot=p, device="eps", dpi = 600, width = 15, height = 15, units = "cm")
ggsave(file=paste(currentTitle, "_ratio_C0", C0, ".pdf", sep=""), plot=p, device=cairo_pdf, dpi = 600, width = 15, height = 15, units = "cm")










######## Now we try to include the data ###############
# The data that we use are those already processed by other scripts, and in particular
# we are interested in the rates of respiration and speed at the reference temperature
# Activation energies and rates at reference temperature can be calculated for each line, or for each 
# experimental conditions: there are 9 experimental conditions, 3 for temperature
# and 3 for medium concentration, but there are 4 replicate lines in each condition
# so the data can be either 9 values of activation energy, or 9x4=36. The files 
# related to respiration data are named as follows:
# allFitResults_metabolic_rate_post-adaptation.csv (36 values of activation energy)
# and 
# allFitResults_metabolic_rate_post-adaptation_all.csv (9 values)
# These files are produced by the script:
# fit_schoolfield_to_respiration_data_post-adaptation.R
#
# The files for movement speed are similar, and we have:
# allFitResults_Schoolfield_on_Speed_acute_response.csv (values of activation energy
# for each line)
# and
# allFitResults_Schoolfield_on_Stokes_Power_acute_response_all.csv (9 values, lines
# are combined together).
# These files are produced by the script:
# fit_schoolfield_equation_to_acute_speed_response.R
# 
# We need to remember that the thermal response curves on acute speed responses 
# were only measured on a subset of replicate lines (typically 2 out of 4)
# For these reasons I here only take the aggregate values (not knowing exactly
# how to pair observations between respiration and movement data
# Body sizes at reference temperature could be obtained from the
# long-term speed and population growth dataset, or from the acute-speed dataset
# in this case it seems more appropriate to use the short-term/acute-speed dataset
# because it is the same that is used for the speed measurements against which 
# we want to compare the effects of cell volume.
# The files for cell volume are:
# estimated_log_volume_adapted_lines.csv (values of cell volume for each line)
# and
# estimated_log_volume_adapted_lines_all.csv (only grouped by condition)
# both files are produced by the script
# analysis_of_body_size_from_tracking_data.R

# the file with the activation energies and rates at reference temperature for movement speed
fileNameSpeed <- "~/Tetrahymena/acute_speed_response/allFitResults_Schoolfield_on_Speed_acute_response_all.csv" # file.choose() # ask the user to select a file name. # file.choose() # ask the user to select a file name.

# the file with the activation energies and rates at reference temperature for metabolic rate measured from respiration
fileNameRespiration <- "~/Tetrahymena/Metabolic_rate_respiration/post-data/allFitResults_metabolic_rate_post-adaptation_all.csv" # file.choose() # ask the user to select a file name.

# the file with measurements of the log10(cell volume) in the 9 adaptation conditions
fileNameCellVolume <- "~/Tetrahymena/body_size/estimated_log_volume_adapted_lines_all.csv" 

dSpeed <- read.table(file = fileNameSpeed, sep = ",", header=TRUE, stringsAsFactors = FALSE, strip.white=TRUE, na.strings = c("NA", " NA"))
dResp <- read.table(file = fileNameRespiration, sep = ",", header=TRUE, stringsAsFactors=FALSE, strip.white=TRUE, na.strings = c("NA", " NA"))
dResp <- subset(dResp, tAdaptation != "control") # exclude the controls from the respiration data
dVol <- read.table(file = fileNameCellVolume, sep = ",", header=TRUE, stringsAsFactors = FALSE, strip.white=TRUE, na.strings = c("NA", " NA"))
names(dResp)
names(dSpeed)
names(dVol)
sapply(dResp, typeof) # sapply(dResp, class)
dResp$tAdaptation <- as.factor(dResp$tAdaptation)
dResp$mediumConcentration <- as.factor(dResp$mediumConcentration)

sapply(dSpeed, typeof) # sapply(dSpeed, class)
dSpeed$tAdaptation <- as.factor(dSpeed$tAdaptation)
dSpeed$mediumConcentration <- as.factor(dSpeed$mediumConcentration)

sapply(dVol, typeof) # sapply(dVol, class)
dVol$tAdaptation <- as.factor(dVol$tAdaptation)
dVol$mediumConcentration <- as.factor(dVol$mediumConcentration)




# merge the data frames
dList <- list(dResp, dSpeed, dVol)
library(tidyverse)
dMerged <- dList %>% reduce(full_join, by = c("tAdaptation", "mediumConcentration"))


originalCondition <- subset(dMerged, tAdaptation == 20 & mediumConcentration == 100)
conditionT25M100 <- subset(dMerged, tAdaptation == 25 & mediumConcentration == 100)
conditionT15M100 <- subset(dMerged, tAdaptation == 15 & mediumConcentration == 100)

conditionT20M200 <- subset(dMerged, tAdaptation == 20 & mediumConcentration == 200)
conditionT25M200 <- subset(dMerged, tAdaptation == 25 & mediumConcentration == 200)
conditionT15M200 <- subset(dMerged, tAdaptation == 15 & mediumConcentration == 200)

conditionT20M50 <- subset(dMerged, tAdaptation == 20 & mediumConcentration == 50)
conditionT25M50 <- subset(dMerged, tAdaptation == 25 & mediumConcentration == 50)
conditionT15M50 <- subset(dMerged, tAdaptation == 15 & mediumConcentration == 50)


activationEnergyT20M100 <- originalCondition$e_from_fit.x
rateAtTAdapt <- originalCondition$r_tref.x
# logCellVolumeT20M100 <- originalCondition$estimatedlogVolume

adaptedTemperature <- 20 # degrees Celsius
Bmeasured <- rateAtTAdapt # nW
Ea <- activationEnergyT20M100 # activation energy
C0 <- 10^-4 #(0.00001 is 10 umol/L) # Oxygen concentration
# according to this reference:
# Wilson, D.F., Rumsey, W.L., Green, T.J. and Vanderkooi, J.M., 1988. The oxygen dependence of mitochondrial oxidative phosphorylation measured by a new optical method for measuring oxygen concentration. Journal of Biological Chemistry, 263(6), pp.2712-2718.
# oxygen concentrations above 0.00001 M or 10 uM would be not limiting anymore at the level of mitochondria


logV_ref <- 4.25
V_ref <- 10^logV_ref # a reference volume in the range of those observed
rs <- (3/8/pi*V_ref)^(1/3) # 25/2 # um
rl <- rs * 2 # 50/2 # um

r_eq <- (rs*rs*rl)^(1/3) # for the purpose of this script I
# approximate the cell to a sphere, because it is easier to
# model diffusion in this way

r_eq_test <- r_eq
# The volume of this cell would be V = 4/3*pi*r_eq^3

# I want to test the effects of changing the volume over the
# range below:
allLogV_test <- seq(3, 6, by=0.05)
allData <- data.frame(logV_test=double(), temp = double(), B = double(), C = double(), D = double(), Energy_supply_nw = double())

for (vvv in 1:length(allLogV_test))
{
  print(vvv)
  r_eq_test <- (10^allLogV_test[vvv]*3/4/pi)^(1/3)
  testedTemperature <- seq(10, 40, by=0.5)# c(10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40)
  
  
  # scale_metabolic_rate_with_temperature(Bmeasured, adaptedTemperature + 273.15, tempH=31+273.15, temp=testedTemperature + 273.15, Ea, 7)
  B <- scale_metabolic_rate_with_body_size_and_temperature(B0=Bmeasured, M0=4/3*pi*r_eq^3, M=4/3*pi*r_eq_test^3, tempRef=adaptedTemperature + 273.15, tempH=31+273.15, temp=testedTemperature + 273.15, E=Ea, ED=7, betaExponent=1)
  # B is the metabolic rate in nW
  
  C <- calculate_gas_concentration_from_partial_pressure(testedTemperature + 273.15, "O2", 0.20948)
  # C is the concentration of O2 in saturated water, in mol/L
  
  D <- scale_diffusion_coefficient_with_temperature(D0 = 2.10e-9, tempRef = 293.15, testedTemperature + 273.15)
  # D is the diffusion coefficient of oxygen in water expressed in m^2/s
  
  Oxygen_supply <- 4*pi*D*((C - C0)*10^9)*(r_eq*10^-6)*60 # this is the current (flux x cell surface area) of Oxygen into the cell
  # This is the oxygen supply in umol/cell/min
  Energy_supply_nW <- respiration_rate_to_W(Oxygen_supply)*10^9
  
  df1 <- data.frame(logV_test = allLogV_test[vvv], temp = testedTemperature, B = B, C = C, D = D, Energy_supply_nw = Energy_supply_nW)
  
  allData <- rbind(allData, df1)
  
}

allData$Energy_balance_nW <- allData$Energy_supply_nw/allData$B


pB <-ggplot(allData,aes(x=temp,y=logV_test,z=B,fill=B)) +
  # geom_tile() + 
  geom_contour(breaks=10^seq(-2,2,by=0.5), aes(colour = after_stat(level)), na.rm=TRUE, size=1.2) +
  scale_fill_viridis(name="Metabolic rate (nW)", option ="C", limits = c(0, 10)) +
  scale_color_viridis(name="Energy supply/demand", option="C", trans = "log10", limits = c(0.01, 100))

pB <-pB + scale_y_continuous(name=expression(paste('log'[10]*'(volume)'," ",  mu, 'm'^3)), limits=c(min(allLogV_test), max(allLogV_test)), expand = c(0, 0))
pB <-pB + scale_x_continuous(name="Temperature (°C)", breaks =seq(10,40,by=5), limits=c(10,35),expand = c(0, 0))
pB <-pB + theme_classic(base_size = 12)
pB <-pB + ggtitle("Metabolic rate (nW)")
pB <-pB + theme(legend.position = "none")

pB #the plot

# save plot pB with the metabolic rate of the cell
ggsave(file=paste(currentTitle, "_MR_with_data_C0", C0, ".png", sep=""), plot=pB, dpi = 600, width = 15, height = 15, units = "cm")
ggsave(file=paste(currentTitle, "_MR_with_data_C0", C0, ".eps", sep=""), plot=pB, device="eps", dpi = 600, width = 15, height = 15, units = "cm")
ggsave(file=paste(currentTitle, "_MR_with_data_C0", C0, ".pdf", sep=""), plot=pB, device=cairo_pdf, dpi = 600, width = 15, height = 15, units = "cm")



p <-ggplot(data=allData,aes(x=temp,y=logV_test,z=Energy_balance_nW)) +
  # geom_tile() + 
  geom_contour(aes(colour = after_stat(level)), breaks=10^seq(0, 11, by=0.5), na.rm=TRUE, size=1.2) +
  geom_contour(breaks=1, colour="black", size=2, na.rm=TRUE) +
  scale_fill_viridis(name="Oxygen supply/demand", option ="C", limits = c(-200, 0)) +
  scale_color_viridis(name="Oxygen supply/demand", option="C", trans = "log10", limits = c(1, 10000))


p <-p + scale_y_continuous(name=expression(paste('log'[10]*'(volume)'," ",  mu, 'm'^3)), limits=c(min(allLogV_test), max(allLogV_test)), expand = c(0, 0))
p <-p + scale_x_continuous(name="Temperature (°C)", breaks =seq(10,40,by=5), limits=c(10,35),expand = c(0, 0))
p <-p + theme_classic(base_size = 12)
p <-p + ggtitle("Oxygen supply/demand")
p <-p + theme(legend.position = "none")
T15 <- conditionT15M100$tAdaptation
V15 <- conditionT15M100$estimatedlogVolume
T20 <- originalCondition$tAdaptation
V20 <- originalCondition$estimatedlogVolume
T25 <- conditionT25M100$tAdaptation
V25 <- conditionT25M100$estimatedlogVolume

T15low <- conditionT15M50$tAdaptation
V15low <- conditionT15M50$estimatedlogVolume
T20low <- conditionT20M50$tAdaptation
V20low <- conditionT20M50$estimatedlogVolume
T25low <- conditionT25M50$tAdaptation
V25low <- conditionT25M50$estimatedlogVolume

T15high <- conditionT15M200$tAdaptation
V15high <- conditionT15M200$estimatedlogVolume
T20high <- conditionT20M200$tAdaptation
V20high <- conditionT20M200$estimatedlogVolume
T25high <- conditionT25M200$tAdaptation
V25high <- conditionT25M200$estimatedlogVolume

p <-p + geom_segment(aes(x = as.numeric(levels(T20))[T20], y = V20, xend = as.numeric(levels(T25))[T25], yend = V25),
                     arrow = arrow(length = unit(0.3, "cm")), colour="darkred", size=1.2)
p <-p + geom_segment(aes(x = as.numeric(levels(T20))[T20], y = V20, xend = as.numeric(levels(T15))[T15], yend = V15),
                     arrow = arrow(length = unit(0.3, "cm")), colour="blue", size=1.2)

# p <-p + geom_segment(aes(x = as.numeric(levels(T20low))[T20low], y = V20low, xend = as.numeric(levels(T25low))[T25low], yend = V25low),
#                      arrow = arrow(length = unit(0.5, "cm")), colour="darkred")
# p <-p + geom_segment(aes(x = as.numeric(levels(T20low))[T20low], y = V20low, xend = as.numeric(levels(T15low))[T15low], yend = V15low),
#                      arrow = arrow(length = unit(0.5, "cm")), colour="blue")
# 
# p <-p + geom_segment(aes(x = as.numeric(levels(T20high))[T20high], y = V20high, xend = as.numeric(levels(T25high))[T25high], yend = V25high),
#                      arrow = arrow(length = unit(0.5, "cm")), colour="darkred")
# p <-p + geom_segment(aes(x = as.numeric(levels(T20high))[T20high], y = V20high, xend = as.numeric(levels(T15high))[T15high], yend = V15high),
#                      arrow = arrow(length = unit(0.5, "cm")), colour="blue")

p <-p + annotate(
  "text", label = expression(paste(Phi, '=1')),
  x = 19, y = 5.9, size = 5, colour = "black"
)

p <-p + annotate(colour = plasma(8)[2],
  "text", label = expression(paste(Phi, '=10')),
  x = 10.5, y = 5.65, size = 5, hjust=0)

p <-p + annotate(colour = plasma(8)[4],
                 "text", label = expression(paste(Phi, '=100')),
                 x = 10.5, y = 4.65, size = 5, hjust=0)

p <-p + annotate(colour = plasma(8)[6],
                 "text", label = expression(paste(Phi, '=1000')),
                 x = 10.5, y = 3.65, size = 5, hjust=0)

p #the plot

# save plot p with the full supply/demand ratio of oxygen to the cell
ggsave(file=paste(currentTitle, "_ratio_with_data_C0", C0, ".png", sep=""), plot=p, dpi = 600, width = 15, height = 15, units = "cm")
ggsave(file=paste(currentTitle, "_ratio_with_data_C0", C0, ".eps", sep=""), plot=p, device="eps", dpi = 600, width = 15, height = 15, units = "cm")
ggsave(file=paste(currentTitle, "_ratio_with_data_C0", C0, ".pdf", sep=""), plot=p, device=cairo_pdf, dpi = 600, width = 15, height = 15, units = "cm")


