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
  mu <- viscA * exp(viscB/temp + viscC*temp + viscD*temp^2) # dynamic viscosity of water, approx. 1e-3 Kg/(m s) Andreade equation with extra parameters
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
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)
setwd("~/Tetrahymena/theoretical_model")


if(!is.null(dev.list())) dev.off()


######################### Predicted speed for short-term exposure #########
# Under short-term exposure, speed increases as the square root of metabolic
# rate, as a consequence of drag being proportional to U^2 in viscous regime

adaptedTemperature <- 20 # degrees Celsius
d <- 1 # food density used 100%
Umeasured <- 500 # um / s
Bmeasured <- 0.7 # nW
rs <- 25/2 # um
rl <- 50/2 # um

Ea <- 0.75 # activation energy

d0 <- 1 # standard food density of 100% medium

# 6 and 250 look ok
fMax0 <- 6 # # nW, maximum feeding rate at standard conditions 
eta0 <- 250 #  # speed of half saturation feeding (um/s)

# If we look at the equation
# B * eta_c =  CD * mu * U^2
# B = sU^2
# Where B is the metabolic rate, eta_c is the efficiency of conversion of 
# metabolic energy to movement (in reality it is a bit more complicate
# as B is the metabolic rate measured from respiration so it has a basal component
# and a movement component (B = Bb + Bm) and the movement component has some
# efficiency of conversion given to dissipation of energy 
# at the level of ciliary movement etc. If the allocation of energy to movement
# and to other metabolic activities does not change, Bb is proportional to Bm
# and also to B, so in this case I can calculate this simplified eta_c that accounts for both 
# the fraction of B allocated to Bm and the efficiency of conversion of Bm to movement)
# CD is the drag coefficient
# mu is the water viscosity
# we can retrieve eta_c
CDmuoveretac <- Bmeasured * 1e-9 / (Umeasured * 1e-6)^2 # s
E <- rl/rs # elongation of the spheroid
CD <- 1.2 * pi * (rs * 10^-6) * (4+E)
mu <- calculate_water_viscosity(adaptedTemperature + 273.15)
eta_c <- CD * mu / CDmuoveretac

testedTemperature <- seq(10, 40, by=1)# c(10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40)
foodDensity <- 1

B <- scale_metabolic_rate_with_temperature(Bmeasured, adaptedTemperature + 273.15, tempH=31+273.15, temp=testedTemperature + 273.15, Ea, 7)
allMu <- calculate_water_viscosity(testedTemperature + 273.15)
predictedSpeed <- sqrt(B * 1e-9 / CDmuoveretac * mu / allMu) * 1e6

df1 <- data.frame(temp = testedTemperature, B = B, rate = predictedSpeed)

labels = c("Speed")
colours = c("#000000")

g0 <- ggplot() +
  geom_line(data = df1, size=2, linetype = 1, aes(x=temp, y = B, color=labels[1])) +
  annotate("text", size=5, x=10.5, y=1.75, label= paste("E=", Ea, "eV", sep=""), hjust = 0, parse=F) +
  # annotate("text", size=5, x=10.5, y=3.15, label= paste("E=", Ea, "eV", sep=""), hjust = 0, parse=F) +
  scale_x_continuous(name=expression(paste("Temperature °C "))) +
  scale_y_continuous(name="metabolic rate (nW)", limits=c(0,2)) +
  # scale_y_continuous(name="metabolic rate (nW)", limits=c(0,3.6)) +
  theme_classic(base_size = 22) + 
  scale_color_manual(name="Contribution", breaks=labels, values=colours) +
  theme(legend.position = "none") # + # This removes the legend above
# ggtitle(paste("Predicted speed, short term exposure"))
g0

ggsave(file="aaaa_estimated_short_term_MR_from_model.pdf", device=cairo_pdf, dpi = 1200, width = 12, height = 10, units = "cm")
ggsave(file="aaaa_estimated_short_term_MR_from_model.png", dpi = 600, width = 12, height = 10, units = "cm")


g1 <- ggplot() + 
  geom_line(data = df1, size=1, linetype = 1, aes(x=temp, y = rate, color=labels[1])) +
  scale_x_continuous(name=expression(paste("Temperature °C "))) +
  scale_y_continuous(name=expression(paste("Speed (", mu, "m/s", ")")),limits=c(0,1100), breaks=seq(0, 1000, by=200)) +
  theme_classic(base_size = 22) + 
  scale_color_manual(name="Contribution", breaks=labels, values=colours) +
  theme(legend.position = "none") # + # This removes the legend above
# ggtitle(paste("Predicted speed, short term exposure"))



# choose model
mod = 'sharpschoolhigh_1981'

# get start vals
start_vals <- get_start_vals(df1$temp, df1$rate, model_name = 'sharpeschoolhigh_1981')

# get limits
low_lims <- get_lower_lims(df1$temp, df1$rate, model_name = 'sharpeschoolhigh_1981')
upper_lims <- get_upper_lims(df1$temp, df1$rate, model_name = 'sharpeschoolhigh_1981')

start_vals
low_lims
upper_lims



# fit model
fit <- nls_multstart(rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 20),
                     data = df1,
                     iter = 500,
                     start_lower = start_vals - 10,
                     start_upper = start_vals + 10,
                     lower = low_lims,
                     upper = upper_lims,
                     supp_errors = 'Y')



# calculate additional traits
calculatedFitParameters <- calc_params(fit) %>%
  # round for easy viewing
  mutate_all(round, 2)

print(fit)

new_data <- data.frame(temp = testedTemperature)
preds <- augment(fit, newdata = new_data)

g1 <- g1 + 
  geom_line(aes(temp, .fitted), preds, col = 'black', size=2) +
  annotate("text", size=5, x=10.5, y=960, label= paste("E=", calculatedFitParameters$e, "eV", sep=""), hjust = 0, parse=F)
g1
ggsave(file="aaaa_estimated_short_term_speed_from_model.pdf", device=cairo_pdf, dpi = 1200, width = 12, height = 10, units = "cm")
ggsave(file="aaaa_estimated_short_term_speed_from_model.png", dpi = 600, width = 12, height = 10, units = "cm")





########### Now use the same strategy to look for optimal speed:
# feeding rate at test temperature, body size and food density
fMax <- fMax0 # * B/Bmeasured # faster metabolic rate increases the maximum feeding rate (shorter handling time)
eta <- eta0 * d0/d # * B/Bmeasured  # CS0/CS # not clear if eta should scale with the cross section
s <- CDmuoveretac * allMu / mu
Uopt <- calculate_optimal_speed(s, eta, fMax)
df2 <- data.frame(temp = testedTemperature[1:21], rate = Uopt[1:21] + rnorm(21, mean=0, sd=80))

# I add an artificial data point for temperature 40 degrees and rate=0
dExtra <- df2[1,]
dExtra$temp <- 40
dExtra$rate <- 0
df2 <- rbind(df2, dExtra)

g2 <- ggplot() + 
  geom_line(data = df2, linetype = 1, size= 1, aes(x=temp, y = rate, color=labels[1])) +
  scale_x_continuous(name=expression(paste("Temperature °C "))) +
  scale_y_continuous(name=expression(paste("Speed (", mu, "m/s", ")")), lim=c(0,800)) +
  theme_classic(base_size = 18) + 
  scale_color_manual(name="Contribution", breaks=labels, values=colours) +
  theme(legend.position = "none") # + # This removes the legend above
# ggtitle(paste("Predicted speed, long term exposure"))

g2


# choose model
mod = 'sharpschoolhigh_1981'

# get start vals
start_vals <- get_start_vals(df2$temp, df2$rate, model_name = 'sharpeschoolhigh_1981')

# get limits
low_lims <- get_lower_lims(df2$temp, df2$rate, model_name = 'sharpeschoolhigh_1981')
upper_lims <- get_upper_lims(df2$temp, df2$rate, model_name = 'sharpeschoolhigh_1981')

start_vals
low_lims
upper_lims



# fit model
fit <- nls_multstart(rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 20),
                     data = df2,
                     iter = 500,
                     start_lower = start_vals - 10,
                     start_upper = start_vals + 10,
                     lower = low_lims,
                     upper = upper_lims,
                     supp_errors = 'Y')



# calculate additional traits
calculatedFitParameters <- calc_params(fit) %>%
  # round for easy viewing
  mutate_all(round, 2)

print(fit)

new_data <- data.frame(temp = testedTemperature)
preds <- augment(fit, newdata = new_data)

g2 <- g2 + 
  geom_line(aes(temp, .fitted), preds, col = 'red', size=2)

g2




variousMu <- calculate_water_viscosity(c(15, 20, 25) + 273.15)
referenceMu <- calculate_water_viscosity(adaptedTemperature + 273.15)

CDmuoveretac * variousMu / referenceMu

allU <- seq(0, 1200, by=1) * 10^-6 # speed in m/s

etaSmall <- eta*10^-6 # convert eta to m/s from um/s
fSmall <- fMax * 10^-9 # convert the feeding rate from nW to W

etaSmallLow <- etaSmall * d0/0.5 # * B/Bmeasured  # CS0/CS # not clear if eta should scale with the cross section
etaSmallMedium <- etaSmall * d0/1 # * B/Bmeasured  # CS0/CS # not clear if eta should scale with the cross section
etaSmallHigh <- etaSmall * d0/2 # * B/Bmeasured  # CS0/CS # not clear if eta should scale with the cross section


der1 <- 1e9 * etaSmall * fSmall / (allU + etaSmall)^2
der1Low <- 1e9 * etaSmallLow * fSmall / (allU + etaSmallLow)^2 # also consider other medium concentrations
der1Medium <- 1e9 * etaSmallMedium * fSmall / (allU + etaSmallMedium)^2 # also consider other medium concentrations
der1High <- 1e9 * etaSmallHigh * fSmall / (allU + etaSmallHigh)^2 # also consider other medium concentrations
der2T15 <- 1e9 * 2 * CDmuoveretac * variousMu[1] / referenceMu * allU
der2T20 <- 1e9 * 2 * CDmuoveretac * variousMu[2] / referenceMu * allU
der2T25 <- 1e9 * 2 * CDmuoveretac * variousMu[3] / referenceMu * allU

derData <- data.frame(allU = allU*10^6, der1 = der1, der2T15=der2T15, der2T20=der2T20, der2T25=der2T25, der1Low = der1Low, der1Medium = der1Medium, der1High = der1High)

# one equation is 2sU
# the other is eta fMax / (U + eta)^2
g3 <- ggplot() + 
  geom_line(data = derData, col="black", linetype = 1, aes(x=allU, y = der1)) +
  geom_line(data = derData, col="#3B9AB2", linetype = 1, aes(x=allU, y = der2T15)) +
  geom_line(data = derData, col="#EBCC2A", linetype = 1, aes(x=allU, y = der2T20)) +
  geom_line(data = derData, col="#F21A00", linetype = 1, aes(x=allU, y = der2T25)) +
  scale_y_continuous(name=expression(paste("terms (nW s / m)")), lim=c(0,5000)) +
  scale_x_continuous(name=expression(paste("Speed U (", mu, "m/s", ")")), lim=c(250, 750)) +
  theme_classic(base_size = 18) + 
  theme(legend.position = "none") # + # This removes the legend above
# ggtitle(paste("Predicted speed, long term exposure"))

g3
ggsave(file="aaaa_estimate_optimal_speed_from_derivatives.pdf", device=cairo_pdf, dpi = 1200, width = 12, height = 10, units = "cm")
ggsave(file="aaaa_estimate_optimal_speed_from_derivatives.png", dpi = 600, width = 12, height = 10, units = "cm")


# one equation is 2sU
# the other is eta fMax / (U + eta)^2
g3Bis <- ggplot() + 
  geom_line(data = derData, col="#A3A3A3", linetype = 1, aes(x=allU, y = der1Low)) +
  geom_line(data = derData, col="#666666", linetype = 1, aes(x=allU, y = der1Medium)) +
  geom_line(data = derData, col="#000000", linetype = 1, aes(x=allU, y = der1High)) +
  geom_line(data = derData, col="#3B9AB2", linetype = 1, aes(x=allU, y = der2T15)) +
  geom_line(data = derData, col="#EBCC2A", linetype = 1, aes(x=allU, y = der2T20)) +
  geom_line(data = derData, col="#F21A00", linetype = 1, aes(x=allU, y = der2T25)) +
  scale_y_continuous(name=expression(paste("terms (nW s / m)")), lim=c(0,5000)) +
  scale_x_continuous(name=expression(paste("Speed U (", mu, "m/s", ")")), lim=c(250, 750)) +
  theme_classic(base_size = 18) + 
  theme(legend.position = "none") # + # This removes the legend above
# ggtitle(paste("Predicted speed, long term exposure"))

g3Bis


ggsave(file="aaaa_estimate_optimal_speed_from_derivatives_all_concentrations.pdf", device=cairo_pdf, dpi = 1200, width = 12, height = 10, units = "cm")
ggsave(file="aaaa_estimate_optimal_speed_from_derivatives_all_concentrations.png", dpi = 600, width = 12, height = 10, units = "cm")





########### Now use the same strategy to look for optimal speed:
# feeding rate at test temperature, body size and food density
fMax <- fMax0 # * B/Bmeasured # faster metabolic rate increases the maximum feeding rate (shorter handling time)

etaLow <- eta0 * d0/0.5 # * B/Bmeasured  # CS0/CS # not clear if eta should scale with the cross section
etaMedium <- eta0 * d0/1 # * B/Bmeasured  # CS0/CS # not clear if eta should scale with the cross section
etaHigh <- eta0 * d0/2 # * B/Bmeasured  # CS0/CS # not clear if eta should scale with the cross section

s <- CDmuoveretac * allMu / mu
# equivalent to: 
# s <- allMu * CD / eta_c
# s <- calculate_stokes_coefficient(eta_c, rs, rl, allMu, "simplifiedSpheroid")
UoptLow <- calculate_optimal_speed(s, etaLow, fMax)
UoptMedium <- calculate_optimal_speed(s, etaMedium, fMax)
UoptHigh <- calculate_optimal_speed(s, etaHigh, fMax)

df3 <- data.frame(temp = testedTemperature[1:21], UoptLow = UoptLow[1:21], UoptMedium = UoptMedium[1:21], UoptHigh = UoptHigh[1:21])

# pick the values that correspond to the experimental conditions
library(reshape2)
selectedValues3 <- melt(subset(df3, temp %in% c(15, 20 ,25)),id.vars=c("temp"))
selectedValues3$mediumConcentration <- 0
selectedValues3$mediumConcentration[selectedValues3$variable == "UoptLow"] <- 50
selectedValues3$mediumConcentration[selectedValues3$variable == "UoptMedium"] <- 100
selectedValues3$mediumConcentration[selectedValues3$variable == "UoptHigh"] <- 200
names(selectedValues3) <- c("temp", "variable", "speed", "mediumConcentration")

lineColours = c("#A3A3A3", "#666666", "#000000", "#FF00FF")

g3 <- ggplot() + 
  geom_line(data = df3, size=0.8, linetype = 1, color=lineColours[1], aes(x=temp, y = UoptLow)) +
  geom_line(data = df3, size=1.4, linetype = 1, color=lineColours[2], aes(x=temp, y = UoptMedium)) +
  geom_line(data = df3, size=2, linetype = 1, color=lineColours[3], aes(x=temp, y = UoptHigh)) +
  geom_point(size=6, data=selectedValues3, aes(x=temp, y=speed, fill=factor(temp * 100 + mediumConcentration),  shape=factor(mediumConcentration))) +
  scale_x_continuous(name=expression(paste("Temperature °C "))) +
  scale_y_continuous(name=expression(paste("Speed (", mu, "m/s", ")")), lim=c(0,600)) +
  theme_classic(base_size = 18) + 
  scale_color_manual(name="Contribution", breaks=labels, values=colours) +
  scale_shape_manual(values=c(21, 24, 22)) + # shapes for the markers
  scale_fill_manual(values=c("#8FCDDC", "#3B9AB2", "#18434E", "#F2DD70", "#EBCC2A", "#463C07", "#FF7C6C", "#F21A00", "#660B00")) +
  theme(legend.position = "none") # + # This removes the legend above
# ggtitle(paste("Predicted speed, long term exposure"))

g3


ggsave(file="aaaa_predicted_long_term_speed_vs_temperature.pdf", device=cairo_pdf, dpi = 1200, width = 12, height = 10, units = "cm")
ggsave(file="aaaa_predicted_long_term_speed_vs_temperature.png", dpi = 600, width = 12, height = 10, units = "cm")


# Now same figure as a function of medium density and for the three temperature
# conditions


allD <- seq(0, 2, by=0.01) # all the range of densities
allEta <- eta0 * d0 / allD
fMax <- fMax0

# change viscosity for each temperature condition
mu15 <- calculate_water_viscosity(zeroCelsius + 15)
mu20 <- calculate_water_viscosity(zeroCelsius + 20)
mu25 <- calculate_water_viscosity(zeroCelsius + 25)

# calculate stokes coefficient (in this simplified model I do not adjust
# cell size or shape according to temperature, and I also keep eta_c the same
# or it would be difficult to understand all the different changes being produced)
s15 <- calculate_stokes_coefficient(eta_c, rs, rl, mu15, "simplifiedSpheroid")
s20 <- calculate_stokes_coefficient(eta_c, rs, rl, mu20, "simplifiedSpheroid")
s25 <- calculate_stokes_coefficient(eta_c, rs, rl, mu25, "simplifiedSpheroid")


Uopt15 <- calculate_optimal_speed(s15, allEta, fMax)
Uopt20 <- calculate_optimal_speed(s20, allEta, fMax)
Uopt25 <- calculate_optimal_speed(s25, allEta, fMax)

df4 <- data.frame(density = allD, Uopt15 = Uopt15, Uopt20 = Uopt20, Uopt25 = Uopt25)


# pick the values that correspond to the experimental conditions
library(reshape2)
selectedValues4 <- melt(subset(df4, density %in% c(0.5, 1, 2)),id.vars=c("density"))
selectedValues4$temp <- 0
selectedValues4$temp[selectedValues4$variable == "Uopt15"] <- 15
selectedValues4$temp[selectedValues4$variable == "Uopt20"] <- 20
selectedValues4$temp[selectedValues4$variable == "Uopt25"] <- 25
names(selectedValues4) <- c("density", "variable", "speed", "temp")



lineColours = c("#3B9AB2", "#EBCC2A", "#F21A00", "#FF00FF")

g4 <- ggplot() + 
  geom_line(data = df4, size=1.4, linetype = 1, color=lineColours[1], aes(x=density*100, y = Uopt15)) +
  geom_line(data = df4, size=1.4, linetype = 1, color=lineColours[2], aes(x=density*100, y = Uopt20)) +
  geom_line(data = df4, size=1.4, linetype = 1, color=lineColours[3], aes(x=density*100, y = Uopt25)) +
  geom_point(size=6, data=selectedValues4, aes(x=density*100, y=speed, fill=factor(temp * 100 + density),  shape=factor(density))) +
  scale_x_continuous(name=expression(paste("Medium concentration (%)"))) +
  scale_y_continuous(name=expression(paste("Speed (", mu, "m/s", ")")), lim=c(0,600)) +
  theme_classic(base_size = 18) + 
  scale_color_manual(name="Contribution", breaks=labels, values=colours) +
  scale_shape_manual(values=c(21, 24, 22)) + # shapes for the markers
  scale_fill_manual(values=c("#8FCDDC", "#3B9AB2", "#18434E", "#F2DD70", "#EBCC2A", "#463C07", "#FF7C6C", "#F21A00", "#660B00")) +
  theme(legend.position = "none") # + # This removes the legend above
# ggtitle(paste("Predicted speed, long term exposure"))

g4
ggsave(file="aaaa_predicted_long_term_speed_vs_concentration.pdf", device=cairo_pdf, dpi = 1200, width = 12, height = 10, units = "cm")
ggsave(file="aaaa_predicted_long_term_speed_vs_concentration.png", dpi = 600, width = 12, height = 10, units = "cm")


library(ggpubr)
# combine together the two plots into a single figure
figure <- ggarrange(g3, g4,
                    # labels = c("A", "B"),
                    ncol = 2, nrow = 1)
figure

ggsave(file="figure_predicted_optimal_speed.eps", device="eps", dpi = 1200, width = 22, height = 8, units = "cm")
ggsave(file="figure_predicted_optimal_speed.png", dpi = 600, width = 22, height = 8, units = "cm")
library(Cairo)
ggsave(file="figure_predicted_optimal_speed.pdf", device=cairo_pdf, dpi = 1200, width = 22, height = 8, units = "cm")


