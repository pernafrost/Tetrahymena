rm(list=ls()) # clean memory
if(!is.null(dev.list())) dev.off()
saveFigures <- TRUE

library(ggplot2)
library(dplyr)

# NOTE that activation energies can be calculated for each line, or for each 
# experimental conditions: there are 9 experimental condition, 3 for temperature
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


# the file with the activation energies for movement speed
fileNameSpeed <- "allFitResults_Schoolfield_on_Speed_acute_response_all.csv" # "/../../acute_speed_response/allFitResults_Schoolfield_on_Speed_acute_response_all.csv" # file.choose() # ask the user to select a file name. # file.choose() # ask the user to select a file name.

# the file with the activation energies for metabolic rate measured from respiration
fileNameRespiration <- "allFitResults_metabolic_rate_post-adaptation_all.csv" # "../Metabolic_rate_respiration/post-data/allFitResults_metabolic_rate_post-adaptation_all.csv" # file.choose() # ask the user to select a file name.
setwd(dirname(fileNameRespiration))

dSpeed <- read.table(file = fileNameSpeed, sep = ",", header=TRUE, stringsAsFactors = FALSE, strip.white=TRUE, na.strings = c("NA", " NA"))
dResp <- read.table(file = fileNameRespiration, sep = ",", header=TRUE, stringsAsFactors=FALSE, strip.white=TRUE, na.strings = c("NA", " NA"))
dResp <- subset(dResp, tAdaptation != "control") # exclude the controls from the respiration data
names(dResp)
names(dSpeed)

# merge the two data frames
dRespAndSpeed <- merge(dResp, dSpeed, by = c("tAdaptation", "mediumConcentration"))

ggplot(dRespAndSpeed, aes(x=e_from_fit.x, y=e_from_fit.y, color=tAdaptation)) +
  geom_point() +
  geom_errorbarh(aes(xmin=e_from_fit.x - sd_e_from_fit.x, xmax = e_from_fit.x + sd_e_from_fit.x)) +
  geom_errorbar(aes(ymin=e_from_fit.y - sd_e_from_fit.y, ymax = e_from_fit.y + sd_e_from_fit.y)) +
  geom_abline(intercept = 0, slope = 0.5, color="black") +
  # geom_jitter(position = position_jitter(height = 0, width = .5)) +
  # geom_smooth(method="lm", se=TRUE, fill="red", formula=y ~ x, na.rm=TRUE) +
  scale_y_continuous(name="Activation Energy (speed)") +
  scale_x_continuous(name="Activation Energy (respiration)") +
  coord_cartesian(xlim=c(0,1), ylim=c(0,1), expand=FALSE) +
  theme_classic(base_size = 15) +
  theme(legend.position = "none") +
  ggtitle("Comparison of activation energies") +
  scale_color_manual(values=c("#3B9AB2", "#EBCC2A", "#F21A00")) # this is the zissou1 palette



if (saveFigures){
    ggsave(file="activation_energy_speed_respiration.png", dpi = 600, width = 12, height = 10, units = "cm")
}
  
  
  # merge the two data frames. This time I invert respiration and speed, given that
  # the error is mainly on respiration, and in this way I can run the ordinary
  # least squares fit.
  
  dSpeedAndResp <- merge(dSpeed, dResp, by = c("tAdaptation", "mediumConcentration"))
  
  ggplot(dSpeedAndResp, aes(x=e_from_fit.x, y=e_from_fit.y)) +
    geom_point() +
    geom_errorbarh(aes(xmin=e_from_fit.x - sd_e_from_fit.x, xmax = e_from_fit.x + sd_e_from_fit.x)) +
    geom_errorbar(aes(ymin=e_from_fit.y - sd_e_from_fit.y, ymax = e_from_fit.y + sd_e_from_fit.y)) +
    geom_smooth(method="lm", se=FALSE, fill="red", formula=y ~ x + 0, na.rm=TRUE) +
    scale_x_continuous(name="Activation Energy (speed)") +
    scale_y_continuous(name="Activation Energy (respiration)") +
    coord_cartesian(xlim=c(0,1), ylim=c(0,1), expand=FALSE) +
    theme_classic(base_size = 15) +
    theme(legend.position = "none") +
    ggtitle("Comparison of activation energies") +
    scale_color_manual(values=c("#3B9AB2", "#EBCC2A", "#F21A00")) # this is the zissou1 palette
  
  # I can fit a straight line
  fitSpeedVsResp <- lm(formula=e_from_fit.y ~ e_from_fit.x + 0, data = dSpeedAndResp)
  print(fitSpeedVsResp)
  
  # or I can simply compare the ratios between the two measures of activation energy
  meanRatioESpeedOverEResp <- mean(dSpeedAndResp$e_from_fit.x/dSpeedAndResp$e_from_fit.y)
  print(meanRatioESpeedOverEResp)
  
  # Here I create a variable with the difference of activation energies:
  dSpeedAndResp$differenceActivationEnergy <- dSpeedAndResp$e_from_fit.x * 2 - dSpeedAndResp$e_from_fit.y
  
  # I check that the data are normally distributed (check with Shapiro test)
  # given that we focus on the difference between the two conditions,
  # we only need to test that the difference is normally distributed;
  # we do not need to test the normality of the individual distributions
  
  shapiro.test(dSpeedAndResp$differenceActivationEnergy)
  # Shapiro-Wilk normality test
  # 
  # data:  dSpeedAndResp$differenceActivationEnergy
  # W = 0.92434, p-value = 0.4294
  # Good! (p is greater than 0.05)
  
  # As everything is good, I can compute the one sample t-test on the difference
  t.test(dSpeedAndResp$differenceActivationEnergy, mu=0)
  
  # In R I can also obtain the same result by calling the t.test with the option "paired"
  t.test(dSpeedAndResp$e_from_fit.x * 2, dSpeedAndResp$e_from_fit.y, paired=TRUE)
  # And similarly with the other value of activation energy returned by the rTPC library
  t.test(dSpeedAndResp$e.x * 2, dSpeedAndResp$e.y, paired=TRUE)

  
  
