# Changes in speed and respiration rate in the cells that had been adapted 
# to different temperature conditions and to different medium densities are 
# partly explained by changes in cell volume
# In this script we consider how much the change in volume alone is sufficient
# to account for the other observed changes in respiration rate and movement speed

rm(list=ls()) # clean memory
if(!is.null(dev.list())) dev.off()
saveFigures <- TRUE
if (saveFigures==TRUE) library(Cairo)
library(ggpubr) # for combining figures

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

setwd("~/Tetrahymena/relationship_between_volume_respiration_and_speed/")

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

tentativeSlopeSpeedvsRespiration <- mean(dMerged$e_from_fit.y/dMerged$e_from_fit.x)


ggplot(dMerged, aes(x=e_from_fit.x, y=e_from_fit.y, fill=tAdaptation, colour=tAdaptation, shape=mediumConcentration, size=mediumConcentration)) +
  geom_errorbarh(aes(xmin=e_from_fit.x - sd_e_from_fit.x, xmax = e_from_fit.x + sd_e_from_fit.x)) +
  geom_errorbar(aes(ymin=e_from_fit.y - sd_e_from_fit.y, ymax = e_from_fit.y + sd_e_from_fit.y)) +
  geom_point(size=3, colour="black") +
  geom_abline(intercept = 0, slope = tentativeSlopeSpeedvsRespiration, color="black") +
  annotate("text", size=5, x=0.1, y=0.8, label= paste("slope: ", round(tentativeSlopeSpeedvsRespiration, 2), sep=""), hjust = 0, parse=F) +
  scale_y_continuous(name="Activation Energy (speed)") +
  scale_x_continuous(name="Activation Energy (respiration)") +
  coord_cartesian(xlim=c(0,1), ylim=c(0,1), expand=FALSE) +
  theme_classic(base_size = 15) +
  theme(legend.position = "none") +
  scale_colour_manual(values=c("#3B9AB2", "#EBCC2A", "#F21A00", "#FF00FF")) +
  scale_fill_manual(values=c("#3B9AB2", "#EBCC2A", "#F21A00", "#FF00FF")) +
  scale_size_manual( values = c(0.4, 0.6, 0.8) ) +
  # ggtitle("Comparison of activation energies") +
  scale_shape_manual(values=c(21, 24, 22))  # shapes for the markers

if (saveFigures){
  ggsave(file="activation_energy_for_speed_and_respiration.png", dpi = 600, width = 12, height = 10, units = "cm")
  ggsave(file="activation_energy_for_speed_and_respiration.pdf", device=cairo_pdf, dpi = 1200, width = 24, height = 20, units = "cm")
}


########### Look at body size and metabolic rate
dMerged$logMR <- log10(dMerged$r_tref.x)
dMerged$logMR_0025 <- log10(dMerged$r_tref_0025.x)
dMerged$logMR_0975 <- log10(dMerged$r_tref_0975.x)
tentativeSlopeMRvsVolume <- lm(formula=logMR ~ estimatedlogVolume, data = dMerged)
print(tentativeSlopeMRvsVolume)
# if logMR = estimatedLogVolume x 1.234 -5.371
# then the metabolic rate would increase more than the volume alone

imposedSlope <- 1
tentativeInterceptMRvsVolumeWithImposedSlope <- mean(dMerged$logMR - dMerged$estimatedlogVolume * imposedSlope) 

gMRVol <- ggplot(dMerged, aes(x=estimatedlogVolume, y=logMR, fill=tAdaptation, colour=tAdaptation, shape=mediumConcentration, size=mediumConcentration)) +
  geom_errorbarh(aes(xmin=estimatedlogVolume_0025, xmax = estimatedlogVolume_0975)) +
  geom_errorbar(aes(ymin=logMR_0025, ymax = logMR_0975)) +
  geom_point(size=3, colour="black") +
  geom_abline(intercept = tentativeInterceptMRvsVolumeWithImposedSlope, slope = imposedSlope, color="black") +
  geom_abline(intercept = tentativeSlopeMRvsVolume$coefficients[1], slope = tentativeSlopeMRvsVolume$coefficients[2], color="darkgrey") +
  annotate("text", size=5, x=3.80, y=0.3, label= paste("best fitting slope: ", round(tentativeSlopeMRvsVolume$coefficients[2], 2), sep=""), hjust = 0, parse=F, colour="darkgrey") +
  annotate("text", size=5, x=3.80, y=0.2, label= paste("theoretical slope: ", round(imposedSlope, 2), sep=""), hjust = 0, parse=F, colour="black") +
  scale_x_continuous(name=expression(paste('log'[10]*'(volume)'," ",  mu, 'm'^3))) +
  scale_y_continuous(name=expression(paste('log'[10]*'(M.R.)'," nW"))) +
  coord_fixed(ratio=1,xlim=c(3.75,4.75), ylim=c(-0.6,0.4), expand=FALSE) +
  theme_classic(base_size = 15) +
  theme(legend.position = "none") +
  scale_colour_manual(values=c("#3B9AB2", "#EBCC2A", "#F21A00", "#FF00FF")) +
  scale_fill_manual(values=c("#3B9AB2", "#EBCC2A", "#F21A00", "#FF00FF")) +
  scale_size_manual( values = c(0.4, 0.6, 0.8) ) +
  scale_shape_manual(values=c(21, 24, 22))  # shapes for the markers
# the superposed line could be compatible with a linear scaling of 
# metabolic rate with cell volume, except in the low-nutrients conditions
# where perhaps the metabolic rate is limited by nutrients
gMRVol

if (saveFigures){
  ggsave(file="log10MR_vs_log10Vol.png", dpi = 600, width = 12, height = 10, units = "cm")
  ggsave(file="log10MR_vs_log10Vol.pdf", device=cairo_pdf, dpi = 1200, width = 24, height = 20, units = "cm")
}





########### Look at body size and movement speed
dMerged$logSpeed <- log10(dMerged$r_tref.y)
dMerged$logSpeed_0025 <- log10(dMerged$r_tref_0025.y)
dMerged$logSpeed_0975 <- log10(dMerged$r_tref_0975.y)
tentativeSlopeSpeedvsVolume <- lm(formula=logSpeed ~ estimatedlogVolume, data = dMerged)
print(tentativeSlopeSpeedvsVolume)
# if logSpeed = estimatedLogVolume x 0.5063 + 0.6121
# then the speed would increase more than the volume tothe power 1/3

imposedSlope <- 1/3
tentativeInterceptSpeedvsVolumeWithImposedSlope <- mean(dMerged$logSpeed - dMerged$estimatedlogVolume * imposedSlope) 


gSpeedVol <- ggplot(dMerged, aes(x=estimatedlogVolume, y=logSpeed, fill=tAdaptation, colour=tAdaptation, shape=mediumConcentration, size=mediumConcentration)) +
  geom_errorbarh(aes(xmin=estimatedlogVolume_0025, xmax = estimatedlogVolume_0975)) +
  geom_errorbar(aes(ymin=logSpeed_0025, ymax = logSpeed_0975)) +
  geom_point(size=3, colour="black") +
  geom_abline(intercept = tentativeInterceptSpeedvsVolumeWithImposedSlope, slope = imposedSlope, color="black") +
  geom_abline(intercept = tentativeSlopeSpeedvsVolume$coefficients[1], slope = tentativeSlopeSpeedvsVolume$coefficients[2], color="darkgrey") +
  annotate("text", size=5, x=3.80, y=3.3, label= paste("best fitting slope: ", round(tentativeSlopeSpeedvsVolume$coefficients[2], 2), sep=""), hjust = 0, parse=F, colour="darkgrey") +
  annotate("text", size=5, x=3.80, y=3.2, label= paste("theoretical slope: ", round(imposedSlope, 2), sep=""), hjust = 0, parse=F, colour="black") +
  scale_x_continuous(name=expression(paste('log'[10]*'(volume)'," ",  mu, 'm'^3))) +
  scale_y_continuous(name=expression(paste('log'[10]*'(speed)'," ",  mu, 'm/s'))) +
  coord_fixed(ratio=1, xlim=c(3.75,4.75), ylim=c(2.4,3.4), expand=FALSE) +
  theme_classic(base_size = 15) +
  theme(legend.position = "none") +
  scale_colour_manual(values=c("#3B9AB2", "#EBCC2A", "#F21A00", "#FF00FF")) +
  scale_fill_manual(values=c("#3B9AB2", "#EBCC2A", "#F21A00", "#FF00FF")) +
  scale_size_manual( values = c(0.4, 0.6, 0.8) ) +
  scale_shape_manual(values=c(21, 24, 22))  # shapes for the markers
# the superposed line could be compatible with a linear scaling of 
# metabolic rate with cell volume, except in the low-nutrients conditions
# where perhaps the metabolic rate is limited by nutrients
gSpeedVol
if (saveFigures){
  ggsave(file="log10Speed_vs_log10Vol.png", dpi = 600, width = 12, height = 10, units = "cm")
  ggsave(file="log10Speed_vs_log10Vol.pdf", device=cairo_pdf, dpi = 1200, width = 24, height = 20, units = "cm")
}



# combine together the two plots into a single figure
fullFigure <- ggarrange(gMRVol, gSpeedVol,
                    ncol = 2, nrow = 1, align="v", common.legend=FALSE)
fullFigure

if (saveFigures){
  ggsave(file="figure_speed_and_MR_based_on_volume.png", dpi = 600, width = 24, height = 10, units = "cm")
  ggsave(file="figure_speed_and_MR_based_on_volume.pdf", device=cairo_pdf, dpi = 1200, width = 24, height = 10, units = "cm")
}

