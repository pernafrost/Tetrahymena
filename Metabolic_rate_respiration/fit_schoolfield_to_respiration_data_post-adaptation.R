rm(list=ls()) # clean memory







verbose = FALSE # This tells the programme to print some information for each function



#################################################################################################################
# Function to calculate the mean and the standard deviation
# for each group
#+++++++++++++++++++++++++
# data : a data frame
# varname : the name of a column containing the variable
#to be summariezed
# groupnames : vector of column names to be used as
# grouping variables
data_summary <- function(data, varname, groupnames){
  require(plyr)
  library("plotrix") # for standard error
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE),
      se = std.error(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- plyr::rename(data_sum, c("mean" = varname))
  return(data_sum)
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




# Consider changing this to true, also consider
# removing the filter on size and speed
# check why allFitResults$e_from_fit and allFitResults$e are different
saveFigures <- TRUE
combineLinesTogether <- TRUE # whether to analyse all experimental lines together or each independently
includeBootstrap <- TRUE # whether to run a bootstrap on the fitted thermal response curve

# dev.off()

# load packages
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)
library(ggplot2)

referenceTemperature <- 20 # note that this is not used; check in the code

fileName <- "~/Tetrahymena/Metabolic_rate_respiration/summary_MR_post_adaptation.csv" # file.choose() # ask the user to select a file name.

setwd(dirname(fileName))

# read the formatted and pre-processed respiration data
allExperimentResults <- read.table(file = fileName, sep = ",", header=TRUE, na.strings = c("NA", " NA"))



allTreatments <- unique(allExperimentResults$treatment)
plotList = list() # To save the plots in a list
iteration = 0
for (aaa in 1:length(allTreatments))
{
  # print(paste("Treatment: ", allTreatments[aaa]))
  allMediumConcentrations <- sort(unique(allExperimentResults$concentration[allExperimentResults$treatment == allTreatments[aaa]]))
  for (mmm in 1:length(allMediumConcentrations))
  {
    # print(paste("Medium concentration: ", allMediumConcentrations[mmm]))
    if (combineLinesTogether == FALSE)
    {
      allLines <- unique(allExperimentResults$replicate[allExperimentResults$treatment == allTreatments[aaa] & allExperimentResults$concentration == allMediumConcentrations[mmm]])
    } else {
      allLines <- "all"
    }
    for (lll in 1:length(allLines))
    {
      print(paste("TAdapt:", allTreatments[aaa], "; C:", allMediumConcentrations[mmm], "; Line: ", allLines[lll]))
      currentTitle <- paste("aaaaa_", allTreatments[aaa], "_", allMediumConcentrations[mmm], "_", allLines[lll], sep="")
      
      iteration <- iteration + 1
      
      
      if (combineLinesTogether == FALSE)
      {
        currentCondition <-subset(allExperimentResults, treatment == allTreatments[aaa] & concentration == allMediumConcentrations[mmm] & replicate ==  allLines[lll])
      } else {
        currentCondition <-subset(allExperimentResults, treatment == allTreatments[aaa] & concentration == allMediumConcentrations[mmm])
      }
      
      
      
      
      ggplot(currentCondition, aes(x=temperature, y=rate.per.cell.nW)) +
        geom_point() + 
        # geom_jitter(position = position_jitter(height = 0, width = .3)) +
        # geom_smooth(method="lm", se=TRUE, fill=NA, formula=y ~ poly(x, 3, raw=TRUE)) + 
        scale_y_continuous(name="metabolic rate (nW)") +
        scale_x_continuous(name="tested temperature", limits=c(8, 42)) +
        theme_classic(base_size = 15) +
        theme(legend.position = "none") + 
        ggtitle(paste("TAdapt:", allTreatments[aaa], "; C:", allMediumConcentrations[mmm], "; Line: ", allLines[lll])) +
        scale_color_manual(values=c("#3B9AB2", "#EBCC2A", "#F21A00")) # this is the zissou1 palette
      if (saveFigures){
        ggsave(file=paste(currentTitle, "_metabolic_rate_nW_vs_temp.png", sep=""), dpi = 600, width = 12, height = 10, units = "cm")
      }
      
      
      ###### Fit schoolfield to metabolic rate        
      
      
      if(!require(dplyr)){install.packages('dplyr')}
      d <- currentCondition %>% 
        select("replicate", "treatment", "rate.per.cell.nW", "temperature")%>%
        dplyr::rename(
          curve_id = replicate,
          growth_temp = treatment,
          rate = rate.per.cell.nW,
          temp = temperature
        )
      
      # convert Watts to eV/s
      # eVperW <- 6241506363094027800
      # d$rate <- d$rate * eVperW
      
      # convert to nanoWatt
      # d$rate <- d$rate * 10^9
      
      
      # convert Watts to J/hour
      # d$rate <- d$rate * 3600
      
      d <- d[complete.cases(d), ]
      # I added this because the get_start_vals function does not seem to work with na values
      # and when using w per cell as units the numbers are outside of the default boundaries
      
      
      # I add an artificial data point for temperature 40 degrees and rate=0
      dExtra <- d[1,]
      dExtra$temp <- 40
      dExtra$rate <- 0
      d <- rbind(d, dExtra)
      
      
      # choose model
      mod = 'sharpschoolhigh_1981'
      
      # get start vals
      start_vals <- get_start_vals(d$temp, d$rate, model_name = 'sharpeschoolhigh_1981')
      
      # get limits
      low_lims <- get_lower_lims(d$temp, d$rate, model_name = 'sharpeschoolhigh_1981')
      upper_lims <- get_upper_lims(d$temp, d$rate, model_name = 'sharpeschoolhigh_1981')
      
      start_vals
      low_lims
      upper_lims
      
      
      
      
      
      
      
      
      
      
      
      # fit model
      fit <- nls_multstart(rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 20),
                           data = d,
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
      
      calculatedFitParameters$tAdaptation <- allTreatments[aaa]
      calculatedFitParameters$line <- allLines[lll]
      calculatedFitParameters$mediumConcentration <- allMediumConcentrations[mmm]
      
      # rate at the adaptation temperature:
      r_tadapt <- (augment(fit, newdata=data.frame(temp=as.numeric(allTreatments[aaa]))))$.fitted
      calculatedFitParameters <- cbind(calculatedFitParameters, data.frame(r_tref=coef(fit)[1], e_from_fit=coef(fit)[2], eh_from_fit=coef(fit)[3], th=coef(fit)[4], r_tadapt = r_tadapt))
      
      print(paste("TAdapt:", allTreatments[aaa], "; C:", allMediumConcentrations[mmm], "; Line: ", allLines[lll]))
      print(calculatedFitParameters)
      
      print(fit)
      
      
      if (includeBootstrap)
      {
        library(minpack.lm)
        # refit model using nlsLM
        fit_nlsLM <- minpack.lm::nlsLM(rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 20),
                                       data = d,
                                       start = coef(fit),
                                       lower = get_lower_lims(d$temp, d$rate, model_name = 'sharpeschoolhigh_1981'),
                                       upper = get_upper_lims(d$temp, d$rate, model_name = 'sharpeschoolhigh_1981'),
                                       weights = rep(1, times = nrow(d)))        
        
        library(car)
        # bootstrap using case resampling
        # boot1 <- Boot(fit_nlsLM, method = 'case')
        boot1 <- Boot(fit_nlsLM, method = 'residual')
        
        # create predictions of each bootstrapped model
        boot1_preds <- boot1$t %>%
          as.data.frame() %>%
          drop_na() %>%
          dplyr::mutate(iter = 1:n()) %>%
          group_by_all() %>%
          do(data.frame(temp = seq(12.5, 30, 0.5))) %>% 
          ungroup() %>%
          dplyr::mutate(pred = sharpeschoolhigh_1981(temp, r_tref, e, eh, th, tref = 20))
        
        # calculate bootstrapped confidence intervals
        boot1_conf_preds <- group_by(boot1_preds, temp) %>%
          dplyr::summarise(conf_lower = quantile(pred, 0.025), conf_upper = quantile(pred, 0.975)) %>%
          ungroup()
        
        # exclude improbable numbers in which the t_max is less than 25 degrees
        # as when the max is too low, then very few data points contribute to the actual
        # activation energy (most of the curve is deactivation)
        boot1$t <- boot1$t[boot1$t[,4]>25,]
        
        CI_r_tref <- quantile(boot1$t[,1], c(0.025,0.975), na.rm=TRUE)
        CI_e_from_fit <- quantile(boot1$t[,2], c(0.025,0.975), na.rm=TRUE)
        se_e_from_fit <- sd(boot1$t[,2], na.rm=TRUE) / sqrt(sum(is.finite(boot1$t[,2])))
        sd_e_from_fit <- sd(boot1$t[,2], na.rm=TRUE)
        sd_r_tref <- sd(boot1$t[,1], na.rm=TRUE)
        
        calculatedFitParameters <- cbind(calculatedFitParameters, data.frame(r_tref_0025=CI_r_tref[1], r_tref_0975=CI_r_tref[2], e_from_fit_0025=CI_e_from_fit[1], e_from_fit_0975=CI_e_from_fit[2], sd_e_from_fit, sd_r_tref))
        
      }
      
      
      if (iteration == 1)
      {
        allFitResults <- calculatedFitParameters[0,]
      }
      
      allFitResults <- rbind(allFitResults, calculatedFitParameters)
      
      # predict new data
      # new_data <- data.frame(temp = seq(min(d$temp), max(d$temp), 0.5))
      new_data <- data.frame(temp = seq(12.5, 30, 0.5))
      preds <- augment(fit, newdata = new_data)
      
      
      # plot data and model fit
      ggplot(d, aes(temp, rate)) +
        geom_line(aes(temp, .fitted), preds, col = 'blue', size=2) +
        geom_point(size=6) +
        theme_classic(base_size = 15) +
        scale_x_continuous(name="Temp. (°C)") +
        scale_y_continuous(name="Rate (nW)", sec.axis = sec_axis( trans=~.*W_to_respiration_rate(1e-9), name=expression(paste("Rate (",  mu, "mol[O2]/cell/min", ")")))) +
        ggtitle(paste("TAdapt:", allTreatments[aaa], "; C:", allMediumConcentrations[mmm], "; Line: ", allLines[lll]))
      if (saveFigures){
        ggsave(file=paste(currentTitle, "_schoolfield.png", sep=""), dpi = 600, width = 12, height = 10, units = "cm")
      }
      
      
      # Summarize the data in order to plot the data with error bar
      d2 <- data_summary(d, varname="rate", 
                         groupnames=c("growth_temp", "temp"))
      # head(d2)
      
      lineColours = c("#3B9AB2", "#EBCC2A", "#F21A00", "#FF00FF")
      markerColours = c("#A3A3A3", "#666666", "#000000", "#FF00FF")
      markerShapes = c(16,17,15,4)
      
      # plot data and model fit
      thisPlot <- ggplot(d2, aes(temp, rate)) +
        geom_line(aes(temp, .fitted), preds, col = lineColours[aaa], size=2) +
        geom_point(size=3, col=markerColours[mmm], shape=markerShapes[mmm]) +
        annotate("text", size=5, x=10.5, y=3.15, label= paste("E=", round(calculatedFitParameters$e_from_fit, 2), "eV", sep=""), hjust = 0, parse=F) +
        geom_errorbar(aes(ymin=rate-sd, ymax=rate+sd), width=1.3, col=markerColours[mmm]) + 
        theme_classic(base_size = 14) +
        scale_x_continuous(name="Temp. (°C)",  limits=c(10, 32.5)) +
        scale_y_continuous(name="Rate (nW)",  limits=c(0, 3.6), sec.axis = sec_axis( trans=~.*W_to_respiration_rate(1e-9), name=expression(paste("Rate (",  mu, "mol[O2]/cell/min", ")")))) # +
      # ggtitle(paste("TAdapt:", allTreatments[aaa], "; C:", allMediumConcentrations[mmm], "; Line: ", allLines[lll])) +
      # scale_fill_manual(values= alpha(c("#3B9AB2", "#EBCC2A", "#F21A00", "#FF00FF")), 0.9) + # this is the zissou1 palette
      # scale_color_manual(values=c("#3B9AB2", "#EBCC2A", "#F21A00", "#FF00FF")) # this is the zissou1 palette
      thisPlot
      if (saveFigures){
        ggsave(file=paste(currentTitle, "_schoolfield_errorbar.png", sep=""), dpi = 600, width = 12, height = 10, units = "cm")
      }
      
      
      
      if (includeBootstrap)
      {
        # plot bootstrapped CIs
        ggplot(d2, aes(temp, rate)) +
          geom_line(aes(temp, .fitted), preds, col = lineColours[aaa], size=2) +
          geom_line(aes(temp, pred, group = iter), boot1_preds, col = lineColours[aaa], alpha = 0.007) +
          geom_point(size=3, col=markerColours[mmm], shape=markerShapes[mmm]) +
          annotate("text", size=5, x=10.5, y=3.15, label= paste("E=", round(calculatedFitParameters$e_from_fit, 2), "eV", sep=""), hjust = 0, parse=F) +
          geom_errorbar(aes(ymin=rate-sd, ymax=rate+sd), width=1.3, col=markerColours[mmm]) + 
          theme_classic(base_size = 14) +
          scale_x_continuous(name="Temp. (°C)",  limits=c(10, 32.5)) +
          scale_y_continuous(name="Rate (nW)",  limits=c(0, 3.6), sec.axis = sec_axis( trans=~.*W_to_respiration_rate(1e-9), name=expression(paste("Rate (",  mu, "mol[O2]/cell/min", ")"))))
        if (saveFigures){
          ggsave(file=paste(currentTitle, "_schoolfield_errorbar_and_bootstrap.png", sep=""), dpi = 600, width = 12, height = 10, units = "cm")
        }
        
        # plot bootstrapped CIs
        ggplot() +
          geom_line(aes(temp, .fitted), preds, col = lineColours[aaa], size=2) +
          geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot1_conf_preds, fill = lineColours[aaa], alpha = 0.3) +
          geom_point(data=d2, aes(temp, rate), size=3, col=markerColours[mmm], shape=markerShapes[mmm]) +
          annotate("text", size=5, x=10.5, y=3.15, label= paste("E=", round(calculatedFitParameters$e_from_fit, 2), "eV", sep=""), hjust = 0, parse=F) +
          geom_errorbar(data=d2, aes(x=temp, ymin=rate-sd, ymax=rate+sd), width=1.3, col=markerColours[mmm]) + 
          theme_classic(base_size = 14) +
          scale_x_continuous(name="Temp. (°C)",  limits=c(10, 32.5))  +
          scale_y_continuous(name="Rate (nW)", sec.axis = sec_axis( trans=~.*W_to_respiration_rate(1e-9), name=expression(paste("Rate (",  mu, "mol[O2]/cell/min", ")"))))
        if (saveFigures){
          ggsave(file=paste(currentTitle, "_schoolfield_errorbar_and_bootstrap2.png", sep=""), dpi = 600, width = 12, height = 10, units = "cm")
        }
        
      }
      
      
      
      plotList[[iteration]] <- thisPlot
      
    }
  }
}

print(allFitResults)

if (combineLinesTogether)
{
  write.table(allFitResults, "allFitResults_metabolic_rate_post-adaptation_all.csv", append = FALSE, sep = ", ", row.names = FALSE)
} else {
  write.table(allFitResults, "allFitResults_metabolic_rate_post-adaptation.csv", append = FALSE, sep = ", ", row.names = FALSE)
}



# Make a figure with all the plots
library(ggpubr)
# combine together the two plots into a single figure
fullFigure <- ggarrange(plotList[[1]], plotList[[2]], plotList[[3]],
                        plotList[[4]], plotList[[5]], plotList[[6]],
                        plotList[[7]], plotList[[8]], plotList[[9]],
                        # labels = c("15-50", "15-100", "15-200", "20-50", "20-100", "20-200", "25-50", "25-100", "25-200"),
                        ncol = 3, nrow = 3, common.legend=TRUE)

# annotate_figure(fullFigure, 
#                 top = "50%                    100%                    200%",
#                 left = "15°C                   20°C                   25°C")
fullFigure
ggsave(file="figure_metabolic_thermal_response_post_adaptation.png", dpi = 600, width = 24, height = 20, units = "cm")
ggsave(file="figure_metabolic_thermal_response_post_adaptation.eps", device="eps", dpi = 1200, width = 24, height = 20, units = "cm")
library(Cairo)
ggsave(file="figure_metabolic_thermal_response_post_adaptation.pdf", device=cairo_pdf, dpi = 1200, width = 24, height = 20, units = "cm")


allFitResults <- subset(allFitResults, tAdaptation != "control")
allFitResults$tAdaptation <- as.numeric(allFitResults$tAdaptation)

# plot data and model fit
ggplot(allFitResults, aes(x=tAdaptation, y=e_from_fit, fill=factor(mediumConcentration + tAdaptation*100), color= factor(mediumConcentration), shape=factor(mediumConcentration), label=factor(line))) +
  # geom_jitter(size=6, position = position_jitter(height = 0, width = .7)) +
  geom_point(size=6) + 
  # geom_text(color="red", size=4) + 
  scale_color_manual(values=c("#A3A3A3", "#666666", "#000000")) +
  scale_fill_manual(values=c("#8FCDDC", "#3B9AB2", "#18434E", "#F2DD70", "#EBCC2A", "#463C07", "#FF7C6C", "#F21A00", "#660B00")) +
  scale_y_continuous(name="activation energy (eV)", limits=c(0, 1)) +
  scale_x_continuous(name="Adaptation temperature (ºC)", limits=c(12.5, 27.5), breaks=c(15, 20, 25)) +
  scale_shape_manual(values=c(21, 24, 22)) + # shapes for the markers
  theme_classic(base_size=18) +
  theme(legend.position = "none") +
  labs(color = "Conc.") #  + # this specifies a custom legend
# ggtitle("Activation energy (eV)")

if (saveFigures){
  if (combineLinesTogether)
  {
    ggsave(file="figure_activation_energy_all.png", dpi = 600, width = 12, height = 10, units = "cm")
    ggsave(file="figure_activation_energy_all.eps", device="eps", dpi = 1200, width = 12, height = 10, units = "cm")
    library(Cairo)
    ggsave(file="figure_activation_energy_all.pdf", device=cairo_pdf, dpi = 1200, width = 12, height = 10, units = "cm")
  } else {
    ggsave(file="aaaa_activation_energy.png", dpi = 600, width = 12, height = 10, units = "cm")
  }
}




if (includeBootstrap)
{
  # plot data and model fit
  ggplot(allFitResults, aes(x=tAdaptation, y=e_from_fit, fill=factor(mediumConcentration + tAdaptation*100), color= factor(mediumConcentration), shape=factor(mediumConcentration), label=factor(line))) +
    # geom_jitter(size=6, position = position_jitter(height = 0, width = .7)) +
    geom_errorbar(aes(ymin=e_from_fit - sd_e_from_fit, ymax=e_from_fit + sd_e_from_fit), width=0, position=position_dodge(width=2)) + 
    geom_point(size=6, position=position_dodge(width=2)) + 
    # geom_errorbar(aes(ymin=e_from_fit_0025, ymax=e_from_fit_0975), width=0, position=position_dodge(width=2)) + 
    scale_color_manual(values=c("#A3A3A3", "#666666", "#000000")) +
    scale_fill_manual(values=c("#8FCDDC", "#3B9AB2", "#18434E", "#F2DD70", "#EBCC2A", "#463C07", "#FF7C6C", "#F21A00", "#660B00")) +
    scale_shape_manual(values=c(21, 24, 22)) + # shapes for the markers
    scale_y_continuous(name="Activation energy (eV)") + #
    scale_x_continuous(name="Adaptation temperature (ºC)", limits=c(12.5, 27.5), breaks=c(15, 20, 25)) +
    coord_cartesian(ylim = c(0, 1)) +
    theme_classic(base_size=18) +
    theme(legend.position = "none") +
    labs(color = "Conc.") #  + # this specifies a custom legend
  # ggtitle("Activation energy (eV)")
  
  if (saveFigures){
    if (combineLinesTogether)
    {
      ggsave(file="figure_activation_energy_all_bootstrap.png", dpi = 600, width = 12, height = 10, units = "cm")
      ggsave(file="figure_activation_energy_all_bootstrap.eps", device="eps", dpi = 1200, width = 12, height = 10, units = "cm")
      library(Cairo)
      ggsave(file="figure_activation_energy_all_bootstrap.pdf", device=cairo_pdf, dpi = 1200, width = 12, height = 10, units = "cm")
    } else {
      ggsave(file="aaaa_activation_energy_bootstrap.png", dpi = 600, width = 12, height = 10, units = "cm")
    }
  }
  
}

# 
# # plot data and model fit
# ggplot(allFitResults, aes(x=tAdaptation, y=topt, color=factor(mediumConcentration), shape=factor(mediumConcentration), label=factor(line))) +
#   geom_point(size=4) +
#   # geom_text(color="red", size=4) + 
#   scale_color_manual(values=c("#A3A3A3", "#666666", "#000000")) +
#   scale_fill_manual(values=c("#A3A3A3", "#666666", "#000000")) +
#   scale_y_continuous(name="Optimal temperature (ºC)") +
#   scale_x_continuous(name="Adaptation temperature (ºC)", limits=c(12.5, 27.5), breaks=c(15, 20, 25)) +
#   theme_classic(base_size=12) +
#   theme(legend.position = "none") +
#   # labs(color = "Conc.") + # this specifies a custom legend
#   ggtitle("Optimal temperature")
# 
# if (saveFigures){
#   if (combineLinesTogether)
#   {
#     ggsave(file="aaaa_optimal_temp_vs_temp_all.png", dpi = 600, width = 12, height = 10, units = "cm")
#   } else {
#     ggsave(file="aaaa_optimal_temp_vs_temp.png", dpi = 600, width = 12, height = 10, units = "cm")
#   }
# }







# plot data and model fit
ggplot(allFitResults, aes(x=tAdaptation, y=r_tref, fill=factor(mediumConcentration + tAdaptation*100), color=factor(mediumConcentration), shape=factor(mediumConcentration), label=factor(line))) +
  geom_point(size=4) +
  # geom_text(color="red", size=4) + 
  scale_color_manual(values=c("#A3A3A3", "#666666", "#000000")) +
  scale_fill_manual(values=c("#8FCDDC", "#3B9AB2", "#18434E", "#F2DD70", "#EBCC2A", "#463C07", "#FF7C6C", "#F21A00", "#660B00")) +
  scale_shape_manual(values=c(21, 24, 22)) + # shapes for the markers
  scale_y_continuous(name=paste("M. R. (nW) at T=", referenceTemperature, "°C"), sec.axis = sec_axis( trans=~.*W_to_respiration_rate(1e-9), name=expression(paste("Rate (",  mu, "mol[O2]/cell/min", ")")))) +
  # scale_y_continuous(name=paste("M. R. (nW) at T=", referenceTemperature, "°C")) +
  scale_x_continuous(name="Adaptation temperature (ºC)", limits=c(12.5, 27.5), breaks=c(15, 20, 25)) +
  theme_classic(base_size=18) +
  theme(legend.position = "none") # +
# labs(color = "Conc.") + # this specifies a custom legend
# ggtitle(paste("M.R. at reference temperature T= ", referenceTemperature, "°C", sep=""))

if (saveFigures){
  if (combineLinesTogether)
  {
    ggsave(file="figure_metabolic_rate_vs_tref_all.png", dpi = 600, width = 12, height = 10, units = "cm")
    ggsave(file="figure_metabolic_rate_vs_tref_all.eps", device="eps", dpi = 1200, width = 12, height = 10, units = "cm")
    # quartz(type = 'pdf', file="aaaa_metabolic_thermal_response_post_adaptation.pdf")
    library(Cairo)
    ggsave(file="figure_metabolic_rate_vs_tref_all.pdf", device=cairo_pdf, dpi = 1200, width = 12, height = 10, units = "cm")
  } else {
    ggsave(file="aaaa_metabolic_rate_vs_tref.png", dpi = 600, width = 12, height = 10, units = "cm")
  }
}



if (includeBootstrap)
{
  
  # plot data and model fit
  ggplot(allFitResults, aes(x=tAdaptation, y=r_tref, color=factor(mediumConcentration), fill=factor(mediumConcentration + tAdaptation*100),  shape=factor(mediumConcentration), label=factor(line))) +
    geom_errorbar(aes(ymin=r_tref - sd_r_tref, ymax=r_tref + sd_r_tref), width=0, position=position_dodge(width=2)) + 
    geom_point(size=6, position=position_dodge(width=2)) + 
    # geom_errorbar(aes(ymin=e_from_fit_0025, ymax=e_from_fit_0975), width=0, position=position_dodge(width=2)) + 
    # geom_text(color="red", size=4) + 
    scale_color_manual(values=c("#A3A3A3", "#666666", "#000000")) +
    scale_fill_manual(values=c("#8FCDDC", "#3B9AB2", "#18434E", "#F2DD70", "#EBCC2A", "#463C07", "#FF7C6C", "#F21A00", "#660B00")) +
    scale_shape_manual(values=c(21, 24, 22)) + # shapes for the markers
    scale_y_continuous(name="Rate at T=20°C (nW)") +
    # scale_y_continuous(name="Rate (nW)", sec.axis = sec_axis( trans=~.*W_to_respiration_rate(1e-9), name=expression(paste("Rate (",  mu, "mol[O2]/cell/min", ")")))) +
    scale_x_continuous(name="Adaptation temperature (ºC)", limits=c(12.5, 27.5), breaks=c(15, 20, 25)) +
    theme_classic(base_size=18) +
    theme(legend.position = "none")#  +
  # labs(color = "Conc.") + # this specifies a custom legend
  # ggtitle(paste("M.R. at reference temperature T= ", referenceTemperature, "°C", sep=""))
  
  if (saveFigures){
    if (combineLinesTogether)
    {
      ggsave(file="figure_metabolic_rate_vs_tref_all_bootstrap.png", dpi = 600, width = 12, height = 10, units = "cm")
      ggsave(file="figure_metabolic_rate_vs_tref_all_bootstrap.eps", device="eps", dpi = 1200, width = 12, height = 10, units = "cm")
      library(Cairo)
      ggsave(file="figure_metabolic_rate_vs_tref_all_bootstrap.pdf", device=cairo_pdf, dpi = 1200, width = 12, height = 10, units = "cm")
    } else {
      ggsave(file="aaaa_metabolic_rate_vs_tref_bootstrap.png", dpi = 600, width = 12, height = 10, units = "cm")
    }
  }
}

# plot data and model fit
ggplot(allFitResults, aes(x=tAdaptation, y=r_tadapt, color=factor(mediumConcentration), fill=factor(mediumConcentration + tAdaptation*100),  shape=factor(mediumConcentration), label=factor(line))) +
  geom_point(size=6, position=position_dodge(width=2)) + 
  # geom_text(color="red", size=4) + 
  scale_color_manual(values=c("#A3A3A3", "#666666", "#000000")) +
  scale_fill_manual(values=c("#8FCDDC", "#3B9AB2", "#18434E", "#F2DD70", "#EBCC2A", "#463C07", "#FF7C6C", "#F21A00", "#660B00")) +
  scale_shape_manual(values=c(21, 24, 22)) + # shapes for the markers
  scale_y_continuous(name="Rate (nW)", sec.axis = sec_axis( trans=~.*W_to_respiration_rate(1e-9), name=expression(paste("Rate (",  mu, "mol[O2]/cell/min", ")")))) +
  scale_x_continuous(name="Adaptation temperature (ºC)", limits=c(12.5, 27.5), breaks=c(15, 20, 25)) +
  theme_classic(base_size=18) +
  theme(legend.position = "none") +
  # labs(color = "Conc.") + # this specifies a custom legend
  ggtitle("M.R. at adaptation temperature")

if (saveFigures){
  if (combineLinesTogether)
  {
    ggsave(file="aaaa_metabolic_rate_vs_tadapt_all.png", dpi = 600, width = 12, height = 10, units = "cm")
  } else {
    ggsave(file="aaaa_metabolic_rate_vs_tadapt.png", dpi = 600, width = 12, height = 10, units = "cm")
  }
}



lineColours = c("#3B9AB2", "#EBCC2A", "#F21A00", "#FF00FF")
lineTypes <- c("twodash", "dashed", "solid", "longdash")
# Now try to plot the fit only
allFitData <- data.frame()
for (iii in 1:dim(allFitResults)[1])
{
  predictedData <- data.frame(temp = seq(10, 30, 0.5))
  predictedData$rate <- sharpeschoolhigh_1981(temp = predictedData$temp, r_tref = allFitResults$r_tref[iii],e = allFitResults$e_from_fit[iii], eh = allFitResults$eh_from_fit[iii], th = allFitResults$th[iii], tref = 20)
  predictedData$tAdaptation <- allFitResults$tAdaptation[iii]
  predictedData$mediumConcentration <- allFitResults$mediumConcentration[iii]
  predictedData$line <- allFitResults$line[iii]
  allFitData <- rbind(allFitData, predictedData)
}

plotG1 <- ggplot(allFitData, aes(x=temp, y=rate, color=as.factor(tAdaptation), size=as.factor(mediumConcentration), alpha=mediumConcentration)) +
  geom_line() +
  scale_x_continuous(name="Temp. (°C)",  limits=c(10, 25)) +
  # scale_y_continuous(name="Rate (nW)",  limits=c(0.1, 3), trans = 'log10', sec.axis = sec_axis( trans=~.*W_to_respiration_rate(1e-9), name=expression(paste("Rate (",  mu, "mol[O2]/cell/min", ")"))))  +
  scale_y_continuous(name="Rate (nW)",  limits=c(0.1, 1), trans = 'log10')  +
  theme_classic(base_size=18) +
  scale_color_manual(values=lineColours) +
  scale_linetype_manual(values=lineTypes) + 
  scale_alpha_continuous(range=c(0.5, 1)) +
  scale_size_manual( values = c(0.8, 1.4, 2) ) +
  theme(legend.position = "none") +
  labs(color = "Conc.") 
plotG1

if (saveFigures){
  if (combineLinesTogether)
  {
    ggsave(file="figure_fitted_rate_vs_temp_all.png", dpi = 600, width = 12, height = 10, units = "cm")
    ggsave(file="figure_fitted_rate_vs_temp_all.eps", device="eps", dpi = 1200, width = 12, height = 10, units = "cm")
    library(Cairo)
    ggsave(file="figure_fitted_rate_vs_temp_all.pdf", device=cairo_pdf, dpi = 1200, width = 12, height = 10, units = "cm")
  } else {
    ggsave(file="figure_fitted_rate_vs_temp.png", dpi = 600, width = 12, height = 10, units = "cm")
  }
}


plotG1 <- ggplot(allFitData, aes(x=temp, y=rate, color=as.factor(tAdaptation), size=as.factor(mediumConcentration))) +
  geom_line() +
  geom_point(size=6, data=allFitResults, color="black", aes(x=tAdaptation, y=r_tadapt, fill=factor(tAdaptation),  shape=factor(mediumConcentration))) +
  scale_x_continuous(name="Temp. (°C)",  limits=c(10, 25)) +
  # scale_y_continuous(name="Rate (nW)",  limits=c(0.1, 3), trans = 'log10', sec.axis = sec_axis( trans=~.*W_to_respiration_rate(1e-9), name=expression(paste("Rate (",  mu, "mol[O2]/cell/min", ")"))))  +
  scale_y_continuous(name="Rate (nW)",  limits=c(0.3, 3), trans = 'log10')  +
  theme_classic(base_size=18) +
  scale_color_manual(values=lineColours) +
  scale_shape_manual(values=c(21, 24, 22)) + # shapes for the markers
  scale_fill_manual(values=lineColours) +
  scale_linetype_manual(values=lineTypes) + 
  scale_alpha_continuous(range=c(0.3, 0.5)) +
  scale_size_manual( values = c(0.4, 0.7, 1) ) +
  theme(legend.position = "none") +
  labs(color = "Conc.") 
plotG1

if (saveFigures){
  if (combineLinesTogether)
  {
    ggsave(file="figure_fitted_rate_vs_tadapt_all.png", dpi = 600, width = 12, height = 10, units = "cm")
    ggsave(file="figure_fitted_rate_vs_tadapt_all.eps", device="eps", dpi = 1200, width = 12, height = 10, units = "cm")
    library(Cairo)
    ggsave(file="figure_fitted_rate_vs_tadapt_all.pdf", device=cairo_pdf, dpi = 1200, width = 12, height = 10, units = "cm")
  } else {
    ggsave(file="figure_fitted_rate_vs_tadapt.png", dpi = 600, width = 12, height = 10, units = "cm")
  }
}





experimentalResultsAtTAdapt <- subset(allExperimentResults, treatment == temperature)
experimentalResultsAtTAdapt$treatment <- as.numeric(experimentalResultsAtTAdapt$treatment)

if (combineLinesTogether)
{
  
  # Summarize the data in order to plot the data with error bar
  dSum <- data_summary(experimentalResultsAtTAdapt, varname="rate.per.cell.nW", 
                       groupnames=c("treatment", "concentration"))
  
  
  # plot data and model fit
  ggplot(dSum, aes(x=treatment, y=rate.per.cell.nW, color=factor(concentration), shape=factor(concentration))) +
    geom_point(size=6, position=position_dodge(width=2)) + 
    # geom_errorbar(aes(ymin=e_from_fit_0025, ymax=e_from_fit_0975), width=0, position=position_dodge(width=2)) + 
    geom_errorbar(aes(ymin=rate.per.cell.nW - sd, ymax=rate.per.cell.nW + sd), width=0, position=position_dodge(width=2)) + 
    # geom_text(color="red", size=4) + 
    scale_color_manual(values=c("#A3A3A3", "#666666", "#000000")) +
    scale_fill_manual(values=c("#A3A3A3", "#666666", "#000000")) +
    scale_y_continuous(name="Rate (nW)", sec.axis = sec_axis( trans=~.*W_to_respiration_rate(1e-9), name=expression(paste("Rate (",  mu, "mol[O2]/cell/min", ")")))) +
    scale_x_continuous(name="Adaptation temperature (ºC)", limits=c(12.5, 27.5), breaks=c(15, 20, 25)) +
    theme_classic(base_size=18) +
    theme(legend.position = "none") +
    # labs(color = "Conc.") + # this specifies a custom legend
    ggtitle(paste("M.R. at adaptation temperature T"))
  
  # plot data and model fit
  ggplot(dSum, aes(x=treatment, y=rate.per.cell.nW, color=factor(concentration), fill=factor(concentration + treatment*100), shape=factor(concentration))) +
    geom_errorbar(aes(ymin=rate.per.cell.nW - sd, ymax=rate.per.cell.nW + sd), width=0, position=position_dodge(width=2)) + 
    geom_point(size=6, position=position_dodge(width=2)) + 
    # geom_errorbar(aes(ymin=e_from_fit_0025, ymax=e_from_fit_0975), width=0, position=position_dodge(width=2)) + 
    # geom_text(color="red", size=4) + 
    scale_color_manual(values=c("#A3A3A3", "#666666", "#000000")) +
    scale_fill_manual(values=c("#8FCDDC", "#3B9AB2", "#18434E", "#F2DD70", "#EBCC2A", "#463C07", "#FF7C6C", "#F21A00", "#660B00")) +
    scale_shape_manual(values=c(21, 24, 22)) + # shapes for the markers
    scale_y_continuous(name="Rate (nW)", sec.axis = sec_axis( trans=~.*W_to_respiration_rate(1e-9), name=expression(paste("Rate (",  mu, "mol[O2]/cell/min", ")")))) +
    scale_x_continuous(name="Adaptation temperature (ºC)", limits=c(12.5, 27.5), breaks=c(15, 20, 25)) +
    theme_classic(base_size=18) +
    theme(legend.position = "none") +
    # labs(color = "Conc.") + # this specifies a custom legend
    ggtitle(paste("M.R. at adaptation temperature T"))
  
  
  if (saveFigures){
    ggsave(file="aaaa_metabolic_rate_at_tadapt_only_all.png", dpi = 600, width = 12, height = 10, units = "cm")
  }
  
} else {
  # Metabolic rate at the adaptation temperature without Schoolfield fitting (only the original data)
  ggplot(experimentalResultsAtTAdapt, aes(x=treatment, y=rate.per.cell.nW, color=factor(concentration), fill=factor(concentration + treatment*100),  shape=factor(concentration), label=factor(replicate))) +
    geom_jitter(size=4, position = position_jitter(height = 0, width = .3)) +
    scale_color_manual(values=c("#A3A3A3", "#666666", "#000000")) +
    scale_fill_manual(values=c("#8FCDDC", "#3B9AB2", "#18434E", "#F2DD70", "#EBCC2A", "#463C07", "#FF7C6C", "#F21A00", "#660B00")) +
    scale_shape_manual(values=c(21, 24, 22)) + # shapes for the markers
    scale_y_continuous(name="Rate (nW)", sec.axis = sec_axis( trans=~.*W_to_respiration_rate(1e-9), name=expression(paste("Rate (",  mu, "mol[O2]/cell/min", ")")))) +
    scale_x_continuous(name="Adaptation temperature (ºC)", limits=c(12.5, 27.5), breaks=c(15, 20, 25)) +
    theme_classic(base_size=18) +
    theme(legend.position = "none") +
    # labs(color = "Conc.") + # this specifies a custom legend
    ggtitle("M.R. at adaptation temperature")
  
  if (saveFigures){
    ggsave(file="aaaa_metabolic_rate_at_tadapt_only.png", dpi = 600, width = 12, height = 10, units = "cm")
  }
  
  
  
