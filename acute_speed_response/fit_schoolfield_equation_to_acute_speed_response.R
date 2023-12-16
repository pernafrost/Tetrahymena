rm(list=ls()) # clean memory

# detach all packages
# lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)
# invisible(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only=TRUE, unload=TRUE))

#################################################################################################################
# Function to calculate the mean and the standard deviation
# for each group (useful to combine lines together)
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


# Consider changing this to true, also consider
# removing the filter on size and speed
# check why allFitResults$e_from_fit and allFitResults$e are different
saveFigures <- TRUE
combineLinesTogether <- TRUE # whether to analyse all experimental lines together or each independently
includeBootstrap <- TRUE # whether to run a bootstrap on the fitted thermal response curve

# dev.off()

# load packages
# remotes::install_github("padpadpadpad/rTPC")
library(rTPC)
# library(plyr)
# library(dplyr)
library(nls.multstart)
library(broom)
library(tidyverse)
library(ggplot2)
if (saveFigures == TRUE)
{        
  library(Cairo)
}

referenceTemperature <- 20 # note that this is not used; check in the code
skipMotherCulture <- TRUE # whether to skip the mother culture from this analysis

fileName <- "~/Tetrahymena/acute_speed_response/track_analysis_results_individual_particles.csv" # file.choose() # ask the user to select a file name.

setwd(dirname(fileName))
allExperimentResults <- read.table(file = fileName, sep = ",", header=TRUE, na.strings = c("NA", " NA"))

library("stringr")
# Extract information from the string in the "fileName" field
# newColumns <- str_split_fixed(as.character(str_trim(tolower(allExperimentResults$fileName))), "/adapted|_tested|_group|diluted", 5)[,2:4]

newColumns <- str_split_fixed(as.character(str_trim(tolower(basename(allExperimentResults$fileName)))), "_|_tested_", n=Inf)[,c(2,3,4,6)]

# remove c for degrees celsius and remove bis from the line
newColumns[,4] <- str_replace(newColumns[,4], "c", "")
newColumns[,3] <- str_replace(newColumns[,3], "bis", "")
newColumns[,4] <- str_replace(newColumns[,4], "bis", "")

# class(newColumns) <- "numeric"
newColumns <- data.frame(tAdaptation = as.numeric(newColumns[,1]), 
                         mediumConcentration = as.numeric(newColumns[,2]),
                         tTest = as.numeric(newColumns[,4]), 
                         line = newColumns[,3])
# colnames(newColumns) <- c("tAdaptation", "tTest", "line")

# add the new columns to the original table
allExperimentResults <- cbind(allExperimentResults, newColumns, deparse.level = 1, stringsAsFactors = default.stringsAsFactors())
# sapply(allExperimentResults, class)

allExperimentResults$logSpeed <- log10(allExperimentResults$medianSpeed)

# allExperimentResults$sqrtSpeed <- sqrt(allExperimentResults$medianSpeed)

# calculate drag power
# Here I assume the approximate equation for the drag coefficient:
# C_D = 1.2 * pi * rs * (4 + E)
# where rs is the small radius and E is the elongation
# and then the Stokes power is C_D * mu * U^2 
# where mu is the viscosity and U is speed
# I convert from um to m by multiplying times 10^-6
allExperimentResults$dragCoefficient <- 1.2 * pi * allExperimentResults$medianAEllipse * 10^-6 /2 * (4 + allExperimentResults$medianElongation)

allExperimentResults$tTestK <- allExperimentResults$tTest + 273.15

viscA = 1.856* 10^-14 # for a measure in Pascal * s
viscB = 4209 # K 
viscC = 0.04527 # K^-1
viscD = -3.366 * 10^-5 # K^-2

allExperimentResults$waterViscosity <- viscA * exp(viscB/allExperimentResults$tTestK + viscC * allExperimentResults$tTestK + viscD *allExperimentResults$tTestK^2)

allExperimentResults$dragPower <- allExperimentResults$dragCoefficient * allExperimentResults$waterViscosity * (allExperimentResults$medianSpeed*10^-6)^2

allExperimentResults$lnDragPower <- log(allExperimentResults$dragPower)

# HERE to continue the Arrhenius plot
kBoltzmann <- 8.617333262145e-5 # Boltzmann constant in eV
allExperimentResults$oneOverkT <- (kBoltzmann * allExperimentResults$tTestK)^-1

oneOverkT293 <- (kBoltzmann * 293)^-1
allExperimentResults$oneOverkTrel <- oneOverkT293 - allExperimentResults$oneOverkT


names(allExperimentResults)

allTAdaptation <- sort(unique(allExperimentResults$tAdaptation))

# coefficients of fitted log(speed) vs. log(cell volume)
allFittedLogSpeedVsLogVolume <- list()

plotList = list() # To save the plots in a list
plotListExtra = list() # I have two plot lists, one for the first nine plots, to make 
# a figure panel, the remaining plots are also kept and arranged in a different list
plotCounter <- 0
plotCounterExtra <- 0
iteration = 0
for (aaa in 1:length(allTAdaptation))
{
  # print(paste("Adaptation temp.: ", allTAdaptation[aaa]))
  allMediumConcentrations <- sort(unique(allExperimentResults$mediumConcentration[allExperimentResults$tAdaptation == allTAdaptation[aaa]]))
  for (mmm in 1:length(allMediumConcentrations))
  {
    # print(paste("Medium concentration: ", allMediumConcentrations[mmm]))
    if (combineLinesTogether == FALSE)
    {
      allLines <- sort(unique(allExperimentResults$line[allExperimentResults$tAdaptation == allTAdaptation[aaa] & allExperimentResults$mediumConcentration == allMediumConcentrations[mmm]]), decreasing=TRUE)
    } else {
      allLines <- "all"
    }
    for (lll in 1:length(allLines))
    {
      print(paste("TAdapt:", allTAdaptation[aaa], "; C:", allMediumConcentrations[mmm], "; Line: ", allLines[lll]))
      currentTitle <- paste("aaaaa_", allTAdaptation[aaa], "_", allMediumConcentrations[mmm], "_", allLines[lll], sep="")
      
      if (skipMotherCulture & is.finite(as.numeric(allLines[lll])))
      {
        next
      }
      
      iteration = iteration + 1
      
      
      if (combineLinesTogether == FALSE)
      {
        currentCondition <-subset(allExperimentResults, tAdaptation == allTAdaptation[aaa] & mediumConcentration == allMediumConcentrations[mmm] & line ==  allLines[lll])
      } else {
        currentCondition <-subset(allExperimentResults, tAdaptation == allTAdaptation[aaa] & mediumConcentration == allMediumConcentrations[mmm])
        if (skipMotherCulture)
        {
          currentCondition <-subset(currentCondition, is.finite(as.numeric(currentCondition$line)) == FALSE) 
        }
      }
      
      
      
      # Picks the first plot for each temperature and concentration condition
      # except for the 20 degrees, 100% condition, for which the first plot refers
      # to the adaptation phase, while here we want the test phase
      if (lll == 1 & !is.finite(as.numeric(allLines[lll]))) { # first iteration, not mother culture
        plotCounter <- plotCounter + 1
      } else if (is.finite(as.numeric(allLines[1])) & lll==length(allLines)) { # last iteration, the first one was the mother culture
        plotCounter <- plotCounter + 1
      } else {
        plotCounterExtra <- plotCounterExtra + 1
      }
      
      
      
      # If I further filter the data:
      # keep those that are big enough and that move
      currentCondition <- subset(currentCondition, medianArea > 400 & medianSpeed > 100)
      
      # Only keep those that move straight
      # currentCondition <- subset(currentCondition, meanAbsAngle*180/pi < 15)
      
      # Only keep a percentile for each temperature group
      allTTest <- unique(currentCondition$tTest)
      # if there are multiple lines combined, calculate the percentile on each line
      allLinesInsideLoop <- unique(currentCondition$line)
      for (iii in 1:length(allLinesInsideLoop))
      {
        for (ttt in 1:length(allTTest))
        {
          currentConditionThisTemperatureAndLine <- subset(currentCondition, tTest == allTTest[ttt] & line == allLinesInsideLoop[iii])
          speedQuantile <- quantile(currentConditionThisTemperatureAndLine$medianSpeed, probs = 0.80, na.rm = TRUE)
          currentCondition$speedQuantile[currentCondition$tTest == allTTest[ttt] & currentCondition$line==allLinesInsideLoop[iii]] <- speedQuantile
        }
      }
      
      currentCondition$estimatedlogVolume <- log10(4/3 * pi* currentCondition$medianBEllipse * currentCondition$medianAEllipse * currentCondition$medianAEllipse / 8)
      
      
      currentConditionTopSpeed <- subset(currentCondition, medianSpeed >= speedQuantile)
      
      
      
      
      
      
      ggplot(currentConditionTopSpeed, aes(x=tTest, y=medianSpeed, color=line)) +
        geom_jitter(position = position_jitter(height = 0, width = .5)) + 
        # geom_jitter(position = position_jitter(height = 0, width = .3)) +
        # geom_smooth(method="lm", se=TRUE, fill=NA, formula=y ~ poly(x, 3, raw=TRUE)) + 
        scale_y_continuous(name=expression(paste("median speed (", mu, "m/s)")), limits=c(0,1000)) +
        scale_x_continuous(name="tested temperature", limits=c(8, 42)) +
        theme_classic(base_size = 15) +
        # theme(legend.position = "none") + 
        ggtitle(paste("TAdapt:", allTAdaptation[aaa], "; C:", allMediumConcentrations[mmm], "; Line: ", allLines[lll])) +
        labs(color="abs(line):")
      # scale_color_manual(values=c("#3B9AB2", "#EBCC2A", "#F21A00")) # this is the zissou1 palette
      if (saveFigures){
        ggsave(file=paste(currentTitle, "_speed_colour_by_line.png", sep=""), dpi = 600, width = 12, height = 10, units = "cm")
      }
      
      Tlabels = c("T15", "T20", "T25")
      Tcolours = c("#3B9AB2", "#EBCC2A", "#F21A00")
      Dlabels = c("50%", "100%", "200%")
      Dcolours = c("#A3A3A3", "#666666", "#000000")
      
      
      ggplot(currentCondition, aes(x=as.factor(tTest), y=medianSpeed, fill=Dlabels[aaa], colour=Tlabels[mmm])) +
        geom_violin(width=1.5, linewidth=1.5) +
        scale_y_continuous(name=expression(paste("speed (", mu, "m/s)")), limits=c(0,1000)) +
        scale_x_discrete(name="tested temperature") +
        theme_classic(base_size = 18) +
        theme(legend.position = "none") +
        # ggtitle(paste("TAdapt:", allTAdaptation[aaa], "; C:", allMediumConcentrations[mmm], "; Line: ", allLines[lll])) +
        scale_fill_manual(name="Medium conc.", breaks=Dlabels, values=Dcolours) +
        scale_colour_manual(name="T. adapt", breaks=Tlabels, values=Tcolours)
      
      if (saveFigures){
        ggsave(file=paste(currentTitle, "_speed_vs_temp_violin.png", sep=""), dpi = 600, width = 12, height = 10, units = "cm")
        ggsave(file=paste(currentTitle, "_speed_vs_temp_violin.png", sep=""), device="eps", dpi = 1200, width = 24, height = 20, units = "cm")
        ggsave(file=paste(currentTitle, "_speed_vs_temp_violin.png", sep=""), device=cairo_pdf, dpi = 1200, width = 24, height = 20, units = "cm")
      }

      
      ggplot(currentCondition, aes(x=tTest, y=medianSpeed, color=medianArea)) +
        geom_jitter(position = position_jitter(height = 0, width = .5)) + 
        # geom_jitter(position = position_jitter(height = 0, width = .3)) +
        # geom_smooth(method="lm", se=TRUE, fill=NA, formula=y ~ poly(x, 3, raw=TRUE)) + 
        scale_y_continuous(name=expression(paste("median speed (", mu, "m/s)")), limits=c(0,1000)) +
        scale_x_continuous(name="tested temperature", limits=c(8, 42)) +
        theme_classic(base_size = 15) +
        # theme(legend.position = "none") + 
        ggtitle(paste("TAdapt:", allTAdaptation[aaa], "; C:", allMediumConcentrations[mmm], "; Line: ", allLines[lll])) +
        scale_color_continuous() +
        labs(color="Area:")
      # scale_color_manual(values=c("#3B9AB2", "#EBCC2A", "#F21A00")) # this is the zissou1 palette
      if (saveFigures){
        ggsave(file=paste(currentTitle, "_speed_vs_temp_jitter.png", sep=""), dpi = 600, width = 12, height = 10, units = "cm")
      }
      
      # a plot to check how speed changes with cell size
      ggplot(currentConditionTopSpeed, aes(x=estimatedlogVolume, y=logSpeed, color=tTest)) +
        geom_point() + 
        # geom_jitter(position = position_jitter(height = 0, width = .3)) +
        geom_smooth(method="lm", se=FALSE, formula=y ~ x, colour="red", na.rm=TRUE) + 
        scale_y_continuous(name=expression(paste("log10(speed) (", mu, "m/s)"))) +
        scale_x_continuous(name=expression(paste('log'[10]*'(volume)'," ",  mu, 'm'^3))) +
        theme_classic(base_size = 15) +
        theme(legend.position = "none") + 
        ggtitle(paste("TAdapt:", allTAdaptation[aaa], "; C:", allMediumConcentrations[mmm], "; Line: ", allLines[lll])) +
        scale_color_continuous() +
        coord_fixed() +
        facet_wrap(vars(tTest))
        labs(color="T:")
      # scale_color_manual(values=c("#3B9AB2", "#EBCC2A", "#F21A00")) # this is the zissou1 palette
      if (saveFigures){
        ggsave(file=paste(currentTitle, "_speed_vs_cell_volume.png", sep=""), dpi = 600, width = 12, height = 15, units = "cm")
      }
      
        library("lmodel2")
        for (ttt in sort(unique(currentConditionTopSpeed$tTest)))
        {
        # print(ttt)
        fitSpeedVsVolume <- lmodel2(formula=estimatedlogVolume ~ logSpeed, data = subset(currentConditionTopSpeed, tTest == ttt), range.y="interval", range.x = "interval", nperm=99)
        # print(fitSpeedVsVolume$regression.results)
        # fitSpeedVsVolume <- lm(formula=estimatedlogVolume ~ logSpeed, data = subset(currentConditionTopSpeed, tTest == ttt))
        # print(coefficients(fitSpeedVsVolume))
        allFittedLogSpeedVsLogVolume[[length(allFittedLogSpeedVsLogVolume)+1]] <- fitSpeedVsVolume$regression.results$Slope
        # summary(fitSpeedVsVolume)
        }
      
      # I am running this on currentConditionTopSpeed because I think it makes more sense to do it only on the cells
      # that move fast and so are well aligned with the image.
      ggplot(currentConditionTopSpeed, aes(x=tTest, y=medianArea, fill=factor(tTest))) +
        geom_boxplot(width=0.6, outlier.shape = NA) + # add a box plot
        scale_y_continuous("median area") +
        scale_x_continuous(name="tested temperature", limits=c(8, 42)) +
        theme_classic(base_size = 15) +
        theme(legend.position = "none") + 
        ggtitle(paste("TAdapt:", allTAdaptation[aaa], "; C:", allMediumConcentrations[mmm], "; Line: ", allLines[lll])) +
        scale_color_continuous()
      if (saveFigures){
        ggsave(file=paste(currentTitle, "_area_vs_temp.png", sep=""), dpi = 600, width = 12, height = 10, units = "cm")
      }
      
      ggplot(currentConditionTopSpeed, aes(x=tTest, y=medianElongation, fill=factor(tTest))) +
        geom_boxplot(width=0.6, outlier.shape = NA) + # add a box plot
        scale_y_continuous("elongation") +
        scale_x_continuous(name="tested temperature", limits=c(8, 42)) +
        theme_classic(base_size = 15) +
        theme(legend.position = "none") + 
        ggtitle(paste("TAdapt:", allTAdaptation[aaa], "; C:", allMediumConcentrations[mmm], "; Line: ", allLines[lll])) +
        scale_color_continuous()
      if (saveFigures){
        ggsave(file=paste(currentTitle, "_elongation_vs_temp.png", sep=""), dpi = 600, width = 12, height = 10, units = "cm")
      }
      
      ggplot(currentConditionTopSpeed, aes(x=tTest, y=estimatedlogVolume, fill=factor(tTest))) +
        geom_boxplot(width=0.6, outlier.shape = NA) + # add a box plot
        scale_y_continuous(name=expression(paste('log'[10]*'(volume)'," ",  mu, 'm'^3))) +
        scale_x_continuous(name="tested temperature", limits=c(8, 42)) +
        theme_classic(base_size = 15) +
        theme(legend.position = "none") + 
        ggtitle(paste("TAdapt:", allTAdaptation[aaa], "; C:", allMediumConcentrations[mmm], "; Line: ", allLines[lll])) +
        scale_color_continuous()
      if (saveFigures){
        ggsave(file=paste(currentTitle, "_logVolume_vs_temp.png", sep=""), dpi = 600, width = 12, height = 10, units = "cm")
      } 
      
      
      ggplot(currentCondition, aes(x=tTest, y=medianSpeed, color=medianElongation)) +
        geom_jitter(position = position_jitter(height = 0, width = .5)) + 
        # geom_jitter(position = position_jitter(height = 0, width = .3)) +
        # geom_smooth(method="lm", se=TRUE, fill=NA, formula=y ~ poly(x, 3, raw=TRUE)) + 
        scale_y_continuous(name=expression(paste("median speed (", mu, "m/s)")), limits=c(0,1000)) +
        scale_x_continuous(name="tested temperature", limits=c(8, 42)) +
        theme_classic(base_size = 15) +
        # theme(legend.position = "none") + 
        ggtitle(paste("TAdapt:", allTAdaptation[aaa], "; C:", allMediumConcentrations[mmm], "; Line: ", allLines[lll])) +
        scale_color_continuous(type = "viridis", option="magma", limits = c(1.4, 2.6)) +
        labs(color="Elongation:")
      # scale_color_manual(values=c("#3B9AB2", "#EBCC2A", "#F21A00")) # this is the zissou1 palette
      
      
      if (saveFigures){
        ggsave(file=paste(currentTitle, "_speed_vs_temp_jitter_color_elongation.png", sep=""), dpi = 600, width = 12, height = 10, units = "cm")
      }
      
      
      
      ggplot(currentCondition, aes(x=tTest, y=medianSpeed, color=meanAbsAngle*180/pi)) +
        geom_jitter(position = position_jitter(height = 0, width = .5)) + 
        # geom_jitter(position = position_jitter(height = 0, width = .3)) +
        # geom_smooth(method="lm", se=TRUE, fill=NA, formula=y ~ poly(x, 3, raw=TRUE)) + 
        scale_y_continuous(name=expression(paste("median speed (", mu, "m/s)")), limits=c(0,1000)) +
        scale_x_continuous(name="tested temperature", limits=c(8, 42)) +
        theme_classic(base_size = 15) +
        # theme(legend.position = "none") + 
        ggtitle(paste("TAdapt:", allTAdaptation[aaa], "; C:", allMediumConcentrations[mmm], "; Line: ", allLines[lll])) +
        scale_color_continuous(type = "viridis", option="magma", limits = c(0, 30)) +
        labs(color="abs(turning angle):")
      # scale_color_manual(values=c("#3B9AB2", "#EBCC2A", "#F21A00")) # this is the zissou1 palette
      
      
      if (saveFigures){
        ggsave(file=paste(currentTitle, "_speed_vs_temp_jitter_color_angle.png", sep=""), dpi = 600, width = 12, height = 10, units = "cm")
      }
      
      
      
      
      ggplot(currentCondition, aes(x=tTest, y=medianSpeed, color=autocorrelationTime*medianSpeed)) +
        geom_jitter(position = position_jitter(height = 0, width = .5)) + 
        # geom_jitter(position = position_jitter(height = 0, width = .3)) +
        # geom_smooth(method="lm", se=TRUE, fill=NA, formula=y ~ poly(x, 3, raw=TRUE)) + 
        scale_y_continuous(name=expression(paste("median speed (", mu, "m/s)")), limits=c(0,1000)) +
        scale_x_continuous(name="tested temperature", limits=c(8, 42)) +
        theme_classic(base_size = 18) +
        scale_color_continuous(type = "viridis", option="magma", limits=c(0,300)) +
        labs(color="autocorrelation\ndistance:")
      
      if (saveFigures){
        ggsave(file=paste(currentTitle, "_speed_vs_temp_jitter_color_autocorr_distance.png", sep=""), dpi = 600, width = 12, height = 10, units = "cm")
      }
      
      
      
      ggplot(currentCondition, aes(x=tTest, y=dragPower)) +
        geom_jitter(position = position_jitter(height = 0, width = .5)) + 
        # geom_jitter(position = position_jitter(height = 0, width = .3)) +
        # geom_smooth(method="lm", se=TRUE, fill=NA, formula=y ~ poly(x, 3, raw=TRUE)) + 
        scale_y_continuous(name="drag power (W)") +
        scale_x_continuous(name="tested temperature", limits=c(8, 42)) +
        theme_classic(base_size = 15) +
        theme(legend.position = "none") + 
        ggtitle(paste("TAdapt:", allTAdaptation[aaa], "; C:", allMediumConcentrations[mmm], "; Line: ", allLines[lll])) +
        scale_color_manual(values=c("#3B9AB2", "#EBCC2A", "#F21A00")) # this is the zissou1 palette
      if (saveFigures){
        ggsave(file=paste(currentTitle, "_stokes_power_vs_temp.png", sep=""), dpi = 600, width = 12, height = 10, units = "cm")
      }
      
      # ggplot(currentCondition, aes(x=as.factor(tTest), y=dragPower, shape=as.factor(mediumConcentration), color=as.factor(mediumConcentration))) +
      #   geom_violin() + 
      #   # geom_jitter(position = position_jitter(height = 0, width = .3)) +
      #   # geom_smooth(method="lm", se=TRUE, fill=NA, formula=y ~ poly(x, 3, raw=TRUE)) + 
      #   scale_y_continuous(name="Stokes power (W)") +
      #   scale_x_discrete(name="tested temperature") +
      #   theme_classic(base_size = 15) +
      #   theme(legend.position = "none") + 
      #   ggtitle(paste("TAdapt:", allTAdaptation[aaa], "; C:", allMediumConcentrations[mmm], "; Line: ", allLines[lll])) +
      #   scale_color_manual(values=c("#3B9AB2", "#EBCC2A", "#F21A00")) # this is the zissou1 palette
      # if (saveFigures){
      #   ggsave(file=paste(currentTitle, "_stokes_power_vs_temp_violin.png", sep=""), dpi = 600, width = 12, height = 10, units = "cm")
      # }
      
      # the lines below are some tests with only the activation part
      
      currentConditionUnder30 <- subset(currentConditionTopSpeed, tTest < 30)
      
      
      ggplot(currentConditionUnder30, aes(x=oneOverkTrel, y=lnDragPower, fill=line)) +
        # geom_point() +
        geom_jitter(position = position_jitter(height = 0, width = .5)) +
        geom_smooth(method="lm", se=TRUE, fill="red", formula=y ~ x, na.rm=TRUE) +
        scale_y_continuous(name="ln(drag power (W))") +
        scale_x_continuous(name="1/kT0 - 1/kT") +
        theme_classic(base_size = 15) +
        theme(legend.position = "none") +
        ggtitle(paste("TAdapt:", allTAdaptation[aaa], "; C:", allMediumConcentrations[mmm], "; Line: ", allLines[lll])) +
        scale_color_manual(values=c("#3B9AB2", "#EBCC2A", "#F21A00")) # this is the zissou1 palette
      if (saveFigures){
        ggsave(file=paste(currentTitle, "_arrhenius.png", sep=""), dpi = 600, width = 12, height = 10, units = "cm")
      }
      # 
      # ggplot(currentConditionUnder30, aes(x=medianArea, y=lnDragPower, fill=line)) +
      #   # geom_point() +
      #   geom_jitter(position = position_jitter(height = 0, width = .5)) +
      #   geom_smooth(method="lm", se=TRUE, fill="red", formula=y ~ x, na.rm=TRUE) +
      #   scale_y_continuous(name="ln(drag power (W))") +
      #   scale_x_continuous(name="median area") +
      #   theme_classic(base_size = 15) +
      #   theme(legend.position = "none") +
      #   ggtitle(paste("TAdapt:", allTAdaptation[aaa], "; C:", allMediumConcentrations[mmm], "; Line: ", allLines[lll])) +
      #   scale_color_manual(values=c("#3B9AB2", "#EBCC2A", "#F21A00")) # this is the zissou1 palette
      # 
      
      # This below does not make much sense because the speed depends primarily on
      # temperature
      #
      # ggplot(currentConditionUnder30, aes(x=medianArea, y=medianSpeed, color=factor(line))) +
      #   # geom_point() +
      #   geom_jitter(position = position_jitter(height = 0, width = .5)) +
      #   # geom_smooth(method="lm", se=TRUE, fill="red", formula=y ~ x, na.rm=TRUE) +
      #   scale_y_continuous(name=expression(paste("median speed (", mu, "m/s)"))) +
      #   scale_x_continuous(name="median area") +
      #   theme_classic(base_size = 15) +
      #   theme(legend.position = "none") +
      #   ggtitle(paste("TAdapt:", allTAdaptation[aaa], "; C:", allMediumConcentrations[mmm], "; Line: ", allLines[lll])) +
      #   scale_color_manual(values=c("#FF0000", "#00FF00", "#00FFFF", "#FFFF00"))
      #
      # if (saveFigures){
      #   ggsave(file=paste(currentTitle, "_median_speed_vs_area.png", sep=""), dpi = 600, width = 12, height = 10, units = "cm")
      # }
      #
      # # Change line colors by groups
      # ggplot(currentCondition, aes(x=medianArea, color=factor(line), fill=factor(line))) +
      #   geom_histogram(position="identity", alpha=0.5, binwidth=20) +
      #   scale_color_manual(values=c("#A3A3A3", "#666666", "#000000")) +
      #   scale_fill_manual(values=c("#A3A3A3", "#666666", "#000000")) +
      #   scale_x_continuous() +
      #   labs(title = paste("TAdapt:", allTAdaptation[aaa], "; C:", allMediumConcentrations[mmm], "; Line: ", allLines[lll])) +
      #   theme_classic(base_size=15) +
      #   theme(legend.position = "none")
      # if (saveFigures){
      #   ggsave(file=paste(currentTitle, "_area_histogram.png", sep=""), dpi = 600, width = 12, height = 10, units = "cm")
      # }
      
      
      currentConditionUnder30 <- subset(currentConditionUnder30, is.finite(lnDragPower) & is.finite(oneOverkT))
      lm(lnDragPower ~ oneOverkTrel, data=currentConditionUnder30)
      
      
      
      
      
      
      
      ###### Fit schoolfield to drag power        
      
      # load in data
      #data("chlorella_tpc")
      
      # keep just a single curve
      # d <- filter(chlorella_tpc, curve_id == 1)
      
      if(!require(dplyr)){install.packages('dplyr')}
      d <- currentConditionTopSpeed %>% 
        dplyr::select("line", "tAdaptation", "dragPower", "tTest")%>%
        dplyr::rename(
          curve_id = line,
          growth_temp = tAdaptation,
          rate = dragPower,
          temp = tTest
        )
      
      # convert Watts to eV/s
      # eVperW <- 6241506363094027800
      # d$rate <- d$rate * eVperW
      
      # convert to nanoWatt
      d$rate <- d$rate * 10^9
      
      
      # convert Watts to J/hour
      # d$rate <- d$rate * 3600
      
      d <- d[complete.cases(d), ]
      # I added this because the get_start_vals function does not seem to work with na values
      # and when using w per cell as units the numbers are outside of the default boundaries
      
      
      
      
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
      
      calculatedFitParameters$tAdaptation <- allTAdaptation[aaa]
      calculatedFitParameters$line <- allLines[lll]
      calculatedFitParameters$mediumConcentration <- allMediumConcentrations[mmm]
      
      
      
      if (is.null(fit))
      {
        calculatedFitParameters <- cbind(calculatedFitParameters, data.frame(r_tref=NA, e_from_fit=NA, eh_from_fit=NA, th=NA, r_tadapt = NA))
      } else{
        # rate at the adaptation temperature:
        r_tadapt <- (augment(fit, newdata=data.frame(temp=as.numeric(allTAdaptation[aaa]))))$.fitted
        calculatedFitParameters <- cbind(calculatedFitParameters, data.frame(r_tref=coef(fit)[1], e_from_fit=coef(fit)[2], eh_from_fit=coef(fit)[3], th=coef(fit)[4]), r_tadapt=r_tadapt)
      }
      
      print(paste("TAdapt:", allTAdaptation[aaa], "; C:", allMediumConcentrations[mmm], "; Line: ", allLines[lll]))
      print(calculatedFitParameters)
      
      print(fit)
      
      if (iteration == 1)
      {
        allFitResults <- calculatedFitParameters[0,]
      }
      
      allFitResults <- rbind(allFitResults, calculatedFitParameters)
      
      # predict new data
      new_data <- data.frame(temp = seq(min(d$temp), max(d$temp), 0.5))
      preds <- augment(fit, newdata = new_data)
      
      
      lineColours = c("#3B9AB2", "#EBCC2A", "#F21A00", "#FF00FF")
      markerColours = c("#A3A3A3", "#666666", "#000000", "#FF00FF")
      markerShapes = c(16,17,15,4)
      
      
      # plot data and model fit
      thisPlot <- ggplot(d, aes(temp, rate)) +
        # geom_point(col=markerColours[mmm], shape=markerShapes[mmm]) +
        geom_jitter(position = position_jitter(height = 0, width = .5), col=markerColours[mmm], shape=markerShapes[mmm]) +
        annotate("text", size=5, x=10.5, y=3e-04, label= paste("E=", round(calculatedFitParameters$e_from_fit, 2), "eV", sep=""), hjust = 0, parse=F) +
        geom_line(aes(temp, .fitted), preds, col = lineColours[aaa], size=2) +
        theme_classic(base_size = 15) +
        scale_x_continuous(name="Temp. (°C)",  limits=c(7.5, 42.5)) +
        scale_y_continuous(name="Stokes power (nW)") +
        # labs(x = 'Temp (ºC)',  y = 'Stokes power (nW)') +
        ggtitle(paste("TAdapt:", allTAdaptation[aaa], "; C:", allMediumConcentrations[mmm], "; Line: ", allLines[lll]))
      
      thisPlot
      
      if (saveFigures){
        ggsave(file=paste(currentTitle, "_schoolfield.png", sep=""), dpi = 600, width = 12, height = 10, units = "cm")
      }
      
      
      
      
      ## Fit schoolfield to the speed directly
      
      
      
      
      if(!require(dplyr)){install.packages('dplyr')}
      d <- currentConditionTopSpeed %>% 
        dplyr::select("line", "tAdaptation", "medianSpeed", "tTest")%>%
        dplyr::rename(
          curve_id = line,
          growth_temp = tAdaptation,
          rate = medianSpeed,
          temp = tTest
        )
      
      
      
      # convert Watts to J/hour
      # d$rate <- d$rate * 3600
      
      d <- d[complete.cases(d), ]
      # I added this because the get_start_vals function does not seem to work with na values
      # and when using w per cell as units the numbers are outside of the default boundaries
      
      
      
      
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
      
      calculatedFitParameters$tAdaptation <- allTAdaptation[aaa]
      calculatedFitParameters$line <- allLines[lll]
      calculatedFitParameters$mediumConcentration <- allMediumConcentrations[mmm]
      
      if (is.null(fit))
      {
        calculatedFitParameters <- cbind(calculatedFitParameters, data.frame(r_tref=NA, e_from_fit=NA, eh_from_fit=NA, th=NA))
      } else{
        r_tadapt <- (augment(fit, newdata=data.frame(temp=as.numeric(allTAdaptation[aaa]))))$.fitted
        calculatedFitParameters <- cbind(calculatedFitParameters, data.frame(r_tref=coef(fit)[1], e_from_fit=coef(fit)[2], eh_from_fit=coef(fit)[3], th=coef(fit)[4]), r_tadapt = r_tadapt)
      }
      
      print(paste("TAdapt:", allTAdaptation[aaa], "; C:", allMediumConcentrations[mmm], "; Line: ", allLines[lll]))
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
          do(data.frame(temp = seq(10, 40, 0.5))) %>% 
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
        allFitResultsOnSpeed <- calculatedFitParameters[0,]
      }
      
      allFitResultsOnSpeed <- rbind(allFitResultsOnSpeed, calculatedFitParameters)
      
      # predict new data
      new_data <- data.frame(temp = seq(min(d$temp), max(d$temp), 0.5))
      preds <- augment(fit, newdata = new_data)
      
      
      
      lineColours = c("#3B9AB2", "#EBCC2A", "#F21A00", "#FF00FF")
      markerColours = c("#A3A3A3", "#666666", "#000000", "#FF00FF")
      markerShapes = c(16,17,15,4)
      
      # plot data and model fit
      thisPlot <- ggplot(d, aes(temp, rate)) +
        # geom_point(col=markerColours[mmm], shape=markerShapes[mmm]) +
        geom_jitter(position = position_jitter(height = 0, width = .5), col=markerColours[mmm], shape=markerShapes[mmm]) +
        annotate("text", size=5, x=10.5, y=1020, label= paste("E=", round(calculatedFitParameters$e_from_fit, 2), "eV", sep=""), hjust = 0, parse=F) +
        geom_line(aes(temp, .fitted), preds, col = lineColours[aaa], size=2) +
        theme_classic(base_size = 15) +
        scale_x_continuous(name="Temp. (°C)",  limits=c(7.5, 42.5)) +
        scale_y_continuous(name=expression(paste("Speed (", mu, "m/s", ")")), limits=c(0,1100), breaks=seq(0, 1000, by=200)) +
        # labs(x = 'Temperature (ºC)', y = 'Speed (um/s)')+
        ggtitle(paste("T:", allTAdaptation[aaa], "; C:", allMediumConcentrations[mmm], "; L: ", allLines[lll]))
      thisPlot
      if (saveFigures){
        ggsave(file=paste(currentTitle, "_Speed_schoolfield.png", sep=""), dpi = 600, width = 12, height = 10, units = "cm")
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
        annotate("text", size=5, x=10.5, y=1020, label= paste("E=", round(calculatedFitParameters$e_from_fit, 2), "eV", sep=""), hjust = 0, parse=F) +
        geom_errorbar(aes(ymin=rate-sd, ymax=rate+sd), width=1.3, col=markerColours[mmm]) + 
        theme_classic(base_size = 14) +
        scale_x_continuous(name="Temp. (°C)",  limits=c(7.5, 42.5)) +
        scale_y_continuous(name=expression(paste("Speed (", mu, "m/s", ")")), limits=c(0,1100), breaks=seq(0, 1000, by=200))
      # ggtitle(paste("T:", allTAdaptation[aaa], "; C:", allMediumConcentrations[mmm], "; L: ", allLines[lll]))
      
      thisPlot
      if (saveFigures){
        ggsave(file=paste(currentTitle, "_Speed_schoolfield_errorbar.png", sep=""), dpi = 600, width = 12, height = 10, units = "cm")
      }
      
      
      if (lll == 1 & !is.finite(as.numeric(allLines[lll]))){
        plotList[[plotCounter]] <- thisPlot
      } else if (is.finite(as.numeric(allLines[1])) & lll==length(allLines)) {
        plotList[[plotCounter]] <- thisPlot
      } else {
        plotListExtra[[plotCounterExtra]] <- thisPlot
      }
      
      
      if (includeBootstrap)
      {
        # plot bootstrapped CIs
        ggplot(d2, aes(temp, rate)) +
          geom_line(aes(temp, .fitted), preds, col = lineColours[aaa], size=2) +
          geom_line(aes(temp, pred, group = iter), boot1_preds, col = lineColours[aaa], alpha = 0.007) +
          geom_point(size=3, col=markerColours[mmm], shape=markerShapes[mmm]) +
          annotate("text", size=5, x=10.5, y=1020, label= paste("E=", round(calculatedFitParameters$e_from_fit, 2), "eV", sep=""), hjust = 0, parse=F) +
          geom_errorbar(aes(ymin=rate-sd, ymax=rate+sd), width=1.3, col=markerColours[mmm]) + 
          theme_classic(base_size = 14) +
          scale_x_continuous(name="Temp. (°C)",  limits=c(7.5, 42.5)) +
          scale_y_continuous(name=expression(paste("Speed (", mu, "m/s", ")")), limits=c(0,1100), breaks=seq(0, 1000, by=200))
        if (saveFigures){
          ggsave(file=paste(currentTitle, "_Speed_schoolfield_errorbar_and_bootstrap.png", sep=""), dpi = 600, width = 12, height = 10, units = "cm")
        }
        
        # plot bootstrapped CIs
        ggplot() +
          geom_line(aes(temp, .fitted), preds, col = lineColours[aaa], size=2) +
          geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot1_conf_preds, fill = lineColours[aaa], alpha = 0.3) +
          geom_point(data=d2, aes(temp, rate), size=3, col=markerColours[mmm], shape=markerShapes[mmm]) +
          annotate("text", size=5, x=10.5, y=1020, label= paste("E=", round(calculatedFitParameters$e_from_fit, 2), "eV", sep=""), hjust = 0, parse=F) +
          geom_errorbar(data=d2, aes(x=temp, ymin=rate-sd, ymax=rate+sd), width=1.3, col=markerColours[mmm]) + 
          theme_classic(base_size = 14) +
          scale_x_continuous(name="Temp. (°C)",  limits=c(7.5, 42.5)) +
          scale_y_continuous(name=expression(paste("Speed (", mu, "m/s", ")")), limits=c(0,1100), breaks=seq(0, 1000, by=200))
        if (saveFigures){
          ggsave(file=paste(currentTitle, "_Speed_schoolfield_errorbar_and_bootstrap2.png", sep=""), dpi = 600, width = 12, height = 10, units = "cm")
        }
        
        
      }
      
      
      
      
      
    }
  }
}


print(allFitResults)

if (combineLinesTogether)
{
  write.table(allFitResults, "allFitResults_Schoolfield_on_Stokes_Power_acute_response_all.csv", append = FALSE, sep = ", ", row.names = FALSE)
  write.table(allFitResultsOnSpeed, "allFitResults_Schoolfield_on_Speed_acute_response_all.csv", append = FALSE, sep = ", ", row.names = FALSE)
} else {
  write.table(allFitResults, "allFitResults_Schoolfield_on_Stokes_Power_acute_response.csv", append = FALSE, sep = ", ", row.names = FALSE)
  write.table(allFitResultsOnSpeed, "allFitResults_Schoolfield_on_Speed_acute_response.csv", append = FALSE, sep = ", ", row.names = FALSE)
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
library(Cairo)

if (saveFigures){
  if (combineLinesTogether)
  {
    ggsave(file="figure_speed_thermal_response_post_adaptation_all.png", dpi = 600, width = 24, height = 20, units = "cm")
    ggsave(file="figure_speed_thermal_response_post_adaptation_all.eps", device="eps", dpi = 1200, width = 24, height = 20, units = "cm")
    ggsave(file="figure_speed_thermal_response_post_adaptation_all.pdf", device=cairo_pdf, dpi = 1200, width = 24, height = 20, units = "cm")
  }
  else
  {
    ggsave(file="figure_speed_thermal_response_post_adaptation.png", dpi = 600, width = 24, height = 20, units = "cm")
    ggsave(file="figure_speed_thermal_response_post_adaptation.eps", device="eps", dpi = 1200, width = 24, height = 20, units = "cm")
    ggsave(file="figure_speed_thermal_response_post_adaptation.pdf", device=cairo_pdf, dpi = 1200, width = 24, height = 20, units = "cm")
    
  }
}




# If I want to exclude the mother culture:
if (skipMotherCulture == TRUE)
{
  allFitResults <- subset(allFitResults, is.finite(as.numeric(allFitResults$line)) == FALSE)
  allFitResultsOnSpeed <- subset(allFitResultsOnSpeed, is.finite(as.numeric(allFitResultsOnSpeed$line)) == FALSE)
}

# plot data and model fit
ggplot(allFitResults, aes(x=tAdaptation, y=e_from_fit, color=factor(mediumConcentration), shape=factor(mediumConcentration), label=factor(line))) +
  # geom_point(size=6) +
  geom_jitter(size=6, position = position_jitter(height = 0, width = .3)) + 
  # geom_text(color="black", size=4) + 
  scale_color_manual(values=c("#A3A3A3", "#666666", "#000000")) +
  scale_fill_manual(values=c("#A3A3A3", "#666666", "#000000")) +
  scale_y_continuous(name="activation energy (eV)") +
  scale_x_continuous(name="Adaptation temperature (ºC)", limits=c(12.5, 27.5), breaks=c(15, 20, 25)) +
  theme_classic(base_size=15) +
  #theme(legend.position = "none") +
  labs(color = "Conc.") + # this specifies a custom legend
  ggtitle("Activation energy (this includes estimation of shape\nperhaps not very reliable in these data)")

if (saveFigures){
  if (combineLinesTogether){
    ggsave(file="aaaa_activation_energy_all.png", dpi = 600, width = 12, height = 10, units = "cm")
  } else {
    ggsave(file="aaaa_activation_energy.png", dpi = 600, width = 12, height = 10, units = "cm")
    
  }
}




# plot data and model fit
ggplot(allFitResults, aes(x=tAdaptation, y=topt, color=factor(mediumConcentration), shape=factor(mediumConcentration), label=factor(line))) +
  # geom_point(size=6) +
  geom_jitter(size=6, position = position_jitter(height = 0, width = .3)) + 
  # geom_text(color="black", size=4) + 
  scale_color_manual(values=c("#A3A3A3", "#666666", "#000000")) +
  scale_fill_manual(values=c("#A3A3A3", "#666666", "#000000")) +
  scale_y_continuous(name="Optimal Temperature (ºC)") +
  scale_x_continuous(name="Adaptation temperature (ºC)", limits=c(12.5, 27.5), breaks=c(15, 20, 25)) +
  theme_classic(base_size=15) +
  theme(legend.position = "none") +
  ggtitle("Optimal temperature on Stokes power")

if (saveFigures){
  if (combineLinesTogether){
    ggsave(file="aaaa_optimal_temp_vs_temp_all.png", dpi = 600, width = 12, height = 10, units = "cm")
  } else {
    ggsave(file="aaaa_optimal_temp_vs_temp.png", dpi = 600, width = 12, height = 10, units = "cm")
    
  }
}








# plot data and model fit
ggplot(allFitResultsOnSpeed, aes(x=tAdaptation, y=e_from_fit, color=factor(mediumConcentration), fill=factor(mediumConcentration + tAdaptation*100), shape=factor(mediumConcentration), label=factor(line))) +
  # geom_point(size=4) +
  geom_jitter(size=4, position = position_jitter(height = 0, width = .3)) + 
  # geom_text(color="black", size=4) + 
  scale_color_manual(values=c("#A3A3A3", "#666666", "#000000")) +
  scale_fill_manual(values=c("#8FCDDC", "#3B9AB2", "#18434E", "#F2DD70", "#EBCC2A", "#463C07", "#FF7C6C", "#F21A00", "#660B00")) +
  scale_shape_manual(values=c(21, 24, 22)) + # shapes for the markers
  scale_y_continuous(name="Activation Energy (eV)", limits=c(0, 1)) +
  scale_x_continuous(name="Adaptation temperature (ºC)", limits=c(12.5, 27.5), breaks=c(15, 20, 25)) +
  theme_classic(base_size=18) +
  theme(legend.position = "none") +
  labs(color = "Conc. (%)") # + # this specifies a custom legend
# ggtitle("Activation energy on speed")

if (saveFigures){
  if (combineLinesTogether){
    ggsave(file="figure_activation_energy_on_speed_all.png", dpi = 600, width = 12, height = 10, units = "cm")
    ggsave(file="figure_activation_energy_on_speed_all.eps", device="eps", dpi = 1200, width = 12, height = 10, units = "cm")
    library(Cairo)
    ggsave(file="figure_activation_energy_on_speed_all.pdf", device=cairo_pdf, dpi = 1200, width = 12, height = 10, units = "cm")
  } else {
    ggsave(file="figure_activation_energy_on_speed.png", dpi = 600, width = 12, height = 10, units = "cm")
    ggsave(file="figure_activation_energy_on_speed.eps", device="eps", dpi = 1200, width = 12, height = 10, units = "cm")
    library(Cairo)
    ggsave(file="figure_activation_energy_on_speed.pdf", device=cairo_pdf, dpi = 1200, width = 12, height = 10, units = "cm")
    
  }
}




# plot data and model fit
ggplot(allFitResultsOnSpeed, aes(x=tAdaptation, y=topt, color=factor(mediumConcentration), fill=factor(mediumConcentration + tAdaptation*100), shape=factor(mediumConcentration), label=factor(line))) +
  # geom_point(size=4) +
  geom_jitter(size=4, position = position_jitter(height = 0, width = .3)) + 
  # geom_text(color="black", size=4) + 
  scale_color_manual(values=c("#A3A3A3", "#666666", "#000000")) +
  scale_fill_manual(values=c("#8FCDDC", "#3B9AB2", "#18434E", "#F2DD70", "#EBCC2A", "#463C07", "#FF7C6C", "#F21A00", "#660B00")) +
  scale_shape_manual(values=c(21, 24, 22)) + # shapes for the markers
  scale_y_continuous(name="Optimal temperature (ºC)") +
  scale_x_continuous(name="Adaptation temperature (ºC)", limits=c(12.5, 27.5), breaks=c(15, 20, 25)) +
  theme_classic(base_size=12) +
  theme(legend.position = "none") +
  labs(color = "Conc.") # + # this specifies a custom legend
# ggtitle("Optimal temperature")

if (saveFigures){
  if (combineLinesTogether){
    ggsave(file="figure_optimal_temp_vs_temp_on_speed_all.png", dpi = 600, width = 12, height = 10, units = "cm")
    ggsave(file="figure_optimal_temp_vs_temp_on_speed_all.eps", device="eps", dpi = 1200, width = 12, height = 10, units = "cm")
    library(Cairo)
    ggsave(file="figure_optimal_temp_vs_temp_on_speed_all.pdf", device=cairo_pdf, dpi = 1200, width = 12, height = 10, units = "cm")
  } else {
    ggsave(file="figure_optimal_temp_vs_temp_on_speed.png", dpi = 600, width = 12, height = 10, units = "cm")
    ggsave(file="figure_optimal_temp_vs_temp_on_speed.eps", device="eps", dpi = 1200, width = 12, height = 10, units = "cm")
    library(Cairo)
    ggsave(file="figure_optimal_temp_vs_temp_on_speed.pdf", device=cairo_pdf, dpi = 1200, width = 12, height = 10, units = "cm")
  }
}





# plot data and model fit
ggplot(allFitResultsOnSpeed, aes(x=tAdaptation, y=r_tref, color=factor(mediumConcentration), fill=factor(mediumConcentration + tAdaptation*100), shape=factor(mediumConcentration), label=factor(line))) +
  # geom_point(size=4) +
  geom_jitter(size=4, position = position_jitter(height = 0, width = .3)) + 
  # geom_text(color="black", size=4) + 
  scale_color_manual(values=c("#A3A3A3", "#666666", "#000000")) +
  scale_fill_manual(values=c("#8FCDDC", "#3B9AB2", "#18434E", "#F2DD70", "#EBCC2A", "#463C07", "#FF7C6C", "#F21A00", "#660B00")) +
  scale_shape_manual(values=c(21, 24, 22)) + # shapes for the markers
  scale_y_continuous(name=expression(paste("Speed (", mu, "m/s", ")")), limits=c(0, 900)) +
  # scale_y_continuous(name=expression(paste("Speed (", mu, "m/s", ")")), sec.axis = sec_axis( trans=~.*1e-6, name=paste("Speed (m/s)"))) +
  scale_x_continuous(name=expression("Adaptation temperature (ºC)"), limits=c(12.5, 27.5), breaks=c(15, 20, 25)) +
  theme_classic(base_size=12) +
  theme(legend.position = "none") +
  labs(color = "Conc.") + # this specifies a custom legend
  ggtitle(paste("Speed at reference temperature (", referenceTemperature, "°C)", sep=""))

if (saveFigures)
{
  if (combineLinesTogether)
  {
    ggsave(file="figure_speed_at_tref_all.png", dpi = 600, width = 12, height = 10, units = "cm")
    ggsave(file="figure_speed_at_tref_all.eps", device="eps", dpi = 1200, width = 12, height = 10, units = "cm")
    library(Cairo)
    ggsave(file="figure_speed_at_tref_all.pdf", device=cairo_pdf, dpi = 1200, width = 12, height = 10, units = "cm")
  } else {
    ggsave(file="figure_speed_at_tref.png", dpi = 600, width = 12, height = 10, units = "cm")
    ggsave(file="figure_speed_at_tref.eps", device="eps", dpi = 1200, width = 12, height = 10, units = "cm")
    # quartz(type = 'pdf', file="aaaa_metabolic_thermal_response_post_adaptation.pdf")
    library(Cairo)
    ggsave(file="figure_speed_at_tref.pdf", device=cairo_pdf, dpi = 1200, width = 12, height = 10, units = "cm")
    
  }
  
}





# plot data and model fit
ggplot(allFitResultsOnSpeed, aes(x=tAdaptation, y=r_tadapt, color=factor(mediumConcentration), fill=factor(mediumConcentration + tAdaptation*100), shape=factor(mediumConcentration), label=factor(line))) +
  # geom_point(size=4) +
  geom_jitter(size=4, position = position_jitter(height = 0, width = .3)) + 
  # geom_text(color="black", size=4) + 
  scale_color_manual(values=c("#A3A3A3", "#666666", "#000000")) +
  scale_fill_manual(values=c("#8FCDDC", "#3B9AB2", "#18434E", "#F2DD70", "#EBCC2A", "#463C07", "#FF7C6C", "#F21A00", "#660B00")) +
  scale_shape_manual(values=c(21, 24, 22)) + # shapes for the markers
  scale_y_continuous(name=expression(paste("Speed (", mu, "m/s", ")")), limits=c(0, 900)) +
  # scale_y_continuous(name=expression(paste("Speed (", mu, "m/s", ")")), sec.axis = sec_axis( trans=~.*1e-6, name=paste("Speed (m/s)"))) +
  scale_x_continuous(name=expression("Adaptation temperature (ºC)"), limits=c(12.5, 27.5), breaks=c(15, 20, 25)) +
  theme_classic(base_size=12) +
  theme(legend.position = "none") +
  labs(color = "Conc.") + # this specifies a custom legend
  ggtitle("Speed at adaptation temperature")

if (saveFigures)
{
  if (combineLinesTogether) {
    ggsave(file="aaaa_speed_at_tadapt_all.png", dpi = 600, width = 12, height = 10, units = "cm")
  } else{
    ggsave(file="aaaa_speed_at_tadapt.png", dpi = 600, width = 12, height = 10, units = "cm")
    
  }
}




if (includeBootstrap)
{
  
  
  # plot data and model fit
  ggplot(allFitResultsOnSpeed, aes(x=tAdaptation, y=e_from_fit, color=factor(mediumConcentration), fill=factor(mediumConcentration + tAdaptation*100), shape=factor(mediumConcentration), label=factor(line))) +
    geom_point(size=6, position=position_dodge(width=2)) + 
    # geom_errorbar(aes(ymin=e_from_fit_0025, ymax=e_from_fit_0975), width=0, position=position_dodge(width=2)) + 
    geom_errorbar(aes(ymin=e_from_fit - sd_e_from_fit, ymax=e_from_fit + sd_e_from_fit), width=0, position=position_dodge(width=2)) + 
    scale_color_manual(values=c("#A3A3A3", "#666666", "#000000")) +
    scale_fill_manual(values=c("#8FCDDC", "#3B9AB2", "#18434E", "#F2DD70", "#EBCC2A", "#463C07", "#FF7C6C", "#F21A00", "#660B00")) +
    scale_shape_manual(values=c(21, 24, 22)) + # shapes for the markers
    scale_y_continuous(name="Activation Energy (eV)", limits=c(0, 1)) +
    scale_x_continuous(name="Adaptation temperature (ºC)", limits=c(12.5, 27.5), breaks=c(15, 20, 25)) +
    theme_classic(base_size=18) +
    theme(legend.position = "none") +
    labs(color = "Conc. (%)") # + # this specifies a custom legend
  # ggtitle("Activation energy on speed")
  
  if (saveFigures){
    if (combineLinesTogether){
      ggsave(file="figure_activation_energy_on_speed_bootstrap_all.png", dpi = 600, width = 12, height = 10, units = "cm")
      ggsave(file="figure_activation_energy_on_speed_bootstrap_all.eps", device="eps", dpi = 1200, width = 12, height = 10, units = "cm")
      library(Cairo)
      ggsave(file="figure_activation_energy_on_speed_bootstrap_all.pdf", device=cairo_pdf, dpi = 1200, width = 12, height = 10, units = "cm")
    }else{
      ggsave(file="figure_activation_energy_on_speed_bootstrap.png", dpi = 600, width = 12, height = 10, units = "cm")
      ggsave(file="figure_activation_energy_on_speed_bootstrap.eps", device="eps", dpi = 1200, width = 12, height = 10, units = "cm")
      library(Cairo)
      ggsave(file="figure_activation_energy_on_speed_bootstrap.pdf", device=cairo_pdf, dpi = 1200, width = 12, height = 10, units = "cm")
    }
  }
  
  
  
  
  # plot data and model fit
  ggplot(allFitResultsOnSpeed, aes(x=tAdaptation, y=r_tref, color=factor(mediumConcentration), fill=factor(mediumConcentration + tAdaptation*100), shape=factor(mediumConcentration), label=factor(line))) +
    geom_point(size=6, position=position_dodge(width=2)) + 
    # geom_errorbar(aes(ymin=r_from_fit_0025, ymax=r_from_fit_0975), width=0, position=position_dodge(width=2)) + 
    geom_errorbar(aes(ymin=r_tref - sd_r_tref, ymax=r_tref + sd_r_tref), width=0, position=position_dodge(width=2)) + 
    scale_color_manual(values=c("#A3A3A3", "#666666", "#000000")) +
    scale_fill_manual(values=c("#8FCDDC", "#3B9AB2", "#18434E", "#F2DD70", "#EBCC2A", "#463C07", "#FF7C6C", "#F21A00", "#660B00")) +
    scale_shape_manual(values=c(21, 24, 22)) + # shapes for the markers
    scale_y_continuous(name=expression(paste("Speed at T=20°C (", mu, "m/s", ")")), limits=c(0, 900)) +
    # scale_y_continuous(name=expression(paste("Speed (", mu, "m/s", ")")), sec.axis = sec_axis( trans=~.*1e-6, name=paste("Speed (m/s)"))) +
    scale_x_continuous(name=expression("Adaptation temperature (ºC)"), limits=c(12.5, 27.5), breaks=c(15, 20, 25)) +
    theme_classic(base_size=18) +
    theme(legend.position = "none") +
    labs(color = "Conc.") # + # this specifies a custom legend
    # ggtitle(paste("Speed at reference temperature (", referenceTemperature, "°C)", sep=""))
  
  if (saveFigures)
  {
    if (combineLinesTogether)
    {
      ggsave(file="figure_speed_at_tref_bootstrap_all.png", dpi = 600, width = 12, height = 10, units = "cm")
      ggsave(file="figure_speed_at_tref_bootstrap_all.eps", device="eps", dpi = 1200, width = 12, height = 10, units = "cm")
      # quartz(type = 'pdf', file="aaaa_metabolic_thermal_response_post_adaptation.pdf")
      library(Cairo)
      ggsave(file="figure_speed_at_tref_bootstrap_all.pdf", device=cairo_pdf, dpi = 1200, width = 12, height = 10, units = "cm")
    }else{
      ggsave(file="figure_speed_at_tref_bootstrap.png", dpi = 600, width = 12, height = 10, units = "cm")
      ggsave(file="figure_speed_at_tref_bootstrap.eps", device="eps", dpi = 1200, width = 12, height = 10, units = "cm")
      # quartz(type = 'pdf', file="aaaa_metabolic_thermal_response_post_adaptation.pdf")
      library(Cairo)
      ggsave(file="figure_speed_at_tref_bootstrap.pdf", device=cairo_pdf, dpi = 1200, width = 12, height = 10, units = "cm")
      
    }
  }
  
  
}





lineColours = c("#3B9AB2", "#EBCC2A", "#F21A00", "#FF00FF")
lineTypes <- c("twodash", "dashed", "solid", "longdash")
# Now try to plot the fit only
allFitData <- data.frame()
for (iii in 1:dim(allFitResultsOnSpeed)[1])
{
  predictedData <- data.frame(temp = seq(10, 40, 0.5))
  predictedData$rate <- sharpeschoolhigh_1981(temp = predictedData$temp, r_tref = allFitResultsOnSpeed$r_tref[iii],e = allFitResultsOnSpeed$e_from_fit[iii], eh = allFitResultsOnSpeed$eh_from_fit[iii], th = allFitResultsOnSpeed$th[iii], tref = 20)
  predictedData$tAdaptation <- allFitResultsOnSpeed$tAdaptation[iii]
  predictedData$mediumConcentration <- allFitResultsOnSpeed$mediumConcentration[iii]
  predictedData$line <- allFitResultsOnSpeed$line[iii]
  allFitData <- rbind(allFitData, predictedData)
}

plotG1 <- ggplot(allFitData, aes(x=temp, y=rate, color=as.factor(tAdaptation), size=as.factor(mediumConcentration), alpha=mediumConcentration)) + # linetype=as.factor(mediumConcentration), 
  geom_line() +
  scale_x_continuous(name="Temp. (°C)",  limits=c(10, 25)) +
  scale_y_continuous(name=expression(paste("Speed (", mu, "m/s", ")")), limits=c(100, 1000), trans = 'log10')  +
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
    ggsave(file="figure_fitted_speed_vs_temp_all.png", dpi = 600, width = 12, height = 10, units = "cm")
    ggsave(file="figure_fitted_speed_vs_temp_all.eps", device="eps", dpi = 1200, width = 12, height = 10, units = "cm")
    library(Cairo)
    ggsave(file="figure_fitted_speed_vs_temp_all.pdf", device=cairo_pdf, dpi = 1200, width = 12, height = 10, units = "cm")
  } else {
    ggsave(file="figure_fitted_speed_vs_temp.png", dpi = 600, width = 12, height = 10, units = "cm")
  }
}


# Look at the data of slope from the fitted models
a <-unlist(allFittedLogSpeedVsLogVolume, recursive=FALSE)
aOLS<- a[seq(1,length(a), by=4)]# this is the ordinary LS
median(aOLS)
aMA<- a[seq(2,length(a), by=4)]# this is the MA
median(aMA)
aSMA<- a[seq(3,length(a), by=4)]# this is the SMA
median(aSMA)
aRMA<- a[seq(4,length(a), by=4)]# this is the RMA
median(aRMA)




