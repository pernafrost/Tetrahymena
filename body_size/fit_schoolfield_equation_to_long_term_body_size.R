rm(list=ls()) # clean memory
if(!is.null(dev.list())) dev.off()





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
includeStartingData <- TRUE # In some conditions, we did not measure the body size of 
# cultures tested at the same temperature at which they were adapted. We could
# complement these data from other measurements

# dev.off()

# load packages
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)
library(ggplot2)

referenceTemperature <- 20 # careful that it is not used in the schoolfield fitting
skipMotherCulture <- TRUE # whether to skip the mother culture from this analysis


fileName <- "~/Tetrahymena/longer_term_speed_response/track_analysis_individual_particles_long_term_response.csv"
if (includeStartingData == TRUE)
{
  fileNameStartingData <- "Tetrahymena/acute_speed_response/track_analysis_results_individual_particles.csv" # file.choose() # ask the user to select a file name.
  allExperimentResultsStartingData <- read.table(file = fileNameStartingData, sep = ",", header=TRUE, na.strings = c("NA", " NA"))
  
  library("stringr")
  # Extract information from the string in the "fileName" field
  # newColumns <- str_split_fixed(as.character(str_trim(tolower(allExperimentResults$fileName))), "/adapted|_tested|_group|diluted", 5)[,2:4]
  
  newColumns <- str_split_fixed(as.character(str_trim(tolower(basename(allExperimentResultsStartingData$fileName)))), "_|_tested_", n=Inf)[,c(2,3,4,6)]
  
  # remove c for degrees celsius and remove bis from the line
  newColumns[,4] <- str_replace(newColumns[,4], "c", "")
  newColumns[,3] <- str_replace(newColumns[,3], "bis", "")
  newColumns[,4] <- str_replace(newColumns[,4], "bis", "")
  
  # class(newColumns) <- "numeric"
  newColumns <- data.frame(tAdaptation = as.numeric(newColumns[,1]), 
                           mediumConcentration = as.numeric(newColumns[,2]),
                           tTest = as.numeric(newColumns[,4]), 
                           line = newColumns[,3],
                           incubationDuration = rep(0,length(newColumns[,1])),
                           repetition = rep(0,length(newColumns[,1])))
  # colnames(newColumns) <- c("tAdaptation", "tTest", "line")
  
  # add the new columns to the original table
  allExperimentResultsStartingData <- cbind(allExperimentResultsStartingData, newColumns, deparse.level = 1, stringsAsFactors = default.stringsAsFactors())
  # sapply(allExperimentResults, class)
  
  # check if the data are for the mother culture and remove them
  allExperimentResultsStartingData <- subset(allExperimentResultsStartingData, line != 1 & line != 2)
  # also exclude the data in which the tTest is very different from the tAdaptation (even if in
  # this dataset the tTest is for an acute response
  allExperimentResultsStartingData <- subset(allExperimentResultsStartingData, abs(tTest - tAdaptation) <=5 ) # tTest not too different from adaptation temperature
  
  # given that tTest here is just for a very short-term exposure, while in all
  # the other data tTest is for at least a few days, I re-assign the values of 
  # tTest before merging the datasets
  allExperimentResultsStartingData$tTest <- allExperimentResultsStartingData$tAdaptation
  
  
  names(allExperimentResultsStartingData)
}

workingDir <- "~/Tetrahymena/body_size"
setwd(workingDir)

# read the data for long-term exposure to temperature
allExperimentResults <- read.table(file = fileName, sep = ",", header=TRUE, na.strings = c("NA", " NA"))


# filter some values with high density
allExperimentResults <- subset(allExperimentResults, meanDensity <= 200000)

library("stringr")
# Extract information from the string in the "fileName" field
newColumns <- str_split_fixed(as.character(str_trim(tolower(basename(allExperimentResults$fileName)))), "_|_growth_|hours", n=Inf)[,c(2,3,4,6,7,9)]

# remove c for degrees celsius and remove bis from the line
newColumns[,4] <- str_replace(newColumns[,4], "c", "")
newColumns[,3] <- str_replace(newColumns[,3], "bis", "")

# class(newColumns) <- "numeric"
newColumns <- data.frame(tAdaptation = as.numeric(newColumns[,1]), 
                         mediumConcentration = as.numeric(newColumns[,2]),
                         tTest = as.numeric(newColumns[,4]), 
                         line = newColumns[,3],
                         incubationDuration = as.numeric(newColumns[,5]),
                         repetition = newColumns[,6])


# add the new columns to the original table
allExperimentResults <- cbind(allExperimentResults, newColumns, deparse.level = 1, stringsAsFactors = default.stringsAsFactors())
# sapply(allExperimentResults, class)

# merge the starting data with the new data
if (includeStartingData == TRUE)
{
  allExperimentResults <- rbind(allExperimentResults, allExperimentResultsStartingData)
}


allExperimentResults$logSpeed <- log10(allExperimentResults$medianSpeed)


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

allExperimentResults$estimatedVolume <- (4/3 * pi* allExperimentResults$medianBEllipse * allExperimentResults$medianAEllipse * allExperimentResults$medianAEllipse / 8)
allExperimentResults$estimatedLogVolume <- log10(allExperimentResults$estimatedVolume)
allExperimentResults$oneOverVolume <- (4/3 * pi* allExperimentResults$medianBEllipse * allExperimentResults$medianAEllipse * allExperimentResults$medianAEllipse / 8)^-1


names(allExperimentResults)

allTAdaptation <- sort(unique(allExperimentResults$tAdaptation))
plotList = list() # To save the plots in a list
plotList2 = list() # there are four plots for each condition and I save them all
plotList3 = list()
plotList4 = list()
plotCounter1 <- 1
plotCounter2 <- 1
plotCounter3 <- 1
plotCounter4 <- 1
plotListOverTime = list()
plotCounterOverTime <- 1
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
      currentTitle <- paste("bbbbb_", allTAdaptation[aaa], "_", allMediumConcentrations[mmm], "_", allLines[lll], sep="")
      
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
      
      # If I further filter the data:
      # keep those that are big enough and that move
      currentCondition <- subset(currentCondition, medianArea > 400 & medianSpeed > 100)
      
      # If also filtering by time (less than two days)
      # currentCondition <- subset(currentCondition,incubationDuration <= 48) 
      # currentTitle <- paste("bbbbb_less_than_48_", allTAdaptation[aaa], "_", allMediumConcentrations[mmm], "_", allLines[lll], sep="")
      
      
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
      currentConditionTopSpeed <- subset(currentCondition, medianSpeed >= speedQuantile)
      
      
      
      ggplot(currentCondition, aes(x=tTest, y=medianArea, color=medianSpeed)) +
        geom_jitter(position = position_jitter(height = 0, width = .5)) + 
        scale_y_continuous(name=expression(paste("median area (", mu, "m^3)"))) +
        scale_x_continuous(name="tested temperature", limits=c(8, 42)) +
        theme_classic(base_size = 15) +
        # theme(legend.position = "none") + 
        ggtitle(paste("TAdapt:", allTAdaptation[aaa], "; C:", allMediumConcentrations[mmm], "; Line: ", allLines[lll])) +
        scale_color_continuous() +
        labs(color="Speed:")
      # scale_color_manual(values=c("#3B9AB2", "#EBCC2A", "#F21A00")) # this is the zissou1 palette
      if (saveFigures){
        ggsave(file=paste(currentTitle, "_area_vs_temp_jitter.png", sep=""), dpi = 600, width = 12, height = 10, units = "cm")
      }
      
      
      
      # I am running this on currentConditionTopSpeed because I think it makes more sense to do it only on the cells
      # that move fast and so are well aligned with the image.
      ggplot(currentConditionTopSpeed, aes(x=tTest, y=medianArea, fill=factor(tTest))) +
        geom_boxplot(width=0.8, outlier.shape = NA) + # add a box plot
        scale_y_continuous("median area") +
        scale_x_continuous(name="tested temperature", limits=c(10, 32.5)) +
        theme_classic(base_size = 15) +
        theme(legend.position = "none") + 
        ggtitle(paste("TAdapt:", allTAdaptation[aaa], "; C:", allMediumConcentrations[mmm], "; Line: ", allLines[lll])) +
        scale_color_continuous()
      if (saveFigures){
        ggsave(file=paste(currentTitle, "_area_vs_temp.png", sep=""), dpi = 600, width = 12, height = 10, units = "cm")
      }
      
      ggplot(currentConditionTopSpeed, aes(x=tTest, y=medianElongation, fill=factor(tTest))) +
        geom_boxplot(width=0.8, outlier.shape = NA) + # add a box plot
        scale_y_continuous("elongation") +
        scale_x_continuous(name="tested temperature", limits=c(10, 32.5)) +
        theme_classic(base_size = 15) +
        theme(legend.position = "none") + 
        ggtitle(paste("TAdapt:", allTAdaptation[aaa], "; C:", allMediumConcentrations[mmm], "; Line: ", allLines[lll])) +
        scale_color_continuous()
      if (saveFigures){
        ggsave(file=paste(currentTitle, "_elongation_vs_temp.png", sep=""), dpi = 600, width = 12, height = 10, units = "cm")
      }
      
      # currentConditionTopSpeed$estimatedlogVolume <- log10(4/3 * pi* currentConditionTopSpeed$medianBEllipse * currentConditionTopSpeed$medianAEllipse * currentConditionTopSpeed$medianAEllipse / 8)
      ggplot(currentConditionTopSpeed, aes(x=tTest, y=estimatedLogVolume, fill=factor(tTest))) +
        geom_boxplot(width=0.8, outlier.shape = NA) + # add a box plot
        scale_y_continuous(name=expression(paste('log'[10]*'(volume)'," ",  mu, 'm'^3))) +
        scale_x_continuous(name="tested temperature", limits=c(10, 32.5)) +
        theme_classic(base_size = 15) +
        theme(legend.position = "none") + 
        ggtitle(paste("TAdapt:", allTAdaptation[aaa], "; C:", allMediumConcentrations[mmm], "; Line: ", allLines[lll])) +
        scale_color_continuous()
      if (saveFigures){
        ggsave(file=paste(currentTitle, "_logVolume_vs_temp.png", sep=""), dpi = 600, width = 12, height = 10, units = "cm")
      }
      
      ggplot(currentConditionTopSpeed, aes(x=tTest, y=estimatedLogVolume, color=medianSpeed)) +
        geom_jitter(position = position_jitter(height = 0, width = .5)) + 
        scale_y_continuous(name=expression(paste('log'[10]*'(volume)'," ",  mu, 'm'^3))) +
        scale_x_continuous(name="tested temperature", limits=c(10, 32.5)) +
        theme_classic(base_size = 15) +
        # theme(legend.position = "none") + 
        ggtitle(paste("TAdapt:", allTAdaptation[aaa], "; C:", allMediumConcentrations[mmm], "; Line: ", allLines[lll])) +
        scale_color_continuous() +
        labs(color="speed:")
      
      if (saveFigures){
        ggsave(file=paste(currentTitle, "_logVolume_vs_temp_jitter.png", sep=""), dpi = 600, width = 12, height = 10, units = "cm")
      }
      
      ### Look at changes of volume over time
      
      # I first group together times that are on the same day
      currentConditionTopSpeed$incubationDurationInDays <- round((currentConditionTopSpeed$incubationDuration)/24)
      # currentConditionTopSpeed$tTest_as_factor <- factor(currentConditionTopSpeed$tTest, levels=paste(as.character(sort(unique(currentConditionTopSpeed$tTest)))))
      currentConditionTopSpeed$tTest_as_factor <- factor(currentConditionTopSpeed$tTest, levels=as.character(seq(12.5, 30, by=2.5)))
      
      currentConditionTopSpeedSummary <- data_summary(currentConditionTopSpeed, varname="estimatedLogVolume", 
                                                      groupnames=c("incubationDurationInDays", "tTest_as_factor", "tAdaptation", "mediumConcentration"))
      
      # if (plotCounterOverTime == 1)
      # {
      #   combinedCurrentConditionTopSpeedSummary <- currentConditionTopSpeedSummary
      # } else{
      #   combinedCurrentConditionTopSpeedSummary <- rbind(combinedCurrentConditionTopSpeedSummary, currentConditionTopSpeedSummary)
      # }
      
      plotColours <- rep(c("#3B9AB2", "#EBCC2A", "#F21A00"),5)
      markerShapes = c(21, 24, 22, 4)
      
      plotOverTime <- ggplot(currentConditionTopSpeedSummary, aes(x=incubationDurationInDays, y=estimatedLogVolume)) + # , color=tTest_as_factor)) +
        geom_hline(yintercept=currentConditionTopSpeedSummary$estimatedLogVolume[currentConditionTopSpeedSummary$tTest_as_factor == allTAdaptation[aaa] & currentConditionTopSpeedSummary$incubationDurationInDays == 0], colour=plotColours[aaa]) +
        geom_errorbar(aes(ymin=estimatedLogVolume-sd, ymax=estimatedLogVolume+sd), width=1.3) + 
        geom_point(size=3, color="black", fill=plotColours[aaa], shape=markerShapes[mmm]) + 
        # geom_line() +
        # geom_smooth(method=lm, formula='y~x', se=FALSE) + # geom_smooth should be done on the original data rather than on these data with errorbar
        # scale_y_continuous(name=expression(paste('log'[10]*'(volume)'," ",  mu, 'm'^3))) +
        scale_y_continuous(name=expression(paste('log'[10]*'(volume)'," ",  mu, 'm'^3)), limits=c(3.8, 4.8), sec.axis = sec_axis(trans=~.*1, name=paste("T:", allTAdaptation[aaa], "°C [", allMediumConcentrations[mmm], "%]", sep="") )) +
        scale_x_continuous(name="Incubation duration (days)", limits=c(-2,10), breaks=seq(0,10,by=5)) +
        theme_classic(base_size = 15) +
        theme(axis.ticks.y.right = element_blank(), axis.text.y.right = element_blank(), axis.line.y.right=element_blank()) +
        # theme(legend.position = "none") + 
        # ggtitle(paste("TAdapt:", allTAdaptation[aaa], "; C:", allMediumConcentrations[mmm], "; Line: ", allLines[lll])) +
        # scale_colour_manual(values=c("#3B9AB2", "#5DAABC", "#88BAAE", "#CAC656", "#E8C31E", "#E2B306", "#E86F00", "#F21A00")) +
        scale_color_manual(values=plotColours) + # this is the zissou1 palette
        scale_shape_manual(values=markerShapes) + # shapes for the markers
        facet_grid(cols=vars(tTest_as_factor),drop=FALSE) +
        labs(color="T test:")
      
      plotOverTime
      plotListOverTime[[plotCounterOverTime]] <- plotOverTime
      plotCounterOverTime <- plotCounterOverTime + 1
      
      if (saveFigures){
        library(Cairo)
        ggsave(file=paste(currentTitle, "_logVolume_vs_time.png", sep=""), dpi = 600, width = 20, height = 8, units = "cm")
        ggsave(file=paste(currentTitle, "_logVolume_vs_time.eps", sep=""), device="eps", dpi = 600, width = 20, height = 8, units = "cm")
        ggsave(file=paste(currentTitle, "_logVolume_vs_time.pdf", sep=""), device=cairo_pdf, dpi = 600, width = 20, height = 8, units = "cm")
      }
      
      ggplot(currentConditionTopSpeed, aes(x=tTest, y=medianElongation, color=medianSpeed)) +
        geom_jitter(position = position_jitter(height = 0, width = .5)) + 
        # geom_jitter(position = position_jitter(height = 0, width = .3)) +
        # geom_smooth(method="lm", se=TRUE, fill=NA, formula=y ~ poly(x, 3, raw=TRUE)) + 
        scale_y_continuous(name="median elongation") +
        scale_x_continuous(name="tested temperature", limits=c(10, 32.5)) +
        theme_classic(base_size = 15) +
        # theme(legend.position = "none") + 
        ggtitle(paste("TAdapt:", allTAdaptation[aaa], "; C:", allMediumConcentrations[mmm], "; Line: ", allLines[lll])) +
        scale_color_continuous(type = "viridis", option="magma") +
        labs(color="Speed:")
      # scale_color_manual(values=c("#3B9AB2", "#EBCC2A", "#F21A00")) # this is the zissou1 palette
      
      
      if (saveFigures){
        ggsave(file=paste(currentTitle, "_elongation_vs_temp_jitter_color_speed.png", sep=""), dpi = 600, width = 12, height = 10, units = "cm")
      }
      
      
      # currentConditionTopSpeed$estimatedVolume <- (4/3 * pi* currentConditionTopSpeed$medianBEllipse * currentConditionTopSpeed$medianAEllipse * currentConditionTopSpeed$medianAEllipse / 8)
      ggplot(currentConditionTopSpeed, aes(x=tTest, y=estimatedVolume, color=medianSpeed)) +
        geom_jitter(position = position_jitter(height = 0, width = .5)) + 
        # geom_jitter(position = position_jitter(height = 0, width = .3)) +
        # geom_smooth(method="lm", se=TRUE, fill=NA, formula=y ~ poly(x, 3, raw=TRUE)) + 
        scale_y_continuous(name=expression(paste("volume (",  mu, 'm'^3,")"))) +
        scale_x_continuous(name="tested temperature", limits=c(10, 32.5)) +
        theme_classic(base_size = 15) +
        # theme(legend.position = "none") + 
        ggtitle(paste("TAdapt:", allTAdaptation[aaa], "; C:", allMediumConcentrations[mmm], "; Line: ", allLines[lll])) +
        scale_color_continuous(type = "viridis", option="magma") +
        labs(color="Speed:")
      # scale_color_manual(values=c("#3B9AB2", "#EBCC2A", "#F21A00")) # this is the zissou1 palette
      
      
      if (saveFigures){
        ggsave(file=paste(currentTitle, "_volume_vs_temp_jitter_color_speed.png", sep=""), dpi = 600, width = 12, height = 10, units = "cm")
      }
      
      # library("wesanderson")
      # colourPalette <- wes_palette("Zissou1", 8, "continuous")
      
      ggplot(currentConditionTopSpeed, aes(x=tTest, y=estimatedLogVolume, fill=factor(tTest))) +
        geom_violin(width=2) + # add a violin plot
        geom_boxplot(width=0.4, outlier.shape = NA, color="black") + # add a violin plot
        scale_y_continuous(name=expression(paste('log'[10]*'(volume)'," ",  mu, 'm'^3))) +
        scale_x_continuous(name="Temp. (°C)", limits=c(10, 32.5)) +
        theme_classic(base_size = 15) +
        theme(legend.position = "none") +
        scale_fill_manual(values=c("#3B9AB2", "#5DAABC", "#88BAAE", "#CAC656", "#E8C31E", "#E2B306", "#E86F00", "#F21A00"))
      
      if (saveFigures){
        ggsave(file=paste(currentTitle, "_logvolume_vs_temp_violin.png", sep=""), dpi = 600, width = 12, height = 10, units = "cm")
        ggsave(file=paste(currentTitle, "_logvolume_vs_temp_violin.eps", sep=""), device="eps", dpi = 600, width = 12, height = 10, units = "cm")
        library(Cairo)
        ggsave(file=paste(currentTitle, "_logvolume_vs_temp_violin.pdf", sep=""), device=cairo_pdf, dpi = 1200, width = 12, height = 10, units = "cm")
      }
      
      
      # currentConditionTopSpeed$estimatedVolume <- (4/3 * pi* currentConditionTopSpeed$medianBEllipse * currentConditionTopSpeed$medianAEllipse * currentConditionTopSpeed$medianAEllipse / 8)
      ggplot(currentConditionTopSpeed, aes(x=tTest, y=oneOverVolume, color=medianSpeed)) +
        geom_jitter(position = position_jitter(height = 0, width = .5)) + 
        # geom_jitter(position = position_jitter(height = 0, width = .3)) +
        # geom_smooth(method="lm", se=TRUE, fill=NA, formula=y ~ poly(x, 3, raw=TRUE)) + 
        scale_y_continuous(name=expression(paste("1/V (",  mu, 'm'^-3,")"))) +
        scale_x_continuous(name="tested temperature", limits=c(10, 32.5)) +
        theme_classic(base_size = 15) +
        # theme(legend.position = "none") + 
        ggtitle(paste("TAdapt:", allTAdaptation[aaa], "; C:", allMediumConcentrations[mmm], "; Line: ", allLines[lll])) +
        scale_color_continuous(type = "viridis", option="magma") +
        labs(color="Speed:")
      # scale_color_manual(values=c("#3B9AB2", "#EBCC2A", "#F21A00")) # this is the zissou1 palette
      
      
      if (saveFigures){
        ggsave(file=paste(currentTitle, "_oneovervolume_vs_temp_jitter_color_speed.png", sep=""), dpi = 600, width = 12, height = 10, units = "cm")
      }
      
      
      
      
      
      currentConditionUnder30 <- subset(currentConditionTopSpeed, tTest < 30)
      
      
      ggplot(currentConditionUnder30, aes(x=oneOverkTrel, y=log(oneOverVolume), fill=line)) +
        # geom_point() + 
        geom_jitter(position = position_jitter(height = 0, width = .5)) +
        geom_smooth(method="lm", se=TRUE, fill="red", formula=y ~ x, na.rm=TRUE) + 
        scale_y_continuous(name=expression(paste("ln(1/V) (",  mu, 'm'^-3,")", sep=""))) +
        scale_x_continuous(name="1/kT0 - 1/kT") +
        theme_classic(base_size = 15) +
        theme(legend.position = "none") + 
        ggtitle(paste("TAdapt:", allTAdaptation[aaa], "; C:", allMediumConcentrations[mmm], "; Line: ", allLines[lll])) +
        scale_color_manual(values=c("#3B9AB2", "#EBCC2A", "#F21A00")) # this is the zissou1 palette
      if (saveFigures){
        ggsave(file=paste(currentTitle, "_arrhenius.png", sep=""), dpi = 600, width = 12, height = 10, units = "cm")
      }
      
      currentConditionUnder30 <- subset(currentConditionUnder30, is.finite(log(oneOverVolume)) & is.finite(oneOverkT))
      lm(log(oneOverVolume) ~ oneOverkTrel, data=currentConditionUnder30)
      
      
      
      
      
      
      
      
      ###### Fit schoolfield to 1/volume         
      
      
      if(!require(dplyr)){install.packages('dplyr')}
      d <- currentConditionTopSpeed %>% 
        dplyr::select("line", "tAdaptation", "oneOverVolume", "tTest", "estimatedVolume")%>%
        dplyr::rename(
          curve_id = line,
          growth_temp = tAdaptation,
          rate = oneOverVolume,
          temp = tTest,
          estimatedVolume = estimatedVolume
        )
      
      
      
      d <- d[complete.cases(d), ]
      # I added this because the get_start_vals function does not seem to work with na values
      # and when using w per cell as units the numbers are outside of the default boundaries
      
      # I add an artificial data point for temperature 40 degrees and rate=0
      dExtra <- d[1,]
      dExtra$temp <- 40
      dExtra$rate <- 0
      dExtra$estimatedVolume <- Inf
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
      
      calculatedFitParameters$tAdaptation <- allTAdaptation[aaa]
      calculatedFitParameters$line <- allLines[lll]
      calculatedFitParameters$mediumConcentration <- allMediumConcentrations[mmm]
      
      
      if (is.null(fit))
      {
        calculatedFitParameters <- cbind(calculatedFitParameters, data.frame(r_tref=NA, e_from_fit=NA, eh_from_fit=NA, th=NA))
      } else{
        calculatedFitParameters <- cbind(calculatedFitParameters, data.frame(r_tref=coef(fit)[1], e_from_fit=coef(fit)[2], eh_from_fit=coef(fit)[3], th=coef(fit)[4]))
      }
      
      print(paste("TAdapt:", allTAdaptation[aaa], "; C:", allMediumConcentrations[mmm], "; Line: ", allLines[lll]))
      print(calculatedFitParameters)
      
      print(fit)
      
      
      
      
      if (includeBootstrap)
      {
        library(minpack.lm)
        # refit model using nlsLM
        if (is.null(fit))
        {
          startingValues <- c(0.00013, 0.3, 20, 37) #
          # Hopefully this should not happen here!!!
          print(paste("Be careful: the fit is null!!!"))
          names(startingValues) <- c("r_tref", "e", "eh", "th")
        } else {
          startingValues <- coef(fit)
        }
        
        fit_nlsLM <- minpack.lm::nlsLM(rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 20),
                                       data = d,
                                       start = startingValues,
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
        allFitResults <- calculatedFitParameters[0,]
      }
      
      allFitResults <- rbind(allFitResults, calculatedFitParameters)
      
      if (!is.null(fit))
      {
        # predict new data
        new_data <- data.frame(temp = seq(min(d$temp), max(d$temp), 0.5))
        preds <- augment(fit, newdata = new_data)
        
        if (!is_empty(preds))
        {
          
          lineColours = c("#3B9AB2", "#EBCC2A", "#F21A00", "#FF00FF")
          markerColours = c("#A3A3A3", "#666666", "#000000", "#FF00FF")
          markerShapes = c(16,17,15,4)
          
          
          # plot data and model fit
          thisPlot <- ggplot(d, aes(temp, rate)) +
            # geom_point(col=markerColours[mmm], shape=markerShapes[mmm]) +
            geom_jitter(position = position_jitter(height = 0, width = .5), col=markerColours[mmm], shape=markerShapes[mmm]) +
            # annotate("text", size=5, x=10.5, y=920, label= paste("E=", round(calculatedFitParameters$e_from_fit, 2), "eV", sep=""), hjust = 0, parse=F) +
            geom_line(aes(temp, .fitted), preds, col = lineColours[aaa], size=2) +
            theme_classic(base_size = 15) +
            scale_x_continuous(name="Temp. (°C)",  limits=c(7.5, 42.5)) +
            scale_y_continuous(name=expression(paste("1/V (",  mu, 'm'^-3,")", sep=""))) +
            # labs(x = 'Temp (ºC)',  y = 'Stokes power (nW)') +
            ggtitle(paste("TAdapt:", allTAdaptation[aaa], "; C:", allMediumConcentrations[mmm], "; Line: ", allLines[lll]))
          
          thisPlot
          
          if (saveFigures){
            ggsave(file=paste(currentTitle, "_schoolfield.png", sep=""), dpi = 600, width = 12, height = 10, units = "cm")
          }
        }
      } else {
        # plot data only
        
        lineColours = c("#3B9AB2", "#EBCC2A", "#F21A00", "#FF00FF")
        markerColours = c("#A3A3A3", "#666666", "#000000", "#FF00FF")
        markerShapes = c(16,17,15,4)
        
        
        # plot data and model fit
        thisPlot <- ggplot(d, aes(temp, rate)) +
          # geom_point(col=markerColours[mmm], shape=markerShapes[mmm]) +
          geom_jitter(position = position_jitter(height = 0, width = .5), col=markerColours[mmm], shape=markerShapes[mmm]) +
          # annotate("text", size=5, x=10.5, y=920, label= paste("E=", round(calculatedFitParameters$e_from_fit, 2), "eV", sep=""), hjust = 0, parse=F) +
          # geom_line(aes(temp, .fitted), preds, col = lineColours[aaa], size=2) +
          theme_classic(base_size = 15) +
          scale_x_continuous(name="Temp. (°C)",  limits=c(7.5, 42.5)) +
          scale_y_continuous(name=expression(paste("1/V (",  mu, 'm'^-3,")", sep=""))) +
          # labs(x = 'Temp (ºC)',  y = 'Stokes power (nW)') +
          ggtitle(paste("TAdapt:", allTAdaptation[aaa], "; C:", allMediumConcentrations[mmm], "; Line: ", allLines[lll]))
        
        thisPlot
        if (saveFigures){
          ggsave(file=paste(currentTitle, "_schoolfield.png", sep=""), dpi = 600, width = 12, height = 10, units = "cm")
        }
      } 
      
      
      
      
      
      
      
      # Summarize the data in order to plot the data with error bar
      d2 <- data_summary(d, varname="rate", 
                         groupnames=c("growth_temp", "temp"))
      # head(d2)
      
      lineColours = c("#3B9AB2", "#EBCC2A", "#F21A00", "#FF00FF")
      markerColours = c("#A3A3A3", "#666666", "#000000", "#FF00FF")
      markerShapes = c(16,17,15,4)
      
      # plot data and model fit
      thisPlot <- ggplot(d2, aes(temp, rate))
      if (!is.null(fit))
      { thisPlot <- thisPlot + 
        geom_line(aes(temp, .fitted), preds, col = lineColours[aaa], size=2)
      }
      thisPlot <- thisPlot +
        geom_point(size=3, col=markerColours[mmm], shape=markerShapes[mmm]) +
        annotate("text", size=5, x=10.5, y=12e-05, label= paste("E=", round(calculatedFitParameters$e_from_fit, 2), "eV", sep=""), hjust = 0, parse=F) +
        geom_errorbar(aes(ymin=rate-sd, ymax=rate+sd), width=1.3, col=markerColours[mmm]) + 
        theme_classic(base_size = 14) +
        scale_x_continuous(name="Temp. (°C)",  limits=c(10, 32.5)) +
        scale_y_continuous(name=expression(paste("1/V (",  mu, 'm'^-3,")", sep="")), lim=c(0, 1.5e-4))
      # ggtitle(paste("T:", allTAdaptation[aaa], "; C:", allMediumConcentrations[mmm], "; L: ", allLines[lll]))
      
      thisPlot
      if (saveFigures){
        ggsave(file=paste(currentTitle, "_one_over_Volume_schoolfield_errorbar.png", sep=""), dpi = 600, width = 12, height = 10, units = "cm")
      }
      
      
      
      
      if (lll == 1){
        plotList[[plotCounter1]] <- thisPlot # add the first plot for each condition to a figure
        plotCounter1 <- plotCounter1 + 1
      } else if (lll == 2)
      {
        plotList2[[plotCounter2]] <- thisPlot # add the second plot for each condition to a figure
        plotCounter2 <- plotCounter2 + 1
      }else if (lll == 3)
      {
        plotList3[[plotCounter3]] <- thisPlot # add the third plot for each condition to a figure
        plotCounter3 <- plotCounter3 + 1
      }else if (lll == 4)
      {
        plotList4[[plotCounter4]] <- thisPlot # add the fourth plot for each condition to a figure
        plotCounter4 <- plotCounter4 + 1
      }
      
      
      
      
      
      # Summarize the data in order to plot the data with error bar
      d3 <- data_summary(d, varname="estimatedVolume", 
                         groupnames=c("growth_temp", "temp"))
      # head(d3)
      
      lineColours = c("#3B9AB2", "#EBCC2A", "#F21A00", "#FF00FF")
      markerColours = c("#A3A3A3", "#666666", "#000000", "#FF00FF")
      markerShapes = c(16,17,15,4)
      
      # plot data and model fit
      thisPlot <- ggplot(d3, aes(temp, estimatedVolume))
      if (!is.null(fit))
      { thisPlot <- thisPlot + 
        geom_line(aes(temp, 1/.fitted), preds, col = lineColours[aaa], size=2)
      }
      thisPlot <- thisPlot +
        geom_point(size=3, col=markerColours[mmm], shape=markerShapes[mmm]) +
        annotate("text", size=5, x=22.5, y=40000, label= paste("E=", round(calculatedFitParameters$e_from_fit, 2), "eV", sep=""), hjust = 0, parse=F) +
        geom_errorbar(aes(ymin=estimatedVolume - sd, ymax=estimatedVolume+sd), width=1.3, col=markerColours[mmm]) + 
        theme_classic(base_size = 14) +
        scale_x_continuous(name="Temp. (°C)",  limits=c(10, 32.5)) +
        scale_y_continuous(name=expression(paste("Volume (",  mu, 'm'^3,")", sep="")), lim=c(0, 50000))
      # ggtitle(paste("T:", allTAdaptation[aaa], "; C:", allMediumConcentrations[mmm], "; L: ", allLines[lll]))
      
      thisPlot
      # NOTE: the fit here goes below the points. This is because the points are the arithmetic mean of the data (mean(V))
      # while the fit is calculated on the inverse, so in a sense on the harmonic mean of the data 1/(mean(1/V)) and the 
      # harmonic mean is always smaller than the aritmetic mean.
      if (saveFigures){
        ggsave(file=paste(currentTitle, "_Volume_schoolfield_errorbar.png", sep=""), dpi = 600, width = 12, height = 10, units = "cm")
      }
      
      
      
      
      if (includeBootstrap)
      {
        # plot bootstrapped CIs
        thisPlot <- ggplot(d2, aes(temp, rate))
        if (!is.null(fit))
        { thisPlot <- thisPlot + 
          geom_line(aes(temp, .fitted), preds, col = lineColours[aaa], size=2)
        }
        thisPlot <- thisPlot +
          geom_line(aes(temp, pred, group = iter), boot1_preds, col = lineColours[aaa], alpha = 0.007) +
          geom_point(size=3, col=markerColours[mmm], shape=markerShapes[mmm]) +
          annotate("text", size=5, x=10.5, y=9e-05, label= paste("E=", round(calculatedFitParameters$e_from_fit, 2), "eV", sep=""), hjust = 0, parse=F) +
          geom_errorbar(aes(ymin=rate-sd, ymax=rate+sd), width=1.3, col=markerColours[mmm]) + 
          theme_classic(base_size = 14) +
          scale_x_continuous(name="Temp. (°C)",  limits=c(10, 42.5)) +
          scale_y_continuous(name=expression(paste("1/V (",  mu, 'm'^-3,")", sep="")), lim=c(1e-5, 3e-4))
        thisPlot
        if (saveFigures){
          ggsave(file=paste(currentTitle, "_Volume_schoolfield_errorbar_and_bootstrap.png", sep=""), dpi = 600, width = 12, height = 10, units = "cm")
        }
        
        # plot bootstrapped CIs
        thisPlot <- ggplot()
        if (!is.null(fit))
        { thisPlot <- thisPlot + 
          geom_line(aes(temp, .fitted), preds, col = lineColours[aaa], size=2)
        }
        thisPlot <- thisPlot +
          geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot1_conf_preds, fill = lineColours[aaa], alpha = 0.3) +
          geom_point(data=d2, aes(temp, rate), size=3, col=markerColours[mmm], shape=markerShapes[mmm]) +
          annotate("text", size=5, x=10.5, y=9e-05, label= paste("E=", round(calculatedFitParameters$e_from_fit, 2), "eV", sep=""), hjust = 0, parse=F) +
          geom_errorbar(data=d2, aes(x=temp, ymin=rate-sd, ymax=rate+sd), width=1.3, col=markerColours[mmm]) + 
          theme_classic(base_size = 14) +
          scale_x_continuous(name="Temp. (°C)",  limits=c(7.5, 42.5)) +
          scale_y_continuous(name=expression(paste("1/V (",  mu, 'm'^-3,")", sep="")), lim=c(1e-5, 3e-4))
        thisPlot
        if (saveFigures){
          ggsave(file=paste(currentTitle, "_Volume_schoolfield_errorbar_and_bootstrap2.png", sep=""), dpi = 600, width = 12, height = 10, units = "cm")
        }
        
        
      }
      
      
      
    }
  }
}



print(allFitResults)

if (combineLinesTogether)
{
  write.table(allFitResults, "allFitResults_Schoolfield_on_one_over_Volume_long_term_response_all.csv", append = FALSE, sep = ", ", row.names = FALSE)
} else {
  write.table(allFitResults, "allFitResults_Schoolfield_on_one_over_Volume_long_term_response.csv", append = FALSE, sep = ", ", row.names = FALSE)
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

if(saveFigures)
{
  if (combineLinesTogether)
  {
    ggsave(file="figure_one_over_Volume_long_term_thermal_response_post_adaptation_all.png", dpi = 600, width = 24, height = 20, units = "cm")
    ggsave(file="figure_one_over_Volume_long_term_thermal_response_post_adaptation_all.eps", device="eps", dpi = 1200, width = 24, height = 20, units = "cm")
    ggsave(file="figure_one_over_Volume_thermal_response_post_adaptation_all.pdf", device=cairo_pdf, dpi = 1200, width = 24, height = 20, units = "cm")
    
  }  else {
    ggsave(file="figure_one_over_Volume_long_term_thermal_response_post_adaptation.png", dpi = 600, width = 24, height = 20, units = "cm")
    ggsave(file="figure_one_over_Volume_long_term_thermal_response_post_adaptation.eps", device="eps", dpi = 1200, width = 24, height = 20, units = "cm")
    ggsave(file="figure_one_over_Volume_thermal_response_post_adaptation.pdf", device=cairo_pdf, dpi = 1200, width = 24, height = 20, units = "cm")
  }
}

if (combineLinesTogether == FALSE)
{
  # combine together the two plots into a single figure
  fullFigure <- ggarrange(plotList2[[1]], plotList2[[2]], plotList2[[3]],
                          plotList2[[4]], plotList2[[5]], plotList2[[6]],
                          plotList2[[7]], plotList2[[8]], plotList2[[9]],
                          # labels = c("15-50", "15-100", "15-200", "20-50", "20-100", "20-200", "25-50", "25-100", "25-200"),
                          ncol = 3, nrow = 3, common.legend=TRUE)
  
  fullFigure
  if(saveFigures)
  {
    ggsave(file="figure_one_over_Volume_long_term_thermal_response_post_adaptation2.png", dpi = 600, width = 24, height = 20, units = "cm")
    ggsave(file="figure_one_over_Volume_long_term_thermal_response_post_adaptation2.eps", device="eps", dpi = 1200, width = 24, height = 20, units = "cm")
    library(Cairo)
    ggsave(file="figure_one_over_Volume_thermal_response_post_adaptation2.pdf", device=cairo_pdf, dpi = 1200, width = 24, height = 20, units = "cm")
  }
  
  # combine together the two plots into a single figure
  fullFigure <- ggarrange(plotList3[[1]], plotList3[[2]], plotList3[[3]],
                          plotList3[[4]], plotList3[[5]], plotList3[[6]],
                          plotList3[[7]], plotList3[[8]], plotList3[[9]],
                          # labels = c("15-50", "15-100", "15-200", "20-50", "20-100", "20-200", "25-50", "25-100", "25-200"),
                          ncol = 3, nrow = 3, common.legend=TRUE)
  
  fullFigure
  if (saveFigures)
  {
    ggsave(file="figure_one_over_Volume_long_term_thermal_response_post_adaptation3.png", dpi = 600, width = 24, height = 20, units = "cm")
    ggsave(file="figure_one_over_Volume_long_term_thermal_response_post_adaptation3.eps", device="eps", dpi = 1200, width = 24, height = 20, units = "cm")
    library(Cairo)
    ggsave(file="figure_one_over_Volume_thermal_response_post_adaptation3.pdf", device=cairo_pdf, dpi = 1200, width = 24, height = 20, units = "cm")
  }
  
  # combine together the two plots into a single figure
  fullFigure <- ggarrange(plotList4[[1]], plotList4[[2]], plotList4[[3]],
                          plotList4[[4]], plotList4[[5]], plotList4[[6]],
                          plotList4[[7]], plotList4[[8]], plotList4[[9]],
                          # labels = c("15-50", "15-100", "15-200", "20-50", "20-100", "20-200", "25-50", "25-100", "25-200"),
                          ncol = 3, nrow = 3, common.legend=TRUE)
  
  fullFigure
  if (saveFigures)
  {
    ggsave(file="figure_one_over_Volume_long_term_thermal_response_post_adaptation4.png", dpi = 600, width = 24, height = 20, units = "cm")
    ggsave(file="figure_one_over_Volume_long_term_thermal_response_post_adaptation4.eps", device="eps", dpi = 1200, width = 24, height = 20, units = "cm")
    library(Cairo)
    ggsave(file="figure_one_over_Volume_thermal_response_post_adaptation4.pdf", device=cairo_pdf, dpi = 1200, width = 24, height = 20, units = "cm")
  }
  
  
}


fullFigureOverTime <- ggarrange(plotListOverTime[[1]]+theme(axis.title.x = element_blank(), strip.background=element_blank(), strip.text=element_blank()), plotListOverTime[[2]]+theme(axis.title.x = element_blank(), strip.background=element_blank(), strip.text=element_blank()), plotListOverTime[[3]]+theme(axis.title.x = element_blank(), strip.background=element_blank(), strip.text=element_blank()),
                                plotListOverTime[[4]]+theme(axis.title.x = element_blank(), strip.background=element_blank(), strip.text=element_blank()), plotListOverTime[[5]]+theme(axis.title.x = element_blank(), strip.background=element_blank(), strip.text=element_blank()), plotListOverTime[[6]]+theme(axis.title.x = element_blank(), strip.background=element_blank(), strip.text=element_blank()),
                                plotListOverTime[[7]]+theme(axis.title.x = element_blank(), strip.background=element_blank(), strip.text=element_blank()), plotListOverTime[[8]]+theme(axis.title.x = element_blank(), strip.background=element_blank(), strip.text=element_blank()), plotListOverTime[[9]]+theme(axis.title.x = element_blank(), strip.background=element_blank(), strip.text=element_blank()),
                                ncol = 1, nrow = 9, common.legend=TRUE)

fullFigureOverTime

if(saveFigures)
{
  if (combineLinesTogether)
  {
    ggsave(file="figure_log_volume_over_time_all.png", dpi = 600, width = 20, height = 72, units = "cm")
    ggsave(file="figure_log_volume_over_time_all.eps", device="eps", dpi = 1200, width = 20, height = 72, units = "cm")
    ggsave(file="figure_log_volume_over_time_all.pdf", device=cairo_pdf, dpi = 1200, width = 20, height = 72, units = "cm")
  }
}

# # Alternative approach to get the same figure:
# 
# plotColours <- rep(c("#3B9AB2", "#EBCC2A", "#F21A00"),5)
# markerShapes = c(21, 24, 22, 4)
# 
# plotOverTime <- ggplot(combinedCurrentConditionTopSpeedSummary, aes(x=incubationDurationInDays, y=estimatedLogVolume)) + # , color=tTest_as_factor)) +
#   geom_hline(yintercept=combinedCurrentConditionTopSpeedSummary$estimatedLogVolume[combinedCurrentConditionTopSpeedSummary$tTest_as_factor == allTAdaptation[aaa] & combinedCurrentConditionTopSpeedSummary$incubationDurationInDays == 0], colour=plotColours[aaa]) +
#   geom_errorbar(aes(ymin=estimatedLogVolume-sd, ymax=estimatedLogVolume+sd), width=1.3) + 
#   geom_point(size=3, color="black", fill=plotColours[aaa], shape=markerShapes[mmm]) + 
#   # geom_line() +
#   # geom_smooth(method=lm, formula='y~x', se=FALSE) + # geom_smooth should be done on the original data rather than on these data with errorbar
#   # scale_y_continuous(name=expression(paste('log'[10]*'(volume)'," ",  mu, 'm'^3))) +
#   scale_y_continuous(name=expression(paste('log'[10]*'(volume)'," ",  mu, 'm'^3)), limits=c(3.8, 4.8), sec.axis = sec_axis(trans=~.*1, name=paste("T:", allTAdaptation[aaa], "°C [", allMediumConcentrations[mmm], "%]", sep="") )) +
#   scale_x_continuous(name="Incubation duration (days)", limits=c(-2,10), breaks=seq(0,10,by=5)) +
#   theme_classic(base_size = 15) +
#   theme(axis.ticks.y.right = element_blank(), axis.text.y.right = element_blank(), axis.line.y.right=element_blank()) +
#   # theme(legend.position = "none") + 
#   # ggtitle(paste("TAdapt:", allTAdaptation[aaa], "; C:", allMediumConcentrations[mmm], "; Line: ", allLines[lll])) +
#   # scale_colour_manual(values=c("#3B9AB2", "#5DAABC", "#88BAAE", "#CAC656", "#E8C31E", "#E2B306", "#E86F00", "#F21A00")) +
#   scale_color_manual(values=plotColours) + # this is the zissou1 palette
#   scale_shape_manual(values=markerShapes) + # shapes for the markers
#   facet_grid(cols=vars(tTest_as_factor), rows=vars(mediumConcentration), drop=FALSE) +
#   labs(color="T test:")
# plotOverTime

# plot data and model fit
ggplot(allFitResults, aes(x=tAdaptation, y=e_from_fit, color=factor(mediumConcentration), shape=factor(mediumConcentration), label=factor(line))) +
  # geom_point(size=6) +
  # geom_jitter(size=3, width=0.4) +
  geom_jitter(size=4, position = position_jitter(height = 0, width = .3)) + 
  # geom_text(color="black", size=4) + 
  scale_color_manual(values=c("#A3A3A3", "#666666", "#000000")) +
  scale_fill_manual(values=c("#A3A3A3", "#666666", "#000000")) +
  scale_y_continuous(name="activation energy (eV)", limits=c(0, 1)) +
  scale_x_continuous(name="Adaptation temperature (ºC)", breaks=c(15,20,25), limits=c(12.5, 27.5)) + 
  theme_classic(base_size=15) +
  theme(legend.position = "none") +
  labs(color = "Conc. (%)") # + # this specifies a custom legend
# ggtitle("Activation energy on one over volume")

if (saveFigures){
  if (combineLinesTogether){
    ggsave(file="figure_activation_energy_on_one_over_Volume_long_term_all.png", dpi = 600, width = 12, height = 10, units = "cm")
    ggsave(file="figure_activation_energy_on_one_over_Volume_long_term_all.eps", device="eps", dpi = 1200, width = 12, height = 10, units = "cm")
    library(Cairo)
    ggsave(file="figure_activation_energy_on_one_over_Volume_long_term_all.pdf", device=cairo_pdf, dpi = 1200, width = 12, height = 10, units = "cm")
  } else {
    ggsave(file="figure_activation_energy_on_one_over_Volume_long_term.png", dpi = 600, width = 12, height = 10, units = "cm")
    ggsave(file="figure_activation_energy_on_one_over_Volume_long_term.eps", device="eps", dpi = 1200, width = 12, height = 10, units = "cm")
    library(Cairo)
    ggsave(file="figure_activation_energy_on_one_over_Volume_long_term.pdf", device=cairo_pdf, dpi = 1200, width = 12, height = 10, units = "cm")
    
  }
  
}




# plot data and model fit
ggplot(allFitResults, aes(x=tAdaptation, y=topt, color=factor(mediumConcentration), shape=factor(mediumConcentration), label=factor(line))) +
  # geom_point(size=6) +
  geom_jitter(size=4, position = position_jitter(height = 0, width = .3)) + 
  # geom_text(color="black", size=4) + 
  scale_color_manual(values=c("#A3A3A3", "#666666", "#000000")) +
  scale_fill_manual(values=c("#A3A3A3", "#666666", "#000000")) +
  scale_y_continuous(name="Optimal temperature (ºC)") +
  scale_x_continuous(name="Adaptation temperature (ºC)", limits=c(12.5, 27.5), breaks=c(15, 20, 25)) +
  theme_classic(base_size=15) +
  theme(legend.position = "none") +
  labs(color = "Conc. (%)") + # this specifies a custom legend
  ggtitle("Optimal temperature")

if (saveFigures){
  if (combineLinesTogether){
    ggsave(file="figure_optimal_temp_vs_temp_on_one_over_Volume_long_term_all.png", dpi = 600, width = 12, height = 10, units = "cm")
    ggsave(file="figure_optimal_temp_vs_temp_on_one_over_Volume_long_term_all.eps", device="eps", dpi = 1200, width = 12, height = 10, units = "cm")
    library(Cairo)
    ggsave(file="figure_optimal_temp_vs_temp_on_one_over_Volume_long_term_all.pdf", device=cairo_pdf, dpi = 1200, width = 12, height = 10, units = "cm")
  } else {
    ggsave(file="figure_optimal_temp_vs_temp_on_one_over_Volume_long_term.png", dpi = 600, width = 12, height = 10, units = "cm")
    ggsave(file="figure_optimal_temp_vs_temp_on_one_over_Volume_long_term.eps", device="eps", dpi = 1200, width = 12, height = 10, units = "cm")
    library(Cairo)
    ggsave(file="figure_optimal_temp_vs_temp_on_one_over_Volume_long_term.pdf", device=cairo_pdf, dpi = 1200, width = 12, height = 10, units = "cm")
    
  }
}




# plot data and model fit
ggplot(allFitResults, aes(x=tAdaptation, y=r_tref, color=factor(mediumConcentration), shape=factor(mediumConcentration), label=factor(line))) +
  # geom_point(size=4) +
  geom_jitter(size=4, position = position_jitter(height = 0, width = .3)) + 
  # geom_text(color="black", size=4) + 
  scale_color_manual(values=c("#A3A3A3", "#666666", "#000000")) +
  scale_fill_manual(values=c("#A3A3A3", "#666666", "#000000")) +
  scale_y_continuous(name=expression(paste("1/V (",  mu, 'm'^-3,")", sep=""))) +
  scale_x_continuous(name=expression("Adaptation temperature (ºC)"), limits=c(12.5, 27.5), breaks=c(15, 20, 25)) +
  theme_classic(base_size=12) +
  theme(legend.position = "none") +
  labs(color = "Conc.") + # this specifies a custom legend
  ggtitle(paste("1/V at reference temperature (", referenceTemperature, "°C)", sep=""))

if (saveFigures)
{
  if (combineLinesTogether){
    ggsave(file="figure_one_over_Volume_at_tref_long_term_all.png", dpi = 600, width = 12, height = 10, units = "cm")
    ggsave(file="figure_one_over_Volume_at_tref_long_term_all.eps", device="eps", dpi = 1200, width = 12, height = 10, units = "cm")
    # quartz(type = 'pdf', file="aaaa_metabolic_thermal_response_post_adaptation.pdf")
    library(Cairo)
    ggsave(file="figure_one_over_Volume_at_tref_long_term_all.pdf", device=cairo_pdf, dpi = 1200, width = 12, height = 10, units = "cm")
  } else {
    ggsave(file="figure_one_over_Volume_at_tref_long_term.png", dpi = 600, width = 12, height = 10, units = "cm")
    ggsave(file="figure_one_over_Volume_at_tref_long_term.eps", device="eps", dpi = 1200, width = 12, height = 10, units = "cm")
    library(Cairo)
    ggsave(file="figure_one_over_Volume_at_tref_long_term.pdf", device=cairo_pdf, dpi = 1200, width = 12, height = 10, units = "cm")
    
  }
}








if (includeBootstrap)
  
{
  # plot data and model fit
  ggplot(allFitResults, aes(x=tAdaptation, y=e_from_fit, color=factor(mediumConcentration), fill=factor(mediumConcentration + tAdaptation*100), shape=factor(mediumConcentration), label=factor(line))) +
    geom_errorbar(aes(ymin=e_from_fit - sd_e_from_fit, ymax=e_from_fit + sd_e_from_fit), width=0, position=position_dodge(width=2)) + 
    geom_point(size=6, position=position_dodge(width=2)) + 
    # geom_errorbar(aes(ymin=e_from_fit_0025, ymax=e_from_fit_0975), width=0, position=position_dodge(width=2)) + 
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
      ggsave(file="figure_activation_energy_on_one_over_Volume_bootstrap_all.png", dpi = 600, width = 12, height = 10, units = "cm")
      ggsave(file="figure_activation_energy_on_one_over_Volume_bootstrap_all.eps", device="eps", dpi = 1200, width = 12, height = 10, units = "cm")
      library(Cairo)
      ggsave(file="figure_activation_energy_on_one_over_Volume_bootstrap_all.pdf", device=cairo_pdf, dpi = 1200, width = 12, height = 10, units = "cm")
    } else {
      ggsave(file="figure_activation_energy_on_one_over_Volume_bootstrap.png", dpi = 600, width = 12, height = 10, units = "cm")
      ggsave(file="figure_activation_energy_on_one_over_Volume_bootstrap.eps", device="eps", dpi = 1200, width = 12, height = 10, units = "cm")
      library(Cairo)
      ggsave(file="figure_activation_energy_on_one_over_Volume_bootstrap.pdf", device=cairo_pdf, dpi = 1200, width = 12, height = 10, units = "cm")
      
    }
  }
  
  
  
  # plot data and model fit
  ggplot(allFitResults, aes(x=tAdaptation, y=r_tref, color=factor(mediumConcentration), fill=factor(mediumConcentration + tAdaptation*100), shape=factor(mediumConcentration), label=factor(line))) +
    geom_errorbar(aes(ymin=r_tref - sd_r_tref, ymax=r_tref + sd_r_tref), width=0, position=position_dodge(width=2)) + 
    geom_point(size=6, position=position_dodge(width=2)) + 
    # geom_errorbar(aes(ymin=r_from_fit_0025, ymax=r_from_fit_0975), width=0, position=position_dodge(width=2)) + 
    scale_color_manual(values=c("#A3A3A3", "#666666", "#000000")) +
    scale_fill_manual(values=c("#8FCDDC", "#3B9AB2", "#18434E", "#F2DD70", "#EBCC2A", "#463C07", "#FF7C6C", "#F21A00", "#660B00")) +
    scale_shape_manual(values=c(21, 24, 22)) + # shapes for the markers
    scale_y_continuous(name=expression(paste("1/V at T=20°C (",  mu, 'm'^-3,")", sep="")), limits=c(1e-05, 1e-04), sec.axis = sec_axis( trans=~log10(1/.), name=expression(paste("Log10(V) at T=20°C (",  mu, 'm'^3,")")))) +
    scale_x_continuous(name=expression("Adaptation temperature (ºC)"), limits=c(12.5, 27.5), breaks=c(15, 20, 25)) +
    theme_classic(base_size=18) +
    theme(legend.position = "none") +
    labs(color = "Conc.") # + # this specifies a custom legend
  
  
  # plot data and model fit
  ggplot(allFitResults, aes(x=tAdaptation, y=r_tref, color=factor(mediumConcentration), fill=factor(mediumConcentration + tAdaptation*100), shape=factor(mediumConcentration), label=factor(line))) +
    geom_errorbar(aes(ymin=r_tref - sd_r_tref, ymax=r_tref + sd_r_tref), width=0, position=position_dodge(width=2)) + 
    geom_point(size=6, position=position_dodge(width=2)) + 
    # geom_errorbar(aes(ymin=r_from_fit_0025, ymax=r_from_fit_0975), width=0, position=position_dodge(width=2)) + 
    scale_color_manual(values=c("#A3A3A3", "#666666", "#000000")) +
    scale_fill_manual(values=c("#8FCDDC", "#3B9AB2", "#18434E", "#F2DD70", "#EBCC2A", "#463C07", "#FF7C6C", "#F21A00", "#660B00")) +
    scale_shape_manual(values=c(21, 24, 22)) + # shapes for the markers
    scale_y_continuous(name=expression(paste("1/V at T=20°C (",  mu, 'm'^-3,")", sep="")), limits=c(1e-05, 1e-04)) +
    scale_x_continuous(name=expression("Adaptation temperature (ºC)"), limits=c(12.5, 27.5), breaks=c(15, 20, 25)) +
    theme_classic(base_size=18) +
    theme(legend.position = "none") +
    labs(color = "Conc.") # + # this specifies a custom legend
  # ggtitle(paste("1/V at reference temperature (", referenceTemperature, "°C)", sep=""))
  
  if (saveFigures)
  {
    if (combineLinesTogether){
      ggsave(file="figure_one_over_Volume_at_tref_bootstrap_all.png", dpi = 600, width = 12, height = 10, units = "cm")
      ggsave(file="figure_one_over_Volume_at_tref_bootstrap_all.eps", device="eps", dpi = 1200, width = 12, height = 10, units = "cm")
      library(Cairo)
      ggsave(file="figure_one_over_Volume_at_tref_bootstrap_all.pdf", device=cairo_pdf, dpi = 1200, width = 12, height = 10, units = "cm")
    } else {
      ggsave(file="figure_one_over_Volume_at_tref_bootstrap.png", dpi = 600, width = 12, height = 10, units = "cm")
      ggsave(file="figure_one_over_Volume_at_tref_bootstrap.eps", device="eps", dpi = 1200, width = 12, height = 10, units = "cm")
      library(Cairo)
      ggsave(file="figure_one_over_Volume_at_tref_bootstrap.pdf", device=cairo_pdf, dpi = 1200, width = 12, height = 10, units = "cm")
      
    }
  }
  
  
}





lineColours = c("#3B9AB2", "#EBCC2A", "#F21A00", "#FF00FF")
lineTypes <- c("twodash", "dashed", "solid", "longdash")
# Now try to plot the fit only
allFitData <- data.frame()
for (iii in 1:dim(allFitResults)[1])
{
  predictedData <- data.frame(temp = seq(10, 40, 0.5))
  predictedData$rate <- sharpeschoolhigh_1981(temp = predictedData$temp, r_tref = allFitResults$r_tref[iii],e = allFitResults$e_from_fit[iii], eh = allFitResults$eh_from_fit[iii], th = allFitResults$th[iii], tref = 20)
  predictedData$tAdaptation <- allFitResults$tAdaptation[iii]
  predictedData$mediumConcentration <- allFitResults$mediumConcentration[iii]
  predictedData$line <- allFitResults$line[iii]
  allFitData <- rbind(allFitData, predictedData)
}

plotG1 <- ggplot(allFitData, aes(x=temp, y=rate, color=as.factor(tAdaptation), size=as.factor(mediumConcentration), alpha=mediumConcentration)) +
  geom_line() +
  scale_x_continuous(name="Temp. (°C)",  limits=c(10, 25)) +
  scale_y_continuous(name=expression(paste("1/V (",  mu, 'm'^-3,")", sep="")), lim=c(1e-5, 1e-4), trans = 'log10') +
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
    ggsave(file="figure_fitted_one_over_Volume_vs_temp_all.png", dpi = 600, width = 12, height = 10, units = "cm")
    ggsave(file="figure_fitted_one_over_Volume_vs_temp_all.eps", device="eps", dpi = 1200, width = 12, height = 10, units = "cm")
    library(Cairo)
    ggsave(file="figure_fitted_one_over_Volume_vs_temp_all.pdf", device=cairo_pdf, dpi = 1200, width = 12, height = 10, units = "cm")
  } else {
    ggsave(file="figure_fitted_one_over_Volume_vs_temp.png", dpi = 600, width = 12, height = 10, units = "cm")
  }
}



plotG1 <- ggplot(allFitData, aes(x=temp, y=1/rate, color=as.factor(tAdaptation), size=as.factor(mediumConcentration), alpha=mediumConcentration)) + # linetype=as.factor(mediumConcentration), 
  geom_line() +
  scale_x_continuous(name="Temp. (°C)",  limits=c(10, 25)) +
  scale_y_continuous(name=expression(paste("volume (",  mu, 'm'^3,")")), trans = 'log10', lim=c(5000, 50000)) +
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
    ggsave(file="figure_fitted_Volume_vs_temp_all.png", dpi = 600, width = 12, height = 10, units = "cm")
    ggsave(file="figure_fitted_Volume_vs_temp_all.eps", device="eps", dpi = 1200, width = 12, height = 10, units = "cm")
    library(Cairo)
    ggsave(file="figure_fitted_Volume_vs_temp_all.pdf", device=cairo_pdf, dpi = 1200, width = 12, height = 10, units = "cm")
  } else {
    ggsave(file="figure_fitted_Volume_vs_temp.png", dpi = 600, width = 12, height = 10, units = "cm")
  }
}

# sessionInfo()
