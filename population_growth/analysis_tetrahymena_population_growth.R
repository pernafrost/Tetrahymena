rm(list=ls()) # clean memory


library(rTPC)
library(nls.multstart)
library(dplyr)
library(broom)
library(ggplot2)
library(stringr)



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
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}




conversionPixelPerMicrometre <- sqrt(779441.5)/2/1000
saveFigures <- TRUE
skipMotherCulture <- TRUE # whether to skip the mother culture from this analysis
includeBootstrap <- TRUE # whether to run a bootstrap on the fitted thermal response curve
dev.off() # sometimes R studio stops showing plots (possibly because you call other devices for pdf)

conversionMillilitresPerPixel <-  4 * 10^-4 / 779441.5
# the volume in which the cells are counted depends on the area selected under the microscope in pixels (d1$AreaPixels). Each square is 0.1 microliters, or 10^-4 milliliters and 4 squares are together 779441 pixels square: (X px) * (4 sq) * (10^-4 ml / sq) / (779441 px) = (Y ml)


fileName <- "/Tetrahymena/population_growth/Tetrahymena_pop_growth.csv"
setwd(dirname(fileName))
allExperimentResults <- read.table(file = fileName, sep = ",", header=TRUE, na.strings = c("NA", " NA"))



# some cleaning out of comments
allExperimentResults$Spores[is.na(allExperimentResults$Spores)] <- 0

allExperimentResults$filteredFileName <- str_replace(allExperimentResults$fileName, "_careful_measured22.5_", "_")
allExperimentResults$filteredFileName <- str_replace(allExperimentResults$filteredFileName, "_subculture_14_04_", "_")
allExperimentResults$filteredFileName <- str_replace(allExperimentResults$filteredFileName, "_contaminated_", "_")
allExperimentResults$filteredFileName <- str_replace(allExperimentResults$filteredFileName, "_Starting_14_04_", "_")
allExperimentResults$filteredFileName <- str_replace(allExperimentResults$filteredFileName, "_contamination_or_unusual_", "_")

# allExperimentResults$fileName[ which(is.na(str_match(allExperimentResults$fileName, "Starting"))==FALSE)]

# rename lines 1 and 2 in the initial folders as they do not have the same initial density
myInd <- which(is.na(str_match(allExperimentResults$Folder_Label, "2021_03_01")) == FALSE)
allExperimentResults$line[myInd]
# Line 1 is split into 1 and 3 for the growth experiment;
# Line 2 is split into 2 and 4 for the growth experiment.
allExperimentResults$line[myInd] <- str_replace(allExperimentResults$line[myInd], "1", "3")
allExperimentResults$line[myInd] <- str_replace(allExperimentResults$line[myInd], "2", "4")
allExperimentResults$filteredFileName[myInd] <- str_replace(allExperimentResults$filteredFileName[myInd], "_1_", "_3_")
allExperimentResults$filteredFileName[myInd] <- str_replace(allExperimentResults$filteredFileName[myInd], "_2_", "_4_")
allExperimentResults$filteredFileName[myInd] <- str_replace(allExperimentResults$filteredFileName[myInd], "_1bis_", "_3bis_")
allExperimentResults$filteredFileName[myInd] <- str_replace(allExperimentResults$filteredFileName[myInd], "_2bis_", "_4bis_")





# Extract information from the string in the "fileName" field
newColumns <- str_split_fixed(as.character(str_trim(tolower(basename(allExperimentResults$filteredFileName)))), "_|hours", n=Inf)[,c(2,3,4,6,7,9,11,13,16)]

# remove c for degrees celsius and remove bis from the line
newColumns[,4] <- str_replace(newColumns[,4], "c", "")
newColumns[,3] <- str_replace(newColumns[,3], "bis", "")

# class(newColumns) <- "numeric"
newColumns <- data.frame(tAdaptation = as.numeric(newColumns[,1]), 
                         mediumConcentration = as.numeric(newColumns[,2]),
                         tTest = as.numeric(newColumns[,4]), 
                         line = newColumns[,3],
                         incubationDuration = as.numeric(newColumns[,5]),
                         repetition = newColumns[,6],
                         frameNum = as.numeric(newColumns[,7]),
                         automaticCount = as.numeric(newColumns[,8]),
                         areaPx = as.numeric(newColumns[,9]))


# add the new columns to the original table
allExperimentResults <- cbind(allExperimentResults, newColumns, deparse.level = 1, stringsAsFactors = default.stringsAsFactors())
# sapply(allExperimentResults, class)

allExperimentResults$volMl <- allExperimentResults$areaPx * conversionMillilitresPerPixel

fileNameStartingDensities <- "/Tetrahymena/population_growth/Starting_Densities_for_Growth_Experiment.csv"
startingDensities <- read.table(file = fileNameStartingDensities, sep = ",", header=TRUE, na.strings = c("NA", " NA"))
fileNameStartingDensities30deg <- "/Tetrahymena/population_growth/Starting_Densities_for_Growth_Experiment_30deg.csv"
startingDensities30deg <- read.table(file = fileNameStartingDensities30deg, sep = ",", header=TRUE, na.strings = c("NA", " NA"))



# Now for each tested condition (e.g. a_15_50_c_growth25) I look at the initial density (which is the same
# across all those adapted in the same conditions)
# adaptedCondition <- paste("a_", allExperimentResults$tAdaptation, "_", allExperimentResults$mediumConcentration, "_", allExperimentResults$line,  sep="")
allExperimentResults$testedCondition <- paste("a_", allExperimentResults$tAdaptation, "_", allExperimentResults$mediumConcentration, "_", allExperimentResults$line, "_growth_", allExperimentResults$tTest, sep="")
testedCondition <- unique(allExperimentResults$testedCondition)
# testedCondition <- unique(paste("a_", allExperimentResults$tAdaptation, "_", allExperimentResults$mediumConcentration, "_", allExperimentResults$line, "_growth_", allExperimentResults$tTest, sep=""))



initialDensities <- allExperimentResults[0, ]
initialDensities [ nrow(initialDensities) + length(testedCondition), ] <- NA
initialDensities$fileName <-  paste(testedCondition, "_initialDensity", sep="") 
initialDensities$filteredFileName <- initialDensities$fileName
initialDensities$tAdaptation <- as.numeric(str_split_fixed(as.character(str_trim(tolower(initialDensities$filteredFileName))), "_", n=Inf)[,2])
initialDensities$mediumConcentration <- as.numeric(str_split_fixed(as.character(str_trim(tolower(initialDensities$filteredFileName))), "_", n=Inf)[,3])
initialDensities$line <- str_split_fixed(as.character(str_trim(tolower(initialDensities$filteredFileName))), "_", n=Inf)[,4]
initialDensities$tTest <- as.numeric(str_split_fixed(as.character(str_trim(tolower(initialDensities$filteredFileName))), "_", n=Inf)[,6])
initialDensities$testedCondition <- testedCondition

# now for each initialDensities line complete with more information
for (iii in 1:dim(initialDensities)[1])
{
  if (is.na(str_match(initialDensities$line[iii], "1")) == FALSE)
  {
    initialDensities$Particle_manual_count[iii] <- 1826
    initialDensities$incubationDuration[iii] <- 0
    initialDensities$volMl[iii] <- 1
  } else if (is.na(str_match(initialDensities$line[iii], "2")) == FALSE)
  {
    initialDensities$Particle_manual_count[iii] <- 3370
    initialDensities$incubationDuration[iii] <- 0
    initialDensities$volMl[iii] <- 1
  } else if (is.na(str_match(initialDensities$line[iii], "3")) == FALSE)
  {
    initialDensities$Particle_manual_count[iii] <- 1958
    initialDensities$incubationDuration[iii] <- 0
    initialDensities$volMl[iii] <- 1
  } else if (is.na(str_match(initialDensities$line[iii], "4")) == FALSE)
  {
    initialDensities$Particle_manual_count[iii] <- 2805
    initialDensities$incubationDuration[iii] <- 0
    initialDensities$volMl[iii] <- 1
  } else if (initialDensities$tTest[iii] == 30)   # check if the starting values come from 30 degrees
  {
    myInd <- which(startingDensities30deg$Adapted_Temperature_C == initialDensities$tAdaptation[iii] & 
                     is.na(str_match(tolower(startingDensities30deg$Replicate_label), initialDensities$line[iii])) == FALSE &
                     startingDensities30deg$Resource_Level_percent == initialDensities$mediumConcentration[iii])
    if (length(myInd) > 0)
    {
      initialDensities$Particle_manual_count[iii] <- startingDensities30deg$initial_density_cells_ml[myInd]
      initialDensities$incubationDuration[iii] <- 0
      initialDensities$volMl[iii] <- 1
    }
  } else {
    myInd <- which(startingDensities$Adapted_Temperature_C == initialDensities$tAdaptation[iii] & 
                     is.na(str_match(tolower(startingDensities$Replicate_label), initialDensities$line[iii])) == FALSE &
                     startingDensities$Resource_Level_percent == initialDensities$mediumConcentration[iii])
    if (length(myInd) > 0)
    {
      initialDensities$Particle_manual_count[iii] <- startingDensities$initial_density_cells_ml[myInd]
      initialDensities$incubationDuration[iii] <- 0
      initialDensities$volMl <- 1
    }
  }
}

allExperimentResults <- rbind(allExperimentResults, initialDensities)

# Now save the population growth table with all the processed additional information, standardized format etc.
write.table(allExperimentResults, file="Tetrahymena_pop_growth_complete.csv", append = FALSE, sep = ", ", row.names = FALSE)

allTAdaptation <- sort(unique(allExperimentResults$tAdaptation))

# Initialise a data frame for the results
populationGrowthRates <- data.frame(line=character(), tAdaptation=double(), mediumConcentration=double(), tTest = double(), growthRate = double(), standardErrorOnSlope = double(), fitIntercept = double())

iteration = 0
iterationCombinedLines = 0
plotList = list() # To save the plots in a list
plotList2 = list() # there are four plots for each condition and I save them all
plotList3 = list()
plotList4 = list()
plotListCombined = list()
plotCounter1 <- 1
plotCounter2 <- 1
plotCounter3 <- 1
plotCounter4 <- 1
plotCounterCombined <- 1
for (aaa in 1:length(allTAdaptation))
{
  allMediumConcentrations <- sort(unique(allExperimentResults$mediumConcentration[allExperimentResults$tAdaptation == allTAdaptation[aaa]]))
  for (mmm in 1:length(allMediumConcentrations))
  {
    allLines <- sort(unique(allExperimentResults$line[allExperimentResults$tAdaptation == allTAdaptation[aaa] & allExperimentResults$mediumConcentration == allMediumConcentrations[mmm]]), decreasing=TRUE)
    iterationCombinedLines = iterationCombinedLines + 1
    for (lll in 1:length(allLines))
    { 
      
      if (skipMotherCulture & is.finite(as.numeric(allLines[lll])))
      {
        next
      }
      
      iteration = iteration + 1
      
      tempDataset <- subset(allExperimentResults, line == allLines[lll] & tAdaptation == allTAdaptation[aaa] & mediumConcentration == allMediumConcentrations[mmm])
      
      
      
      if (dim(tempDataset)[1] == 0)
      {
        print("check as this should never happen!")
        next
      }
      print(paste("Population growth line:", allLines[lll], ", conc:", allMediumConcentrations[mmm], "%, T adapt.:", allTAdaptation[aaa], "C", sep=""))
      
      tempDataset$log2Density <- log2(tempDataset$Particle_manual_count/tempDataset$volMl)
      tempDataset <- tempDataset[is.finite(tempDataset$log2Density), ]
      myFit <- tempDataset %>% group_by(tTest) %>% do(tidy(lm(formula = log2Density ~ incubationDuration, data=.)))
      myFitIntercepts <- myFit[myFit$term == "(Intercept)",]
      myFitSlopes <- myFit[myFit$term == "incubationDuration",]
      
      # maxDensity
      maxDensity <- tempDataset %>% group_by(tTest) %>% dplyr::summarise(log2DensityMax = max(log2Density))
      # This can give us an idea of the maximum number of generations in culture
      # mean or max(populationGrowthRates$maxDensity - populationGrowthRates$fitIntercept)
      
      for (iii in 1:dim(myFitSlopes)[1])
      {
        growthRateInfo <- data.frame(line=allLines[lll], tAdaptation=allTAdaptation[aaa], mediumConcentration=allMediumConcentrations[mmm], tTest = myFitSlopes$tTest[iii], growthRate = myFitSlopes$estimate[iii], standardErrorOnSlope = myFitSlopes$std.error[lll], fitIntercept = myFitIntercepts$estimate[iii], maxDensity=maxDensity$log2DensityMax[iii])
        populationGrowthRates <- rbind(populationGrowthRates, growthRateInfo)
      }
      
      
      
      
      print(
         ggplot(tempDataset, aes(x=incubationDuration, y=log10(Particle_manual_count/volMl), colour=as.factor(tTest))) +
          geom_point(size=4) + 
          # geom_jitter() + 
          geom_smooth(method=lm, formula=' y ~ x', se=FALSE, fullrange=TRUE) + 
          theme_classic(base_size = 18) +
          scale_x_continuous(name="time (hours)") +
          scale_y_continuous(name="log10 density (cell/ml)") +
          scale_color_manual(values=c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")) + 
          # theme(legend.position = "none") + 
          labs(colour = "Growth\ntemperature:") +
          ggtitle(paste("Population growth\nline:", allLines[lll], ", conc:", allMediumConcentrations[mmm], "%, T adapt.:", allTAdaptation[aaa], "C", sep=""))
      )
      if (saveFigures){
        ggsave(file=paste("aaaa_short_term_population_growth_line_", allLines[lll], "_conc_", allMediumConcentrations[mmm], "_tadapt_", allTAdaptation[aaa], "_experiment.eps", sep=""), device="eps", dpi = 1200, width = 14, height = 10, units = "cm")
        ggsave(file=paste("aaaa_short_term_population_growth_line_", allLines[lll], "_conc_", allMediumConcentrations[mmm], "_tadapt_", allTAdaptation[aaa], "_experiment.png", sep=""), device="png", dpi = 600, width = 14, height = 10, units = "cm")
        
      }
      
      
      
      
      
      
      # FITTING SHARPE-SCHOOLFIELD
      
      currentCondition <-subset(populationGrowthRates, tAdaptation == allTAdaptation[aaa] & mediumConcentration == allMediumConcentrations[mmm] & line == allLines[lll])
      
      
      if(!require(dplyr)){install.packages('dplyr')}
      d <- currentCondition %>%
        select("line", "tAdaptation", "growthRate", "tTest")%>%
        dplyr::rename(
          curve_id = line,
          growth_temp = tAdaptation,
          rate = growthRate,
          temp = tTest
        )
      
      # add one element with zero growth to avoid 
      # the failing of the fitting with the deactivation part
      extraTemp <- 40
      d <- rbind(d, data.frame(curve_id = allLines[lll], growth_temp = allTAdaptation[aaa], rate = 0, temp = extraTemp))
      
      
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
      
      # # If the fit failed
      # if (is.null(coef(fit)))
      # {
      #   calculatedFitParameters <- cbind(calculatedFitParameters, data.frame(r_tref=NA, e_from_fit=NA, eh_from_fit=NA, th=NA))
      # }
      # else
      # {
      calculatedFitParameters <- cbind(calculatedFitParameters, data.frame(r_tref=coef(fit)[1], e_from_fit=coef(fit)[2], eh_from_fit=coef(fit)[3], th=coef(fit)[4]))
      #}
      
      
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
      
      
      
      # plot data and model fit
      lineColours = c("#3B9AB2", "#EBCC2A", "#F21A00", "#FF00FF")
      markerColours = c("#A3A3A3", "#666666", "#000000", "#FF00FF")
      markerShapes = c(16,17,15,4)
      
      
      # plot data and model fit
      thisPlot <- ggplot(d, aes(temp, rate*24)) +
        # geom_jitter(position = position_jitter(height = 0, width = .5), col=markerColours[mmm], shape=markerShapes[mmm]) +
        annotate("text", size=5, x=10.5, y=3, label= paste("E=", round(calculatedFitParameters$e_from_fit, 2), "eV", sep=""), hjust = 0, parse=F) +
        geom_line(aes(temp, .fitted*24), preds, col = lineColours[aaa], size=2) +
        geom_point(col=markerColours[mmm], shape=markerShapes[mmm]) +
        theme_classic(base_size = 15) +
        scale_x_continuous(name="Temp. (°C)",   limits=c(10, 32.5)) +
        scale_y_continuous(name="Growth (gen/day)", limits=c(0,5), breaks=seq(0,5,by=1)) # +
      # ggtitle(paste("T:", allTAdaptation[aaa], "; C:", allMediumConcentrations[mmm], "; L: ", allLines[lll]))
      thisPlot
      if (saveFigures){
        currentTitle <- paste("aaaa_short_term_population_growth_line_", allLines[lll], "_conc_", allMediumConcentrations[mmm], "_tadapt_", allTAdaptation[aaa], sep="")
        ggsave(file=paste(currentTitle, "_schoolfield.png", sep=""), dpi = 600, width = 12, height = 10, units = "cm")
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
      
      
    }
    
    
    # Repeat the fit on the combined lines:
    
    
    
    # FITTING SHARPE-SCHOOLFIELD
    
    currentCondition <-subset(populationGrowthRates, tAdaptation == allTAdaptation[aaa] & mediumConcentration == allMediumConcentrations[mmm])
    
    
    if(!require(dplyr)){install.packages('dplyr')}
    d <- currentCondition %>%
      select("line", "tAdaptation", "growthRate", "tTest")%>%
      dplyr::rename(
        curve_id = line,
        growth_temp = tAdaptation,
        rate = growthRate,
        temp = tTest
      )
    
    # add one element with zero growth to avoid 
    # the failing of the fitting with the deactivation part
    extraTemp <- 40
    d <- rbind(d, data.frame(curve_id = "all", growth_temp = allTAdaptation[aaa], rate = 0, temp = extraTemp))
    
    
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
    
    # # If the fit failed
    # if (is.null(coef(fit)))
    # {
    #   calculatedFitParameters <- cbind(calculatedFitParameters, data.frame(r_tref=NA, e_from_fit=NA, eh_from_fit=NA, th=NA))
    # }
    # else
    # {
    calculatedFitParameters <- cbind(calculatedFitParameters, data.frame(r_tref=coef(fit)[1], e_from_fit=coef(fit)[2], eh_from_fit=coef(fit)[3], th=coef(fit)[4]))
    #}
    
    
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
    
    
    
    
    if (iterationCombinedLines == 1)
    {
      allFitResultsLinesCombined <- calculatedFitParameters[0,]
    }
    
    allFitResultsLinesCombined <- rbind(allFitResultsLinesCombined, calculatedFitParameters)
    
    # predict new data
    new_data <- data.frame(temp = seq(min(d$temp), max(d$temp), 0.5))
    preds <- augment(fit, newdata = new_data)
    
    
    
    # plot data and model fit
    lineColours = c("#3B9AB2", "#EBCC2A", "#F21A00", "#FF00FF")
    markerColours = c("#A3A3A3", "#666666", "#000000", "#FF00FF")
    markerShapes = c(16,17,15,4)
    
    
    # plot data and model fit
    thisPlot <- ggplot(d, aes(temp, rate*24)) +
      # geom_jitter(position = position_jitter(height = 0, width = .5), col=markerColours[mmm], shape=markerShapes[mmm]) +
      annotate("text", size=5, x=10.5, y=3, label= paste("E=", round(calculatedFitParameters$e_from_fit, 2), "eV", sep=""), hjust = 0, parse=F) +
      geom_line(aes(temp, .fitted*24), preds, col = lineColours[aaa], size=2) +
      geom_point(col=markerColours[mmm], shape=markerShapes[mmm]) +
      theme_classic(base_size = 15) +
      scale_x_continuous(name="Temp. (°C)", limits=c(10, 32.5)) +
      scale_y_continuous(name="Growth (gen/day)", limits=c(0,5), breaks=seq(0,5,by=1)) # +
    # ggtitle(paste("T:", allTAdaptation[aaa], "; C:", allMediumConcentrations[mmm], "; L: ", allLines[lll]))
    thisPlot
    if (saveFigures){
      currentTitle <- paste("aaaa_short_term_population_growth_line_all_conc_", allMediumConcentrations[mmm], "_tadapt_", allTAdaptation[aaa], sep="")
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
    thisPlot <- ggplot(d2, aes(temp, rate*24)) +
      geom_line(aes(temp, .fitted*24), preds, col = lineColours[aaa], size=2) +
      geom_point(size=3, col=markerColours[mmm], shape=markerShapes[mmm]) +
      annotate("text", size=5, x=10.5, y=3, label= paste("E=", round(calculatedFitParameters$e_from_fit, 2), "eV", sep=""), hjust = 0, parse=F) +
      geom_errorbar(aes(ymin=(rate-sd)*24, ymax=(rate+sd)*24), width=1.3, col=markerColours[mmm]) + 
      theme_classic(base_size = 14) +
      scale_x_continuous(name="Temp. (°C)", limits=c(10, 32.5)) +
      scale_y_continuous(name="Growth (gen/day)", limits=c(0,5), breaks=seq(0,5,by=1)) # +
    
    thisPlot
    if (saveFigures){
      ggsave(file=paste(currentTitle, "_schoolfield_errorbar.png", sep=""), dpi = 600, width = 12, height = 10, units = "cm")
    }
    
    plotListCombined[[plotCounterCombined]] <- thisPlot
    plotCounterCombined <- plotCounterCombined + 1
    
    
    if (includeBootstrap)
    {
      # plot bootstrapped CIs
      thisPlot <- ggplot(d2, aes(temp, rate*24))
      if (!is.null(fit))
      { thisPlot <- thisPlot + 
        geom_line(aes(temp, .fitted*24), preds, col = lineColours[aaa], size=2)
      }
      thisPlot <- thisPlot +
        geom_line(aes(temp, pred*24, group = iter), boot1_preds, col = lineColours[aaa], alpha = 0.007) +
        geom_point(size=3, col=markerColours[mmm], shape=markerShapes[mmm]) +
        annotate("text", size=5, x=10.5, y=3, label= paste("E=", round(calculatedFitParameters$e_from_fit, 2), "eV", sep=""), hjust = 0, parse=F) +
        geom_errorbar(aes(ymin=(rate-sd)*24, ymax=(rate+sd)*24), width=1.3, col=markerColours[mmm]) + 
        theme_classic(base_size = 14) +
        scale_x_continuous(name="Temp. (°C)", limits=c(10, 32.5)) +
        scale_y_continuous(name="Growth (gen/day)", limits=c(0,5), breaks=seq(0,5,by=1)) # +
      thisPlot
      if (saveFigures){
        ggsave(file=paste(currentTitle, "_schoolfield_errorbar_and_bootstrap.png", sep=""), dpi = 600, width = 12, height = 10, units = "cm")
      }
      
      # plot bootstrapped CIs
      thisPlot <- ggplot()
      if (!is.null(fit))
      { thisPlot <- thisPlot + 
        geom_line(aes(temp, .fitted*24), preds, col = lineColours[aaa], size=2)
      }
      thisPlot <- thisPlot +
        geom_ribbon(aes(temp, ymin = conf_lower*24, ymax = conf_upper*24), boot1_conf_preds, fill = lineColours[aaa], alpha = 0.3) +
        geom_point(data=d2, aes(temp, rate*24), size=3, col=markerColours[mmm], shape=markerShapes[mmm]) +
        annotate("text", size=5, x=10.5, y=3, label= paste("E=", round(calculatedFitParameters$e_from_fit, 2), "eV", sep=""), hjust = 0, parse=F) +
        geom_errorbar(data=d2, aes(x=temp, ymin=(rate-sd)*24, ymax=(rate+sd)*24), width=1.3, col=markerColours[mmm]) + 
        theme_classic(base_size = 14) +
        scale_x_continuous(name="Temp. (°C)", limits=c(10, 32.5)) +
        scale_y_continuous(name="Growth (gen/day)", limits=c(0,5), breaks=seq(0,5,by=1)) # +
      thisPlot
      if (saveFigures){
        ggsave(file=paste(currentTitle, "_schoolfield_errorbar_and_bootstrap2.png", sep=""), dpi = 600, width = 12, height = 10, units = "cm")
      }
      
      
    }
    
  }
}



print(allFitResults)


write.table(allFitResults, "allFitResults_Schoolfield_on_short_term_population_growth.csv", append = FALSE, sep = ", ", row.names = FALSE)
write.table(allFitResultsLinesCombined, "allFitResults_Schoolfield_on_short_term_population_growth_all.csv", append = FALSE, sep = ", ", row.names = FALSE)


# Make a figure with all the plots
library(ggpubr)
# combine together the two plots into a single figure
fullFigure <- ggarrange(plotListCombined[[1]], plotListCombined[[2]], plotListCombined[[3]],
                        plotListCombined[[4]], plotListCombined[[5]], plotListCombined[[6]],
                        plotListCombined[[7]], plotListCombined[[8]], plotListCombined[[9]],
                        # labels = c("15-50", "15-100", "15-200", "20-50", "20-100", "20-200", "25-50", "25-100", "25-200"),
                        ncol = 3, nrow = 3, common.legend=TRUE)

# annotate_figure(fullFigure, 
#                 top = "50%                    100%                    200%",
#                 left = "15°C                   20°C                   25°C")
fullFigure

if(saveFigures)
{
  ggsave(file="figure_population_growth_thermal_response_post_adaptation_all.png", dpi = 600, width = 24, height = 20, units = "cm")
  ggsave(file="figure_population_growth_thermal_response_post_adaptation_all.eps", device="eps", dpi = 1200, width = 24, height = 20, units = "cm")
  library(Cairo)
  ggsave(file="figure_population_growth_thermal_response_post_adaptation_all.pdf", device=cairo_pdf, dpi = 1200, width = 24, height = 20, units = "cm")
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

if(saveFigures)
{
  ggsave(file="figure_population_growth_thermal_response_post_adaptation.png", dpi = 600, width = 24, height = 20, units = "cm")
  ggsave(file="figure_population_growth_thermal_response_post_adaptation.eps", device="eps", dpi = 1200, width = 24, height = 20, units = "cm")
  library(Cairo)
  ggsave(file="figure_population_growth_thermal_response_post_adaptation.pdf", device=cairo_pdf, dpi = 1200, width = 24, height = 20, units = "cm")
}

# combine together the two plots into a single figure
fullFigure <- ggarrange(plotList2[[1]], plotList2[[2]], plotList2[[3]],
                        plotList2[[4]], plotList2[[5]], plotList2[[6]],
                        plotList2[[7]], plotList2[[8]], plotList2[[9]],
                        # labels = c("15-50", "15-100", "15-200", "20-50", "20-100", "20-200", "25-50", "25-100", "25-200"),
                        ncol = 3, nrow = 3, common.legend=TRUE)

# annotate_figure(fullFigure, 
#                 top = "50%                    100%                    200%",
#                 left = "15°C                   20°C                   25°C")
fullFigure

if(saveFigures)
{
  ggsave(file="figure_population_growth_thermal_response_post_adaptation2.png", dpi = 600, width = 24, height = 20, units = "cm")
  ggsave(file="figure_population_growth_thermal_response_post_adaptation2.eps", device="eps", dpi = 1200, width = 24, height = 20, units = "cm")
  library(Cairo)
  ggsave(file="figure_population_growth_thermal_response_post_adaptation2.pdf", device=cairo_pdf, dpi = 1200, width = 24, height = 20, units = "cm")
}

# combine together the two plots into a single figure
fullFigure <- ggarrange(plotList3[[1]], plotList3[[2]], plotList3[[3]],
                        plotList3[[4]], plotList3[[5]], plotList3[[6]],
                        plotList3[[7]], plotList3[[8]], plotList3[[9]],
                        # labels = c("15-50", "15-100", "15-200", "20-50", "20-100", "20-200", "25-50", "25-100", "25-200"),
                        ncol = 3, nrow = 3, common.legend=TRUE)

# annotate_figure(fullFigure, 
#                 top = "50%                    100%                    200%",
#                 left = "15°C                   20°C                   25°C")
fullFigure

if(saveFigures)
{
  ggsave(file="figure_population_growth_thermal_response_post_adaptation3.png", dpi = 600, width = 24, height = 20, units = "cm")
  ggsave(file="figure_population_growth_thermal_response_post_adaptation3.eps", device="eps", dpi = 1200, width = 24, height = 20, units = "cm")
  library(Cairo)
  ggsave(file="figure_population_growth_thermal_response_post_adaptation3.pdf", device=cairo_pdf, dpi = 1200, width = 24, height = 20, units = "cm")
}

# combine together the two plots into a single figure
fullFigure <- ggarrange(plotList4[[1]], plotList4[[2]], plotList4[[3]],
                        plotList4[[4]], plotList4[[5]], plotList4[[6]],
                        plotList4[[7]], plotList4[[8]], plotList4[[9]],
                        # labels = c("15-50", "15-100", "15-200", "20-50", "20-100", "20-200", "25-50", "25-100", "25-200"),
                        ncol = 3, nrow = 3, common.legend=TRUE)

# annotate_figure(fullFigure, 
#                 top = "50%                    100%                    200%",
#                 left = "15°C                   20°C                   25°C")
fullFigure

if(saveFigures)
{
  ggsave(file="figure_population_growth_thermal_response_post_adaptation4.png", dpi = 600, width = 24, height = 20, units = "cm")
  ggsave(file="figure_population_growth_thermal_response_post_adaptation4.eps", device="eps", dpi = 1200, width = 24, height = 20, units = "cm")
  library(Cairo)
  ggsave(file="figure_population_growth_thermal_response_post_adaptation4.pdf", device=cairo_pdf, dpi = 1200, width = 24, height = 20, units = "cm")
}



ggplot(populationGrowthRates, aes(x=tTest, y=growthRate*24)) +
  geom_point(size=4, shape="+", aes(colour=as.factor(line))) + 
  # geom_jitter() + 
  geom_smooth(method=lm, formula=' y ~ x', se=FALSE, fullrange=TRUE) + 
  theme_classic(base_size = 18) +
  scale_x_continuous(name="temperature (°C)") +
  scale_y_continuous(name="Growth (gen/day)") +
  scale_color_manual(values=c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")) + 
  # theme(legend.position = "none") + 
  labs(colour = "Line:") +
  ggtitle("Population growth") +
  facet_grid(mediumConcentration ~ tAdaptation)
if (saveFigures){
  ggsave(file="aaaa_short_term_population_growth_experiment.eps", device="eps", dpi = 1200, width = 14, height = 10, units = "cm")
}


populationGrowthRates$tNormalized = 1 / populationGrowthRates$tTest - 1/(20+273.15)
populationGrowthRates$logGrowth = log(populationGrowthRates$growthRate)





# plot data and model fit
ggplot(allFitResults, aes(x=tAdaptation, y=e_from_fit, color=factor(mediumConcentration), shape=factor(mediumConcentration), label=factor(line))) +
  # geom_point(size=6) +
  # geom_jitter(size=3, width=0.4) +
  geom_jitter(size=4, position = position_jitter(height = 0, width = .3)) + 
  # geom_text(color="black", size=4) + 
  scale_color_manual(values=c("#A3A3A3", "#666666", "#000000")) +
  scale_fill_manual(values=c("#A3A3A3", "#666666", "#000000")) +
  scale_y_continuous(name="activation energy (eV)", limits=c(0, 1.5)) +
  scale_x_continuous(name="Adaptation temperature (ºC)", breaks=c(15,20,25), limits=c(12.5, 27.5)) +
  theme_classic(base_size=15) +
  theme(legend.position = "none") +
  labs(color = "Conc. (%)") # + # this specifies a custom legend
# ggtitle("Activation energy on population growth")

if (saveFigures){
  ggsave(file="figure_activation_energy_on_short_term_population_growth.png", dpi = 600, width = 12, height = 14, units = "cm")
  ggsave(file="figure_activation_energy_on_short_term_population_growth.eps", device="eps", dpi = 1200, width = 12, height = 14, units = "cm")
  library(Cairo)
  ggsave(file="figure_activation_energy_on_short_term_population_growth.pdf", device=cairo_pdf, dpi = 1200, width = 12, height = 14, units = "cm")
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
  ggsave(file="aaaa_optimal_temp_vs_adaptation_temp_on_short_term_pop_growth.png", dpi = 600, width = 12, height = 10, units = "cm")
}




# plot data and model fit
referenceTemperature <- 20
ggplot(allFitResults, aes(x=tAdaptation, y=r_tref*24, color=factor(mediumConcentration), shape=factor(mediumConcentration), label=factor(line))) +
  # geom_point(size=4) +
  geom_jitter(size=4, position = position_jitter(height = 0, width = .3)) + 
  # geom_text(color="black", size=4) + 
  scale_color_manual(values=c("#A3A3A3", "#666666", "#000000")) +
  scale_fill_manual(values=c("#A3A3A3", "#666666", "#000000")) +
  scale_y_continuous(name="Growth (gen/day)") +
  scale_x_continuous(name=expression("Adaptation temperature (ºC)"), limits=c(12.5, 27.5), breaks=c(15, 20, 25)) +
  theme_classic(base_size=12) +
  theme(legend.position = "none") +
  labs(color = "Conc.") + # this specifies a custom legend
  ggtitle(paste("Population growth at reference temperature (", referenceTemperature, "°C)", sep=""))


if (saveFigures)
{
  ggsave(file="figure_short_term_population_growth_at_tRef.png", dpi = 600, width = 12, height = 10, units = "cm")
  ggsave(file="figure_short_term_population_growth_at_tRef.eps", device="eps", dpi = 1200, width = 12, height = 10, units = "cm")
  library(Cairo)
  ggsave(file="figure_short_term_population_growth_at_tRef.pdf", device=cairo_pdf, dpi = 1200, width = 12, height = 10, units = "cm")
}







# plot data and model fit
ggplot(allFitResultsLinesCombined, aes(x=tAdaptation, y=e_from_fit, color=factor(mediumConcentration), shape=factor(mediumConcentration), label=factor(line))) +
  # geom_point(size=6) +
  # geom_jitter(size=3, width=0.4) +
  geom_jitter(size=4, position = position_jitter(height = 0, width = .3)) + 
  # geom_text(color="black", size=4) + 
  scale_color_manual(values=c("#A3A3A3", "#666666", "#000000")) +
  scale_fill_manual(values=c("#A3A3A3", "#666666", "#000000")) +
  scale_y_continuous(name="activation energy (eV)", limits=c(0, 1.5)) +
  scale_x_continuous(name="Adaptation temperature (ºC)", breaks=c(15,20,25), limits=c(12.5, 27.5)) +
  theme_classic(base_size=15) +
  theme(legend.position = "none") +
  labs(color = "Conc. (%)") # + # this specifies a custom legend
# ggtitle("Activation energy on population growth")

if (saveFigures){
  ggsave(file="figure_activation_energy_on_short_term_population_growth_all.png", dpi = 600, width = 12, height = 14, units = "cm")
  ggsave(file="figure_activation_energy_on_short_term_population_growth_all.eps", device="eps", dpi = 1200, width = 12, height = 14, units = "cm")
  library(Cairo)
  ggsave(file="figure_activation_energy_on_short_term_population_growth_all.pdf", device=cairo_pdf, dpi = 1200, width = 12, height = 14, units = "cm")
}





# plot data and model fit
ggplot(allFitResultsLinesCombined, aes(x=tAdaptation, y=topt, color=factor(mediumConcentration), shape=factor(mediumConcentration), label=factor(line))) +
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
  ggsave(file="aaaa_optimal_temp_vs_adaptation_temp_on_short_term_pop_growth_all.png", dpi = 600, width = 12, height = 10, units = "cm")
}




# plot data and model fit
referenceTemperature <- 20
ggplot(allFitResultsLinesCombined, aes(x=tAdaptation, y=r_tref*24, color=factor(mediumConcentration), shape=factor(mediumConcentration), label=factor(line))) +
  # geom_point(size=4) +
  geom_jitter(size=4, position = position_jitter(height = 0, width = .3)) + 
  # geom_text(color="black", size=4) + 
  scale_color_manual(values=c("#A3A3A3", "#666666", "#000000")) +
  scale_fill_manual(values=c("#A3A3A3", "#666666", "#000000")) +
  scale_y_continuous(name="Growth (gen/day)") +
  scale_x_continuous(name=expression("Adaptation temperature (ºC)"), limits=c(12.5, 27.5), breaks=c(15, 20, 25)) +
  theme_classic(base_size=12) +
  theme(legend.position = "none") +
  labs(color = "Conc.") + # this specifies a custom legend
  ggtitle(paste("Population growth at reference temperature (", referenceTemperature, "°C)", sep=""))


if (saveFigures)
{
  ggsave(file="figure_short_term_population_growth_at_tRef_all.png", dpi = 600, width = 12, height = 10, units = "cm")
  ggsave(file="figure_short_term_population_growth_at_tRef_all.eps", device="eps", dpi = 1200, width = 12, height = 10, units = "cm")
  library(Cairo)
  ggsave(file="figure_short_term_population_growth_at_tRef_all.pdf", device=cairo_pdf, dpi = 1200, width = 12, height = 10, units = "cm")
}


if (includeBootstrap)
{
  # plot data and model fit
  ggplot(allFitResultsLinesCombined, aes(x=tAdaptation, y=e_from_fit, color=factor(mediumConcentration), fill=factor(mediumConcentration + tAdaptation*100), shape=factor(mediumConcentration), label=factor(line))) +
    geom_errorbar(aes(ymin=e_from_fit - sd_e_from_fit, ymax=e_from_fit + sd_e_from_fit), width=0, position=position_dodge(width=2)) + 
    geom_point(size=6, position=position_dodge(width=2)) + 
    # geom_errorbar(aes(ymin=e_from_fit_0025, ymax=e_from_fit_0975), width=0, position=position_dodge(width=2)) + 
    scale_color_manual(values=c("#A3A3A3", "#666666", "#000000")) +
    scale_fill_manual(values=c("#8FCDDC", "#3B9AB2", "#18434E", "#F2DD70", "#EBCC2A", "#463C07", "#FF7C6C", "#F21A00", "#660B00")) +
    scale_shape_manual(values=c(21, 24, 22)) + # shapes for the markers
    scale_y_continuous(name="activation energy (eV)") +
    scale_x_continuous(name="Adaptation temperature (ºC)", breaks=c(15,20,25), limits=c(12.5, 27.5)) +
    coord_cartesian(ylim = c(0, 1)) +
    theme_classic(base_size=18) +
    theme(legend.position = "none") +
    labs(color = "Conc. (%)") # + # this specifies a custom legend
  # ggtitle("Activation energy on population growth")
  
  if (saveFigures){
    ggsave(file="figure_activation_energy_on_short_term_population_growth_all_bootstrap.png", dpi = 600, width = 12, height = 10, units = "cm")
    ggsave(file="figure_activation_energy_on_short_term_population_growth_all_bootstrap.eps", device="eps", dpi = 1200, width = 12, height = 10, units = "cm")
    library(Cairo)
    ggsave(file="figure_activation_energy_on_short_term_population_growth_all_bootstrap.pdf", device=cairo_pdf, dpi = 1200, width = 12, height = 10, units = "cm")
  }
  
  
  
  
  
  # plot data and model fit
  referenceTemperature <- 20
  ggplot(allFitResultsLinesCombined, aes(x=tAdaptation, y=r_tref*24, color=factor(mediumConcentration), fill=factor(mediumConcentration + tAdaptation*100), shape=factor(mediumConcentration), label=factor(line))) +
    geom_errorbar(aes(ymin= (r_tref - sd_r_tref)*24, ymax=(r_tref + sd_r_tref)*24), width=0, position=position_dodge(width=2)) + 
    geom_point(size=6, position=position_dodge(width=2)) + 
    scale_color_manual(values=c("#A3A3A3", "#666666", "#000000")) +
    scale_fill_manual(values=c("#8FCDDC", "#3B9AB2", "#18434E", "#F2DD70", "#EBCC2A", "#463C07", "#FF7C6C", "#F21A00", "#660B00")) +
    scale_shape_manual(values=c(21, 24, 22)) + # shapes for the markers
    scale_y_continuous(name="Growth at T=20°C (gen/day)") +
    scale_x_continuous(name=expression("Adaptation temperature (ºC)"), limits=c(12.5, 27.5), breaks=c(15, 20, 25)) +
    theme_classic(base_size=18) +
    theme(legend.position = "none") +
    labs(color = "Conc.")  # this specifies a custom legend
    # ggtitle(paste("Population growth at reference temperature (", referenceTemperature, "°C)", sep=""))
  
  
  if (saveFigures)
  {
    ggsave(file="figure_short_term_population_growth_at_tRef_all_bootstrap.png", dpi = 600, width = 12, height = 10, units = "cm")
    ggsave(file="figure_short_term_population_growth_at_tRef_all_bootstrap.eps", device="eps", dpi = 1200, width = 12, height = 10, units = "cm")
    library(Cairo)
    ggsave(file="figure_short_term_population_growth_at_tRef_all_bootstrap.pdf", device=cairo_pdf, dpi = 1200, width = 12, height = 10, units = "cm")
  }
}








lineColours = c("#3B9AB2", "#EBCC2A", "#F21A00", "#FF00FF")
lineTypes <- c("twodash", "dashed", "solid", "longdash")
# Now try to plot the fit only
allFitData <- data.frame()
for (iii in 1:dim(allFitResultsLinesCombined)[1])
{
  predictedData <- data.frame(temp = seq(10, 40, 0.5))
  predictedData$rate <- sharpeschoolhigh_1981(temp = predictedData$temp, r_tref = allFitResultsLinesCombined$r_tref[iii],e = allFitResultsLinesCombined$e_from_fit[iii], eh = allFitResultsLinesCombined$eh_from_fit[iii], th = allFitResultsLinesCombined$th[iii], tref = 20)
  predictedData$tAdaptation <- allFitResultsLinesCombined$tAdaptation[iii]
  predictedData$mediumConcentration <- allFitResultsLinesCombined$mediumConcentration[iii]
  predictedData$line <- allFitResultsLinesCombined$line[iii]
  allFitData <- rbind(allFitData, predictedData)
}

plotG1 <- ggplot(allFitData, aes(x=temp, y=rate*24, color=as.factor(tAdaptation), size=as.factor(mediumConcentration), alpha=mediumConcentration)) + # linetype=as.factor(mediumConcentration), 
  geom_line() +
  scale_x_continuous(name="Temp. (°C)",  limits=c(10, 25)) +
  scale_y_continuous(name="Growth (gen/day)", limits=c(0.4, 4), trans = 'log10')  +
  theme_classic(base_size=18) +
  scale_color_manual(values=lineColours) +
  scale_linetype_manual(values=lineTypes) + 
  scale_alpha_continuous(range=c(0.5, 1)) +
  scale_size_manual( values = c(0.8, 1.4, 2) ) +
  theme(legend.position = "none") +
  labs(color = "Conc.") 
plotG1


if (saveFigures){
    ggsave(file="figure_fitted_growth_vs_temp_all.png", dpi = 600, width = 12, height = 10, units = "cm")
    ggsave(file="figure_fitted_growth_vs_temp_all.eps", device="eps", dpi = 1200, width = 12, height = 10, units = "cm")
    library(Cairo)
    ggsave(file="figure_fitted_growth_vs_temp_all.pdf", device=cairo_pdf, dpi = 1200, width = 12, height = 10, units = "cm")
}

