rm(list=ls()) # clean memory
# read the results again


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
saveFigures <- TRUE

# dev.off()

# load packages
library(ggplot2)
library("stringr")

# Now read the data from 2020
fileName <- "~/Tetrahymena/tetrahymena_density_speed/track_analysis_results_individual_particles.csv" # file.choose() # ask the user to select a file name.
setwd(dirname(fileName))


allExperimentResults <- read.table(file = fileName, sep = ",", header=TRUE, na.strings = c("NA", " NA"))

# Extract information from the string in the "fileName" field
# Extract information from the string in the "fileName" field
newColumns <- str_split_fixed(as.character(str_trim(tolower(allExperimentResults$fileName))), "tracking/|ctracking_diluted/|ctracking_dilute/|xdilute|erx|er|x|cr", 5)[,2:3]
# class(newColumns) <- "numeric"
newColumns <- data.frame(
  tAdaptation = as.numeric(newColumns[,2]),
  mediumConcentration = 100*as.numeric(newColumns[,1]),
  tTest = as.numeric(newColumns[,2]) )
# replace water with zero mediumConcentration (but careful that it is replacing anything that is na)
newColumns$mediumConcentration[is.na(newColumns$mediumConcentration)] = 0
# colnames(newColumns) <- c("tAdaptation", "tTest", "line")

# This below is to split data for which the medium was
# diluted with fresh medium at the same concentration
# just before measuring, so as to get lower density of
# tetrahymena
# newColumns$line <- "fDiluted"
# newColumns$line[which(is.na(str_match(allExperimentResults$fileName, "dilut")))] <- "fUndiluted"

newColumns$line <- "f2020"

# add the new columns to the original table
allExperimentResults <- cbind(allExperimentResults, newColumns, deparse.level = 1, stringsAsFactors = default.stringsAsFactors())


allExperimentResults$logSpeed <- log10(allExperimentResults$medianSpeed)

names(allExperimentResults)

allTAdaptation <- unique(allExperimentResults$tAdaptation)

for (aaa in 1:length(allTAdaptation))
{
  allTTest <-  unique(allExperimentResults$tTest[allExperimentResults$tAdaptation == allTAdaptation[aaa]])
  for (ttt in 1:length(allTTest))
  {
    allLines <- unique(allExperimentResults$line[allExperimentResults$tAdaptation == allTAdaptation[aaa] & allExperimentResults$tTest == allTTest[ttt]])
    for (lll in 1:length(allLines))
    {
      print(paste("TAdapt: ", allTAdaptation[aaa], "; tTest: ", allTTest[ttt]," Line: ", allLines[lll], sep=""))
      currentTitle <- paste("aaaaa_", allTAdaptation[aaa], "_", allTTest[ttt], "_", allLines[lll], sep="")
      
      
      
      currentCondition <- subset(allExperimentResults, tAdaptation == allTAdaptation[aaa] & tTest == allTTest[ttt] & line==allLines[lll])
      
      
      
      # If I further filter the data:
      # keep those that are big enough and that move
      currentCondition <- subset(currentCondition, medianArea > 400 & medianSpeed > 100)
      
      # Only keep those that move straight
      # currentCondition <- subset(currentCondition, meanAbsAngle*180/pi < 15)
      
      
      
      
      # Only keep a percentile for each medium concentration
      allMediumConcentrations <- unique(currentCondition$mediumConcentration)
      
      for (mmm in 1:length(allMediumConcentrations))
      {
        currentConditionThisConcentration <- subset(currentCondition, mediumConcentration == allMediumConcentrations[mmm])
        speedQuantile <- quantile(currentConditionThisConcentration$medianSpeed, probs = 0.80, na.rm = TRUE)
        currentCondition$speedQuantile[currentCondition$mediumConcentration == allMediumConcentrations[mmm]] <- speedQuantile
      }
      currentConditionTopSpeed <- subset(currentCondition, medianSpeed >= speedQuantile)
      
      ggplot(currentCondition, aes(x=mediumConcentration, y=medianSpeed, color=medianElongation)) +
        geom_jitter(position = position_jitter(height = 0, width = .5)) + 
        # geom_jitter(position = position_jitter(height = 0, width = .3)) +
        # geom_smooth(method="lm", se=TRUE, fill=NA, formula=y ~ poly(x, 3, raw=TRUE)) + 
        scale_y_continuous(name=expression(paste("median speed (", mu, "m/s)")), limits=c(0,1000)) +
        scale_x_continuous(name="Medium Concentration (%)") +
        theme_classic(base_size = 15) +
        # theme(legend.position = "none") + 
        ggtitle(paste("TAdapt:", allTAdaptation[aaa], "; tTest:", allTTest[ttt], "; Line: ", allLines[lll])) +
        scale_color_continuous(type = "viridis", option="magma", limits = c(1.4, 3.6)) +
        labs(color="Elongation:")
      # scale_color_manual(values=c("#3B9AB2", "#EBCC2A", "#F21A00")) # this is the zissou1 palette
      # if (saveFigures){
      #   ggsave(file=paste(currentTitle, "_speed_vs_concentration_jitter.png", sep=""), dpi = 600, width = 12, height = 10, units = "cm")
      # }
      
      
      
      
      ggplot(currentConditionTopSpeed, aes(x=mediumConcentration, y=medianSpeed, color=medianElongation)) +
        geom_jitter(position = position_jitter(height = 0, width = .5)) + 
        # geom_jitter(position = position_jitter(height = 0, width = .3)) +
        # geom_smooth(method="lm", se=TRUE, fill=NA, formula=y ~ poly(x, 3, raw=TRUE)) + 
        scale_y_continuous(name=expression(paste("median speed (", mu, "m/s)")), limits=c(0,1000)) +
        scale_x_continuous(name="Medium Concentration (%)") +
        theme_classic(base_size = 15) +
        # theme(legend.position = "none") + 
        ggtitle(paste("TAdapt:", allTAdaptation[aaa], "; tTest:", allTTest[ttt], "; Line: ", allLines[lll])) +
        scale_color_continuous(type = "viridis", option="magma", limits = c(1.4, 3.6)) +
        labs(color="Elongation:")
      # scale_color_manual(values=c("#3B9AB2", "#EBCC2A", "#F21A00")) # this is the zissou1 palette
      # if (saveFigures){
      #   ggsave(file=paste(currentTitle, "_top_speed_vs_concentration_jitter.png", sep=""), dpi = 600, width = 12, height = 10, units = "cm")
      # }
      
      
      
      # Summarize the data in order to plot the data with error bar
      d <- data_summary(currentConditionTopSpeed, varname="medianSpeed", 
                         groupnames=c("mediumConcentration", "tAdaptation", "tTest"))
      # head(d)
      
      markerColours = c("#0022AA", "#3B9AB2", "#EBCC2A", "#F21A00", "#FF44AA")
      markerShapes = c(23, 21, 24, 22, 25) # shapes for the markers
      
      # plot data and model fit
      thisPlot <- ggplot(d, aes(x=mediumConcentration, y=medianSpeed)) +
        geom_errorbar(aes(ymin=medianSpeed-sd, ymax=medianSpeed+sd), width=1.3, col="black") + 
        geom_point(size=6, fill=markerColours[aaa], color="black", shape=markerShapes[aaa]) +
        annotate("text", size=6, x=140, y=900, label= paste("T=", allTAdaptation[aaa], "°C", sep=""), hjust = 0, parse=F) +
        theme_classic(base_size = 18) +
        scale_x_continuous(name="Medium Concentration (%)",  limits=c(0, 200)) +
        scale_y_continuous(name=expression(paste("Speed (", mu, "m/s", ")")), limits=c(0,1000), breaks=seq(0, 1000, by=200))

      # print(thisPlot)
      thisPlot
      
      if (saveFigures){
        ggsave(file=paste(currentTitle, "_top_speed_vs_concentration_errorbar.png", sep=""), dpi = 600, width = 12, height = 10, units = "cm")
        ggsave(file=paste(currentTitle, "_top_speed_vs_concentration_errorbar.eps", sep=""), device="eps", dpi = 1200, width = 12, height = 10, units = "cm")
        library(Cairo)
        ggsave(file=paste(currentTitle, "_top_speed_vs_concentration_errorbar.pdf", sep=""), device=cairo_pdf, dpi = 1200, width = 12, height = 10, units = "cm")
     }
      

      
      
    }
    
    
  }
}

