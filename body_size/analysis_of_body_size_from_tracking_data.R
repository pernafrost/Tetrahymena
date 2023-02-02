rm(list=ls()) # clean memory

skipMotherCulture <- TRUE
saveFigures <- TRUE
dev.off()


library(ggplot2)

fileName <- "~/speed_response_2021/track_analysis_results_individual_particles.csv" # file.choose() # ask the user to select a file name.
workingDir <- "~/body_size"
setwd(workingDir)
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


names(allExperimentResults)

allTAdaptation <- sort(unique(allExperimentResults$tAdaptation))

iteration <- 0
for (aaa in 1:length(allTAdaptation))
{
  # print(paste("Adaptation temp.: ", allTAdaptation[aaa]))
  allMediumConcentrations <- sort(unique(allExperimentResults$mediumConcentration[allExperimentResults$tAdaptation == allTAdaptation[aaa]]))
  for (mmm in 1:length(allMediumConcentrations))
  {
    allLines <- sort(unique(allExperimentResults$line[allExperimentResults$tAdaptation == allTAdaptation[aaa] & allExperimentResults$mediumConcentration == allMediumConcentrations[mmm]]), decreasing=TRUE)
    for (lll in 1:length(allLines))
    {
      if (skipMotherCulture & is.finite(as.numeric(allLines[lll])))
      {
        next
      }
      
      print(paste("TAdapt:", allTAdaptation[aaa], "; C:", allMediumConcentrations[mmm], "; Line: ", allLines[lll]))
      currentTitle <- paste("aaaaa_", allTAdaptation[aaa], "_", allMediumConcentrations[mmm], "_", allLines[lll], sep="")
      
      
      iteration <- iteration + 1
      

      currentCondition <-subset(allExperimentResults, tAdaptation == allTAdaptation[aaa] & mediumConcentration == allMediumConcentrations[mmm] & line ==  allLines[lll])

      
      
      
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
      currentConditionTopSpeed <- subset(currentCondition, medianSpeed >= speedQuantile)
      
      
      
      
      
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
      
      currentConditionTopSpeed$estimatedlogVolume <- log10(4/3 * pi* currentConditionTopSpeed$medianBEllipse * currentConditionTopSpeed$medianAEllipse * currentConditionTopSpeed$medianAEllipse / 8)
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
      
      
      # currentConditionTopSpeedReasonableTemperature <- subset(currentConditionTopSpeed, tTest < 30) # tTest not too high as above 30 degrees cells change shape etc.
      currentConditionTopSpeedReasonableTemperature <- subset(currentConditionTopSpeed, abs(tTest - tAdaptation) <=5 ) # tTest not too different from adaptation temperature
      
      
      
      if (iteration == 1)
      {
        d <- currentConditionTopSpeedReasonableTemperature[0,]
      }
      
      d <- rbind(d, currentConditionTopSpeedReasonableTemperature)

      
      
      
      
    }
  }
}


currentTitle <- "figure_cell"

markerColours = c("#A3A3A3", "#666666", "#000000", "#FF00FF")
markerShapes = c(16,17,15,4)

ggplot(d, aes(x=factor(mediumConcentration), y=medianBEllipse, col=factor(tAdaptation))) +
  geom_boxplot(width=0.6, outlier.shape = NA, colour="black", aes(fill = factor(tAdaptation))) + # add a box plot
  theme_classic(base_size = 22) +
  theme(legend.position = "none") + 
  scale_y_continuous(name=expression(paste("Major axis length"," (",  mu, "m)")), lim=c(0,100)) + 
  scale_x_discrete(name="Adaptation conditions", labels=c("50%", "100%", "200%")) +
  scale_fill_manual(values= alpha(c("#3B9AB2", "#EBCC2A", "#F21A00", "#FF00FF")), 0.9) + # this is the zissou1 palette
  scale_color_manual(values=c("#3B9AB2", "#EBCC2A", "#F21A00", "#FF00FF")) # this is the zissou1 palette
if (saveFigures){
  ggsave(file=paste(currentTitle, "_length.png", sep=""), dpi = 600, width = 12, height = 10, units = "cm")
}


ggplot(d, aes(x=factor(mediumConcentration, levels=sort(unique(as.numeric(mediumConcentration)))), y=medianAEllipse, col=factor(tAdaptation))) +
  geom_boxplot(width=0.6, outlier.shape = NA, colour="black", aes(fill = factor(tAdaptation))) + # add a box plot
  theme_classic(base_size = 22) +
  theme(legend.position = "none") + 
  scale_y_continuous(name=expression(paste("Minor axis length"," (",  mu, "m)"))) +
  scale_x_discrete(name="Adaptation conditions", labels=c("50%", "100%", "200%")) +
  scale_fill_manual(values= alpha(c("#3B9AB2", "#EBCC2A", "#F21A00", "#FF00FF")), 0.9) + # this is the zissou1 palette
  scale_color_manual(values=c("#3B9AB2", "#EBCC2A", "#F21A00", "#FF00FF")) # this is the zissou1 palette
if (saveFigures){
  ggsave(file=paste(currentTitle, "_width.png", sep=""), dpi = 600, width = 12, height = 10, units = "cm")
}




ggplot(d, aes(x=factor(mediumConcentration, levels=sort(unique(as.numeric(mediumConcentration)))), y=medianElongation, col=factor(tAdaptation))) +
  geom_boxplot(width=0.6, outlier.shape = NA, colour="black", aes(fill = factor(tAdaptation))) + # add a box plot
  theme_classic(base_size = 22) +
  theme(legend.position = "none") + 
  scale_y_continuous(name=expression(paste("Aspect ratio"))) +
  scale_x_discrete(name="Adaptation conditions", labels=c("50%", "100%", "200%")) +
  scale_fill_manual(values= alpha(c("#3B9AB2", "#EBCC2A", "#F21A00", "#FF00FF")), 0.9) + # this is the zissou1 palette
  scale_color_manual(values=c("#3B9AB2", "#EBCC2A", "#F21A00", "#FF00FF")) # this is the zissou1 palette
if (saveFigures){
  ggsave(file=paste(currentTitle, "_elongation.png", sep=""), dpi = 600, width = 12, height = 10, units = "cm")
}


d$estimatedVolume <- 4/3 * pi * d$medianBEllipse * d$medianAEllipse * d$medianAEllipse / 8
ggplot(d, aes(x=factor(mediumConcentration, levels=sort(unique(as.numeric(mediumConcentration)))), y=estimatedVolume, col=factor(tAdaptation))) +
  geom_boxplot(width=0.6, outlier.shape = NA, colour="black", aes(fill = factor(tAdaptation))) + # add a box plot
  theme_classic(base_size = 22) +
  theme(legend.position = "none") + 
  scale_y_continuous(name=expression(paste('estimated volume '," ",  mu, 'm'^3))) +
  scale_x_discrete(name="Adaptation conditions", labels=c("50%", "100%", "200%")) +
  scale_fill_manual(values= alpha(c("#3B9AB2", "#EBCC2A", "#F21A00", "#FF00FF")), 0.9) + # this is the zissou1 palette
  scale_color_manual(values=c("#3B9AB2", "#EBCC2A", "#F21A00", "#FF00FF")) # this is the zissou1 palette
if (saveFigures){
  ggsave(file=paste(currentTitle, "_volume.png", sep=""), dpi = 600, width = 12, height = 10, units = "cm")
}


d$estimatedlogVolume <- log10(d$estimatedVolume)
ggplot(d, aes(x=factor(mediumConcentration, levels=sort(unique(as.numeric(mediumConcentration)))), y=estimatedlogVolume, col=factor(tAdaptation))) +
  # geom_violin(width=0.5, aes(fill = factor(tAdaptation))) +
  geom_boxplot(width=0.5, outlier.shape = NA, colour="black", aes(fill = factor(tAdaptation))) + # add a box plot
  theme_classic(base_size = 22) +
  theme(legend.position = "none") + 
  scale_y_continuous(name=expression(paste('log'[10]*'(volume)'," ",  mu, 'm'^3))) +
  scale_x_discrete(name="Adaptation conditions", labels=c("50%", "100%", "200%")) +
  scale_fill_manual(values= alpha(c("#3B9AB2", "#EBCC2A", "#F21A00", "#FF00FF")), 0.9) + # this is the zissou1 palette
  scale_color_manual(values=c("#3B9AB2", "#EBCC2A", "#F21A00", "#FF00FF")) # this is the zissou1 palette
if (saveFigures){
  ggsave(file=paste(currentTitle, "_logvolume.eps", sep=""), device="eps", dpi = 1200, width = 12, height = 10, units = "cm")
  library(Cairo)
  ggsave(file=paste(currentTitle, "_logvolume.pdf", sep=""), device=cairo_pdf, dpi = 1200, width = 12, height = 10, units = "cm")
  ggsave(file=paste(currentTitle, "_logvolume.png", sep=""), dpi = 600, width = 12, height = 10, units = "cm")
}



