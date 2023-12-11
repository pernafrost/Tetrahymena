rm(list=ls()) # clean memory

skipMotherCulture <- TRUE
saveFigures <- TRUE
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


library(ggplot2)

fileName <- "~/Tetrahymena/acute_speed_response/track_analysis_results_individual_particles.csv" # file.choose() # ask the user to select a file name.
workingDir <- "~/Tetrahymena/body_size"
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

# check if the data are for the mother culture
allExperimentResults$isMotherCulture <- (allExperimentResults$line == 1 | allExperimentResults$line == 2)


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
      # if (skipMotherCulture & is.finite(as.numeric(allLines[lll])))
      # {
      #   next
      # }
      
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

d1 <- d # make a copy of the data frame including both 
# mother culture and subcultures

if (skipMotherCulture)
{
  d <- subset(d, isMotherCulture == FALSE)
}


ggplot(d, aes(x=factor(mediumConcentration), y=medianBEllipse, col=factor(tAdaptation))) +
  geom_boxplot(width=0.6, outlier.shape = NA, colour="black", aes(fill = factor(tAdaptation))) + # add a box plot
  theme_classic(base_size = 22) +
  theme(legend.position = "none") + 
  scale_y_continuous(name=expression(paste("Major axis length"," (",  mu, "m)")), lim=c(0,100)) + 
  scale_x_discrete(name="Adaptation conditions", labels=c("50%", "100%", "200%")) +
  scale_fill_manual(values= alpha(c("#3B9AB2", "#EBCC2A", "#F21A00", "#FF00FF")), 0.9) + # this is the zissou1 palette
  scale_color_manual(values=c("#3B9AB2", "#EBCC2A", "#F21A00", "#FF00FF")) # this is the zissou1 palette
if (saveFigures){
  ggsave(file=paste(currentTitle, "_length.eps", sep=""), device="eps", dpi = 1200, width = 12, height = 10, units = "cm")
  ggsave(file=paste(currentTitle, "_length.png", sep=""), dpi = 600, width = 12, height = 10, units = "cm")
  library(Cairo)
  ggsave(file=paste(currentTitle, "_length.pdf", sep=""), device=cairo_pdf, dpi = 1200, width = 12, height = 10, units = "cm")
  
}


ggplot(d, aes(x=factor(mediumConcentration, levels=sort(unique(as.numeric(mediumConcentration)))), y=medianAEllipse, col=factor(tAdaptation))) +
  geom_boxplot(width=0.6, outlier.shape = NA, colour="black", aes(fill = factor(tAdaptation))) + # add a box plot
  theme_classic(base_size = 22) +
  theme(legend.position = "none") + 
  scale_y_continuous(name=expression(paste("Minor axis length"," (",  mu, "m)")), lim=c(0,50)) +
  scale_x_discrete(name="Adaptation conditions", labels=c("50%", "100%", "200%")) +
  scale_fill_manual(values= alpha(c("#3B9AB2", "#EBCC2A", "#F21A00", "#FF00FF")), 0.9) + # this is the zissou1 palette
  scale_color_manual(values=c("#3B9AB2", "#EBCC2A", "#F21A00", "#FF00FF")) # this is the zissou1 palette
if (saveFigures){
  ggsave(file=paste(currentTitle, "_width.eps", sep=""), device="eps", dpi = 1200, width = 12, height = 10, units = "cm")
  ggsave(file=paste(currentTitle, "_width.png", sep=""), dpi = 600, width = 12, height = 10, units = "cm")
  library(Cairo)
  ggsave(file=paste(currentTitle, "_width.pdf", sep=""), device=cairo_pdf, dpi = 1200, width = 12, height = 10, units = "cm")
  
}




ggplot(d, aes(x=factor(mediumConcentration, levels=sort(unique(as.numeric(mediumConcentration)))), y=medianElongation, col=factor(tAdaptation))) +
  geom_boxplot(width=0.6, outlier.shape = NA, colour="black", aes(fill = factor(tAdaptation))) + # add a box plot
  theme_classic(base_size = 22) +
  theme(legend.position = "none") + 
  scale_y_continuous(name=expression(paste("Aspect ratio")), lim=c(0,4)) +
  scale_x_discrete(name="Adaptation conditions", labels=c("50%", "100%", "200%")) +
  scale_fill_manual(values= alpha(c("#3B9AB2", "#EBCC2A", "#F21A00", "#FF00FF")), 0.9) + # this is the zissou1 palette
  scale_color_manual(values=c("#3B9AB2", "#EBCC2A", "#F21A00", "#FF00FF")) # this is the zissou1 palette
if (saveFigures){
  ggsave(file=paste(currentTitle, "_elongation.eps", sep=""), device="eps", dpi = 1200, width = 12, height = 10, units = "cm")
  ggsave(file=paste(currentTitle, "_elongation.png", sep=""), dpi = 600, width = 12, height = 10, units = "cm")
  library(Cairo)
  ggsave(file=paste(currentTitle, "_elongation.pdf", sep=""), device=cairo_pdf, dpi = 1200, width = 12, height = 10, units = "cm")
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



# here plot with error bars
# Summarize the data in order to plot the data with error bar
d2 <- data_summary(d1, varname="estimatedlogVolume", 
                   groupnames=c("tAdaptation", "mediumConcentration", "isMotherCulture"))

d2subculture <- subset(d2, isMotherCulture == FALSE)
d2mother <- subset(d2, isMotherCulture == TRUE)

if (skipMotherCulture)
{
  d2 <- d2subculture
}

plotG2bis <- ggplot(d2, aes(x=factor(mediumConcentration, levels=sort(unique(as.numeric(mediumConcentration)))), estimatedlogVolume, colour=factor(tAdaptation))) +
  # geom_violin(aes(fill = factor(TreatmentNumber))) +
  # coord_cartesian(ylim = c(0,150000)) + # limit y axis to 150000
  # stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="red") +
  # scale_fill_grey() + 
  # geom_boxplot(width=0.6, outlier.shape = NA, colour="black", aes(fill = factor(tAdapt))) + # add a box plot
  geom_errorbar(aes(ymin=estimatedlogVolume-sd, ymax=estimatedlogVolume+sd), width=0.1, position=position_dodge(width=0.2)) + 
  geom_point(size=6, shape=21, position=position_dodge(width=0.2), aes(col=factor(tAdaptation), fill=factor(tAdaptation))) + 
  theme_classic(base_size = 22) +
  theme(legend.position = "none") + 
  scale_y_continuous(name=expression(paste('log'[10]*'(volume)'," ",  mu, 'm'^3))) +
  scale_x_discrete(name="Adaptation conditions", labels=c("50%", "100%", "200%")) +
  scale_fill_manual(values= alpha(c("#3B9AB2", "#EBCC2A", "#F21A00", "#FF00FF")), 0.9) + # this is the zissou1 palette
  scale_color_manual(values=c("#000000", "#000000", "#000000", "#FF00FF")) # this is the zissou1 palette
plotG2bis

if (saveFigures){
  ggsave(file=paste(currentTitle, "_logvolume_errorbar.eps", sep=""), device="eps", dpi = 1200, width = 12, height = 10, units = "cm")
  library(Cairo)
  ggsave(file=paste(currentTitle, "_logvolume_errorbar.pdf", sep=""), device=cairo_pdf, dpi = 1200, width = 12, height = 10, units = "cm")
  ggsave(file=paste(currentTitle, "_logvolume_errorbar.png", sep=""), dpi = 600, width = 12, height = 10, units = "cm")
}


# Let's have a look at the changes of size relative to the mother culture
# assuming that cells have the same size as in the mother culture at time 0,
# when they first are subcultured into the new medium, and they have their own
# measured size 30 days later.
# There are some assumptions here, that these start and end dates are correct (while in
# practice there is some variation on the dates when the cells are measured
# and that the trend of log volume vs. time is linear (which probably is not
# the case as the volume changes a lot at the beginning. However, similar
# assumptions were also made when fitting linear curves to the population growth
# rates during adaptation and subcultures.
# assigne time zero to the mother culture and 30 to the subcultures
d1$time <- (1 - as.numeric(d1$isMotherCulture))*30
d1subculture <- subset(d1, isMotherCulture == FALSE)
d1mother <- subset(d1, isMotherCulture == TRUE)

#define standard error of mean function
std.error <- function(x) sd(x)/sqrt(length(x))

iteration<-0
for (aaa in 1:length(allTAdaptation))
{
  for (mmm in 1:length(allMediumConcentrations))
  { 
    d1current0 <- subset(d1subculture, tAdaptation == allTAdaptation[aaa] & mediumConcentration == allMediumConcentrations[mmm])
    d1current <- bind_rows(d1current0, d1mother)
    fittedSizeChange <- lm(formula = estimatedlogVolume ~ time, data=d1current)
    print(paste("T:", allTAdaptation[aaa], "; C=", allMediumConcentrations[mmm], "; Volume=10^(", round(fittedSizeChange$coefficients[1], 5) , "+", round(fittedSizeChange$coefficients[2],5), "t)", sep=""))
    deltaBodySize <- data.frame(tAdaptation=allTAdaptation[aaa], mediumConcentration=allMediumConcentrations[mmm], sizeT0=mean(d1mother$estimatedlogVolume), sizeT0sd=sd(d1mother$estimatedlogVolume), sizeT0se=std.error(d1mother$estimatedlogVolume), sizeT30=mean(d1current0$estimatedlogVolume), sizeT30sd=sd(d1current0$estimatedlogVolume), sizeT30se=std.error(d1current0$estimatedlogVolume), fitIntercept=fittedSizeChange$coefficients[1], fitSlope=fittedSizeChange$coefficients[2], row.names=NULL)
    
    if (iteration == 0)
    {
      bodySizeFitResults <- deltaBodySize
    } else
    {
      bodySizeFitResults <- rbind(bodySizeFitResults, deltaBodySize)
    }
    iteration <- iteration + 1
  }
}

write.table(bodySizeFitResults, "estimated_body_size_change_during_adaptation.csv", append = FALSE, sep = ",", row.names = FALSE)


################# DP CODE ######################
# need to specify that variables are categorical (rather than continuous)
d$tAdaptation <- as.factor(d$tAdaptation)
d$mediumConcentration <- as.factor(d$mediumConcentration)
str(d)

# start with ANOVA analysis with interaction between adaptation temp and concentration
m <- lm(estimatedlogVolume ~ tAdaptation * mediumConcentration, data = d)
plot(m) # model is nicely behaved and assumption are met
summary(m) #R-squared values indicates approx 34% of variation in data is explained

# because we have multiple measurements from each culture (line) we should use linear-mixed effects 
# models to account for this

# fit using lme4 package but this does not give p-values
library(lme4)
m <- lmer(estimatedlogVolume ~ tAdaptation * mediumConcentration + (1 | line), data = d)
summary(m)

# check model performace 
library(performance)
mod.check <- check_model(m)
mod.check # all looks good. High multi-collinearity but this is not a major issue as we are not interested in estimating parameters

library(lmerTest)
anova(m)

library('emmeans')
emmeans(m, ~ tAdaptation)
emmeans(m, ~ mediumConcentration)
emmeans(m, ~ tAdaptation * mediumConcentration)

# cell shrinkage at low nutrients
(4.45 - 4.12) / 4.45 *100

# cell shrinkage at intermediate nutrients
(4.46 - 4.23) / 4.46 *100

# cell shrinkage at high nutrients
(4.30 - 4.23) / 4.30 *100

