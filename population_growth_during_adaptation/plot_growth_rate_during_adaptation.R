rm(list=ls()) # clean memory
if(!is.null(dev.list())) dev.off()


#################################################################################################################
# Function to calculate the mean and the standard deviation
# for each group (to plot as errorbars, rather than as a boxplot)
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
  data_sum<-plyr::ddply(data, groupnames, .fun=summary_func,
                        varname)
  data_sum <- plyr::rename(data_sum, c("mean" = varname))
  return(data_sum)
}



# file name file.choose() # 
fileName <- "~/Tetrahymena/population_growth_during_adaptation/growth_data_during_adaptation.csv"
# The following file is produced by the script analysis_of_body_size_from_tracking_data.R
fileNameBodySizeData <- "/Tetrahymena/body_size/estimated_body_size_change_during_adaptation.csv"

# change directory
setwd(dirname(fileName))

# read data from file
populationGrowthData <- read.table(file = fileName, sep = ",", header=TRUE, na.strings = c("NA", " NA"))

# header names
names(populationGrowthData)

library("stringr")

# create new columns in which splits information given by file name: on the population, treatment and date of experiment
newColumns <- str_split_fixed(as.character(str_trim(tolower(populationGrowthData$name))), "\\.|-|_", 3)
colnames(newColumns) <- c("tAdapt", "density", "line")


# convert strings to numbers
library(readr)
newColumns[,1:2] = (parse_number(newColumns[,1:2] ))

# add the new columns to the original table
populationGrowthData <- cbind(populationGrowthData, newColumns, deparse.level = 1, stringsAsFactors = FALSE)

# Find the data for the mother culture (as it has same density and temperature, but I prefer to consider it as a separate culture)
populationGrowthData$isMotherCulture <- as.numeric(populationGrowthData$line == 1 | populationGrowthData$line == 2)

# also create a column with both density and temperature
populationGrowthData$tAdapt_and_density <- paste("T", populationGrowthData$tAdapt, "_C", populationGrowthData$density, "_M", populationGrowthData$isMotherCulture, sep="")

plotColours <- rep(c("#3B9AB2", "#EBCC2A", "#F21A00"),5)
# plotColoursFill <- c("#18434E", "#45ABC4", "#8FCDDC", "#463C07", "#E4C315", "#F2DD70", "#660B00", "#FF250A", "#FF7C6C", "#3B9AB2", "#EBCC2A", "#F21A00", "black")
markerShapes = c(21, 24, 22, 4)
markerFillColours = rep(c("#666666", "#A3A3A3", "#000000", "#FF00FF"),4)


# Convert date to numeric
# These columns are here to extract the date as a single number, so first makes a new Date column and then converts it with as.numeric
populationGrowthData$date <- as.POSIXlt(populationGrowthData$current_date, format="%Y-%m-%d")
populationGrowthData$date_number <- as.numeric(populationGrowthData$date)/86400 # convert to units of days
populationGrowthData$date_number <- populationGrowthData$date_number - min(populationGrowthData$date_number)

# isolate the data for the mother culture (before splitting)
motherCulture <- subset(populationGrowthData, isMotherCulture == 1)

library(ggplot2)
ggplot(motherCulture, aes(x=date_number, y=growth_rate_gen_per_day, shape=line)) +
  geom_point(size=4, alpha=0.9, aes(colour=factor(line))) + 
  geom_smooth(method=lm, formula='y~x', se=TRUE, fullrange=TRUE, colour="black", aes(fill=factor(tAdapt))) + 
  theme_classic(base_size = 22) +
  scale_y_continuous(name="generations per day") +
  scale_x_continuous(name="adaptation time (days)") + 
  theme(legend.position = "none") + 
  scale_color_manual(values=c("#EBCC2A", "#3B9AB2", "#F21A00")) + # this is the zissou1 palette
  scale_fill_manual(values= c("#EBCC2A", "#3B9AB2", "#F21A00"))# this is the zissou1 palette

ggsave(file="tetrahymena_mother_culture_growth_rate.png", dpi = 600, width = 15, height = 12, units = "cm")
ggsave(file="tetrahymena_mother_culture_growth_rate.eps", device="eps", dpi = 1200, width = 15, height = 10, units = "cm")

# The cell count is not very reliable if the density is below a certain value
# Say that we need at least 3 cells in 6 squares of the grid
# As each square is 0.1 uL, this is a density of 3/0.6*1000=5000 cells/ml
# populationGrowthData <- subset(populationGrowthData, current_count >= 5000 & previous_count >= 5000)

plotG1 <- ggplot(populationGrowthData, aes(x=date_number, y=growth_rate_gen_per_day, shape=factor(as.numeric(density)), color=tAdapt)) +
  geom_point(size=2) + 
  geom_smooth(method=lm, formula='y~x', se=FALSE, fullrange=FALSE, aes(fill=factor(tAdapt_and_density), linetype=factor(density))) + 
  theme_classic(base_size = 22) +
  scale_y_continuous(name="generations per day", limits=c(0,5.7), breaks=seq(0,5, by=1)) +
  scale_x_continuous(name="adaptation time (days)") + 
  theme(legend.position = "none") + 
  scale_color_manual(values=plotColours) + # this is the zissou1 palette
  scale_fill_manual(values= alpha(plotColours), 0.9)# this is the zissou1 palette
plotG1


plotG1 <- ggplot(populationGrowthData, aes(x=date_number, y=growth_rate_gen_per_day, shape=factor(as.numeric(density)), color=tAdapt, fill=tAdapt)) +
  # geom_point(size=2, alpha=0.4, colour="black") +
  # geom_jitter(size=2, alpha=0.4, colour="black", position = position_jitter(height = 0, width = .4)) +
  geom_smooth(method=lm, formula='y~x', se=FALSE, fullrange=FALSE, aes(fill=factor(tAdapt_and_density), size=factor(as.numeric(density)))) + 
  theme_classic(base_size = 22) +
  scale_y_continuous(name="generations per day", limits=c(0,5.7), breaks=seq(0,5, by=1)) +
  scale_x_continuous(name="adaptation time (days)") + 
  theme(legend.position = "none") + 
  scale_color_manual(values=plotColours) + # this is the zissou1 palette
  scale_shape_manual(values=markerShapes) + # shapes for the markers
  scale_size_manual( values = c(0.8, 1.4, 2) ) +
  # scale_alpha_continuous(range=c(0.5, 1)) + 
  scale_fill_manual(values= plotColours)# this is the zissou1 palette
plotG1

ggsave(file="tetrahymena_growth_rate_during_adaptation_no_markers.png", dpi = 600, width = 15, height = 10, units = "cm")
ggsave(file="tetrahymena_growth_rate_during_adaptation_no_markers.eps", device="eps", dpi = 1200, width = 15, height = 10, units = "cm")
library(Cairo)
ggsave(file="tetrahymena_growth_rate_during_adaptation_no_markers.pdf", device="pdf", dpi = 1200, width = 15, height = 10, units = "cm")

# find the end of the subculture
endOfSubcultures <- max(populationGrowthData$date_number)
nLastDays <- 25
# and isolate the last 15 days (to exclude the transient when everything happens)
populationGrowthDataStableSubcultures <- subset(populationGrowthData, date_number > endOfSubcultures - nLastDays)
populationGrowthDataStableSubcultures$date_number_from_subculture <- populationGrowthDataStableSubcultures$date_number -54

# adding new columns with values as factors to simplify facet titles,
# sorting of plots in a figure, etc.
populationGrowthDataStableSubcultures$density_as_factor = factor(paste(populationGrowthDataStableSubcultures$density, "%", sep=""), levels=paste(as.character(sort(as.numeric(unique(populationGrowthDataStableSubcultures$density)))), "%", sep=""))
populationGrowthDataStableSubcultures$tAdapt_as_factor = factor(paste(populationGrowthDataStableSubcultures$tAdapt, "°C", sep=""), levels=paste(as.character(sort(unique(populationGrowthDataStableSubcultures$tAdapt))), "°C", sep=""))

growthRateStableSubculturesFit <- as.data.frame(populationGrowthDataStableSubcultures 
                               %>% group_by(tAdapt_and_density, tAdapt, density) 
                               %>% do(
                                 popfitIntercept = (coefficients(summary(lm(growth_rate_gen_per_day ~ date_number_from_subculture, data = .)))[1,1]), 
                                 popfitSlope = (coefficients(summary(lm(growth_rate_gen_per_day ~ date_number_from_subculture, data = .)))[2,1]),
                                 popfitInterceptSe = (coefficients(summary(lm(growth_rate_gen_per_day ~ date_number_from_subculture, data = .)))[1,2]), 
                                 popfitSlopeSe = (coefficients(summary(lm(growth_rate_gen_per_day ~ date_number_from_subculture, data = .)))[2,2]))
)
sapply(growthRateStableSubculturesFit, typeof) # sapply(growthRateStableSubculturesFit, class)
growthRateStableSubculturesFit$popfitIntercept <- as.double(growthRateStableSubculturesFit$popfitIntercept)
growthRateStableSubculturesFit$popfitSlope <- as.double(growthRateStableSubculturesFit$popfitSlope)
growthRateStableSubculturesFit$popfitInterceptSe <- as.double(growthRateStableSubculturesFit$popfitInterceptSe)
growthRateStableSubculturesFit$popfitSlopeSe <- as.double(growthRateStableSubculturesFit$popfitSlopeSe)
sapply(growthRateStableSubculturesFit, typeof) # sapply(growthRateStableSubculturesFit, class)

# This shows which fits are decreasing or increasing, or stable over time:
popGrowthRateFitStableSubculture <- cbind(growthRateStableSubculturesFit$tAdapt_and_density, growthRateStableSubculturesFit$popfitSlope, growthRateStableSubculturesFit$popfitSlope - growthRateStableSubculturesFit$popfitSlopeSe, growthRateStableSubculturesFit$popfitSlope + growthRateStableSubculturesFit$popfitSlopeSe)
popGrowthRateFitStableSubculture
write.table(popGrowthRateFitStableSubculture, "changes_of_pop_growth_rate_late_in_adaptation.csv", append = FALSE, sep = ",", row.names = FALSE)


plotG1err <- ggplot(populationGrowthDataStableSubcultures, aes(x=date_number - endOfSubcultures+nLastDays, y=growth_rate_gen_per_day, shape=factor(as.numeric(density)), color=tAdapt, fill=tAdapt)) +
  geom_point(size=2, alpha=0.4, colour="black") +
  # geom_jitter(size=2, alpha=0.4, colour="black", position = position_jitter(height = 0, width = .4)) +
  geom_smooth(method=lm, formula='y~x', se=TRUE, fullrange=FALSE, aes(fill=factor(tAdapt), size=factor(as.numeric(density)))) + 
  theme_classic(base_size = 14) +
  scale_y_continuous(name="generations per day", limits=c(0,5), breaks=seq(0,5, by=2)) +
  scale_x_continuous(name="adaptation time (days)", limits=c(0,nLastDays), breaks=seq(0,nLastDays, by=10)) + 
  theme(legend.position = "none") + 
  scale_color_manual(values=plotColours) + # this is the zissou1 palette
  scale_shape_manual(values=markerShapes) + # shapes for the markers
  scale_size_manual( values = c(0.8, 1.4, 2) ) +
  # scale_alpha_continuous(range=c(0.5, 1)) + 
  scale_fill_manual(values= plotColours) # this is the zissou1 palette


plotG1err <- plotG1err + facet_grid(rows=vars(tAdapt_as_factor), cols=vars(density_as_factor))
plotG1err

ggsave(file="tetrahymena_growth_rate_during_adaptation_no_markers_err.png", dpi = 600, width = 15, height = 12, units = "cm")
ggsave(file="tetrahymena_growth_rate_during_adaptation_no_markers_err.eps", device="eps", dpi = 1200, width = 15, height = 12, units = "cm")
library(Cairo)
ggsave(file="tetrahymena_growth_rate_during_adaptation_no_markers_err.pdf", device="pdf", dpi = 1200, width = 15, height = 12, units = "cm")



plotG1 <- ggplot(populationGrowthData, aes(x=date_number, y=growth_rate_gen_per_day, shape=factor(as.numeric(density)), color=tAdapt, fill=tAdapt)) +
  geom_point(size=2, alpha=0.4, colour="black") +
  # geom_jitter(size=2, alpha=0.4, colour="black", position = position_jitter(height = 0, width = .4)) +
  geom_smooth(method=lm, formula='y~x', se=FALSE, fullrange=FALSE, aes(fill=factor(tAdapt_and_density), size=factor(as.numeric(density)))) + 
  theme_classic(base_size = 22) +
  scale_y_continuous(name="generations per day", limits=c(0,5.7), breaks=seq(0,5, by=1)) +
  scale_x_continuous(name="adaptation time (days)") + 
  theme(legend.position = "none") + 
  scale_color_manual(values=plotColours) + # this is the zissou1 palette
  scale_shape_manual(values=markerShapes) + # shapes for the markers
  scale_size_manual( values = c(0.8, 1.4, 2) ) +
  # scale_alpha_continuous(range=c(0.5, 1)) + 
  scale_fill_manual(values= plotColours)# this is the zissou1 palette
plotG1

ggsave(file="tetrahymena_growth_rate_during_adaptation.png", dpi = 600, width = 15, height = 10, units = "cm")
ggsave(file="tetrahymena_growth_rate_during_adaptation.eps", device="eps", dpi = 1200, width = 15, height = 10, units = "cm")
library(Cairo)
ggsave(file="tetrahymena_growth_rate_during_adaptation.pdf", device="pdf", dpi = 1200, width = 15, height = 10, units = "cm")









# Isolate the data for the 9 experimental cultures after splitting the mother culture

experimentalCultures <- subset(populationGrowthData, isMotherCulture == 0)






plotG2 <- ggplot(experimentalCultures, aes(x=factor(density, levels=unique(as.numeric(density))), growth_rate_gen_per_day, col=tAdapt_and_density)) +
  # geom_violin(aes(fill = factor(TreatmentNumber))) +
  #coord_cartesian(ylim = c(0,150000)) + # limit y axis to 150000
  # stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="red") +
  # scale_fill_grey() + 
  geom_boxplot(width=0.6, outlier.shape = NA, colour="black", aes(fill = factor(tAdapt))) + # add a box plot
  theme_classic(base_size = 22) +
  theme(legend.position = "none") + 
  scale_y_continuous(name="generations per day", limits=c(0,5.7), breaks=seq(0,5, by=1)) +
  scale_x_discrete(name="Adaptation conditions", labels=c("50%", "100%", "200%")) +
  scale_fill_manual(values= alpha(c("#3B9AB2", "#EBCC2A", "#F21A00", "#FF00FF")), 0.9) + # this is the zissou1 palette
  scale_color_manual(values=c("#3B9AB2", "#EBCC2A", "#F21A00", "#FF00FF")) # this is the zissou1 palette
plotG2

ggsave(file="avg_growth_rate_during_adaptation.eps", device="eps", dpi = 1200, width = 12, height = 10, units = "cm")
ggsave(file="avg_growth_rate_during_adaptation.png", dpi = 600, width = 12, height = 10, units = "cm")
library(Cairo)
ggsave(file="avg_growth_rate_during_adaptation.pdf", device=cairo_pdf, dpi = 1200, width = 12, height = 10, units = "cm")



# here plot with error bars
# Summarize the data in order to plot the data with error bar
d2 <- data_summary(experimentalCultures, varname="growth_rate_gen_per_day", 
                   groupnames=c("tAdapt", "density", "tAdapt_and_density"))

plotG2bis <- ggplot(d2, aes(x=factor(density, levels=sort(unique(as.numeric(density)))), growth_rate_gen_per_day, colour=factor(tAdapt))) +
  # geom_violin(aes(fill = factor(TreatmentNumber))) +
  # coord_cartesian(ylim = c(0,150000)) + # limit y axis to 150000
  # stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="red") +
  # scale_fill_grey() + 
  # geom_boxplot(width=0.6, outlier.shape = NA, colour="black", aes(fill = factor(tAdapt))) + # add a box plot
  geom_errorbar(aes(ymin=growth_rate_gen_per_day-sd, ymax=growth_rate_gen_per_day+sd), width=0.1, position=position_dodge(width=0.2)) + 
  geom_point(size=6, shape=21, position=position_dodge(width=0.2), aes(col=factor(tAdapt), fill=factor(tAdapt))) + 
  theme_classic(base_size = 22) +
  theme(legend.position = "none") + 
  scale_y_continuous(name="generations per day", limits=c(0,5.7), breaks=seq(0,5, by=1)) +
  scale_x_discrete(name="Adaptation conditions", labels=c("50%", "100%", "200%")) +
  scale_fill_manual(values= alpha(c("#3B9AB2", "#EBCC2A", "#F21A00", "#FF00FF")), 0.9) + # this is the zissou1 palette
  scale_color_manual(values=c("#000000", "#000000", "#000000", "#FF00FF")) # this is the zissou1 palette
plotG2bis

ggsave(file="avg_growth_rate_during_adaptation_with_sd_errorbar.eps", device="eps", dpi = 1200, width = 12, height = 10, units = "cm")
ggsave(file="avg_growth_rate_during_adaptation_with_sd_errorbar.png", dpi = 600, width = 12, height = 10, units = "cm")
library(Cairo)
ggsave(file="avg_growth_rate_during_adaptation_with_sd_errorbar.pdf", device=cairo_pdf, dpi = 1200, width = 12, height = 10, units = "cm")



library(ggpubr)
# combine together the two plots into a single figure
figure <- ggarrange(plotG1, plotG2,
                    # labels = c("A", "B"),
                    ncol = 2, nrow = 1)
figure

ggsave(file="figure_growth_rate_during_adaptation.eps", device="eps", dpi = 1200, width = 22, height = 10, units = "cm")
ggsave(file="figure_growth_rate_during_adaptation.png", dpi = 600, width = 22, height = 10, units = "cm")
library(Cairo)
ggsave(file="figure_growth_rate_during_adaptation.pdf", device=cairo_pdf, dpi = 1200, width = 22, height = 10, units = "cm")



# Look at the density of cells in each culture over time
# we can first check that the data below are the same:
# cbind(populationGrowthData$starting_count*2^(populationGrowthData$culture_time_days*populationGrowthData$growth_rate_gen_per_day), populationGrowthData$current_count)
# cbind(log2(populationGrowthData$current_count/populationGrowthData$starting_count)/populationGrowthData$culture_time_days, populationGrowthData$growth_rate_gen_per_day)

# Here are the densities hour by hour:
# densities_by_hour <- populationGrowthData$starting_count[jjj]*2^(populationGrowthData$culture_time_days[jjj]*seq(0,populationGrowthData$growth_rate_gen_per_day[jjj],by=1/24))
bigList <- list()
for (jjj in 1:dim(populationGrowthData)[1])
{
  print(jjj)
  if (is.finite(populationGrowthData$growth_rate_gen_per_day[jjj]))
  {
    culture_time_by_hour <- seq(0,populationGrowthData$growth_rate_gen_per_day[jjj],by=1/24)
    # important: note that the calculation below assumes exponential growth at all densities
    densities_by_hour <- populationGrowthData$starting_count[jjj]*2^(populationGrowthData$culture_time_days[jjj]*culture_time_by_hour)
    nElements <- length(culture_time_by_hour)
    repLine <- rep(populationGrowthData$line[jjj], nElements)
    repTAdapt <- rep(populationGrowthData$tAdapt[jjj], nElements)
    repDensity <- rep(populationGrowthData$density[jjj], nElements)
    repIsMotherCulture <- rep(populationGrowthData$isMotherCulture[jjj], nElements)
    repTAdapt_and_density <- rep(populationGrowthData$tAdapt_and_density[jjj], nElements)
    
    
    bigList[[jjj]] <- cbind(repLine, repTAdapt, repTAdapt_and_density, repDensity, culture_time_by_hour, densities_by_hour, repIsMotherCulture)
  }
}

bigMatrix <- do.call("rbind", bigList)
colnames(bigMatrix) <- c("line", "tAdapt", "tAdapt_and_density", "density", "culture_time_by_hour", "densities_by_hour", "isMotherCulture")

populationGrowthDataByHour <- as.data.frame(bigMatrix)

#changing columns to factors/numeric and date variables
populationGrowthDataByHour <- populationGrowthDataByHour %>% dplyr::mutate_at(c('tAdapt', 'density', 'culture_time_by_hour', 'densities_by_hour'), as.numeric)

# adding new columns with values as factors to simplify facet titles,
# sorting of plots in a figure, etc.
populationGrowthDataByHour$density_as_factor = factor(paste(populationGrowthDataByHour$density, "%", sep=""), levels=paste(as.character(sort(as.numeric(unique(populationGrowthDataByHour$density)))), "%", sep=""))
populationGrowthDataByHour$tAdapt_as_factor = factor(paste(populationGrowthDataByHour$tAdapt, "°C", sep=""), levels=paste(as.character(sort(unique(populationGrowthDataByHour$tAdapt))), "°C", sep=""))



experimentalCulturesByHour <- subset(populationGrowthDataByHour, isMotherCulture == 0)

# check below whether the plot is done on populationGrowthDataByHour or on experimentalCulturesByHour
plotDensity <- ggplot(data=experimentalCulturesByHour, aes(x=log10(densities_by_hour), fill=as.factor(density_as_factor))) +
  geom_histogram(aes(y= ..density..), breaks=seq(1, 7, by = 1), closed="left") +
  geom_density(lwd = 1.5, alpha = 0.25, aes(colour=as.factor(tAdapt), fill=as.factor(density_as_factor))) +
  # geom_density(lwd = 1, alpha = 0, aes(colour=as.factor(tAdapt), fill=as.factor(density_as_factor))) +
  geom_vline(xintercept=rep(6,dim(populationGrowthDataByHour)[1]), color="green", linetype="dashed") +
  theme_classic() +
  theme(legend.position="none") +
  scale_color_manual(values=c("#3B9AB2", "#EBCC2A", "#F21A00", "#FF00FF")) +
  scale_fill_manual(values=c("#A3A3A3", "#666666", "#000000")) +
  scale_y_continuous(name="fraction of time") +
  scale_x_continuous(name="density (cells/ml)", breaks=seq(2,6, by=2), labels=c('100', '10k', '1M'))
plotDensity <- plotDensity + facet_grid(rows=vars(tAdapt_as_factor), cols=vars(density_as_factor))

plotDensity

meanDensity <- mean(experimentalCulturesByHour$densities_by_hour)
medianDensity <- median(experimentalCulturesByHour$densities_by_hour)
meanLogDensity <- mean(log10(experimentalCulturesByHour$densities_by_hour))
print(paste("The mean density across all conditions is ", meanDensity, "; mean log(density)=", meanLogDensity, sep=""))
print(paste("The median density across all conditions is ", medianDensity, sep=""))

# look at the proportion of extreme values in the reconstructed data (with the caveat
# that we assume exponential growth
sum(experimentalCulturesByHour$densities_by_hour > 500000)
length(experimentalCulturesByHour$densities_by_hour)
# or in the orginal data (now the caveat is the opposite, that we look at extreme values
# while these were probably reached less than one generation before, because of the 
# exponential growth.)
sum(experimentalCultures$current_count > 500000)
length(experimentalCultures$current_count)

ggsave(file="density_over_subculture_time.eps", device="eps", dpi = 1200, width = 12, height = 10, units = "cm")
ggsave(file="density_over_subculture_time.png", dpi = 600, width = 12, height = 10, units = "cm")
library(Cairo)
ggsave(file="density_over_subculture_time.pdf", device=cairo_pdf, dpi = 1200, width = 12, height = 10, units = "cm")



# adding new columns with values as factors to simplify facet titles,
# sorting of plots in a figure, etc.
experimentalCultures$density_as_factor = factor(paste(experimentalCultures$density, "%", sep=""), levels=paste(as.character(sort(as.numeric(unique(experimentalCultures$density)))), "%", sep=""))
experimentalCultures$tAdapt_as_factor = factor(paste(experimentalCultures$tAdapt, "°C", sep=""), levels=paste(as.character(sort(unique(experimentalCultures$tAdapt))), "°C", sep=""))


# fit logistic growth curve (all data combined):
# tempData <- experimentalCultures[1,]
# tempData$growth_rate_gen_per_day[1] <- 0 # add an extra data point with no growth
# tempData$current_count[1] <- 2000000
# tempData <- rbind(experimentalCultures, tempData)
# fit <- nls(data=tempData[c('growth_rate_gen_per_day', 'current_count')], formula=growth_rate_gen_per_day ~ a * (1 - current_count/b), start=list(a=2, b=1000000))

allFittedCurves <- data.frame(current_count=double(), .fitted=double(), tAdapt=character(), density=character())
# fit logistic growth curve (one for each subplot)
for (ddd in unique(experimentalCultures$density))
{
  for (ttt in unique(experimentalCultures$tAdapt))
  {
    experimentalCultures1 <- subset(experimentalCultures, tAdapt == ttt & density == ddd)
    # for 15 degrees and 200% medium, the fit fails because of a high value of 
    # generations per day when the current_count is very low. However, if the final count is low
    # the estimation of population density is less accurate (fewer cells under the microscope
    # a density of 10000 cells per ml means 1 cell per 0.1 uL square)
    # and so I simply remove it to force the fit.
    if (ttt == "15" & ddd == "200")
    {
      # which_min_current_count <- which.min(experimentalCultures1$current_count)
      # experimentalCultures1 <- experimentalCultures1[-which_min_current_count,]
      experimentalCultures1 <- subset(experimentalCultures1, current_count >= 5000)
    }
    # tempData <- experimentalCultures1[1,]
    # tempData$growth_rate_gen_per_day[1] <- 0 # add an extra data point with no growth
    # tempData$current_count[1] <- 2000000
    # tempData <- rbind(experimentalCultures1, tempData)
    fit1 <- nls(data=experimentalCultures1[c('growth_rate_gen_per_day', 'current_count')], formula=growth_rate_gen_per_day ~ a * (1 - current_count/b), start=list(a=2, b=1000000))
    fit2 <- nls(data=experimentalCultures1[c('growth_rate_gen_per_day', 'current_count')], formula=growth_rate_gen_per_day ~ a , start=list(a=2))
    
    aic_on_models <- AIC(fit1,fit2)
    bestModel <- which.min(aic_on_models$AIC)
    if (bestModel == 1)
    {
      fit <- fit1
      print(paste("density=", ddd, "; temperature=", ttt, "; best fit is saturating growth - dAIC=", round(diff(aic_on_models$AIC), 4), sep=""))
      print(fit)
    } else 
    {
      fit <- fit2
      print(paste("density=", ddd, "; temperature=", ttt, "; best fit is growth independent of density - dAIC=", round(diff(aic_on_models$AIC), 4), sep=""))
    }
    
    # compare the two fits:
    # library(AICcmodavg)
    # models <- list(fit, fit2) #define list of models
    # mod.names <- c('logistic.growth', 'flat') #specify model names
    # aictab(cand.set = models, modnames = mod.names) #calculate AIC of each model
    
    # draw fitted curve
    new_data <- data.frame(current_count = 10^(seq(2, 6, 0.1)))
    library(broom)
    preds <- broom::augment(fit, newdata = new_data)
    preds$tAdapt <- ttt
    preds$density <- ddd
    preds$tAdapt_as_factor <- experimentalCultures1$tAdapt_as_factor[1]
    preds$density_as_factor <- experimentalCultures1$density_as_factor[1]
    allFittedCurves <- rbind(allFittedCurves, preds)
  }
}

# Look at the effect of final density on the growth rate, to check for 
# deviations from the exponential growth
plotDensityEffects <- ggplot(experimentalCultures, aes(x=log10(current_count), y=growth_rate_gen_per_day, shape=density_as_factor, color=tAdapt_as_factor, fill=density_as_factor)) +
  geom_point(size=2, alpha=1, colour="black") +
  # geom_line(aes(current_count, .fitted), data=preds, col="blue") +
  geom_line(aes(log10(current_count), .fitted, size=density_as_factor), data=allFittedCurves) +
  # geom_smooth(method='lm', formula = y ~ x, se=FALSE, fullrange=FALSE, aes(fill=factor(tAdapt_and_density), size=factor(as.numeric(density)))) + 
  # geom_jitter(size=2, alpha=0.4, colour="black", position = position_jitter(height = 0, width = .4)) +
  theme_classic() +
  scale_y_continuous(name="generations per day", limits=c(0,5.7), breaks=seq(0,5, by=1)) +
  scale_x_continuous(name="maximum subculture density (cells/ml)", breaks=seq(2,6, by=2), labels=c('100', '10k', '1M'), limits=c(1,7)) + 
  theme(legend.position = "none") + 
  scale_color_manual(values=plotColours) + # this is the zissou1 palette
  scale_shape_manual(values=markerShapes) + # shapes for the markers
  scale_size_manual( values = c(0.8, 1.4, 2) ) +
  # scale_alpha_continuous(range=c(0.5, 1)) + 
  scale_fill_manual(values=c("#A3A3A3", "#666666", "#000000"))# this is the zissou1 palette
plotDensityEffects <- plotDensityEffects + facet_grid(rows=vars(tAdapt_as_factor), cols=vars(density_as_factor))

plotDensityEffects

ggsave(file="growth_vs_maximum_density.eps", device="eps", dpi = 1200, width = 12, height = 10, units = "cm")
ggsave(file="growth_vs_maximum_density.png", dpi = 600, width = 12, height = 10, units = "cm")
library(Cairo)
ggsave(file="growth_vs_maximum_density.pdf", device=cairo_pdf, dpi = 1200, width = 12, height = 10, units = "cm")


  
  ## Let's have a look at the cumulative number of generations.
  populationGrowthData$n_generations_in_interval <- populationGrowthData$growth_rate_gen_per_day * populationGrowthData$culture_time_days
  populationGrowthData$n_generations_in_interval[is.na(populationGrowthData$n_generations_in_interval)] <- 0
  populationGrowthData$culture_time_days[is.na(populationGrowthData$culture_time_days)] <- 0
  populationGrowthData$cumulative_generations <- ave(populationGrowthData$n_generations_in_interval, populationGrowthData$name, FUN=cumsum)
  populationGrowthData$cumulative_time <- ave(populationGrowthData$culture_time_days, populationGrowthData$name, FUN=cumsum) + ave(populationGrowthData$date_number, populationGrowthData$name, FUN=min) - min(populationGrowthData$date_number[populationGrowthData$tAdapt == 25]) # I align the zero with the start of the cultures at 25 degrees (the cultures at 15 degrees start a bit earlier)
  # Note that in the data there are some missing counts at the very end, which means that the cumulative time
  # is not the same as the adaptation time. In other words, in the plot below I shouldn't
  # use date_number but the cumulative_time. 
  
  experimentalCultures <- subset(populationGrowthData, isMotherCulture == 0)
  
  
  plotCumG <- ggplot(experimentalCultures, aes(x=cumulative_time, y=cumulative_generations, shape=factor(as.numeric(density)), color=tAdapt, fill=tAdapt)) +
    # geom_point(size=2, alpha=0.4, colour="black") +
    geom_smooth(method=lm, formula='y~x', se=FALSE, fullrange=TRUE, aes(fill=factor(tAdapt_and_density), size=factor(as.numeric(density)))) + 
    geom_jitter(size=2, alpha=0.4, colour="black", position = position_jitter(height = 0, width = .4)) +
    theme_classic(base_size = 22) +
    scale_y_continuous(name=expression(paste("Num. generations")), limits=c(0,75)) +
    scale_x_continuous(name="adaptation time (days)", limits=c(-5, 21), breaks=seq(0, 21, by=7)) +
    coord_cartesian(xlim=c(-5,25)) + 
    theme(legend.position = "none") + 
    scale_color_manual(values=plotColours) + # this is the zissou1 palette
    scale_shape_manual(values=markerShapes) + # shapes for the markers
    scale_size_manual( values = c(0.8, 1.4, 2) ) +
    # scale_alpha_continuous(range=c(0.5, 1)) + 
    scale_fill_manual(values= plotColours)# this is the zissou1 palette
  plotCumG
  # note that the xlimits remove some data points, but mainly at the end when
  # we were no longer monitoring densities, and the experiment had started
  
  ggsave(file="tetrahymena_cumulative_generations.eps", device="eps", dpi = 1200, width = 12, height = 10, units = "cm")
  ggsave(file="tetrahymena_cumulative_generations.png", dpi = 600, width = 12, height = 10, units = "cm")
  library(Cairo)
  ggsave(file="tetrahymena_cumulative_generations.pdf", device=cairo_pdf, dpi = 1200, width = 12, height = 10, units = "cm")
  

  
  
  
  
  # run Anova
  # check for similar variances
  bartlett.test(growth_rate_gen_per_day ~ tAdapt_and_density, data = experimentalCultures)

   # Shapiro test to check if data in each group are normally distributed
  for (condition in unique(experimentalCultures$tAdapt_and_density))
  {
    print(condition)
  print(with(experimentalCultures, shapiro.test(growth_rate_gen_per_day[tAdapt_and_density == condition])))
  }

 # If the Shapiro-Wilk test indicates that the data are normally distributed, then we can go ahead and calculate the anova:
anova_result = aov(growth_rate_gen_per_day ~ tAdapt + density, data = experimentalCultures)
anova_result
summary(anova_result) # this calculates the p value of the test

#  If not, we can calculate a Kruskal-Wallis test:
kruskal.test(growth_rate_gen_per_day ~ density, data = experimentalCultures)



################# DP CODE ######################
# start with ANOVA analysis with interaction between adaptation temp and concentration
m <- lm(growth_rate_gen_per_day ~ tAdapt * density, data = experimentalCultures)
plot(m) # model is nicely behaved and assumption are met
summary(m) #R-squared values indicates approx 58% of variation in data is explained

# because we have multiple measurements from each culture (line) we should use linear-mixed effects 
# models to account for this

# fit using lme4 package but this does not give p-values
library(lme4)
m <- lmer(growth_rate_gen_per_day ~ tAdapt * density + (1 | line), data = experimentalCultures)
summary(m)

# check model performace 
library(performance)
mod.check <- check_model(m)
mod.check # all looks good. High multi-collinearity but this is not a major issue as we are not interested in estimating parameters
#  Could not compute standard errors from random effects, but this is fine.

library(lmerTest)
anova(m)

library('emmeans')
emmeans(m, ~ tAdapt)
emmeans(m, ~ density)

print("Now I try to estimate biomass production.")
print("In order to do this, I need to read an estimate of the change in cell volume")
print(paste("in the file", fileNameBodySizeData))
print("which is produced by the script analysis_of_body_size_from_tracking_data.R")
bodySizeFitResults <- read.table(file = fileNameBodySizeData, sep = ",", header=TRUE, na.strings = c("NA", " NA"))

bodySizeFitResults$tAdapt_and_density <- paste("T", bodySizeFitResults$tAdaptation, "_C", bodySizeFitResults$mediumConcentration, "_M0", sep="")

# fit straight lines to population growth data
populationGrowthData$date_number_from_subculture <- populationGrowthData$date_number -54
growthRateFit <- as.data.frame(populationGrowthData 
              %>% group_by(tAdapt_and_density, tAdapt, density) 
              %>% do(
              popfitIntercept = (coefficients(summary(lm(growth_rate_gen_per_day ~ date_number_from_subculture, data = .)))[1,1]), 
              popfitSlope = (coefficients(summary(lm(growth_rate_gen_per_day ~ date_number_from_subculture, data = .)))[2,1]),
              popfitInterceptSe = (coefficients(summary(lm(growth_rate_gen_per_day ~ date_number_from_subculture, data = .)))[1,2]), 
              popfitSlopeSe = (coefficients(summary(lm(growth_rate_gen_per_day ~ date_number_from_subculture, data = .)))[2,2]))
              )
sapply(growthRateFit, typeof) # sapply(growthRateFit, class)
growthRateFit$popfitIntercept <- as.double(growthRateFit$popfitIntercept)
growthRateFit$popfitSlope <- as.double(growthRateFit$popfitSlope)
growthRateFit$popfitInterceptSe <- as.double(growthRateFit$popfitInterceptSe)
growthRateFit$popfitSlopeSe <- as.double(growthRateFit$popfitSlopeSe)
sapply(growthRateFit, typeof) # sapply(growthRateFit, class)


dfChangesDuringAdaptation <- merge(growthRateFit, bodySizeFitResults, by="tAdapt_and_density")

# If the cells were keeping the same volume, the biomass production per unit biomass
# would be the same as the population growth rate.
# I have two equations for the change over time of Tetrahymena:
# I fitted a linear trendline to the log of cell volume
# log(V(t)) = log(V(0)) + a*t
# and the population growth is described by:
# dN/dt = r*N --> dN/N = rdt --> log(N(t)) = r*t + log(N0)
# these give:
# a = log(V(t)/V(0))/t
# and
# r = log(N(t)/N(0))/t
# I used log2 for the population growth rate and log10 for the body size
# If I wanted to use the natural log I should take r/log(2) and a/log(10)
# If I consider that r changes over time:
# dN/dt = r*(1+bt)*N --> dN/N = r(1+bt)dt --> log(N) = r*t + 1/2*b*t^2 + log(N0)
# r + 1/2bt = log(N(t)/N(0))/t
# I could calculate the biomass production per cell, or the biomass production per unit biomass
# the biomass production per unit biomass is the same as the population growth rate, except that at the end of one day the size 
# of individual cells has changed:
dfChangesDuringAdaptation$biomassProductionperDayT0 <- 10^dfChangesDuringAdaptation$fitSlope * 2^dfChangesDuringAdaptation$popfitIntercept
dfChangesDuringAdaptation$biomassProductionperDayT30 <- 10^dfChangesDuringAdaptation$fitSlope * 2^(dfChangesDuringAdaptation$popfitIntercept + dfChangesDuringAdaptation$popfitSlope * 30)

dfChangesDuringAdaptation$biomassProductionum3perCellPerDayT0 <- 10^dfChangesDuringAdaptation$sizeT0 * 2^dfChangesDuringAdaptation$popfitIntercept
dfChangesDuringAdaptation$biomassProductionum3perCellPerDayT30 <- 10^dfChangesDuringAdaptation$sizeT30 * 2^(dfChangesDuringAdaptation$popfitIntercept + dfChangesDuringAdaptation$popfitSlope * 30)

dfChangesDuringAdaptation$growthRateT0 <- dfChangesDuringAdaptation$popfitIntercept
dfChangesDuringAdaptation$growthRateT30 <- dfChangesDuringAdaptation$popfitIntercept + dfChangesDuringAdaptation$popfitSlope * 30


# Here I need to reorganise the dataframe because I want the data at time 0 and time
# 30 in the same column. I use the library reshape2 for this
if (!require(reshape2))
{install.packages("reshape2")}


# In order to use ggplot I need to reshape the dataframe
library(reshape2)
dR <- reshape(dfChangesDuringAdaptation, direction='long', 
              varying=c('biomassProductionperDayT0', 'biomassProductionperDayT30'), 
              timevar='t0ort30',
              times=c('t0', 't30'),
              v.names='biomassProdRate',
              idvar='newID')


plotBiomass <- ggplot(dR, aes(t0ort30, biomassProdRate, colour=factor(mediumConcentration + tAdaptation*100), fill=factor(tAdaptation), shape=factor(mediumConcentration), size=factor(mediumConcentration))) +
  # geom_boxplot(width=0.6) + # add a box plot
  geom_line(aes(group = newID), position=position_dodge(width=0)) +
  geom_point(size=3, position=position_dodge(width=0)) +
  theme_classic(base_size = 20) + # this defines the "theme" and the font size, other options for theme are theme_minimal, theme_dark, etc.
  # labs(fill = "Condition:") +# this specifies a custom legend
  theme(legend.position="none") + 
  scale_y_continuous(name=expression(paste(Delta, ' Biomass'," (1/day)")), limits=c(0,20)) +
  # scale_y_continuous(name=expression(paste('Biomass prod.'," (",  mu, 'm'^3, "/", mu, 'm'^3, "/day)")), limits=c(0,20)) +
  scale_x_discrete(name="Time (days)", labels=c("0", 30), expand=c(0.1,0.1)) + # in this case the x scale is not numerical
  scale_colour_manual(values=c("#8FCDDC", "#3B9AB2", "#18434E", "#F2DD70", "#EBCC2A", "#463C07", "#FF7C6C", "#F21A00", "#660B00")) +
  scale_fill_manual(values=c("#3B9AB2", "#EBCC2A", "#F21A00", "#FF00FF")) +
  scale_size_manual( values = c(0.4, 0.6, 0.8) ) +
  scale_shape_manual(values=c(21, 24, 22)) # shapes for the markers
plotBiomass


dR <- reshape(dfChangesDuringAdaptation, direction='long', 
              varying=c('biomassProductionum3perCellPerDayT0', 'biomassProductionum3perCellPerDayT30'), 
              timevar='t0ort30',
              times=c('t0', 't30'),
              v.names='biomassProdRatePerCell',
              idvar='newID')


plotBiomassPerCell <- ggplot(dR, aes(t0ort30, biomassProdRatePerCell, colour=factor(mediumConcentration + tAdaptation*100), fill=factor(tAdaptation), shape=factor(mediumConcentration), size=factor(mediumConcentration))) +
  # geom_boxplot(width=0.6) + # add a box plot
  geom_line(aes(group = newID), position=position_dodge(width=0)) +
  geom_point(size=3, position=position_dodge(width=0)) +
  theme_classic(base_size = 20) + # this defines the "theme" and the font size, other options for theme are theme_minimal, theme_dark, etc.
  # labs(fill = "Condition:") +# this specifies a custom legend
  theme(legend.position="none") + 
  scale_y_continuous(name=expression(paste(Delta, ' Biomass'," (",  mu, 'm'^3, "/cell/day)")), trans="log10", limits=c(1e4, 1e6)) +
  scale_x_discrete(name="Time (days)", labels=c("0", 30), expand=c(0.1,0.1)) + # in this case the x scale is not numerical
  scale_colour_manual(values=c("#8FCDDC", "#3B9AB2", "#18434E", "#F2DD70", "#EBCC2A", "#463C07", "#FF7C6C", "#F21A00", "#660B00")) +
  scale_fill_manual(values=c("#3B9AB2", "#EBCC2A", "#F21A00", "#FF00FF")) +
  scale_size_manual( values = c(0.4, 0.6, 0.8) ) +
  scale_shape_manual(values=c(21, 24, 22)) # shapes for the markers
plotBiomassPerCell

    
dR <- reshape(dfChangesDuringAdaptation, direction='long', 
              varying=c('growthRateT0', 'growthRateT30'), 
              timevar='t0ort30',
              times=c('t0', 't30'),
              v.names='growthRate',
              idvar='newID')

plotGrowth <- ggplot(dR, aes(t0ort30, growthRate, colour=factor(mediumConcentration + tAdaptation*100), fill=factor(tAdaptation), shape=factor(mediumConcentration), size=factor(mediumConcentration))) +
  # geom_boxplot(width=0.6) + # add a box plot
  geom_line(aes(group = newID), position=position_dodge(width=0)) +
  geom_point(size=3, position=position_dodge(width=0)) +
  theme_classic(base_size = 20) + # this defines the "theme" and the font size, other options for theme are theme_minimal, theme_dark, etc.
  # labs(fill = "Condition:") +# this specifies a custom legend
  theme(legend.position="none") + 
  scale_y_continuous(name=expression(paste('Growth'," (gen/day)")), limits=c(0,5)) +
  scale_x_discrete(name="Time (days)", labels=c("0", 30), expand=c(0.1,0.1)) + # in this case the x scale is not numerical
  scale_colour_manual(values=c("#8FCDDC", "#3B9AB2", "#18434E", "#F2DD70", "#EBCC2A", "#463C07", "#FF7C6C", "#F21A00", "#660B00")) +
  scale_fill_manual(values=c("#3B9AB2", "#EBCC2A", "#F21A00", "#FF00FF")) +
  scale_size_manual( values = c(0.4, 0.6, 0.8) ) +
  scale_shape_manual(values=c(21, 24, 22)) # shapes for the markers
plotGrowth

dR <- reshape(dfChangesDuringAdaptation, direction='long', 
              varying=c('sizeT0', 'sizeT30'), 
              timevar='t0ort30',
              times=c('t0', 't30'),
              v.names='bodySize',
              idvar='newID')

plotBodySize <- ggplot(dR, aes(t0ort30, 10^bodySize, colour=factor(mediumConcentration + tAdaptation*100), fill=factor(tAdaptation), shape=factor(mediumConcentration), size=factor(mediumConcentration))) +
  # geom_boxplot(width=0.6) + # add a box plot
  geom_line(aes(group = newID), position=position_dodge(width=0)) +
  geom_point(size=3, position=position_dodge(width=0)) +
  theme_classic(base_size = 20) + # this defines the "theme" and the font size, other options for theme are theme_minimal, theme_dark, etc.
  # labs(fill = "Condition:") +# this specifies a custom legend
  theme(legend.position="none") + 
  scale_y_continuous(name=expression(paste('Volume '," ",  mu, 'm'^3)), limits=c(0,30000)) +
  scale_x_discrete(name="Time (days)", labels=c("0", 30), expand=c(0.1,0.1)) + # in this case the x scale is not numerical
  scale_colour_manual(values=c("#8FCDDC", "#3B9AB2", "#18434E", "#F2DD70", "#EBCC2A", "#463C07", "#FF7C6C", "#F21A00", "#660B00")) +
  scale_fill_manual(values=c("#3B9AB2", "#EBCC2A", "#F21A00", "#FF00FF")) +
  scale_size_manual( values = c(0.4, 0.6, 0.8) ) +
  scale_shape_manual(values=c(21, 24, 22)) # shapes for the markers
plotBodySize

if (!require(ggpubr))
{install.packages("ggpubr")}

# Make a figure with all the plots
library(ggpubr)
# combine together the two plots into a single figure
fullFigure <- ggarrange(plotGrowth, plotBodySize, plotBiomass, plotBiomassPerCell, ncol = 2, nrow = 2, common.legend=FALSE, align = "v")
fullFigure
ggsave(file="figure_growth_and_biomass_over_time.png", dpi = 600, width = 24, height = 20, units = "cm")
ggsave(file="figure_growth_and_biomass_over_time.eps", device="eps", dpi = 1200, width = 24, height = 20, units = "cm")
library(Cairo)
ggsave(file="figure_growth_and_biomass_over_time.pdf", device=cairo_pdf, dpi = 1200, width = 24, height = 20, units = "cm")

# This shows the values of fitted slope for the population growth rate
# to see when the growth rate is decreasing or increasing over time
cbind(dfChangesDuringAdaptation$tAdapt_and_density, dfChangesDuringAdaptation$popfitSlope - dfChangesDuringAdaptation$popfitSlopeSe, dfChangesDuringAdaptation$popfitSlope + dfChangesDuringAdaptation$popfitSlopeSe)

write.table(dfChangesDuringAdaptation, "changes_of_size_and_pop_growth_during_adaptation.csv", append = FALSE, sep = ",", row.names = FALSE)

