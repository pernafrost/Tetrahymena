rm(list=ls()) # clean memory

# file name file.choose() # 
fileName <- "~/Tetrahymena/population_growth_during_adaptation/growth_data_during_adaptation.csv"

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
populationGrowthData <- cbind(populationGrowthData, newColumns, deparse.level = 1, stringsAsFactors = default.stringsAsFactors())

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
populationGrowthData$dateNumber <- as.numeric(populationGrowthData$date)/86400 # convert to units of days
populationGrowthData$dateNumber <- populationGrowthData$dateNumber - min(populationGrowthData$dateNumber)

# isolate the data for the mother culture (before splitting)
motherCulture <- subset(populationGrowthData, isMotherCulture == 1)

library(ggplot2)
ggplot(motherCulture, aes(x=dateNumber, y=growth_rate_gen_per_day, shape=line)) +
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
  
  plotG1 <- ggplot(populationGrowthData, aes(x=dateNumber, y=growth_rate_gen_per_day, shape=factor(as.numeric(density)), color=tAdapt)) +
    geom_point(size=2) + 
    geom_smooth(method=lm, formula='y~x', se=FALSE, fullrange=FALSE, aes(fill=factor(tAdapt_and_density), linetype=factor(density))) + 
    theme_classic(base_size = 22) +
    scale_y_continuous(name="generations per day", limits=c(0,5.7), breaks=seq(0,5, by=1)) +
    scale_x_continuous(name="adaptation time (days)") + 
    theme(legend.position = "none") + 
    scale_color_manual(values=plotColours) + # this is the zissou1 palette
    scale_fill_manual(values= alpha(plotColours), 0.9)# this is the zissou1 palette
  plotG1
  
  plotG1 <- ggplot(populationGrowthData, aes(x=dateNumber, y=growth_rate_gen_per_day, shape=factor(as.numeric(density)), color=tAdapt, fill=tAdapt)) +
    geom_point(size=2, alpha=0.4, colour="black") + 
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
  
  ggsave(file="tetrahymena_growth_rate_during_adaptation.png", dpi = 600, width = 15, height = 12, units = "cm")
  ggsave(file="tetrahymena_growth_rate_during_adaptation.eps", device="eps", dpi = 1200, width = 15, height = 10, units = "cm")
  
  
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
  
  
