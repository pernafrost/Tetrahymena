rm(list=ls()) # clean memory
if(!is.null(dev.list())) dev.off()
saveFigures <- TRUE

library(ggplot2)
library(dplyr)

# NOTE that activation energies can be calculated for each line, or for each 
# experimental condition: there are 9 experimental conditions, 3 for temperature
# and 3 for medium concentration, but there are 4 replicate lines in each condition
# so the data can be either 9 values of activation energy, or 9x4=36. The files 
# related to population growth are named as follows:
# allFitResults_Schoolfield_on_short_term_population_growth.csv (36 values of activation energy)
# and 
# allFitResults_Schoolfield_on_short_term_population_growth_all.csv (9 values)
# These files are produced by the script:
# analysis_tetrahymena_population_growth.R
#
# The files for body size, or more precisely for 1/V where V is the volume of 
# the cell are the following ones:
# allFitResults_Schoolfield_on_one_over_Volume_long_term_response.csv (values of activation energy
# for each line, 36 rows in the dataset)
# and
# allFitResults_Schoolfield_on_one_over_Volume_long_term_response_all.csv (9 values, lines
# are combined together).
# These files are produced by the script:
# fit_schoolfield_equation_to_long_term_body_size.R


# the file with the activation energies for 1/(cell volume)
# fileNameCellVolume <- "~/Tetrahymena/comparison_of_activation_energy_for_cell_size_and_population_growth/allFitResults_Schoolfield_on_one_over_Volume_long_term_response_all.csv" # file.choose() # ask the user to select a file name.
fileNameCellVolume <- "~/Tetrahymena/comparison_of_activation_energy_for_cell_size_and_population_growth/allFitResults_Schoolfield_on_one_over_Volume_long_term_response.csv" # file.choose() # ask the user to select a file name.


# the file with the activation energies for population growth
# fileNamePopulationGrowth <- "~/Tetrahymena/comparison_of_activation_energy_for_cell_size_and_population_growth/allFitResults_Schoolfield_on_short_term_population_growth_all.csv" # file.choose() # ask the user to select a file name.
fileNamePopulationGrowth <- "~/Tetrahymena/comparison_of_activation_energy_for_cell_size_and_population_growth/allFitResults_Schoolfield_on_short_term_population_growth.csv" # file.choose() # ask the user to select a file name.
setwd(dirname(fileNamePopulationGrowth))

dVolume <- read.table(file = fileNameCellVolume, sep = ",", header=TRUE, stringsAsFactors = FALSE, strip.white=TRUE, na.strings = c("NA", " NA"))
dPopGrowth <- read.table(file = fileNamePopulationGrowth, sep = ",", header=TRUE, stringsAsFactors=FALSE, strip.white=TRUE, na.strings = c("NA", " NA"))
names(dPopGrowth)
names(dVolume)

# merge the two data frames
dPopGrowthAndVolume <- merge(dPopGrowth, dVolume, by = c("tAdaptation", "mediumConcentration", "line"))

ggplot(dPopGrowthAndVolume, aes(x=e_from_fit.x, y=e_from_fit.y, color=as.factor(tAdaptation))) +
  geom_point() +
  geom_errorbarh(aes(xmin=e_from_fit.x - sd_e_from_fit.x, xmax = e_from_fit.x + sd_e_from_fit.x)) +
  geom_errorbar(aes(ymin=e_from_fit.y - sd_e_from_fit.y, ymax = e_from_fit.y + sd_e_from_fit.y)) +
  geom_abline(intercept = 0, slope = 1, color="black") +
  # geom_jitter(position = position_jitter(height = 0, width = .5)) +
  # geom_smooth(method="lm", se=TRUE, fill="red", formula=y ~ x, na.rm=TRUE) +
  scale_y_continuous(name="Activation Energy (1/V)") +
  scale_x_continuous(name="Activation Energy (population growth)") +
  coord_cartesian(xlim=c(0,1.2), ylim=c(0,1.2), expand=FALSE) +
  theme_classic(base_size = 15) +
  theme(legend.position = "none") +
  ggtitle("Comparison of activation energies") +
  scale_color_manual(values=c("#3B9AB2", "#EBCC2A", "#F21A00")) # this is the zissou1 palette
# The black line is not a fit, but it is the line that we would have if all the 
# population growth was explained by cell size only


if (saveFigures){
    ggsave(file="activation_energy_population_growth_cell_volume.png", dpi = 600, width = 12, height = 10, units = "cm")
}
  
  
  ggplot(dPopGrowthAndVolume, aes(x=e_from_fit.x, y=e_from_fit.y)) +
    geom_point() +
    geom_errorbarh(aes(xmin=e_from_fit.x - sd_e_from_fit.x, xmax = e_from_fit.x + sd_e_from_fit.x)) +
    geom_errorbar(aes(ymin=e_from_fit.y - sd_e_from_fit.y, ymax = e_from_fit.y + sd_e_from_fit.y)) +
    geom_smooth(method="lm", se=FALSE, color="red", formula=y ~ x + 0, na.rm=TRUE) +
    scale_y_continuous(name="Activation Energy (1/V)") +
    scale_x_continuous(name="Activation Energy (population growth)") +
    coord_cartesian(xlim=c(0,1), ylim=c(0,1), expand=FALSE) +
    theme_classic(base_size = 15) +
    theme(legend.position = "none") +
    ggtitle("Comparison of activation energies") +
    scale_color_manual(values=c("#3B9AB2", "#EBCC2A", "#F21A00")) # this is the zissou1 palette
  
  
  
  if (saveFigures){
    ggsave(file="activation_energy_population_growth_cell_volume_with_fit.png", dpi = 600, width = 12, height = 10, units = "cm")
  }
  
  
  # I can fit a straight line
  fitVolumeVsPopGrowth <- lm(formula=e_from_fit.y ~ e_from_fit.x + 0, data = dPopGrowthAndVolume)
  print(fitVolumeVsPopGrowth)
  
  # also try to avoid regression dilution because of error on x as well as on y
  # library("lmodel2")
  # fitVolumeVsPopGrowth <- lmodel2(formula=e_from_fit.y ~ e_from_fit.x , data = dPopGrowthAndVolume, range.y="interval", range.x = "interval", nperm=99)
  # print(fitVolumeVsPopGrowth)
  
  # or I can simply compare the ratios between the two measures of activation energy
  meanRatioEPopGrowthOverEVolume <- mean(dPopGrowthAndVolume$e_from_fit.x/dPopGrowthAndVolume$e_from_fit.y)
  print(meanRatioEPopGrowthOverEVolume)
  meanRatioEVolumeOverEPopGrowth <- 1/meanRatioEPopGrowthOverEVolume
  print(meanRatioEVolumeOverEPopGrowth)
  
  # Here I create a variable with the difference of activation energies:
  dPopGrowthAndVolume$differenceActivationEnergy <- dPopGrowthAndVolume$e_from_fit.x - dPopGrowthAndVolume$e_from_fit.y
  
  # I check that the data are normally distributed (check with Shapiro test)
  # given that we focus on the difference between the two conditions,
  # we only need to test that the difference is normally distributed;
  # we do not need to test the normality of the individual distributions
  
  shapiro.test(dPopGrowthAndVolume$differenceActivationEnergy)
  # Shapiro-Wilk normality test
  #
  # data:  dPopGrowthAndVolume$differenceActivationEnergy
  # W = 0.96807, p-value = 0.3752
  
  
  # As everything is good, I can compute the one sample t-test on the difference
  t.test(dPopGrowthAndVolume$differenceActivationEnergy, mu=0)
  
  # In R I can also obtain the same result by calling the t.test with the option "paired"
  t.test(dPopGrowthAndVolume$e_from_fit.x, dPopGrowthAndVolume$e_from_fit.y, paired=TRUE)
  # And similarly with the other value of activation energy returned by the rTPC library
  t.test(dPopGrowthAndVolume$e.x, dPopGrowthAndVolume$e.y, paired=TRUE)

  
  
