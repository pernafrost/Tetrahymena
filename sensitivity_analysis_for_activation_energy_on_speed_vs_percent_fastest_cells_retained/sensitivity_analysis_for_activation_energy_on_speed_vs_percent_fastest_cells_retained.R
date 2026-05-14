rm(list=ls()) # clean memory
if(!is.null(dev.list())) dev.off()

library(ggplot2)

setwd("/Users/perna/Library/CloudStorage/Dropbox/tetrahymena_results/sensitivity_analysis_for_activation_energy_on_speed_vs_percent_fastest_cells_retained/")


allAcuteSpeedQuantilesTogether <- read.table("allAcuteSpeedQuantilesTogether.csv", header=TRUE, sep=',', strip.white = TRUE)
allSpeedQuantilesTogether <- read.table("allSpeedQuantilesTogether.csv", header=TRUE, sep=',', strip.white = TRUE)

plotColourLevels <- paste(as.character(sort(unique(allAcuteSpeedQuantilesTogether$tAdaptation))))
plotColours <- rep(c("#3B9AB2", "#EBCC2A", "#F21A00"),5)


# plot data and model fit
ggplot(allAcuteSpeedQuantilesTogether, aes(x=speed_quantile, y=e_from_fit, color=factor(tAdaptation), size=factor(mediumConcentration), shape=factor(mediumConcentration))) +
  geom_line(linetype = "solid") +
  # geom_errorbar(aes(ymin=e_from_fit_0025, ymax=e_from_fit_0975), width=1.3) + 
  scale_color_manual(values=plotColours, limits=plotColourLevels) +
  scale_y_continuous(name="Activation energy (eV)", limits=c(0, 0.5)) +
  scale_x_reverse(name="% fastest cells retained", breaks=seq(0,100,by=20), limits=c(0,100)) +
  scale_size_manual( values = c(0.8, 1.4, 2) ) +
  theme_classic(base_size=15) + # Inverts the x-axis
  theme(legend.position = "none")


# plot data and model fit
ggplot(allSpeedQuantilesTogether, aes(x=speed_quantile, y=e_from_fit, color=factor(tAdaptation), size=factor(mediumConcentration), shape=factor(mediumConcentration))) +
  geom_line(linetype = "dashed") +
  # geom_errorbar(aes(ymin=e_from_fit_0025, ymax=e_from_fit_0975), width=1.3) + 
  scale_color_manual(values=plotColours, limits=plotColourLevels) +
  scale_y_continuous(name="Activation energy (eV)", limits=c(0, 0.5)) +
  scale_x_reverse(name="% fastest cells retained", breaks=seq(0,100,by=20), limits=c(0,100)) +
  scale_size_manual( values = c(0.8, 1.4, 2) ) +
  theme_classic(base_size=15) + # Inverts the x-axis
  theme(legend.position = "none")


# plot data and model fit
figSpeed <- ggplot(allAcuteSpeedQuantilesTogether, aes(x=speed_quantile, y=r_tref, color=factor(tAdaptation), size=factor(mediumConcentration), shape=factor(mediumConcentration))) +
  geom_line(linetype = "solid") +
  # geom_errorbar(aes(ymin=r_tref_0025, ymax=r_tref_0975), width=1.3) + 
  scale_color_manual(values=plotColours, limits=plotColourLevels) +
  scale_y_continuous(name=expression(paste("Speed at T=20°C (", mu, "m/s", ")")), limits=c(0, 900)) +
    scale_x_reverse(name="% fastest cells retained", breaks=seq(0,100,by=20), limits=c(0,100)) +
  scale_size_manual( values = c(0.8, 1.4, 2) ) +
  theme_classic(base_size=15) + # Inverts the x-axis
  theme(legend.position = "none")
figSpeed

ggsave(file="effect_of_percentile_choice_on_speed.pdf", device=cairo_pdf, dpi = 1200, width = 12, height = 12, units = "cm")


# plot data and model fit
ggplot(allSpeedQuantilesTogether, aes(x=speed_quantile, y=r_tref, color=factor(tAdaptation), size=factor(mediumConcentration), shape=factor(mediumConcentration))) +
  geom_line(linetype = "dashed") +
  # geom_errorbar(aes(ymin=r_tref_0025, ymax=r_tref_0975), width=1.3) + 
  scale_color_manual(values=plotColours, limits=plotColourLevels) +
  scale_y_continuous(name=expression(paste("Speed (", mu, "m/s", ")")), limits=c(0, 900)) +
  scale_x_reverse(name="% fastest cells retained", breaks=seq(0,100,by=20), limits=c(0,100)) +
  scale_size_manual( values = c(0.8, 1.4, 2) ) +
  theme_classic(base_size=15) + # Inverts the x-axis
  theme(legend.position = "none")




plotCounter <- 1
plotList = list()
allTAdaptation <- sort(unique(allSpeedQuantilesTogether$tAdaptation))
for (aaa in 1:length(allTAdaptation))
{
  allMediumConcentrations <- sort(unique(allSpeedQuantilesTogether$mediumConcentration[allSpeedQuantilesTogether$tAdaptation == allTAdaptation[aaa]]))
  for (mmm in 1:length(allMediumConcentrations))
  {
    print(paste("aaa=", aaa, "; mmm=", mmm))
    currentConditionAcute <-subset(allAcuteSpeedQuantilesTogether, tAdaptation == allTAdaptation[aaa] & mediumConcentration == allMediumConcentrations[mmm])
    currentConditionLongTerm <- subset(allSpeedQuantilesTogether, tAdaptation == allTAdaptation[aaa] & mediumConcentration == allMediumConcentrations[mmm])

    thisPlot <- ggplot() +
      geom_line(data=currentConditionAcute, aes(x=speed_quantile, y=e_from_fit, color=factor(tAdaptation)), linetype = "solid", size=2) +
      geom_line(data=currentConditionLongTerm, aes(x=speed_quantile, y=e_from_fit, color=factor(tAdaptation)), size=2, lty="11") +
      # geom_errorbar(data=currentConditionLongTerm, aes(x=speed_quantile, ymin=e_from_fit_0025, ymax=e_from_fit_0975, color=factor(tAdaptation)), width=1.3) + 
      # geom_errorbar(data=currentConditionAcute, aes(x=speed_quantile, ymin=e_from_fit_0025, ymax=e_from_fit_0975, color=factor(tAdaptation)), width=1.3) + 
      scale_color_manual(values=plotColours, limits=plotColourLevels) +
      scale_y_continuous(name="Activation energy (eV)", limits=c(0, 0.5)) +
      scale_x_reverse(name="% fastest cells retained", breaks=seq(0,100,by=20), limits=c(0,100)) +
      theme_classic(base_size=15) + # Inverts the x-axis
      theme(legend.position = "none")
    thisPlot
    plotList[[plotCounter]] <- thisPlot
    plotCounter <- plotCounter + 1
  }
}





# Make a figure with all the plots
library(ggpubr)
# combine together the two plots into a single figure
fullFigure <- ggarrange(plotList[[1]], plotList[[2]], plotList[[3]],
                        plotList[[4]], plotList[[5]], plotList[[6]],
                        plotList[[7]], plotList[[8]], plotList[[9]],
                        # labels = c("15-50", "15-100", "15-200", "20-50", "20-100", "20-200", "25-50", "25-100", "25-200"),
                        ncol = 3, nrow = 3, common.legend=FALSE)

fullFigure

ggsave(file="effect_of_percentile_choice_on_activation_energy.pdf", device=cairo_pdf, dpi = 1200, width = 24, height = 20, units = "cm")



#fullFullFigure <- ggarrange(figSpeed, fullFigure, ncol=1, nrow=2, common.legend=FALSE)
#fullFullFigure
