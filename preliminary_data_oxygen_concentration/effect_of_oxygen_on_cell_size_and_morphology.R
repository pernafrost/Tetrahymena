
# clean memory
rm(list=ls())
library(Cairo)
library(ggplot2)
saveFigures <- TRUE

# Read data on cell shape and size (obtained from imageJ macro)
fileName <- "~/Tetrahymena/preliminary_data_oxygen_concentration/data_cell_volume_vs_oxygen.csv"# file.choose()
setwd(dirname(fileName))

# read the file
myData <- read.table(fileName, header=TRUE, sep=",")

#look for different fields (column names) in the file
names(myData)

myData$treatment  = factor(myData$treatment,c("hypoxia","normoxia","hyperoxia"))
myData$tt = factor(myData$tt,c("15/hypo","15/norm","15/hyper","20/hypo","20/norm","20/hyper","25/hypo","25/norm","25/hyper"))

myData$elongation = myData$Major / myData$Minor
myData$estimatedlogVolume <- log10(4/3 * pi* myData$Major * myData$Minor * myData$Minor / 8)




plotG <- ggplot(myData, aes(x=factor(adaptationTemperature, levels=unique(as.numeric(adaptationTemperature))), estimatedlogVolume, colour=tt)) +
  geom_boxplot(width=0.6, outlier.shape = NA, aes(fill = factor(treatment), colour = factor(adaptationTemperature, levels=unique(as.numeric(adaptationTemperature))))) + # add a box plot
  theme_classic(base_size = 22) +
  theme(legend.position = "none") + 
  scale_y_continuous(name=expression(paste('log'[10]*'(volume)'," ",  mu, 'm'^3)), limits=c(3.75, 4.75)) +
  scale_x_discrete(name="Adaptation conditions", labels=c("15°C", "20°C", "25°C")) +
  scale_fill_manual(values=c("#DDDDDD","#A3A3A3", "#666666", "#000000")) +
  scale_color_manual(values=c("#3B9AB2", "#EBCC2A", "#F21A00", "#FF00FF")) # this is the zissou1 palette
plotG

if (saveFigures){
ggsave(file="cell_volume_vs_oxygen_concentration.eps", device="eps", dpi = 1200, width = 12, height = 10, units = "cm")
ggsave(file="cell_volume_vs_oxygen_concentration.png", dpi = 600, width = 12, height = 10, units = "cm")
ggsave(file="cell_volume_vs_oxygen_concentration.pdf", device=cairo_pdf, dpi = 1200, width = 12, height = 10, units = "cm")
}


ggplot(myData, aes(x=factor(adaptationTemperature, levels=unique(as.numeric(adaptationTemperature))), y=Major, colour=tt)) +
  geom_boxplot(width=0.6, outlier.shape = NA, aes(fill = factor(treatment), colour = factor(adaptationTemperature, levels=unique(as.numeric(adaptationTemperature))))) + # add a box plot
  theme_classic(base_size = 22) +
  theme(legend.position = "none") + 
  scale_y_continuous(name=expression(paste("Major axis length"," (",  mu, "m)")), lim=c(0,100)) + 
  scale_x_discrete(name="Adaptation conditions", labels=c("15°C", "20°C", "25°C")) +
  scale_fill_manual(values=c("#DDDDDD","#A3A3A3", "#666666", "#000000")) +
  scale_color_manual(values=c("#3B9AB2", "#EBCC2A", "#F21A00", "#FF00FF")) # this is the zissou1 palette

if (saveFigures){
  ggsave(file="cell_length_vs_oxygen_concentration.eps", device="eps", dpi = 1200, width = 12, height = 10, units = "cm")
  ggsave(file="cell_length_vs_oxygen_concentration.png", dpi = 600, width = 12, height = 10, units = "cm")
  ggsave(file="cell_length_vs_oxygen_concentration.pdf", device=cairo_pdf, dpi = 1200, width = 12, height = 10, units = "cm")
}




ggplot(myData, aes(x=factor(adaptationTemperature, levels=unique(as.numeric(adaptationTemperature))), y=Minor, colour=tt)) +
  geom_boxplot(width=0.6, outlier.shape = NA, aes(fill = factor(treatment), colour = factor(adaptationTemperature, levels=unique(as.numeric(adaptationTemperature))))) + # add a box plot
  theme_classic(base_size = 22) +
  theme(legend.position = "none") + 
  scale_y_continuous(name=expression(paste("Minor axis length"," (",  mu, "m)")), lim=c(0,50)) + 
  scale_x_discrete(name="Adaptation conditions", labels=c("15°C", "20°C", "25°C")) +
  scale_fill_manual(values=c("#DDDDDD","#A3A3A3", "#666666", "#000000")) +
  scale_color_manual(values=c("#3B9AB2", "#EBCC2A", "#F21A00", "#FF00FF")) # this is the zissou1 palette

if (saveFigures){
  ggsave(file="cell_width_vs_oxygen_concentration.eps", device="eps", dpi = 1200, width = 12, height = 10, units = "cm")
  ggsave(file="cell_width_vs_oxygen_concentration.png", dpi = 600, width = 12, height = 10, units = "cm")
  ggsave(file="cell_width_vs_oxygen_concentration.pdf", device=cairo_pdf, dpi = 1200, width = 12, height = 10, units = "cm")
}



ggplot(myData, aes(x=factor(adaptationTemperature, levels=unique(as.numeric(adaptationTemperature))), y=elongation, colour=tt)) +
  geom_boxplot(width=0.6, outlier.shape = NA, aes(fill = factor(treatment), colour = factor(adaptationTemperature, levels=unique(as.numeric(adaptationTemperature))))) + # add a box plot
  theme_classic(base_size = 22) +
  theme(legend.position = "none") + 
  scale_y_continuous(name=expression(paste("Aspect ratio")), lim=c(0,4)) +
  scale_x_discrete(name="Adaptation conditions", labels=c("15°C", "20°C", "25°C")) +
  scale_fill_manual(values=c("#DDDDDD","#A3A3A3", "#666666", "#000000")) +
  scale_color_manual(values=c("#3B9AB2", "#EBCC2A", "#F21A00", "#FF00FF")) # this is the zissou1 palette

if (saveFigures){
  ggsave(file="cell_elongation_vs_oxygen_concentration.eps", device="eps", dpi = 1200, width = 12, height = 10, units = "cm")
  ggsave(file="cell_elongation_vs_oxygen_concentration.png", dpi = 600, width = 12, height = 10, units = "cm")
  ggsave(file="cell_elongation_vs_oxygen_concentration.pdf", device=cairo_pdf, dpi = 1200, width = 12, height = 10, units = "cm")
}


