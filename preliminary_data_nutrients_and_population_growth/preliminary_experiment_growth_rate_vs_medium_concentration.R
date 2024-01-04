rm(list=ls()) # clean memory
if(!is.null(dev.list())) dev.off()

library(ggplot2)
library(Cairo)

fileName <- "~/Tetrahymena/preliminary_data_nutrients_and_population_growth/data_growth_rate_vs_medium_concentration.csv" # file.choose() # ask the user to select a file name.

dataFileFolder <- dirname(fileName) # here we ask the user to select a file with the data that we want to plot.
setwd(dirname(fileName))

# read the table of data inside the file. Note the additional parameters
# indicating that the table has a header line and that the different fields of the table
# are separated by a comma
dTot <- read.table(fileName, header=TRUE, sep=',')

# show header names
names(dTot)
allFitResults <- data.frame(temperature=integer(), v_max=double(), Km=double(), B=double())
allTemperatures <- unique(dTot$Temperature_C)
for (ttt in 1:length(allTemperatures))
{
currentTemperature = allTemperatures[ttt]
# select the data at the current temperature
d <- subset(dTot, Temperature_C ==currentTemperature & Concentration_solution < 400)

# the function that I want to fit is a Michaelis-Menten equation with some constant loss B
fMM <- function(S) v_max*S/(Km + S) - B

# I first rewrite the equation in the form of a double-reciprocal equation (Lineweaver-Burke)


nlsfit <- nls(Growth_rate ~ v_max * Concentration_solution / (Km + Concentration_solution) -B, data = d, start=list(v_max=max(d$Growth_rate),Km=30,B=1))
v_max <- coefficients(nlsfit)[1]
Km <- coefficients(nlsfit)[2]
B <- coefficients(nlsfit)[3]
allFitResults <- rbind(allFitResults, data.frame(temperature=currentTemperature, v_max=v_max, Km=Km, B=B))

lineColours = c("#6ACFE4", "#3B9AB2", "#EBCC2A", "#F21A00", "#BB2200")


# ggplot the data frame d.
g <- ggplot(d, aes(x= Concentration_solution, y= Growth_rate)) +
  # geom_vline(aes(xintercept=Km)) +
  geom_function(fun = fMM, colour = "black", size=1) +
  geom_jitter(position = position_jitter(height = 0, width = 3), size=3, fill=lineColours[ttt], colour="black", alpha=0.5, shape=21) + 
  theme_classic(base_size = 15) +
  geom_text(aes(x=150, y=4.2, label=paste("T =", currentTemperature, "°C", sep="")), size=5)+
  scale_x_continuous( name="Concentration", breaks=c(0, 50, 100, 150, 200), labels=c("0", "50%", "100%", "150%", "200%")) +
  scale_y_continuous(name="Growth rate (gen/day)") +
  coord_cartesian(ylim=c(-1, 5))
  # scale_colour_grey() + # an option to draw in greyscale
  # ggtitle(paste("Growth rate of Tetrahymena at", currentTemperature, "degrees")) +
  # geom_text(aes(x=150, y=v_max - B+ 0.2, label=paste("V max. =", round(v_max, 3), "gen / day")))+
  # geom_text(aes(x=Km+10, y=v_max+1, label=paste("Km = ", round(Km, 3))))

print(g)

ggsave(file=paste("Michaelis_Menten_T", currentTemperature, ".png", sep=""), dpi = 600, width = 10, height = 10, units = "cm")
ggsave(file=paste("Michaelis_Menten_T", currentTemperature, ".eps", sep=""), device="eps", dpi = 600, width = 10, height = 10, units = "cm")
ggsave(file=paste("Michaelis_Menten_T", currentTemperature, ".pdf", sep=""), device=cairo_pdf, dpi = 600, width = 10, height = 10, units = "cm")


}


write.table(allFitResults, "allFitResults_growth_vs_medium_concentration.csv", append = FALSE, sep = ", ", row.names = FALSE)

