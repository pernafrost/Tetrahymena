rm(list = ls()) # clear history

verbose = FALSE # This tells the programme to print some information for each function


#################################################################################################################
# FUNCTION respiration_rate_to_W
# Converts the respiration rate of the culture or of the cell to metabolic rate (in Watts)
# input:
# respiration rate of the organism or of the culture in micromol / L / min
# dPredator is the density of the predator (tetrahymena) in cells / ml
# output:
# metabolic power consumption (W) of one single cell
# To make a test let's use some values from Katsu-Kimura, who
# report an oxygen consumption for a single cell of paramecium
# of 0.19-4.4 * 10^-9 L/cell/hour and transform the number to get
# 0.38-8.8 * 10^-5 J/cell/hour. I convert L to micromol by considering
# 1L * 1atm = n * R * T with T = 298 and R = 0.082057366080960
# oxygenConsumption_in_mol_per_min <- 0.19e-9 * 1/(293*0.082057366080960)*60
# [1] 4.741543e-10 mol/L/min
# oxygenConsumption <- oxygenConsumption_in_mol_per_min*10^6
# > respiration_rate_to_W(oxygenConsumption)
# [1] 3.781983e-06 # W/cell
# > respiration_rate_to_W(oxygenConsumption) * 3600
# [1] 0.01361514 # J/cell/h
respiration_rate_to_W <- function(oxygenConsumption, dPredator) {
  if (missing(dPredator))
  {
    dPredator <- 1e-3  # if dPredator is omitted, we pretend that there is only one cell per L
    # so that the function calculates the conversion per cell, rather than per L
    if (verbose)
    {
      print(paste("In function respiration_rate_to_W, calculating conversion of respiration rate per cell to W /cell"))
    }
  }
  # Convert rate of respiration in umol /min / L to rate of energy consumption in W/cell using the following assumptions:
  # Aerobic respiration goes as 1G + 6O2 = 6CO2 + 6H2O (1G is one mol of glucose) and produces 2871458 J (about 686 kCal)
  JoulesPerMol <- 2871458/6 # (J/mol) there are 2871458 J in one mol of glucose and 2871458/6 J in one mol of produced O2
  JoulesPerMicroMol <- JoulesPerMol * 1e-6 # (J/umol)
  metabolicRateJoulesPerMinPerL <- oxygenConsumption * JoulesPerMicroMol # Joules / L / min
  # I divide this metabolic rate by the number of cells of predator per litre of suspension
  # to get the metabolic rate per cell
  metabolicRatePerCellPerMin <- metabolicRateJoulesPerMinPerL / (dPredator * 1e3) #  Joules / min /cell
  metabolicRateW <- metabolicRatePerCellPerMin/60 # (Joules / s / cell) metabolicRatePerCellPerMin was in Joules per min, 
  return(metabolicRateW)
}



########################################################
## Packages and libraries                              #
########################################################

# Function to Install and Load R Packages
Install_And_Load <- function(Required_Packages)
{
  Remaining_Packages <- Required_Packages[!(Required_Packages %in% installed.packages()[,"Package"])];
  
  if(length(Remaining_Packages)) 
  {
    install.packages(Remaining_Packages);
  }
  for(package_name in Required_Packages)
  {
    library(package_name,character.only=TRUE,quietly=TRUE);
  }
}

# Specify the list of required packages to be installed and load    
Required_Packages=c("readxl", "lme4", "ggplot2",  "plyr", "dplyr", "reshape", "reshape2", "ggforce", "stringr");

# Call the Function
Install_And_Load(Required_Packages);

# function for calculating regression coefficents 
sumfun <- function(x) c(coef(x),summary(x)$r.squared, pVal <- anova(x)$'Pr(>F)'[1]) 

########################################################
# SET PATH TO THE DATA ON YOUR DEVICE     
########################################################

# aFile<-file.choose() # CHOOSE ONE FILE FROM FOLDER IN ORDER TO SET YOUR DATA PATH 
aFile <- "~/post-data/Respiration_20c_25c_tested_at_12.5c_SDR856_Oxygen.xlsx"
data.path = dirname(aFile)
setwd(data.path) 

#################################################
# ADD LIST OF ALL FILE NAMES TO BE ANALYISID HERE
#################################################


all.metabolism.files<-list.files(path = data.path,pattern = "*Oxygen.xlsx",recursive = FALSE) # LIST SENSORDISH DATA FILES 
# all.sampleID.files<-list.files(path = data.path,pattern = "*Sample.csv",recursive = FALSE)
all.sampleID.files <- str_replace(all.metabolism.files, "Oxygen.xlsx", "Sample.csv")


#################################################
# run analysis
#################################################

datalist = list()   

for (fff in 1:length(all.metabolism.files)){
  print(paste("file running ", fff, " of ", length(all.metabolism.files)))
  current.metabolism.file <- all.metabolism.files[fff]
  df.wide <- read_excel(paste(data.path, current.metabolism.file,  sep="/"), skip = 12)
  
  #convert file to long format
  df.wide <- as.data.frame(df.wide)
  df.long <- melt.data.frame(df.wide, id.vars = c(2,28,31), measure.vars = c(3:26)) # C
  names(df.long)[c(1,2, 3,4, 5)] <- c("time", "temp", "internal.temp", "position", "o2.umol.l") 
  
  # merge with sampleID file
  sampleID <- read.csv(paste(data.path, all.sampleID.files[fff], sep="/"))
  df.long <- merge(df.long, sampleID, by = "position")

  # merge with metafile
  df.long$file <- rep(current.metabolism.file, length(df.long$sample))

  df.long$time.m <- (df.long$time)# convert time to minutes if needed
  df.long$o2.umol.l <- as.numeric(df.long$o2.umol.l)
   
  # add treatment info
  df.long <- cbind(df.long, colsplit(df.long$sample, "_", c("concentration", "treatment", "replicate"))) 
  df.long$treatment <- as.character(df.long$treatment)
  
  datalist[[fff]] <- df.long 
  
}

#################################################
# combine all data frames
#################################################

df.all <- dplyr::bind_rows(datalist)
  
#################################################
# subtract control values from each datum
#################################################
  
df.all <- subset(df.all, o2.umol.l != "NA")

df.all<- ddply(df.all, .(temp, concentration, time.m), transform,
                  control = mean(o2.umol.l[treatment == "control"]))
  
df.all$o2.umol.l.corrected <- df.all$o2.umol.l - df.all$control

#################################################
# analyse each temperature seperately
#################################################

datalist = list()  

uniq <- unique(unlist(df.all$temp))

for (iii in 1:length(uniq)){
  print(paste("running ", iii, " of ", length(uniq)))
  df <- subset(df.all, temp == uniq[iii])

  
  #################################################
  # SET VALUES FOR ANALYSIS
  #################################################
  
  use.imposed.start = T # CHOOSE WHETHER TO USE A FIXED START (=F) POINT OR VARY BY RUN (=T)
  start.time <- 20 # CHOOSE FIXED START POINT
  incubation.time <- 180 # CHOOSE LENGTH OF INCUBATIOSN TO QUANTIFY METABOLISM
  
   if (use.imposed.start == TRUE)
   {    start.time <- df$start 
          end.time <-  df$end
   }
  
   df.sub <- subset (df, time.m > start.time & time.m < end.time) 

   # plot raw data and analysis subset
   pdf(paste(data.path, uniq[iii], '.rates.pdf', sep = ''), width = 15, height = 20)   
  p <- ggplot(data = df, aes(y=o2.umol.l.corrected, x= time.m)) + geom_point(colour = "grey") +
    facet_wrap(~sample, scales = "free", ncol = 6, nrow = 8) +
    geom_smooth(data = df.sub, aes(y=o2.umol.l.corrected, x= time.m), method = "lm", colour = "red", se = FALSE, formula="y~x")
  plot(p)
  dev.off()
  
  pdf(paste(data.path, uniq[iii], '.internal.temp.pdf', sep = ''), width = 5, height = 5)   
  p <- ggplot(data = df, aes(y=internal.temp, x= time.m)) + geom_point(colour = "grey") 
  plot(p)
  dev.off()
  
 

  # # get rates from linear model
  df.sub <- subset(df.sub, o2.umol.l != "No Sensor") # remove any empty locations
  df.sub$code <- paste(df.sub$sample,  df.sub$density, sep = "_") # add density information
  
  df.list <- lmList(o2.umol.l.corrected ~ time.m  |code, data = df.sub)
  is.na(df.list) <- 0
  df.sum <- as.data.frame(t(sapply(df.list,sumfun)))
  df.sum$code <- as.factor(rownames(df.sum))
  names(df.sum)[c(5,1,2,3,4)] <- c("code", "intercept", "slope", "r2", "p.val")
  df.sum$temperature <- df.sub$temp[iii] #   Add back coding
  
  # add interation to the list
  datalist[[iii]] <- df.sum 

}

#################################################
# combine all data frames
#################################################

df.rates <- dplyr::bind_rows(datalist)

df.rates <- cbind(df.rates, colsplit(df.rates$code, "_", c("concentration", "treatment", "replicate", "density")))

df.rates$rate.per.cell <- df.rates$slope / (df.rates$density *1000) # density corrected rates


# Now convert the respiration rate from umol / L / min to nanoWatts
df.rates$rate.per.cell.nW <- -respiration_rate_to_W(df.rates$rate.per.cell) * 10^9

################################################
# EXPORT DATAFRAME 
################################################

# exports data file
write.csv(df.rates, file = "summary_MR_post_adaptation.csv", row.names = FALSE)



################################################
# END
################################################


# ANCOVA

df.rates$treatment <- as.factor(df.rates$treatment)
df.rates$concentration <- as.factor(df.rates$concentration)
df.rates$logRate <- log10(-(df.rates$rate.per.cell*60))

df.rates$kT <- 1 / (0.00008617343*(273.15+df.rates$temperature))
df.rates$ktc <-  (1/(0.0000862*293.15))- df.rates$kT # 1.ktc-1.kt set at 20c

m <- lmer(logRate ~  ktc + treatment + concentration + ktc:treatment + ktc:concentration+ (1|replicate), data = df.rates)
anova(m)
visreg(m, "ktc", by = "treatment")


ggplot(df.rates,aes(x=ktc,y=log10(-rate.per.cell),color=replicate))+geom_point(size = 3) +
  stat_smooth(method="lm", se=F, colour = "black", size = 1) +
  facet_grid(concentration~treatment, scales = "free_y")  #theme(legend.position = "none") 

df.rates.100 <- subset(df.rates, concentration == "100")
m <- lm(logRate ~  ktc + treatment + ktc:treatment,  data = df.rates.100)
anova(m)


m <- lmer(logRate ~  ktc + treatment + ktc:treatment + (1|replicate),  data = df.rates.100)
anova(m)
visreg(m, "ktc", by = "treatment")


#####################
# plot aesthetics 
#####################

txt.a = element_text (size = 18)
txt.b = element_text (size = 14)

pub.tweeks <- theme_bw() + theme(panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank(),
                                 strip.background = element_blank(),
                                 panel.border = element_rect(colour = "black"),
                                 strip.text.x = element_text(size = 16),
                                 axis.title.x = txt.a, axis.title.y = txt.a,
                                 axis.text.x = txt.b, axis.text.y = txt.b) 


########################################################

df.plot <- subset(df.rates, logRate >=-14 & is.na(str_match(treatment, "control")))

pdf("AE.pdf", width = 10, height = 5) 
ggplot(df.plot,aes(x=ktc,y=logRate))+ geom_point(aes(fill = factor(treatment)), pch=21, size=5) + 
  stat_smooth(method="lm", se=F, colour = "black", size = 1) +
  facet_grid(~treatment, scales = "free_y")  + pub.tweeks +
  scale_fill_manual(values=c("#3B9AB2", "#EBCC2A", "#F21A00")) + # this is the zissou1 palette
  theme(legend.position = "none") +
  xlab(expression(paste(Standardised~temperature:~1/kT[c] - 1/kT~(1/eV)))) + 
  # ylab(expression(paste(Respiration~(log~(O[2]~cell^-1~hour^-1)))))
  ylab(expression(paste("Respiration log10(",  mu, "mol[O2]/cell/hour", ")")))
dev.off()


