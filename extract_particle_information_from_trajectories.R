# These are the columns in a trajectory file: frame,time,frame,xCircle,yCircle,x,y,orientation,mass,particle
rm(list=ls()) # clean memory

if (FALSE)
{
  require(zoo) # I use the zoo library for the lagged series, but I have to convert back and forth
  require('tidyr') # for adding NA to the missing frames
  
  showFig <- 0
  
  conversionPixelPerMicrometre <- sqrt(779441.5)/2/1000
  
  
  conversionMillilitresPerPixel <-  4 * 10^-4 / 779441.5
  # the volume in which the cells are counted depends on the area selected under the microscope in pixels (d1$AreaPixels). Each square is 0.1 microliters, or 10^-4 milliliters and 4 squares are together 779441 squared pixels: (X px) * (4 sq) * (10^-4 ml / sq) / (779441 px) = (Y ml)
  
  
  # for the framesPerSecond <- 54 used timeScaleForSpeed <- 5 and timeScaleForOrientation <- 5
  framesPerSecond <- 30
  timeScaleForSpeed <- 3 # this is the interval for calculating speed, in frames
  timeScaleForOrientation <- 3 # similar as above but for orientation. It is important to be a bit careful in this case as the orientation is used to calculate the autocorrelation time, so this time scale should not be too long
  
  
  fileName <- "/Users/perna/Dropbox/tetrahymena_results/tracking_results/speed_response_2021/test.csv" #only to pick the directory
  # file.choose() # ask the user to select a file name.
  # fileName <- "/Users/perna/Dropbox/tetrahymena_results/tracking_results/population_growth_2021/population_growth_2021_04_26/test.csv" #only to pick the directory
  # fileName <- "/Users/perna/Dropbox/tetrahymena_results/tetrahymena_density_speed/density_speed2021_54fps/test.csv"
  # fileName <- "/Users/perna/Dropbox/tetrahymena_results/tetrahymena_density_speed/density_speed2021_30fps/test.csv"
  # fileName <- "/Users/perna/Dropbox/tetrahymena_results/tetrahymena_density_speed/concentration_experiment2020/test.csv"
  
  
  setwd(dirname(fileName))
  getwd()
  
  if (file.exists("track_analysis_results_individual_particles.csv")) 
  {
    #Delete file if it exists
    file.remove("track_analysis_results_individual_particles.csv")
  }
  
  allFiles <- list.files(path = ".", pattern = "*trajectories.csv", recursive = TRUE, all.files =FALSE, ignore.case = TRUE, include.dirs = FALSE, full.names = TRUE)
  
  # data frame with the main results for all files!
  # allExperimentResults <- data.frame(fileName=allFiles, medianParticleArea = NA, medianParticleSpeed = NA, estimatedDiffusionCoefficient = NA, meanNParticlesPerFrame = NA, meanDensity = NA, nFrames = NA)
  
  firstFileIteration = 1
  for (fc in 1:length(allFiles))
  {
    # fc = 179
    fileName = allFiles[fc]
    print(paste("processing:", basename(fileName), "(n.", fc, "of", length(allFiles), ").", Sys.time()))
    # read the table of data inside the file. Notice the additional parameters
    # indicating that the table has a header line and that the different fields of the table
    # are separated by a comma
    d <- read.table(fileName, header=TRUE, sep=',')
    
    if (dim(d)[1] == 0)
    {
      print(paste("skipping ", fileName, "because there were no trajectories of sufficient length"))
      next
    }
    
    
    # # limit to the first N frames
    # if (tail(d$frame,1) > 3000)
    # {
    #   d <- d[which(d$frame <= 3000),]
    # }
    # 
    # d <- d[which(d$frame >= 20 * framesPerSecond),] # exclude first 20 seconds
    # # if (length(d) > 5000) {d <- d[1:50000,]}
    #     #{ d <- sample(d, 10000, replace = FALSE) }
    
    if (dim(d)[1] == 0)
    {
      print(paste("skipping ", fileName, "probably because video is too short"))
      next
    }
    # show header names
    names(d)
    
    # add new columns for the orientation from movement and for the autocorrelation time
    # these can be useful when then saving back the trajectories with this additional information
    d$orientationFromMovement <- d$orientation
    d$autocorrelationTimeFrame <- rep(NA, length(d$orientation))
    
    uniqueParticles <- unique(unlist(d$particle))
    # num.reps <- length(uniqueParticles)
    
    # if (length(uniqueParticles) >= 1000)
    # {
    #   print(paste("skipping ", fileName, "probably because too high density"))
    #   next
    # }
    
    # for each individual video there are two main types of data that are collected: statistics on the movement and characteristics of particles and statistics on the number of particles per frame
    allParticleResults <- data.frame(fileName = fileName, particle = uniqueParticles, medianFrame= NA, medianArea = NA, medianSpeed = NA, trajectoryLength = NA, autocorrelationTime = NA, autocorrelationTimeFrame = NA, autocorrelationTimeFrame2 = NA, autocorrelationValue = NA, meanAbsAngle = NA, minSpeed = NA, maxSpeed = NA, percentile75Speed = NA, percentile25Speed = NA, meanSpeed = NA, meanArea = NA, minArea = NA, maxArea = NA, percentile75Area = NA, percentile25Area = NA, minAEllipse = NA, maxAEllipse = NA, meanAEllipse = NA, medianAEllipse = NA, minBEllipse = NA, maxBEllipse = NA, meanBEllipse = NA, medianBEllipse = NA, medianElongation = NA)
    
    
    
    for (iii in 1:length(uniqueParticles))
    {
      # print(paste("running ", iii, " of ", length(uniqueParticles)))
      particleX <- subset(d, particle == uniqueParticles[iii])
      
      # I complete the data frame mainly because of the calculation of speed
      # If I have missing data points then I change the time lag over which I 
      # calculate speed and direction of movement
      particleX <- particleX %>% complete(frame = full_seq(frame, 1))
      
      allParticleResults$medianArea[iii] <- median(particleX$mass, na.rm=TRUE) / conversionPixelPerMicrometre^2
      allParticleResults$meanArea[iii] <- mean(particleX$mass, na.rm=TRUE) / conversionPixelPerMicrometre^2
      allParticleResults$minArea[iii] <- min(particleX$mass, na.rm=TRUE) / conversionPixelPerMicrometre^2
      allParticleResults$maxArea[iii] <- max(particleX$mass, na.rm=TRUE) / conversionPixelPerMicrometre^2
      allParticleResults$percentile75Area[iii] <- quantile(particleX$mass, probs = 0.75, na.rm=TRUE) / conversionPixelPerMicrometre^2
      allParticleResults$percentile25Area[iii] <- quantile(particleX$mass, probs = 0.25, na.rm=TRUE) / conversionPixelPerMicrometre^2
      
      allParticleResults$minAEllipse[iii] <- min(particleX$aEllipse, na.rm=TRUE) / conversionPixelPerMicrometre
      allParticleResults$maxAEllipse[iii] <- max(particleX$aEllipse, na.rm=TRUE) / conversionPixelPerMicrometre
      allParticleResults$meanAEllipse[iii] <- mean(particleX$aEllipse, na.rm=TRUE) / conversionPixelPerMicrometre
      allParticleResults$medianAEllipse[iii] <- median(particleX$aEllipse, na.rm=TRUE) / conversionPixelPerMicrometre
      
      allParticleResults$minBEllipse[iii] <- min(particleX$bEllipse, na.rm=TRUE) / conversionPixelPerMicrometre
      allParticleResults$maxBEllipse[iii] <- max(particleX$bEllipse, na.rm=TRUE) / conversionPixelPerMicrometre
      allParticleResults$meanBEllipse[iii] <- mean(particleX$bEllipse, na.rm=TRUE) / conversionPixelPerMicrometre
      allParticleResults$medianBEllipse[iii] <- median(particleX$bEllipse, na.rm=TRUE) / conversionPixelPerMicrometre
      
      allParticleResults$medianElongation[iii] <- median(particleX$bEllipse / particleX$aEllipse, na.rm=TRUE)
      
      
      particleSpeed <- sqrt((diff(particleX$x, timeScaleForSpeed))^2 + (diff(particleX$y, timeScaleForSpeed))^2) * framesPerSecond / timeScaleForSpeed / conversionPixelPerMicrometre
      
      particleOrientationFromMovement <- atan2(- diff(particleX$y, timeScaleForOrientation),diff(particleX$x, timeScaleForOrientation)) # This is the orientation estimated from the movement itself. Another orientation can be measured from the axes of the ellipse directly in each image. I change the sign of the y axis because in the image the y axis goes from top to bottom
      
      allParticleResults$medianSpeed[iii] <- median(particleSpeed, na.rm=TRUE)
      allParticleResults$minSpeed[iii] <- min(particleSpeed, na.rm=TRUE)
      allParticleResults$maxSpeed[iii] <- max(particleSpeed, na.rm=TRUE)
      allParticleResults$percentile75Speed[iii] <- quantile(particleSpeed, probs = 0.75, na.rm=TRUE)
      allParticleResults$percentile25Speed[iii] <- quantile(particleSpeed, probs = 0.25, na.rm=TRUE)
      allParticleResults$meanSpeed[iii] <- mean(particleSpeed, na.rm=TRUE)
      
      # particleOrientationFromMovement is the orientation
      # I first calculate the change of orientation OVER TIME SCALE FOR SPEED,
      # then I unwrap, I take the absolute value
      allParticleResults$meanAbsAngle[iii] <- mean( abs((pi + diff(particleOrientationFromMovement, timeScaleForSpeed)) %% (2*pi) - pi), na.rm=TRUE)
      
      allParticleResults$trajectoryLength[iii] <- nrow(particleX)
      # look at the frame in which the particles occur
      allParticleResults$medianFrame[iii]<-median(particleX$frame)
      
      
      # Calculate the directional autocorrelation
      
      allLags <- seq(0, framesPerSecond*2, by=1) # these are the lags for the autocorrelation. Explores two seconds, does not care about doing it symmetric because it is an autocorrelation
      
      
      # prepare an empty vector
      autocorrelation = rep(NA, length(allLags))
      
      # for each lag, calculate the autocorrelation
      for (lll in 1:length(allLags))
      {
        # converts all the angles to the zoo format
        orientation <- zoo(particleOrientationFromMovement)
        
        # apply the lag, fill with na the unmatched values at the beginning or at the end
        orientation_n <- stats::lag(orientation, k=allLags[lll], na.pad=T)
        
        # now converts back to numeric
        particleOrientationFromMovement_n <- as.numeric(orientation_n)
        
        # the directional autocorrelation is the cos of the difference in direction
        cosAngle <- cos(particleOrientationFromMovement[!is.na(particleOrientationFromMovement_n)] - particleOrientationFromMovement_n[!is.na(particleOrientationFromMovement_n)])
        
        # I want to focus on points when the trajectory is relatively straight
        # that is, angles smaller than about 37 degrees
        # cos(37/180*pi) = 0.799
        autocorrelation[lll] <- mean(cosAngle[cosAngle>= 0.8], na.rm=T)
        
      }
      
      # calculation of the typical time after which tetrahymena has completed one cycle of body oscillation. 
      # This is measured by first looking at the difference of orientation (cosine of the angle) between time T and time T plus or minus a lag. 
      # When tetrahymena has completed one cycle, it should have roughly the same orientation and so the cos(angle) 
      # should be close to 1. Now we want to find all the lags for which the autocorrelation is maximum again and
      # see how distant they are from each other, or better we find the lag that corresponds to the first local maximum

      localMaxima =which(diff(sign(diff(autocorrelation)))==-2)+1
      allParticleResults$autocorrelationTimeFrame2[iii] <- median(diff(allLags[localMaxima]))
      allParticleResults$autocorrelationTimeFrame[iii] <- allLags[localMaxima[1]]
      allParticleResults$autocorrelationValue[iii] <- autocorrelation[localMaxima[1]]
      allParticleResults$autocorrelationTime[iii] <- allLags[localMaxima[1]]/framesPerSecond
      
      
      # In the next few lines try to put back some of the average information into the trajectory file
      # this is mostly useful for plotting them on the video with a little bit less noise
      particleIndexInOriginalTable <- which(d$particle == uniqueParticles[iii])
      d$mass[particleIndexInOriginalTable] <- allParticleResults$medianArea[iii] * conversionPixelPerMicrometre^2
      d$aEllipse[particleIndexInOriginalTable] <- allParticleResults$medianAEllipse[iii] * conversionPixelPerMicrometre
      d$bEllipse[particleIndexInOriginalTable] <- allParticleResults$medianBEllipse[iii] * conversionPixelPerMicrometre
      d$autocorrelationTimeFrame[particleIndexInOriginalTable] <- allParticleResults$autocorrelationTimeFrame[iii]
      
      # replaces the orientation with the orientation from movement
      calculatedOrientation <- - 180/pi*c(rep(NA, ceiling(timeScaleForOrientation/2)), particleOrientationFromMovement, rep(NA, floor(timeScaleForOrientation/2))) # I minus again to go back to image coordinates y starting at the top
      # I also have to remove the NA introduced above when I filled missing particles
      d$orientationFromMovement[particleIndexInOriginalTable] <- calculatedOrientation[!is.na(particleX$particle)]
      
    }
    
    # if there is no orientation from movement pick the orientation from shape
    d$orientationFromMovement[is.na(d$orientationFromMovement)] <- d$orientation[is.na(d$orientationFromMovement)]
    
    # number of frames
    allParticleResults$nFrames <- max(d$frame)
    
    # Here compute statistics on the frames
    allParticleResults$meanNParticlesPerFrame <- nrow(d) / tail(d$frame,1) # there are n lines in the file (one for each particle) and a smaller number of frames in the video. This gives me the average number of particles per frame
    
    # read ROI file:
    getwd()
    roiFileName <- sub("_trajectories.csv$", "_ROI.txt", fileName)
    roiCoords <- read.table(roiFileName, col.names = c('x', 'y'), header=FALSE, sep=',')
    # install.packages("geometry")
    library(geometry)
    roiArea <- polyarea(roiCoords$x, roiCoords$y)
    roiVolml <- roiArea * conversionMillilitresPerPixel
    allParticleResults$meanDensity <- allParticleResults$meanNParticlesPerFrame / roiVolml
    
    # write down the results
    write.table(allParticleResults, file = sub("_trajectories.csv$", "_particles.csv", fileName), append = FALSE, sep = ", ", row.names = FALSE)
    
    if (firstFileIteration == 1)
    {
      firstFileIteration = 0
      write.table(allParticleResults, file = "track_analysis_results_individual_particles.csv", append = FALSE, sep = ", ", row.names = FALSE)
    } else{
      write.table(allParticleResults, file = "track_analysis_results_individual_particles.csv", append = TRUE, sep = ", ", row.names = FALSE, col.names=FALSE)
    }
    
    # also saves the original trajectories with associated particle information
    # write down the results
    write.table(d, file = sub("_trajectories.csv$", "_trajectories_with_aggregate_info.csv", fileName), append = FALSE, sep = ", ", row.names = FALSE)
    
    
  }
  
  plot(allParticleResults$medianSpeed, allParticleResults$autocorrelationTime)
  plot(allParticleResults$medianArea, allParticleResults$medianSpeed)
  plot(allParticleResults$medianSpeed, allParticleResults$meanAbsAngle)
  plot(sqrt(allParticleResults$medianArea), allParticleResults$medianElongation)
  plot(allParticleResults$medianArea, allParticleResults$maxSpeed)
  plot(allParticleResults$medianElongation, allParticleResults$maxSpeed)
  plot(allParticleResults$medianElongation, allParticleResults$autocorrelationTime)
  hist(allParticleResults$meanSpeed,20)
  plot(allParticleResults$medianArea, allParticleResults$percentile75Area)
  hist(allParticleResults$meanArea,20)
  
} else
{
  print("Skipping this part because it processes a lot of trajectory files but you probably don't want to")
}


