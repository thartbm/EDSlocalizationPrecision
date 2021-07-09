
# Participants and Styles -----

getGroupParticipants <- function(group) {
  
  
  pdf <- read.csv('data/participants_files.csv',stringsAsFactors=FALSE)
  participants <- pdf$ID[pdf$folder == group]
  
  print(group)
  print(participants)
  
  return(participants)
  
}


getStyle <- function() {
  
  
  # Set styles for the two groups (colors, labels, max rotation)
  groups    =  c('sEDS', 
                 'zEDS')
  rotations =  c(30,
                 30)
  solidcolors =  c(rgb(229, 22,  54,  255, max = 255), 
                   rgb(136, 153, 255, 255, max = 255))
  
  transcolors =  c(rgb(229, 22,  54,  47,  max = 255), 
                   rgb(136, 153, 255, 47,  max = 255))
  
  linestyles = c(1,
                 1)
  labels <-    c('control',
                 'EDS')
  
  # these are not really different for the groups:
  figwidth <- c(6.5,6.5)
  pointsize <- c(1, 0.83) # the second one is used in fig 5
  fontsize  <- c(1,1)
  
  
  styles <- data.frame(groups,rotations,solidcolors,transcolors,linestyles,labels,figwidth,fontsize,pointsize)
  colnames(styles) <- c('group','rotation','color_solid','color_trans','linestyle','label','figwidth','fontsize','pointsize')
  
  return(styles)
  
}

# Descriptives and statistics -----

t.interval = function(data, variance = var(data), conf.level = 0.95) {
  
  z = qt((1 - conf.level)/2, df = length(data) - 1, lower.tail = FALSE)
  
  xbar = mean(data)
  sdx = sqrt(variance/length(data))
  
  return(c(xbar - z * sdx, xbar + z * sdx))
  
}

getConfidenceInterval <- function(data, variance = var(data), conf.level = 0.95, method='t-distr', resamples=1000, FUN=mean, returndist=FALSE) {
  
  if (method %in% c('t-distr','t')) {
    
    z = qt((1 - conf.level)/2, df = length(data) - 1, lower.tail = FALSE)
    
    xbar = mean(data)
    sdx = sqrt(variance/length(data))
    
    return(c(xbar - z * sdx, xbar + z * sdx))
    
  }
  
  # add sample z-distribution?
  
  # for bootstrapping:
  
  if (method %in% c('bootstrap','b')) {
    
    data <- data[which(is.finite(data))] #need is.finite due to NA values
    
    samplematrix <- matrix(sample(data, size = resamples*length(data), replace = TRUE), nrow = resamples)
    BS <- apply(samplematrix, c(1), FUN=FUN) 
    
    lo <- (1-conf.level)/2.
    hi <- 1 - lo
    
    if (returndist) {
      percentiles <- data.frame(percentile=seq(.01,.99,.01),value=quantile(BS, probs=seq(.01,.99,.01)))
      densdist <- density(BS, bw='SJ', from=min(percentiles$value), to=max(percentiles$value))  
      return(list('percentiles'=percentiles, 'density'=densdist, 'CI95'=quantile(BS, probs = c(lo,hi))))
    } else {
      return(quantile(BS, probs = c(lo,hi)))
    }
    
  }
  
}


etaSquaredTtest <- function(g1,g2=NA,na.rm=TRUE,mu=0) {
  
  doOneSample <- FALSE
  doTwoSample <- FALSE
  
  if (length(g2) == 1) {
    if (is.na(g2)) {
      doOneSample <- TRUE
    } else {
      # set mu to the single value in g2 and do a one sample one anyway?
    }
  } else {
    doTwoSample <- TRUE
  }
  
  if (doOneSample) {
    
    # compare group 1 mean with mu as explanation
    SStotal <- sum((g1-mean(g1,na.rm=na.rm))^2)
    SSeffect <- sum(((mean(g1, na.rm=na.rm) - mu)^2)*length(g1))
    # 
    # 
    return(SSeffect / SStotal)
    
  }
  
  if (doTwoSample) {
    
    overallmean <- mean(c(g1,g2),na.rm=na.rm)
    # compare overall mean with group means as explanation
    SStotal <- sum((c(g1,g2) - overallmean)^2, na.rm=na.rm)
    SSeffect <- sum(length(g1)*(mean(g1,na.rm=na.rm)-overallmean)^2, length(g2)*(mean(g2,na.rm=na.rm)-overallmean)^2)
    return(SSeffect / SStotal)
    
  }
  
}

# Reach processing functions -----

rotateTrajectory <- function(X,Y,angle) {
  
  # create rotation matrix to rotate the X,Y coordinates
  th <- (angle/180) * pi
  R <- t(matrix(data=c(cos(th),sin(th),-sin(th),cos(th)),nrow=2,ncol=2))
  
  # put coordinates in a matrix as well
  coordinates <- matrix(data=c(X,Y),ncol=2)
  
  # rotate the coordinates
  Rcoordinates <- coordinates %*% R
  
  # return the rotated reach
  return(Rcoordinates)
  
}

getAlignedBiases <- function(reachAngles) {
  
  # # remove outliers:
  # sets <- unique(reachAngles$set)
  # idxOK <- c()
  # 
  # # collect indices of non-outliers, separately per set:
  # for (setno in c(1:length(sets))) {
  #   
  #   set <- sets[setno]
  #   setrows <- which(reachAngles$set == set)
  #   angdevs <- reachAngles$angular_deviation[setrows]
  #   adm <- mean(angdevs, na.rm=TRUE)
  #   std3 <- sd(angdevs, na.rm=TRUE) * 3
  #   idxOK <- c(idxOK, setrows[intersect(which(angdevs > (adm-std3)),which(angdevs < (adm+std3)))])
  #   
  # }
  # 
  # # keep only the non-outliers:
  # reachAngles <- reachAngles[idxOK,]
  
  # calculate the median angular deviation per target angle
  reachAngles <- aggregate(as.formula('angular_deviation ~ target_angle'), reachAngles, median, na.rm=TRUE)
  
  # return the result
  return(reachAngles)
  
}



getTrialReachAngleAt <- function(trialdf, location='maxvel') {
  
  
  # location (string) determines where the angle of thereach is determines, it is one of:
  # maxvel: maximum velocity (default)
  # endpoint: end of the reach
  # cmX: the last sample before this distance from home, where X is replaced by a numeral
  
  # return a matrix of two numbers:
  reachangle = matrix(data=NA,nrow=1,ncol=2)
  
  # if the trial was rejected, return empty matrix now
  if (trialdf[1,'trialselected'] == 0) {
    
    return(reachangle);
    
  }
  
  # extract the relevant reach information
  X <- trialdf[trialdf$sampleselected == 1,'Xrobot']
  Y <- trialdf[trialdf$sampleselected == 1,'Yrobot']
  MV <- trialdf[trialdf$sampleselected == 1,'maxvelocity']
  angle <- trialdf[1,'target']
  
  # print(X)
  
  # rotate the trajectory
  # (this avoids problems in the output of atan2 for large angles)
  trajectory <- rotateTrajectory(X,Y,-1*angle)
  X <- trajectory[,1]
  Y <- trajectory[,2]
  
  # now try find the specified location in this reach:
  # if we can't find it, we need to know
  invalidlocation <- TRUE
  
  # maximum velocity, should be in the data
  if (location == 'maxvel') {
    rown <- which(MV == 1)
    if (length(rown) > 1) {
      rown <- rown[1]
    }
    if (length(rown) == 0) {
      # no maximum velocity defined!
      return(reachangle)
    }
    invalidlocation <- FALSE
  }
  # end point, just the last point in the selected stretch of the reach
  if (location == 'endpoint') {
    rown <- length(X)
    invalidlocation <- FALSE
  }
  # cutoff in centimers, the last sample before this cutoff distance is reached
  # this assumes that people don't go back, or that there is only one movement from home to target
  if (substring(location,1,2) == 'cm') {
    distance <- as.numeric(substring(location, 3))
    
    # get the distance from home:
    dist <- sqrt(X^2 + Y^2)
    
    # if there are no selected samples below 3 cm: return NAs
    if (length(which(dist < distance)) == 0) {
      return(reachangle)
    }
    
    # find the last sample, where dist < 3
    rown <- max(which(dist < distance))
    invalidlocation <- FALSE
  }
  
  # if we don't have a valid location, we can't calculate an angle to return
  if (invalidlocation) {
    return(reachangle)
  }
  
  # calculate the angle at that point for the rotated trajectory
  # this is the angular deviation we are looking for
  angulardeviation <- (atan2(Y[rown],X[rown]) / pi) * 180
  
  # put the result in the little matrix:
  reachangle[1,1] <- angulardeviation
  reachangle[1,2] <- angle
  
  return(reachangle)
  
}


# Convert Methods Figures -----

library('rsvg')
library(magick)

rePlotMethodsFigX <- function(target='inline', fig=1) {
  
  if (target=='svg') {
    # nothing to do here:
    return()
  }

  dpi <- c(300,1200)[fig]
  
  # read in the figure from the SVG source file:
  FigX <- magick::image_read_svg(sprintf('doc/Fig%d.svg',fig))

  # extract width and height properties:
  width <- magick::image_info(FigX)[['width']]
  height <- magick::image_info(FigX)[['height']]

  # read style settings for the project to get the standard width:
  styles <- getStyle()
  # get the output figure width
  fw <- styles$figwidth[1]
  # determine what the output figure height should be:
  fh <- fw * (height/width)


  # create the graphics device (if not inline):
  if (target == 'tiff') {
    tiff( filename=sprintf('doc/Fig%d.tiff',fig),
          width=fw*dpi, height=fh*dpi, units='px',
          type='cairo',compression='lzw',res=dpi)
  }
  if (target == 'pdf') {
    pdf( file = sprintf('doc/Fig%d.pdf',fig),
         width=fw*dpi, height=fh*dpi)
  }

  # output the figure:
  par(mai=c(0,0,0,0))
  plot(FigX)

  # close non-inline graphics device:
  if (target %in% c('tiff','pdf')) {
    dev.off()
  }

}

# Relative log-likelihoods -----

relativeLikelihood <- function(crit) {
  return( exp( ( min( crit  ) - crit  ) / 2 ) )
}