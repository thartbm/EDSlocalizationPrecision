source('R/common.R')


# Data Handling -----

getLocalization <- function(groups=c('sEDS', 'zEDS')) {
  
  # use all groups:
  # groups <- c('30implicit', '30explicit', 'cursorjump', '60implicit', '60explicit')
  
  # loop through groups:
  for (group in groups) {
    
    # we store stuff temporarily:
    grouplocalization <- getGroupLocalization(group)
    
    #save(grouplocalization,file=sprintf('../data/%s_localization.dat',group))
    write.csv(grouplocalization,file=sprintf('data/%s_localization.csv',group),row.names=FALSE,quote=FALSE)
  }
  
  # success?
  return(TRUE)
  
}



getGroupLocalization <- function(group) {
  
  # This should return a dataframe with columns:
  # 1) participant
  # 2) aligned / rotated
  # 3) active / passive
  # 4) block
  # 5) trial
  # 6) arc position
  # 7) hand angle
  # 8) tap angle
  
  # get all participants in this group:
  participants <- getGroupParticipants(group)
  
  # create empty dataframe to store all info:
  alldata <- data.frame()
  
  # loop through participants:
  for (pp.id in participants) {
    
    # get localization data for participant:
    df <- getParticipantLocalization(pp.id)
    
    # concatenate all dataframes into one:
    if (prod(dim(alldata)) == 0) {
      alldata <- df
    } else {
      alldata <- rbind(alldata, df)
    }
    
  }
  
  # return dataframe with localization info for all participants in the group:
  return(alldata)
  
}



getParticipantLocalization <- function(participant) {
  
  # get info on filenames from this list:
  pdf <- read.csv('data/participants_files.csv',stringsAsFactors=FALSE)
  
  # all the participants data will be somewhere in this folder
  folder <- pdf[pdf$ID==participant, 'folder']
  
  # we'll loop through combinations of these two, to get all files:
  conditions <- c('aligned','rotated')
  reachtypes <- c('active','passive')
  
  # create empty data frame:
  participantdf <- NA
  
  # start looping:
  for (rotated in c(1:2)) {
    
    condition <- conditions[rotated]
    
    for (passive in c(1:2)) {
      
      reachtype <- reachtypes[passive]
      
      # get variables necessary to open the data
      column <- sprintf('%s_%stap',condition,reachtype)
      filename <- pdf[pdf$ID==participant, column]
      path2file <- sprintf('data/%s/%s/%s', folder, participant, filename)
      
      # read the data, but skip the header line... too many weird characters:
      sourcedf <- read.table(path2file,header=TRUE,stringsAsFactors=FALSE)
      sourcedf <- sourcedf[which(sourcedf$selected == 1),]

      # get hand coordinates at the end of the reaches, relative to the home position:
      handX <- sourcedf$handx_cm
      handY <- sourcedf$handy_cm
      
      # rotate by -90 degrees, to avoid atan2 problems around 180/-180
      coords <- rotateTrajectory(handX, handY, -90)
      handX <- coords[,1]
      handY <- coords[,2]
      
      # calculate angles (then add 90 degrees again):
      handangles <- (atan2(handY, handX) / pi) * 180
      handangles <- handangles + 90
      
      # get the tap coordinates, relative to the perceived home position
      tapX <- sourcedf$tapx_cm
      tapY <- sourcedf$tapy_cm
      
      # correct with circle fitting...
      cfit <- circleFit(tapX, tapY)
      #print(cfit)
      tapX <- tapX - cfit$x
      tapY <- tapY - cfit$y
      
      # rotate by -90 degrees to avoid atan2 problesm around 180/-180
      coords <- rotateTrajectory(tapX, tapY, -90)
      tapX <- coords[,1]
      tapY <- coords[,2]
      
      # get angles and add the 90 degrees back in:
      tapangles <- (atan2(tapY, tapX) / pi) * 180
      tapangles <- tapangles + 90
      
      # prepare a dataframe with all collected data
      # with the same numbr of rows as the source data
      nrows <- nrow(sourcedf)
      
      # get other sequences ready:
      participant_col <- rep(participant,nrows)
      #print(length(participant_col))
      rotated_col <- seq(rotated-1,rotated-1,length.out=nrows)
      #print(length(rotated_col))
      passive_col <- seq(passive-1,passive-1,length.out=nrows)
      #print(length(passive_col))
      block_col <- sourcedf$block
      #print(length(block_col))
      trial_col <- sourcedf$trial
      #print(length(trial_col))
      
      # put everything in dataframe
      targetdf <- data.frame(participant_col,rotated_col,passive_col,block_col,trial_col,sourcedf$arcangle_deg,handangles,tapangles)
      # give the right names to the columns
      colnames(targetdf) <- c('participant','rotated','passive','block','trial','arc','hand_angle','tap_angle')
      
      # concatenate all dataframes into one for this participant:
      if (is.data.frame(participantdf)) {
        participantdf <- rbind(participantdf, targetdf)
      } else {
        participantdf <- targetdf
      }
      
    }
    
  }
  
  #print(str(participantdf))
  return(participantdf)
  
}

# Circle fitting -----
# copied from SMCL package (Oct 14th, 2019)

circleFit <- function(x, y, radius=12, verbosity=0) {
  
  coords  <- data.frame(x, y)
  
  if ('optimx' %in% installed.packages()) {
    
    library(optimx)
    
    lower <- c(min(coords$x)-radius,min(coords$x)-radius)
    upper <- c(max(coords$x)+radius,max(coords$y)+radius)
    
    circlefit <- optimx(par=c('x'=0, 'y'=0), circleFitError, gr = NULL, method='L-BFGS-B', lower=lower, upper=upper, coords=coords, radius=radius)
    
    return(list('x'=circlefit$x, 'y'=circlefit$y))
    
  } else {
    
    if (verbosity > 0) cat('optimx not installed, falling back on optim\n')
    
    circlefit <-  optim(par = c('x'=0, 'y'=0), circleFitError, gr = NULL, coords=coords, radius=radius)
    
    return(list('x'=circlefit$par['x'],'y'=circlefit$par['y']))
    
  }
  
}


circleFitError <- function(par, coords, radius){
  
  return(mean((sqrt((coords$x - par['x'])^2 + (coords$y - par['y'])^2) - radius)^2, na.rm=TRUE))
  
}

# Data Handling (two) -----

getLocalizationForANOVA <- function(groups=c('sEDS','zEDS')) {
  
  # exploratory function, just to see what the data looks like
  
  # interpoints <- seq(30,135,15) + (15/2)
  #interpoints <- seq(40,140,20)
  interpoints <- c(45,90,135)
  nreps <- length(interpoints)
  
  for (group in groups) {
    
    df <- read.csv(sprintf('data/%s_localization.csv',group),stringsAsFactors=FALSE)

    # throw out rows where either hand or tap is missing
    df <- df[is.finite(df$hand_angle) & is.finite(df$tap_angle),]
    
    # collect all data for the group in one dataframe:
    outdf <- data.frame()
    
    participants <- c(unique(df$participant))
    
    
    for (reachtype in c(0,1)) {
      
      for (condition in c(0,1)) {
        
        subdf <- df[(df$rotated == condition) & (df$passive == reachtype),]
        
        # and get a list of named lists of hand positions and localization biases premade
        # we can use this to bootstrap over participants
        
        bias <- list()
        hand <- list()
        for (p in participants) {
          
          tap  <- subdf[subdf$participant == p, 'tap_angle']
          hand <- subdf[subdf$participant == p, 'hand_angle']
          bias <- tap - hand
          
          # select non-outliers:
          bias.sd <- sd(bias)
          idx <- intersect(which(bias > (mean(bias) - (3 * bias.sd))), which(bias < (mean(bias) + (3 * bias.sd))))
          
          # remove outliers for every participant?
          bias <- bias[idx]
          hand <- hand[idx]
          
          kreg <- kernelRegression(x=hand,y=bias,width=15,interpoints=interpoints)
          
          tempdf <- data.frame(rep(group,nreps), rep(p,nreps), rep(condition,nreps), rep(reachtype,nreps), interpoints, kreg)
          colnames(tempdf) <- c('group', 'participant', 'rotated_b', 'passive_b', 'handangle_deg', 'bias_deg')
          
          # add this to group data frame
          outdf <- rbind(outdf, tempdf)
          
        } # end loop over participants
        
      } # end condition loop
      
    } # end reachtype loop
    
    
    
    # save the dataframe for the group to a CSV:
    
    filename = sprintf('data/%s_loc_AOV.csv',group)
    write.csv(outdf, file=filename, row.names=FALSE, quote=FALSE)
    
  } # end group loop
  
} # end function


kernelRegression <- function(x,y,width,interpoints) {
  
  output <- c()
  
  for (interpoint in interpoints) {
    
    w <- exp(-1 * ((x-interpoint)^2/(2 * width^2)))
    output <- c(output, sum(w*y)/sum(w))
    
  }
  
  return(output)
  
}

getLocalizationTdistributionConfidenceIntervals <- function(groups=c('sEDS','zEDS')) {
  
  # interpoints <- seq(0,180,1)
  interpoints <- seq(30,150,1)
  # iterations <- 10000 # ideally 10000 or more? no more iterations: use t-distribution
  percentiles <- c(0.025, 0.5, 0.975) # this is not actually used...
  
  for (group in groups) {
    
    df <- read.csv(sprintf('data/%s_localization.csv',group),stringsAsFactors=FALSE)
    
    # correct using 3rd order polynomial
    # and removing taps where the hand location is further than 20 degrees from the arc centre
    #df <- correctLocalizationsP3(df)
    # throw out rows where either hand or tap is missing
    df <- df[is.finite(df$hand_angle) & is.finite(df$tap_angle),]
    
    participants <- c(unique(df$participant))
    
    KRbiases <- array(NA, dim=c(2, 2, length(interpoints), length(participants)))
    
    for (reachtype in c(0,1)) {
      
      for (condition in c(0,1)) {
        
        print(sprintf('%s %d %d',group,reachtype,condition))
        
        # let's subset the dataframe to speed things up
        subdf <- df[(df$rotated == condition) & (df$passive == reachtype),]
        
        # and get a list of named lists of hand positions and localization biases premade
        # we can use this to bootstrap over participants
        # participants <- c(unique(subdf$participant))
        # pp.bias <- list()
        # pp.hand <- list()
        
        for (ppno in c(1:length(participants))) {
          
          p <- participants[ppno]
          
          tap  <- subdf[subdf$participant == p, 'tap_angle']
          hand <- subdf[subdf$participant == p, 'hand_angle']
          bias <- tap - hand
          
          # select non-outliers:
          bias.sd <- sd(bias)
          idx <- intersect(which(bias > (mean(bias) - (3 * bias.sd))), which(bias < (mean(bias) + (3 * bias.sd))))
          
          # remove outliers for every participant?
          bias <- bias[idx]
          hand <- hand[idx]
          
          # collect all the kernel regression interpolated points to later do confidence intervals on:
          KRbiases[reachtype+1,condition+1,,ppno] <- kernelRegression(x=hand,y=bias,width=10,interpoints=interpoints)
          
        }
        
        confint_kernelregression <- matrix(data=NA,nrow=length(interpoints),ncol=length(percentiles))
        
        for (intpt in c(1:length(interpoints))) {
          
          confint_kernelregression[intpt,2] <- mean(KRbiases[reachtype+1,condition+1,intpt,])
          confint_kernelregression[intpt,c(1,3)] <- t.interval(KRbiases[reachtype+1,condition+1,intpt,])
          
        }
        
        # reachtype (0,1)
        # condition (0,1)
        
        # put the current combination in data frame:
        if (reachtype == 0 & condition == 0) {
          localizationCI <- data.frame(interpoints,confint_kernelregression[,1],confint_kernelregression[,2],confint_kernelregression[,3])
          colnames(localizationCI) <- c('angle',sprintf(c('%s_p2.5','%s_p50','%s_p97.5'),sprintf('%s%s',c('AL','RO')[condition+1],c('act','pas')[reachtype+1])))
        } else {
          tempLocalizationCI <- data.frame(confint_kernelregression[,1],confint_kernelregression[,2],confint_kernelregression[,3])
          colnames(tempLocalizationCI) <- sprintf(c('%s_p2.5','%s_p50','%s_p97.5'),sprintf('%s%s',c('AL','RO')[condition+1],c('act','pas')[reachtype+1]))
          localizationCI <- cbind(localizationCI, tempLocalizationCI)
        }
        
      }
      
      # now within reachtype, get a confidence interval on the difference between rotated and aligned:
      
      subdata <- KRbiases[reachtype+1,2,,] - KRbiases[reachtype+1,1,,]
      
      confint_kernelregression <- matrix(data=NA,nrow=length(interpoints),ncol=length(percentiles))
      
      for (intpt in c(1:length(interpoints))) {
        
        confint_kernelregression[intpt,2] <- mean(subdata[intpt,])
        confint_kernelregression[intpt,c(1,3)] <- t.interval(subdata[intpt,])
        
      }
      
      tempLocalizationCI <- data.frame(confint_kernelregression[,1],confint_kernelregression[,2],confint_kernelregression[,3])
      colnames(tempLocalizationCI) <- sprintf(c('%s_p2.5','%s_p50','%s_p97.5'),sprintf('%s',c('act','pas')[reachtype+1]))
      localizationCI <- cbind(localizationCI, tempLocalizationCI)
      
    }
    
    
    # now within reachtype, get a confidence interval on the difference between rotated and aligned:
    
    subdata <- (KRbiases[1,2,,] - KRbiases[1,1,,]) - (KRbiases[2,2,,] - KRbiases[2,1,,])
    
    confint_kernelregression <- matrix(data=NA,nrow=length(interpoints),ncol=length(percentiles))
    
    for (intpt in c(1:length(interpoints))) {
      
      confint_kernelregression[intpt,2] <- mean(subdata[intpt,])
      confint_kernelregression[intpt,c(1,3)] <- t.interval(subdata[intpt,])
      
    }
    
    tempLocalizationCI <- data.frame(confint_kernelregression[,1],confint_kernelregression[,2],confint_kernelregression[,3])
    colnames(tempLocalizationCI) <- sprintf(c('%s_p2.5','%s_p50','%s_p97.5'),'PredCons')
    localizationCI <- cbind(localizationCI, tempLocalizationCI)
    
    write.csv(localizationCI, sprintf('data/%s_localization_tCI.csv',group), row.names=FALSE)
    
  }
  
}



# Figures -----

oldPlotLocalizationShifts <- function(target='inline') {
  
  styles <- getStyle()
  
  fw <- styles$figwidth[1]
  fh <- fw * (8/7.5)
  fs <- styles$fontsize[1]
  fs <- 0.85
  
  if (target == 'svg') {
    svglite::svglite(file='doc/Fig5.svg', width=fw, height=fh, system_fonts=list(sans='Arial'))
  }
  if (target == 'tiff') {
    tiff(filename='doc/Fig5.tiff',width=fw*1200,height=fh*1200,units='px',type='cairo',compression='lzw',res=1200)
  }
  if (target == 'pdf') {
    pdf(file='doc/Fig5.pdf',width=fw,height=fh)
  }
  
  
  #par(mar=c(4,4,2,0.1))
  
  # 1   2
  # 3 4 5
  # 6 7 8
  # layout(matrix(c(1,1,1,2,2,2,3,3,4,4,5,5,6,6,7,7,8,8), nrow=3, ncol=6, byrow = TRUE))
  
  # 1   2
  # 3   4
  # 5 6 7
  layout(matrix(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,6,6,7,7), nrow=3, ncol=6, byrow = TRUE))
  
  par(mar=c(3,3,2,0.1), cex=1, cex.main=1, cex.axis=0.85, cex.lab=0.85, mgp=c(3,.75,0))
  
  
  # # # # # # # # # #
  # panels A & B: aligned active and passive localization biases
  
  for (reachtype.idx in 1:2) {
    
    reachtype <- c('active','passive')[reachtype.idx]
    
    passive <- 2 - reachtype.idx
    
    #plot(x=c(30,150),y=c(0,0),type='l',lty=2,col=rgb(127, 127, 127, 255, max = 255),xlim=c(30,150),ylim=c(-17,17),main=sprintf('%s localization bias',reachtype),xlab='hand angle [°]',ylab='localization bias [°]',axes=FALSE,font.main=1)
    plot(x=c(30,150),y=c(0,0),type='l',lty=2,col=rgb(127, 127, 127, 255, max = 255),xlim=c(30,150),ylim=c(-17,17),main=sprintf('%s localization bias',reachtype),xlab='',ylab='',axes=FALSE,font.main=1)
    
    title(xlab='hand angle [°]', line=2)
    title(ylab='localization bias [°]', line=2)
  
    mtext(c('A','B')[reachtype.idx], outer=FALSE, side=3, las=1, line=1, adj=0, padj=1)
    axis(1, at=c(45,90,135),cex.axis=fs)
    axis(2, at=c(-12,0,12),cex.axis=fs)
    
    for (groupno in c(1:length(styles$group))) {
      
      group <- styles$group[groupno]
      
      localization <- read.csv(sprintf('data/%s_localizationbias.csv',group))
      localization <- localization[which(localization$passive == passive & localization$rotated == 0),]
      
      participants <- unique(localization$participant)
      
      for (participant in participants) {
        
        partbias <- localization[which(localization$participant == participant),]
        lines(partbias$handlocation_deg, partbias$localizationbias_deg, col=as.character(styles$color_trans[groupno]))
        
      }
      
      groupbias <- aggregate(localizationbias_deg ~ handlocation_deg, data=localization, FUN=mean, na.rm=F)
      idx <- which(groupbias$handlocation_deg >= 50 & groupbias$handlocation_deg <= 130)
      
      pX <- c()
      pY <- c()
      
      for (handangle in seq(130,50,-1)) {
        CI <- t.interval(localization$localizationbias_deg[which(localization$handlocation_deg == handangle)])
        pX <- c(handangle, pX, handangle)
        pY <- c(CI[1], pY, CI[2])
      }
      
      polygon(pX, pY, col=as.character(styles$color_trans[groupno]), border=NA)
      
      lines(groupbias$handlocation_deg[idx], groupbias$localizationbias_deg[idx], col=as.character(styles$color_solid[groupno]), lw=2, lty=1)
      
    }
  
  }
  
  # # # # # # # # # #
  # panels C & D: active and passive localization
  
  for (reachtype.idx in 1:2) {
    
    reachtype <- c('active','passive')[reachtype.idx]
    
    # create plot panel with the right properties
    #plot(c(25,175),c(0,0),type='l',main=sprintf('%s shift',reachtype),xlim=c(25,175),ylim=c(2,-17),axes=FALSE,xlab='hand angle [°]', ylab='localization shift [°]',lty=2,col=rgb(.5,.5,.5),font.main=1)
    plot(c(25,175),c(0,0),type='l',main=sprintf('%s shift',reachtype),xlim=c(25,175),ylim=c(2,-17),axes=FALSE,xlab='', ylab='',lty=2,col=rgb(.5,.5,.5),font.main=1)
    
    title(xlab='hand angle [°]', line=2)
    title(ylab='localization shift [°]', line=2)
    
    #mtext(c('A','B')[reachtype.idx], side=3, outer=TRUE, at=c(c(0,1/3)[reachtype.idx],1), line=-1, adj=0, padj=1)
    mtext(c('C','D')[reachtype.idx], outer=FALSE, side=3, las=1, line=1, adj=0, padj=1)
    
    axis(1, at=c(45,90,135),cex.axis=fs)
    axis(2, at=c(0,-6,-12),cex.axis=fs)
    
    
    for (groupno in c(1:length(styles$group))) {
      
      group <- styles$group[groupno]
      
      localization <- read.csv(sprintf('data/%s_localization_tCI.csv',group))
      
      angles <- localization$angle
      lo <- localization[,sprintf('%s_p2.5',c('act','pas')[reachtype.idx])]
      hi <- localization[,sprintf('%s_p97.5',c('act','pas')[reachtype.idx])]
      
      idx <- which((angles >= 30) & (angles <= 150))
      
      coord.x = c(angles[idx],rev(angles[idx]));
      coord.y = c(lo[idx],rev(hi[idx]))
      polygon(coord.x,coord.y,col=as.character(styles$color_trans[groupno]),border=NA)
      
    }
    
    for (groupno in c(1:length(styles$group))) {
      
      group <- styles$group[groupno]
      
      localization <- read.csv(sprintf('data/%s_localization_tCI.csv',group))
      
      angles <- localization$angle
      centre <- localization[,sprintf('%s_p50',c('act','pas')[reachtype.idx])]
      
      idx <- which((angles >= 30) & (angles <= 150))
      
      lines(angles[idx],centre[idx],col=as.character(styles$color_solid[groupno]),lw=2,lty=styles$linestyle[groupno])
      
    }
    
    # ADD AVERAGE DOTS
    
    for (groupno in c(1:length(styles$group))) {
      
      group <- styles$group[groupno]
      
      localization <- read.csv(sprintf('data/%s_loc_AOV.csv',group))
      localization <- localization[which(localization$passive_b == (reachtype.idx-1)),]
      localization <- aggregate(bias_deg ~ participant*rotated_b, data=localization, FUN=mean)
      shift <- localization$bias_deg[which(localization$rotated_b == 1)] - localization$bias_deg[which(localization$rotated_b == 0)]
      
      xloc <- 150 + (groupno*6)
      CI <- t.interval(shift)
      #arrows(xloc, CI[2], xloc, CI[1], length=0.05, angle=90, code=3, col=as.character(styles$color_solid[groupno]), lty=styles$linestyle[groupno])
      lines(c(xloc, xloc), c(CI[2], CI[1]), col=as.character(styles$color_solid[groupno]), lty=styles$linestyle[groupno])
      #lines(c(xloc-1.5, xloc+1.5), c(CI[1], CI[1]), col=as.character(styles$color_solid[groupno]), lty=1)
      #lines(c(xloc-1.5, xloc+1.5), c(CI[2], CI[2]), col=as.character(styles$color_solid[groupno]), lty=1)
      points(xloc, mean(shift), col=as.character(styles$color_solid[groupno]), pch=16, cex=styles$pointsize[2])
      
    }
    
  }
  
  legend(90,-17.5,as.character(styles$label),col=as.character(styles$color_solid),lty=styles$linestyle,bty='n',lw=2,cex=fs, seg.len=3)
  
  # # # # # # # # # #
  # panel E: predicted sensory consequences
  
  # plot(c(25,175),c(0,0),type='l',main='predicted consequences',xlim=c(25,175),ylim=c(2,-17),axes=FALSE,xlab='hand angle [°]', ylab='update [°]',lty=2,col=rgb(.5,.5,.5),font.main=1)
  # 
  # mtext('E', outer=FALSE, side=3, las=1, line=1, adj=0, padj=1)
  # 
  # axis(1, at=c(45,90,135),cex.axis=0.85)
  # axis(2, at=c(0,-5,-10,-15),cex.axis=0.85)
  # 
  # for (groupno in c(1:length(styles$group))) {
  #   
  #   group <- styles$group[groupno]
  #   
  #   localization <- read.csv(sprintf('data/%s_localization_tCI.csv',group))
  #   
  #   angles <- localization$angle
  #   lo <- localization$PredCons_p2.5
  #   hi <- localization$PredCons_p97.5
  #   
  #   idx <- which((angles >= 30) & (angles <= 150))
  #   
  #   coord.x = c(angles[idx],rev(angles[idx]));
  #   coord.y = c(lo[idx],rev(hi[idx]))
  #   polygon(coord.x,coord.y,col=as.character(styles$color_trans[groupno]),border=NA)
  #   
  # }
  # 
  # for (groupno in c(1:length(styles$group))) {
  #   
  #   group <- styles$group[groupno]
  #   
  #   localization <- read.csv(sprintf('data/%s_localization_tCI.csv',group))
  #   
  #   angles <- localization$angle
  #   PSQ <- localization$PredCons_p50
  #   
  #   idx <- which(is.finite(PSQ) & (angles >= 30) & (angles <= 150))
  #   
  #   lines(angles[idx],PSQ[idx],col=as.character(styles$color_solid[groupno]),lw=2,lty=styles$linestyle[groupno])
  #   
  # }
  # 
  # for (groupno in c(1:length(styles$group))) {
  #   
  #   group <- styles$group[groupno]
  #   
  #   shifts <- list()
  #   
  #   for (reachtype.idx in c(1,2)) {
  #     localization <- read.csv(sprintf('data/%s_loc_AOV.csv',group))
  #     localization <- localization[which(localization$passive_b == (reachtype.idx-1)),]
  #     localization <- aggregate(bias_deg ~ participant*rotated_b, data=localization, FUN=mean)
  #     shift <- localization$bias_deg[which(localization$rotated_b == 1)] - localization$bias_deg[which(localization$rotated_b == 0)]
  #     shifts[[reachtype.idx]] <- shift
  #   }
  #   
  #   shift <- shifts[[1]] - shifts[[2]]
  #   
  #   xloc <- 155 + (groupno*4)
  #   CI <- t.interval(shift)
  #   lines(c(xloc, xloc), c(CI[2], CI[1]), col=as.character(styles$color_solid[groupno]), lty=styles$linestyle[groupno])
  #   lines(c(xloc-1.5, xloc+1.5), c(CI[1], CI[1]), col=as.character(styles$color_solid[groupno]), lty=1)
  #   lines(c(xloc-1.5, xloc+1.5), c(CI[2], CI[2]), col=as.character(styles$color_solid[groupno]), lty=1)
  #   points(xloc, mean(shift), col=as.character(styles$color_solid[groupno]), pch=19)
  #   
  # }
  # 
  # 
  # legend(60,-15,as.character(styles$label),col=as.character(styles$color_solid),lty=styles$linestyle,bty='n',lw=2,cex=0.85, seg.len=3)
  
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  #
  #             NOW THE LOCALIZATION PRECISION
  #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  
  # get the data:
  
  SlocVar <- read.csv('data/sEDS_localization_var.csv', stringsAsFactors = FALSE)
  SlocVar$group <- 'sEDS'
  ZlocVar <- read.csv('data/zEDS_localization_var.csv', stringsAsFactors = FALSE)
  ZlocVar$group <- 'zEDS'
  
  locVar <- rbind(SlocVar, ZlocVar)
  locVar$std <- sqrt(locVar$var)
  
  # get the y-axis equal for all panels:
  
  ylims <- c(0,12) # one individual point now falls of the figure
  
  # # # # # # # # # #
  # panel F/E: localization variance as in ANOVA
  
  #plot(c(0.5,2.5),c(0,0),col=rgb(0.5,0.5,0.5),type='l',lty=2,xlim=c(0.5,2.5),ylim=ylims,main='localization\nprecision',xlab='condition',ylab='localization precision, SD [°]',xaxt='n',yaxt='n',bty='n',font.main=1)
  plot(c(0.5,2.5),c(0,0),col=rgb(0.5,0.5,0.5),type='l',lty=2,xlim=c(0.5,2.5),ylim=ylims,main='localization\nprecision',xlab='',ylab='',xaxt='n',yaxt='n',bty='n',font.main=1)
  
  title(xlab='condition', line=2)
  title(ylab='localization precision, SD [°]', line=2)
  
  #mtext('D', side=3, outer=TRUE, at=c(0,1), line=-1, adj=0, padj=1)
  mtext('E', outer=FALSE, side=3, las=1, line=1, adj=0, padj=1)
  avgLocCI <- aggregate(std ~ group + passive + rotated, data=locVar, FUN=t.interval)
  
  for (group in styles$group) {
    
    col <- as.character(styles$color_trans[which(styles$group == group)])
    
    for (passive in c(0,1)) {
      
      loLV <- avgLocCI$std[which(avgLocCI$group == group & avgLocCI$passive == passive),1]
      hiLV <- avgLocCI$std[which(avgLocCI$group == group & avgLocCI$passive == passive),2]
      
      Y <- c(loLV, rev(hiLV))
      X <- c(1,2,2,1)
      
      polygon(X,Y,col=col,border=NA)
      
    }
    
  }
  
  avgLocStd <- aggregate(std ~ group + passive + rotated, data=locVar, FUN=mean, na.rm=TRUE)
  
  for (group in styles$group) {
    
    col <- as.character(styles$color_solid[which(styles$group == group)])
    
    for (passive in c(0,1)) {
      
      linestyle <- passive + 1
      
      LV <- avgLocStd$std[which(avgLocStd$group == group & avgLocStd$passive == passive)]
      
      lines(x=c(1,2),y=LV,lty=linestyle,col=col)
      
    }
    
  }
  
  axis(side=1,at=c(1,2),labels=c('aligned','rotated'))
  axis(side=2,at=seq(0,12,4))
  
  
  # # # # # # # # # #
  # panel G/F: aligned localization variance descriptives
  
  
  #plot(c(0,5),c(0,0),col=rgb(0.5,0.5,0.5),type='l',lty=2,xlim=c(0.5,4.5),ylim=ylims,main='aligned',xlab='condition',ylab='individual precision, SD [°]',xaxt='n',yaxt='n',bty='n',font.main=1)
  plot(c(0,5),c(0,0),col=rgb(0.5,0.5,0.5),type='l',lty=2,xlim=c(0.5,4.5),ylim=ylims,main='aligned',xlab='',ylab='',xaxt='n',yaxt='n',bty='n',font.main=1)
  
  title(xlab='condition', line=2)
  title(ylab='precision, SD [°]', line=1)
  
  mtext('F', outer=FALSE, side=3, las=1, line=1, adj=0, padj=1)
  
  ####
  
  for (passive in c(0,1)) {
    
    for (groupno in c(1,2)) {
      
      group <- styles$group[groupno]
      
      #cols <- styles$color_solid[groupno]
      #colt <- styles$color_trans[groupno]
      
      locvars <- locVar$std[which(locVar$group == group & locVar$rotated == 0 & locVar$passive == passive)]
      
      #print(str(nocursors))
      X <- rep( (passive*2)+(groupno)-.33, length(locvars) )
      #X <- rep(((passive*2)+((groupno-1.7)/2.5))+2,length(locvars))
      
      Y <- locvars
      
      points(x=X,y=Y,pch=16,cex=styles$pointsize[2],col=as.character(styles$color_trans[groupno]))
      
      meandist <- getConfidenceInterval(data=locvars, method='bootstrap', resamples=5000, FUN=mean, returndist=TRUE)
      
      DX <- meandist$density$x
      DY <- meandist$density$y / max(meandist$density$y) / 3
      
      DX <- c(DX[1], DX, DX[length(DX)])
      DY <- c(0,     DY, 0)
      
      Xoffset <- (passive*2)+(groupno)
      
      polygon(x=DY+Xoffset, y=DX, border=FALSE, col=as.character(styles$color_trans[groupno]))
      
      lines(x=rep(Xoffset,2),y=meandist$CI95,col=as.character(styles$color_solid[groupno]))
      
      points(x=Xoffset,y=mean(locvars),pch=16,cex=styles$pointsize[2],col=as.character(styles$color_solid[groupno]))
      
    }
    
  }
  
  ####
  
  axis(side=1, at=c(1.5,3.5), labels=c('active','passive'))
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  #
  # panel F/G: rotated localization variance descriptives
  #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  
  #plot(c(0,5),c(0,0),col=rgb(0.5,0.5,0.5),type='l',lty=2,xlim=c(0.5,4.5),ylim=ylims,main='rotated',xlab='condition',ylab='individual precision, SD [°]',xaxt='n',yaxt='n',bty='n',font.main=1)
  plot(c(0,5),c(0,0),col=rgb(0.5,0.5,0.5),type='l',lty=2,xlim=c(0.5,4.5),ylim=ylims,main='rotated',xlab='',ylab='',xaxt='n',yaxt='n',bty='n',font.main=1)
  
  title(xlab='condition', line=2)
  title(ylab='precision, SD [°]', line=1)
  
  mtext('G', outer=FALSE, side=3, las=1, line=1, adj=0, padj=1)
  
  for (passive in c(0,1)) {
    
    for (groupno in c(1,2)) {
      
      group <- styles$group[groupno]
      
      #cols <- styles$color_solid[groupno]
      #colt <- styles$color_trans[groupno]
      
      locvars <- locVar$std[which(locVar$group == group & locVar$rotated == 1 & locVar$passive == passive)]
      
      #print(str(nocursors))
      X <- rep( (passive*2)+(groupno)-.33, length(locvars) )
      #X <- rep(((passive*2)+((groupno-1.7)/2.5))+2,length(locvars))
      
      Y <- locvars
      
      points(x=X,y=Y,pch=16,cex=styles$pointsize[2],col=as.character(styles$color_trans[groupno]))
      
      meandist <- getConfidenceInterval(data=locvars, method='bootstrap', resamples=5000, FUN=mean, returndist=TRUE)
      
      DX <- meandist$density$x
      DY <- meandist$density$y / max(meandist$density$y) / 3
      
      DX <- c(DX[1], DX, DX[length(DX)])
      DY <- c(0,     DY, 0)
      
      Xoffset <- (passive*2)+(groupno)
      
      polygon(x=DY+Xoffset, y=DX, border=FALSE, col=as.character(styles$color_trans[groupno]))
      
      lines(x=rep(Xoffset,2),y=meandist$CI95,col=as.character(styles$color_solid[groupno]))
      
      points(x=Xoffset,y=mean(locvars),pch=16,cex=styles$pointsize[2],col=as.character(styles$color_solid[groupno]))
      
    }
    
  }
  
  axis(side=1, at=c(1.5,3.5), labels=c('active','passive'))
  
  if (target %in% c('tiff','svg','pdf')) {
    dev.off()
  }
  
  
}



plotLocalizationShifts <- function(target='inline') {
  
  styles <- getStyle()
  
  fw <- styles$figwidth[1]
  fh <- fw * (5.7/8)
  fs <- styles$fontsize[1]
  fs <- 0.85
  
  
  # 40000000 pixles in total only...
  
  if (target == 'svg') {
    svglite::svglite(file='doc/Fig5.svg', width=fw, height=fh, scaling=1, fix_text_size = FALSE)
  }
  if (target == 'tiff') {
    tiff(filename='doc/Fig5.tiff',width=fw*1150,height=fh*1150,units='px',type='cairo',compression='lzw',res=1200)
  }
  if (target == 'pdf') {
    pdf(file='doc/Fig5.pdf',width=fw,height=fh)
  }
  
  
  layout(matrix(c(1,2,3,4), nrow=2, ncol=2, byrow = TRUE))
  
  par(mar=c(3,3,2,0.1), cex=1, cex.main=1, cex.axis=0.85, cex.lab=0.85, mgp=c(3,.75,0))
  
  
  # # # # # # # # # #
  # panels A & B: aligned active and passive localization biases
  
  for (reachtype.idx in 1:2) {
    
    reachtype <- c('active','passive')[reachtype.idx]
    
    passive <- 2 - reachtype.idx
    
    #plot(x=c(30,150),y=c(0,0),type='l',lty=2,col=rgb(127, 127, 127, 255, max = 255),xlim=c(30,150),ylim=c(-17,17),main=sprintf('%s localization bias',reachtype),xlab='hand angle [°]',ylab='localization bias [°]',axes=FALSE,font.main=1)
    plot(x=c(30,150),y=c(0,0),type='l',lty=2,col=rgb(127, 127, 127, 255, max = 255),xlim=c(30,150),ylim=c(-17,17),main=sprintf('%s localization bias',reachtype),xlab='',ylab='',axes=FALSE,font.main=1)
    
    title(xlab='hand angle [°]', line=2)
    title(ylab='localization bias [°]', line=2)
    
    mtext(c('A','B')[reachtype.idx], outer=FALSE, side=3, las=1, line=1, adj=0, padj=1)
    axis(1, at=c(45,90,135),cex.axis=fs)
    axis(2, at=c(-12,0,12),cex.axis=fs)
    
    for (groupno in c(1:length(styles$group))) {
      
      group <- styles$group[groupno]
      
      localization <- read.csv(sprintf('data/%s_localizationbias.csv',group))
      localization <- localization[which(localization$passive == passive & localization$rotated == 0),]
      
      participants <- unique(localization$participant)
      
      for (participant in participants) {
        
        partbias <- localization[which(localization$participant == participant),]
        lines(partbias$handlocation_deg, partbias$localizationbias_deg, col=as.character(styles$color_trans[groupno]))
        
      }
      
      groupbias <- aggregate(localizationbias_deg ~ handlocation_deg, data=localization, FUN=mean, na.rm=F)
      idx <- which(groupbias$handlocation_deg >= 50 & groupbias$handlocation_deg <= 130)
      
      pX <- c()
      pY <- c()
      
      for (handangle in seq(130,50,-1)) {
        CI <- t.interval(localization$localizationbias_deg[which(localization$handlocation_deg == handangle)])
        pX <- c(handangle, pX, handangle)
        pY <- c(CI[1], pY, CI[2])
      }
      
      polygon(pX, pY, col=as.character(styles$color_trans[groupno]), border=NA)
      
      lines(groupbias$handlocation_deg[idx], groupbias$localizationbias_deg[idx], col=as.character(styles$color_solid[groupno]), lw=2, lty=1)
      
    }
    
  }
  
  # # # # # # # # # #
  # panels C & D: active and passive localization
  
  for (reachtype.idx in 1:2) {
    
    reachtype <- c('active','passive')[reachtype.idx]
    
    # create plot panel with the right properties
    #plot(c(25,175),c(0,0),type='l',main=sprintf('%s shift',reachtype),xlim=c(25,175),ylim=c(2,-17),axes=FALSE,xlab='hand angle [°]', ylab='localization shift [°]',lty=2,col=rgb(.5,.5,.5),font.main=1)
    plot(c(25,175),c(0,0),type='l',main=sprintf('%s shift',reachtype),xlim=c(25,175),ylim=c(2,-17),axes=FALSE,xlab='', ylab='',lty=2,col=rgb(.5,.5,.5),font.main=1)
    
    title(xlab='hand angle [°]', line=2)
    title(ylab='localization shift [°]', line=2)
    
    #mtext(c('A','B')[reachtype.idx], side=3, outer=TRUE, at=c(c(0,1/3)[reachtype.idx],1), line=-1, adj=0, padj=1)
    mtext(c('C','D')[reachtype.idx], outer=FALSE, side=3, las=1, line=1, adj=0, padj=1)
    
    axis(1, at=c(45,90,135),cex.axis=fs)
    axis(2, at=c(0,-6,-12),cex.axis=fs)
    
    
    for (groupno in c(1:length(styles$group))) {
      
      group <- styles$group[groupno]
      
      localization <- read.csv(sprintf('data/%s_localization_tCI.csv',group))
      
      angles <- localization$angle
      lo <- localization[,sprintf('%s_p2.5',c('act','pas')[reachtype.idx])]
      hi <- localization[,sprintf('%s_p97.5',c('act','pas')[reachtype.idx])]
      
      idx <- which((angles >= 30) & (angles <= 150))
      
      coord.x = c(angles[idx],rev(angles[idx]));
      coord.y = c(lo[idx],rev(hi[idx]))
      polygon(coord.x,coord.y,col=as.character(styles$color_trans[groupno]),border=NA)
      
    }
    
    for (groupno in c(1:length(styles$group))) {
      
      group <- styles$group[groupno]
      
      localization <- read.csv(sprintf('data/%s_localization_tCI.csv',group))
      
      angles <- localization$angle
      centre <- localization[,sprintf('%s_p50',c('act','pas')[reachtype.idx])]
      
      idx <- which((angles >= 30) & (angles <= 150))
      
      lines(angles[idx],centre[idx],col=as.character(styles$color_solid[groupno]),lw=2,lty=styles$linestyle[groupno])
      
    }
    
    # ADD AVERAGE DOTS
    
    for (groupno in c(1:length(styles$group))) {
      
      group <- styles$group[groupno]
      
      localization <- read.csv(sprintf('data/%s_loc_AOV.csv',group))
      localization <- localization[which(localization$passive_b == (reachtype.idx-1)),]
      localization <- aggregate(bias_deg ~ participant*rotated_b, data=localization, FUN=mean)
      shift <- localization$bias_deg[which(localization$rotated_b == 1)] - localization$bias_deg[which(localization$rotated_b == 0)]
      
      # xlim: 175
      xloc <- 145 + (groupno*10)
      xlocs <- xloc + c(0,7.5)
      CI <- t.interval(shift)
      
      # first version:
      #arrows(xloc, CI[2], xloc, CI[1], length=0.05, angle=90, code=3, col=as.character(styles$color_solid[groupno]), lty=styles$linestyle[groupno])
      #lines(c(xloc-1.5, xloc+1.5), c(CI[1], CI[1]), col=as.character(styles$color_solid[groupno]), lty=1)
      #lines(c(xloc-1.5, xloc+1.5), c(CI[2], CI[2]), col=as.character(styles$color_solid[groupno]), lty=1)
      
      # second version:
      #lines(c(xloc, xloc), c(CI[2], CI[1]), col=as.character(styles$color_solid[groupno]), lty=styles$linestyle[groupno])
      #points(xloc, mean(shift), col=as.character(styles$color_solid[groupno]), pch=16, cex=styles$pointsize[2])
      
      # "error bar" version:
      polygon(x=c(xlocs,rev(xlocs)),
              y=rep(CI,each=2),
              col=as.character(styles$color_trans[groupno]),border=NA)
      lines(x=xlocs,y=rep(mean(shift),2),
            col=as.character(styles$color_solid[groupno]))
      
      
      
    }
    
  }
  
  legend(90,-17.5,as.character(styles$label),col=as.character(styles$color_solid),lty=styles$linestyle,bty='n',lw=2,cex=fs, seg.len=3)
  
  if (target %in% c('tiff','svg','pdf')) {
    dev.off()
  }
  
  
}

# OLDplotLocalizationPrecision <- function(target='inline') {
#   
#   styles <- getStyle()
#   
#   fw <- styles$figwidth[1]
#   fh <- fw * (2.85/8)
#   fs <- styles$fontsize[1]
#   fs <- 0.85
#   
#   if (target == 'svg') {
#     svglite::svglite(file='doc/Fig6.svg', width=fw, height=fh, scaling=1, fix_text_size = FALSE)
#   }
#   if (target == 'tiff') {
#     tiff(filename='doc/Fig6.tiff',width=fw*1200,height=fh*1200,units='px',type='cairo',compression='lzw',res=1200)
#   }
#   if (target == 'pdf') {
#     pdf(file='doc/Fig6.pdf',width=fw,height=fh)
#   }
#   
#   
#   layout(matrix(c(1,2,3), nrow=1, ncol=3, byrow = TRUE))
#   
#   par(mar=c(3,3,2,0.1), cex=1, cex.main=1, cex.axis=0.85, cex.lab=0.85, mgp=c(3,.75,0))
#   
# 
#   # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#   #
#   #             NOW THE LOCALIZATION PRECISION
#   #
#   # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#   
#   # get the data:
#   
#   SlocVar <- read.csv('data/sEDS_localization_var.csv', stringsAsFactors = FALSE)
#   SlocVar$group <- 'sEDS'
#   ZlocVar <- read.csv('data/zEDS_localization_var.csv', stringsAsFactors = FALSE)
#   ZlocVar$group <- 'zEDS'
#   
#   locVar <- rbind(SlocVar, ZlocVar)
#   locVar$std <- sqrt(locVar$var)
#   
#   # get the y-axis equal for all panels:
#   
#   ylims <- c(0,12) # one individual point now falls of the figure
#   
#   # # # # # # # # # #
#   # panel F/E: localization variance as in ANOVA
#   
#   #plot(c(0.5,2.5),c(0,0),col=rgb(0.5,0.5,0.5),type='l',lty=2,xlim=c(0.5,2.5),ylim=ylims,main='localization\nprecision',xlab='condition',ylab='localization precision, SD [°]',xaxt='n',yaxt='n',bty='n',font.main=1)
#   plot(c(0.5,2.5),c(0,0),col=rgb(0.5,0.5,0.5),type='l',lty=2,xlim=c(0.5,2.5),ylim=ylims,main='localization\nprecision',xlab='',ylab='',xaxt='n',yaxt='n',bty='n',font.main=1)
#   
#   title(xlab='condition', line=2)
#   title(ylab='localization precision, SD [°]', line=2)
#   
#   #mtext('D', side=3, outer=TRUE, at=c(0,1), line=-1, adj=0, padj=1)
#   mtext('A', outer=FALSE, side=3, las=1, line=1, adj=0, padj=1)
#   avgLocCI <- aggregate(std ~ group + passive + rotated, data=locVar, FUN=t.interval)
#   
#   for (group in styles$group) {
#     
#     col <- as.character(styles$color_trans[which(styles$group == group)])
#     
#     for (passive in c(0,1)) {
#       
#       loLV <- avgLocCI$std[which(avgLocCI$group == group & avgLocCI$passive == passive),1]
#       hiLV <- avgLocCI$std[which(avgLocCI$group == group & avgLocCI$passive == passive),2]
#       
#       Y <- c(loLV, rev(hiLV))
#       X <- c(1,2,2,1)
#       
#       polygon(X,Y,col=col,border=NA)
#       
#     }
#     
#   }
#   
#   avgLocStd <- aggregate(std ~ group + passive + rotated, data=locVar, FUN=mean, na.rm=TRUE)
#   
#   for (group in styles$group) {
#     
#     col <- as.character(styles$color_solid[which(styles$group == group)])
#     
#     for (passive in c(0,1)) {
#       
#       linestyle <- passive + 1
#       
#       LV <- avgLocStd$std[which(avgLocStd$group == group & avgLocStd$passive == passive)]
#       
#       lines(x=c(1,2),y=LV,lty=linestyle,col=col)
#       
#     }
#     
#   }
#   
#   legend(0.7,4.5,as.character(styles$label),col=as.character(styles$color_solid),lty=styles$linestyle,bty='n',lw=2,cex=fs, seg.len=3)
#   
#   axis(side=1,at=c(1,2),labels=c('aligned','rotated'))
#   axis(side=2,at=seq(0,12,4))
#   
#   
#   # # # # # # # # # #
#   # panel G/F: aligned localization variance descriptives
#   
#   
#   #plot(c(0,5),c(0,0),col=rgb(0.5,0.5,0.5),type='l',lty=2,xlim=c(0.5,4.5),ylim=ylims,main='aligned',xlab='condition',ylab='individual precision, SD [°]',xaxt='n',yaxt='n',bty='n',font.main=1)
#   plot(c(0,5),c(0,0),col=rgb(0.5,0.5,0.5),type='l',lty=2,xlim=c(0.5,4.5),ylim=ylims,main='aligned',xlab='',ylab='',xaxt='n',yaxt='n',bty='n',font.main=1)
#   
#   title(xlab='condition', line=2)
#   title(ylab='precision, SD [°]', line=1)
#   
#   mtext('B', outer=FALSE, side=3, las=1, line=1, adj=0, padj=1)
#   
#   # # # #
#   
#   for (passive in c(0,1)) {
#     
#     for (groupno in c(1,2)) {
#       
#       group <- styles$group[groupno]
#       
#       #cols <- styles$color_solid[groupno]
#       #colt <- styles$color_trans[groupno]
#       
#       locvars <- locVar$std[which(locVar$group == group & locVar$rotated == 0 & locVar$passive == passive)]
#       
#       #print(str(nocursors))
#       X <- rep( (passive*2)+(groupno)-.33, length(locvars) )
#       #X <- rep(((passive*2)+((groupno-1.7)/2.5))+2,length(locvars))
#       
#       Y <- locvars
#       
#       points(x=X,y=Y,pch=16,cex=styles$pointsize[2],col=as.character(styles$color_trans[groupno]))
#       
#       meandist <- getConfidenceInterval(data=locvars, method='bootstrap', resamples=5000, FUN=mean, returndist=TRUE)
#       
#       DX <- meandist$density$x
#       DY <- meandist$density$y / max(meandist$density$y) / 3
#       
#       DX <- c(DX[1], DX, DX[length(DX)])
#       DY <- c(0,     DY, 0)
#       
#       Xoffset <- (passive*2)+(groupno)
#       
#       polygon(x=DY+Xoffset, y=DX, border=FALSE, col=as.character(styles$color_trans[groupno]))
#       
#       lines(x=rep(Xoffset,2),y=meandist$CI95,col=as.character(styles$color_solid[groupno]))
#       
#       points(x=Xoffset,y=mean(locvars),pch=16,cex=styles$pointsize[2],col=as.character(styles$color_solid[groupno]))
#       
#     }
#     
#   }
#   
#   # # # #
#   
#   axis(side=1, at=c(1.5,3.5), labels=c('active','passive'))
#   
#   # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#   #
#   # panel F/G: rotated localization variance descriptives
#   #
#   # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#   
#   #plot(c(0,5),c(0,0),col=rgb(0.5,0.5,0.5),type='l',lty=2,xlim=c(0.5,4.5),ylim=ylims,main='rotated',xlab='condition',ylab='individual precision, SD [°]',xaxt='n',yaxt='n',bty='n',font.main=1)
#   plot(c(0,5),c(0,0),col=rgb(0.5,0.5,0.5),type='l',lty=2,xlim=c(0.5,4.5),ylim=ylims,main='rotated',xlab='',ylab='',xaxt='n',yaxt='n',bty='n',font.main=1)
#   
#   title(xlab='condition', line=2)
#   title(ylab='precision, SD [°]', line=1)
#   
#   mtext('C', outer=FALSE, side=3, las=1, line=1, adj=0, padj=1)
#   
#   for (passive in c(0,1)) {
#     
#     for (groupno in c(1,2)) {
#       
#       group <- styles$group[groupno]
#       
#       #cols <- styles$color_solid[groupno]
#       #colt <- styles$color_trans[groupno]
#       
#       locvars <- locVar$std[which(locVar$group == group & locVar$rotated == 1 & locVar$passive == passive)]
#       
#       #print(str(nocursors))
#       X <- rep( (passive*2)+(groupno)-.33, length(locvars) )
#       #X <- rep(((passive*2)+((groupno-1.7)/2.5))+2,length(locvars))
#       
#       Y <- locvars
#       
#       points(x=X,y=Y,pch=16,cex=styles$pointsize[2],col=as.character(styles$color_trans[groupno]))
#       
#       meandist <- getConfidenceInterval(data=locvars, method='bootstrap', resamples=5000, FUN=mean, returndist=TRUE)
#       
#       DX <- meandist$density$x
#       DY <- meandist$density$y / max(meandist$density$y) / 3
#       
#       DX <- c(DX[1], DX, DX[length(DX)])
#       DY <- c(0,     DY, 0)
#       
#       Xoffset <- (passive*2)+(groupno)
#       
#       polygon(x=DY+Xoffset, y=DX, border=FALSE, col=as.character(styles$color_trans[groupno]))
#       
#       lines(x=rep(Xoffset,2),y=meandist$CI95,col=as.character(styles$color_solid[groupno]))
#       
#       points(x=Xoffset,y=mean(locvars),pch=16,cex=styles$pointsize[2],col=as.character(styles$color_solid[groupno]))
#       
#     }
#     
#   }
#   
#   axis(side=1, at=c(1.5,3.5), labels=c('active','passive'))
#   
#   if (target %in% c('tiff','svg','pdf')) {
#     dev.off()
#   }
#   
#   
# }


plotLocalizationPrecision <- function(target='inline') {
  
  styles <- getStyle()
  
  fw <- styles$figwidth[1]
  fh <- fw * (5.7/8)
  fs <- styles$fontsize[1]
  fs <- 0.85
  
  if (target == 'svg') {
    svglite::svglite(file='doc/Fig6.svg', width=fw, height=fh, scaling=1, fix_text_size = FALSE)
  }
  if (target == 'tiff') {
    tiff(filename='doc/Fig6.tiff',width=fw*1200,height=fh*1200,units='px',type='cairo',compression='lzw',res=1200)
  }
  if (target == 'pdf') {
    pdf(file='doc/Fig6.pdf',width=fw,height=fh)
  }
  
  layout(matrix(c(1,2,3,4), nrow=2, ncol=2, byrow = TRUE))
  
  par(mar=c(3,3,2,0.1), cex=1, cex.main=1, cex.axis=0.85, cex.lab=0.85, mgp=c(3,.75,0))
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  #
  #             NOW THE LOCALIZATION PRECISION
  #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  
  # get the data:
  
  SlocVar <- read.csv('data/sEDS_localization_var.csv', stringsAsFactors = FALSE)
  SlocVar$group <- 'sEDS'
  ZlocVar <- read.csv('data/zEDS_localization_var.csv', stringsAsFactors = FALSE)
  ZlocVar$group <- 'zEDS'
  
  locVar <- rbind(SlocVar, ZlocVar)
  locVar$std <- sqrt(locVar$var)
  
  # get the y-axis equal for all panels:
  
  ylims <- c(0,12) # one individual point now falls of the figure
  
  # # # # # # # # # #
  # panel A & B:  aligned localization "variance" as in ANOVA
  
  #plot(c(0.5,2.5),c(0,0),col=rgb(0.5,0.5,0.5),type='l',lty=2,xlim=c(0.5,2.5),ylim=ylims,main='localization\nprecision',xlab='condition',ylab='localization precision, SD [°]',xaxt='n',yaxt='n',bty='n',font.main=1)
  plot(c(0.5,2.5),c(0,0),col=rgb(0.5,0.5,0.5),type='l',lty=2,
       xlim=c(0.5,2.5),ylim=ylims,
       main='aligned\nlocalization precision',xlab='',ylab='',
       xaxt='n',yaxt='n',bty='n',font.main=1)
  
  title(xlab='condition', line=2)
  title(ylab='localization precision, SD [°]', line=2)
  
  #mtext('D', side=3, outer=TRUE, at=c(0,1), line=-1, adj=0, padj=1)
  mtext('A', outer=FALSE, side=3, las=1, line=1, adj=0, padj=1)
  
  avgLocCI <- aggregate(std ~ group + passive + rotated, data=locVar, FUN=t.interval)
  avgLocStd <- aggregate(std ~ group + passive + rotated, data=locVar, FUN=mean, na.rm=TRUE)
  
  for (groupno in c(1:length(styles$group))) {
    
    group <- styles$group[groupno]
    
    col <- as.character(styles$color_trans[which(styles$group == group)])
    
    loLV <- avgLocCI$std[which(avgLocCI$group == group & avgLocCI$rotated == 0),1]
    hiLV <- avgLocCI$std[which(avgLocCI$group == group & avgLocCI$rotated == 0),2]
    
    for (cond_no in c(1,2)) {
      
      Xoffset <- cond_no + (groupno / 2) - 0.875
      Xoffsets <- Xoffset + c(0,.25)
      
      Xpol <- c( Xoffsets, rev(Xoffsets) )
      Ypol <- rep( c(loLV[cond_no], hiLV[cond_no]), each=2 )
      
      polygon(Xpol, Ypol, col=as.character(styles$color_trans[groupno]), border=NA)
      meLV <- avgLocStd$std[which(avgLocCI$group == group & avgLocCI$rotated == 0)]
      lines(Xoffsets, rep(meLV[cond_no], 2), col=as.character(styles$color_solid[groupno]),lw=2)
      
    }
    
  }
  
  # only passive:
  xpos <- 2 + (c(1,2)/2) - 0.75
  lines(xpos,rep(10,2),col='#999999')
  text(mean(xpos),11,'*',cex=2,col='#999999')
  
  axis(side=1,at=c(1,2),labels=c('active','passive'))
  axis(side=2,at=seq(0,12,4))

  # # # # # # # # # #
  # panel B:  rotated localization "variance" as in ANOVA
  
  #plot(c(0.5,2.5),c(0,0),col=rgb(0.5,0.5,0.5),type='l',lty=2,xlim=c(0.5,2.5),ylim=ylims,main='localization\nprecision',xlab='condition',ylab='localization precision, SD [°]',xaxt='n',yaxt='n',bty='n',font.main=1)
  plot(c(0.5,2.5),c(0,0),col=rgb(0.5,0.5,0.5),type='l',lty=2,
       xlim=c(0.5,2.5),ylim=ylims,
       main='rotated\nlocalization precision',xlab='',ylab='',
       xaxt='n',yaxt='n',bty='n',font.main=1)
  
  title(xlab='condition', line=2)
  title(ylab='localization precision, SD [°]', line=2)
  
  #mtext('D', side=3, outer=TRUE, at=c(0,1), line=-1, adj=0, padj=1)
  mtext('B', outer=FALSE, side=3, las=1, line=1, adj=0, padj=1)
  
  #avgLocCI <- aggregate(std ~ group + passive + rotated, data=locVar, FUN=t.interval)
  
  for (groupno in c(1:length(styles$group))) {
    
    group <- styles$group[groupno]
    
    col <- as.character(styles$color_trans[which(styles$group == group)])
    
    loLV <- avgLocCI$std[which(avgLocCI$group == group & avgLocCI$rotated == 1),1]
    hiLV <- avgLocCI$std[which(avgLocCI$group == group & avgLocCI$rotated == 1),2]
    
    for (cond_no in c(1,2)) {
      
      Xoffset <- cond_no + (groupno / 2) - 0.875
      Xoffsets <- Xoffset + c(0,.25)
      
      Xpol <- c( Xoffsets, rev(Xoffsets) )
      Ypol <- rep( c(loLV[cond_no], hiLV[cond_no]), each=2 )
      
      polygon(Xpol, Ypol, col=as.character(styles$color_trans[groupno]), border=NA)
      meLV <- avgLocStd$std[which(avgLocCI$group == group & avgLocCI$rotated == 1)]
      lines(Xoffsets, rep(meLV[cond_no], 2), col=as.character(styles$color_solid[groupno]),lw=2)
      
    }
    
  }
  
  # avgLocStd <- aggregate(std ~ group + passive + rotated, data=locVar, FUN=mean, na.rm=TRUE)
  # 
  # for (group in styles$group) {
  #   
  #   col <- as.character(styles$color_solid[which(styles$group == group)])
  #   
  #   linestyle <- 2
  #   
  #   LV <- avgLocStd$std[which(avgLocStd$group == group & avgLocStd$rotated == 1)]
  #   
  #   lines(x=c(1,2),y=LV,lty=linestyle,col=col)
  #   
  # }
  
  for (condno in c(1,2)) {
    
    xpos <- condno + (c(1,2)/2) - 0.75
    lines(xpos,rep(10,2),col='#999999')
    text(mean(xpos),11,'*',cex=2,col='#999999')
    
  }
  
  axis(side=1,at=c(1,2),labels=c('active','passive'))
  #axis(side=2,at=seq(0,12,4))
  
  legend(0.7,4.5,as.character(styles$label),col=as.character(styles$color_solid),lty=styles$linestyle,bty='n',lw=2,cex=fs, seg.len=3)
  
  # # # # # # # # # #
  # panel C: aligned localization variance descriptives
  
  
  #plot(c(0,5),c(0,0),col=rgb(0.5,0.5,0.5),type='l',lty=2,xlim=c(0.5,4.5),ylim=ylims,main='aligned',xlab='condition',ylab='individual precision, SD [°]',xaxt='n',yaxt='n',bty='n',font.main=1)
  plot(c(0,5),c(0,0),col=rgb(0.5,0.5,0.5),type='l',lty=2,xlim=c(0.5,4.5),ylim=ylims,main='aligned',xlab='',ylab='',xaxt='n',yaxt='n',bty='n',font.main=1)
  
  title(xlab='condition', line=2)
  title(ylab='precision, SD [°]', line=2)
  
  mtext('C', outer=FALSE, side=3, las=1, line=1, adj=0, padj=1)
  
  # # # #
  
  for (passive in c(0,1)) {
    
    for (groupno in c(1,2)) {
      
      group <- styles$group[groupno]
      
      #cols <- styles$color_solid[groupno]
      #colt <- styles$color_trans[groupno]
      
      locvars <- locVar$std[which(locVar$group == group & locVar$rotated == 0 & locVar$passive == passive)]
      
      #print(str(nocursors))
      X <- rep( (passive*2)+(groupno)-.33, length(locvars) )
      #X <- rep(((passive*2)+((groupno-1.7)/2.5))+2,length(locvars))
      
      Y <- locvars
      
      points(x=X,y=Y,pch=16,cex=styles$pointsize[2],col=as.character(styles$color_trans[groupno]))
      
      meandist <- getConfidenceInterval(data=locvars, method='bootstrap', resamples=5000, FUN=mean, returndist=TRUE)
      
      DX <- meandist$density$x
      DY <- meandist$density$y / max(meandist$density$y) / 3
      
      DX <- c(DX[1], DX, DX[length(DX)])
      DY <- c(0,     DY, 0)
      
      Xoffset <- (passive*2)+(groupno)
      
      polygon(x=DY+Xoffset, y=DX, border=FALSE, col=as.character(styles$color_trans[groupno]))
      
      lines(x=rep(Xoffset,2),y=meandist$CI95,col=as.character(styles$color_solid[groupno]))
      
      points(x=Xoffset,y=mean(locvars),pch=16,cex=styles$pointsize[2],col=as.character(styles$color_solid[groupno]))
      
    }
    
  }
  
  ####
  
  axis(side=1, at=c(1.5,3.5), labels=c('active','passive'))
  axis(side=2,at=seq(0,12,4))
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  #
  # panel D: rotated localization variance descriptives
  #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  
  #plot(c(0,5),c(0,0),col=rgb(0.5,0.5,0.5),type='l',lty=2,xlim=c(0.5,4.5),ylim=ylims,main='rotated',xlab='condition',ylab='individual precision, SD [°]',xaxt='n',yaxt='n',bty='n',font.main=1)
  plot(c(0,5),c(0,0),col=rgb(0.5,0.5,0.5),type='l',lty=2,xlim=c(0.5,4.5),ylim=ylims,main='rotated',xlab='',ylab='',xaxt='n',yaxt='n',bty='n',font.main=1)
  
  title(xlab='condition', line=2)
  title(ylab='precision, SD [°]', line=2)
  
  mtext('D', outer=FALSE, side=3, las=1, line=1, adj=0, padj=1)
  
  for (passive in c(0,1)) {
    
    for (groupno in c(1,2)) {
      
      group <- styles$group[groupno]
      
      #cols <- styles$color_solid[groupno]
      #colt <- styles$color_trans[groupno]
      
      locvars <- locVar$std[which(locVar$group == group & locVar$rotated == 1 & locVar$passive == passive)]
      
      #print(str(nocursors))
      X <- rep( (passive*2)+(groupno)-.33, length(locvars) )
      #X <- rep(((passive*2)+((groupno-1.7)/2.5))+2,length(locvars))
      
      Y <- locvars
      
      points(x=X,y=Y,pch=16,cex=styles$pointsize[2],col=as.character(styles$color_trans[groupno]))
      
      meandist <- getConfidenceInterval(data=locvars, method='bootstrap', resamples=5000, FUN=mean, returndist=TRUE)
      
      DX <- meandist$density$x
      DY <- meandist$density$y / max(meandist$density$y) / 3
      
      DX <- c(DX[1], DX, DX[length(DX)])
      DY <- c(0,     DY, 0)
      
      Xoffset <- (passive*2)+(groupno)
      
      polygon(x=DY+Xoffset, y=DX, border=FALSE, col=as.character(styles$color_trans[groupno]))
      
      lines(x=rep(Xoffset,2),y=meandist$CI95,col=as.character(styles$color_solid[groupno]))
      
      points(x=Xoffset,y=mean(locvars),pch=16,cex=styles$pointsize[2],col=as.character(styles$color_solid[groupno]))
      
    }
    
  }
  
  axis(side=1, at=c(1.5,3.5), labels=c('active','passive'))
  #axis(side=2,at=seq(0,12,4))
  
  if (target %in% c('tiff','svg','pdf')) {
    dev.off()
  }
  
  
}


# Statistics -----

localizationAOV <- function(groups=c('sEDS','zEDS')) {
  
  df <- data.frame()
  
  for (group in groups) {
    
    filename <- sprintf('data/%s_loc_AOV.csv',group)
    
    if (prod(dim(df)) == 0) {
      
      df <- read.csv(filename, stringsAsFactors = TRUE)
      
    } else {
      
      df <- rbind(df, read.csv(filename, stringsAsFactors = TRUE))
      
    }
    
  }
  
  library(data.table)
  setnames(df, old = c('participant','rotated_b','passive_b','handangle_deg'), new = c('subject','rotated','passive','angle'))
  
  
  dfs <- aggregate(cbind(group, bias_deg) ~ rotated + passive + subject, data=df, FUN=mean)
  dfs$group <- factor(dfs$group, levels=c(1,2), labels=c('control', 'EDS'))
  
  
  factorvars <- c('subject', 'rotated', 'passive')
  
  for (fv in factorvars) {
    
    dfs[,fv] <- as.factor(dfs[,fv])
    
  }
  
  # full ANOVA:
  print('*** full ANOVA')
#  print(ezANOVA(data=df,dv=bias_deg,wid=subject,within=c(rotated,passive,angle),between=group,type=3))
  print(ez::ezANOVA(data=dfs,dv=bias_deg,wid=subject,within=c(rotated,passive),between=group,type=3))
  
  # ANOVA on change in localization:
  # for (passive in c(0,1)) {
  #   print(c('*** active','*** passive')[passive+1])
  #   print(ezANOVA(data=df[which(df$passive == passive),],dv=bias_deg,wid=subject,diff=rotated,within=c(angle),between=group,type=3))
  # }
  # 
  # group_col <- c()
  # subject_col <- c()
  # angle_col <- c()
  # bias_col <- c()
  # 
  # for (participant in unique(df$subject)) {
  #   
  #   group <- df$group[df$subject == participant][1]
  #   
  #   for (angle in unique(df$angle)) {
  #     
  #     act_al <- df$bias[(df$subject == participant) & (df$angle == angle) & (df$rotated == 0) & (df$passive == 0)]
  #     act_ro <- df$bias[(df$subject == participant) & (df$angle == angle) & (df$rotated == 1) & (df$passive == 0)]
  #     pas_al <- df$bias[(df$subject == participant) & (df$angle == angle) & (df$rotated == 0) & (df$passive == 1)]
  #     pas_ro <- df$bias[(df$subject == participant) & (df$angle == angle) & (df$rotated == 1) & (df$passive == 1)]
  #     
  #     act_bias_shift <- act_ro - act_al
  #     pas_bias_shift <- pas_ro - pas_al
  #     
  #     bias <- act_bias_shift - pas_bias_shift
  #     
  #     # add to lists:
  #     group_col <- c(group_col, group)
  #     subject_col <- c(subject_col, participant)
  #     angle_col <- c(angle_col, angle)
  #     bias_col <- c(bias_col, bias)
  #     
  #   }
  #   
  # }
  # 
  # df2 <- data.frame(group_col,subject_col,angle_col,bias_col)
  # colnames(df2) <- c('group', 'subject', 'angle', 'bias_deg')
  # factorvars <- c('group', 'subject', 'angle')
  # 
  # for (fv in factorvars) {
  #   
  #   df2[,fv] <- as.factor(df2[,fv])
  #   
  # }
  # 
  # 
  # # ANOVA on change in predicted sensory consequences
  # print('*** predicted consequences')
  # print(ezANOVA(data=df2,dv=bias_deg,wid=subject,within=angle,between=group,type=3))
  
  dfd <- aggregate(bias_deg ~ passive + subject, data=dfs, FUN=diff)
  dfd$group <- dfs$group[which(dfs$passive == 0)]
  
  # ANOVA on localization shifts
  print('*** localization shifts')
  print(ez::ezANOVA(data=dfd,dv=bias_deg,wid=subject,between=group,within=passive,type=3))
  
  
}