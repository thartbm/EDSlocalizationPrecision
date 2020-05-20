source('R/common.R')

# Data handling -----

getLearningCurves <- function(groups = c('sEDS','zEDS')) {
  
  for (groupno in 1:length(groups)) {
    
    grouplearningcurves = getGroupLearningCurves(groups[groupno])
    
    filename <- sprintf('data/%s_curves.csv',groups[groupno])
    
    write.csv(filename,x=grouplearningcurves,row.names=FALSE,quote=FALSE)
    
  }

}

getGroupLearningCurves<- function(group) {
  
  # determine which participants are in this group, based on the main participants.csv
  participants <- getGroupParticipants(group)
  
  # make a dataframe that holds the relevant info
  GroupLearningCurves <- data.frame(matrix(data=NA,nrow=90,ncol=length(participants)))
  colnames(GroupLearningCurves) <- participants
  
  # loop through participants:
  for (pp in participants) {
    
    # get each participants' learning curve and store in dataframe:
    GroupLearningCurves[,pp] <- getParticipantLearningCurve(pp)
    
  }
  
  # return the learning curve dataframe:
  return(GroupLearningCurves)
  
}

getParticipantLearningCurve <- function(participant) {
  
  pdf <- read.csv('data/participants_files.csv',stringsAsFactors=FALSE)
  
  folder <- pdf[pdf$ID==participant, "folder"]
  
  alignedreachfile <- pdf[pdf$ID==participant,"aligned_train"]
  alignedreachfile <- sprintf('data/%s/%s/%s',folder,participant,alignedreachfile)
  alignedReachAngles <- getReachAngles(alignedreachfile)
  
  # aligned reaches need to be boiled down to simple biases, one for each target angle
  # for now we will use the final one third of each set
  alignedBiases <- getAlignedBiases(alignedReachAngles)
  
  #str(alignedReachAngles$trial)
  
  rotatedreachfile <- pdf[pdf$ID==participant,"rotated_train"]
  rotatedreachfile <- sprintf('data/%s/%s/%s',folder,participant,rotatedreachfile)
  rotatedReachAngles <- getReachAngles(rotatedreachfile)
  
  # for now, we are only interested in the initial training block of 90 trials:
  rotatedReachAngles <- rotatedReachAngles[rotatedReachAngles$trial < 91,]
  
  # this seems unnecessary:
  AD_idx <- which(is.finite(alignedBiases$angular_deviation))
  
  for (target_angle in alignedBiases$target_angle) {
    
    # aligned bias for current target angle:
    bias <- alignedBiases$angular_deviation[alignedBiases$target_angle == target_angle]
    
    # get rows with current target angle
    row_idx <- which(rotatedReachAngles$target_angle == target_angle)
    row_idx <- intersect(row_idx, AD_idx)
    
    rotatedReachAngles$angular_deviation[row_idx] <- rotatedReachAngles$angular_deviation[row_idx] - bias
  }
  
  ParticipantLearningCurve = matrix(data=rotatedReachAngles$angular_deviation,nrow=90,ncol=1)
  
  return(ParticipantLearningCurve)
  
}

getReachAngles <- function(filename) {
  print(filename)
  # for each reach, this should return the angular deviation, as well as the target angle
  
  # read the data file:
  columnnames <- c("trial", "target", "Xcursor","Ycursor","Xrobot","Yrobot","XscreenO","YscreenO","Xhome","Yhome","XTarget","Ytarget","blocknumber","rotation","time","trialselected","sampleselected","sampleinterpolated","maxvelocity","unsure")
  reachdf <- read.table(filename,stringsAsFactors=FALSE, col.names=columnnames, row.names = NULL);
  
  
  # fix, so that origin is (0,0), at least for the robot coordinates that we use here:
  reachdf$Yrobot = reachdf$Yrobot + 8.5
  
  # decide if it is aligned or rotated:
  
  if (regexpr('align',filename)[1] == -1) {
    condition <- 'rotated'
    # rotated file, use one set of trial numbers
    trials <- c(1:90, 109:138, 157:186, 205:234, 253:282, 301:330, 349:387, 397:426, 445:474, 493:522, 541:570, 589:618)
    #print(length(unique(reachdf$trial[which(reachdf$trial %in% trials)])))
    # Block 1
    # 1 – 90: trainingc
    # 109 – 138: training
    # 157 – 186: training
    # 
    # Block 2
    # 205 – 234: training
    # 253 – 282: training
    # 301 – 330: training
    # 
    # Block 3
    # 349 – 378: training
    # 397 – 426: training
    # 445 – 474: training
    # 
    # Block 4
    # 493 – 522: training
    # 541 – 570: training
    # 589 – 618: training
  } else {
    condition <- 'aligned'
    # aligned file, use other trial numbers
    trials <- c(1:45, 64:72, 91:99, 109:117, 136:144, 163:171, 181:189, 208:216, 235:243, 253:261, 280:288, 307:315)
    
    # Block 2
    # 1 – 45: training
    # 64 – 72: training
    # 91 – 99: training
    # 
    # Block 2
    # 109 – 117: training
    # 136 – 144: training
    # 163 – 171: training
    # 
    # Block 3
    # 181 – 189: training
    # 208 – 216: training
    # 235 – 243: training
    # 
    # Block 4
    # 253 – 261: training
    # 280 – 288: training
    # 307 – 315: training
  }
  
  # provide set indexing, so we can look at the last 1/3 of each set...
  # there has to be a set number for every trial:
  set <- trials
  
  # sets are defined by gaps in trial numbers, larger than 1
  setgaps <- which(diff(trials) > 1)
  # we need to know the previous setgap:
  previous <- 0
  
  # loop through setgaps:
  for (setgapn in 1:length(setgaps)) {
    
    setgap <- setgaps[setgapn]
    
    # set numbers are constant in between gaps
    set[previous+1:setgap] <- setgapn
    # update the previous setgap for the next loop iteration
    previous <- setgap
    
  }
  
  # final set is after last gap:
  set[previous+1:length(trials)] <- length(setgaps) + 1
  
  # sometimes there are extra trials for some reason, cut off at appropriate length:
  set <- set[1:length(trials)]
  
  # now, if stuff is aligned, we want only the final third of the trials
  # ans for rotated, we only want the first set of 90
  if (condition == 'aligned') {
    
    newtrials <- c()
    newset <- c()
    
    for (set.id in unique(set)) {
      # which trials are in the set?
      settrials <- trials[which(set == set.id)]
      # how many is one third of those trials?
      Ntrials <- length(settrials)
      Ntrials <- Ntrials / 3
      # get those trial numbers:
      settrials <- settrials[length(settrials)-Ntrials+1:length(settrials)]
      settrials <- settrials[1:Ntrials]
      # append them to our list:
      newtrials <- c(newtrials,settrials)
      
      # also make new set list...
      setset <- 1:Ntrials
      setset[] <- set.id
      newset <- c(newset,setset)
    }
    
    trials <- newtrials
    set <- newset
  }
  
  if (condition == 'rotated') {
    
    trials <- trials[which(set == 1)]
    set <- 1:length(trials)
    set[] <- 1
    
  }
  # now start collecting the reach angles
  
  # first make a matrix to store stuff
  ReachAngles = matrix(data=NA,nrow=length(trials),ncol=4)
  # first column is set:
  ReachAngles[,1] <- set
  # second columns is trial number:
  ReachAngles[,2] <- trials
  
  # loop through all trials to get 
  for (trial.idx in 1:length(trials)) {
    
    trialno = trials[trial.idx]
    
    # we need to check if the trial is actually in the data
    if (length(which(reachdf$trial == trialno)) > 0) {
      
      trialdf <- reachdf[reachdf$trial == trialno,]
      ReachAngles[trial.idx,3:4] <- getTrialReachAngleAt(trialdf, location='maxvel')
      
    }
    # no else: if the trial is not there, the values remain NA
    
  }
  
  # convert to data frame
  ReachAngles <- data.frame(ReachAngles)
  
  # set columns names so other functions know what to do whith this info:
  colnames(ReachAngles) <- c('set','trial','angular_deviation','target_angle')
  
  return(ReachAngles)
  
}



getBlockedLearningCurves <- function(group, blockdefs) {
  
  curves <- read.csv(sprintf('data/%s_curves.csv',group), stringsAsFactors=FALSE)  
  
  N <- dim(curves)[2]
  
  blocked <- array(NA, dim=c(N,length(blockdefs)))
  
  for (ppno in c(1:N)) {
    
    for (blockno in c(1:length(blockdefs))) {
      
      blockdef <- blockdefs[[blockno]]
      blockstart <- blockdef[1]
      blockend <- blockstart + blockdef[2] - 1
      samples <- curves[blockstart:blockend,ppno]
      blocked[ppno,blockno] <- mean(samples, na.rm=TRUE)

    }
    
  }
  
  return(blocked)
  
}

getReachPrecision <- function(groups = c('sEDS','zEDS')) {
  
  for (groupno in 1:length(groups)) {
    
    groupReachPrecision = getGroupReachPrecision(groups[groupno])
    
    filename <- sprintf('data/%s_training_var.csv',groups[groupno])
    
    write.csv(filename,x=groupReachPrecision,row.names=FALSE,quote=FALSE)
    
  }
  
}

getGroupReachPrecision <- function(group) {
  
  # determine which participants are in this group, based on the main participants.csv
  participants <- getGroupParticipants(group)
  
  # for every participant, we'll get the reach precision in the last 15 trials 
  # of the first aligned and rotated training session
  aligned <- c()
  rotated <- c()
  
  for (participant in participants) {
    
    pdf <- read.csv('data/participants_files.csv',stringsAsFactors=FALSE)
    
    folder <- pdf[pdf$ID==participant, "folder"]
    
    columnnames <- c("trial", "target", "Xcursor","Ycursor","Xrobot","Yrobot","XscreenO","YscreenO","Xhome","Yhome","XTarget","Ytarget","blocknumber","rotation","time","trialselected","sampleselected","sampleinterpolated","maxvelocity","unsure")
    
    alignedreachfile <- pdf[pdf$ID==participant,"aligned_train"]
    alignedreachfile <- sprintf('data/%s/%s/%s',folder,participant,alignedreachfile)
    
    alignedreachdf <- read.table(alignedreachfile,stringsAsFactors=FALSE, col.names=columnnames, row.names = NULL);
    alignedreachdf <- alignedreachdf[which(alignedreachdf$trial > 30 & alignedreachdf$trial < 46),]
    
    alignedreachdeviations <- getDFreachAngles(alignedreachdf)
    alignedreachdeviations <- removeTargetBias(alignedreachdeviations)
    
    aligned <- c(aligned, sd(alignedreachdeviations$angular_deviation, na.rm=TRUE))
    
    #getTrialReachAngleAt(trialdf, location='maxvel')
    
    rotatedreachfile <- pdf[pdf$ID==participant,"rotated_train"]
    rotatedreachfile <- sprintf('data/%s/%s/%s',folder,participant,rotatedreachfile)
    
    rotatedreachdf <- read.table(rotatedreachfile,stringsAsFactors=FALSE, col.names=columnnames, row.names = NULL);
    rotatedreachdf <- rotatedreachdf[which(rotatedreachdf$trial > 75 & rotatedreachdf$trial < 91),]
    
    rotatedreachdeviations <- getDFreachAngles(rotatedreachdf)
    rotatedreachdeviations <- removeTargetBias(rotatedreachdeviations)
    
    rotated <- c(rotated, sd(rotatedreachdeviations$angular_deviation, na.rm=TRUE))
    
  }
  
  reachPrecision <- data.frame('participant'=participants, aligned, rotated)
  
  return(reachPrecision)
  
}

getDFreachAngles <- function(reachdf) {
  
  reachdf$Yrobot = reachdf$Yrobot + 8.5
  
  trials <- unique(reachdf$trial)
  
  ReachAngles = matrix(data=NA,nrow=length(trials),ncol=2)
  
  # loop through all trials to get 
  for (trial.idx in 1:length(trials)) {
    
    trialno <- trials[trial.idx]
    
    # we need to check if the trial is actually in the data
    if (length(which(reachdf$trial == trialno)) > 0) {
      
      trialdf <- reachdf[reachdf$trial == trialno,]
      ReachAngles[trial.idx,1:2] <- getTrialReachAngleAt(trialdf, location='maxvel')
      
    }
    # no else: if the trial is not there, the values remain NA
    
  }
  
  # convert to data frame
  ReachAngles <- data.frame(ReachAngles)
  
  # set columns names so other functions know what to do whith this info:
  colnames(ReachAngles) <- c('angular_deviation','target_angle')
  
  return(ReachAngles)
  
}

removeTargetBias <- function(reachdf) {
  
  targets <- unique(reachdf$target_angle)
  
  for (target in targets) {
    
    idx <- which(reachdf$target_angle == target)
    reachdf$angular_deviation[idx] <- reachdf$angular_deviation[idx] - mean(reachdf$angular_deviation[idx])
    
  }
  
  return(reachdf)
  
}

# Figures -----

plotLearningCurves <- function(target='inline') {
  
  styles <- getStyle()
  
  if (target == 'svg') {
    svglite::svglite(file='doc/Fig3.svg', width=7.5, height=5.5, system_fonts=list(sans='Arial'))
  }
  
  #par(mfrow=c(1,2), mar=c(4,4,2,0.1))
  par(mar=c(4,4,2,0.1))
  
  
  layout(matrix(c(1,2,3,4), nrow=2, ncol=2, byrow = TRUE), widths=c(1,1), heights=c(1,1))
  
  # # # # # # # # # #
  # panel A: actual learning curves
  
  ylims=c(-.1*max(styles$rotation),max(styles$rotation)+(.2*max(styles$rotation)))
  plot(c(-1,36),c(0,0),col=rgb(0.5,0.5,0.5),type='l',lty=2,xlim=c(-1,36),ylim=ylims,xlab='trial',ylab='reach deviation [°]',xaxt='n',yaxt='n',bty='n',main='',font.main=1)
  
  #mtext('A', side=3, outer=TRUE, at=c(0,1), line=-1, adj=0, padj=1)
  mtext('A', outer=FALSE, side=3, las=1, line=1, adj=0, padj=1)
  
  for (groupno in c(1:length(styles$group))) {
    
    group  <- styles$group[groupno]
    curves <- read.csv(sprintf('data/%s_curves.csv',group), stringsAsFactors=FALSE)  
    curve  <- apply(curves, c(1), mean, na.rm=T)
    
    lines(c(1:15),curve[1:15],col=as.character(styles$color_solid[groupno]),lty=styles$linestyle[groupno],lw=2)
    
    lines(c(21:35),curve[76:90],col=as.character(styles$color_solid[groupno]),lty=styles$linestyle[groupno],lw=2)
  }
  
  # axis(side=1, at=c(1,10,20,30))
  axis(side=1, at=c(1,5,10,25,30,35), labels=c('1','5','10','80','85','90'),cex.axis=1.00)
  axis(side=2, at=c(0,10,20,30),cex.axis=1.00)
  
  legend(20,15,styles$label,col=as.character(styles$color_solid),lty=styles$linestyle,bty='n',lw=2,cex=1.00,seg.len=3)
  
  
  # # # # # # # # # #
  # panel B: blocked learning curves
  
  plot(c(0,5),c(0,0),col=rgb(0.5,0.5,0.5),type='l',lty=2,xlim=c(0.5,4.5),ylim=ylims,xlab='trial set',ylab='reach deviation [°]',xaxt='n',yaxt='n',bty='n',main='',font.main=1)
  
  #mtext('B', side=3, outer=TRUE, at=c(3/7,1), line=-1, adj=0, padj=1)
  mtext('B', outer=FALSE, side=3, las=1, line=1, adj=0, padj=1)
  
  blockdefs <- list(c(1,3),c(4,3),c(76,15))
  
  alllines <- array(NA,dim=c(length(styles$group),length(blockdefs)))
  allpolys <- array(NA,dim=c(length(styles$group),2*length(blockdefs)))
  
  for (groupno in c(1:length(styles$group))) {
    
    group <- styles$group[groupno]
    
    blocked <- getBlockedLearningCurves(group, blockdefs)
    
    alllines[groupno,] <- apply(blocked, c(2), mean, na.rm=T)
    
    blockedCI <- apply(blocked, c(2), t.interval)
    allpolys[groupno,] <- c(blockedCI[1,], rev(blockedCI[2,]))
    
  }
  
  # first plot all the polygons representing confidence intervals, so those are in the background
  
  for (groupno in c(1:length(styles$group))) {
    
    polX <- c(c(1,2,4),rev(c(1,2,4)))
    polY <- allpolys[groupno,]
    polygon(polX,polY,col=as.character(styles$color_trans[groupno]),border=NA)
    
  }
  
  # then plot the lines representing the means, so those are in the foreground
  
  for (groupno in c(1:length(styles$group))) {
    
    lines(c(1,2,4),alllines[groupno,],col=as.character(styles$color_solid[groupno]),lty=styles$linestyle[groupno],lw=2)
    
  }
  
  # legend(2,0.45,styles$label,col=as.character(styles$color),lty=styles$linestyle,bty='n',cex=0.7)
  
  axis(side=1, at=c(1,2,4), labels=c('1-3','4-6','76-90'),cex.axis=1.00)
  axis(side=2, at=c(0,10,20,30),labels=c('0','10','20','30'),cex.axis=1.00)
  
  
  # # # # # # # # # #
  # panel C: individual participants in the first trial set
  
  plot(c(0,5),c(0,0),col=rgb(0.5,0.5,0.5),type='l',lty=2,xlim=c(0.5,6.5),ylim=ylims,xlab='trial set',ylab='individual reach deviation [°]',xaxt='n',yaxt='n',bty='n',main='',font.main=1)
  
  #mtext('C', side=3, outer=TRUE, at=c(5/7,1), line=-1, adj=0, padj=1)
  mtext('C', outer=FALSE, side=3, las=1, line=1, adj=0, padj=1)
  
  blockdefs <- list(c(1,3),c(4,3),c(76,15))
  
  for (groupno in c(1:length(styles$group))) {
    
    group <- styles$group[groupno]
    
    blocked <- getBlockedLearningCurves(group, blockdefs)
    
    for (blockno in c(1,2,3)) {
      
      X <- rep((groupno-(1/3)+((blockno-1)*2)),dim(blocked)[1])
      Y <- c(blocked[,blockno])
      #print(Y)
      points(x=X,y=Y,pch=16,cex=1.5,col=as.character(styles$color_trans[groupno]))
      
      meandist <- getConfidenceInterval(data=c(blocked[,blockno]), method='bootstrap', resamples=5000, FUN=mean, returndist=TRUE)
      
      DX <- meandist$density$x
      DY <- meandist$density$y / max(meandist$density$y) / 2.5
      
      DX <- c(DX[1], DX, DX[length(DX)])
      DY <- c(0,     DY, 0)
      
      polygon(x=DY+groupno+((blockno-1)*2), y=DX, border=FALSE, col=as.character(styles$color_trans[groupno]))
      
      lines(x=rep(groupno+((blockno-1)*2),2),y=meandist$CI95,col=as.character(styles$color_solid[groupno]))
      #print(meandist$CI95)
      points(x=groupno+((blockno-1)*2),y=mean(c(blocked[,blockno])),pch=16,cex=1.5,col=as.character(styles$color_solid[groupno]))
      
    }
  }
  
  
  axis(side=1, at=c(1.5,3.5,5.5),labels=c('1-3','4-6','76-90'),cex.axis=1.00)
  axis(side=2, at=c(0,10,20,30),labels=c('0','10','20','30'),cex.axis=1.00)
  
  # # # # # # # # # #
  # panel D: aligned and rotated reach precision
  
  ylims=c(-.1*max(styles$rotation),(max(styles$rotation)*(2/3))+(.2*(max(styles$rotation))))
  
  plot(c(0.5,4.5),c(0,0),col=rgb(0.5,0.5,0.5),type='l',lty=2,xlim=c(0,5),ylim=ylims,xlab='session',ylab='reach precision, SD [°]',xaxt='n',yaxt='n',bty='n',main='reach precision',font.main=1)
  
  #mtext('A', side=3, outer=TRUE, at=c(0,1), line=-1, adj=0, padj=1)
  mtext('D', outer=FALSE, side=3, las=1, line=1, adj=0, padj=1)
  
  for (groupno in c(1:length(styles$group))) {
    
    group  <- styles$group[groupno]
    reachprecision <- read.csv(sprintf('data/%s_training_var.csv',group), stringsAsFactors=FALSE)  
    
    conditions <- c('aligned','rotated')
    
    for (conditionno in c(1:length(conditions))) {
      
      condition <- conditions[conditionno]
      
      X <- rep((groupno-(1/3)+((conditionno-1)*2)),dim(reachprecision)[1])
      Y <- c(reachprecision[,condition])
      
      points(x=X,y=Y,pch=16,cex=1.5,col=as.character(styles$color_trans[groupno]))
      
      meandist <- getConfidenceInterval(data=Y, method='bootstrap', resamples=5000, FUN=mean, returndist=TRUE)
      
      DX <- meandist$density$x
      DY <- meandist$density$y / max(meandist$density$y) / 2.5
      
      DX <- c(DX[1], DX, DX[length(DX)])
      DY <- c(0,     DY, 0)
      
      polygon(x=DY+groupno+((conditionno-1)*2), y=DX, border=FALSE, col=as.character(styles$color_trans[groupno]))
      
      lines(x=rep(groupno+((conditionno-1)*2),2),y=meandist$CI95,col=as.character(styles$color_solid[groupno]))
      points(x=groupno+((conditionno-1)*2),y=mean(Y),pch=16,cex=1.5,col=as.character(styles$color_solid[groupno]))
      
    }
    
  }
  
  axis(side=1, at=c(1.5,3.5), labels=c('aligned','rotated'),cex.axis=1.00)
  axis(side=2, at=c(0,10,20),cex.axis=1.00)
  
  
  if (target == 'svg') {
    dev.off()
  }
  
}



# Statistics ----


learningCurveANOVA <- function() {
  
  styles <- getStyle()
  blockdefs <- list(c(1,3),c(4,3),c(76,15))
  
  LC4aov <- getLearningCurves4ANOVA(styles, blockdefs)                      
  
  #str(LC4aov)
  
  #learning curve ANOVA's
  # for ez, case ID should be a factor:
  LC4aov$participant <- as.factor(LC4aov$participant)
  print(ez::ezANOVA(data=LC4aov, wid=participant, dv=reachdeviation, within=block,between=c(group),type=3))
  
}

getLearningCurves4ANOVA <- function(styles, blockdefs) {
  
  # set up vectors that will form a data frame for the ANOVA(s):
  group          <- c()
  participant    <- c()
  block          <- c()
  reachdeviation <- c()
  
  # keeping count of unique participants:
  startingID <- 0
  
  for (groupno in c(1:length(styles$group))) {
    
    groupname <- styles$group[groupno]

    # block the data:
    blocked <- getBlockedLearningCurves(groupname, blockdefs)
    # this is now the exact same data as the data that is plotted!
    
    # we need to know the number of participants to replicate some values:
    N <- dim(blocked)[1]
    
    for (blockno in c(1:length(blockdefs))) {
      
      group           <- c(group, rep(groupname, N))
      participant     <- c(participant, c(startingID : (startingID + N - 1)))
      block           <- c(block, rep(blockno, N))
      reachdeviation  <- c(reachdeviation, blocked[,blockno])
      
    }
    
    startingID <- startingID + N
    
  }
  
  # put it in a data frame:
  LCaov <- data.frame(group, participant, block, reachdeviation)
  
  # set relevant columns as factors:
  LCaov$group <- as.factor(LCaov$group)
  LCaov$block <- as.factor(LCaov$block)
  
  return(LCaov)
  
}

# The function below could be used to investigate any block * group interaction found in the first ANOVA, but those do not exist.

blockLearningANOVA <- function(block=1) {
  
  styles <- getStyle()
  blockdefs <- list(c(1,3),c(4,3),c(76,15))
  
  LC4aov <- getLearningCurves4ANOVA(styles, blockdefs)
  
  LC4aov <- LC4aov[which(LC4aov$block == block),]
  
  #learning curve ANOVA's
  # for ez, case ID should be a factor:
  LC4aov$participant <- as.factor(LC4aov$participant)
  print(ez::ezANOVA(data=LC4aov, wid=participant, dv=reachdeviation, between=c(group),type=3))
  
}

getReachPrecision4aov <- function() {
  
  styles <- getStyle()
  
  group          <- c()
  participant    <- c()
  training       <- c()
  reachprecision <- c()
  
  for (groupno in c(1:length(styles$group))) {
    
    groupname  <- styles$group[groupno]
    RP <- read.csv(sprintf('data/%s_training_var.csv',groupname), stringsAsFactors=FALSE)  
    
    N <- dim(RP)[1]
    
    for (session in c('aligned','rotated')) {
      
      group          <- c(group,          rep(groupname,N))
      participant    <- c(participant,    RP$participant)
      training       <- c(training,       rep(session,N))
      reachprecision <- c(reachprecision, RP[,session])
      
    }
    
  }
  
  RP4A <- data.frame(group,participant,training,reachprecision)
  RP4A$group       <- as.factor(RP4A$group)
  RP4A$participant <- as.factor(RP4A$participant)
  RP4A$training    <- as.factor(RP4A$training)
  
  return(RP4A)
  
}

reachPrecisionANOVA <- function() {
  
  RP4aov <- getReachPrecision4aov()
  
  print(ez::ezANOVA(data=RP4aov, wid=participant, dv=reachprecision, between=c(group), within=c(training), type=3))
  
}