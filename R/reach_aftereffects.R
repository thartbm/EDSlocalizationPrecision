
source('R/common.R')

# Data handling no-cursors -----

getAllNoCursors <- function(groups=c('sEDS','zEDS')) {
  
  # loop through groups:
  for (group in groups) {
    
    # we store the dataframes in this list:
    groupdf <- getGroupNoCursors(group)
    filename <- sprintf('data/%s_nocursors.csv',group)
    write.csv(filename,x=groupdf,row.names=FALSE,quote=FALSE)
  }
  
}


getGroupNoCursors <- function(group) {
  
  participants <- getGroupParticipants(group)
  
  # make a dataframe that holds the relevant info
  GroupNoCursors <- data.frame(matrix(data=NA,nrow=length(participants),ncol=4))
  colnames(GroupNoCursors) <- c('participant','aligned','exclusive','inclusive')
  GroupNoCursors$participant <- participants
  
  # loop through participants:
  for (pp in participants) {
    
    # get each participants' after effects and store in dataframe:
    GroupNoCursors[which(GroupNoCursors$participant == pp), c('aligned','exclusive','inclusive')] <- getParticipantNoCursors(pp)
    
  }
  
  # return the after effects dataframe:
  return(GroupNoCursors)
  
}

getParticipantNoCursors <- function(participant) {
  
  pdf <- read.csv('data/participants_files.csv',stringsAsFactors=FALSE)
  
  folder <- pdf[pdf$ID==participant, 'folder']
  version <- pdf[pdf$ID==participant, 'version']
  
  alignednocursor <- pdf[pdf$ID==participant,'aligned_nocursor']
  alignednocursor <- sprintf('data/%s/%s/%s',folder,participant,alignednocursor)
  alignedNoCursor <- getNoCursors(alignednocursor, version='aligned')
  
  # aligned reaches need to be boiled down to simple biases, one for each target angle
  # for now we will use the final one third of each set
  alignedBiases <- getAlignedBiases(alignedNoCursor)
  # the getAlignedBiases function should work, if the dataframe has the same column names...
  # otherwise, it might actually be a one-liner, with aggregate()
  
  rotatednocursor <- pdf[pdf$ID==participant,'rotated_nocursor']
  rotatednocursor <- sprintf('data/%s/%s/%s',folder,participant,rotatednocursor)
  rotatedNoCursor <- getNoCursors(rotatednocursor,version=version)
  
  # get three numbers: aligned / exclusive / inclusive
  reachAngles <- aggregate(as.formula('angular_deviation ~ set'), rotatedNoCursor, median, na.rm=TRUE)
  reachAngles <- c(mean(alignedBiases$angular_deviation), c(reachAngles$angular_deviation))
  
  return(reachAngles)
  
}

getNoCursors <- function(filename,version='aligned') {
  
  # returns a dataframe with columns:
  # 1) instruction (0=excluding/aligned, 1=including)
  # 2) trial number
  # 3) angular deviation
  # 4) target angle
  
  print(filename)
  # for each reach, this should return the angular deviation, as well as the target angle
  
  # read the data file:
  columnnames <- c("trial", "target", "Xcursor","Ycursor","Xrobot","Yrobot","XscreenO","YscreenO","Xhome","Yhome","XTarget","Ytarget","blocknumber","rotation","time","trialselected","sampleselected","sampleinterpolated","maxvelocity")
  reachdf <- read.table(filename,stringsAsFactors=FALSE);
  
  if (dim(reachdf)[2] == 19) {
    colnames(reachdf) <- columnnames
  }
  if (dim(reachdf)[2] == 20) {
    colnames(reachdf) <- c(columnnames, 'unsure')
  }
  
  # fix, so that origin is (0,0), at least for the robot coordinates that we use here:
  reachdf$Yrobot = reachdf$Yrobot + 8.5
  
  # we go through all trials in the file
  trials <- unique(reachdf$trial)
  
  # we need to know what kind of trials they are
  instruction <- 1:length(trials)
  
  if (version == 'aligned') {
    
    instruction[] <- 0
    
  } else {
    
    if (version == 'IE') {
      
      instructions <- c(1,0)
      
    }
    
    if (version == 'EI') {
      
      instructions <- c(0,1)
      
    }
    
    instruction <- c()
    trials <- c()
    
    for (block in 1:4) {
      
      starttrial = c(187,331,475,619)[block] - 1
      
      for (set in 1:2) {
        
        newinstruction <- 1:9
        newinstruction[] <- instructions[set]
        instruction <- c(instruction, newinstruction)
        
        newtrial <- (1:9) + starttrial
        starttrial <- starttrial + 9
        trials <- c(trials, newtrial)
        
      }
      
      instructions <- rev(instructions)
      
    }
    
  }
  
  # make matrix/dataframe to put the data in:
  ReachAngles = matrix(data=NA,nrow=length(trials),ncol=4)
  # first column is set:
  ReachAngles[,1] <- instruction
  # second columns is trial number:
  ReachAngles[,2] <- trials
  
  for (trial.id in 1:length(trials)) {
    
    trialno <- trials[trial.id]
    
    # only for existing trials:
    if (length(which(reachdf$trial == trialno)) > 0) {
      
      trialdf <- reachdf[reachdf$trial == trialno,]
      
      # ReachAngles[trial.id,3:4] <- getTrialReachAngleAt3cm(trialdf)
      # ReachAngles[trial.id,3:4] <- getTrialReachAngleAtEndPoint(trialdf)
      # ReachAngles[trial.id,3:4] <- getTrialReachAngleAt(trialdf, location='cm4')
      ReachAngles[trial.id,3:4] <- getTrialReachAngleAt(trialdf, location='maxvel')
      
    }
    
  }
  
  ReachAngles <- data.frame(ReachAngles)
  colnames(ReachAngles) <- c('set','trial','angular_deviation','target_angle')
  
  return(ReachAngles)
  
}

# Figures -----

plotReachAftereffects <- function(target='inline') {
  
  if (target == 'svg') {
    svglite::svglite(file='doc/Fig4.svg', width=8, height=3, system_fonts=list(sans='Arial'))
  }
  
  styles <- getStyle()
  
  #par(mfrow=c(1,1), mar=c(4,4,2,0.1))
  
  par(mar=c(4,4,2,0.1))
  
  
  layout(matrix(c(1,2), nrow=1, ncol=2, byrow = TRUE), widths=c(2,3), heights=c(1,1))
  
  
  ylims=c(-.1*max(styles$rotation),max(styles$rotation)+(.2*max(styles$rotation)))
  plot(c(0.8,2.2),c(0,0),type='l',lty=2,col=rgb(.5,.5,.5),xlim=c(0.75,2.25),ylim=ylims,bty='n',
       xaxt='n',yaxt='n',xlab='strategy use',ylab='reach deviation [°]',main='',font.main=1)
  
  mtext('A', side=3, outer=TRUE, at=c(0,1), line=-1, adj=0, padj=1)
  
  for (groupno in c(1:length(styles$group))) {
    
    group <- styles$group[groupno]
    
    reachaftereffects <- read.csv(sprintf('data/%s_nocursors.csv',group), stringsAsFactors=FALSE)
    
    reachaftereffects$exclusive <- reachaftereffects$exclusive - reachaftereffects$aligned
    reachaftereffects$inclusive <- reachaftereffects$inclusive - reachaftereffects$aligned
    
    meanExc <- mean(reachaftereffects$exclusive)
    meanInc <- mean(reachaftereffects$inclusive)
    
    coord.x <- c(1,1,2,2)
    coord.y <- c(t.interval(reachaftereffects$exclusive),rev(t.interval(reachaftereffects$inclusive)))
    polygon(coord.x, coord.y, col=as.character(styles$color_trans[groupno]), border=NA)
    
  }
  
  for (groupno in c(1:length(styles$group))) {
    
    group <- styles$group[groupno]
    offset <- (groupno - ((length(styles$group) - 1) / 2)) * .035
    
    reachaftereffects <- read.csv(sprintf('data/%s_nocursors.csv',group), stringsAsFactors=FALSE)
    

    reachaftereffects$exclusive <- reachaftereffects$exclusive - reachaftereffects$aligned
    reachaftereffects$inclusive <- reachaftereffects$inclusive - reachaftereffects$aligned
    
    meanExc <- mean(reachaftereffects$exclusive)
    meanInc <- mean(reachaftereffects$inclusive)
    
    lines(c(1,2),c(meanExc,meanInc),col=as.character(styles$color_solid[groupno]),lty=styles$linestyle[groupno],lw=2)
    
  }
  
  axis(side=1, at=c(1,2), labels=c('without strategy','with strategy'),cex.axis=1.00)
  if (max(styles$rotation) == 30) {
    axis(side=2, at=c(0,10,20,30),cex.axis=1.00)
  }
  
  # legend(0.5,max(styles$rotation)*(7/6),styles$label,col=as.character(styles$color),lty=styles$linestyle,bty='n',cex=0.85)
  legend(0.8,30,styles$label,col=as.character(styles$color_solid),lw=2,lty=styles$linestyle,bty='n',cex=1.00,seg.len = 3)
  
  
  plot(c(0.4,3.6),c(0,0),type='l',lty=2,col=rgb(.5,.5,.5),xlim=c(0.5,3.5),ylim=ylims,bty='n',
       xaxt='n',yaxt='n',xlab='strategy use',ylab='reach deviation [°]',main='',font.main=1)
  
  
  mtext('B', side=3, outer=TRUE, at=c(2/5,1), line=-1, adj=0, padj=1)
  
  for (groupno in c(1:length(styles$group))) {
    
    group <- styles$group[groupno]
    #print(group)
    reachaftereffects <- read.csv(sprintf('data/%s_nocursors.csv',group), stringsAsFactors=FALSE)
    #print(str(reachaftereffects))
    conditions <- c('aligned','exclusive','inclusive')
    for (conditionno in c(1:length(conditions))) {
      
      #print(conditionno)
      
      condition <- conditions[conditionno]
      
      nocursors <- as.numeric(reachaftereffects[,condition])
      #print(str(nocursors))
      X <- rep((conditionno+((groupno-1.7)/2.5)),length(nocursors))
      Y <- as.numeric(nocursors)
      points(x=X,y=Y,pch=16,cex=1.5,col=as.character(styles$color_trans[groupno]))
      
      meandist <- getConfidenceInterval(data=c(nocursors), method='bootstrap', resamples=5000, FUN=mean, returndist=TRUE)
      
      DX <- meandist$density$x
      DY <- meandist$density$y / max(meandist$density$y) / 6
      
      DX <- c(DX[1], DX, DX[length(DX)])
      DY <- c(0,     DY, 0)
      
      Xoffset <- (conditionno+((groupno-1.4)/2.5))
      
      polygon(x=DY+Xoffset, y=DX, border=FALSE, col=as.character(styles$color_trans[groupno]))
      
      lines(x=rep(Xoffset,2),y=meandist$CI95,col=as.character(styles$color_solid[groupno]))
      #print(meandist$CI95)
      points(x=Xoffset,y=mean(c(nocursors)),pch=16,cex=1.5,col=as.character(styles$color_solid[groupno]))
    
    }
    
  }
  
  axis(side=1, at=c(1,2,3),labels=conditions)
  axis(side=2, at=c(0,10,20,30),labels=c('0','10','20','30'),cex.axis=1.00)
  
  if (target == 'svg') {
    dev.off()
  }
  
}

# Statistics -----

getRAE4ANOVA <- function(styles) {
  
  group       <- c()
  participant    <- c()
  strategy       <- c()
  reachdeviation <- c()
  
  # keeping count of unique participants:
  startingID <- 0
  
  for (groupno in c(1:length(styles$group))) {
    
    groupname <- styles$group[groupno]
    
    df <- read.csv(sprintf('data/%s_nocursors.csv',groupname),stringsAsFactors=F)
    
    df$exclusive <- df$exclusive - df$aligned
    df$inclusive <- df$inclusive - df$aligned
    
    # we need to know the number of participants to replicate some values:
    N <- dim(df)[1]
    
    for (thisstrategy in c('exclusive','inclusive')) {
      
      group           <- c(group, rep(as.character(styles$label[groupno]), N))
      participant     <- c(participant, c(startingID : (startingID + N - 1)))
      strategy        <- c(strategy, rep(thisstrategy, N))
      reachdeviation  <- c(reachdeviation, df[,thisstrategy])
      
    }
    
    startingID <- startingID + N
    
  }
  
  # put it in a data frame:
  RAEaov <- data.frame(group, participant, strategy, reachdeviation)
  
  # set relevant columns as factors:
  RAEaov$group <- as.factor(RAEaov$group)
  RAEaov$strategy <- as.factor(RAEaov$strategy)
  
  return(RAEaov)
  
}

RAE.ANOVA <- function() {
  
  styles <- getStyle()
  
  RAE4aov <- getRAE4ANOVA(styles)                      
  
  #learning curve ANOVA's
  # for ez, case ID should be a factor:
  RAE4aov$participant <- as.factor(RAE4aov$participant)
  print(ezANOVA(data=RAE4aov, wid=participant, dv=reachdeviation, within=strategy, between=c(group),type=3))
  
  # using afex:
  #RAE4aov$group <- as.factor(RAE4aov$group)
  #RAE4aov$strategy <- as.factor(RAE4aov$strategy)
  #afex::aov_car(reachdeviation ~ group * strategy + Error(participant/strategy), data=RAE4aov, type=3)
  
  # this seems to work too:
  #library('car')
  #mod <- lm(reachdeviation ~ group + strategy - participant, data=RAE4aov)
  #Anova(mod, type=3)
  
}

NoCursorANOVA <- function() {
  
  styles <- getStyle()
  
  NC4aov <- getNoCursors4ANOVA(styles)
  
  NC4aov$participant <- as.factor(NC4aov$participant)
  print(ezANOVA(data=NC4aov, wid=participant, dv=reachdeviation, within=training, between=c(group),type=3))
  
}

getNoCursors4ANOVA <- function(styles) {
  
  # placeholder for data frame:
  NC4aov <- NA
  
  # loop through groups to collect their data:
  for (groupno in c(1:length(styles$group))) {
    
    groupname <- styles$group[groupno]
    
    
    df <- read.csv(sprintf('data/%s_nocursors.csv',groupname),stringsAsFactors=F)
    
    AL.NC <- df[,c('participant','aligned')]
    colnames(AL.NC)[2] <- 'reachdeviation'
    AL.NC$training <- 'aligned'
    
    RO.NC <- df[,c('participant','exclusive')]
    colnames(RO.NC)[2] <- 'reachdeviation'
    RO.NC$training <- 'rotated'
    
    df <- rbind(AL.NC, RO.NC)
    df$group <- groupname
    
    if (is.data.frame(NC4aov)) {
      NC4aov <- rbind(NC4aov, df)
    } else {
      NC4aov <- df
    }
    
  }
  
  NC4aov$group <- as.factor(NC4aov$group)
  NC4aov$training <- as.factor(NC4aov$training)
  
  return(NC4aov)
  
}


NoCursorTtests <- function() {
  
  styles <- getStyle()
  
  RAE4aov <- getRAE4ANOVA(styles)
  
  RAE4aov$participant <- as.factor(RAE4aov$participant)
  
  EDSwis <- RAE4aov$reachdeviation[which(RAE4aov$group == 'EDS' & RAE4aov$strategy == 'inclusive')]
  EDSwos <- RAE4aov$reachdeviation[which(RAE4aov$group == 'EDS' & RAE4aov$strategy == 'exclusive')]
  CTRwis <- RAE4aov$reachdeviation[which(RAE4aov$group == 'control' & RAE4aov$strategy == 'inclusive')]
  CTRwos <- RAE4aov$reachdeviation[which(RAE4aov$group == 'control' & RAE4aov$strategy == 'exclusive')]
  
  cat('WITH strategy in both groups:\n')
  
  print(t.test(EDSwis, CTRwis, alternative = "less"))
  
  cat(sprintf('eta-squared: %0.5f\n\n', etaSquaredTtest(EDSwis, EDSwos)))
  
  
  cat('Control group:\n')

  print(t.test(CTRwos, CTRwis, alternative = "less"))

  cat(sprintf('eta-squared: %0.5f\n\n', etaSquaredTtest(CTRwos, CTRwis)))

  cat('EDS group:\n')

  print(t.test(EDSwos, EDSwis, alternative = "greater"))

  cat(sprintf('eta-squared: %0.5f\n', etaSquaredTtest(EDSwos, EDSwis)))
  
  
}
