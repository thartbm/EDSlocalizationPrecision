
# Data Handling -----


detrendLocalization <- function(groups=c('sEDS','zEDS')) {
  
  detrendedLocalization <- list()
  
  for (group in groups) {
    
    loc <- read.csv(sprintf('data/%s_localization.csv',group), stringsAsFactors = FALSE)
    loc$tap_error <- loc$tap_angle - loc$hand_angle
    loc$predicted_error <- NA
    
    participants <- unique(loc$participant)
    
    for (participant in participants) {
      
      for (rotated in c(0,1)) {
        
        for (passive in c(0,1)) {
          
          idx <- which(loc$participant == participant & loc$rotated == rotated & loc$passive == passive)
          
          subdf <- loc[idx,]
          hand_angle <- subdf$hand_angle
          tap_error  <- subdf$tap_error
          
          predictedError <- c()
          
          for (trialidx in c(1:length(idx))) {
            
            spl <- smooth.spline(x=hand_angle[-trialidx], y=tap_error[-trialidx], spar=0.90, keep.data=F )
            
            thissampleprediction <- predict(spl, x=hand_angle[trialidx])$y
            
            predictedError <- c(predictedError, thissampleprediction)
            
          }
          
          loc$predicted_error[idx] <- predictedError
          
        }
        
      }
      
    }
    
    # end of group loop
    
    detrendedLocalization[[group]] <- loc
    
  }
  
  return(detrendedLocalization)
  
}

getLocalizationVariance <- function() {
  
  DetLoc <- detrendLocalization()
  
  for (group in names(DetLoc)) {
    
    loc <- DetLoc[[group]]
    
    participant <- c()
    rotated <- c()
    passive <- c()
    variance <- c()
    
    participants <- unique(loc$participant)
      
    for (pp in participants) {
        
      for (rot in c(0,1)) {
          
        for (pas in c(0,1)) {
            
          idx <- which(loc$participant == pp & loc$rotated == rot & loc$passive == pas)
            
          subdf <- loc[idx,]
          pcVar <- var(subdf$tap_error - subdf$predicted_error)
            
          participant <- c(participant, pp)
          rotated     <- c(rotated,     rot)
          passive     <- c(passive,     pas)
          variance    <- c(variance,    pcVar)
            
        }
        
      }
        
    }

    groupVar <- data.frame(participant, rotated, passive, variance)
    
    # WRITE CSVs!
    write.csv(groupVar,file=sprintf('data/%s_localization_var.csv',group),row.names=FALSE,quote=FALSE)
    
  }
  
}

# Plots -----

plotLocalizationVariance <- function(target='inline') {
  
  
  styles <- getStyle()
  
  if (target == 'svg') {
    svglite::svglite(file='doc/Fig6.svg', width=8, height=3, system_fonts=list(sans='Arial'))
  }
  
  par(mar=c(4,4,2,0.1))

  layout(matrix(c(1,2,3), nrow=1, ncol=3, byrow = TRUE), widths=c(3,2,2), heights=c(1,1))
  
  # get the data:
  
  SlocVar <- read.csv('data/sEDS_localization_var.csv', stringsAsFactors = FALSE)
  SlocVar$group <- 'sEDS'
  ZlocVar <- read.csv('data/zEDS_localization_var.csv', stringsAsFactors = FALSE)
  ZlocVar$group <- 'zEDS'
  
  locVar <- rbind(SlocVar, ZlocVar)
  locVar$std <- sqrt(locVar$var)
  
  # get the y-axis equal for all panels:
  
  ylims <- c(0,16)
  
  # # # # # # # # # #
  # panel A: localization variance as in ANOVA
  
  plot(c(0.5,2.5),c(0,0),col=rgb(0.5,0.5,0.5),type='l',lty=2,xlim=c(0.5,2.5),ylim=ylims,xlab='condition',ylab='standard deviation',xaxt='n',yaxt='n',bty='n',main='',font.main=1)
  
  mtext('A', side=3, outer=TRUE, at=c(0,1), line=-1, adj=0, padj=1)
  
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
  axis(side=2,at=seq(0,16,4))
  
  
  # # # # # # # # # #
  # panel B: aligned localization variance descriptives
  
  
  plot(c(0,5),c(0,0),col=rgb(0.5,0.5,0.5),type='l',lty=2,xlim=c(0.5,4.5),ylim=ylims,xlab='condition',ylab='',xaxt='n',yaxt='n',bty='n',main='',font.main=1)
  
  mtext('B', side=3, outer=TRUE, at=c(3/7,1), line=-1, adj=0, padj=1)
  
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
      
      points(x=X,y=Y,pch=16,cex=1.5,col=as.character(styles$color_trans[groupno]))
      
      meandist <- getConfidenceInterval(data=locvars, method='bootstrap', resamples=5000, FUN=mean, returndist=TRUE)
       
      DX <- meandist$density$x
      DY <- meandist$density$y / max(meandist$density$y) / 3
       
      DX <- c(DX[1], DX, DX[length(DX)])
      DY <- c(0,     DY, 0)
      
      Xoffset <- (passive*2)+(groupno)
       
      polygon(x=DY+Xoffset, y=DX, border=FALSE, col=as.character(styles$color_trans[groupno]))
       
      lines(x=rep(Xoffset,2),y=meandist$CI95,col=as.character(styles$color_solid[groupno]))
      
      points(x=Xoffset,y=mean(locvars),pch=16,cex=1.5,col=as.character(styles$color_solid[groupno]))
      
    }
    
  }
  
  ####
  
  axis(side=1, at=c(1.5,3.5), labels=c('active','passive'))

  # # # # # # # # # #
  # panel C: rotated localization variance descriptives
  
  
  plot(c(0,5),c(0,0),col=rgb(0.5,0.5,0.5),type='l',lty=2,xlim=c(0.5,4.5),ylim=ylims,xlab='condition',ylab='',xaxt='n',yaxt='n',bty='n',main='',font.main=1)
  
  mtext('C', side=3, outer=TRUE, at=c(5/7,1), line=-1, adj=0, padj=1)
  
  
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
      
      points(x=X,y=Y,pch=16,cex=1.5,col=as.character(styles$color_trans[groupno]))
      
      meandist <- getConfidenceInterval(data=locvars, method='bootstrap', resamples=5000, FUN=mean, returndist=TRUE)
      
      DX <- meandist$density$x
      DY <- meandist$density$y / max(meandist$density$y) / 3
      
      DX <- c(DX[1], DX, DX[length(DX)])
      DY <- c(0,     DY, 0)
      
      Xoffset <- (passive*2)+(groupno)
      
      polygon(x=DY+Xoffset, y=DX, border=FALSE, col=as.character(styles$color_trans[groupno]))
      
      lines(x=rep(Xoffset,2),y=meandist$CI95,col=as.character(styles$color_solid[groupno]))
      
      points(x=Xoffset,y=mean(locvars),pch=16,cex=1.5,col=as.character(styles$color_solid[groupno]))
      
    }
    
  }
  
  axis(side=1, at=c(1.5,3.5), labels=c('active','passive'))
  
    
}

# Statistics -----

localizationSTD.ANOVA <- function() {

  SlocVar <- read.csv('data/sEDS_localization_var.csv', stringsAsFactors = FALSE)
  SlocVar$group <- 'control'
  ZlocVar <- read.csv('data/zEDS_localization_var.csv', stringsAsFactors = FALSE)
  ZlocVar$group <- 'EDS'
  
  locVar <- rbind(SlocVar, ZlocVar)
  locVar$std <- sqrt(locVar$var)
  
  locVar$rotated <- as.factor(locVar$rotated)
  locVar$passive <- as.factor(locVar$passive)
  
  print(ez::ezANOVA(data=locVar,dv=std,wid=participant,between=group,within=c(rotated,passive),type=3))

}

localicalizationSTDt.tests <- function() {
  
  groups <- c('sEDS', 'zEDS')
  
  locVar <- list()
  
  for (group in groups) {
    groupLocVar <- read.csv(sprintf('data/%s_localization_var.csv', group), stringsAsFactors = FALSE)
    groupLocVar$std <- sqrt(groupLocVar$variance)
    locVar[[group]] <- groupLocVar
  }
  
  for (rotated in c(0,1)) {
    for (passive in c(0,1)) {
      stdS <- locVar[['sEDS']]$std[which(locVar[['sEDS']]$rotated == rotated & locVar[['sEDS']]$passive == passive)]
      stdZ <- locVar[['zEDS']]$std[which(locVar[['zEDS']]$rotated == rotated & locVar[['zEDS']]$passive == passive)]
      cat(sprintf('\n%s %s\n\n',c('ALIGNED','ROTATED')[rotated+1],c('ACTIVE','PASSIVE')[passive+1] ))
      print(t.test(stdS,stdZ,alternative='less'))
      cat(sprintf('eta-squared: %0.5f\n\n', etaSquaredTtest(stdS, stdZ)))
    }
  }
  
}

localizationSTDlogreg <- function() {
  
  SlocVar <- read.csv('data/sEDS_localization_var.csv', stringsAsFactors = FALSE)
  SlocVar$group <- 'control'
  ZlocVar <- read.csv('data/zEDS_localization_var.csv', stringsAsFactors = FALSE)
  ZlocVar$group <- 'EDS'
  
  locVar <- rbind(SlocVar, ZlocVar)
  locVar$std <- sqrt(locVar$var)
  
  locVar <- locVar[which(locVar$rotated == 0 & locVar$passive == 1),]
  
  locVar$group <- as.factor(locVar$group)
  
  mylogit <- glm(group ~ std, data = locVar, family = "binomial")
  summary(mylogit)
}