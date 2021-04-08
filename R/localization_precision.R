
# Data Handling -----

getLocalizationTrends <- function(groups=c('sEDS','zEDS')) {
  
  for (group in groups) {
    
    # set up vectors that will become data frame columns:
    participant <- c()
    rotated <- c()
    passive <- c()
    handlocation_deg <- c()
    localizationbias_deg <- c()
    
    # hand locations where we estimate bias:
    X <- seq(35,145)
    
    # read in the screened and compiled data:
    loc <- read.csv(sprintf('data/%s_localization.csv',group), stringsAsFactors = FALSE)
    loc$tap_error <- loc$tap_angle - loc$hand_angle
    loc$predicted_error <- NA
    
    # loop through participants:
    participants <- unique(loc$participant)
    
    for (pp in participants) {
      
      # look at rotated/aligned and passive/active seprately:
      for (rot in c(0,1)) {
        
        for (pas in c(0,1)) {
          
          # index for all relevant rows:
          idx <- which(loc$participant == pp & loc$rotated == rot & loc$passive == pas)
          
          # get relevant data:
          subdf <- loc[idx,]
          hand_angle <- subdf$hand_angle
          tap_error  <- subdf$tap_error

          # spline interpolate through the data:
          spl <- smooth.spline(x=hand_angle, y=tap_error, spar=0.90, keep.data=F )
          splTrend <- predict(spl, x=X)$y
          
          # make sure we don't extrapolate:
          nanidx <- c(which(X < min(hand_angle)), which(X > max(hand_angle)))
          splTrend[nanidx] <- NA
          
          # add data to vectors:
          participant <- c(participant, rep(pp, length(X)))
          rotated <- c(rotated, rep(rot, length(X)))
          passive <- c(passive, rep(pas, length(X)))
          handlocation_deg <- c(handlocation_deg, X)
          localizationbias_deg <- c(localizationbias_deg, splTrend)
          
        }
        
      }
      
    }
    
    # end of group loop
    groupLocBias <- data.frame(participant, rotated, passive, handlocation_deg, localizationbias_deg)
    
    write.csv(groupLocBias, file=sprintf('data/%s_localizationbias.csv',group), row.names=FALSE, quote=FALSE)

  }
  
  # return nothing: data is stored in files (or downloaded from OSF)
  
}

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
  
  # detrend the localization responses, so we get rid of biases
  # and only have variance around this bias:
  DetLoc <- detrendLocalization()
  
  allVar <- NA
  
  for (group in names(DetLoc)) {
    
    loc <- DetLoc[[group]]
    
    # get variance split up by condition:
    loc$variance <- loc$tap_error - loc$predicted_error
    groupVar <- aggregate(variance ~ rotated + passive + participant, data = loc, FUN=var)
    
    # write this to csv:
    write.csv(groupVar,file=sprintf('data/%s_localization_var.csv',group),row.names=FALSE,quote=FALSE)
    
    # also get generic variance per participant:
    if (is.data.frame(allVar)) {
      allVar <- rbind(allVar, aggregate(variance ~ participant, data = loc, FUN=var))
    } else {
      allVar <- aggregate(variance ~ participant, data = loc, FUN=var)
    }
    
  }
  
  participants <- read.csv('data/participants.csv', stringsAsFactors = FALSE)
  participants <- participants[order(participants$ID),]
  
  allVar <- cbind(allVar[order(allVar$participant),], participants[,c('Beighton','folder','age','sex','type')])
  
  write.csv(allVar,file=sprintf('data/all_localization_var.csv',group),row.names=FALSE,quote=FALSE)
  
}




# Plots -----

# plotLocalizationVariance <- function(target='inline') {
#   
#   
#   styles <- getStyle()
#   
#   if (target == 'svg') {
#     svglite::svglite(file='doc/Fig6.svg', width=8, height=3, system_fonts=list(sans='Arial'))
#   }
#   
#   par(mar=c(4,4,2,0.1))
# 
#   layout(matrix(c(1,2,3), nrow=1, ncol=3, byrow = TRUE), widths=c(3,2,2), heights=c(1,1))
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
#   ylims <- c(0,16)
#   
#   # # # # # # # # # #
#   # panel A: localization variance as in ANOVA
#   
#   plot(c(0.5,2.5),c(0,0),col=rgb(0.5,0.5,0.5),type='l',lty=2,xlim=c(0.5,2.5),ylim=ylims,xlab='condition',ylab='standard deviation',xaxt='n',yaxt='n',bty='n',main='',font.main=1)
#   
#   mtext('A', side=3, outer=TRUE, at=c(0,1), line=-1, adj=0, padj=1)
#   
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
#   axis(side=1,at=c(1,2),labels=c('aligned','rotated'))
#   axis(side=2,at=seq(0,16,4))
#   
#   
#   # # # # # # # # # #
#   # panel B: aligned localization variance descriptives
#   
#   
#   plot(c(0,5),c(0,0),col=rgb(0.5,0.5,0.5),type='l',lty=2,xlim=c(0.5,4.5),ylim=ylims,xlab='condition',ylab='',xaxt='n',yaxt='n',bty='n',main='',font.main=1)
#   
#   mtext('B', side=3, outer=TRUE, at=c(3/7,1), line=-1, adj=0, padj=1)
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
#       points(x=X,y=Y,pch=16,cex=1.5,col=as.character(styles$color_trans[groupno]))
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
#       points(x=Xoffset,y=mean(locvars),pch=16,cex=1.5,col=as.character(styles$color_solid[groupno]))
#       
#     }
#     
#   }
#   
#   # # # # 
#   
#   axis(side=1, at=c(1.5,3.5), labels=c('active','passive'))
# 
#   # # # # # # # # # #
#   # panel C: rotated localization variance descriptives
#   
#   
#   plot(c(0,5),c(0,0),col=rgb(0.5,0.5,0.5),type='l',lty=2,xlim=c(0.5,4.5),ylim=ylims,xlab='condition',ylab='',xaxt='n',yaxt='n',bty='n',main='',font.main=1)
#   
#   mtext('C', side=3, outer=TRUE, at=c(5/7,1), line=-1, adj=0, padj=1)
#   
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
#       points(x=X,y=Y,pch=16,cex=1.5,col=as.character(styles$color_trans[groupno]))
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
#       points(x=Xoffset,y=mean(locvars),pch=16,cex=1.5,col=as.character(styles$color_solid[groupno]))
#       
#     }
#     
#   }
#   
#   axis(side=1, at=c(1.5,3.5), labels=c('active','passive'))
#   
#   if (target == 'svg') {
#     dev.off()
#   }
#   
#     
# }

plotBeightonLocSTD <- function(target='inline') {
  
  styles <- getStyle()
  
  locSTD <- read.csv('data/all_localization_var.csv', stringsAsFactors = FALSE)
  locSTD$std <- sqrt(locSTD$variance)
  
  fw<-4
  fh<-4
  
  if (target == 'svg') {
    svglite::svglite(file='doc/Fig7.svg', width=fw, height=fh, scaling=1, fix_text_size = FALSE)
  }
  if (target == 'tiff') {
    tiff(filename='doc/Fig7.tiff',width=fw*1200,height=fh*1200,units='px',type='cairo',compression='lzw',res=1200)
  }
  if (target == 'pdf') {
    pdf(file='doc/Fig7.pdf',width=fw,height=fh)
  }
  
  
  
  plot(-1000,-1000,
       main='localization precision and hypermobility',
       xlab='Beighton score',ylab='localization SD [°]',
       xlim=c(-1,10),ylim=c(0,10),
       bty='n',ax=F,font.main=1)
  
  myLinReg <- lm(std ~ Beighton, data=locSTD)
  
  Xh <- seq(0,9,.01)
  LRfit <- predict(myLinReg, newdata=data.frame(Beighton=Xh), interval='confidence')
  
  polygon(c(Xh,rev(Xh)), c(LRfit[,'lwr'], rev(LRfit[,'upr'])), col='#E9E9E9', border=NA)
  
  lines(x=Xh,y=LRfit[,'fit'], col='#000000')
  
  for (group_no in c(1:length(styles$group))) {
    
    group <- styles$group[group_no]
    
    subdf <- locSTD[which(locSTD$folder == group),]
    
    col <- as.character(styles$color_solid[group_no])
    
    points(subdf$Beighton, subdf$std, col=col, cex=1, pch=1)
    
  }
  
  labels <- c(as.character(styles$label))
  colors <- c(as.character(styles$color_solid))
  legend(0, 3, labels, col=colors, bty='n', cex=1, pch=c(1,1))
  
  labels <- c('regression')
  colors <- c('#000000')
  legend(4.5, 3, labels, col=colors, bty='n', cex=1, lw=c(1), seg.len = c(1), pch=c(NA))
  
  axis(1,seq(0,9,3))
  axis(2,seq(0,10,5))
  
  if (target %in% c('tiff','svg','pdf')) {
    dev.off()
  }
  
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


correlateBeightonLocSTD <- function() {
  
  locSTD <- read.csv('data/all_localization_var.csv', stringsAsFactors = FALSE)
  locSTD$std <- sqrt(locSTD$variance)
    
  #print(cor.test(locSTD$std, locSTD$Beighton))
  
  print(summary(lm(std ~ Beighton, data=locSTD)))
  
  
}

predictGroupBYlocvar <- function() {
  
  locSTD <- read.csv('data/all_localization_var.csv', stringsAsFactors = FALSE)
  
  locSTD$std <- sqrt(locSTD$variance)
  
  locSTD$group <- locSTD$type
  locSTD$group[which(locSTD$group != 'control')] <- 'EDS'
  locSTD$group <- as.factor(locSTD$group)
  
  mylogit <- glm(group ~ std, data = locSTD, family = "binomial")
  summary(mylogit)
  
}

localizationSTDlogreg <- function() {
  
  SlocVar <- read.csv('data/sEDS_localization_var.csv', stringsAsFactors = FALSE)
  SlocVar$group <- 'control'
  ZlocVar <- read.csv('data/zEDS_localization_var.csv', stringsAsFactors = FALSE)
  ZlocVar$group <- 'EDS'
  
  locVar <- rbind(SlocVar, ZlocVar)
  locVar$std <- sqrt(locVar$var)
  
  locVar$group <- as.factor(locVar$group)
  
  for (rotated in c(0,1)) {
    
    for (passive in c(0,1)) {
      
      mylogit <- glm(group ~ std, data = locVar[which(locVar$rotated == rotated & locVar$passive == passive),], family = "binomial")
      
      cat(sprintf('Predict group from localization variance in %s, %s localization:\n', c('ALIGNED','ROTATED')[rotated+1], c('ACTIVE','PASSIVE')[passive+1]))
      
      print(summary(mylogit))
      
    }
    
  }
  
}
