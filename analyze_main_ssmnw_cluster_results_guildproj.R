##############################################
###### ANALYSE GUILDPROJ CLUSTER RESULTS ######
##############################################
#
#   @DESCRIPTION
#   Analysis code.
#   Plot and maybe statistically analyze the
#   stability and network structure results 
#   for the guild sub-project of SSMNW.
#   Written for main_ssmnw results run and collated
#   on the NSC cluster in Linkoping.
#
#   @AUTHOR
#   Alva Curtsdotter, postdoc @ Emory University, 
#   2018-09-26
#   
##############################################

#### Setup ####
{
  
  # clearing the workspace. 
  rm(list=ls(all=TRUE))
  
  # Read in libraries:
  # tidyverse reads all packages in the tidyverse, in which ggplot and dplyr are found
  require(tidyverse)
  require(gridExtra)
  # Binomial CIs
  require(binom)
  
  # Data paths and such
  timestamp   <-  '2018Oct26_2123'
  data.dir    <- paste('C:/Users/icurtsd/Documents/Project Atlanta/Project Structure and stability of mutualistic networks/data/clusterdata/analysis/', timestamp, '/', sep='')
  data.dir.2  <- 'C:/Users/icurtsd/Documents/Project Atlanta/Project Structure and stability of mutualistic networks/data/clusterdata/analysis/bla/jacobianConnectance/'
  data.file.1 <- 'collatedStabResults.csv'
  data.file.2 <- 'stabilitySummary.csv'

  save.dir <- data.dir
  if ( !(dir.exists(save.dir)) ){ 
    dir.create(save.dir)
  } 
  
  # Some variable definition
  design <- 'substitutive'
}


### Load and prep data ###
{
  # Load data frame with replicate level data of stability and network structure
  setwd(data.dir)
  rep.frame <- read.csv(data.file.1, sep=' ')
  str(rep.frame)
  
  # Load data frame with scenario level data of stability
  setwd(data.dir)
  scen.frame <- read.csv(data.file.2, sep=' ')
  scen.frame <- scen.frame[rowSums(is.na(scen.frame)) != ncol(scen.frame), ]
  str(scen.frame)
  
  # scen.frame$rc <- NA
  # wgcIgnoreLevels <- sort(unique(scen.frame$wgcIgnore))
  # wgcStrings <- c('wgC=0', 'wgC=btwgC')
  # weightLevels <- unique(scen.frame$intWeightStd)
  # for (j in weightLevels){
  #   for (i in 1:length(wgcIgnoreLevels)){
  #     scen.frame$rc[scen.frame$wgcIgnore==wgcIgnoreLevels[i] & scen.frame$intWeightStd==j] <- paste(wgcStrings[i], ', IS=', round(j,digits=2) , sep='')
  #   }
  # }
  # # scen.frame$rc <- as.factor(scen.frame$rc)
  # scen.frame$rc = factor(scen.frame$rc,levels(factor(scen.frame$rc))[c(1,4,2,5,3,6)])
  # # levels(scen.frame$rc) <-  str(as.character(levels(scen.frame$rc)[c(1,4,2,5,3,6)]))
  
  scen.frame$colCode <- NA
  scen.frame$colCode[scen.frame$wgcIgnore==0 & scen.frame$guildStruct==0] <- 'pure rand'
  scen.frame$colCode[scen.frame$wgcIgnore==0 & scen.frame$guildStruct==1] <- 'pure guild'
  scen.frame$colCode[scen.frame$wgcIgnore==1 & scen.frame$guildStruct==0] <- 'mix rand'
  scen.frame$colCode[scen.frame$wgcIgnore==1 & scen.frame$guildStruct==1] <- 'mix guild'
  scen.frame$colCode <- factor(scen.frame$colCode)
  levels(scen.frame$colCode)
  scen.frame$colCode = factor(scen.frame$colCode,levels(factor(scen.frame$colCode))[c(4,3,2,1)])
  str(scen.frame)
    
  scen.frame$btwGuildC <- round(scen.frame$btwGuildC, digits = 2)
  
  # scen.frame$guildStruct <- as.factor(scen.frame$guildStruct)
  
  scen.frame$locStabCount <- scen.frame$localStability/100*unique(scen.frame$trials, na.rm=T)
  
  bob <- binom.confint(as.integer(scen.frame$locStabCount), rep(max(scen.frame$trials), nrow(scen.frame)), methods = "agresti-coull")
  scen.frame <- cbind(scen.frame, bob[,4:6])
  rm(bob) #cleanup
  scen.frame$lower[scen.frame$lower < 0] <- 0
  scen.frame$upper[scen.frame$upper > 1] <- 1
  str(scen.frame)
  
  scen.frame$jacC <- scen.frame$btwGuildC/2 #true for substitutive case and no wGuildC case)
  scen.frame$complexity <- scen.frame$intWeightStd*sqrt(scen.frame$S*scen.frame$jacC)
}


### Load and check some stuff ###
{
  setwd(data.dir.2)
  
  allVerdicts <- as.data.frame(matrix(NA, nrow=nrow(scen.frame), ncol=2))
  colnames(allVerdicts) <- c('scenID', 'verdict')
  for (scid in scen.frame$scenID){
    load(paste('jacobianCheck_scenario_', scid, '.RData', sep=''))
    allVerdicts$scenID[scid]  <- scid
    allVerdicts$verdict[scid] <- verdict
  }
  if ( any(allVerdicts$verdict==0) ){
    print('Buggy scenarios found')
    buggy <- which(allVerdicts$verdict==0)
    print(buggy)
  }else{
    print('All scenarios seem OK.')
  }
  
  allCvals <- as.data.frame(matrix(NA, nrow=nrow(scen.frame), ncol=2))
  colnames(allVerdicts) <- c('scenID', 'verdict')
  
}

### Plot stuff ###
{
  # ggplot: Stability bars (based on scen.frame)
  {
    # Plot options
    scaleOp =  'fixed'# 'free'# 
    facetRow <- "guildStruct"
    facetCol <- "wgcIgnore"
    
    # initiate plot
    setwd(save.dir)
    # windows(width=8.3, height=11.7)
    pdf(file=paste('LocalStabilityBars_onintWeightStd_GuildProj_ggplot_', scaleOp, 'Y.pdf', sep=''), width=8.3, height=11.7)
    
    for ( i in sort(unique(scen.frame$intTypes)) ){
      for ( j in sort(unique(scen.frame$S[scen.frame$intTypes==i])) ){
        
        data <- scen.frame[scen.frame$intTypes==i & scen.frame$S==j ,]
        
        p1 <- ggplot(data, aes(x=intWeightStd, y = localStability)) + geom_bar(stat="identity") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          facet_wrap(c(facetRow, facetCol), nrow = length(unique(eval(parse(text=paste('data$',facetRow,sep=''))))), 
                     ncol = length(unique(eval(parse(text=paste('data$',facetCol,sep=''))))), scales = scaleOp,
                     shrink = TRUE, labeller = "label_both", as.table = TRUE,
                     switch = NULL, drop = TRUE, dir = "h", strip.position = "top")
        
        grid.arrange(
          p1,
          nrow = 1,
          top = paste( i , 'with', j, 'species, req jC =', scen.frame$jacC, sep=' '),#paste( i , 'with', j, 'species, ', design, ' design' , sep=' '),#paste( i, 'species', wGuildCCode$name[wGuildCCode$value==j], 'within guild competition', sep=' '),
          newpage = T
        )
        
        
      } # Loop over j
    } # Loop over i
    dev.off()
  }
  
  # ggplot: Stability bars (based on rep.frame) # Was just for silly debugging
  {
    # # Plot options
    # scaleOp =  'fixed'# 'free'# 
    # 
    # # initiate plot
    # setwd(save.dir)
    # # windows(width=8.3, height=11.7)
    # pdf(file=paste('LocalStabilityBars_NestProj_repFrame_ggplot_', scaleOp, 'Y.pdf', sep=''), width=8.3, height=11.7)
    # 
    # for ( i in sort(unique(rep.frame$intTypes)) ){
    #   for ( j in sort(unique(rep.frame$S[rep.frame$intTypes==i])) ){
    #     
    #     data <- rep.frame[rep.frame$intTypes==i & rep.frame$S==j ,]
    #     
    #     p1 <- ggplot(data, aes(x=nestCode, y = localStability)) + geom_col() +
    #       theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    #       facet_wrap(c("wGcode", "btwGuildC"), nrow = length(unique(data$wGcode)), ncol = length(unique(data$btwGuildC)), scales = scaleOp,
    #                  shrink = TRUE, labeller = "label_both", as.table = TRUE,
    #                  switch = NULL, drop = TRUE, dir = "h", strip.position = "top")
    #     
    #     grid.arrange(
    #       p1,
    #       nrow = 1,
    #       top = paste( i , 'with', j, 'species, ', design, ' design' , sep=' '),#paste( i, 'species', wGuildCCode$name[wGuildCCode$value==j], 'within guild competition', sep=' '),
    #       newpage = T
    #     )
    #     
    #     
    #   } # Loop over j
    # } # Loop over i
    # dev.off()
  }
  
  # ggplot: Stability points with error bars (based on scen.frame)
  {
    # Plot options
    scaleOp =  'fixed'# 'free'# 
    limits.all <- aes(ymax = upper, ymin = lower)
    
    # initiate plot
    setwd(save.dir)
    # windows(width=8.3, height=11.7)
    pdf(file=paste('LocalStabilityPointsAndErrorbars_GuildProj_ggplot_', scaleOp, 'Y.pdf', sep=''), width=8.3, height=11.7)
    
    for ( i in sort(unique(scen.frame$intTypes)) ){
      for ( j in sort(unique(scen.frame$S[scen.frame$intTypes==i])) ){
        
        data <- scen.frame[scen.frame$intTypes==i & scen.frame$S==j ,]
        
        p1 <- ggplot(data, aes(x=guildStruct, y = mean)) + 
          geom_point(size = 2) +
          geom_errorbar(limits.all, width=0.2, linetype = 1) + 
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          facet_wrap(c("rc", "btwGuildC"), nrow = length(unique(data$rc)), ncol = length(unique(data$btwGuildC)), scales = scaleOp,
                     shrink = TRUE, labeller = "label_both", as.table = TRUE,
                     switch = NULL, drop = TRUE, dir = "h", strip.position = "top")
        
        grid.arrange(
          p1,
          nrow = 1,
          top = paste( i , 'with', j, 'species, ', design, ' design' , sep=' '),#paste( i, 'species', wGuildCCode$name[wGuildCCode$value==j], 'within guild competition', sep=' '),
          newpage = T
        )
        
        
      } # Loop over j
    } # Loop over i
    dev.off()
  }
  
  # ggplot: Stability vs sigma*sqrt(SC)
  {
    # Plot options
    scaleOp =  'free'# 'fixed'# 
    limits.all <- aes(ymax = upper, ymin = lower)
    
    facetRow <- "guildStruct"
    facetCol <- "wgcIgnore"
    
    # initiate plot
    setwd(save.dir)
    windows(width=8.3, height=11.7)
    # pdf(file=paste('LocalStability_vs_Complexity_GuildProj_ggplot_', scaleOp, 'Y.pdf', sep=''), width=8.3, height=11.7)
    
    
    for ( i in sort(unique(scen.frame$intTypes)) ){
      
      
      data <- scen.frame[scen.frame$intTypes==i ,]
      
      p1 <- ggplot(data, aes(x=complexity, y = mean)) + 
        geom_point(size = 2) +
        geom_errorbar(limits.all, width=0.2, linetype = 1) + 
        geom_line() + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        facet_wrap(c(facetRow, facetCol), nrow = length(unique(paste('data$',facetRow,sep=''))), 
                   ncol = length(unique(paste('data$',facetCol,sep=''))), scales = scaleOp,
                   shrink = TRUE, labeller = "label_both", as.table = TRUE,
                   switch = NULL, drop = TRUE, dir = "h", strip.position = "top")
      
      grid.arrange(
        p1,
        nrow = 1,
        top = paste( i , 'with', unique(scen.frame$S), 'species, and jacC=', unique(scen.frame$jacC), ', ', design, ' design' , sep=' '),#paste( i, 'species', wGuildCCode$name[wGuildCCode$value==j], 'within guild competition', sep=' '),
        newpage = T
      )
      
    } # Loop over i
    dev.off()
  }
  
  # base plot: Stability vs sigma*sqrt(SC)
  # Link to add logistic regression
  {
    # Plot options
    scaleOp =  'free'# 'fixed'# 
    limits.all <- aes(ymax = upper, ymin = lower)
    
     # initiate plot
    setwd(save.dir)
    # windows(width=8.3, height=11.7)
    pdf(file=paste('LocalStability_vs_Complexity_GuildProj_baseplot_', scaleOp, 'Y.pdf', sep=''), width=8.3, height=11.7)
    par(mfrow=(c(length(unique(scen.frame$intTypeSymmetry)), length(unique(scen.frame$colCode)))))
    
    for ( i in sort(unique(scen.frame$intTypes)) ){
      data <- scen.frame[scen.frame$intTypes==i ,]
      for ( j in sort(unique(data$intTypeSymmetry)) ){
        for ( k in sort(unique(data$colCode)) ){
          print(j)
          plotdata <- data[data$colCode==k & data$intTypeSymmetry==j ,]
          plotdata <- plotdata[order(plotdata$intWeightStd),];
          rownames(plotdata) <- seq(1,nrow(plotdata),1)
          
          plot(plotdata$complexity, plotdata$mean)
          abline(h=0.5, col='red')
          title(main=paste(substr(i,1,4), k, ', cor:', j, sep=' '), cex.main=0.9)
          
          closest  <- plotdata[abs(plotdata$mean - 0.5)==min(abs(plotdata$mean - 0.5)),]
          if (nrow(closest)==1){
            rows     <- c(as.numeric(rownames(closest))-1, as.numeric(rownames(closest)), as.numeric(rownames(closest))+1)
            funcdata <- plotdata[rows,]
            
            fit     <- lm(data=funcdata, mean ~ complexity)
            newx    <- data.frame(complexity=seq(min(funcdata$complexity),max(funcdata$complexity), 0.0001))
            fitline <- predict(fit,newdata=newx)
            est     <- data.frame(newx, fitline)
            cross   <- est[which.min(abs(est$fitline - 0.5)),]
            
            lines(est, col='blue', lwd=2)
            abline(v=cross$complexity, col='green')
            text(x=1.15*cross$complexity, y=0.55, labels=paste(round(cross$complexity, digits=2) ))
          }
        }
      }
    }
    dev.off()
  }
}

