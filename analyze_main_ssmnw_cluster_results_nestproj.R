##############################################
###### ANALYSE NESTPROJ CLUSTER RESULTS ######
##############################################
#
#   @DESCRIPTION
#   Analysis code.
#   Plot and maybe statistically analyze the
#   stability and network structure results 
#   for the nestedness sub-project of SSMNW.
#   Written for main_ssmnw results run and collated
#   on the NSC cluster in Linkoping.
#
#   @AUTHOR
#   Alva Curtsdotter, postdoc @ Emory University, 
#   2018-08-20
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
  data.dir <- 'C:/Users/icurtsd/Documents/Project Atlanta/Project Structure and stability of mutualistic networks/data/clusterdata/analysis/2018Jul31_2035/'
  
  filename1 <- 'collatedStabAndNetworkResults.csv'
  filename2 <- 'collatedRowandColSumMetrics.csv'
  filename3 <- 'collatedJacobianCheckMetrics.csv'
  filename4 <- 'stabilitySummary.csv'
  
  save.dir <- paste(data.dir, 'new181102/', sep='')
  if ( !(dir.exists(save.dir)) ){ 
    dir.create(save.dir)
  } 
  
  # Some variable definition
  checkMetrics  <- c('Cj.tot', 'Cj.off', 'btwgC', 'wgC', 'corrcoef' ) 
  
  structMetrics <- c('bNestNODF', 'bNestTemp', 'bNestSpecRad', 'qNestNODF', 'qNestWine', 'qNestSpecRad',
                      'bFcomp', 'qFcomp', 'bFcomp2', 'qFcomp2', 'fc.hl', 'fc.ll', 'no.hl', 'no.ll')
  # structMetrics  <- c('bNestNODF', 'bNestTemp', 'bNestSpecRad', 'qNestNODF', 'qNestWine', 'qNestSpecRad', 'bFcomp', 'qFcomp')
  
  rcSumMetrics  <- c( 'max.rsJ','max.csJ','max.rsBbip','max.csBbip','max.rsQbip','max.csQbip', 
                       'mean.rsJ','mean.csJ','mean.rsBbip','mean.csBbip','mean.rsQbip','mean.csQbip',
                       'std.rsJ','std.csJ','std.rsBbip','std.csBbip','std.rsQbip','std.csQbip',
                       'lm.rsJ','lm.csJ','lm.rsBbip','lm.csBbip','lm.rsQbip','lm.csQbip',
                       'PJ.rsJ','PJ.csJ', 'PJ.tsJ','PJ.rsBbip','PJ.csBbip', 'PJ.tsBbip', 'PJ.rsQbip','PJ.csQbip', 'PJ.tsQbip')
  
  qNestMetrics  <- c('qNestNODF', 'qNestWine', 'qNestSpecRad')
  
  qCompMetrics  <- c('qFcomp', 'qFcomp2', 'fc.hl', 'fc.ll', 'no.hl', 'no.ll')
}

### Load and prep data ###
{
  # Load data frame with replicate level data of stability and network structure
  setwd(data.dir)
  rep.frame <- read.csv(filename1, sep=' ')
  str(rep.frame)
  
  qFcomp.factorial <- as.factor(ceiling(rep.frame$qFcomp))
  rep.frame <- cbind(rep.frame, qFcomp.factorial)
  rep.frame$qFcomp.factorial <- factor(qFcomp.factorial, levels = c(as.character(sort(unique(ceiling(rep.frame$qFcomp))))))

  rep.frame$binNestLevel.factorial   <- as.factor(rep.frame$binNestLevel)
  rep.frame$quantNestLevel.factorial <- as.factor(rep.frame$quantNestLevel)

  rep.frame$wGcode = NA
  rep.frame$wGcode[rep.frame$wGuildC == 0 ]                                     = 'no'
  rep.frame$wGcode[rep.frame$wGuildC == 1 & rep.frame$wGuildWeightStd == 0.25]  = 'intraGuildC=1 and weak'
  rep.frame$wGcode[rep.frame$wGuildC == 1 & rep.frame$wGuildWeightStd == 0.50]  = 'intraGuildC=1 and strong'
  rep.frame$wGcode[rep.frame$wgcIgnore == 1 & rep.frame$wGuildWeightStd == 0.25]  = 'intraGuildC=btwGuildC and weak'
  rep.frame$wGcode[rep.frame$wgcIgnore == 1 & rep.frame$wGuildWeightStd == 0.50]  = 'intraGuildC=btwGuildC and strong'
  rep.frame$wGcode <- factor(rep.frame$wGcode, levels = c('no', 'intraGuildC=btwGuildC and weak', 'intraGuildC=1 and weak', 'intraGuildC=btwGuildC and strong', 'intraGuildC=1 and strong'))

  rep.frame$nestCode <- factor(rep.frame$nestCode, levels = c('bNqN', 'bNqR', 'bNqA', 'bRqN', 'bRqR', 'bRqA','bAqN', 'bAqR', 'bAqA'))
  
  str(rep.frame)
  
  # Load some more network metrics if present
  
  # The row and col sum metrics
  setwd(data.dir)
  file2Present <- length(dir(pattern=filename2))>0
  if ( file2Present ){
    rc.frame <- read.csv(filename2, sep=' ')
    str(rc.frame)
  } 
  
  # The actual jacobian connectance and pairwise IS correlations
  setwd(data.dir)
  file3Present <- length(dir(pattern=filename3))>0
  if ( file3Present ){
    c.frame <- read.csv(filename3, sep=' ')
    str(c.frame)
  } 
  
  # Load data frame with scenario level data of stability
  setwd(data.dir)
  scen.frame <- read.csv(filename4, sep=' ')
  str(scen.frame)
  
  scen.frame$binNestLevel.factorial   <- as.factor(scen.frame$binNestLevel)
  scen.frame$quantNestLevel.factorial <- as.factor(scen.frame$quantNestLevel)
  
  scen.frame$locStabCount <- scen.frame$localStability/100*unique(scen.frame$trials)
  
  bob <- binom.confint(as.integer(scen.frame$locStabCount), rep(max(scen.frame$trials), nrow(scen.frame)), methods = "agresti-coull")
  scen.frame <- cbind(scen.frame, bob[,4:6])
  rm(bob) #cleanup
  scen.frame$lower[scen.frame$lower < 0] <- 0
  scen.frame$upper[scen.frame$upper > 1] <- 1
  
  scen.frame$wGcode <- factor(scen.frame$wGcode, levels = c('no', 'intraGuildC=btwGuildC and weak', 'intraGuildC=1 and weak', 'intraGuildC=btwGuildC and strong', 'intraGuildC=1 and strong'))
  scen.frame$nestCode <- factor(scen.frame$nestCode, levels = c('bNqN', 'bNqR', 'bNqA', 'bRqN', 'bRqR', 'bRqA','bAqN', 'bAqR', 'bAqA'))
  
  str(scen.frame)
  
}

### Plot stuff ###
{
  # ggplot: Stability bars (based on scen.frame)
  {
    # Plot options
    scaleOp =  'fixed'# 'free'# 
    
    # initiate plot
    setwd(save.dir)
    # windows(width=8.3, height=11.7)
    pdf(file=paste('LocalStabilityBars_NestProj_ggplot_', scaleOp, 'Y.pdf', sep=''), width=8.3, height=11.7)
    
    for ( i in sort(unique(scen.frame$intTypes)) ){
      for ( j in sort(unique(scen.frame$S[scen.frame$intTypes==i])) ){
        
        data <- scen.frame[scen.frame$intTypes==i & scen.frame$S==j ,]
        
        p1 <- ggplot(data, aes(x=nestCode, y = localStability)) + geom_bar(stat="identity") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          facet_wrap(c("wGcode", "btwGuildC"), nrow = length(unique(data$wGcode)), ncol = length(unique(data$btwGuildC)), scales = scaleOp,
                     shrink = TRUE, labeller = "label_both", as.table = TRUE,
                     switch = NULL, drop = TRUE, dir = "h", strip.position = "top")
        
        grid.arrange(
          p1,
          nrow = 1,
          top = paste( i , 'with', j, 'species' , sep=' '),#paste( i, 'species', wGuildCCode$name[wGuildCCode$value==j], 'within guild competition', sep=' '),
          newpage = T
        )
        
        
      } # Loop over j
    } # Loop over i
    dev.off()
  }
  
  # ggplot: Stability bars (based on rep.frame) # Was just for debugging
  {
    # Plot options
    scaleOp =  'fixed'# 'free'# 
    
    # initiate plot
    setwd(save.dir)
    # windows(width=8.3, height=11.7)
    pdf(file=paste('LocalStabilityBars_NestProj_repFrame_ggplot_', scaleOp, 'Y.pdf', sep=''), width=8.3, height=11.7)
    
    for ( i in sort(unique(rep.frame$intTypes)) ){
      for ( j in sort(unique(rep.frame$S[rep.frame$intTypes==i])) ){
        
        data <- rep.frame[rep.frame$intTypes==i & rep.frame$S==j ,]
        
        p1 <- ggplot(data, aes(x=nestCode, y = localStability)) + geom_col() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          facet_wrap(c("wGcode", "btwGuildC"), nrow = length(unique(data$wGcode)), ncol = length(unique(data$btwGuildC)), scales = scaleOp,
                     shrink = TRUE, labeller = "label_both", as.table = TRUE,
                     switch = NULL, drop = TRUE, dir = "h", strip.position = "top")
        
        grid.arrange(
          p1,
          nrow = 1,
          top = paste( i , 'with', j, 'species' , sep=' '),#paste( i, 'species', wGuildCCode$name[wGuildCCode$value==j], 'within guild competition', sep=' '),
          newpage = T
        )
        
        
      } # Loop over j
    } # Loop over i
    dev.off()
  }
  
  # ggplot: Stability points with error bars (based on scen.frame)
  {
    # Plot options
    scaleOp =  'fixed'# 'free'# 
    limits.all <- aes(ymax = upper, ymin = lower)
    
    # initiate plot
    setwd(save.dir)
    # windows(width=8.3, height=11.7)
    pdf(file=paste('LocalStabilityPointsAndErrorbars_NestProj_ggplot_', scaleOp, 'Y.pdf', sep=''), width=8.3, height=11.7)
    
    for ( i in sort(unique(scen.frame$intTypes)) ){
      for ( j in sort(unique(scen.frame$S[scen.frame$intTypes==i])) ){
        
        data <- scen.frame[scen.frame$intTypes==i & scen.frame$S==j ,]
        
        p1 <- ggplot(data, aes(x=nestCode, y = mean, colour = binNestLevel.factorial, shape = quantNestLevel.factorial)) + 
          geom_point(size = 2) +
          geom_errorbar(limits.all, width=0.2, linetype = 1) + 
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          facet_wrap(c("wGcode", "btwGuildC"), nrow = length(unique(data$wGcode)), ncol = length(unique(data$btwGuildC)), scales = scaleOp,
                     shrink = TRUE, labeller = "label_both", as.table = TRUE,
                     switch = NULL, drop = TRUE, dir = "h", strip.position = "top")
        
        grid.arrange(
          p1,
          nrow = 1,
          top = paste( i , 'with', j, 'species' , sep=' '),#paste( i, 'species', wGuildCCode$name[wGuildCCode$value==j], 'within guild competition', sep=' '),
          newpage = T
        )
        
        
      } # Loop over j
    } # Loop over i
    dev.off()
  }
  
  # ggplot: Network metric boxplots (based on rep.frame)
  {
    # Plot options
    scaleOp =  'fixed'# 'free'# 
    
    # initiate plot
    setwd(save.dir)
    pdf(file= paste('NetworkMetricBoxes_NestProj_ggplot_', scaleOp, 'Y.pdf', sep=''), width=8.3, height=11.7)
    
    for ( i in sort(unique(rep.frame$S)) ){
      
      for ( j in 1:length(structMetrics) ){
        
        data <- rep.frame[rep.frame$S==i ,]  
        
        # median.data <- aggregate(data[, structMetrics], list(data$nestCode), median)
        # bnl <- as.factor(c(1,1,1,0,0,0,-1,-1,-1))
        # qnl <- as.factor(c(1,0,-1,1,0,-1,1,0,-1))
        # median.data <- cbind(median.data, bnl, qnl )
        
        eval(parse(text=paste('p',j,  '<- ggplot(data, aes(x=nestCode, y =', structMetrics[j], ', colour=binNestLevel.factorial, shape=quantNestLevel.factorial)) + 
                              geom_boxplot() + 
                              #geom_point(median.data, aes(x=Group.1, y=', structMetrics[j], '), colour=bnl, shape=qnl) +
                              theme(legend.position="none") +                              
                              theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
                              facet_wrap(c("btwGuildC"), nrow = 1, ncol = length(unique(data$btwGuildC)), scales = "', scaleOp, '",
                              shrink = TRUE, labeller = "label_both", as.table = TRUE,
                              switch = NULL, drop = TRUE, dir = "h", strip.position = "top")', sep='')))
        
        
        
      } # Loop over j
      
      
      grid.arrange(
        p1,p2,p3,p4,p5,p6,p7,p8,
        nrow = length(structMetrics),
        top = paste( 'All network metrics for', i, 'species networks' , sep=' '),
        newpage = T
      )
      
    } # Loop over i
    dev.off()
    }
  
  # ggplot: Nestedness vs complementarity scatter plot (based on rep.frame)
  {
    # Plot options
    scaleOp =  'fixed'# 'free'# 
    
    # initiate plot
    setwd(save.dir)
    # windows(width=8.3, height=11.7)
    pdf(file=paste('qNestMetrics_vs_qFcomp_ScatterPlot_NestProj_ggplot_', scaleOp, 'Y.pdf', sep=''), width=8.3, height=11.7)
    
    for ( i in sort(unique(rep.frame$S)) ){
      
      for ( j in 1:length(qNestMetrics) ){
        
        data <- rep.frame[rep.frame$S==i ,]  

        eval(parse(text=paste('p',j,  '<- ggplot(data, aes(x=qFcomp, y =', qNestMetrics[j], ', colour=binNestLevel.factorial, shape=quantNestLevel.factorial)) + 
                              geom_point(alpha=0.3) + 
                              theme(legend.position="none") +                              
                              theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                              facet_wrap(c("btwGuildC"), nrow = 1, ncol = length(unique(data$btwGuildC)), scales = "', scaleOp, '",
                              shrink = TRUE, labeller = "label_both", as.table = TRUE,
                              switch = NULL, drop = TRUE, dir = "h", strip.position = "top")', sep='')))
      } # Loop over j
      
      grid.arrange(
        p1,p2,p3,
        nrow = length(qNestMetrics),
        top = paste( 'Quantitative nestedness metrics vs qFcomp for', i, 'species networks' , sep=' '),
        newpage = T
      )

    } # Loop over i
    dev.off()
  }
  
  # ggplot: Nestedness vs complementarity geom_smooth plot (based on rep.frame), grouped by nestCode
  {
    # Plot options
    scaleOp =  'fixed'# 'free'# 
    
    # initiate plot
    setwd(save.dir)
    # windows(width=8.3, height=11.7)
    pdf(file=paste('qNestMetrics_vs_qFcomp_SmoothPlot_GroupedByNestCode_NestProj_ggplot_', scaleOp, 'Y.pdf', sep=''), width=8.3, height=11.7)
    
    for ( i in sort(unique(rep.frame$S)) ){
      
      for ( j in 1:length(qNestMetrics) ){
        
        data <- rep.frame[rep.frame$S==i ,]  
        
        eval(parse(text=paste('p',j,  '<- ggplot(data, aes(x=qFcomp, y =', qNestMetrics[j], ', colour=nestCode)) + 
                              geom_smooth(method="auto") + 
                              theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                              facet_wrap(c("btwGuildC"), nrow = 1, ncol = length(unique(data$btwGuildC)), scales = "', scaleOp, '",
                              shrink = TRUE, labeller = "label_both", as.table = TRUE,
                              switch = NULL, drop = TRUE, dir = "h", strip.position = "top")', sep='')))
      } # Loop over j
      
      grid.arrange(
        p1,p2,p3,
        nrow = length(qNestMetrics),
        top = paste( 'Quantitative nestedness metrics vs qFcomp for', i, 'species networks' , sep=' '),
        newpage = T
      )
      
    } # Loop over i
    dev.off()
  }
  
  # ggplot: stab vs comp
  {
    # Plot options
    scaleOp =  'free'#'fixed'#  
    
    # initiate plot
    setwd(save.dir)
    # windows(width=8.3, height=11.7)
    pdf(file=paste('LocalStability_vs_qFcomp_Boxplot_NestProj_ggplot_', scaleOp, 'Y.pdf', sep=''), width=8.3, height=11.7)
    
    for ( i in sort(unique(rep.frame$intTypes)) ){
      for ( j in sort(unique(rep.frame$S[rep.frame$intTypes==i])) ){
        
        data <- rep.frame[rep.frame$intTypes==i & rep.frame$S==j ,]
        
        p1 <- ggplot(data, aes(x=qFcomp.factorial, y = localStability)) + 
          geom_boxplot() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          facet_wrap(c("wGcode", "btwGuildC"), nrow = length(unique(data$wGcode)), ncol = length(unique(data$btwGuildC)), scales = scaleOp,
                     shrink = TRUE, labeller = "label_both", as.table = TRUE,
                     switch = NULL, drop = TRUE, dir = "h", strip.position = "top")
        
        grid.arrange(
          p1,
          nrow = 1,
          top = paste( i , 'with', j, 'species' , sep=' '),
          newpage = T
        )
        
        
      } # Loop over j
    } # Loop over i
    dev.off()
  }
  
  # ggplot: comp vs stab 
  {
    # Plot options
    scaleOp =  'free'#'fixed'#  
    
    # initiate plot
    setwd(save.dir)
    # windows(width=8.3, height=11.7)
    pdf(file=paste('qFcomp_vs_LocalStability_Boxplot_NestProj_ggplot_', scaleOp, 'Y.pdf', sep=''), width=8.3, height=11.7)
    
    for ( i in sort(unique(rep.frame$intTypes)) ){
      for ( j in sort(unique(rep.frame$S[rep.frame$intTypes==i])) ){
        
        data <- rep.frame[rep.frame$intTypes==i & rep.frame$S==j ,]
        
        p1 <- ggplot(data, aes(x=as.factor(localStability), y = qFcomp)) + 
          geom_boxplot() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          facet_wrap(c("wGcode", "btwGuildC"), nrow = length(unique(data$wGcode)), ncol = length(unique(data$btwGuildC)), scales = scaleOp,
                     shrink = TRUE, labeller = "label_both", as.table = TRUE,
                     switch = NULL, drop = TRUE, dir = "h", strip.position = "top")
        
        grid.arrange(
          p1,
          nrow = 1,
          top = paste( i , 'with', j, 'species' , sep=' '),
          newpage = T
        )
        
        
      } # Loop over j
    } # Loop over i
    dev.off()
  }
  
  # ggplot: comp vs stab, grouped by nestCode
  {
    # Plot options
    scaleOp =  'free'#'fixed'#  
    
    # initiate plot
    setwd(save.dir)
    # windows(width=8.3, height=11.7)
    pdf(file=paste('qFcomp_vs_LocalStability_Boxplot_GroupedByNestCode_NestProj_ggplot_', scaleOp, 'Y.pdf', sep=''), width=8.3, height=11.7)
    
    for ( i in sort(unique(rep.frame$intTypes)) ){
      for ( j in sort(unique(rep.frame$S[rep.frame$intTypes==i])) ){
        
        data <- rep.frame[rep.frame$intTypes==i & rep.frame$S==j ,]
        
        p1 <- ggplot(data, aes(x=as.factor(localStability), y = qFcomp, colour = nestCode)) + 
          geom_boxplot() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          facet_wrap(c("wGcode", "btwGuildC"), nrow = length(unique(data$wGcode)), ncol = length(unique(data$btwGuildC)), scales = scaleOp,
                     shrink = TRUE, labeller = "label_both", as.table = TRUE,
                     switch = NULL, drop = TRUE, dir = "h", strip.position = "top")
        
        grid.arrange(
          p1,
          nrow = 1,
          top = paste( i , 'with', j, 'species' , sep=' '),
          newpage = T
        )
        
        
      } # Loop over j
    } # Loop over i
    dev.off()
  }
  
  # ggplot: Stability vs row and col sum metrics
  {
    # Plot options
    scaleOp =  'fixed'# 'free'# 
    
    # initiate plot
    setwd(save.dir)
    # windows(width=8.3, height=11.7)
    pdf(file=paste('LocalStability_vs_maxrsJ_NestProj_ggplot_', scaleOp, 'Y.pdf', sep=''), width=8.3, height=11.7)
    
    for ( i in sort(unique(scen.frame$intTypes)) ){
      for ( j in sort(unique(scen.frame$S[scen.frame$intTypes==i])) ){
        
        data1 <- rep.frame[rep.frame$intTypes==i & rep.frame$S==j ,]
        data2 <- rc.frame[rc.frame$intTypes==i & rc.frame$S==j ,]
        data  <- cbind(data1,data2)
        
        p1 <- ggplot(data, aes(x=max.rsJ, y = localStability)) + 
          geom_point() +
          geom_smooth(method = "glm", method.args = list(family = "binomial"),se = FALSE)  +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          facet_wrap(c("wGcode", "btwGuildC"), nrow = length(unique(data$wGcode)), ncol = length(unique(data$btwGuildC)), scales = scaleOp,
                     shrink = TRUE, labeller = "label_both", as.table = TRUE,
                     switch = NULL, drop = TRUE, dir = "h", strip.position = "top")
        
        grid.arrange(
          p1,
          nrow = 1,
          top = paste( i , 'with', j, 'species' , sep=' '),#paste( i, 'species', wGuildCCode$name[wGuildCCode$value==j], 'within guild competition', sep=' '),
          newpage = T
        )
        
        
      } # Loop over j
    } # Loop over i
    dev.off()
  }
  
  # ggplot: Stability vs row and col sum metrics
  {
    # Plot options
    scaleOp =  'fixed'# 'free'# 
    
    # initiate plot
    setwd(save.dir)
    # windows(width=8.3, height=11.7)
    pdf(file=paste('LocalStability_vs_maxrsJ_histogram_NestProj_ggplot_', scaleOp, 'Y.pdf', sep=''), width=8.3, height=11.7)
    
    for ( i in sort(unique(scen.frame$intTypes)) ){
      for ( j in sort(unique(scen.frame$S[scen.frame$intTypes==i])) ){
        
        data1 <- rep.frame[rep.frame$intTypes==i & rep.frame$S==j ,]
        data2 <- rc.frame[rc.frame$intTypes==i & rc.frame$S==j ,]
        data  <- cbind(data1,data2)
        # data  <- data[data$max.rsJ==1,] 
        
        p1 <- ggplot(data, aes(x=max.rsJ)) + 
          geom_histogram() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          facet_wrap(c("wGcode", "btwGuildC"), nrow = 5, ncol = 3, scales = scaleOp,
                     shrink = TRUE, labeller = "label_both", as.table = TRUE,
                     switch = NULL, drop = TRUE, dir = "h", strip.position = "top")
        
        grid.arrange(
          p1,
          nrow = 1,
          top = paste( i , 'with', j, 'species' , sep=' '),#paste( i, 'species', wGuildCCode$name[wGuildCCode$value==j], 'within guild competition', sep=' '),
          newpage = T
        )
        
        
      } # Loop over j
    } # Loop over i
    dev.off()
  }
  
  # histograms of correlation coefficients of pairwise IS
  {
    # Plot options
    scaleOp =  'fixed'# 'free'# 
    
    # initiate plot
    setwd(save.dir)
    # windows(width=8.3, height=11.7)
    pdf(file=paste('CorrCoef_histograms_groupedByNestCode_NestProj_ggplot_', scaleOp, 'Y.pdf', sep=''), width=8.3, height=11.7)
    
    for ( i in sort(unique(scen.frame$intTypes)) ){
      # for ( j in sort(unique(scen.frame$S[scen.frame$intTypes==i])) ){
      #   
        
        data <- c.frame[c.frame$intTypes==i  ,]
        
        p1 <- ggplot(data, aes(x=corrcoef)) + 
          geom_histogram() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          facet_wrap(c("binNestLevel", "quantNestLevel"), nrow = 3, ncol = 3, scales = scaleOp,
                     shrink = TRUE, labeller = "label_both", as.table = TRUE,
                     switch = NULL, drop = TRUE, dir = "h", strip.position = "top")
        
        grid.arrange(
          p1,
          nrow = 1,
          top = paste( i , sep=' '),#paste( i, 'species', wGuildCCode$name[wGuildCCode$value==j], 'within guild competition', sep=' '),
          newpage = T
        )
        
    } # Loop over i
    dev.off()
  }

}

