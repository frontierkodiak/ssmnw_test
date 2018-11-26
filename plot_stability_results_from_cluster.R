##############################################
##########  PLOT main_ssmnw RESULTS  #########
##############################################
#
#   @DESCRIPTION
#   Plots data from the cluster
#
#   @AUTHOR
#   Alva Curtsdotter, postdoc @ Emory University, 
#   2018-08-05
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
  data.dir <- 'C:/Users/icurtsd/Documents/Project Atlanta/Project Structure and stability of mutualistic networks/data/clusterdata/analysis/2018Aug04_2228/'
  
  save.dir <- data.dir
  if ( !(dir.exists(save.dir)) ){ 
    dir.create(save.dir)
  } 
  
}

### Load and prep data ###
{
  setwd(data.dir)
  res.frame <- read.csv("stabilitySummary.csv", sep=' ')
  str(res.frame)
  
  res.frame$wGcode = NA                                   
  res.frame$wGcode[res.frame$wgcIgnore == 0 & res.frame$wGuildWeightStd == 0.25]  = paste('intraGuildC=',  signif(res.frame$wGuildC[res.frame$wgcIgnore == 0 & res.frame$wGuildWeightStd == 0.25],digits=2), ' and weak wgIS', sep='')
  res.frame$wGcode[res.frame$wgcIgnore == 0 & res.frame$wGuildWeightStd == 0.50]  = paste('intraGuildC=',  signif(res.frame$wGuildC[res.frame$wgcIgnore == 0 & res.frame$wGuildWeightStd == 0.50],digits=2), ' and strong wgIS', sep='')
  res.frame$wGcode[res.frame$wgcIgnore == 1 & res.frame$wGuildWeightStd == 0.25]  = paste('intraGuildC=',  signif(res.frame$btwGuildC[res.frame$wgcIgnore == 1 & res.frame$wGuildWeightStd == 0.25],digits=2), ' and weak wgIS', sep='')
  res.frame$wGcode[res.frame$wgcIgnore == 1 & res.frame$wGuildWeightStd == 0.50]  = paste('intraGuildC=',  signif(res.frame$btwGuildC[res.frame$wgcIgnore == 1 & res.frame$wGuildWeightStd == 0.50],digits=2), ' and strong wgIS', sep='')
  res.frame$wGcode[res.frame$wGuildC == 0 ]                                       = 'no' 
  res.frame$wGcode <- factor(res.frame$wGcode, levels = c('no', sort(unique(res.frame$wGcode))[1:(length(unique(res.frame$wGcode))-1)] ) )
  
  res.frame$wGuildC[res.frame$wgcIgnore == 1] = res.frame$btwGuildC[res.frame$wgcIgnore == 1]
  res.frame$wGuildC <- signif(res.frame$wGuildC, digits=2)
  res.frame$btwGuildC <- signif(res.frame$btwGuildC, digits=2)
  
  res.frame$guildStructCode = NA
  res.frame$guildStructCode[res.frame$guildStruct==0] <- 'no'
  res.frame$guildStructCode[res.frame$guildStruct==1] <- 'yes'
  res.frame$guildStructCode <- factor(res.frame$guildStructCode, levels=c('no','yes'))
  
  #res.frame$intTypes[res.frame$intTypes == 1] = 'competition'
  res.frame$intTypes[res.frame$intTypes == 2] = 'mutualism'
  res.frame$intTypes[res.frame$intTypes == 3] = 'predator-prey'
  
  res.frame$locStabCount <- res.frame$localStability/100*unique(res.frame$trials)
  
  bob <- binom.confint(as.integer(res.frame$locStabCount), rep(unique(res.frame$trials), nrow(res.frame)), methods = "agresti-coull")
  res.frame <- cbind(res.frame, bob[,4:6])
  rm(bob) #cleanup
  res.frame$lower[res.frame$lower < 0] <- 0
  res.frame$upper[res.frame$upper > 1] <- 1
  
}

### PLOTS ###

### Local stability bars ###
{
  
  # Local stability - barplots 
  
  #ggplot for wgcIgnore == 0
  {
    # initiate plot
    setwd(save.dir)
    # windows(width=8.3, height=11.7)
    pdf(file='LocalStabilityBars_GuildProj_wgcIgnoreIsFalse_ggplot_freeY.pdf', width=8.3, height=11.7)
    
    foo <- res.frame[res.frame$wgcIgnore==0, ]
    
    for ( i in sort(unique(foo$intTypes)) ){
      
      for ( j in sort(unique(foo$wGuildWeightStd)) ){
      
        for ( k in sort(unique(foo$S[foo$intTypes==i])) ){
        
        data <- foo[foo$intTypes==i & foo$wGuildWeightStd==j & foo$S==k ,]
        
        limits.all <- aes(ymax = upper, ymin = lower)
        p1 <- ggplot(data, aes(x=guildStruct, y = mean)) +
          geom_bar(stat="identity") + geom_errorbar(limits.all, width=0.2, linetype = 1) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          facet_wrap(c("wGuildC", "btwGuildC"), nrow = length(unique(data$wGuildC)), ncol = length(unique(data$btwGuildC)), scales = "free",#  scales = "fixed",# 
                     shrink = TRUE, labeller = "label_both", as.table = TRUE,
                     switch = NULL, drop = TRUE, dir = "h", strip.position = "top")
        
        grid.arrange(
          p1,
          nrow = 1,
          top = paste( i , 'with', k, 'species and wG IS =', j , sep=' '),
          newpage = T
        )
        
        } # Loop over k
      } # Loop over j
    } # Loop over i
    dev.off()
  }
  
  #ggplot for wgcIgnore == 1
  {
    # initiate plot
    setwd(save.dir)
    # windows(width=8.3, height=11.7)
    pdf(file='LocalStabilityBars_GuildProj_wgcIgnoreIsTrue_ggplot_freeY.pdf', width=8.3, height=11.7)
    
    foo <- res.frame[res.frame$wgcIgnore==1, ]
    
    for ( i in sort(unique(foo$intTypes)) ){
      
      for ( j in sort(unique(foo$wGuildWeightStd)) ){
          
          data <- foo[foo$intTypes==i & foo$wGuildWeightStd==j,]
          
          limits.all <- aes(ymax = upper, ymin = lower)
          p1 <- ggplot(data, aes(x=guildStruct, y = mean)) +
            geom_bar(stat="identity") + geom_errorbar(limits.all, width=0.2, linetype = 1) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            facet_wrap(c("S", "btwGuildC"), nrow = length(unique(data$S)), ncol = length(unique(data$btwGuildC[data$S==unique(data$S)[1]])), scales = "free",#    scales = "fixed",#  
                       shrink = TRUE, labeller = "label_both", as.table = TRUE,
                       switch = NULL, drop = TRUE, dir = "h", strip.position = "top")
          
          grid.arrange(
            p1,
            nrow = 1,
            top = paste( i , 'with wG IS =', j , sep=' '),#paste( i, 'species', wGuildCCode$name[wGuildCCode$value==j], 'within guild competition', sep=' '),
            newpage = T
          )
          
      } # Loop over j
    } # Loop over i
    dev.off()
  }
}

### Interaction plots ###
{
  #ggplot for wgcIgnore == 1
  {
    # initiate plot
    setwd(save.dir)
    # windows(width=8.3, height=11.7)
    pdf(file='InteractionPlots_GuildProj_wgcIgnoreIsTrue_ggplot_freeY.pdf', width=8.3, height=11.7)
    
    foo <- res.frame[res.frame$wgcIgnore==1, ]
    
    for ( i in sort(unique(foo$intTypes)) ){
          
          data <- foo[foo$intTypes==i ,]
          
          limits.all <- aes(ymax = upper, ymin = lower)
          p1 <- ggplot(data, aes(x=guildStructCode, y = mean, colour = wGuildC)) + 
            geom_point(size = 3) +
            geom_line(aes(group = interaction(wGuildC))) +
            geom_errorbar(limits.all, width=0.2, linetype = 1) + 
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            xlab("Guild Structured") + 
            ylab("Proportion of stable networks") + 
            labs(color = "within-guild connectance") +
            facet_wrap(c("S", "wGuildWeightStd"), nrow = length(unique(data$S)), ncol = length(unique(data$wGuildWeightStd)), scales = "free",#   scales = "fixed",# 
                       shrink = TRUE, labeller = "label_both", as.table = TRUE,
                       switch = NULL, drop = TRUE, dir = "h", strip.position = "top")
          
          grid.arrange(
            p1,
            nrow = 1,
            top = paste( i , sep=' '),
            newpage = T
          )
        
    } # Loop over i
    dev.off()
  }
  
  #ggplot for wGuildC == 0
  {
    # initiate plot
    setwd(save.dir)
    # windows(width=8.3, height=11.7)
    pdf(file='InteractionPlots_GuildProj_wGuildC0_ggplot_freeY.pdf', width=8.3, height=11.7)
    
    data <- res.frame[res.frame$wGuildC==0, ]
      
      limits.all <- aes(ymax = upper, ymin = lower)
      p1 <- ggplot(data, aes(x=guildStructCode, y = mean, colour = btwGuildC)) + 
        geom_point(size = 3) +
        geom_line(aes(group = interaction(btwGuildC))) +
        geom_errorbar(limits.all, width=0.2, linetype = 1) + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        xlab("Guild Structured") + 
        ylab("Proportion of stable networks") + 
        labs(color = "between-guild connectance") +
        facet_wrap(c("S", "intTypes"), nrow = length(unique(data$S)), ncol = length(unique(data$intTypes)), scales = "free",#   scales = "fixed",# 
                   shrink = TRUE, labeller = "label_both", as.table = TRUE,
                   switch = NULL, drop = TRUE, dir = "h", strip.position = "top")
      
      grid.arrange(
        p1,
        nrow = 1,
        newpage = T
      )

    dev.off()
  }
  
  ##ggplot ESA presentation Fig: Effect of guildStruct Try 1
  {
    # initiate plot
    setwd(save.dir)
    # windows(width=8.3, height=11.7)
    pdf(file='InteractionPlots_GuildProj_testESA_MutScens_ggplot_freeY.pdf', width=8.3, height=11.7)
    
    foo1a <- res.frame[res.frame$wGuildC==0 & res.frame$intTypes=='mutualism' & is.element(res.frame$btwGuildC, c(0.2, 0.56)) & res.frame$wGuildWeightStd==0.25, ]
    foo1a <- res.frame[res.frame$wGuildC==0 & is.element(res.frame$btwGuildC, c(0.2, 0.56)), ]
    foo1b <- foo1a
    foo1b$wGuildWeightStd <- 0.5
    foo1 <- rbind(foo1a,foo1b)
    foo2 <- res.frame[res.frame$wgcIgnore==1 & is.element(res.frame$btwGuildC, c(0.2, 0.56)), ]
    data <- rbind(foo1,foo2)
    data$wgCompCode <- NA
    data$wgCompCode[data$wGuildC==0] <- 'no within-guild competition'
    data$wgCompCode[data$wGuildC!=0] <- 'within-guild competition'
    
    
    limits.all <- aes(ymax = upper, ymin = lower)
      p1 <- ggplot(data, aes(x=guildStructCode, y = mean, colour = btwGuildC, linetype = wgCompCode, shape = wgCompCode)) + 
        geom_point(size = 3) +
        geom_line(aes(group = interaction(wgCompCode))) +
        geom_errorbar(limits.all, width=0.2, linetype = 1) + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        xlab("Guild Structured") + 
        ylab("Proportion of stable networks") + 
        labs(color = "between-guild connectance") +
        facet_wrap(c("btwGuildC", "wGuildWeightStd"), nrow = length(unique(data$btwGuildC)), ncol = length(unique(data$wGuildWeightStd)), scales = "free",#   scales = "fixed",# 
                   shrink = TRUE, labeller = "label_both", as.table = TRUE,
                   switch = NULL, drop = TRUE, dir = "h", strip.position = "top")
      
      grid.arrange(
        p1,
        nrow = 1,
        top = paste( i , sep=' '),
        newpage = T
      )
    dev.off()
  }
  
  ##ggplot ESA presentation Fig: Effect of guildStruct Try 2
  {
    # initiate plot
    setwd(save.dir)
    # windows(width=8.3, height=11.7)
    pdf(file='InteractionPlots_GuildProj_testESA_guildStructEffect_ggplot_freeY.pdf', width=8.3, height=11.7)
    
    foo1a <- res.frame[res.frame$wGuildC==0 & res.frame$intTypes=='mutualism' & is.element(res.frame$btwGuildC, c(0.2, 0.56)) & res.frame$wGuildWeightStd==0.25, ]
    foo1b <- res.frame[res.frame$wGuildC==0 & res.frame$intTypes=='predator-prey' & is.element(res.frame$btwGuildC, c(0.016, 0.60)) & res.frame$wGuildWeightStd==0.25, ]
    foo1 <- rbind(foo1a,foo1b)
    foo2 <- res.frame[res.frame$wgcIgnore==1 & is.element(res.frame$btwGuildC, c(0.2, 0.56,0.016, 0.60)) & res.frame$wGuildWeightStd==0.25, ]
    data <- rbind(foo1,foo2)
    data$wgCompCode <- NA
    data$wgCompCode[data$wGuildC==0] <- 'no'
    data$wgCompCode[data$wGuildC!=0] <- 'yes'
    
    
    limits.all <- aes(ymax = upper, ymin = lower)
    p1 <- ggplot(data, aes(x=guildStructCode, y = mean, colour = wgCompCode, shape = wgCompCode)) + 
      geom_point(size = 3) +
      geom_line(aes(group = interaction(wgCompCode))) +
      geom_errorbar(limits.all, width=0.2, linetype = 1) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      xlab("Guild Structured") + 
      ylab("Proportion of stable networks") + 
      labs(color = "presence of within-guild competition") +
      facet_wrap(c('intTypes',"btwGuildC"), nrow = 2, ncol = 2, scales = "free",#   scales = "fixed",# 
                 shrink = TRUE, labeller = "label_both", as.table = TRUE,
                 switch = NULL, drop = TRUE, dir = "h", strip.position = "top")
    
    grid.arrange(
      p1,
      nrow = 1,
      newpage = T
    )
    dev.off()
  }
  
  ##ggplot ESA presentation Fig 1.
  {
    # initiate plot
    setwd(save.dir)
    #windows(width=8.3, height=11.7/2)
    pdf(file='InteractionPlots_GuildProj_ESAfig1.pdf', width=2*8.3/3, height=11.7/2)
    
    foo1a <- res.frame[res.frame$wGuildC==0 & res.frame$intTypes=='mutualism' & is.element(res.frame$btwGuildC, c(0.2, 0.56)) & res.frame$wGuildWeightStd==0.25, ]
    foo1b <- res.frame[res.frame$wGuildC==0 & res.frame$intTypes=='predator-prey' & is.element(res.frame$btwGuildC, c(0.016, 0.60)) & res.frame$wGuildWeightStd==0.25, ]
    foo1 <- rbind(foo1a,foo1b)
    foo2 <- res.frame[res.frame$wgcIgnore==1 & is.element(res.frame$btwGuildC, c(0.2, 0.56,0.016, 0.60)) & res.frame$wGuildWeightStd==0.25, ]
    data <- rbind(foo1,foo2)
    data$wgCompCode <- NA
    data$wgCompCode[data$wGuildC==0] <- 'no'
    data$wgCompCode[data$wGuildC!=0] <- 'yes'
    
    limits.all <- aes(ymax = upper, ymin = lower)
    
    p1 <- ggplot(data[data$intTypes=='mutualism' & data$btwGuildC==0.2,], aes(x=guildStructCode, y = mean, colour = wgCompCode)) + 
      geom_point(size = 3) +
      theme(legend.position="none") +
      geom_line(aes(group = interaction(wgCompCode))) +
      geom_errorbar(limits.all, width=0.2, linetype = 1) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      xlab("Guild Structured") + 
      ylab("Proportion of stable networks") + 
      ylim(0.6,1)
    
    p2 <- ggplot(data[data$intTypes=='mutualism' & data$btwGuildC==0.56,], aes(x=guildStructCode, y = mean, colour = wgCompCode)) + 
      geom_point(size = 3) +
      theme(legend.position="none") +
      geom_line(aes(group = interaction(wgCompCode))) +
      geom_errorbar(limits.all, width=0.2, linetype = 1) + 
      theme( axis.text.x = element_text(angle = 45, hjust = 1)) +
      xlab("Guild Structured") + 
      ylab("Proportion of stable networks") + 
      ylim(0,0.1)
    
    p3 <- ggplot(data[data$intTypes=='predator-prey' & data$btwGuildC==0.016,], aes(x=guildStructCode, y = mean, colour = wgCompCode)) + 
      geom_point(size = 3) +
      theme(legend.position="none") +
      geom_line(aes(group = interaction(wgCompCode))) +
      geom_errorbar(limits.all, width=0.2, linetype = 1) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      xlab("Guild Structured") + 
      ylab("Proportion of stable networks") + 
      ylim(0.9,1)
      
    p4 <- ggplot(data[data$intTypes=='predator-prey' & data$btwGuildC==0.6,], aes(x=guildStructCode, y = mean, colour = wgCompCode)) + 
      geom_point(size = 3) +
      theme(legend.position="none") +
      geom_line(aes(group = interaction(wgCompCode))) +
      geom_errorbar(limits.all, width=0.2, linetype = 1) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      xlab("Guild Structured") + 
      ylab("Proportion of stable networks") + 
      ylim(0,1)
    
    grid.arrange(
      p1,p2,p3,p4,
      nrow = 2, ncol =2 ,
      newpage = T
    )
    dev.off()
  }
}

### Heatmaps ###
{}