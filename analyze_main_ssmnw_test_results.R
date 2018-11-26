##############################################
########## ANALYSE main_ssmnw TEST RESULTS #########
##############################################
#
#   @DESCRIPTION
#   Analysis code.
#   Summarize the results for local stability,
#   reactivity, and resilience, for a large 
#   number of replicates for a set of network
#   structures.
#   Written for main_ssmnw. First time using ggplot.
#
#   @AUTHOR
#   Alva Curtsdotter, postdoc @ Emory University, 
#   2018-06-20
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
  # To read .mat files and such
  require('R.matlab')
  # Network structure analysis
  require(bipartite)
  require(vegan)
  
  # Data paths and such
  data.dir <- 'C:/Users/icurtsd/Documents/Project Atlanta/Project Structure and stability of mutualistic networks/data/main_test/testresults_2018Jul13T0105/'  
  
  resultVarNames <- c('localStability', 'reactivity', 'resilience')
  structMetrics  <- c('bNestNODF', 'bNestTemp', 'bNestSpecRad', 'qNestNODF', 'qNestWine', 'qNestSpecRad', 'bFcomp', 'qFcomp')
  
  save.dir <- paste('C:/Users/icurtsd/Documents/Project Atlanta/Project Structure and stability of mutualistic networks/data/main_test/testanalysis_', substr(data.dir, regexpr('2018', data.dir), nchar(data.dir)), sep='')
  
  if ( length(dir(save.dir)) == 0 ){ #this for some reason doesn't work...but it's fine anyway.
    dir.create(save.dir)
  } 
  
  # Options
  analyseStability <- 1
  analyseBtwGuildStructure <- 1 # Note that if 1, there are some option to define at the start of this section.
  analyseJacobianStructure <- 0 # Not implemented yet.
  
  # Some specific stuff
  intTypesCode <- data.frame( value=c(1,2,3), name=c('competition','mutualism','trophic'), stringsAsFactors = FALSE)
  guildStructCode <- data.frame( value=c(0,1), name=c('without','with'), stringsAsFactors = FALSE)
  wGuildCCode <- data.frame( value=c(0,1), name=c('without','with'), stringsAsFactors = FALSE)
}

### Load and analyse data ###
{
  setwd(data.dir)
  
  # Read parameterStruct for the simulation
  file1  <- dir(pattern = 'paramStruct')
  params <- readMat(file1)
  
  # Get the paramStruct into shape.
  params.data <- lapply(params$params[2], unlist, use.names=FALSE)
  params.data <- as.data.frame(params.data)
  names(params.data) <- c('scenID', lapply(params$params[1], unlist, use.names=FALSE)[[1]])
  str(params.data)
  # params.data$scenID <- 1:nrow(params.data)
  
  params.data$nestCode = NA
  params.data$nestCode[params.data$binNestLevel == 1 & params.data$quantNestLevel == 1]  = 'bNqN'
  params.data$nestCode[params.data$binNestLevel == 1 & params.data$quantNestLevel == 0]  = 'bNqR'
  params.data$nestCode[params.data$binNestLevel == 1 & params.data$quantNestLevel == -1] = 'bNqA'
  params.data$nestCode[params.data$binNestLevel == 0 & params.data$quantNestLevel == 1]  = 'bRqN'
  params.data$nestCode[params.data$binNestLevel == 0 & params.data$quantNestLevel == 0]  = 'bRqR'
  params.data$nestCode[params.data$binNestLevel == 0 & params.data$quantNestLevel == -1] = 'bRqA'
  params.data$nestCode[params.data$binNestLevel == -1 & params.data$quantNestLevel == 1]  = 'bAqN'
  params.data$nestCode[params.data$binNestLevel == -1 & params.data$quantNestLevel == 0]  = 'bAqR'
  params.data$nestCode[params.data$binNestLevel == -1 & params.data$quantNestLevel == -1] = 'bAqA'
  params.data$nestCode <- factor(params.data$nestCode, levels = c('bNqN', 'bNqR', 'bNqA', 'bRqN', 'bRqR', 'bRqA','bAqN', 'bAqR', 'bAqA'))

  
  if ( analyseStability == 1 ){
    
    # List all result files
    files <- dir(pattern='ResMat*')
    str(files)
    
    # Check that number of result files equals the number of scenarios.
    if ( length(files) != nrow(params.data)  ) print('Dimension error of result data!')
    
    # Create the data frame to store the results for all scenarios and replicates!
    foo <- readMat(files[1])$resMat
    res.frame <- as.data.frame(matrix(NA, ncol=ncol(foo)+ncol(params.data), nrow=nrow(foo)*nrow(params.data)))
    names(res.frame) <- c(names(params.data), resultVarNames)
    
    # Load all results into res.frame
    for ( i in 1:length(files)){
      
      file <- files[i]
      
      foo <- readMat(file)$resMat
      str(foo)
      
      scenID <- eval(parse(text=(substr(file, nchar("ResMat_scenID_")+1, nchar(file)-nchar(".mat")))))
      
      res.frame[seq(i*nrow(foo)-nrow(foo)+1, i*nrow(foo), 1), 1:ncol(params.data)] <- params.data[params.data$scenID==scenID,]
      res.frame$trials[seq(i*nrow(foo)-nrow(foo)+1, i*nrow(foo), 1)] <- 1:nrow(foo)
      res.frame[seq(i*nrow(foo)-nrow(foo)+1, i*nrow(foo), 1), seq(1+ncol(params.data),ncol(foo)+ncol(params.data),1)] <- foo
      
    }
    # res.frame <- res.frame[order(res.frame$scenID),]
    str(res.frame)
    
    # setwd(save.dir)
    # write.table(res.frame, 'stabilityAnalysis.csv')
  }
  
  if ( analyseBtwGuildStructure == 1 ){
    
    # Define non-redundant scenario list and number of replicates to run
    scenList <- NA
    sLevels <- sort(unique(params.data$S))
    nestLevels <- sort(unique(params.data$nestCode))
    n <-1
    for (i in sLevels){
      sub <- params.data
      sub <- sub[sub$S==i,]
      bgcLevels <-  sort(unique(sub$btwGuildC))
      for (j in bgcLevels){
        for (k in nestLevels){
          sub <- params.data
          sub <- sub[sub$S==i,]
          sub <- sub[sub$btwGuildC==j,]
          sub <- sub[sub$nestCode==k,]
          scenList[n] <- sub$scenID[1]
          n <- n+1
        }
      }
    }
    scenList <- sort(scenList)
    reps <- 200
    
    # Create the data frame to store the results for all scenarios and replicates!
    bipStruct.frame <- as.data.frame(matrix(NA, ncol=1+length(structMetrics), nrow=length(scenList)*reps))
    names(bipStruct.frame) <- c(names(params.data)[1], structMetrics)
    str(bipStruct.frame)
    
    # Load all results into bipStruct.frame
    n <- 1
    for ( i in scenList){
      for (j in 1:reps){
        
        file <- paste('Jacobian_scenID_', i, '_replicate_', j, '.csv', sep='')
        
        bip <- as.matrix(read.table(file, sep=","))
        bip <- bip[(1+nrow(bip)/2):nrow(bip),1:(ncol(bip)/2)]
        # str(bip)
        
        bbip <- bip
        bbip[bbip!=0] <- 1
        
        qbip <- abs(bip)
        
        bipStruct.frame$scenID[n] <- i   #eval(parse(text=(substr(file, nchar("Jacobian_scenID_")+1, regexpr('_replicate_', file)[1]-1))))
        
        # Binary Nestedness metrics
        bipStruct.frame$bNestTemp[n] <- nested(bbip, method='binmatnest2') 
        bipStruct.frame$bNestNODF[n] <- nested(bbip, method='NODF2')
        
        zeroblock    <- matrix(0,nrow=nrow(bbip),ncol=ncol(bbip))
        A            <- cbind(rbind(zeroblock,bbip), rbind(t(bbip),zeroblock))
        bipStruct.frame$bNestSpecRad[n] <- max(eigen(A)$values)
        
        # Quantitative nestedness metrics
        bipStruct.frame$qNestNODF[n] <- nested(qbip, method='weighted NODF')
        bipStruct.frame$qNestWine[n] <- wine(qbip, nreps = 100)$wine
        
        zeroblock    <- matrix(0,nrow=nrow(qbip),ncol=ncol(qbip))
        A            <- cbind(rbind(zeroblock,qbip), rbind(t(qbip),zeroblock))
        bipStruct.frame$qNestSpecRad[n] <- max(eigen(A)$values)
        
        # Binary complementarity
        ppm <- bbip
        distance = vegdist(t(ppm), method="euclidean", bin=F)
        distance[is.nan(distance)] <- 0
        tree <- hclust(distance)
        bipStruct.frame$bFcomp[n] <- treeheight(tree)
        
        # Quantitative complementarity
        ppm <- qbip
        distance = vegdist(t(ppm), method="euclidean", bin=F)
        distance[is.nan(distance)] <- 0
        tree <- hclust(distance)
        bipStruct.frame$qFcomp[n] <- treeheight(tree)
        
        n <- n+1
      }
    }
    
    # Save bipStruct.frame since it takes so long to do this analysis
    setwd(save.dir)
    write.table(bipStruct.frame, paste('bipNwStructAnalysis_', length(scenList), '_scenarios_', reps, '_replicates.csv', sep=''))
  }
  
}

# Temp ditty to salvage some scenarios. All commented away.
{
  # afo <- bipStruct.frame
  # afo <- afo[,!(names(afo) %in% c("bNestNODF", "qNestNODF"))]
  # afo <- na.omit(afo)
  # str(afo)
  # 
  # run <- unique(afo$scenID)
  # subOrg <- params.data[is.element(params.data$scenID, run ),]
  # 
  # subOrg$nestCode = NA
  # subOrg$nestCode[subOrg$binNestLevel == 1 & subOrg$quantNestLevel == 1]  = 'bNqN'
  # subOrg$nestCode[subOrg$binNestLevel == 1 & subOrg$quantNestLevel == 0]  = 'bNqR'
  # subOrg$nestCode[subOrg$binNestLevel == 1 & subOrg$quantNestLevel == -1] = 'bNqA'
  # subOrg$nestCode[subOrg$binNestLevel == 0 & subOrg$quantNestLevel == 1]  = 'bRqN'
  # subOrg$nestCode[subOrg$binNestLevel == 0 & subOrg$quantNestLevel == 0]  = 'bRqR'
  # subOrg$nestCode[subOrg$binNestLevel == 0 & subOrg$quantNestLevel == -1] = 'bRqA'
  # subOrg$nestCode[subOrg$binNestLevel == -1 & subOrg$quantNestLevel == 1]  = 'bAqN'
  # subOrg$nestCode[subOrg$binNestLevel == -1 & subOrg$quantNestLevel == 0]  = 'bAqR'
  # subOrg$nestCode[subOrg$binNestLevel == -1 & subOrg$quantNestLevel == -1] = 'bAqA'
  # subOrg$nestCode <- factor(subOrg$nestCode, levels = c('bNqN', 'bNqR', 'bNqA', 'bRqN', 'bRqR', 'bRqA','bAqN', 'bAqR', 'bAqA'))
  # 
  # save.frame <- as.data.frame(matrix(NA, ncol=ncol(afo), nrow=135000))
  # names(save.frame) <- names(afo)
  # str(save.frame)
  # 
  # counter <- 1
  # sLevels <- sort(unique(subOrg$S))
  # for (i in sLevels){
  #   sub <- subOrg
  #   sub <- sub[sub$S==i,]
  #   bgcLevels <-  sort(unique(sub$btwGuildC))
  #   for (j in bgcLevels){
  #     sub <- sub[sub$btwGuildC==j,]
  #     nestLevels <- sort(unique(sub$nestCode))
  #     for (k in nestLevels){
  #       sub <- sub[sub$nestCode==k,]
  #       scens <- sort(unique(sub$scenID))
  #       for (l in scens){
  #         if (sum(afo$scenID==l, na.rm = T)==1000){
  #           tosave <- afo[ afo$scenID == l,]
  #           save.frame[seq(counter, counter+nrow(tosave)-1, 1), 1:ncol(tosave)] <- tosave
  #           counter <- counter + nrow(tosave)
  #         }
  #       }
  #     }
  #   }
  # }
}

### Plot results ###
{
  # Plot stability results only
  if ( analyseStability == 1 ){
  
  if ( res.frame$project == 1 ){
    
    ### Plots for guild structure evaluation ###
    {
      
      # Local stability - barplots
      # ggplot
      {
        # initiate plot
        setwd(save.dir)
        # windows(width=8.3, height=11.7)
        pdf(file='LocalStabilityBars_GuildProj_ggplot.pdf', width=8.3, height=11.7)
        
        for ( i in unique(res.frame$intTypes) ){
          
          for ( j in unique(res.frame$guildStruct) ){
            
            data <- res.frame[res.frame$intTypes==i & res.frame$guildStruct==j ,]
            
            p1 <- ggplot(data, aes(x=btwGuildC, y = localStability)) + geom_bar(stat="identity") +
              scale_x_continuous( breaks = unique(data$btwGuildC), labels = scales::percent ) + 
              facet_wrap(c("S", "wGuildC"), nrow = 3, ncol = 3, scales = "fixed",
                         shrink = TRUE, labeller = "label_both", as.table = TRUE,
                         switch = NULL, drop = TRUE, dir = "h", strip.position = "top")
            
            grid.arrange(
              p1,
              nrow = 1,
              top = paste('interaction type = ', i, ' and guild structure = ', j, sep=''),
              newpage = T
            )
            
            
          }
        }
        dev.off()
        
        
        # Basic style
        
        # initiate plot
        setwd(save.dir)
        # windows(width=8.3, height=11.7)
        pdf(file='LocalStabilityBars_GuildProj_basicplot.pdf', width=8.3, height=11.7)
        
        # Get number of levels for main variables
        numComplexityLevels <- c(length(unique(res.frame$S)), length(unique(res.frame$wGuildC)), length(unique(res.frame$btwGuildC)))
        names(numComplexityLevels) <- c("S", "wGuildC", "btwGuildC")
        numComplexityLevels <-  numComplexityLevels[order( numComplexityLevels, decreasing = F)]
        
        # Use the number of levels above to structure plot
        plotRows <- numComplexityLevels[1]
        plotCols <- numComplexityLevels[2]
        par(mfrow=c(plotRows,plotCols))
        
        data <- res.frame %>%
          group_by(intTypes, guildStruct, S, wGuildC, btwGuildC) %>%
          summarize(localStability=signif(sum(localStability)/n(), digits=2))
        
        ymaxMat <- data %>%
          group_by(intTypes, guildStruct, S) %>%
          summarize(value=max(localStability))
        
        
        for ( i in unique(res.frame$intTypes) ){
          
          for ( j in unique(res.frame$guildStruct) ){
            
            for ( k in eval(parse(text=paste( 'unique(res.frame$', names(numComplexityLevels)[1], ')', sep=''))) ){
              
              for ( l in eval(parse(text=paste( 'unique(res.frame$', names(numComplexityLevels)[2], ')', sep='')))  ){
                
                plotdata <- data[data$intTypes==i & data$guildStruct ==j & eval(parse(text=paste( 'data$', names(numComplexityLevels)[1], '==k', sep=''))) & eval(parse(text=paste( 'data$', names(numComplexityLevels)[2], '==l', sep=''))),]
                
                ymaxUse <- ymaxMat$value[ymaxMat$intTypes==i & ymaxMat$guildStruct ==j & eval(parse(text=paste( 'ymaxMat$', names(numComplexityLevels)[1], '==k', sep='')))]
                
                bph <- barplot(plotdata$localStability, ylim=c(0,ymaxUse),las=2, names.arg=unique(plotdata$btwGuildC), 
                               xlab = "btwGuildC", ylab='local stability')
                title(main=paste( intTypesCode$name[intTypesCode$value==i], guildStructCode$name[guildStructCode==j], 'guild structure', sep=' '),
                      sub=paste(names(numComplexityLevels)[1], '=', k, ' and ', names(numComplexityLevels)[2], '=', l, sep='' ))
                
                text(bph, ymaxUse, paste(signif(100*plotdata$localStability,digit=2), '%',sep=''), pos=1 )
                
              }
            }
          }
        }
        dev.off()
      }
      
      ### REACTIVITY - BOX PLOTS commented away
      
      {
        # setwd(save.dir)
        # #windows(width=8.3, height=11.7)
        # pdf(file='ReactivityBoxes_GuildProj.pdf', width=8.3, height=11.7)
        # plotRows <- 3
        # plotCols <- 3
        # par(mfrow=c(plotRows,plotCols))
        # 
        # setwd(data.dir)
        # for (i in 1:length(unique(res.frame$intSym))){ # 
        #   
        #   v1 <- unique(res.frame$intSym)[i]
        #   
        #   for (j in 1:length(unique(res.frame$intType))){
        #     
        #     v2 <- unique(res.frame$intType)[j]
        #     
        #     for (k in 1:length(unique(res.frame$compType))){
        #       
        #       v3 <- unique(res.frame$compType)[k]
        #       
        #       data <- as.data.frame(matrix(NA, ncol=2, nrow=data.nrow*length(unique(res.frame$network)) ))
        #       colnames(data) <- c('network', 'reactivity')
        #       
        #       
        #       for (l in 1:length(unique(res.frame$network))){
        #         
        #         v4 <- unique(res.frame$network)[l]
        #         
        #         file <- paste("ToyResMat_IntType_", v2, "_IntSym_", v1, "_CompType_", v3 ,"_Network_", v4, ".mat", sep="" )
        #         data$reactivity[((l-1)*data.nrow+1):(l*data.nrow)] <- readMat(file)$resMat[,2]
        #         data$network[((l-1)*data.nrow+1):(l*data.nrow)] <- v4
        #         
        #       }
        #       
        #       ymin <- min(res.frame$minReactivity[res.frame$intSym==v1 & res.frame$intType==v2])
        #       ymax <- max(res.frame$maxReactivity[res.frame$intSym==v1 & res.frame$intType==v2])
        #       
        #       
        #       bxph <- boxplot(reactivity ~ network, data, ylim=c(ymin, ymax), las=2 )
        #       title(paste(substr(v2,1,4),
        #                   substr(v1,1,3),
        #                   v3, sep=' ') )
        #       
        #     }
        #   }
        # }
        # setwd(save.dir)
        # dev.off()
      }
      
      
      ### RESILIENCE - BOX PLOTS commented away
      
      {
        # setwd(save.dir)
        # #windows(width=8.3, height=11.7)
        # pdf(file='ResilienceBoxes_GuildProj.pdf', width=8.3, height=11.7)
        # plotRows <- 3
        # plotCols <- 3
        # par(mfrow=c(plotRows,plotCols))
        # 
        # setwd(data.dir)
        # for (i in 1:length(unique(res.frame$intSym))){ # 
        #   
        #   v1 <- unique(res.frame$intSym)[i]
        #   
        #   for (j in 1:length(unique(res.frame$intType))){
        #     
        #     v2 <- unique(res.frame$intType)[j]
        #     
        #     for (k in 1:length(unique(res.frame$compType))){
        #       
        #       v3 <- unique(res.frame$compType)[k]
        #       
        #       data <- as.data.frame(matrix(NA, ncol=2, nrow=data.nrow*length(unique(res.frame$network)) ))
        #       colnames(data) <- c('network', 'resilience')
        #       
        #       
        #       for (l in 1:length(unique(res.frame$network))){
        #         
        #         v4 <- unique(res.frame$network)[l]
        #         
        #         file <- paste("ToyResMat_IntType_", v2, "_IntSym_", v1, "_CompType_", v3 ,"_Network_", v4, ".mat", sep="" )
        #         data$resilience[((l-1)*data.nrow+1):(l*data.nrow)] <- readMat(file)$resMat[,3]
        #         data$network[((l-1)*data.nrow+1):(l*data.nrow)] <- v4
        #         
        #       }
        #       
        #       ymin <- min(res.frame$minResilience[res.frame$intSym==v1 & res.frame$intType==v2])
        #       ymax <- max(res.frame$maxResilience[res.frame$intSym==v1 & res.frame$intType==v2])
        #       
        #       
        #       bxph <- boxplot(resilience ~ network, data, ylim=c(ymin, ymax), las=2 )
        #       title(paste(substr(v2,1,4),
        #                   substr(v1,1,3),
        #                   v3, sep=' ') )
        #       
        #     }
        #   }
        # }
        # setwd(save.dir)
        # dev.off()
      }
      
    }
    
  } else if (res.frame$project == 2 ){
    ### Plots for nestedness exploration ###
    {
      # Some res.frame prep
      {
      res.frame$nestCode = NA
      
      res.frame$nestCode[res.frame$binNestLevel == 1 & res.frame$quantNestLevel == 1]  = 'bNqN'
      res.frame$nestCode[res.frame$binNestLevel == 1 & res.frame$quantNestLevel == 0]  = 'bNqR'
      res.frame$nestCode[res.frame$binNestLevel == 1 & res.frame$quantNestLevel == -1] = 'bNqA'
      
      res.frame$nestCode[res.frame$binNestLevel == 0 & res.frame$quantNestLevel == 1]  = 'bRqN'
      res.frame$nestCode[res.frame$binNestLevel == 0 & res.frame$quantNestLevel == 0]  = 'bRqR'
      res.frame$nestCode[res.frame$binNestLevel == 0 & res.frame$quantNestLevel == -1] = 'bRqA'
      
      res.frame$nestCode[res.frame$binNestLevel == -1 & res.frame$quantNestLevel == 1]  = 'bAqN'
      res.frame$nestCode[res.frame$binNestLevel == -1 & res.frame$quantNestLevel == 0]  = 'bAqR'
      res.frame$nestCode[res.frame$binNestLevel == -1 & res.frame$quantNestLevel == -1] = 'bAqA'
      res.frame$nestCode <- factor(res.frame$nestCode, levels = c('bNqN', 'bNqR', 'bNqA', 'bRqN', 'bRqR', 'bRqA','bAqN', 'bAqR', 'bAqA'))
      
      res.frame$wGcode = NA
      res.frame$wGcode[res.frame$wGuildC == 0 ]                                     = 'no'
      res.frame$wGcode[res.frame$wGuildC == 1 & res.frame$wGuildWeightStd == 0.25]  = 'intraGuildC=1 and weak'
      res.frame$wGcode[res.frame$wGuildC == 1 & res.frame$wGuildWeightStd == 0.50]  = 'intraGuildC=1 and strong'
      res.frame$wGcode[res.frame$wgcIgnore == 1 & res.frame$wGuildWeightStd == 0.25]  = 'intraGuildC=btwGuildC and weak'
      res.frame$wGcode[res.frame$wgcIgnore == 1 & res.frame$wGuildWeightStd == 0.50]  = 'intraGuildC=btwGuildC and strong'
      res.frame$wGcode <- factor(res.frame$wGcode, levels = c('no', 'intraGuildC=btwGuildC and weak', 'intraGuildC=1 and weak', 'intraGuildC=btwGuildC and strong', 'intraGuildC=1 and strong'))
      
      # res.frame$wGuildC[res.frame$wgcIgnore == 1] = res.frame$btwGuildC[res.frame$wgcIgnore == 1]
      
      res.frame$intTypes[res.frame$intTypes == 1] = 'competition'
      res.frame$intTypes[res.frame$intTypes == 2] = 'mutualism'
      res.frame$intTypes[res.frame$intTypes == 3] = 'predator-prey'
      
    }
    
      # Local stability - barplots 
      
      #ggplot
      {
        # initiate plot
        setwd(save.dir)
        # windows(width=8.3, height=11.7)
        # pdf(file='LocalStabilityBars_NestProj_weakerWithinGuildIS_ggplot.pdf', width=8.3, height=11.7)
        pdf(file='LocalStabilityBars_NestProj_ggplot_freeY.pdf', width=8.3, height=11.7)
        
        for ( i in sort(unique(res.frame$intTypes)) ){
          # for ( i in sort(unique(res.frame$S)) ){
          
          # for ( j in unique(res.frame$wGuildC) ){
          # for ( j in unique(res.frame$wGcode) ){ 
          for ( j in sort(unique(res.frame$S[res.frame$intTypes==i])) ){
            
            # data <- res.frame[res.frame$intTypes==i & res.frame$wGuildC==j ,]
            # data <- res.frame[res.frame$S==i & res.frame$wGuildC==j ,]
            # data <- res.frame[res.frame$S==i & res.frame$wGcode==j ,]
            # data <- res.frame[res.frame$S==i ,]  
            data <- res.frame[res.frame$intTypes==i & res.frame$S==j ,]
            
            # p1 <- ggplot(data, aes(x=nestCode, y = localStability)) + geom_bar(stat="identity") +
            #   facet_wrap(c("S", "btwGuildC"), nrow = length(unique(data$S)), ncol = length(unique(data$btwGuildC)), scales = "fixed",
            #              shrink = TRUE, labeller = "label_both", as.table = TRUE,
            #              switch = NULL, drop = TRUE, dir = "h", strip.position = "top")
            # p1 <- ggplot(data, aes(x=nestCode, y = localStability)) + geom_bar(stat="identity") +
            #   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            #   facet_wrap(c("intTypes", "btwGuildC"), nrow = length(unique(data$intTypes)), ncol = length(unique(data$btwGuildC)), scales = "fixed",
            #              shrink = TRUE, labeller = "label_both", as.table = TRUE,
            #              switch = NULL, drop = TRUE, dir = "h", strip.position = "top")
            p1 <- ggplot(data, aes(x=nestCode, y = localStability)) + geom_bar(stat="identity") +
              theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
              facet_wrap(c("wGcode", "btwGuildC"), nrow = length(unique(data$wGcode)), ncol = length(unique(data$btwGuildC)), scales = "free", #scales = "fixed", # 
                         shrink = TRUE, labeller = "label_both", as.table = TRUE,
                         switch = NULL, drop = TRUE, dir = "h", strip.position = "top")
            
            
            # grid.arrange(
            #   p1,
            #   nrow = 1,
            #   top = paste( intTypesCode$name[intTypesCode$value==i], wGuildCCode$name[wGuildCCode$value==j], 'within guild competition', sep=' '),
            #   newpage = T
            # )
            # grid.arrange(
            #   p1,
            #   nrow = 1,
            #   top = paste( i, 'species with', j , 'intra-guild competition', sep=' '),#paste( i, 'species', wGuildCCode$name[wGuildCCode$value==j], 'within guild competition', sep=' '),
            #   newpage = T
            # )
            # grid.arrange(
            #   p1,
            #   nrow = 1,
            #   top = paste( i, 'species with predator-prey inter-guild interactions', sep=' '),#paste( i, 'species', wGuildCCode$name[wGuildCCode$value==j], 'within guild competition', sep=' '),
            #   newpage = T
            # )
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
      
      # Basic style commented away
      {
      #   # initiate plot
      #   setwd(save.dir)
      #   # windows(width=8.3, height=11.7)
      #   pdf(file='LocalStabilityBars_NestProj_basicplot.pdf', width=8.3, height=11.7)
      #   
      #   # Get number of levels for main variables
      #   numComplexityLevels <- c(length(unique(res.frame$S)), length(unique(res.frame$btwGuildC)), length(unique(res.frame$nestCode)))
      #   names(numComplexityLevels) <- c("S", "btwGuildC", "nestCode")
      #   numComplexityLevels <-  numComplexityLevels[order( numComplexityLevels, decreasing = F)]
      #   
      #   # Use the number of levels above to structure plot
      #   plotRows <- numComplexityLevels[1]
      #   plotCols <- numComplexityLevels[2]
      #   par(mfrow=c(plotRows,plotCols))
      #   
      #   data <- res.frame %>%
      #     group_by(intTypes, guildStruct, S, wGuildC, btwGuildC) %>%
      #     summarize(localStability=signif(sum(localStability)/n(), digits=2))
      #   
      #   ymaxMat <- data %>%
      #     group_by(intTypes, guildStruct, S) %>%
      #     summarize(value=max(localStability))
      #   
      #   
      #   for ( i in unique(res.frame$intTypes) ){
      #     
      #     for ( j in unique(res.frame$guildStruct) ){
      #       
      #       for ( k in eval(parse(text=paste( 'unique(res.frame$', names(numComplexityLevels)[1], ')', sep=''))) ){
      #         
      #         for ( l in eval(parse(text=paste( 'unique(res.frame$', names(numComplexityLevels)[2], ')', sep='')))  ){
      #           
      #           plotdata <- data[data$intTypes==i & data$guildStruct ==j & eval(parse(text=paste( 'data$', names(numComplexityLevels)[1], '==k', sep=''))) & eval(parse(text=paste( 'data$', names(numComplexityLevels)[2], '==l', sep=''))),]
      #           
      #           ymaxUse <- ymaxMat$value[ymaxMat$intTypes==i & ymaxMat$guildStruct ==j & eval(parse(text=paste( 'ymaxMat$', names(numComplexityLevels)[1], '==k', sep='')))]
      #           
      #           bph <- barplot(plotdata$localStability, ylim=c(0,ymaxUse),las=2, names.arg=unique(plotdata$btwGuildC), 
      #                          xlab = "btwGuildC", ylab='local stability')
      #           title(main=paste( intTypesCode$name[intTypesCode$value==i], guildStructCode$name[guildStructCode==j], 'guild structure', sep=' '),
      #                 sub=paste(names(numComplexityLevels)[1], '=', k, ' and ', names(numComplexityLevels)[2], '=', l, sep='' ))
      #           
      #           text(bph, ymaxUse, paste(signif(100*plotdata$localStability,digit=2), '%',sep=''), pos=1 )
      #           
      #         }
      #       }
      #     }
      #   }
      #   dev.off()
      }
      
      # Local stability - barplots OLD commented away
      {
        # setwd(save.dir)
        # #windows(width=8.3, height=11.7)
        # pdf(file='LocalStabilityBars.pdf', width=8.3, height=11.7)
        # par(mfrow=c(3,4))
        # 
        # num.nwt <- length(unique(res.frame$network))
        # grp.nwt <- c(1,2,5,8)
        # scenNums <- sort(as.vector(kronecker(matrix(1,length(files)/num.nwt,1),t(grp.nwt))+matrix(seq(0,length(files)-1,num.nwt),ncol=length(grp.nwt), nrow=length(files)/num.nwt)))
        # 
        # for (i in scenNums  ){ # (i in seq(1,nrow(res.frame)-1,3)){ # 
        #   
        #   if ( is.element( i, seq(1,nrow(res.frame)-1,num.nwt) ) ){
        #     
        #     ymax        <- max(res.frame$propStable[i:(i+num.nwt-1)], na.rm = T)
        #     data        <- res.frame$propStable[i]
        #     names(data) <- res.frame$network[i]
        #     
        #   }else{
        #     
        #     data        <- res.frame$propStable[seq(i,i+2,1)]
        #     names(data) <- res.frame$network[seq(i,i+2,1)]
        #     
        #   }
        #   
        #   bph <- barplot(data, ylim=c(0,ymax),las=2)
        #   title(paste(substr(res.frame$intType[i],1,4),
        #               substr(res.frame$intSym[i],1,3),
        #               res.frame$compType[i], sep=' ') )
        #   text(bph, ymax, paste(signif(100*data,digit=2), '%',sep=''), pos=1 )
        # }
        # dev.off()
      }
       
      ### REACTIVITY - BOX PLOTS commented away
      {
        # setwd(save.dir)
        # #windows(width=8.3, height=11.7)
        # pdf(file='ReactivityBoxes.pdf', width=8.3, height=11.7)
        # par(mfrow=c(3,2))
        # 
        # setwd(data.dir)
        # for (i in 1:length(unique(res.frame$intType))){ # 
        #   
        #   v1 <- unique(res.frame$intType)[i]
        #   
        #   for (j in 1:length(unique(res.frame$compType))){
        #     
        #     v2 <- unique(res.frame$compType)[j]
        #     
        #     for (k in 1:length(unique(res.frame$intSym))){
        #       
        #       v3 <- unique(res.frame$intSym)[k]
        #       
        #       data <- as.data.frame(matrix(NA, ncol=2, nrow=data.nrow*length(unique(res.frame$network)) ))
        #       colnames(data) <- c('network', 'reactivity')
        #       
        #       for (l in 1:length(unique(res.frame$network))){
        #         
        #         v4 <- unique(res.frame$network)[l]
        #         
        #         file <- paste("ToyResMat_IntType_", v1, "_IntSym_", v3, "_CompType_", v2 ,"_Network_", v4, ".mat", sep="" )
        #         data$reactivity[((l-1)*data.nrow+1):(l*data.nrow)] <- readMat(file)$resMat[,2]
        #         data$network[((l-1)*data.nrow+1):(l*data.nrow)] <- v4
        #       }
        #       
        #       bxph <- boxplot(reactivity ~ network, data, las=2 )
        #       title(paste(substr(v1,1,4),
        #                   substr(v3,1,3),
        #                   v2, sep=' ') )
        #     }
        #   }
        # }
        # setwd(save.dir)
        # dev.off()
      }
      
      
      ### RESILIENCE - BOX PLOTS commented away
      {
        # setwd(save.dir)
        # #windows(width=8.3, height=11.7)
        # pdf(file='ResilienceBoxes.pdf', width=8.3, height=11.7)
        # par(mfrow=c(3,2))
        # 
        # setwd(data.dir)
        # for (i in 1:length(unique(res.frame$intType))){ # 
        #   
        #   v1 <- unique(res.frame$intType)[i]
        #   
        #   for (j in 1:length(unique(res.frame$compType))){
        #     
        #     v2 <- unique(res.frame$compType)[j]
        #     
        #     for (k in 1:length(unique(res.frame$intSym))){
        #       
        #       v3 <- unique(res.frame$intSym)[k]
        #       
        #       data <- as.data.frame(matrix(NA, ncol=2, nrow=data.nrow*length(unique(res.frame$network)) ))
        #       colnames(data) <- c('network', 'resilience')
        #       
        #       for (l in 1:length(unique(res.frame$network))){
        #         
        #         v4 <- unique(res.frame$network)[l]
        #         
        #         file <- paste("ToyResMat_IntType_", v1, "_IntSym_", v3, "_CompType_", v2 ,"_Network_", v4, ".mat", sep="" )
        #         data$reactivity[((l-1)*data.nrow+1):(l*data.nrow)] <- readMat(file)$resMat[,3]
        #         data$network[((l-1)*data.nrow+1):(l*data.nrow)] <- v4
        #       }
        #       
        #       bxph <- boxplot(reactivity ~ network, data, las=2 )
        #       title(paste(substr(v1,1,4),
        #                   substr(v3,1,3),
        #                   v2, sep=' ') )
        #     }
        #   }
        # }
        # setwd(save.dir)
        # dev.off()
      }
  
    }
  }
  
  ### Reactivity and resilience histograms ### Commented away.
  {
    #   
    #   
    #   for (i in 1:3){ # Loop over comptypes
    #     
    #     #windows(width=8.3, height=11.7)
    #     pdf(file=paste('Reactivity_Comptype_', res.frame$comptype[9*(i-1)+1], '.pdf', sep=''), width=8.3, height=11.7)
    #     par(mfrow=c(3,3))
    #     
    #     for (j in 1:9){ # Loop over networks
    #       
    #       r <- 9*(i-1)+j
    #       
    #       file <- paste("ToyResMat_Comptype_", res.frame$comptype[r], "_Network_", res.frame$network[r], ".mat", sep="" )
    #       data <- readMat(file)$resMat
    #       reactivity <- data[,2]
    #       
    #       if  ( sum(data[,1]==1) == 0 ){
    #         plot.new()
    #       }else{
    #         hist(reactivity, main=paste(res.frame$comptype[r], ' and ', res.frame$network[r], sep='') )
    #       }
    #       
    #     }
    #     dev.off()
    #   }
    #   
    #   
    #   # Resilience
    #   
    #   
    #   for (i in 1:3){ # Loop pver comptypes
    #     
    #     #windows(width=8.3, height=11.7)
    #     pdf(file=paste('Resilience_Comptype_', res.frame$comptype[9*(i-1)+1], '.pdf', sep=''), width=8.3, height=11.7)
    #     par(mfrow=c(3,3))
    #     
    #     for (j in 1:9){ # Loop over networks
    #       
    #       r <- 9*(i-1)+j
    #       
    #       file <- paste("ToyResMat_Comptype_", res.frame$comptype[r], "_Network_", res.frame$network[r], ".mat", sep="" )
    #       data <- readMat(file)$resMat
    #       resilience <- data[,3]
    #       
    #       if  ( sum(data[,1]==1) == 0 ){
    #         plot.new()
    #       }else{
    #         hist(resilience, main=paste(res.frame$comptype[r], ' and ', res.frame$network[r], sep='') )
    #       }
    #       
    #     }
    #     dev.off()
    #   }
  }
}
  
  # Plot network results only
  if ( analyseBtwGuildStructure == 1 ){
    
    ### Plots for guild structure evaluation ###
    if ( params.data$project == 1 ){
      
      ### Plots for nestedness exploration ###  
    } else if ( params.data$project == 2 ){
      
      # Get scenario data
      scens    <- sort(unique(bipStruct.frame$scenID))
      scenData <- do.call("rbind", replicate(reps, params.data[params.data$scenID==scens[1],], simplify = FALSE))
      for (i in scens[2:length(scens)]){
        tmp <- do.call("rbind", replicate(reps, params.data[params.data$scenID==i,], simplify = FALSE))
        scenData <- rbind(scenData, tmp)
      }
      
      # Just bug precaution. If scenarios don't match up, the analysis will ground to a halt here.
      if (!(all(bipStruct.frame$scenID == scenData$scenID))){
        rm(scenData)
      }
      
      # Combine network results with scenario info
      bipStruct.frame <- cbind(bipStruct.frame, scenData[, !(names(scenData) %in% c("scenID"))])
      rm(scenData)
      
      #ggplot
      {
        # initiate plot
        setwd(save.dir)
        pdf(file= 'NetworkMetricBars_NestProj_ggplot_fixedY.pdf', width=8.3, height=11.7)
        
        for ( i in sort(unique(res.frame$S)) ){
           
          for ( j in 1:length(structMetrics) ){
             
            data <- bipStruct.frame[bipStruct.frame$S==i ,]  
           
            eval(parse(text=paste('p',j,  '<- ggplot(data, aes(x=nestCode, y =', structMetrics[j] ,')) + geom_boxplot() +
              theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
              facet_wrap(c("btwGuildC"), nrow = 1, ncol = length(unique(data$btwGuildC)), scales = "fixed",
              shrink = TRUE, labeller = "label_both", as.table = TRUE,
                         switch = NULL, drop = TRUE, dir = "h", strip.position = "top")', sep='')))
            
          } # Loop over j
          
          
          grid.arrange(
            p1,p2,p3,p4,p5,p6,p7,p8,
            nrow = length(handleNames),
            top = paste( 'All network metrics for', i, 'species networks' , sep=' '),
            newpage = T
          )
          
        } # Loop over i
        dev.off()
      }
      
    }
 
  }
 
   # Plot stability and network results combined
  if ( analyseStability == 1 & analyseBtwGuildStructure == 1 ){
    
    ### Plots for guild structure evaluation ###
    if (params.data$project == 1){
      
      # Nothing develop for project 1 so far.
      
      ### Plots for nestedness exploration ###  
    } else if (params.data$project == 2){
      
      # Combined local stability and qFcomp plot
      #ggplot
      {
        # initiate plot
        setwd(save.dir)
        # windows(width=8.3, height=11.7)
        pdf(file='LocalStabilityBars_and_qFcompBoxPlots_NestProj_ggplot_fixedY.pdf', width=8.3, height=11.7)
        
        for ( i in sort(unique(res.frame$intTypes)) ){
          
          for ( j in sort(unique(res.frame$S[res.frame$intTypes==i])) ){
            
              data <- res.frame[res.frame$intTypes==i & res.frame$S==j ,]
              p1   <- ggplot(data, aes(x=nestCode, y = localStability)) + geom_bar(stat="identity") +
              theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
              facet_wrap(c("wGcode", "btwGuildC"), nrow = length(unique(data$wGcode)), ncol = length(unique(data$btwGuildC)), scales = "free", #scales = "fixed", # 
                         shrink = TRUE, labeller = "label_both", as.table = TRUE,
                         switch = NULL, drop = TRUE, dir = "h", strip.position = "top")
              
            data <- bipStruct.frame[bipStruct.frame$S==j ,]  
            p2   <- ggplot(data, aes(x=nestCode, y = qFcomp)) + geom_boxplot() +
                              theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                                  facet_wrap(c("btwGuildC"), nrow = 1, ncol = length(unique(data$btwGuildC)), scales = "fixed",
                                  shrink = TRUE, labeller = "label_both", as.table = TRUE,
                                  switch = NULL, drop = TRUE, dir = "h", strip.position = "top")
            
            grid.arrange(
              p1, p2,
              nrow = 5,
              layout_matrix = rbind(c(1),c(1),c(1),c(1),c(2)),
              top = paste( i , 'with', j, 'species' , sep=' '),
              newpage = T
            )
            
            
          } # Loop over j
        } # Loop over i
        dev.off()
      }
      
      # Plotting regression style plots, with qFcomp as predictor of local stability
      
      # This version makes no sense - since I don't have the actual nw metrics for all scenarios!
      {
        # Get stability data
        scens    <- sort(unique(res.frame$scenID))
        stabData <- res.frame[res.frame$scenID==scens[1], ][1:reps,]
        nw.tmp   <- bipStruct.frame[ bipStruct.frame$S == res.frame$S[res.frame$scenID==scens[1]][1] &  
                                       bipStruct.frame$btwGuildC == res.frame$btwGuildC[res.frame$scenID==scens[1]][1] &
                                       bipStruct.frame$nestCode == res.frame$nestCode[res.frame$scenID==scens[1]][1], 
                                     names(bipStruct.frame) %in% c( structMetrics)]
        stabData <- cbind(stabData, nw.tmp)
        for (i in scens[2:length(scens)]){
          stab.tmp <- res.frame[res.frame$scenID==i, ][1:reps,]
          nw.tmp   <- bipStruct.frame[ bipStruct.frame$S == res.frame$S[res.frame$scenID==i][1] &  
                                         bipStruct.frame$btwGuildC == res.frame$btwGuildC[res.frame$scenID==i][1] &
                                         bipStruct.frame$nestCode == res.frame$nestCode[res.frame$scenID==i][1], 
                                       names(bipStruct.frame) %in% c( structMetrics)]
          stab.tmp <- cbind(stab.tmp, nw.tmp)
          stabData <- rbind(stabData, stab.tmp)
        }
        
        # # Just bug precaution. If scenarios don't match up, the analysis will ground to a halt here.
        # if (!(all(bipStruct.frame$scenID == stabData$scenID))){
        #   rm(stabData)
        # }

        
        # Plot local stability as a function of qFcomp
        
      }
      
      # This version is also weird.. I actually need the nw metrics for all scenarios for this anlysis 
      # to make sense. But for now this is what makes most sense...at least the local stability values I plot 
      # are for those actual networks and their metrics...
       {
        # Get stability data
        scens    <- sort(unique(bipStruct.frame$scenID))
        stabData <- res.frame[res.frame$scenID==scens[1], names(res.frame) %in% c("scenID","localStability", "reactivity", "resilience")][1:reps,]
        for (i in scens[2:length(scens)]){
          tmp <- res.frame[res.frame$scenID==i, names(res.frame) %in% c("scenID","localStability", "reactivity", "resilience")][1:reps,]
          stabData <- rbind(stabData, tmp)
        }
        
        # Just bug precaution. If scenarios don't match up, the analysis will ground to a halt here.
        if (!(all(bipStruct.frame$scenID == stabData$scenID))){
          rm(stabData)
        }
        
        # Combine network results with scenario info
        bipStruct.frame <- cbind(bipStruct.frame, stabData[, !(names(stabData) %in% c("scenID"))])
        rm(stabData)
        
        # Plot local stability as a function of qFcomp
        
        # initiate plot
        setwd(save.dir)
        # windows(width=8.3, height=11.7)
        pdf(file='LocalStability_vs_qFcomp_NestProj_ggplot.pdf', width=8.3, height=11.7)
        
        data <- bipStruct.frame
        
        ggplot(data, aes(x=qFcomp, y = localStability)) + geom_point() + #geom_jitter() #geom_bin2d() #geom_count 
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          facet_wrap(c("S","btwGuildC"), nrow = length(unique(data$S)), ncol = 3, scales = "fixed",
                     shrink = TRUE, labeller = "label_both", as.table = TRUE,
                     switch = NULL, drop = TRUE, dir = "h", strip.position = "top")
        
        dev.off()
        
      }
      
    }
  }
}
