##############################################
### analyze jacobian connectance on cluster ###
##############################################
#
#   @DESCRIPTION
#   Analysis code for project SSMNW, guildProj.
#   For a given scenario, analyze jacobian connectance
#   for all replicates.
#
#   @AUTHOR
#   Alva Curtsdotter, postdoc @ Emory University, 
#   2018-08-01
#   
##############################################

#### Setup ####
{
  
  # clearing the workspace. 
  rm(list=ls(all=TRUE))
  
  # Read in libraries:
  require(R.matlab)   # To read .mat files and such
  
  # # Reading command line input (when on cluster)
  args=(commandArgs())
  print(args)
  eval(parse(text = args[[10]])) # gets the variable scid, i.e. scenario id.
  eval(parse(text = args[[11]])) # gets the variable timestamp
  #scid<-1 #if testrunning locally
  
  # Data paths and such
  data.dir <- paste('/proj/ecoservice/users/x_alvcu/SSMNW/Results/', timestamp,'/', sep='')  
  
  #data.dir <- '/proj/ecoservice/users/x_alvcu/SSMNW/Results/2018Oct11_1746/'  
  #data.dir <- 'C:/Users/icurtsd/Documents/Project Atlanta/Project Structure and stability of mutualistic networks/data/clusterdata/analysis/2018Oct11_1746/'  
  
  save.dir <- paste('/proj/ecoservice/users/x_alvcu/SSMNW/analysis/', timestamp, '/networkStructure/', sep='')  
  #save.dir <- paste('/proj/ecoservice/users/x_alvcu/SSMNW/analysis/', substr(data.dir, regexpr('2018', data.dir), nchar(data.dir)), 'jacobianConnectance/', sep='')
  #save.dir <- data.dir
  if ( !(dir.exists(save.dir)) ){ 
    dir.create(save.dir)
  } 
  
  # Definitions
  structMetrics  <- c('Cj.tot', 'Cj.off', 'propNegTot', 'propNegOff', 'propPosOff')
  corrMetrics    <- c('pcorr', 'pcorr.abs', 'scorr', 'scorr.abs' )
  
  
}

### Load parameter matrix ###
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
}

### Load and Analyse networks ###
{
  setwd(data.dir)
  
  # Load the 3D matrix of Jacobian matrices
  file        <- dir(pattern=paste('allJacobians_scenID_', scid, '.mat', sep=''))
  filePresent <- length(file)>0
  
  if ( filePresent ){
    
    allJ <- readMat(file)$repMat
    if ( length(allJ)==0 ) {
      allJ <- readMat(file)$x
    }
    reps <- unique(params.data$saveJac)
    
    if ( dim(allJ)[3] != reps  ) print('Dimension error of result data!')
    
    # Create the data frame to store the summarized results for all replicates!
    c.frame <- as.data.frame(matrix(NA, ncol=1+length(structMetrics), nrow=params.data$saveJac[scid]))
    names(c.frame) <- c('scenID', structMetrics)
    str(c.frame)
    
    e.frame <- c.frame
    
    corr.frame <- as.data.frame(matrix(NA, ncol=1+length(corrMetrics), nrow=params.data$saveJac[scid]))
    names(corr.frame) <- c('scenID', corrMetrics)
    str(corr.frame)
    
    for (r in 1:params.data$saveJac[scid]){
      
      #define jacobian
      J <- allJ[,,r]
      J.off <- J * !diag(as.vector(matrix(1,ncol=1,nrow=nrow(J))))
      S <- dim(J)[1]
      s <- S/2
      
      #calculate connectance
      c.frame$scenID[r]  <- scid
      c.frame$Cj.tot[r]  <- sum(J!=0)/S^2
      c.frame$Cj.off[r]  <- sum(J.off!=0)/S^2
      c.frame$propNegTot[r] <- sum(J<0)/sum(J!=0)
      c.frame$propNegOff[r] <- sum(J.off<0)/sum(J.off!=0) 
      c.frame$propPosOff[r] <- sum(J.off>0)/sum(J.off!=0)
      
      # expected values
      e.frame$scenID[r]  <- scid
      e.frame$Cj.tot[r]  <- params.data$btwGuildC[scid]/2 + 1/params.data$S[scid]
      e.frame$Cj.off[r]  <- params.data$btwGuildC[scid]/2
      if ( params.data$intTypes[scid] == 2 ){
        if ( params.data$wgcIgnore[scid] == 0 & params.data$wGuildC[scid] == 0){
          e.frame$propNegOff[r] <- 0
          e.frame$propPosOff[r] <- 1
        }else if ( params.data$wgcIgnore[scid] == 1 ){
          e.frame$propNegOff[r] <- 0.5
          e.frame$propPosOff[r] <- 0.5
        }
      } else if ( params.data$intTypes[scid] == 3 ){
        if ( params.data$wgcIgnore[scid] == 0 & params.data$wGuildC[scid] == 0){
          e.frame$propNegTot[r] <-  (0.5*e.frame$Cj.off[r]*params.data$S[scid]^2 + params.data$S[scid])/(e.frame$Cj.off[r]*params.data$S[scid]^2 + params.data$S[scid])
          e.frame$propNegOff[r] <- 0.5
          e.frame$propPosOff[r] <- 0.5
        }else if ( params.data$wgcIgnore[scid] == 1 ){
          e.frame$propNegTot[r] <-  (0.75*e.frame$Cj.off[r]*params.data$S[scid]^2 + params.data$S[scid])/(e.frame$Cj.off[r]*params.data$S[scid]^2 + params.data$S[scid])
          e.frame$propNegOff[r] <- 0.75
          e.frame$propPosOff[r] <- 0.25
        }
      }
      
      # Calculate correlations
      
      corr.frame$scenID[r]  <- scid
      
      if ( params.data$guildStruct[scid] == 0 ){
        
        if ( params.data$intTypes[scid] == 1 ) {
          
          corr.frame$scorr[r]     <- NA
          corr.frame$scorr.abs[r] <- NA
          corr.frame$pcorr[r]     <- NA
          corr.frame$pcorr.abs[r] <- NA
          
        } else { # else of if ( params.data$intTypes[scid] == 1 ) {
          
          # Get pairwise bg interactions, in a way that works for guildstructured and random webs.
          triu <- J.off[upper.tri(J.off)]
          tJ.off <- t(J.off)
          tril <- tJ.off[upper.tri(tJ.off)]
          
          vals <- cbind(triu, tril)
          vals <- vals[!(vals[,1]<0 & vals[,2]<0),] # note that this only works if the intType is not competition. It's not possible to disentangle in a scrambled web which pairwise interactions used to be wg and which used to be bg.
          vals <- vals[!(vals[,1]==0 & vals[,2]==0),]
          
          if ( params.data$intTypes[scid] == 2 ) {
            
            corr.frame$scorr[r]     <- NA
            corr.frame$scorr.abs[r] <- NA
            
          } else if ( params.data$intTypes[scid] == 3 ) {
            
            tmp1 <- vals[,1]
            tmp2 <- vals[,2]
            tmp1[tmp1>0] <- vals[ vals[,2]<0 ,2] # this only works for pred-prey
            tmp2[tmp2<0] <- vals[ vals[,1]>0 ,1]
            vals <- cbind(tmp1, tmp2)
            
            corr.frame$scorr[r]     <- cor(vals[,1], vals[,2], method='spearman')
            corr.frame$scorr.abs[r] <- cor(abs(vals[,1]), abs(vals[,2]), method='spearman')
            
          }
          
          corr.frame$pcorr[r]     <- cor(vals[,1], vals[,2], method='pearson')
          corr.frame$pcorr.abs[r] <- cor(abs(vals[,1]), abs(vals[,2]), method='pearson')
          
        } #end of if ( params.data$intTypes[scid] == 1 ) {
        
        
      } else { #else of if ( params.data$guildStruct[scid] == 0 ){
        
        APchunk <- J.off[(1+s):S, 1:s]
        APvals  <- APchunk[APchunk!=0]
        PAchunk <- t(J.off[1:s, (1+s):S])
        PAvals  <- PAchunk[PAchunk!=0]
        vals <- cbind(APvals, PAvals)
        
        corr.frame$scorr[r]     <- cor(vals[,1], vals[,2], method='spearman')
        corr.frame$scorr.abs[r] <- cor(abs(vals[,1]), abs(vals[,2]), method='spearman')
        
        corr.frame$pcorr[r]     <- cor(vals[,1], vals[,2], method='pearson')
        corr.frame$pcorr.abs[r] <- cor(abs(vals[,1]), abs(vals[,2]), method='pearson')
        
      }#end of if ( params.data$guildStruct[scid] == 0 ){
      
    }#end of replicates for-loop
    
    str(corr.frame)
    str(c.frame)
    str(e.frame)
    diff <- c.frame-e.frame
    tolerance <- 4/S^2 # 4 links off [if 1 rounded up link per block]
    if (any( abs(diff[,]) > abs(tolerance), na.rm=T ) ) {
      verdict <- 0
      print('Jacobian not as expected. Investigate!')
    } else {
      verdict <- 1
    }
    
  }#end of if filePresent
  
}

### Save res.frame of network metrics for all replicates ###
{
  setwd(save.dir)
  
  savename <- paste('jacobianCheck_scenario_', scid, '.RData', sep='')
  save(c.frame, e.frame, diff, tolerance, verdict, corr.frame, file=savename)

}  
