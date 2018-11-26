##############################################
### analyze bipartite network structure on cluster ###
##############################################
#
#   @DESCRIPTION
#   Analysis code for project SSMNW.
#   For a given scenario, analyze network structure 
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
  require(bipartite)  # Network structure analysis
  require(vegan)      # Network structure analysis
  
  # find out which type of machine you're running on
  os <- .Platform$OS.type
  
  # Reading arguments, defining paths
  if ( os == 'windows' ){ #i.e. my local pc. only for development. 
    
    # define the scenario id to work with
    scid <- 73
    
    # Data paths and such
    
    main.dir  <- 'C:/Users/icurtsd/Documents/Project Atlanta/Project Structure and stability of mutualistic networks/data/'
    sub.dir   <- 'main_test/' # 'clusterdata/analysis/' # 
    stupid1   <- 'testresults_' # '' # 
    stupid2   <- 'testanalysis_' # '' # 
    timestamp <-'2018Jul11T1636/' #  '2018Oct11_1746/' # '2018Oct04T1505/' # 
    
    data.dir <- paste(main.dir, sub.dir, stupid1, timestamp, sep='')
    save.dir <- paste(main.dir, sub.dir, stupid2, substr(data.dir, regexpr('201', data.dir), nchar(data.dir)), 'networkStructure/',  sep='')
    
  } else if ( os == 'unix'){ #i.e. cluster [or a Mac unfortunately]
    
    # Reading command line input (when on cluster)
    args=(commandArgs())
    print(args) 
    eval(parse(text = args[[10]])) # gets the variable scid, i.e. scenario id.
    eval(parse(text = args[[11]])) # gets the variable timestamp
    
    # Data paths and such
    data.dir <- paste('/proj/ecoservice/users/x_alvcu/SSMNW/Results/', timestamp, '/', sep='')  
    save.dir <- paste('/proj/ecoservice/users/x_alvcu/SSMNW/analysis/', timestamp, '/networkStructure/', sep='')
    
  }
  
  # Create save.dir it it doesn't already exist.
  if ( !(dir.exists(save.dir)) ){ 
    dir.create(save.dir)
  } 
  
  # Some definitions
  checkMetrics  <- c('Cj.tot', 'Cj.off', 'btwgC', 'wgC', 'corrcoef' ) 
  
  structMetrics  <- c('bNestNODF', 'bNestTemp', 'bNestSpecRad', 'qNestNODF', 'qNestWine', 'qNestSpecRad',
                      'bFcomp', 'qFcomp', 'bFcomp2', 'qFcomp2', 'fc.hl', 'fc.ll', 'no.hl', 'no.ll')

  rcSumMetrics   <- c( 'max.rsJ','max.csJ','max.rsBbip','max.csBbip','max.rsQbip','max.csQbip', 
                       'mean.rsJ','mean.csJ','mean.rsBbip','mean.csBbip','mean.rsQbip','mean.csQbip',
                       'std.rsJ','std.csJ','std.rsBbip','std.csBbip','std.rsQbip','std.csQbip',
                       'lm.rsJ','lm.csJ','lm.rsBbip','lm.csBbip','lm.rsQbip','lm.csQbip',
                       'PJ.rsJ','PJ.csJ', 'PJ.tsJ','PJ.rsBbip','PJ.csBbip', 'PJ.tsBbip', 'PJ.rsQbip','PJ.csQbip', 'PJ.tsQbip')
  
  old <- 0 # (as.Date(substring(timestamp, 1, nchar(timestamp)-6), "%Y%b%d")-as.Date("2018-09-01")) < 0
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
  
  if ( old ) {
    
    # Check that number of result files equals the number of replicates.
    files      <- dir(pattern=paste('Jacobian_scenID_', scid, '_replicate_*', sep=''))
    replicates <- unique(params.data$trials)
    if ( length(files) != replicates  ) print('Dimension error of result data!')
    
  } else {
    file <- dir(pattern=paste('allJacobians_scenID_', scid, '.mat', sep=''))
    allJ <- eval(parse(text=paste('readMat(file)$', names(readMat(file)), sep='')))
    
    if (length(intersect(names(params.data), 'saveJac')) > 0){
      replicates <- min(unique(params.data$saveJac), dim(allJ)[3])
    } else {
      replicates <- min(unique(params.data$trials), dim(allJ)[3])
    }
    if ( dim(allJ)[3] != replicates  ) print('Potential Dimension error of result data!')
  }
  
  # Create the data frame to store the summarized results for all replicates!
  c.frame <- as.data.frame(matrix(NA, ncol=2+length(checkMetrics), nrow=replicates))
  names(c.frame) <- c('scenID', 'repID', checkMetrics)
  str(c.frame)
  
  res.frame <- as.data.frame(matrix(NA, ncol=2+length(structMetrics), nrow=replicates))
  names(res.frame) <- c('scenID', 'repID', structMetrics)
  str(res.frame)
  
  rcSum.frame <- as.data.frame(matrix(NA, ncol=2+length(rcSumMetrics), nrow=replicates))
  names(rcSum.frame) <- c('scenID', 'repID', rcSumMetrics)
  str(rcSum.frame)
  
  
  for (r in 1:replicates){
    
    # Get Jacobian for this replicate
    if ( old ) {
      # read in data file
      file <- paste('Jacobian_scenID_', scid, '_replicate_', r, '.csv', sep='')
      J  <- as.matrix(read.table(file, sep=","))
    } else {
      J <- allJ[,,r]
    }
    
    # Extract (one of) the bipartite chunks (the APchunk, i.e. plant effect on animals)  
    bip  <- J[(1+nrow(J)/2):nrow(J),1:(ncol(J)/2)]
    bip2 <- J[1:(nrow(J)/2), (1+ncol(J)/2):ncol(J)]
    ppchunk <- J[1:(nrow(J)/2), 1:(ncol(J)/2) ]
    # str(bip)
    
    # Binary bipartite structure
    bbip  <- bip
    bbip[bbip!=0]   <- 1
    bbip2 <- bip2
    bbip2[bbip2!=0] <- 1
    
    # Quantitative bipartite structure ABSOLUTE VALUES!!!
    qbip  <- abs(bip)
    qbip2 <- abs(bip2)
    
    # Save scenario id
    c.frame$scenID[r]      <- scid
    res.frame$scenID[r]    <- scid
    rcSum.frame$scenID[r]  <- scid
    
    # Save replicate id
    c.frame$repID[r]      <- r
    res.frame$repID[r]    <- r
    rcSum.frame$repID[r]  <- r
    
    # Do some checks of the Jacobian structure
    J.off <- J * !diag(as.vector(matrix(1,ncol=1,nrow=nrow(J))))
    ppchunk.off <- ppchunk * !diag(as.vector(matrix(1,ncol=1,nrow=nrow(ppchunk))))
    S <- dim(J)[1]
    s <- S/2
    
    c.frame$Cj.tot[r] <- sum(J!=0)/S^2           #calculate total connectance
    c.frame$Cj.off[r] <- sum(J.off!=0)/S^2       #calculate off-diagonal connectance
    c.frame$btwgC[r]  <- sum(bbip)/s^2
    c.frame$wgC[r]    <- sum(ppchunk.off!=0)/s^2
    
    corr.x <- bip[bip!=0]                      #calculate correlation btw pairwise interacion strenghts
    corr.y <- t(bip2)
    corr.y <- bip2[bip2!=0]
    corrcoef <- signif(cor(cbind(corr.x,corr.y))[2],2)
    c.frame$corrcoef[r] <- corrcoef
    
    # Binary Nestedness metrics
    res.frame$bNestTemp[r] <- nested(bbip, method='binmatnest2') 
    res.frame$bNestNODF[r] <- nested(bbip, method='NODF2')
    
    zeroblock    <- matrix(0,nrow=nrow(bbip),ncol=ncol(bbip))
    A            <- cbind(rbind(zeroblock,bbip), rbind(t(bbip),zeroblock))
    res.frame$bNestSpecRad[r] <- max(eigen(A)$values)
    
    # Quantitative nestedness metrics
    res.frame$qNestNODF[r] <- nested(qbip, method='weighted NODF')
    res.frame$qNestWine[r] <- wine(qbip, nreps = 100)$wine
    
    zeroblock    <- matrix(0,nrow=nrow(qbip),ncol=ncol(qbip))
    A            <- cbind(rbind(zeroblock,qbip), rbind(t(qbip),zeroblock))
    res.frame$qNestSpecRad[r] <- max(eigen(A)$values)
    
    # Binary complementarity
    ppm <- bbip
    distance = vegdist(t(ppm), method="euclidean", bin=T)
    distance[is.nan(distance)] <- 0
    tree <- hclust(distance)
    res.frame$bFcomp[r] <- treeheight(tree)
    
    ppm <- bbip2
    distance = vegdist(t(ppm), method="euclidean", bin=T)
    distance[is.nan(distance)] <- 0
    tree <- hclust(distance)
    res.frame$bFcomp2[r] <- treeheight(tree)
    
    # Quantitative complementarity
    ppm <- qbip
    distance = vegdist(t(ppm), method="euclidean", bin=F)
    distance[is.nan(distance)] <- 0
    tree <- hclust(distance)
    res.frame$qFcomp[r] <- treeheight(tree)
    
    ppm <- qbip2
    distance = vegdist(t(ppm), method="euclidean", bin=F)
    distance[is.nan(distance)] <- 0
    tree <- hclust(distance)
    res.frame$qFcomp2[r] <- treeheight(tree)
    
    # Use bipartites functions for complementarity calculations
    bpc <- networklevel(bbip, index=c("functional complementarity", "niche overlap"), level="both")
    
    res.frame$fc.hl[r] <- bpc["functional.complementarity.HL"]
    res.frame$fc.ll[r] <- bpc["functional.complementarity.LL"]
    res.frame$no.hl[r] <- bpc["niche.overlap.HL"]
    res.frame$no.ll[r] <- bpc["niche.overlap.LL"]
    
    # Calculate Row and Column Sums
    rsJ <- rowSums(J)
    csJ <- colSums(J)
    rsBbip <- rowSums(bbip)
    csBbip <- colSums(bbip)
    rsQbip <- rowSums(qbip)
    csQbip <- colSums(qbip)
    
    # Max, mean, std, and linear trend in row and column sums
    rcSum.frame$max.rsJ[r]  <- max(rsJ)
    rcSum.frame$mean.rsJ[r] <- mean(rsJ)
    rcSum.frame$std.rsJ[r]  <- sd(rsJ)
    tempdata <- as.data.frame(cbind(x=1:length(rsJ),y=rsJ))
    rcSum.frame$lm.rsJ[r]   <- lm(tempdata$y ~ tempdata$x)$coefficients[2]
    
    rcSum.frame$max.csJ[r]  <- max(csJ)
    rcSum.frame$mean.csJ[r] <- mean(csJ)
    rcSum.frame$std.csJ[r]  <- sd(csJ)
    tempdata <- as.data.frame(cbind(x=1:length(csJ),y=csJ))
    rcSum.frame$lm.csJ[r]   <- lm(tempdata$y ~ tempdata$x)$coefficients[2]
    
    rcSum.frame$max.rsBbip[r]  <- max(rsBbip)
    rcSum.frame$mean.rsBbip[r] <- mean(rsBbip)
    rcSum.frame$std.rsBbip[r]  <- sd(rsBbip)
    tempdata <- as.data.frame(cbind(x=1:length(rsBbip),y=rsBbip))
    rcSum.frame$lm.rsBbip[r]   <- lm(tempdata$y ~ tempdata$x)$coefficients[2]
    
    rcSum.frame$max.csBbip[r]  <- max(csBbip)
    rcSum.frame$mean.csBbip[r] <- mean(csBbip)
    rcSum.frame$std.csBbip[r]  <- sd(csBbip)
    tempdata <- as.data.frame(cbind(x=1:length(csBbip),y=csBbip))
    rcSum.frame$lm.csBbip[r]   <- lm(tempdata$y ~ tempdata$x)$coefficients[2]
    
    rcSum.frame$max.rsQbip[r]  <- max(rsQbip)
    rcSum.frame$mean.rsQbip[r] <- mean(rsQbip)
    rcSum.frame$std.rsQbip[r]  <- sd(rsQbip)
    tempdata <- as.data.frame(cbind(x=1:length(rsQbip),y=rsQbip))
    rcSum.frame$lm.rsQbip[r]   <- lm(tempdata$y ~ tempdata$x)$coefficients[2]
    
    rcSum.frame$max.csQbip[r]  <- max(csQbip)
    rcSum.frame$mean.csQbip[r] <- mean(csQbip)
    rcSum.frame$std.csQbip[r]  <- sd(csQbip)
    tempdata <- as.data.frame(cbind(x=1:length(csQbip),y=csQbip))
    rcSum.frame$lm.csQbip[r]   <- lm(tempdata$y ~ tempdata$x)$coefficients[2]
    
    # Evenness of row and column sums (PJ stands for Pielou's J)
    rcSum.frame$PJ.rsJ[r] <- diversity(rsJ+abs(min(rsJ)-1))/log(length(rsJ))
    rcSum.frame$PJ.csJ[r] <- diversity(csJ+abs(min(csJ)-1))/log(length(csJ))
    rcSum.frame$PJ.tsJ[r] <- diversity(c(csJ+abs(min(csJ)-1), rsJ+abs(min(rsJ)-1)) )/log(length(rsJ)+length(csJ))
    
    rcSum.frame$PJ.rsBbip[r] <- diversity(rsBbip+abs(min(rsBbip)-1))/log(length(rsBbip))
    rcSum.frame$PJ.csBbip[r] <- diversity(csBbip+abs(min(csBbip)-1))/log(length(csBbip))
    rcSum.frame$PJ.tsBbip[r] <- diversity(c(csBbip+abs(min(csBbip)-1), rsBbip+abs(min(rsBbip)-1)) )/log(length(rsBbip)+length(csBbip))
    
    rcSum.frame$PJ.rsQbip[r] <- diversity(rsQbip+abs(min(rsQbip)-1))/log(length(rsQbip))
    rcSum.frame$PJ.csQbip[r] <- diversity(csQbip+abs(min(csQbip)-1))/log(length(csQbip))
    rcSum.frame$PJ.tsQbip[r] <- diversity(c(csQbip+abs(min(csQbip)-1), rsQbip+abs(min(rsQbip)-1)) )/log(length(rsQbip)+length(csQbip))
    
    
  }#end of replicates for-loop
  
  str(c.frame)
  str(res.frame)
  str(rcSum.frame)
  
}

### Save res.frame of network metrics for all replicates ###
{
  setwd(save.dir)
  
  savename <- paste('checkMetrics_scenario_', scid, '.csv', sep='')
  write.table(c.frame, savename)
  
  savename <- paste('bipNetworkMetrics_scenario_', scid, '.csv', sep='')
  write.table(res.frame, savename)
  
  savename <- paste('rowAndColSumMetrics_scenario_', scid, '.csv', sep='')
  write.table(rcSum.frame, savename)
}  
