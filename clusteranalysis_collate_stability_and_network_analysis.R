##############################################
### collate stability and network analysis results on cluster ###
##############################################
#
#   @DESCRIPTION
#   Analysis code for project SSMNW.
#   Collect and collate all results for 
#   stability and network metrics from all 
#   scenarios and replicates, into one humungous
#   data frame (ncol = ncol(params.data)+ 
#   ncol(stabilityRestults))+ncol(networkMterics), and
#   nrow=scenarios*replicates.
#   Subset the data frame if and as wanted.                
#
#   @AUTHOR
#   Alva Curtsdotter, postdoc @ Emory University, 
#   2018-08-02
#   
##############################################

#### Setup ####
{
  
  # clearing the workspace. 
  rm(list=ls(all=TRUE))
  
  # Read in libraries:
  require('R.matlab') # To read .mat files and such
  
  # find out which type of machine you're running on
  os <- .Platform$OS.type
  
  # Reading arguments, defining paths
  if ( os == 'windows' ){ #i.e. my local pc. only for development. 
    
    # define the simulation folder to work with
    timestamp <- '2018Oct20_0125'
    # Collating options
    collate.stab <- 1
    collate.nwm  <- 0
    
    # Data paths and such
    main.dir  <- 'C:/Users/icurtsd/Documents/Project Atlanta/Project Structure and stability of mutualistic networks/data/'
    sub.dir   <- 'clusterdata/analysis/' # 'main_test/' # 
    stupid1   <- '' #'testresults_' #  
    
    data.dir.1 <- paste(main.dir, sub.dir, stupid1, timestamp, sep='')
    save.dir <- data.dir.1 
    
  } else if ( os == 'unix'){ #i.e. cluster [or a Mac unfortunately]
    
    # Reading command line input (when on cluster)
    args=(commandArgs())
    print(args) 
    eval(parse(text = args[[10]])) # gets the variable collate.stab # make this a for-loop if u want greater flexibility
    eval(parse(text = args[[11]])) # gets the variable collate.nwm
    eval(parse(text = args[[12]])) # gets the variable timestamp

    # Data paths and such
    data.dir.1 <- paste('/proj/ecoservice/users/x_alvcu/SSMNW/Results/', timestamp,'/', sep='')  
    save.dir <- paste('/proj/ecoservice/users/x_alvcu/SSMNW/analysis/', timestamp,'/', sep='')  
    
    if ( collate.nwm ){
      data.dir.2 <- paste('/proj/ecoservice/users/x_alvcu/SSMNW/analysis/', timestamp, '/networkStructure/', sep='')
    }
    
  }
  print(data.dir.1)
  print(save.dir)
  
  if ( !(dir.exists(save.dir)) ){ 
    dir.create(save.dir)
  } 
  
  resultVarNames <- c('localStability', 'reactivity', 'resilience')
  
 
  # Subsetting options
  do.subset <- 0 # to subset or not to subset? that is the question.
  
  # Columns to drop
  # colDrop <- c('project', 'numGuilds', 'guildStruct', 'wgcIgnore', 'wGuildCompType',
  #              'numIntTypes', 'intTypeSymmetry', 'JacobianSum', 'saveJac')
  colDrop <- c('')
  
  # Scenarios, i.e. rows, to drop:
  # Hardcoded in the subsetting section
}

### Load parameter matrix ###
{
  setwd(data.dir.1)
  
  # Read parameterStruct for the simulation
  file1  <- dir(pattern = 'paramStruct')
  params <- readMat(file1)
  
  # Get the paramStruct into shape.
  params.data <- lapply(params$params[2], unlist, use.names=FALSE)
  params.data <- as.data.frame(params.data)
  names(params.data) <- c('scenID', lapply(params$params[1], unlist, use.names=FALSE)[[1]])
  #str(params.data)
  
  if (params.data$project==2){
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
  }
  
  # params.data$wGcode = NA
  # params.data$wGcode[params.data$wGuildC == 0 ]                                     = 'no'
  # params.data$wGcode[params.data$wGuildC == 1 & params.data$wGuildWeightStd == 0.25]  = 'intraGuildC=1 and weak'
  # params.data$wGcode[params.data$wGuildC == 1 & params.data$wGuildWeightStd == 0.50]  = 'intraGuildC=1 and strong'
  # params.data$wGcode[params.data$wgcIgnore == 1 & params.data$wGuildWeightStd == 0.25]  = 'intraGuildC=btwGuildC and weak'
  # params.data$wGcode[params.data$wgcIgnore == 1 & params.data$wGuildWeightStd == 0.50]  = 'intraGuildC=btwGuildC and strong'
  # 
  #params.data$wGuildC[params.data$wgcIgnore == 1] = params.data$btwGuildC[params.data$wgcIgnore == 1]
  
  params.data$intTypes[params.data$intTypes == 1] = 'competition'
  params.data$intTypes[params.data$intTypes == 2] = 'mutualism'
  params.data$intTypes[params.data$intTypes == 3] = 'predator-prey'
  
  str(params.data)
}

### Prepare saving of stability and network results ###
{
  if ( collate.stab == 1 ){
    # List all stability result files
    setwd(data.dir.1)  
    stab.files <- dir(pattern='ResMat*')
    #str(stab.files)
    
    # Check that number of stability result files equals the number of scenarios.
    if ( length(stab.files) != nrow(params.data)  ) print('Dimension error of stability result data!')
    
    # Load first file
    stab.foo <- readMat(stab.files[1])$resMat
    #str(stab.foo)
  }
  
  if ( collate.nwm == 1 ){
    
    setwd(data.dir.2)
    
    if (params.data$project==1){
      
      # List all network result files
      c.files <- dir(pattern='jacobianCheck*')
      
      # Check that number of network result files equals the number of scenarios.
      if ( length(c.files) != nrow(params.data)  ) print('Dimension error of check network result data!')
      
      # Load first file
      load(c.files[1])
      
    } else if (params.data$project==2){
      
      # List all network result files
      c.files <- dir(pattern='checkMetrics*')
      nwm.files <- dir(pattern='bipNetworkMetrics*')
      rcSum.files <- dir(pattern='rowAndColSumMetrics*')
      #str(nwm.files)
      
      # Check that number of network result files equals the number of scenarios.
      if ( length(c.files) != nrow(params.data)  ) print('Dimension error of check network result data!')
      if ( length(nwm.files) != nrow(params.data)  ) print('Dimension error of network result data!')
      if ( length(rcSum.files) != nrow(params.data)  ) print('Dimension error of network row and col sum result data!')
      
      # Load first file
      c.foo <- read.csv(c.files[1], sep=' ')
      nwm.foo <- read.csv(nwm.files[1], sep=' ')
      rcSum.foo <- read.csv(rcSum.files[1], sep=' ')
      #str(nwm.foo)
    }
  }
  
  # Create the data frame to store the results for all scenarios and replicates!
  if ( collate.stab == 1 && collate.nwm == 0 ){
    
    scenarios <- params.data$scenID
    nscen     <- nrow(params.data)
    nreps     <- nrow(stab.foo)
    
    res.frame <- as.data.frame(matrix(NA, ncol=ncol(stab.foo)+ncol(params.data), nrow=nreps*nscen))
    names(res.frame) <- c(names(params.data),  resultVarNames)
    str(res.frame)
  } 
  
  if ( collate.stab == 0 && collate.nwm == 1 ){
    
    if (params.data$project==1){
      
      scenarios <- params.data$scenID
      nscen     <- nrow(params.data)
      nreps     <- nrow(c.frame)
      
      res.frame <- as.data.frame(matrix(NA, ncol=ncol(c.frame)+ncol(corr.frame)+ncol(params.data)-2, nrow=nreps*nscen ))
      names(res.frame) <- c(names(params.data), names(c.frame)[!(names(c.frame) %in% c('scenID'))], names(corr.frame)[!(names(corr.frame) %in% c('scenID'))])
      str(res.frame)
      
      expC.frame <- as.data.frame(matrix(NA, ncol=ncol(e.frame)+ncol(params.data)-1, nrow=nreps*nscen ))
      names(expC.frame) <- c(names(params.data), names(e.frame)[!(names(e.frame) %in% c('scenID'))])
      str(expC.frame)
      
    } else if (params.data$project==2){
      
      scenarios <- unique(nwm.foo$scenID)
      nscen     <- length(scenarios)
      nreps     <- nrow(nwm.foo)/nscen
      
      check.frame <- as.data.frame(matrix(NA, ncol=ncol(c.foo)+ncol(params.data)-1, nrow=nreps*nscen ))
      names(check.frame) <- c(names(params.data), names(c.foo)[!(names(c.foo) %in% c('scenID'))])
      str(check.frame)
      
      res.frame <- as.data.frame(matrix(NA, ncol=ncol(nwm.foo)+ncol(params.data)-1, nrow=nreps*nscen ))
      names(res.frame) <- c(names(params.data), names(nwm.foo)[!(names(nwm.foo) %in% c('scenID'))])
      str(res.frame)
      
      rc.frame <- as.data.frame(matrix(NA, ncol=ncol(rcSum.foo)+ncol(params.data)-1, nrow=nreps*nscen ))
      names(rc.frame) <- c(names(params.data), names(rcSum.foo)[!(names(rcSum.foo) %in% c('scenID'))])
      str(rc.frame)
    }
  } 
  
  if ( collate.stab == 1 && collate.nwm == 1 ){
    
    if (params.data$project==1){
      
      scenarios <- params.data$scenID
      nscen     <- nrow(params.data)
      nreps     <- nrow(c.frame)
      print(scenarios)
      print(nscen)
      print(nreps)
      
      res.frame <- as.data.frame(matrix(NA, ncol=ncol(stab.foo)+ncol(c.frame)+ncol(corr.frame)+ncol(params.data)-2, nrow=nreps*nscen ))
      names(res.frame) <- c(names(params.data),  resultVarNames, names(c.frame)[!(names(c.frame) %in% c('scenID'))], names(corr.frame)[!(names(corr.frame) %in% c('scenID'))])
      str(res.frame)
      
      expC.frame <- as.data.frame(matrix(NA, ncol=ncol(e.frame)+ncol(params.data)-1, nrow=nreps*nscen ))
      names(expC.frame) <- c(names(params.data), names(e.frame)[!(names(e.frame) %in% c('scenID'))])
      str(expC.frame)
      
    } else if (params.data$project==2){
      
      scenarios <- params.data$scenID
      nscen     <- nrow(params.data)
      # nreps     <- nrow(stab.foo)
      nreps     <- nrow(nwm.foo)
      print(scenarios)
      print(nscen)
      print(nreps)
      
      check.frame <- as.data.frame(matrix(NA, ncol=ncol(c.foo)+ncol(params.data)-1, nrow=nreps*nscen ))
      names(check.frame) <- c(names(params.data), names(c.foo)[!(names(c.foo) %in% c('scenID'))])
      str(check.frame)
      
      res.frame <- as.data.frame(matrix(NA, ncol=ncol(stab.foo)+ncol(nwm.foo)+ncol(params.data)-1, nrow=nreps*nscen ))
      names(res.frame) <- c(names(params.data),  resultVarNames, names(nwm.foo)[!(names(nwm.foo) %in% c('scenID'))])
      str(res.frame)
      
      rc.frame <- as.data.frame(matrix(NA, ncol=ncol(rcSum.foo)+ncol(params.data)-1, nrow=nreps*nscen ))
      names(rc.frame) <- c(names(params.data), names(rcSum.foo)[!(names(rcSum.foo) %in% c('scenID'))])
      str(rc.frame)
    }
  }
}

### Load and Analyse stability results ###
{
  print('Entering Loop')
  # Load all results into res.frame
  for ( i in scenarios){
    # print(i)
    # browser()
    # print('test1')
    
     # print('test2')
    
    if ( collate.stab == 1 ){
      
      # print('collate.stab is True')
      res.frame[seq(i*nreps-nreps+1, i*nreps, 1), 1:ncol(params.data)] <- params.data[params.data$scenID==i,]  # Save params.data data
      res.frame$trials[seq(i*nreps-nreps+1, i*nreps, 1)] <- 1:nreps  # Save replicate numbers
      # print('test3')
      
      stab.file <- paste('ResMat_scenID_', i, '.mat', sep='')
      
      setwd(data.dir.1)
      
      filePresent <- length(dir(pattern=stab.file))>0
      if ( filePresent ){
        stab.foo <- readMat(stab.file)$resMat
        #str(stab.foo)
        res.frame[seq(i*nreps-nreps+1, i*nreps, 1), seq(1+ncol(params.data),ncol(params.data)+ncol(stab.foo), 1)] <- stab.foo[1:nreps,]
      } else{
        res.frame[seq(i*nreps-nreps+1, i*nreps, 1), seq(1+ncol(params.data),ncol(params.data)+ncol(stab.foo), 1)] <- NA
      }
      
    }
    
    # print('test4')
    
    if ( collate.nwm == 1 ){
      
      if (params.data$project==1){
        
        expC.frame[seq(i*nreps-nreps+1, i*nreps, 1), 1:ncol(params.data)] <- params.data[params.data$scenID==i,] 
       
        #print('collate.nwm is True')
        setwd(data.dir.2) 
        
        nwm.file <- paste('jacobianCheck_scenario_', i, '.RData', sep='')
        filePresent <- length(dir(pattern=nwm.file))>0
        
        if ( filePresent ){
          
          load(nwm.file)
          if ( i != unique(c.frame$scenID)){
            print('ID error in data!')
            rm('i') # should throw error it and terminate analysis.
          } 
          
          expC.frame[seq(i*nreps-nreps+1, i*nreps, 1), seq(1+ncol(params.data),ncol(params.data)+ncol(e.frame)-1, 1)] <- e.frame[, !(names(e.frame) %in% c('scenID'))]
          
          if ( collate.stab == 1){
            res.frame[seq(i*nreps-nreps+1, i*nreps, 1), seq(1+ncol(params.data)+ncol(stab.foo), ncol(params.data)+ncol(stab.foo)+ncol(c.frame)-1, 1)] <- c.frame[, !(names(c.frame) %in% c('scenID'))]
            res.frame[seq(i*nreps-nreps+1, i*nreps, 1), seq(1+ncol(params.data)+ncol(stab.foo)+ncol(c.frame), ncol(params.data)+ncol(stab.foo)+ncol(c.frame)+ncol(corr.frame)-2, 1)] <- corr.frame[, !(names(corr.frame) %in% c('scenID'))]
          } else {
            res.frame[seq(i*nreps-nreps+1, i*nreps, 1), 1:ncol(params.data)] <- params.data[params.data$scenID==i,]  # Save params.data data
            res.frame$trials[seq(i*nreps-nreps+1, i*nreps, 1)] <- 1:nreps  # Save replicate numbers
            res.frame[seq(i*nreps-nreps+1, i*nreps, 1), seq(1+ncol(params.data), ncol(params.data)+ncol(c.frame)-1, 1)] <- c.frame[, !(names(c.frame) %in% c('scenID'))]
            res.frame[seq(i*nreps-nreps+1, i*nreps, 1), seq(1+ncol(params.data)+ncol(c.frame), ncol(params.data)+ncol(c.frame)+ncol(corr.frame)-2, 1)] <- corr.frame[, !(names(corr.frame) %in% c('scenID'))]
          }
          
        }else{
          
          expC.frame[seq(i*nreps-nreps+1, i*nreps, 1), seq(1+ncol(params.data),ncol(params.data)+ncol(e.frame)-1, 1)] <- NA
          
          if ( collate.stab == 1){
            res.frame[seq(i*nreps-nreps+1, i*nreps, 1), seq(1+ncol(params.data)+ncol(stab.foo), ncol(params.data)+ncol(stab.foo)+ncol(c.frame)-1, 1)] <- NA
            res.frame[seq(i*nreps-nreps+1, i*nreps, 1), seq(1+ncol(params.data)+ncol(stab.foo)+ncol(c.frame), ncol(params.data)+ncol(stab.foo)+ncol(c.frame)+ncol(corr.frame)-2, 1)] <- NA
          } else {
            res.frame[seq(i*nreps-nreps+1, i*nreps, 1), 1:ncol(params.data)] <- params.data[params.data$scenID==i,]  # Save params.data data
            res.frame$trials[seq(i*nreps-nreps+1, i*nreps, 1)] <- 1:nreps  # Save replicate numbers
            res.frame[seq(i*nreps-nreps+1, i*nreps, 1), seq(1+ncol(params.data), ncol(params.data)+ncol(c.frame)-1, 1)] <- NA
            res.frame[seq(i*nreps-nreps+1, i*nreps, 1), seq(1+ncol(params.data)+ncol(c.frame), ncol(params.data)+ncol(c.frame)+ncol(corr.frame)-2, 1)] <- NA
          }
          
        }

      } else if (params.data$project==2){
        
        check.frame[seq(i*nreps-nreps+1, i*nreps, 1), 1:ncol(params.data)] <- params.data[params.data$scenID==i,] 
        rc.frame[seq(i*nreps-nreps+1, i*nreps, 1), 1:ncol(params.data)] <- params.data[params.data$scenID==i,] 
        
        #print('collate.nwm is True')
        setwd(data.dir.2) 
        
        c.file <- paste('checkMetrics_scenario_', i, '.csv', sep='')
        filePresent <- length(dir(pattern=c.file))>0
        if ( filePresent ){
          c.foo  <- read.csv(c.file, sep=' ')
          check.frame[seq(i*nreps-nreps+1, i*nreps, 1), seq(1+ncol(params.data),ncol(params.data)+ncol(c.foo)-1, 1)] <- c.foo[, !(names(c.foo) %in% c('scenID'))]
        }else{
          check.frame[seq(i*nreps-nreps+1, i*nreps, 1), seq(1+ncol(params.data),ncol(params.data)+ncol(c.foo)-1, 1)] <- NA
        }
        
        rcSum.file  <- paste('rowAndColSumMetrics_scenario_', i, '.csv', sep='')
        filePresent <- length(dir(pattern=rcSum.file))>0
        if ( filePresent ){
          rc.foo  <- read.csv(rcSum.file, sep=' ')
          rc.frame[seq(i*nreps-nreps+1, i*nreps, 1), seq(1+ncol(params.data),ncol(params.data)+ncol(rc.foo)-1, 1)] <- rc.foo[, !(names(rc.foo) %in% c('scenID'))]
        }else{
          rc.frame[seq(i*nreps-nreps+1, i*nreps, 1), seq(1+ncol(params.data),ncol(params.data)+ncol(rc.foo)-1, 1)] <- NA
        }
        
        nwm.file <- paste('bipNetworkMetrics_scenario_', i, '.csv', sep='')
        filePresent <- length(dir(pattern=nwm.file))>0
        if ( filePresent ){
          nwm.foo  <- read.csv(nwm.file, sep=' ')
          
          if ( i != unique(nwm.foo$scenID)){
            print('ID error in data!')
            rm('i') # should throw error it and terminate analysis.
          } 
          
          if ( collate.stab == 1){
            res.frame[seq(i*nreps-nreps+1, i*nreps, 1), seq(1+ncol(params.data)+ncol(stab.foo), ncol(params.data)+ncol(stab.foo)+ncol(nwm.foo)-1, 1)] <- nwm.foo[, !(names(nwm.foo) %in% c('scenID'))]
          } else {
            res.frame[seq(i*nreps-nreps+1, i*nreps, 1), 1:ncol(params.data)] <- params.data[params.data$scenID==i,]  # Save params.data data
            res.frame$trials[seq(i*nreps-nreps+1, i*nreps, 1)] <- 1:nreps  # Save replicate numbers
            res.frame[seq(i*nreps-nreps+1, i*nreps, 1), seq(1+ncol(params.data),ncol(params.data)+ncol(nwm.foo)-1, 1)] <- nwm.foo[, !(names(nwm.foo) %in% c('scenID'))]
          }
          
        }else{
          
          if ( collate.stab == 1){
            res.frame[seq(i*nreps-nreps+1, i*nreps, 1), seq(1+ncol(params.data)+ncol(stab.foo), ncol(params.data)+ncol(stab.foo)+ncol(nwm.foo)-1, 1)] <- NA
          } else {
            res.frame[seq(i*nreps-nreps+1, i*nreps, 1), 1:ncol(params.data)] <- params.data[params.data$scenID==i,]  # Save params.data data
            res.frame$trials[seq(i*nreps-nreps+1, i*nreps, 1)] <- 1:nreps  # Save replicate numbers
            res.frame[seq(i*nreps-nreps+1, i*nreps, 1), seq(1+ncol(params.data),ncol(params.data)+ncol(nwm.foo)-1, 1)] <- NA
          }
          
        } # end of if ( filePresent ) statement for nwm
        
        
      }
      
    } #end of if ( collate.nwm == 1 )
    
    # print('test5')
    
  } # end of loop
  
  print('Done w Loop')
  str(res.frame)
  
}

### Subset the full result frame ###
{
  
  if ( do.subset ){
    
    sub.frame <- res.frame[,!(names(res.frame) %in% colDrop)]
    sub.frame <- sub.frame[sub.frame$wGuildC!=1,]
    sub.frame <- sub.frame[sub.frame$wGuildWeightStd!=0.5,]
    sub.frame <- sub.frame[sub.frame$intTypes!=levels(sub.frame$intTypes)[1],]
    sub.frame <- sub.frame[!(sub.frame$intTypes==levels(sub.frame$intTypes)[2] & sub.frame$S>10),]
    sub.frame <- sub.frame[!(sub.frame$intTypes==levels(sub.frame$intTypes)[3] & sub.frame$S!=128),]
    sub.frame <- sub.frame[!(sub.frame$btwGuildC==min(sub.frame$btwGuildC[sub.frame$S==10]) & sub.frame$S==10),]
    sub.frame <- sub.frame[!(sub.frame$btwGuildC==min(sub.frame$btwGuildC[sub.frame$S==128]) & sub.frame$S==128),]
    sub.frame <- droplevels(sub.frame)
    str(sub.frame)
  }
  
}

### Save summarized results - full and subsetted ###
{
  
  if ( collate.stab == 1 && collate.nwm == 0 ){
    
    savename1 <- 'collatedStabResults.csv'
    savename2 <- 'collatedStabResultsSubsetted.csv'
    
  } else if ( collate.stab == 0 && collate.nwm == 1 ){
    
    if ( params.data$project==1 ){
      
      savename1 <- 'collatedNetworkResults.csv'
      savename2 <- 'collatedNetworkResultsSubsetted.csv'
      savename3 <- 'collatedExpectedJacobianCheckMetrics.csv'

      
    } else if ( params.data$project==2 ){
      
      savename1 <- 'collatedNetworkResults.csv'
      savename2 <- 'collatedNetworkResultsSubsetted.csv'
      savename3 <- 'collatedJacobianCheckMetrics.csv'
      savename4 <- 'collatedRowandColSumMetrics.csv'
      
    }
    
  } else if ( collate.stab == 1 && collate.nwm == 1 ){
    
    if ( params.data$project==1 ){
      
      savename1 <- 'collatedStabAndNetworkResults.csv'
      savename2 <- 'collatedStabAndNetworkResultsSubsetted.csv'
      savename3 <- 'collatedExpectedJacobianCheckMetrics.csv'
      
      
    } else if ( params.data$project==2 ){
      
      savename1 <- 'collatedStabAndNetworkResults.csv'
      savename2 <- 'collatedStabAndNetworkResultsSubsetted.csv'
      savename3 <- 'collatedJacobianCheckMetrics.csv'
      savename4 <- 'collatedRowandColSumMetrics.csv'
      
    }
    
  }
  
  setwd(save.dir)
  write.table(res.frame, savename1)
  
  if ( collate.nwm == 1 ){
    if ( params.data$project==1 ){
      
      write.table(expC.frame, savename3)
      
    } else if ( params.data$project==2 ) {
      
      write.table(check.frame, savename3)
      write.table(rc.frame, savename4)
      
    }
  }
  
  if ( do.subset ){
    write.table(sub.frame, savename2)
  }
  
}  