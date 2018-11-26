##############################################
### summarize stability results on cluster ###
##############################################
#
#   @DESCRIPTION
#   Analysis code for project SSMNW.
#   Summarize local stability for all scenarios and 
#   and save new summary data frame.
#   Subset the data frame if and as wanted.
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
  require('R.matlab') # To read .mat files and such
  
  # find out which type of machine you're running on
  os <- .Platform$OS.type
  
  # Reading arguments, defining paths
  if ( os == 'windows' ){ #i.e. my local pc. only for development. 
    
    # define the simulation folder to work with
    timestamp <- '2018Oct20_0125'
    
    # Data paths and such
    main.dir  <- 'C:/Users/icurtsd/Documents/Project Atlanta/Project Structure and stability of mutualistic networks/data/'
    sub.dir   <- 'clusterdata/analysis/' # 'main_test/' # 
    stupid1   <- '' #'testresults_' #  
   
    data.dir <- paste(main.dir, sub.dir, stupid1, timestamp, sep='')
    save.dir <- data.dir 
    
  } else if ( os == 'unix'){ #i.e. cluster [or a Mac unfortunately]
    
    # Reading command line input (when on cluster)
    args=(commandArgs())
    print(args) 
    eval(parse(text = args[[10]])) # gets the variable timestamp
    
    # Data paths and such
    data.dir <- paste('/proj/ecoservice/users/x_alvcu/SSMNW/Results/', timestamp,'/', sep='')  
    save.dir <- paste('/proj/ecoservice/users/x_alvcu/SSMNW/analysis/', timestamp,'/', sep='')  
    
  }
  print(data.dir)
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
  setwd(data.dir)
  
  # Read parameterStruct for the simulation
  file1  <- dir(pattern = 'paramStruct')
  params <- readMat(file1)
  
  # Get the paramStruct into shape.
  params.data <- lapply(params$params[2], unlist, use.names=FALSE)
  params.data <- as.data.frame(params.data)
  names(params.data) <- c('scenID', lapply(params$params[1], unlist, use.names=FALSE)[[1]])
  str(params.data)
  
  if ( is.element(unique(params.data$project),2)  ){
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
  
  params.data$wGuildC[params.data$wgcIgnore == 1] = params.data$btwGuildC[params.data$wgcIgnore == 1]
  
  params.data$intTypes[params.data$intTypes == 1] = 'competition'
  params.data$intTypes[params.data$intTypes == 2] = 'mutualism'
  params.data$intTypes[params.data$intTypes == 3] = 'predator-prey'
  
  str(params.data)
}

### Load and Analyse stability results ###
{
  
  # List all result files
  setwd(data.dir)  
  files <- dir(pattern='ResMat*')
  str(files)
  
  # Check that number of result files equals the number of scenarios.
  if ( length(files) != nrow(params.data)  ) print('Dimension error of result data!')
  
  # Create the data frame to store the summarized results for all scenarios!
  foo <- readMat(files[1])$resMat
  res.frame <- as.data.frame(matrix(NA, ncol=ncol(foo)+ncol(params.data), nrow=nrow(params.data)))
  names(res.frame) <- c(names(params.data), resultVarNames)
  
  # Load all results into res.frame
  for ( i in 1:length(files)){
    
    file <- files[i]
    
    foo <- readMat(file)$resMat
    #str(foo)
    
    scenID <- eval(parse(text=(substr(file, nchar("ResMat_scenID_")+1, nchar(file)-nchar(".mat")))))
    
    res.frame[i, 1:ncol(params.data)] <- params.data[params.data$scenID==scenID,]
    res.frame[i, 1+ncol(params.data)] <- 100*sum(foo[,1])/nrow(foo)
    res.frame[i, 2+ncol(params.data)] <- mean(foo[,2], na.rm=T)
    res.frame[i, 3+ncol(params.data)] <- mean(foo[,3], na.rm=T)
    
  }
  
  # Make characters into factors
  # res.frame$wGcode   <- factor(res.frame$wGcode, levels = c('no', 'intraGuildC=btwGuildC and weak', 'intraGuildC=1 and weak', 'intraGuildC=btwGuildC and strong', 'intraGuildC=1 and strong'))
  # res.frame$nestCode <- factor(res.frame$nestCode, levels = c('bNqN', 'bNqR', 'bNqA', 'bRqN', 'bRqR', 'bRqA','bAqN', 'bAqR', 'bAqA'))
  res.frame$intTypes <- factor(res.frame$intTypes, levels = c('competition', 'mutualism', 'predator-prey'))
  
  res.frame <- res.frame[order(res.frame$scenID),]
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
  
  setwd(save.dir)
  write.table(res.frame, 'stabilitySummary.csv')
  
  if ( do.subset ){
    write.table(sub.frame, 'stabilitySummarySubsetted.csv')
  }
  
}  
  