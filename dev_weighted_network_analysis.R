##############################################
########## SSMNW WEIGHTED NETWORK ANALYSIS ############
##############################################
#
#   @DESCRIPTION
#   Development code.
#   Calculates network measures..
#   Written for SSMNW, to calculate
#   quantitative nestedness for the weighted 
#   incidence matrices in the pilot project.
#   Might include other metrics as well.
#
#   @AUTHOR
#   Alva Curtsdotter, postdoc @ Emory University, 
#   started during Stanford visit, 2017-12-05
#   
##############################################

# Setup --------------------------------------------------------------------------------------------
{
  require(bipartite)
  require(vegan)
  
  SAVE <- 1 # To save files or not.
  if (SAVE){
    save.dir <- paste('//storage.slu.se/Home$/aacu0001/My Documents/Project Atlanta/Project Structure and stability of mutualistic networks/data/toydata/', format(Sys.time(), format="%Y%b%dT%H%M"), '_dev_network_analysis_output', sep='')
    if ( length(dir(save.dir)) == 0 ){
      dir.create(save.dir)
    }
  }
  
  signTypes <- c('positive', 'predpreylike')# c('positive', 'negative')# Positive and negatvive yield the same result.
  trials   <- 1000
}

# Topologies ---------------------------------------------------------------------------------------------
{
  #% Nested weight order
  T1O1 <- matrix(c(21,19,18,12,11,1, 20,17,13,10,2,0, 16,14,9,3,0,0, 15,8,4,0,0,0, 7,5,0,0,0,0, 6,0,0,0,0,0), nrow=6, ncol=6, byrow=T)
  T2O1 <- matrix(c(21,19,18,13,0,6, 20,17,14,0,7,5, 16,0,12,8,4,0, 15,11,0,3,0,0, 10,0,0,0,1,0, 9,2,0,0,0,0), nrow=6, ncol=6, byrow=T)
  T3O1 <- matrix(c(21,0,19,0,15,7, 20,18,0,0,8,6, 0,16,14,9,5,0, 17,13,10,0,0,0, 12,11,0,3,0,0, 0,0,4,2,0,1), nrow=6, ncol=6, byrow=T)
  
  #% Intermediately nested weight order
  T1O2 <- matrix(c(11,18,4,8,20,5, 13,7,6,15,14,0, 3,17,21,19,0,0, 10,2,16,0,0,0, 9,12,0,0,0,0, 1,0,0,0,0,0), nrow=6, ncol=6, byrow=T)
  T2O2 <- matrix(c(11,18,12,21,0,14, 13,7,4,0,15,5, 3,0,6,16,19,0, 10,17,0,8,0,0, 9,0,0,0,20,0, 1,2,0,0,0,0), nrow=6, ncol=6, byrow=T)
  T3O2 <- matrix(c(11,0,17,0,8,20, 13,9,0,0,15,14, 0,1,2,6,19,0, 3,18,12,0,0,0, 10,7,0,21,0,0,  0,0,4,16,0,5), nrow=6, ncol=6, byrow=T)
  
  #% Antinested weight order
  T1O3 <- matrix(c(1,2,6,7,15,16, 3,5,8,14,17,0, 4,9,13,18,0,0, 10,12,19,0,0,0, 11,20,0,0,0,0, 21,0,0,0,0,0), nrow=6, ncol=6, byrow=T)
  T2O3 <- matrix(c(1,2,6,7,0,13, 3,5,8,0,14,20, 4,0,12,15,19,0, 9,11,0,18,0,0, 10,0,0,0,21,0, 16,17,0,0,0,0), nrow=6, ncol=6, byrow=T)
  T3O3 <- matrix(c(1,0,3,0,7,15, 2,4,0,0,14,16, 0,6,8,13,17,0, 5,9,12,0,0,0, 10,11,0,19,0,0, 0,0,18,20,0,21), nrow=6, ncol=6, byrow=T)
  
  networks <-  c('T1O1', 'T1O2', 'T1O3', 'T2O1', 'T2O2', 'T2O3', 'T3O1', 'T3O2', 'T3O3')#c('T1', 'T2', 'T3')#c('T1q1', 'T1q2',  'T1q3', 'T1aq1', 'T1aq2', 'T1aq3', 'T3q1','T3aq1', 'T3aq2', 'T3aq3') #c('T3q1', 'T3q1t', 'T3aq1', 'T3aq1t', 'T3aq2', 'T3aq2t', 'T3aq3') #c('T1', 'T2', 'T3') #c('q1', 'q2',  'q3', 'aq1', 'aq2', 'aq3') #c('q1', 'q1t', 'q2', 'q2t', 'q3', 'q3t', 'aq1', 'aq1t', 'aq2', 'aq2t', 'aq3', 'aq3t')#c('T1', 'T2a', 'T2b', 'T3a', 'T3b')#
  
  if ( SAVE ){
    setwd(save.dir)
    for (i in 1:length(networks)){
      filename = paste(save.dir, '/network', networks[i], '.csv', sep='')
      write.table(eval(parse(text = networks[i])), filename)
    }
  }
}

# CREATE WEIGHTED TOPOLOGIES ----------------------------------------------

# Setup structures to hold results
{
  numScen <- length(networks)*length(signTypes)
  
  sumStats.qmets.B <- as.data.frame(matrix(NA,ncol=10,nrow=numScen))
  colnames(sumStats.qmets.B) <- c( 'network', 'signType', 'mean.wNODF', 'sd.wNODF', 'mean.wine', 'sd.wine', 'mean.specRad', 'sd.specRad', 'mean.Fcomp', 'sd.Fcomp')
  #str(sumStats.qmets.B)
  
  sumStats.qmets.tB <- as.data.frame(matrix(NA,ncol=10,nrow=numScen))
  colnames(sumStats.qmets.tB) <- c( 'network', 'signType', 'mean.wNODF', 'sd.wNODF', 'mean.wine', 'sd.wine', 'mean.specRad', 'sd.specRad', 'mean.Fcomp', 'sd.Fcomp')
  #str(sumStats.qmets.tB)
  # Note that the transpose only seems to matter to Fcomp.
  # And also a little bit for specRad (not A and tA, but A with B and tB vs A with tB and B).
  
  r <- 1
}

# Generate random weights and distribute them according to order

nlinks <- sum(sum(T1O1>0))
meanW <- 0
stdW <- sqrt(1/4)

for ( i in 1:length(networks)){
  B <- eval(parse(text=c(networks[i])))
  
  for ( j in 1:length(signTypes)){
    signType <- signTypes[j]
    
    quant.metrics.B <- as.data.frame(matrix(NA,ncol=5,nrow=trials))
    colnames(quant.metrics.B) <- c('trial', 'wNODF', 'wine', 'specRad', 'Fcomp')
    #str(quant.metrics.B)
    
    quant.metrics.tB <- as.data.frame(matrix(NA,ncol=5,nrow=trials))
    colnames(quant.metrics.tB) <- c('trial', 'wNODF', 'wine', 'specRad', 'Fcomp')
    #str(quant.metrics.tB)
    
    for ( k in 1:trials){
      
      quant.metrics.B$trial[k]  <- k
      quant.metrics.tB$trial[k] <- k
      
      weights = abs(rnorm(nlinks, mean=meanW, sd=stdW))
      weights = sort(weights)
      
      if (signType == 'positive' || signType == 'predpreylike'){
        
        for ( l in 1:nlinks ){
          B[B==l] <- weights[l]    # Distribute weights according to order.
        }
        
      }  else if (signType == 'negative'){
        
        for ( l in 1:nlinks ){
          B[B==l] <- -weights[l]   # Distribute weights according to order.
        }
        
      } # End of if signType
      
      tB <- t(B)
      
      # Quantitative nestedness -----------------------------------------------------------------------
      
      quant.metrics.B$wNODF[k]  <- nested(B, method='weighted NODF', rescale=F, normalised=F) # doubtful that rescale works for all metrics so better to not use it!
      quant.metrics.tB$wNODF[k] <- nested(tB, method='weighted NODF', rescale=F, normalised=F)
      
      quant.metrics.B$wine[k]  <- nested(B, method='wine', rescale=F, normalised=F) # doubtful that rescale works for all metrics so better to not use it!
      quant.metrics.tB$wine[k] <- nested(tB, method='wine', rescale=F, normalised=F)
      
      # Spectral method to nestedness -------------------------------------------------------------------------
      
      zeroblock  <- matrix(0,nrow=nrow(B),ncol=ncol(B))
      
      if (signType == 'predpreylike'){
        A          <- cbind(rbind(zeroblock,B), rbind(-tB,zeroblock))
        tA         <- cbind(rbind(zeroblock,tB), rbind(-B,zeroblock))
      }else{
        A          <- cbind(rbind(zeroblock,B), rbind(tB,zeroblock))
        tA         <- cbind(rbind(zeroblock,tB), rbind(B,zeroblock))
      }
      
      quant.metrics.B$specRad[k]  <- max(Re(eigen(A)$values))
      quant.metrics.tB$specRad[k] <- max(Re(eigen(tA)$values))
      # Note that I am expecting the specRad for A and tA to be the same.
      # But do a dummycheck for that later on.
      
      # Quantitative complementarity ----------------------------------------------------------------
      
      # Implementing Tylianakis' complementarity measure, Fcomp (is identical to Petchey's Trophic Diversity if made only on the predation matrix rather than community matrix, independent sum of predation and rsource matrix, which is what Petchey et al actually used...)
      
      ppm <- B
      distance = vegdist(t(ppm), method="euclidean", bin=F)
      distance[is.nan(distance)] <- 0
      tree <- hclust(distance)
      quant.metrics.B$Fcomp[k] <- treeheight(tree)
      
      ppm <- tB
      distance = vegdist(t(ppm), method="euclidean", bin=F)
      distance[is.nan(distance)] <- 0
      tree <- hclust(distance)
      quant.metrics.tB$Fcomp[k] <- treeheight(tree)
      
      
      # Saving some weighted networks
      if ( SAVE && k == 1 ){
        setwd(save.dir)
        filename = paste(save.dir, '/network_', networks[i], '_signType_', signType, '_rep_', k, '.csv', sep='')
        write.table(B, filename)
      }
      
    } # End of for trials
    
    # Put summaryStats for this scenario into the saving structure
    {
      Bmeans <- colMeans(quant.metrics.B)
      Bsds   <- sqrt(colSums((quant.metrics.B-kronecker(matrix(1,nrow(quant.metrics.B),1),t(Bmeans)))^2)/(trials-1))
      
      sumStats.qmets.B$network[r]  <- networks[i]
      sumStats.qmets.B$signType[r] <- signType
      
      sumStats.qmets.B$mean.wNODF[r]   <- Bmeans[2]
      sumStats.qmets.B$mean.wine[r]    <- Bmeans[3]
      sumStats.qmets.B$mean.specRad[r] <- Bmeans[4]
      sumStats.qmets.B$mean.Fcomp[r]   <- Bmeans[5]
      
      sumStats.qmets.B$sd.wNODF[r]   <- Bsds[2]
      sumStats.qmets.B$sd.wine[r]    <- Bsds[3]
      sumStats.qmets.B$sd.specRad[r] <- Bsds[4]
      sumStats.qmets.B$sd.Fcomp[r]   <- Bsds[5]
      
      tBmeans <- colMeans(quant.metrics.tB)
      tBsds   <- sqrt(colSums((quant.metrics.tB-kronecker(matrix(1,nrow(quant.metrics.tB),1),t(tBmeans)))^2)/(trials-1))
      
      sumStats.qmets.tB$network[r]  <- networks[i]
      sumStats.qmets.tB$signType[r] <- signType
      
      sumStats.qmets.tB$mean.wNODF[r]   <- tBmeans[2]
      sumStats.qmets.tB$mean.wine[r]    <- tBmeans[3]
      sumStats.qmets.tB$mean.specRad[r] <- tBmeans[4]
      sumStats.qmets.tB$mean.Fcomp[r]   <- tBmeans[5]
      
      sumStats.qmets.tB$sd.wNODF[r]   <- tBsds[2]
      sumStats.qmets.tB$sd.wine[r]    <- tBsds[3]
      sumStats.qmets.tB$sd.specRad[r] <- tBsds[4]
      sumStats.qmets.tB$sd.Fcomp[r]   <- tBsds[5]
      
      r <- r+1
    }
    #head(sumStats.qmets.tB)
    #head(sumStats.qmets.B)
    
  } # End of for signTypes
} # End of for networks

print(sumStats.qmets.B)
print(sumStats.qmets.tB)

tDiff <- sumStats.qmets.tB[,3:ncol(sumStats.qmets.tB)]-sumStats.qmets.B[,3:ncol(sumStats.qmets.B)]


if ( SAVE ){
  setwd(save.dir)
  filename = paste(save.dir, '/summaryStats_weightedIncidenceMats.csv', sep='')
  write.table(sumStats.qmets.B, filename)
  filename = paste(save.dir, '/summaryStats_weightedIncidenceMatsTransposed.csv', sep='')
  write.table(sumStats.qmets.tB, filename)
  filename = paste(save.dir, '/summaryStats_diffWeightedIncidenceMatsAndTranspose.csv', sep='')
  write.table(tDiff, filename)
}

