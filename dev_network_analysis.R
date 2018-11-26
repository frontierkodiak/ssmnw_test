##############################################
########## SSMNW NETWORK ANALYSIS ############
##############################################
#
#   @DESCRIPTION
#   Development code.
#   Calculates network measures..
#   Written for SSMNW, to calculate binary and
#   quantitative nestedness for the pilot project.
#   Might include other metrics as well.
#
#   @AUTHOR
#   Alva Curtsdotter, postdoc @ Emory University, 
#   started during Stanford visit, 2017-12-05
#   
##############################################

# Setup --------------------------------------------------------------------------------------------

require(bipartite)
require(vegan)

SAVE <- 1 # To save files or not.
save.dir <- paste('//storage.slu.se/Home$/aacu0001/My Documents/Project Atlanta/Project Structure and stability of mutualistic networks/data/toydata/', format(Sys.time(), format="%Y%b%dT%H%M"), '_dev_network_analysis_output', sep='')
if ( length(dir(save.dir)) == 0 ){
  dir.create(save.dir)
}

# Topologies ---------------------------------------------------------------------------------------------

# Discarded 18 links webs. Had chosen T1, T2b, and T3b. T1, T2, T3 21 link webs are developments of these three.
#T1  <- matrix(c(1,1,1,1,1,1, 1,1,1,1,1,0, 1,1,1,0,0,0, 1,1,0,0,0,0, 1,0,0,0,0,0, 1,0,0,0,0,0), nrow=6, ncol=6, byrow=T)
#T2a <- matrix(c(1,1,1,1,1,0, 1,1,1,0,0,1, 1,1,0,1,0,0, 1,0,1,0,0,0, 1,0,0,1,0,1, 0,0,1,0,0,0), nrow=6, ncol=6, byrow=T)
#T2b <- matrix(c(1,1,1,1,1,0, 1,0,1,1,0,1, 0,1,0,0,1,1, 1,0,0,0,0,1, 1,0,0,1,0,1, 0,0,1,0,0,0), nrow=6, ncol=6, byrow=T)
#T3a <- matrix(c(1,0,0,0,1,1, 1,1,0,0,0,1, 1,1,1,0,0,0, 0,1,1,1,0,0, 0,0,1,1,1,0, 0,0,0,1,1,1), nrow=6, ncol=6, byrow=T)
#T3b <- matrix(c(1,1,1,0,0,0, 1,1,0,1,0,0, 1,0,0,0,1,1, 0,1,0,0,1,1, 0,0,1,1,1,0, 0,0,1,1,0,1), nrow=6, ncol=6, byrow=T)

T1 <- matrix(c(1,1,1,1,1,1, 1,1,1,1,1,0, 1,1,1,1,0,0, 1,1,1,0,0,0, 1,1,0,0,0,0, 1,0,0,0,0,0), nrow=6, ncol=6, byrow=T)
T2 <- matrix(c(1,1,1,1,0,1, 1,1,1,0,1,1, 1,0,1,1,1,0, 1,1,0,1,0,0, 1,0,0,0,1,0, 1,1,0,0,0,0), nrow=6, ncol=6, byrow=T)
T3 <- matrix(c(1,0,1,0,1,1, 1,1,0,0,1,1, 0,1,1,1,1,0, 1,1,1,0,0,0, 1,1,0,1,0,0, 0,0,1,1,0,1), nrow=6, ncol=6, byrow=T)

T1q1  <- matrix(c(21,20,19,18,17,16, 15,14,13,12,11,0, 10,9,8,7,0,0, 6,5,4,0,0,0, 3,2,0,0,0,0, 1,0,0,0,0,0), nrow=6, ncol=6, byrow=F) # My interp of Phil\s verbal description in ghost paper.
T1q1t <- t(T1q1)
T1q2  <- matrix(c(21,20,18,15,11,6, 19,17,14,10,5,0, 16,13,9,4,0,0, 12,8,3,0,0,0, 7,2,0,0,0,0, 1,0,0,0,0,0), nrow=6, ncol=6, byrow=F) # More diag disp, always decreasin up right diag.
T1q2t <- t(T1q2)
T1q3  <- matrix(c(21,20,16,15,7,6, 19,17,14,8,5,0, 18,13,9,4,0,0, 12,10,3,0,0,0, 11,2,0,0,0,0, 1,0,0,0,0,0), nrow=6, ncol=6, byrow=F) # serpentine 
T1q3t <- t(T1q3)

T1aq1   <- matrix(c(1,3,6,10,15,21, 2,5,9,14,20,0, 4,8,13,19,0,0, 7,12,18,0,0,0, 11,17,0,0,0,0, 16,0,0,0,0,0), nrow=6, ncol=6, byrow=F) # More diag disp, always decreasin up right diag.
T1aq1t  <- t(T1aq1)
T1aq2   <- matrix(c(1,7,12,16,19,21, 2,8,13,17,20,0, 3,9,14,18,0,0, 4,10,15,0,0,0, 5,11,0,0,0,0, 6,0,0,0,0,0), nrow=6, ncol=6, byrow=F) # Row wise up, always starting at diagonal element
T1aq2t  <- t(T1aq2)
T1aq3   <- matrix(c(1,3,4,10,11,21, 2,5,9,12,20,0, 6,8,13,19,0,0, 7,14,18,0,0,0, 15,17,0,0,0,0, 16,0,0,0,0,0), nrow=6, ncol=6, byrow=F)  # serpentine, starting lower left going right diag up
T1aq3t  <- t(T1aq3)

T3q1  <- matrix(c(21,0,19,0,15,7, 20,18,0,0,8,6, 0,16,14,9,5,0, 17,13,10,0,0,0, 12,11,0,3,0,0, 0,0,4,2,0,1), nrow=6, ncol=6, byrow=T)
T3q1t  <- t(T3q1)

T3aq1 <- matrix(c(1,0,4,0,11,17, 2,5,0,0,18,16, 0,7,12,19,15,0, 8,13,20,0,0,0, 14,21,0,10,0,0, 0,0,9,6,0,3), nrow=6, ncol=6, byrow=T)
T3aq1t  <- t(T3aq1)
T3aq2 <- matrix(c(1,0,3,0,7,15, 2,4,0,0,14,16, 0,6,8,13,17,0, 5,9,12,0,0,0, 10,11,0,19,0,0, 0,0,18,20,0,21), nrow=6, ncol=6, byrow=T)
T3aq2t  <- t(T3aq2)
T3aq3 <- matrix(c(16,0,1,0,4,10, 2,17,0,0,11,5, 0,6,3,12,18,0, 7,13,19,0,0,0, 14,8,0,20,0,0, 0,0,9,15,0,21), nrow=6, ncol=6, byrow=T)

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

# Basic calculations and dummy checks -----------------------------------------------------------------

basic <- as.data.frame(matrix(NA,ncol=length(networks),nrow=2))
colnames(basic) <- networks
rownames(basic) <- c('rows', 'linksum')

for ( i in 1:length(networks) ){
  nw <- eval(parse(text=networks[i]))
  
  basic[,i] <- c(nrow(nw), sum(nw) )
}
print(basic)

radsummor   <- as.data.frame(matrix(NA,ncol=length(networks),nrow=nrow(nw)))
kolsummor   <- as.data.frame(matrix(NA,ncol=length(networks),nrow=ncol(nw)))
kolsumslope <- as.data.frame(matrix(NA,ncol=length(networks),nrow=1))

colnames(radsummor) <- colnames(kolsummor) <- colnames(kolsumslope) <- networks
rownames(radsummor) <- rownames(kolsummor) <- seq(1,nrow(nw))

for ( i in 1:length(networks) ){
  nw <- eval(parse(text=networks[i]))
  
  radsummor[,i] <- c(rowSums(nw))
  kolsummor[,i] <- c(colSums(nw))
  linmodkolsums <- lm( kolsummor[,i] ~ seq(1,length(kolsummor[,i])) )
  kolsumslope[i] <- coefficients(linmodkolsums)[2]
}
print(radsummor)
print(kolsummor)
print(kolsumslope)
print('kolsumslope doesnt make sense unless kolumns are ordered by kolumnsum which it currently isnt in the transposed ones!')

# Organize network, if topological - otherwise don't run this section! -----------------------------------------------------------------------------------------

for ( i in 1:length(networks) ){
  nw <- eval(parse(text=networks[i]))
  
  nw <- nw[order(radsummor[,i],decreasing=T),]
  nw <- nw[,order(kolsummor[,i],decreasing=T)]
  
  eval(parse(text=paste(networks[i], '<-nw',sep=''))) 
}

# Topological (binary) nestedness -----------------------------------------------------------------------
# Actually, two metrics are quantitative no binary, wine and weighted NODF.

topo.nest <- as.data.frame(matrix(NA,ncol=length(networks),nrow=10))
colnames(topo.nest) <- networks

for ( i in 1:length(networks) ){
  nw <- eval(parse(text=networks[i]))
  
  topo.nest[,i] <- nested(nw, method='ALL', rescale=F, normalised=F) # doubtful that rescale works for all metrics so better to not use it!
 
  #nesteddisc(nw)  # vegan
  #nestednodf(nw) # vegan
  #nestedtemp(nw) # vegan
  #nestedn0(nw) # vegan
  
}
rownames(topo.nest) <- names(nested(nw, method='ALL', rescale=T, normalised=T))
print(topo.nest)

if ( SAVE ){
  setwd(save.dir)
  filename = paste(save.dir, '/nestedness.csv', sep='')
  write.table(topo.nest, filename)
}

# Spectral method to nestedness -------------------------------------------------------------------------

specrad <- as.data.frame(matrix(NA,ncol=length(networks),nrow=1))
colnames(specrad) <- networks

for ( i in 1:length(networks) ){
  nw <- eval(parse(text=networks[i]))
  
  zeroblock  <- matrix(0,nrow=nrow(nw),ncol=ncol(nw))
  A          <- cbind(rbind(zeroblock,nw), rbind(t(nw),zeroblock))
  specrad[i] <- max(eigen(A)$values)
  
}
print(specrad)

upperbound <- sqrt(sum(nw>0))
print(upperbound)

if ( SAVE ){
  setwd(save.dir)
  filename = paste(save.dir, '/spectral.csv', sep='')
  write.table(specrad, filename)
}

# Topological (binary) complementarity ----------------------------------------------------------------

# Implementing Tylianakis' complementarity measure, Fcomp (is identical to Petchey's Trophic Diversity if made only on the predation matrix rather than community matrix, independent sum of predation and rsource matrix, which is what Petchey et al actually used...)
topo.comp <- as.data.frame(matrix(NA,ncol=length(networks),nrow=1))
colnames(topo.comp) <- networks
rownames(topo.comp) <- 'Fcomp'

for ( i in 1:length(networks) ){
  nw <- eval(parse(text=networks[i]))
  
  ppm <- nw
  
  distance = vegdist(t(ppm), method="euclidean", bin=F)
  distance[is.nan(distance)] <- 0
  
  tree <- hclust(distance)
  #plot(tree)    
  topo.comp[,i] <- treeheight(tree)
}
print(topo.comp)

if ( SAVE ){
  setwd(save.dir)
  filename = paste(save.dir, '/complementarity.csv', sep='')
  write.table(topo.comp, filename)
}

