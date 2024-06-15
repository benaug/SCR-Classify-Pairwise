NimModel <- nimbleCode({
  #--------------------------------------------------------------
  # priors
  #--------------------------------------------------------------
  D ~ dunif(0,10) #expected D
  lam0 ~ dunif(0,5) #baseline detection rate
  sigma ~ dunif(0,10) #detection spatial scale
  for(i in 1:2){
    lambda.match[i] ~ dunif(0,50) #match score|match type (correct,incorrect)
  }
  #--------------------------------------------------------------
  lambda <- D*area #expected N
  N ~ dpois(lambda) #realized N
  for(i in 1:M){
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    lam[i,1:J] <- GetDetectionRate(s = s[i,1:2], X = X[1:J,1:2], J=J,sigma=sigma, lam0=lam0, z=z[i])
    y.true[i,1:J] ~ dPoissonVector(lam=lam[i,1:J]*K1D[1:J],z=z[i]) #vectorized obs mod
  }
  #Pairwise sample-level match scores - upper triangle matrix
  for(i in 1:(n.samples-1)){
    for(j in (i+1):n.samples){
      match.type[i,j] <- 1*equals(ID[i],ID[j]) + 2*(1-equals(ID[i],ID[j]))
      scores[i,j] ~ dpois(lambda.match[match.type[i,j]])
    }
  }
  #latent ID derived variables
  #number of detections per individual
  capcounts[1:M] <- Getcapcounts(y.true=y.true[1:M,1:J])
  #number of detected individuals
  n <- Getncap(capcounts=capcounts[1:M],ID=ID[1:n.samples])
})# end model