NimModel <- nimbleCode({
  #--------------------------------------------------------------
  # priors
  #--------------------------------------------------------------
  D ~ dunif(0,10) #expected D
  lam0 ~ dunif(0,5)
  sigma ~ dunif(0,10)
  #match parameters. Assuming no variation across classifiers, but can modify code below to allow that.
  for(m in 1:2){ #m=1 is correct score, m=2 is incorrect score
    for(c in 1:n.cat){
      alpha.match[m,c] <- 1 #dirichlet prior parameter. Can use prior distribution as well, with support in (0,Inf)
    }
    pi.match[m,1:n.cat] ~ ddirch(alpha.match[m,1:n.cat])
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
      for(o in 1:n.classifiers){
        scores[o,i,j] ~ dcat(pi.match[match.type[i,j],1:n.cat])
      }
    }
  }
  #latent ID derived variables
  #number of detections per individual
  capcounts[1:M] <- Getcapcounts(y.true=y.true[1:M,1:J])
  #number of detected individuals
  n <- Getncap(capcounts=capcounts[1:M],ID=ID[1:n.samples])
})# end model