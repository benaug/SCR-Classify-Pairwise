e2dist = function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

get.area = function (X, buff){
    area <- diff(range(X[,1])+c(-buff,buff))*diff(range(X[,2])+c(-buff,buff))
  return(area)
}

sim.pairwise.poisson <-
  function(D=NA,area=NA,lam0=NA,sigma=NA,K=NA,X=NA,buff=NA,lambda.match=NA,seed=NA){
    if(!is.na(seed)){
      set.seed(seed)
    }
    lambda <- D*area
    N <- rpois(1,lambda)
    # simulate a population of activity centers
    xlim <- c(min(X[,1]),max(X[,1]))+c(-buff, buff)
    ylim <- c(min(X[,2]),max(X[,2]))+c(-buff, buff)
    s <- cbind(runif(N, xlim[1],xlim[2]), runif(N,ylim[1],ylim[2]))
    D <- e2dist(s,X)
    lamd <- lam0*exp(-D*D/(2*sigma*sigma))
    J <- nrow(X)
    
    # Capture individuals
    y <- array(0,dim <- c(N,J,K))
    for(i in 1:N){
      for(j in 1:J){
        for(k in 1:K){
          y[i,j,k] <- rpois(1,lamd[i,j])
        }
      }
    } 
    
    #discard uncaptured inds
    caught <- which(apply(y,c(1),sum)>0)
    y.true <- y
    y <- y[caught,,]
    if(K==1){
      y <- array(y,dim=c(dim(y),1))
    }
    n <- length(caught)
    n.samples <- sum(y)
    
    #disaggregate samples
    ID <- this.j <- this.k <- rep(NA,n.samples)
    idx <- 1
    for(i in 1:n){ #loop through inds (uncaptured already removed)
      for(j in 1:J){ #then traps
        for(k in 1:K){ #then occasions
          if(y[i,j,k]>0){ #is there at least one sample here?
            for(l in 1:y[i,j,k]){ #then samples
              ID[idx] <- i #ID numbers don't count uncaptured guys
              this.j[idx] <- j
              this.k[idx] <- k
              idx <- idx+1
            }
          }
        }
      }
    }
    
    #get pairwise match scores
    match.type <- scores <- matrix(NA,n.samples,n.samples)
    for(i in 1:(n.samples-1)){
      for(j in (i+1):n.samples){
        match.type[i,j] <- 1*(ID[i]==ID[j]) + 2*(ID[i]!=ID[j])
        scores[i,j] <- rpois(1,lambda.match[match.type[i,j]])
      }
    }
    
    out<-list(this.j=this.j,this.k=this.k,scores=scores,match.type=match.type,ID=ID,n.samples=n.samples,
              X=X,K=K,buff=buff,K=K,N=N,n=n,y=y,xlim=xlim,ylim=ylim,seed=seed,area=area,lambda=lambda)
    return(out)
  }
