getCellR = function(u,res,cells,xlim,ylim){
  inout=1*(u[1]>xlim[1]&u[1]<xlim[2]&u[2]>ylim[1]&u[2]<ylim[2])
  if(inout==1){
    this.cell=cells[trunc(u[1]/res)+1,trunc(u[2]/res)+1]
  }else{
    this.cell=0
  }
  return(this.cell)
}

init.data.SCR.pairwise <- function(data=NA,M=NA,inits=NA){
  library(abind)
  this.j <- data$this.j
  this.k <- data$this.k
  X <- as.matrix(data$X)
  J <- nrow(X)
  K <- data$K
  n.samples <- length(this.j)
  
  #state space extent
  buff <- data$buff
  xlim <- c(min(X[,1]),max(X[,1]))+c(-buff, buff)
  ylim <- c(min(X[,2]),max(X[,2]))+c(-buff, buff)
  
  # lam0 <- inits$lam0
  sigma <- inits$sigma
  
  #assign random activity centers
  s <- cbind(runif(M,xlim[1],xlim[2]), runif(M,ylim[1],ylim[2]))
  D <- e2dist(s, X)
  
  # lamd <- lam0*exp(-D*D/(2*sigma*sigma))
  lamd <- exp(-D*D/(2*sigma*sigma))
  #Build y.true
  y.true <- array(0,dim=c(M,J,K))
  ID <- rep(NA,n.samples)
  for(l in 1:n.samples){
    propdist <- lamd[,this.j[l]]
    propdist <- propdist/sum(propdist)
    ID[l] <- sample(1:M,1,replace=FALSE,prob=propdist)
    y.true[ID[l],this.j[l],this.k[l]] <- y.true[ID[l],this.j[l],this.k[l]]+1
  }
  z <- 1*(rowSums(y.true)>0)
  N <- sum(z) #must use these z and N inits
  
  y.true2D <- apply(y.true,c(1,2),sum)
  
  return(list(y.true=y.true2D,y.true3D=y.true,z=z,s=s,ID=ID,n.samples=n.samples,xlim=xlim,ylim=ylim,
              this.k=this.k,N=N))
}