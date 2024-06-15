#This script considers categorical match scores and multiple classifiers (ML/AI/humans/etc.)
source("sim.pairwise.multiCat.R")
library(nimble) #need nimble for data simulator, loading here

#make state space first to get area
buff <- 3 #state space buffer. Should be at least 3 sigma.
X <- as.matrix(expand.grid(3:11,3:11))
area <- get.area(X,buff)

#D and area determine expected N
D <- 0.3
lambda <- D*area #expected N
lambda

lam0 <- 0.25
sigma <- 0.50
K <- 5
n.classifiers <- 3
# 2 category case - e.g., correct/incorrect match classifications from humans
# pi.match <- matrix(c(0.95,0.05, #match parameters, P(observe 1 or 2|match)
#                     0.1,0.9), #nonmatch parameters, P(observe 1 or 2|nonmatch)
#                   byrow=TRUE,ncol=2)

#3 category case - e.g., correct/unsure/incorrect match classifications from humans
pi.match <- matrix(c(0.8,0.15,0.05, #match parameters, P(observe 1, 2, or 3|match)
                    0.1,0.2,0.7), #nonmatch parameters, P(observe 1, 2, or 3|nonmatch)
byrow=TRUE,ncol=3)
n.cat <- ncol(pi.match)

data <- sim.pairwise.multiCat(D=D,area=area,lam0=lam0,sigma=sigma,K=K,X=X,buff=buff,pi.match=pi.match,
                                 n.classifiers=n.classifiers)

data$N #realized N
data$n #number of inds captured
table(rowSums(data$y)) #number of inds captures X times
str(data$scores) #upper triangle matrix of scores
data$scores[,1,] #scores for sample 1 over classifiers and samples.

#the only data are the pairwise sample scores, the trap of each sample, and the occasion of each sample
str(data$scores) #n.classifier x n.samples x n.samples scores
str(data$this.j) #n.samples length vector of trap IDs
str(data$this.k) #n.samples length vector of occasion IDs

#summarize match scores
#just looking at 1st classifier here
matches <- data$scores[1,,][data$match.type==1]
matches <- matches[-which(is.na(matches))]
nonmatches <- data$scores[1,,][data$match.type==2]
nonmatches <- nonmatches[-which(is.na(nonmatches))]
table(matches)
table(nonmatches)

##Fit Model
library(nimble)
library(coda)
nimbleOptions(determinePredictiveNodesInModel = FALSE)
source("init.data.pairwise.multiCat.R")
source("Nimble Functions pairwise multiCat.R") #nimble functions used in data simulator
source("NimbleModel pairwise multiCat.R")
source("sSampler.R")

M <- 150 #data augmentation limit. Must be larger than simulated N. If N posterior hits M, need to raise M and try again.
if(M<=data$N)stop("Raise M to be larger than N.")

inits <- list(sigma=0.5) #needs to be set somewhere in the ballpark of truth

#not using match scores to initialize data, just spatial location and sigma init
nimbuild <- init.data.pairwise.multiCat(data=data,inits=inits,M=M)
n.samples <- nimbuild$n.samples
J <- nrow(data$X)
K1D <- rep(K,J)

Niminits <- list(z=nimbuild$z,s=nimbuild$s,
                 ID=nimbuild$ID,capcounts=rowSums(nimbuild$y.true),
                 y.true=nimbuild$y.true,N=nimbuild$N,
                 sigma=inits$sigma,lam0=lam0)

#constants for Nimble
constants <- list(M=M,n.samples=n.samples,K1D=K1D,J=J,area=area,X=X,n.classifiers=n.classifiers,
                  xlim=data$xlim,ylim=data$ylim,n.cat=n.cat)

#supply data to nimble
Nimdata <- list(y.true=matrix(NA,nrow=M,ncol=J),scores=data$scores)

# set parameters to monitor
parameters <- c('lam0','pi.match','sigma','N','n','lambda',"D")

#can also monitor a different set of parameters with a different thinning rate
nt <- 2 #thinning rate

# Build the model, configure the mcmc, and compile
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,inits=Niminits)
#tell nimble which nodes to configure so we don't waste time for samplers we will replace below
config.nodes <- c("D","lam0","sigma",'pi.match')
# config.nodes <- c()

conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt,useConjugacy = FALSE,
                      nodes=config.nodes) 

##Here, we remove the default sampler for y.true
#and replace it with the custom "IDSampler".
# conf$removeSampler("y.true")
# conf$removeSampler("scores")
calcNodes.y.true <- Rmodel$expandNodeNames("y.true")
calcNodes.scores <- Rmodel$expandNodeNames("scores")
#need to make a map to reference the match nodes all on the upper diagonal
tmp <- matrix(NA,n.samples,n.samples)
idx <- 1
for(j in 2:n.samples){
  for(i in 1:(j-1)){
    tmp[i,j] <- idx
    idx <- idx + 1
  }
}
map.match <- matrix(NA,n.samples,n.samples-1)
map.match[1,] <- tmp[1,-1]
for(i in 2:n.samples){
  tmp2 <- c(tmp[,i],tmp[i,])
  map.match[i,] <- tmp2[!is.na(tmp2)]
}
#
tmp <- array(NA,dim=c(n.classifiers,n.samples,n.samples))
idx <- 1
for(j in 2:n.samples){
  for(i in 1:(j-1)){
    for(o in 1:n.classifiers){
      tmp[o,i,j] <- idx
      idx <- idx + 1
    }
  }
}
map.score <- array(NA,dim=c(n.classifiers,n.samples,n.samples-1))
for(o in 1:n.classifiers){
  map.score[o,1,] <- tmp[o,1,-1]
  for(i in 2:n.samples){
    tmp2 <- c(tmp[o,,i],tmp[o,i,])
    map.score[o,i,] <- tmp2[!is.na(tmp2)]
  }
}

calcNodes.match.type <- Rmodel$expandNodeNames("match.type")
calcNodes.all <- c(calcNodes.y.true,calcNodes.match.type,calcNodes.scores)
conf$addSampler(target = paste0("y.true[1:",M,",1:",K,"]"),
                type = 'IDSampler',control = list(M=M,J=J,n.samples=nimbuild$n.samples,
                                                  map.score=map.score,
                                                  map.match=map.match,
                                                  this.j=data$this.j,
                                                  calcNodes.y.true=calcNodes.y.true,
                                                  calcNodes.match.type=calcNodes.match.type,
                                                  calcNodes.scores=calcNodes.scores,
                                                  calcNodes.all=calcNodes.all),
                silent = TRUE)

###*required* sampler replacement for "alternative data augmentation" N/z update
z.ups <- round(M*0.5) # how many N/z proposals per iteration? Not sure what is optimal, setting to 25% of M here.
# conf$removeSampler("N")
#nodes used for update, calcNodes + z nodes
y.nodes <- Rmodel$expandNodeNames(paste("y.true[1:",M,",1:",J,"]"))
lam.nodes <- Rmodel$expandNodeNames(paste("lam[1:",M,",1:",J,"]"))
N.node <- Rmodel$expandNodeNames(paste("N"))
z.nodes <- Rmodel$expandNodeNames(paste("z[1:",M,"]"))
calcNodes <- c(N.node,lam.nodes,y.nodes)
conf$addSampler(target = c("N"),
                type = 'zSampler',control = list(z.ups=z.ups,M=M,
                                                 y.nodes=y.nodes,lam.nodes=lam.nodes,
                                                 N.node=N.node,z.nodes=z.nodes,
                                                 calcNodes=calcNodes),
                silent = TRUE)

# replace default activity center sampler
# conf$removeSampler(paste("s[1:",M,", 1:2]", sep=""))
for(i in 1:M){
  conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
                  type = 'sSampler',control=list(i=i,xlim=data$xlim,ylim=data$ylim),silent = TRUE)
}

# conf$removeSampler(c("lam0","sigma"))
conf$addSampler(target = c("lam0","sigma"),
                type = 'RW_block',control=list(adaptive=TRUE,tries=1),silent = TRUE)

# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=1) #this will run in R, used for better debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run the model.
start.time2<-Sys.time()
Cmcmc$run(2000,reset=FALSE) #can keep running this line to extend run
end.time<-Sys.time()
end.time-start.time  # total time for compilation, replacing samplers, and fitting
end.time-start.time2 # post-compilation run time

library(coda)
mvSamples = as.matrix(Cmcmc$mvSamples)
plot(mcmc(mvSamples[500:nrow(mvSamples),])) #discarding some burnin here. Can't plot 1st sample which is all NA

lambda #target expected abundance
data$N #target realized abundance
data$n #target number detected (n)
