GetDetectionRate <- nimbleFunction(
  run = function(s = double(1), lam0=double(0), sigma=double(0), 
                 X=double(2), J=double(0), z=double(0)){ 
    returnType(double(1))
    if(z==0) return(rep(0,J))
    if(z==1){
      d2 <- ((s[1]-X[1:J,1])^2 + (s[2]-X[1:J,2])^2)
      ans <- lam0*exp(-d2/(2*sigma^2))
      return(ans)
    }
  }
)

dPoissonVector <- nimbleFunction(
  run = function(x = double(1), lam = double(1), z = double(0),
                 log = integer(0)) {
    returnType(double(0))
    if(z==0){
      if(sum(x)>0){ #need this so z is not turned off if samples allocated to individual
        return(-Inf)
      }else{
        return(0)
      }
    }else{
      logProb <- sum(dpois(x, lambda = lam, log = TRUE))
      return(logProb)
    }
  }
)

#make dummy random vector generator to make nimble happy
rPoissonVector <- nimbleFunction(
  run = function(n = integer(0),lam = double(1),z = double(0)) {
    returnType(double(1))
    J=nimDim(lam)[1]
    out=numeric(J,value=0)
    return(out)
  }
)

Getcapcounts <- nimbleFunction(
  run = function(y.true=double(2)){
    returnType(double(1))
    M <- nimDim(y.true)[1]
    J <- nimDim(y.true)[2]
    capcounts=numeric(M, value = 0)
    for(i in 1:M){
      capcounts[i]=sum(y.true[i,1:J])
    }
    return(capcounts)
  }
)
Getncap <- nimbleFunction(
  run = function(capcounts=double(1),ID=double(1)){ #don't need ID, but nimble requires is it used in a function 
    returnType(double(0))
    M <- nimDim(capcounts)[1]
    nstate <- numeric(M, value = 0)
    for(i in 1:M){
      if(capcounts[i]>0){
        nstate[i] <- 1
      }
    }
    n.cap <- sum(nstate)
    return(n.cap)
  }
)

#Required custom update for N/z
zSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    M <- control$M
    z.ups <- control$z.ups
    y.nodes <- control$y.nodes
    lam.nodes <- control$lam.nodes
    N.node <- control$N.node
    z.nodes <- control$z.nodes
    calcNodes <- control$calcNodes
  },
  run = function() {
    for(up in 1:z.ups){ #how many updates per iteration?
      #propose to add/subtract 1
      updown <- rbinom(1,1,0.5) #p=0.5 is symmetric. If you change this, must account for asymmetric proposal
      reject <- FALSE #we auto reject if you select a detected call
      if(updown==0){#subtract
        # find all z's currently on
        z.on <- which(model$z==1)
        n.z.on <- length(z.on)
        pick <- rcat(1,rep(1/n.z.on,n.z.on)) #select one of these individuals
        pick <- z.on[pick]
        #prereject turning off individuals currently allocated samples
        if(model$capcounts[pick]>0){#is this an individual with samples?
          reject <- TRUE
        }
        if(!reject){
          #get initial logprobs for N and y
          lp.initial.N <- model$getLogProb(N.node)
          lp.initial.y <- model$getLogProb(y.nodes[pick])

          #propose new N/z
          model$N[1] <<-  model$N[1] - 1
          model$z[pick] <<- 0
          
          #turn lam off
          model$calculate(lam.nodes[pick])

          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.node)
          lp.proposed.y <- model$calculate(y.nodes[pick]) #will always be 0

          #MH step
          log_MH_ratio <- (lp.proposed.N + lp.proposed.y) - (lp.initial.N + lp.initial.y)
          accept <- decide(log_MH_ratio)

          if(accept) {
            mvSaved["N",1][1] <<- model[["N"]]
            mvSaved["lam",1][pick,] <<- model[["lam"]][pick,]
            mvSaved["z",1][pick] <<- model[["z"]][pick]
          }else{
            model[["N"]] <<- mvSaved["N",1][1]
            model[["lam"]][pick,] <<- mvSaved["lam",1][pick,]
            model[["z"]][pick] <<- mvSaved["z",1][pick]
            model$calculate(y.nodes[pick])
            model$calculate(N.node)
          }
        }
      }else{#add
        if(model$N[1] < M){ #cannot update if z maxed out. Need to raise M
          z.off <- which(model$z==0)
          n.z.off <- length(z.off)
          pick <- rcat(1,rep(1/n.z.off,n.z.off)) #select one of these individuals
          pick <- z.off[pick]
          #get initial logprobs for N and y
          lp.initial.N <- model$getLogProb(N.node)
          lp.initial.y <- model$getLogProb(y.nodes[pick]) #will always be 0
          
          #propose new N/z
          model$N[1] <<-  model$N[1] + 1
          model$z[pick] <<- 1
          
          #turn lam on
          model$calculate(lam.nodes[pick])
          
          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.node)
          lp.proposed.y <- model$calculate(y.nodes[pick])
          
          #MH step
          log_MH_ratio <- (lp.proposed.N + lp.proposed.y) - (lp.initial.N + lp.initial.y)
          accept <- decide(log_MH_ratio)
          if(accept) {
            mvSaved["N",1][1] <<- model[["N"]]
            mvSaved["lam",1][pick,] <<- model[["lam"]][pick,]
            mvSaved["z",1][pick] <<- model[["z"]][pick]
          }else{
            model[["N"]] <<- mvSaved["N",1][1]
            model[["lam"]][pick,] <<- mvSaved["lam",1][pick,]
            model[["z"]][pick] <<- mvSaved["z",1][pick]
            model$calculate(y.nodes[pick])
            model$calculate(N.node)
          }
        }
      }
    }
    #copy back to mySaved to update logProbs which was not done above
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)

IDSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    M <-control$M
    J <- control$J
    n.samples <- control$n.samples
    this.j <- control$this.j
    map.score <- control$map.score
    calcNodes.y.true <- control$calcNodes.y.true
    calcNodes.match.type <- control$calcNodes.match.type
    calcNodes.scores <- control$calcNodes.scores
    calcNodes.all <- control$calcNodes.all
  },
  run = function() {
    s <- model$s
    z <- model$z
    sigma <- model$sigma[1]
    
    for(l in 1:n.samples){
      propprobs <- model$lam[,this.j[l]]*z
      propprobs <- propprobs/sum(propprobs)
      pick <- rcat(1,prob=propprobs)
      if(model$ID[l]!=pick){ #skip if propose same ID.
        swapped <- c(mvSaved["ID",1][l],pick) #order swap.out then swap.in
        propprob <- propprobs[swapped[2]]
        backprob <- propprobs[swapped[1]]
        #focal column and row. diagonal element counted twice, but doesn't matter because diagonal
        #elements (self match) likelihood never changes.
        these.scores <- map.score[l,]
        lp.initial.y.true <- model$getLogProb(calcNodes.y.true[swapped])
        lp.initial.scores <- model$getLogProb(calcNodes.scores[these.scores])
        lp.initial.total <- lp.initial.y.true + lp.initial.scores
        model[["ID"]][l] <<- pick
        #p(select sample from this individual at this trap on this occasion)
        focalprob <- sum(mvSaved["ID",1]==mvSaved["ID",1][l]&this.j==this.j[l])/n.samples
        focalbackprob <- sum(model[["ID"]]==model[["ID"]][l]&this.j==this.j[l])/n.samples
        
        #forgive me, but this is just easier for me to track here
        ID.curr <- swapped[1]
        ID.cand <- swapped[2]
        
        model[["y.true"]][ID.curr,this.j[l]] <<- model[["y.true"]][ID.curr,this.j[l]] - 1
        model[["y.true"]][ID.cand,this.j[l]] <<- model[["y.true"]][ID.cand,this.j[l]] + 1

        #need to update match.type before scores logprob
        model$calculate(calcNodes.match.type[these.scores])

        lp.proposed.y.true <- model$calculate(calcNodes.y.true[swapped])
        lp.proposed.scores <- model$calculate(calcNodes.scores[these.scores])
        lp.proposed.total <- lp.proposed.y.true + lp.proposed.scores
        
        #MH step
        log_MH_ratio <- (lp.proposed.total + log(backprob) + log(focalbackprob)) -
          (lp.initial.total + log(propprob) + log(focalprob))
        accept <- decide(log_MH_ratio)
        
        if(accept){
          #copy model to mvSaved - dont need to keep up with logprobs
          mvSaved["y.true",1][ID.curr,this.j[l]] <<- model[["y.true"]][ID.curr,this.j[l]]
          mvSaved["y.true",1][ID.cand,this.j[l]] <<- model[["y.true"]][ID.cand,this.j[l]]
          mvSaved["ID",1][l] <<- model[["ID"]][l]
          # mvSaved["match.type",1] <<- model[["match.type"]] #this works, but how to only move elements that changed?
          if(l<n.samples){ #last sample only has a column
            mvSaved["match.type",1][l,] <<- model[["match.type"]][l,]
          }
          mvSaved["match.type",1][,l] <<- model[["match.type"]][,l]
        }else{
          #set model back to mvSaved states
          model[["y.true"]][ID.curr,this.j[l]] <<- mvSaved["y.true",1][ID.curr,this.j[l]]
          model[["y.true"]][ID.cand,this.j[l]] <<- mvSaved["y.true",1][ID.cand,this.j[l]]
          #need to set logProbs back, not updated in mvSaved, so recalculate. not most efficient
          model$calculate(calcNodes.y.true[swapped])
          model[["ID"]][l] <<- mvSaved["ID",1][l]
          model$calculate(calcNodes.match.type[these.scores])
          model$calculate(calcNodes.scores[these.scores])
        }
      }
    }
    capcounts <- Getcapcounts(y.true=model$y.true[1:M,1:J])
    n <- Getncap(capcounts=capcounts,ID=model$ID)
    model$capcounts <<- capcounts
    model$n[1] <<- n
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes.all, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)
