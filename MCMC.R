Trials = 100000

PFunc = function(p){
  if(p == 1){
    return(0.2)
  }else{
    return(0.8)
  }
}

DFunc = function(d){
  if(d == 1){
    return(0.4)
  }else{
    return(0.6)
  }
}

HFunc = function(p,h){
  if(h == 1 && p == 1){
    answer = 0.9
  }else if(h == 1 && p == 0){
    answer = 0.2
  }else if(h == 0 && p == 1){
    answer = 0.1
  }else if(h == 0 && p == 0){
    answer = 0.8
  }
  return(answer)
}

AFunc = function(u,a){
  if(u == 1 && a == 1){
    answer = 0.95
  }else if(u == 1 && a == 0){
    answer = 0.05
  }else if(u == 0 && a == 1){
    answer = 0.5
  }else if(u == 0 && a == 0){
    answer = 0.5
  }
  return(answer)
}

UFunc = function(p,d,u){
  if(p == 1 && d == 1 && u == 1){
    answer = 0.999
  }else if(p == 1 && d == 1 && u == 0){
    answer = 0.001
  }else if(p == 1 && d == 0 && u == 1){
    answer = 0.9
  }else if(p == 1 && d == 0 && u == 0){
    answer = 0.1
  }else if(p == 0 && d == 1 && u == 1){
    answer = 0.9
  }else if(p == 0 && d == 1 && u == 0){
    answer = 0.1
  }else if(p == 0 && d == 0 && u == 1){
    answer = 0.01
  }else if(p == 0 && d == 0 && u == 0){
    answer = 0.99
  }
  return(answer)
}

ProbFunc = function(p,d,u,h,a){
  Answer = PFunc(p)*DFunc(d)*UFunc(p,d,u)*HFunc(p,h)*AFunc(u,a)
  return(Answer)
}

# Doing MCMC
MCMCFunc = function(Start){
  Matrix = matrix(0,nrow = Trials, ncol = 3)

  for(i in 1:Trials){
    rand1  = runif(1)
    rand2  = runif(1)
    if(rand1 < 1/3){
      if(rand2 < min(ProbFunc(1-Start[1],Start[2],Start[3],1,1)/ProbFunc(Start[1],Start[2],Start[3],1,1),1)){
        Start[1] = 1-Start[1]
      }
    }else if(rand1<2/3){
      if(rand2 < min(ProbFunc(Start[1],1-Start[2],Start[3],1,1)/ProbFunc(Start[1],Start[2],Start[3],1,1),1)){
        Start[2] = 1-Start[2]
      }
    }else{
      if(rand2 < min(ProbFunc(Start[1],Start[2],1-Start[3],1,1)/ProbFunc(Start[1],Start[2],Start[3],1,1),1)){
        Start[3] = 1-Start[3]
      }
    }
    Matrix[i,] = Start
  }
  return(Matrix[,1])
}

MCMC     = matrix(0,nrow = Trials,ncol = 8)
MCMC[,1] = MCMCFunc(c(0,0,0))
MCMC[,2] = MCMCFunc(c(0,0,1))
MCMC[,3] = MCMCFunc(c(0,1,0))
MCMC[,4] = MCMCFunc(c(0,1,1))
MCMC[,5] = MCMCFunc(c(1,0,0))
MCMC[,6] = MCMCFunc(c(1,0,1))
MCMC[,7] = MCMCFunc(c(1,1,0))
MCMC[,8] = MCMCFunc(c(1,1,1))

MCMC  = MCMC[((Trials/2)+1):Trials,]
CMean = vector()
Temp  = vector()
for(i in 1:8){
  CMean = c(CMean,mean(MCMC[,i]))
  Temp  = c(Temp,var(MCMC[,i]))
}
MCMCMean = mean(CMean)
B        = 8*var(CMean)
W        = mean(Temp)
Var      = ((Trials - 1)*W + B)/Trials
R        = sqrt(Var/W)

print(MCMCMean)
print(R)

# Doing Gibbs Sampling
GibbsFunc = function(Start){
  Matrix = matrix(0,nrow = Trials, ncol = 3)
  
  for(i in 1:Trials){
    
    rand  = runif(3)
    
    if(rand[1] < ProbFunc(1,Start[2],Start[3],1,1)/(ProbFunc(1,Start[2],Start[3],1,1)+(ProbFunc(0,Start[2],Start[3],1,1)))){
      Start[1] = 1
    }else{
      Start[1] = 0
    }
    
    if(rand[2] < ProbFunc(Start[1],1,Start[3],1,1)/(ProbFunc(Start[1],1,Start[3],1,1)+(ProbFunc(Start[1],0,Start[3],1,1)))){
      Start[2] = 1
    }else{
      Start[2] = 0
    }
    
    if(rand[3] < ProbFunc(Start[1],Start[2],1,1,1)/(ProbFunc(Start[1],Start[2],1,1,1)+(ProbFunc(Start[1],Start[2],0,1,1)))){
      Start[3] = 1
    }else{
      Start[3] = 0
    }
    
    Matrix[i,] = Start
  }
  return(Matrix[,1])
}

Gibbs     = matrix(0,nrow = Trials,ncol = 8)
Gibbs[,1] = GibbsFunc(c(0,0,0))
Gibbs[,2] = GibbsFunc(c(0,0,1))
Gibbs[,3] = GibbsFunc(c(0,1,0))
Gibbs[,4] = GibbsFunc(c(0,1,1))
Gibbs[,5] = GibbsFunc(c(1,0,0))
Gibbs[,6] = GibbsFunc(c(1,0,1))
Gibbs[,7] = GibbsFunc(c(1,1,0))
Gibbs[,8] = GibbsFunc(c(1,1,1))

Gibbs = Gibbs[((Trials/2)+1):Trials,]
CMean = vector()
Temp  = vector()
for(i in 1:8){
  CMean = c(CMean,mean(Gibbs[,i]))
  Temp  = c(Temp,var(Gibbs[,i]))
}
GibbsMean = mean(CMean)
B        = 8*var(CMean)
W        = mean(Temp)
Var      = ((Trials - 1)*W + B)/Trials
R        = sqrt(Var/W)

print(GibbsMean)
print(R)
