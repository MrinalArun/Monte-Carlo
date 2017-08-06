# Parts a and b

Pages = 6
Q = matrix(0,nrow = Pages, ncol = Pages)
Q[1,2] = 0.5
Q[1,3] = 0.5
Q[2,3] = 1
Q[3,1] = 1
Q[4,1] = 0.5
Q[4,3] = 0.5
Q[5,3] = 1
Q[6,2] = 0.5
Q[6,5] = 0.5

Eps = 0.15

B = matrix(1/Pages,nrow = Pages, ncol = Pages)
A = (1-Eps)*Q + Eps*B
Mu  = matrix(1/Pages, nrow = 1, ncol = Pages)

for(i in 1:1000){
  Mu = Mu %*% A
}

print(Mu)

# Part c - Testing it out for different values of Epsilon

RankMat = matrix(0,nrow = 20, ncol = Pages)
for(j in 1:20){
  Eps = j*0.05
  A   = (1-Eps)*Q + Eps*B
  Mu  = matrix(1/Pages, nrow = 1, ncol = Pages)
  for(i in 1:1000){
    Mu = Mu %*% A
  }
  RankMat[j,] = Mu
}

print(RankMat)

Epsilons = c(1:20)*0.05
plot(Epsilons,RankMat[,1],type = "l",col = 1,ylim = c(0,0.5))

for(i in 2:6){
  lines(Epsilons,RankMat[,i],col = i)
}