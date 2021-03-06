# Initialise HMM
hmm = initHMM(c("A","B"), c("L","R"), transProbs=matrix(c(.8,.2,.2,.8),2),
              emissionProbs=matrix(c(.6,.4,.4,.6),2))
print(hmm)
# Sequence of observations
observations = c("L","L","R","R")
# Calculate backward probablities
logBackwardProbabilities = backward(hmm,observations)
print(exp(logBackwardProbabilities))

# Calculate forward probablities
observations = c("L","L","R","R")
logForwardProbabilities = forward(hmm,observations)
print(exp(logForwardProbabilities))

# Calculate posterior probablities of the states
posterior = posterior(hmm,observations)
print(posterior)

# Calculate Viterbi path
viterbi = viterbi(hmm,observations)
print(viterbi)