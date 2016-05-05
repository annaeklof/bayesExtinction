### Code for analyzing extinctions in food webs using Bayesian networks. Code is written by Anna Eklöf, Si Tang and Stefano Allesina. 
  # Checked 2014-04-30 for R version 3.0.2.
  # Checked 2016-05-05 for R version 3.2.4.

### REFERENCE: Eklöf, Anna, Si Tang, and Stefano Allesina. "Secondary extinctions in food webs: a Bayesian network approach." Methods in Ecology and Evolution 4.8 (2013): 760-770.

### Please cite and acknowledge properly!

source("BN-BayesNet.R")
source("BN-Likelihoods.R")
source("BN-DFS.R")

LaunchAnalysis <- function(Adjacency,
                           Extant,
                           ProbExtinction,
                           Functional = "topo",
                           Alpha = 0.0,
                           Beta = 0.0,
                           Label){
    ## Wrapper for the analysis.
    ##############################
    ## INPUT
    ##############################
    # Adjacency: square, binary matrix
    #            Adjacency[i,j] = 1 means that j consumes i
    # Extant: observed presence of the species
    #         binary, Reps x S matrix, where
    #         Reps is the number of observations
    #         and each row is a different replicate
    #         marking whether the species was present (1) or absent (0)
    # ProbExtinction: vector of length S, measuring
    #                 the probability that each species goes extinct due
    #                 to external causes
    # Functional: functional form, as in the paper.
    #             Options are "topo", "linear", and "nonlinear"
    # Alpha, Beta: parameters for nonlinear functional response
    # Label: label for the system
    #
    # Results returned as "results.Rdata" stored in "FW"
    ##############################

    # First, we need to make the adjacency matrix acyclic
    # For this, we use BN-DFS.R
    A <- DFS(A=Adjacency)
    if (sum(A) < sum(Adjacency)){
        print("The network is not acyclic, and thus some links were removed!")
    }
    ## Call the function
    FW <- GetBayesNetMarginals(A,
                              Extant,
                              ProbExtinction,
                              Functional,
                              Alpha,
                              Beta,
                              Label)
    FW$Functional <- NULL
    FW$OriginalAdjacency <- Adjacency
    FW$RemovedLinks <- Adjacency - A
    
    save(FW, file='results.RData')
    return(FW)
}

# Run the code by writing TestCode()
TestCode <- function(){
    # Number of species
    S <- 10
    # Construct an adjacency matrix with zeros on the diagonal
    Adjacency <- (matrix(runif(S*S), S, S) < 0.5) * 1
    diag(Adjacency) <- 0
    # Define number of replicates
    Reps <- 50
    # Extant (see above) 
    Extant <- (matrix(runif(Reps * S), Reps, S) < 0.2) * 1

    print(LaunchAnalysis(Adjacency, Extant, rep(0.1, S), "nonlinear", 0.5, 3, "FW1"))
    
}
