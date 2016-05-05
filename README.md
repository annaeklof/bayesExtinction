# bayesExtinction
Anna Eklöf, Si Tang and Stefano Allesina

# Overview
This repository includes the code for analyzing secondary extinctions in food webs using Bayesian networks. Code written by Anna Eklöf, Si Tang and Stefano Allesina.

REFERENCE: Eklöf, Anna, Si Tang, and Stefano Allesina. "Secondary extinctions in food webs: a Bayesian network approach." Methods in Ecology and Evolution 4.8 (2013): 760-770.

Please cite and acknowledge properly!

Contact: anna.eklof@liu.se

# Run the code
   
   BN-Wrapper.R is where the start function is located. Function to initialize analyse is TestCode
   
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
