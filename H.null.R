#################
# H.null.R
#################

#################
## DESCRIPTION: function to randomize a community matrix of species-level categorical trait frequencies (e.g. habitat, geographic module or genus).

## Methods:
# The function randomizes rows keeping row sums and using weighted probabilities based on either global frequencies or community-level frequencies

#################
## Created by Alexandre Mestre
## Date of creation: 20/03/2020

#################
## Parameters
# matrix: MxN matrix where M (rows) are community entities (e.g. interaction modules), and N (columns) are types of a categorical trait (e.g. habitat types, genera, geographic modules...). Each row stores the frequency of each trait variant across species within the target community.
# pool: If pool = "global", weighted probabilities are based on global pool frequencies. If pool = "community", weighted probabilities are based on community-level frequencies.
# weights: frequencies used as weights in the randomisation of rows. If pool = "global", weights is a vector that applies to all the communities (i.e. rows).  If pool = "community", weights is a matrix with the same dimensions than the input matrix. In this case, each community has its own probability weights (e.g. based on geographic availability of species).

#################
##### Output
# The randomized matrix

H.null <- function ( matrix , pool = "global", weights )
{
   
totals <- apply(matrix, MARGIN = 1, sum)

   if (pool == "community")
      {
          Ntypes <- ncol(weights)
          types <- c(1:Ntypes)
          rand.matrix <- sapply(
                                1:nrow(matrix),
                                function (i)
                                   {
                                       samp <- sample(
                                                        x = types,
                                                        size = totals[i],
                                                        replace = TRUE,
                                                        prob = weights[i,]
                                                     )
                                       freq <- rep(0,Ntypes)
                                       freq <- sapply(1:Ntypes, function (i) {freq[i] <- sum(samp == i)})
                                       return(freq)
                                   }
                               )
         # We transpose the matrix because sapply store vectors by columns
         rand.matrix <- t(rand.matrix)
         colnames(rand.matrix) <- colnames(matrix)
         rownames(rand.matrix) <- rownames(matrix)
      }

   else
      {
          Ntypes <- length(weights)
          types <- c(1:Ntypes)
          rand.matrix <- sapply(
                                1:nrow(matrix),
                                function (i)
                                   {
                                       samp <- sample(
                                                        x = types,
                                                        size = totals[i],
                                                        replace = TRUE,
                                                        prob = weights
                                                     )
                                       freq <- rep(0,Ntypes)
                                       freq <- sapply(1:Ntypes, function (i) {freq[i] <- sum(samp == i)})
                                       return(freq)
                                   }
                               )
         # We transpose the matrix because sapply store vectors by columns
         rand.matrix <- t(rand.matrix)
         colnames(rand.matrix) <- colnames(matrix)
         rownames(rand.matrix) <- rownames(matrix)
      }

   return(rand.matrix)
}

