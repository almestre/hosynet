#################
# null.weights.geo.R
#################

#################
## DESCRIPTION: function to calculate a matrix of weighted probabilities of trait frequencies based on geographic availability of species in the target community.

## Methods:
# The function uses a community matrix and a vector of geographic modules. For each community (i.e. row of the community matrix), it obtains the set of species that belong to the geographic modules represented in the target community. Then, it calculates the frequencies of trait variants in the selected species and store them as a row of the output matrix.

#################
## Created by Alexandre Mestre
## Date of creation: 20/03/2020

#################
## Parameters
# com.matrix: MxN matrix where M (rows) are community entities (e.g. interaction modules), and N (columns) are species.
# trait: vector with trait variants of species composing the global pool (same order than columns of com.matrix).
# geo_partition: vector with the geographic module identities off all the species composing the global pool (same order than columns of com.matrix).

#################
##### Output
# Matrix with the weighted probabilities, where rows are community entities and columns are trait variants.

null.weights.geo <- function ( com.matrix , trait, geo_partition )
{
   traitID <- sort(unique(trait))
   Ntrait <- length(traitID)

   trait.freq <- function (traitSp , traitFreq, trait.id = traitID )
      {
         selec <- trait.id == traitSp
         traitFreq[ selec ] <- traitFreq[ selec ] + 1
         return(traitFreq)
      }
   
   weights <- sapply (
                        1:nrow(com.matrix),
                        function (i)
                           {   # Identify geographic modules of the target community
                               selec <- com.matrix[i,] > 0
                               geo.mod <- geo_partition[selec]

                               # Select trait variant of all the species that belong to the geo. mods
                               selec2 <- geo_partition %in% geo.mod
                               tr <- trait[selec2]

                               # Calculate trait frequencies
                               tr.Freq <- rep(0,Ntrait)
                               names(tr.Freq) <- traitID
                               for (j in 1:length(tr))
                                  {
                                     tr.Freq <- trait.freq(traitSp = tr[j], traitFreq = tr.Freq)
                                  }

                               return(tr.Freq)                               
                           }
                      )

   weights <- t(weights)
   colnames(weights) <-  traitID
   rownames(weights) <- rownames(com.matrix)

   return(weights)
}

