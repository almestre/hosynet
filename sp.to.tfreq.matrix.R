#################
# sp.to.tfreq.matrix
#################

#################
## Created by Alexandre Mestre
## Date of creation: 21/03/2020

#################
## DESCRIPTION: function to convert a community matrix with species occurrences into a community matrix with frequencies of variants a categorigcal trait e.g. (geographic module, genus...).

#################
## Parameters
# com.matrix: community matrix with MxN dimensions where M (rows) are community entities (e.g. interaction modules), and N (cols) are species.
# trait.values: vector storing trait variants of all the species of the global pool (same order than columns of com.matrix).

#################
## Output
# Community matrix with MxN dimensions where M (rows) are community entities (e.g. interaction modules), and N (cols) are trait variants. For each community, it stores the frequency of occurrence of a species trait variant in the set of species composing the community. 


sp.to.tfreq.matrix <- function( com.matrix , trait.values )
   {

    # First, we build a function to fill a given habitat row based on a species habitat preferences in the target community 
    f1 <- function ( traitSp , traitFreq )
       { # traitSp will be trait value of the target species
         # traitFreq will be the vector to store trait frequencies in the target community
         selec <- traitID == traitSp
         traitFreq[ selec ] <- traitFreq[ selec ] + 1
         return(traitFreq)
       }


    f2 <- function (comSp , trait = trait.values)
       {  
          # comSp will be a row in the input com.matrix
          # trait.sp will be the vector of trait values of species in the target community
          selec <- comSp > 0
          trait.sp <- trait[ selec ]
          trait.freq[] <- rep(0,Ntrait) # reset trait-frequency storage vector

          x <- sapply (1:length(trait.sp),
                    function (i)
                           {
                              trait.freq <- f1(
                                                 traitSp = trait.sp[i],
                                                 traitFreq = trait.freq
                                              )
                           }
                       )
          trait.freq <- rowSums(x)
          return(trait.freq)
       }

    traitID <- sort(unique(trait.values))
    Ntrait <- length(traitID)
    trait.freq <- rep(0,Ntrait) # vector to store trait frequencies in the target community
    names(trait.freq) <- traitID
    trait.matrix <- t(apply (com.matrix , MARGIN = 1, f2))
    colnames(trait.matrix) <- traitID
    rownames(trait.matrix) <- rownames(com.matrix)
    return(trait.matrix)
   }


