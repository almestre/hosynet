#################
# H.entropy.R
#################

#################
## DESCRIPTION: function to calculate the H index of entropy of a categorical trait in a set of communities.

## Methods:
# Shannon entropy (H):
# H = -1*sum(pi*ln(pi)) where pi is the proportion of occurrences of the type "i" of the categorical trait in the set of species composing the target community.

#################
## Created by Alexandre Mestre
## Date of creation: 19/03/2020

#################
## Parameters
# matrix: MxN matrix where M (rows) are community entities (e.g. interaction modules), and N (columns) are types of a categorical variable (e.g. habitat types, genera, geographic modules...). Each row stores the frequency of each variable type in species within the target community.

#################
##### Output
# A vector with the H indices

H.entropy <- function(matrix)
{
   totals <- apply(matrix, MARGIN = 1, sum)
   pi <- matrix/totals
   lnpi <- log(pi)
   pi.lnpi <- pi*lnpi
   h.obs <- -1 * apply( pi.lnpi, MARGIN = 1, sum, na.rm = TRUE )   

   return(h.obs)
}
