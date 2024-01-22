#################
# rpd.R
#################

#################
## DESCRIPTION: function to calculate the RPD (relative phylogenetic distinctiveness) for each module of a modular interaction network (for a given partition).

## Methods:
# The RPD is defined in Calatayud (2016; they provide the equation):
# First, we calculate the phylogenetic mean pairwise distances (MPD) between the taxon i and all other taxa belonging to its module (i.e., MPD i intracommunity ) minus the mean pairwise distances of taxon i and all other taxa from different modules (MPD i intercommunity ). Then, we computed relative phylogenetic distinctiveness (RPD) as the inverse of the average between all taxa in a module.

#################
## Created by Alexandre Mestre
## Date of creation: 14/03/2020

#################
## Parameters
# dist_mat: Matrix of pairwise phylogenetic distances between species. The same function can be used to calculate relative functional distinctiveness if we use a distance matrix of interspecific differences in a functional trait.
# partition: vector storing the module indentity for each of the species. The vector should be labeled with the species names, because labels are used to order distance matrix according to the partition vector.

#################
## Output
# RPD: vector with the rpd values for the modules (length equal to the number of modules of the parittion)

rpd <- function(dist_mat, partition)
{
   # First, we order data based on species labels
   dist_mat <- dist_mat [order(rownames(dist_mat)) , order(colnames(dist_mat)) ]
   partition <- partition[ order(names(partition)) ]

   modID <- sort(unique(partition))
   Nmod <- length(modID)
   RPD <- rep(NA,Nmod)
   av <- function(x) {sum(x)/(length(x)-1)} # discount the pairwise with the same species (for mpd_intra)

   for (i in 1:Nmod)
      {
         # Select rows with species belonging to the targent module
         selec <- partition == modID[i]
         Nsp <- sum(selec) # number of species belonging to the target module

         # Calculate mpd_inter
         mpd_inter <- apply(dist_mat[selec, -selec, drop = FALSE], MARGIN = 1, mean)

         # Calculate mpd_intra
         mpd_intra <- apply(dist_mat[selec, selec, drop = FALSE], MARGIN = 1, av)
         
         # Calculate rpd
         RPD[i] <- -1*sum(mpd_intra - mpd_inter, na.rm = TRUE)/Nsp
      }

   return(RPD)
}

