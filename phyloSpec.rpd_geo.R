#################
# phyloSpec.rpd_geo.R
#################

#################
## DESCRIPTION: function to calculate the geographically corrected RPD (relative phylogenetic distinctiveness) for each species of symbiont based on the hosts used by the symbiont.

#NOTE This function uses the version of the geocorrected metric that excludes intramodular species for the intercommunity comparisons.

## Methods:
# The RPD is defined in Calatayud (2016; they provide the equation):
# First, we calculate the phylogenetic mean pairwise distances (MPD) between the taxon i and all other taxa belonging to its module (i.e., MPD i intracommunity ) minus the mean pairwise distances of taxon i and all other taxa from different modules (MPD i intercommunity ). Then, we computed relative phylogenetic distinctiveness (RPD) as the inverse of the average between all taxa in a module.
# For the geographically corrected version, the source of species for the null models is constrained by geographic restrictions based on a geographic partition.  In this case, the mpd_inter for each species is calculate based on all the species that belong to the any of the geographic modules represented in the target community.

#################
## Created by Alexandre Mestre
## Date of creation: 04/04/2020

#################
## Parameters
# dist_mat: Matrix of pairwise phylogenetic distances between species. The same function can be used to calculate relative functional distinctiveness if we use a distance matrix of interspecific differences in a functional trait.
# int_mat: matrix storing the symbiont-host interactions, with host species in rows and symbiont species in columns
# geo_partition: vector storing the identity of the geographic module for each of the species (same species order than "partition" vector).

#################
## Output
# RPD_geo: vector with the geographically corrected rpd values for the modules (length equal to the number of modules of the parittion) 

phyloSpec.rpd_geo <- function(dist_mat, int_mat, geo_partition)
{
   # First, we order data based on species labels
   dist_mat <- dist_mat [order(rownames(dist_mat)) , order(colnames(dist_mat)) ]
   geo_partition <- geo_partition[ order(names(geo_partition)) ]
   
   Nmod <- length(colnames(int_mat))
   RPD_geo <- rep(NA,Nmod)
   av <- function(x) {sum(x)/(length(x)-1)} # discount the pairwise with the same species (for mpd_intra)

   for (i in 1:Nmod)
      {
         # Select rows with species belonging to the targent module
         selec <- int_mat[,i] > 0
         Nsp <- sum(selec) # number of species belonging to the target module

         # Calculate mpd_intra
         mpd_intra <- apply(dist_mat[selec, selec, drop = FALSE], MARGIN = 1, av)
         
         # Calculate mpd_inter
         mpd_inter <- rep(NA,Nsp)
         iselec <- which(int_mat[,i] > 0)
         sp_intra <- rownames(int_mat)[selec]
         for (j in 1:Nsp)
            {
               # select the hosts belonging to any geo-module represented in the set of hosts used by
               # the target symbiont species
               geoselec <- geo_partition %in% geo_partition[iselec]
               sp_geo <- colnames(dist_mat)[geoselec]
               sp_inter <- sp_geo[!(sp_geo %in% sp_intra)]
               interselec <- colnames(dist_mat) %in% sp_inter
               mpd_inter[j] <- av(dist_mat[iselec[j], interselec])
            }
         # Calculate rpd
         RPD_geo[i] <- -1*sum(mpd_intra - mpd_inter, na.rm = TRUE)/Nsp
      }

   return(RPD_geo)
}
