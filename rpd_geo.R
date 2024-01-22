#################
# rpd_geo.R
#################

#################
## DESCRIPTION: function to calculate the geographically corrected RPD (relative phylogenetic distinctiveness) for each module of a modular interaction network (for a given partition).

#NOTE This function calculates the version of the geocorrected metric that excludes intramodular species for the intercommunity comparisons. Moreover, for intercommunity comparisons, it includes all the species that belong to any of the geogrpahic modules represented in the target community (either an interaction module or the set of hosts used by a symbiont)

## Methods:
# The RPD is defined in Calatayud (2016; they provide the equation):
# First, we calculate the phylogenetic mean pairwise distances (MPD) between the taxon i and all other taxa belonging to its module (i.e., MPD i intracommunity ) minus the mean pairwise distances of taxon i and all other taxa from different modules (MPD i intercommunity ). Then, we computed relative phylogenetic distinctiveness (RPD) as the inverse of the average between all taxa in a module.
# For the geographically corrected version, the source of species for the null models is constrained by geographic restrictions based on a geographic partition. In this case, the mpd_inter for each species is calculate based on all the species that belong to the any of the geographic modules represented in the target community.

#################
## Created by Alexandre Mestre
## Date of creation: 07/04/2020

#################
## Parameters
# dist_mat: Matrix of pairwise phylogenetic distances between species. The same function can be used to calculate relative functional distinctiveness if we use a distance matrix of interspecific differences in a functional trait.
# partition: vector storing the module indentity for each of the species. The vector should be labeled with the species names, because labels are used to order distance matrix according to the partition vector.
# geo_partition: vector storing the identity of the geographic module for each of the species (same species order than "partition" vector).

#################
## Output
# RPD_geo: vector with the geographically corrected rpd values for the modules (length equal to the number of modules of the parittion) 

rpd_geo <- function(dist_mat, partition, geo_partition)
{
   # First, we order data based on species labels
   dist_mat <- dist_mat [order(rownames(dist_mat)) , order(colnames(dist_mat)) ]
   partition <- partition[ order(names(partition)) ]
   geo_partition <- geo_partition[ order(names(geo_partition)) ]
   
   modID <- sort(unique(partition))
   Nmod <- length(modID)
   RPD_geo <- rep(NA,Nmod)
   av <- function(x) {sum(x)/(length(x)-1)} # discount the pairwise with the same species (for mpd_intra)

   for (i in 1:Nmod)
      {
         # Select rows with species belonging to the targent module
         selec <- partition == modID[i]
         Nsp <- sum(selec) # number of species belonging to the target module

         # Calculate mpd_intra
         mpd_intra <- apply(dist_mat[selec, selec, drop = FALSE], MARGIN = 1, av)
         
         # Calculate mpd_inter
         mpd_inter <- rep(NA,Nsp)
         iselec <- which(partition == modID[i])
         sp_intra <- colnames(dist_mat)[selec]
         for (j in 1:Nsp)
            {
               # select the species belonging to any geo-module that appear in the target int-module
               geoselec <- geo_partition %in% geo_partition[iselec]
               sp_geo <- colnames(dist_mat)[geoselec]
               sp_inter <- sp_geo[!(sp_geo %in% sp_intra)]
               interselec <- colnames(dist_mat) %in% sp_inter
               mpd_inter[j] <- mean(dist_mat[iselec[j], interselec])
            }
         # Calculate rpd
         RPD_geo[i] <- -1*sum(mpd_intra - mpd_inter, na.rm = TRUE)/Nsp
      }

   return(RPD_geo)
}
