#################
# nri.R
#################

#################
## DESCRIPTION: function to calculate the geographically corrected NRI (net relatedness index, also called nearest relative index) for each module of a modular interaction network (for a given partition). 

## Methods:
# The NRI is calculated using the function ses.mpd() of picante. This functions calculates the standardized effect size of MPD for each community (each module in our case).
# Then, to obtain the NRI, we multiply the sesMPD by -1.
# The geographically corrected version of NRI is based on null models that assume geographic restrictions on availability of taxa. For that, null models are built from species pools composed only of those species belonging to all the geographic modules represented in the target interaction module.
# If we do not provide a geographic partition, then the standard nri is calculated.

# For the geo-corrected version, we use the option 'taxa.labels' for the null models. In that case, the function creates a specific distance matrix for each interaction module, which is composed of the pool of geographically available species. For the standard version of nri, we use the option "taxa.labels"

#################
## Created by Alexandre Mestre
## Date of creation: 16/03/2020

#################
## Parameters
# matrix: MxN incidence matrix where M (i.e. rows) represent species of level one (e.g. hosts), and N (i.e. cols) reperesent species of the level 2 (e.g. symbionts).
# tree: phylogenetic tree of species belonging to the target level
# partition: vector that stores the identities of the interaction modules for the species belonging to the target level (either level 1 or 2)
# geo_partition: vector that stores the identities of the geographic modules for species belonging to the target level (either level 1 or 2). If geo_partition = NA, then standard nri is provided instead of a geo-corrected version.
# level: If level = 1, we will obtain the nri for species of level 1 (i.e rows of the incidence matrix); if level = 2, we obtain nri for species of level 2 (columns the incidence matrix). 
# runs: number of randomizations

#################
##### Output (same output structure than ses.mpd()
##### A data frame of results for each interaction module:
#ntaxa Number of taxa in community
#mpd.obs Observed mpd in community
#mpd.rand.mean Mean mpd in null communities
#mpd.rand.sd Standard deviation of mpd in null communities
#mpd.obs.rank Rank of observed mpd vs. null communities
#mpd.obs.z Standardized effect size of mpd vs. null communities (= (mpd.obs - mpd.rand.mean)/ mpd.rand.sd, equivalent to -NRI)
#mpd.obs.p P-value (quantile) of observed mpd vs. null communities (= mpd.obs.rank / runs + 1)
#runs Number of randomizations

nri <- function(matrix, tree, partition, geo_partition = NA, level = 1, runs = 999)
{
# First we build the community matrix for the ses.mpd function (where modules represent communities)
# Rows of the community matrix = modules
# Cols of the community matrix = host species
# Nsp <- tapply(bpartg1,bpartg1,length) # number of species (symbionts + hosts) per module

   if (level == 2)
      {matrix = t(matrix)}

   Nsp <- nrow(matrix)
   partid <- sort(unique(partition))
   Npart <- length(partid)

   com.matrix <- matrix(rep(0,Npart*Nsp),c(Npart,Nsp))
   colnames(com.matrix) <- rownames(matrix)
   rownames(com.matrix) <- partid

   for (i in 1:Nsp)
      {
         selec <- partid == partition[i]
         com.matrix [selec,i] = 1
      }

   if (sum(is.na(geo_partition) == 1))
      {  ### Standard NRI
         NRI <- ses.mpd(com.matrix,cophenetic(tree),null.model="taxa.labels",runs=runs)
      }
   else
      {  ### Geo-corrected NRI

         # containers to store lists
         lcom.mat <- list()
         ltree <- list()

         for (i in 1:Npart)
            { # Build a list of community matrices (with a single row representing the target module)
              # We add a row with zeros to prevent an error in the ses.mpd()
              # ses.mpd() requires at least two communities (i.e. two rows)
              # Then we will keep only the first row
              lcom.mat[[i]] <- com.matrix[ i , ,drop = FALSE]
              lcom.mat[[i]] <- rbind(lcom.mat[[i]], rep(0,ncol(com.matrix)))

              # Build a list of trees with the geographically available species for the target module
              # Select geographic modules represented in the target interaction module
              selec1 <- com.matrix[ i , ] > 0
              geomodID <- geo_partition[selec1]               
              # Select all the species present in the selected geographic modules
              selec2 <- geo_partition %in% geomodID 
              sp_pool <- colnames(com.matrix)[selec2]
              # Build a subtree containing the selected species
              selec3 <- tree$tip.label %in% sp_pool
              ltree[[i]] <- drop.tip (tree, tree$tip.label[!selec3])

            }

          f <- function (x,y)
               {
                  # for cases with only one species in the tree
                  if (length(y$tip.label) == 1)
                     {
                        data <- rep(NA,8)   
                     }
                  else
                     {
                  data <- ses.mpd(
                           samp = x,
                           dis = cophenetic(y),
                           null.model = "taxa.labels",
                           runs = runs
                           )[1,] # we keep the first row
                     }
                  return(data)
               }

          NRItemp <- mapply (f,
                             lcom.mat,
                             ltree,
                             SIMPLIFY = FALSE
                            )

          NRI <- NRItemp[[1]]
          for (i in 2:Npart)
             {
               NRI <- rbind(NRI,NRItemp[[i]])
             }
      
      }
   return(NRI)
}
