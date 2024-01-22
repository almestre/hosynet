#################
# ses.pd_geo.R
#################

#################
## DESCRIPTION: function to calculate the geographically corrected version of the "standardized effect size of phylogenetic diversity" (sesPD).

## Methods:
# The sesPD is calculated using the function ses.pd.fixed() of picante. This function calculates the sesPD for each community (each set of hosts used by a target symbiont species in our case).
# The geographically corrected version of sesPD is based on null models that assume geographic restrictions on availability of taxa. For that, null models are built from species pools composed only of those species belonging to all the geographic modules represented in the set of hosts used by the target symbiont.

# Community matrix as input:
# Rows of the community matrix = communities (symbiont species in our case)
# Cols of the community matrix = species (host species in our case)

# We use the option 'taxa.labels' for the null models. The function creates a specific distance matrix for each symbiont, which is composed of the pool of geographically available hosts to the symbiont.

#################
## Created by Alexandre Mestre
## Date of creation: 16/03/2020

#################
## Parameters
# matrix: MxN incidence matrix where M (i.e. rows) represent species of level one (e.g. hosts), and N (i.e. cols) reperesent species of the level 2 (e.g. symbionts).
# tree: phylogenetic tree of species belonging to the target level
# geopartition: vector that stores the identities of the geographic modules for species belonging to the level 1 (hosts in our case).
# runs: number of randomizations

#################
## A data frame of results for each community
#ntaxa Number of taxa in community
#pd.obs Observed PD in community
#pd.rand.mean Mean PD in null communities
#pd.rand.sd Standard deviation of PD in null communities
#pd.obs.rank Rank of observed PD vs. null communities
#pd.obs.z Standardized effect size of PD vs. null communities (= (pd.obs - pd.rand.mean) / pd.rand.sd)
#pd.obs.p P-value (quantile) of observed PD vs. null communities (= mpd.obs.rank / runs + 1)
#runs Number of randomizations

ses.pd_geo <- function(com.matrix, tree, geopartition, runs = 999)
{
   Nsp <- nrow(com.matrix)

   # containers to store lists
   lcom.mat <- list()
   ltree <- list()
   for (i in 1:Nsp)
      { # Build a list of community matrices (with a single row representing the target module)
        # We add a row with zeros to prevent an error in the ses.pd.fixed()
        # ses.pd.fixed() requires at least two communities (i.e. two rows)
        # Then we will keep only the first row
        lcom.mat[[i]] <- com.matrix[i , ,drop = FALSE]
        lcom.mat[[i]] <- rbind(lcom.mat[[i]], rep(0,ncol(com.matrix)))

        #### Build a list of trees with the geographically available species for the target module
        # Select geographic modules represented in the target interaction module
        selec1 <- com.matrix[ i , ] > 0
        geomodID <- geopartition[selec1]               
        # Select all the species present in the selected geographic modules
        selec2 <- geopartition %in% geomodID
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
               data <- ses.pd.fixed(
                        samp = x,
                        tree = y,
                        null.model = "taxa.labels",
                        runs = runs,
                        include.root = TRUE
                        )[1,] # we keep the first row
                  }
               return(data)
            }

        sesPD_temp <- mapply (f,
                             lcom.mat,
                             ltree,
                             SIMPLIFY = FALSE
                            )

   sesPD <- sesPD_temp[[1]]
   for (i in 2:Nsp)
      {
          sesPD <- rbind(sesPD,sesPD_temp[[i]])
      }

   return(sesPD)
}

