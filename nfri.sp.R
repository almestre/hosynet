#################
# nfri.sp.R
#################

#################
## DESCRIPTION: function to calculate the geographically corrected NFRI (net functional relatedness index, also called nearest relative index) for each set of hosts used by a symbiont species. It can be used to analyse the diversity of host body sizes used by a symbiont species.

## Methods:
# The NFRI is akin to the NRI index of phylogenetic relatedness.
# In this case, the difference is that we use a functional distance matrix instead of a phylogenetic distance matrix.
# The nfri.sp() function is an adjustment of the nri() that directly uses a distance matrix as input instead of a phylo tree. Otherwise, is equivalent to nri().
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
# com.matrix: MxN incidence matrix where M (i.e. rows) represent species of symbionts, and N (i.e. cols) reperesent species of hosts.
# dist_matrix: matrix of pairwise functional distances ("as.matrix" format).
# geo_partition: vector that stores the identities of the geographic modules for species belonging to the target level (i.e. hosts). If geo_partition = NA, then standard nri is provided instead of a geo-corrected version.
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

nfri.sp <- function(com.matrix, dist_matrix, geo_partition = NA, runs = 999)
{
# First we build the community matrix for the ses.mpd function (where modules represent communities)
# Rows of the community matrix = modules
# Cols of the community matrix = host species
# Nsp <- tapply(bpartg1,bpartg1,length) # number of species (symbionts + hosts) per module

   if (sum(is.na(geo_partition) == 1))
      {  ### Standard NRI
         NRI <- ses.mpd(com.matrix,dist_matrix,null.model="taxa.labels",runs=runs)
      }
   else
      {  ### Geo-corrected NRI

         # containers to store lists
         lcom.mat <- list()
         ldist <- list()

         for (i in 1:nrow(com.matrix))
            { # Build a list of community matrices (with a single row representing the target module)
              # We add a row with zeros to prevent an error in the ses.mpd()
              # ses.mpd() requires at least two communities (i.e. two rows)
              # Then we will keep only the first row
              lcom.mat[[i]] <- com.matrix[ i , ,drop = FALSE]
              lcom.mat[[i]] <- rbind(lcom.mat[[i]], rep(0,ncol(com.matrix)))

              # Build a list of distance matrices with the geographically available species for the target module
              # Select geographic modules represented in the target interaction module
              selec1 <- com.matrix[ i , ] > 0
              geomodID <- geo_partition[selec1]
              # Build a distance matrix containing only the species from the selected geographic modules
              selec2 <- geo_partition %in% geomodID
              ldist[[i]] <- dist_matrix[selec2,selec2]
            }

          f <- function (x,y)
               {
                  # for cases with only one species in the tree
                  if (length(labels(y)) == 0)
                     {
                        data <- rep(NA,8)   
                     }
                  else
                     {
                  data <- ses.mpd(
                           samp = x,
                           dis = y,
                           null.model = "taxa.labels",
                           runs = runs
                           )[1,] # we keep the first row
                     }
                  return(data)
               }

          NRItemp <- mapply (f,
                             lcom.mat,
                             ldist,
                             SIMPLIFY = FALSE
                            )

          NRI <- NRItemp[[1]]
          for (i in 2:nrow(com.matrix))
             {
               NRI <- rbind(NRI,NRItemp[[i]])
             }
      
      }
   return(NRI)
}
