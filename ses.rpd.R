#################
# ses.rpd.R
#################

#################
## DESCRIPTION: function to calculate the standardised effects of RPD (relative phylogenetic distinctiveness) for each module of a modular interaction network (for a given partition).

## Methods:
# The RPD is defined in Calatayud (2016; they provide the equation):
# First, we calculate the phylogenetic mean pairwise distances (MPD) between the taxon i and all other taxa belonging to its module (i.e., MPD i intracommunity ) minus the mean pairwise distances of taxon i and all other taxa from different modules (MPD i intercommunity ). Then, we computed relative phylogenetic distinctiveness (RPD) as the inverse of the average between all taxa in a module.

# The function provides random mean and sd obtained from null models based on randomisation of labels in the distance matrix. We also provide the z-score (z.rpd) and the p-value of the observed rpd.  Moreover, for intercommunity comparisons, it includes all the species that belong to any of the geographic modules represented in the target community (either an interaction module or the set of hosts used by a symbiont)

# The z.rpd is calculated as follows: z.rpd = obs_value - mean(null_values)/sd(null_values)
# The p-value is the proportion of null values below (if obs_value - mean(null_values) < 0) or above (if obs_value - mean(null_values) < 0) the observed value.

#################
## Created by Alexandre Mestre
## Date of creation: 14/03/2020

#################
## Parameters
# dist_mat: Matrix of pairwise phylogenetic distances between species. The same function can be used to calculate relative functional distinctiveness if we use a distance matrix of interspecific differences in a functional trait.
# partition: vector storing the module indentity for each of the species (same order than rows and columns of the distance matrix).
# geo_partition: vector storing the identity of the geographic module for each of the species (same species order than "partition" vector). If geo_partition = NA, then output is based on the standard rpd (i.e. without geo-correction). If a geo_partition is provided, then the output is based on the rpd_geo.
# runs: number of null randomisations (i.e. number of null values that will be created to test for significance of the observed value)

#################
##### Output
##### A data frame of results for each interaction module:
#ntaxa Number of taxa in the interaction module
#rpd.obs Observed mpd in the interaction module
#rpd.rand.mean Mean mpd in null communities
#rpd.rand.sd Standard deviation of mpd in null communities
#rpd.obs.rank Rank of observed mpd vs. null communities
#rpd.obs.z Standardized effect size of rpd vs. null modules (= (rpd.obs - rpd.rand.mean)/ rpd.rand.sd)
#rpd.obs.p P-value (quantile) of observed rpd vs. null modules (= rpd.obs.rank / runs + 1)
#runs Number of randomizations

ses.rpd <- function(dist_mat, partition, geo_partition = NA, runs=999)
{
   if (sum(is.na(geo_partition) == 1))
   {
      rpd.obs <- rpd(dist_mat = dist_mat, partition = partition)
   }
   else
   {
      rpd.obs <- rpd_geo2(dist_mat = dist_mat, partition = partition, geo_partition = geo_partition)
   }

   #function to randomise labels of the distance matrix
   rand.taxa.labels <- function(x)
                          {
                              rand.taxa <- sample(colnames(x))
                              colnames(x) <- rand.taxa
                              rownames(x) <- rand.taxa
                              return(x)
                          }

   if (sum(is.na(geo_partition) == 1))
   {   
      rpd.rand <- t(replicate(runs, rpd( dist_mat = rand.taxa.labels(dist_mat), partition = partition )))
   }
   else
   {
      rpd.rand <- t(replicate(runs, rpd_geo2( dist_mat = rand.taxa.labels(dist_mat), partition = partition, geo_partition = geo_partition )))
   }

   ntaxa <- tapply(partition, partition, length)
   modID <- sort(unique(partition))
   rpd.rand.mean <- apply(X = rpd.rand, MARGIN = 2, FUN = mean, na.rm = TRUE)
   rpd.rand.sd <- apply(X = rpd.rand, MARGIN = 2, FUN = sd, na.rm = TRUE)
   rpd.obs.z <- (rpd.obs - rpd.rand.mean)/rpd.rand.sd
   rpd.obs.rank <- apply(X = rbind(rpd.obs, rpd.rand), MARGIN = 2, FUN = rank)[1, ]
   rpd.obs.rank <- ifelse(is.na(rpd.rand.mean), NA, rpd.obs.rank)

#   return(list(ntaxa, rpd.obs, rpd.rand.mean, 
#            rpd.rand.sd, rpd.obs.rank, rpd.obs.z, rpd.obs.rank/(runs + 
#                1), runs, row.names = modID))

   return(data.frame(
                     ntaxa,
                     rpd.obs,
                     rpd.rand.mean,
                     rpd.rand.sd,
                     rpd.obs.rank,
                     rpd.obs.z,
                     rpd.obs.p = ((runs + 1) - rpd.obs.rank)/(runs + 1), # p-valor basat en obs > null
                     runs = runs,
                     row.names = modID
          ))

}

