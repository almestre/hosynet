#################
# H.test.R
#################

#################
## DESCRIPTION: function to test for clustering of a categorical trait within communities.

## Methods:
# The function uses the Shannon entropy index (H). For each community, the function calculates the observed H (Hobs) and compare it against H values estimated based on null models (Hnull).

# Hipotheses:
# Null hipothesis -> Hobs = Hnull (no clustering)
# Alternative hipotehsis -> Hobs < Hnull (clustering)

# The p-value is calculated as the proportion of null values that are below the observed value.

# The function provides the random mean and sd obtained from the null models. It also provide the p-value of the observed H.

# The null models are randomisations of the matrix of trait-type frequencies. The function randomizes rows keeping row sums and using weighted probabilities based on either global frequencies or community-level frequencies.

# Two types of null models are implemented:
# 1. Based on trait frequencies of the global species pool (i.e. same weighted probabilities for all the communities)
# 2. Based on trait frequencies of the geographically available species of the target community (i.e. weighed probabilities vary among communities)


#################
## Created by Alexandre Mestre
## Date of creation: 20/03/2020

#################
## Parameters
# matrix: MxN matrix where M (rows) are community entities (e.g. interaction modules), and N (columns) are types of a categorical variable (e.g. habitat types, genera, geographic modules...). Each row stores the frequency of each variable type in species within the target community.
# pool: If pool = "global", weighted probabilities for null models are based on global pool frequencies. If pool = "community", weighted probabilities for null models are based on community-level frequencies.
# weights: frequencies used as weights in the randomisation of rows. If pool = "global", weights is a vector that applies to all the communities (i.e. rows).  If pool = "community", weights is a matrix with the same dimensions than the input matrix. In this case, each community has its own probability weights (e.g. based on geographic availability of species).
# runs: number of randomizations

#################
##### Output
##### A data frame of results for each community:
#ntaxa Number of taxa in the community
#h.obs Observed H in the community
#h.rand.mean Mean H in null communities
#h.rand.sd Standard deviation of H in null communities
#h.obs.rank Rank of observed H vs. null communities
#rpd.obs.z Standardized effect size of rpd vs. null modules (= (rpd.obs - rpd.rand.mean)/ rpd.rand.sd)
#h.obs.p P-value (quantile) of observed H vs. null modules (= h.obs.rank / runs + 1)
#runs Number of randomizations

H.test <- function ( matrix, pool = "global", weights, runs=999)
{
   h.obs <- H.entropy(matrix)
   h.rand <- t(replicate( runs, H.entropy( H.null( matrix, pool = pool, weights = weights ))))

   ntaxa <- apply(matrix, MARGIN = 1, sum)
   h.rand.mean <- apply(X = h.rand, MARGIN = 2, FUN = mean, na.rm = TRUE)
   h.rand.sd <- apply(X = h.rand, MARGIN = 2, FUN = sd, na.rm = TRUE)
   h.obs.z <- (h.obs - h.rand.mean)/h.rand.sd
   h.obs.rank <- apply(X = rbind(h.obs, h.rand), MARGIN = 2, FUN = rank)[1, ]
   h.obs.rank <- ifelse(is.na(h.rand.mean), NA, h.obs.rank)

   return(data.frame(
                     ntaxa,
                     h.obs,
                     h.rand.mean,
                     h.rand.sd,
                     h.obs.rank,
                     h.obs.z,
                     h.obs.p = (h.obs.rank)/(runs + 1), # p-valor basat en obs < null
                     runs = runs,
                     row.names = rownames(matrix)
          ))

}

