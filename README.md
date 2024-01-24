# hosynet

Matlab, Python and R tools for the analysis of symbiont-host interaction networks

## Author
Alexandre Mestre

## Aims
This repository was developed as part of a research line focused on disentangling the drivers of freshwater host-symbiont interaction communities, using American crayfish-ostracod associations. The results of this research are shown in an ongoing SCI paper:

Mestre A, Mesquita-Joanes F, and Peres-Neto PR (2014). Habitat specialization of crayfish drives specialization in host usage by their ostracod symbionts (In prep.)

## 1_create_null_models.m

### Description
Matlab script that creates the null models for the analysis of the modularity of a bipartite interaction network (e.g., characterising symbiont-host associations).

### Methods
Number of null models = 999
The null models are randomised versions of an MxN adjacency matrix for an undirected bipartite network
Type of null model = fixed null model
Description of the fixed null model:
In this model, each random network has the same number of connections per species as the real network. Random networks (pmatrices) are constructed using an independent swap algorithm in which the original matrix is reshuffled by repeatedly swapping 2 # 2 submatrices that preserve the row and column totals. The null model implemented in our study maintained the bipartite structure of the networks by allowing connections between hosts and symbionts and not allowing connections among hosts or among symbionts.
Matlab function swap.m written by B. Semmens
999 random matrices for each real network


## 2_create_network_partitions.m

### Description
Matlab script that creates a set of network partitions from an MxN adjacency matrix (called A), covering a representative range of gamma values.
The MxN adjacency matrix of symbiont-host associations

## 3_estimate_best_gamma.py

### Description
Python script to estimate the best gamma parameter for the modularity of a bipartite network. 

The estimation is based on a set of network partitions covering a representative range of gamma values, that were previously obtained using the Matlab script '2_create_network_partitions.m'.

The script requires the previous execution of the following Matlab scripts (in the given order):
1. '1_create_null_models.m' (required for 'modularity_test.m')
2. '2_create_network_partitions.m'
3. 'modularity_test.m'

### Methods
It uses the CHAMP package to create an ensemble of partitions that can be used to identify the partition with the largest domain of optimisation based on a graphical approach described in:

Weir et al. (2017). Post-Processing Partitions to Identify Domains of Modularity Optimization. Algorithms, 10:93.

The plot to identify the best partition can be obtained with the following script that processes the ensemble created in this script: 'plot_bestg.py'

## modularity_test.m

### Description
Matlab script that calculates an optimised partition for a bipartite network using the 'iterated_genlouvain' algorithm, and tests the significance of its modularity based on a set of null models (see methods section below).
The null models are obtained from the following Matlab script: 1_nullModels.m

We have to provide a gamma value. The best gamma value can be graphically estimated by applying the CHAMP package, with the following scripts:
1. '1_create_null_models.m' (required for 'modularity_test.m')
2. '2_create_network_partitions.m'
3. '3_estimate_best_gamma.py'
4. 'plot_bestg.py'

### Methods (modularity test)
Number of null models = 999
The null models are randomised versions of an MxN adjacency matrix for an undirected bipartite network
Type of null model = fixed null model
Description of the fixed null model:
In this model, each random network has the same number of connections per species as the real network. Random networks (pmatrices) are constructed using an independent swap algorithm in which the original matrix is reshuffled by repeatedly swapping 2 # 2 submatrices that preserve the row and column totals. The null model implemented in our study maintained the bipartite structure of the networks by allowing connections between hosts and symbionts and not allowing connections among hosts or among symbionts.
For each of the 999 random matrices, we calculate the modularity (Q)
Then, we compare the distribution of Q from null models with the Q of our optimised partition (i.e. the observed network).
The statistical significance of the modularity of the observed network (p-value) is the fraction of the 999 random matrices with a modularity value equal to or higher than the observed one.

## plot_bestg.py

### Description
Python script to obtain the domain-of-optimisation plot that can be used to estimate the best gamma parameter for the modularity of a bipartite network.

The script requires the previous execution of the following Matlab and Python scripts (in the listed order):
1. '1_create_null_models.m' (required for 'modularity_test.m')
2. '2_create_network_partitions.m'
3. 'modularity_test.m'
4. '3_estimate_best_gamma.py'

The interpretation of the resulting plot is explained in:
# Weir et al. (2017). Post-Processing Partitions to Identify Domains of Modularity Optimization. Algorithms, 10:93.

## plot_AMI.py

### Description
Python script to obtain a heatmap plot representing the pairwise adjusted mutual information (AMI) between the optimal partitions identified by the Champ algorithm. 

The script requires the previous execution of the following Matlab and Python scripts (in the listed order):
1. '1_create_null_models.m' (required for 'modularity_test.m')
2. '2_create_network_partitions.m'
3. 'modularity_test.m'
4. '3_estimate_best_gamma.py'

The interpretation of the resulting plot is explained in:
Weir et al. (2017). Post-Processing Partitions to Identify Domains of Modularity Optimization. Algorithms, 10:93.

## plot_graph.py

### Description
Python script to obtain a graphical representation of a bipartite network (i.e. a bipartite graph) that was created with the Matlab script 'modularity_test.m'.

## H.entropy.R

### Description
R function to calculate the H index of entropy of a categorical trait in a set of communities.

### Methods
Shannon entropy (H):
H = -1*sum(pi*ln(pi)) where pi is the proportion of occurrences of the type "i" of the categorical trait in the set of species composing the target community.

### Parameters
matrix: MxN matrix where M (rows) are community entities (e.g. interaction modules), and N (columns) are types of a categorical variable (e.g. habitat types, genera, geographic modules...). Each row stores the frequency of each variable type in species within the target community.

### Output
A vector with the H indices

## H.null.R

### Description
R function to randomize a community matrix of species-level categorical trait frequencies (e.g. habitat, geographic module or genus).

### Methods
The function randomizes rows keeping row sums and using weighted probabilities based on either global frequencies or community-level frequencies.

### Parameters

matrix: MxN matrix where M (rows) are community entities (e.g. interaction modules), and N (columns) are types of a categorical trait (e.g. habitat types, genera, geographic modules...). Each row stores the frequency of each trait variant across species within the target community.

pool: If pool = "global", weighted probabilities are based on global pool frequencies. If pool = "community", weighted probabilities are based on community-level frequencies.

weights: frequencies used as weights in the randomisation of rows. If pool = "global", weights is a vector that applies to all the communities (i.e. rows).  If pool = "community", weights is a matrix with the same dimensions as the input matrix. In this case, each community has its own probability weights (e.g. based on geographic availability of species).

### Output
The randomized matrix

## H.test.R

### Description
R function to test for clustering of a categorical trait within communities.

### Methods:
The function uses the Shannon entropy index (H). For each community, the function calculates the observed H (Hobs) and compare it against H values estimated based on null models (Hnull).

Hypotheses:
Null hypothesis -> Hobs = Hnull (no clustering)
Alternative hypotehsis -> Hobs < Hnull (clustering)

The p-value is calculated as the proportion of null values that are below the observed value.

The function provides the random mean and SD obtained from the null models. It also provides the p-value of the observed H.

The null models are randomisations of the matrix of trait-type frequencies. The function randomizes rows keeping row sums and using weighted probabilities based on either global frequencies or community-level frequencies.

Two types of null models are implemented:
1. Based on trait frequencies of the global species pool (i.e. same weighted probabilities for all the communities)
2. Based on trait frequencies of the geographically available species of the target community (i.e. weighed probabilities vary among communities)

### Parameters

matrix: MxN matrix where M (rows) are community entities (e.g. interaction modules), and N (columns) are types of a categorical variable (e.g. habitat types, genera, geographic modules...). Each row stores the frequency of each variable type in species within the target community.

pool: If pool = "global", weighted probabilities for null models are based on global pool frequencies. If pool = "community", weighted probabilities for null models are based on community-level frequencies.

weights: frequencies used as weights in the randomisation of rows. If pool = "global", weights is a vector that applies to all the communities (i.e. rows).  If pool = "community", weights is a matrix with the same dimensions as the input matrix. In this case, each community has its own probability weights (e.g. based on geographic availability of species).
runs: number of randomizations

### Output
A data frame of results for each community:
ntaxa Number of taxa in the community
h.obs Observed H in the community
h.rand.mean Mean H in null communities
h.rand.sd Standard deviation of H in null communities
h.obs.rank Rank of observed H vs. null communities
rpd.obs.z Standardized effect size of rpd vs. null modules (= (rpd.obs - rpd.rand.mean)/ rpd.rand.sd)
h.obs.p P-value (quantile) of observed H vs. null modules (= h.obs.rank / runs + 1)
runs Number of randomizations

## cast_to_matrix

### Description
R function to build a matrix of pairwise associations of two "character" variables stored in a data frame. Examples of applications:
1. Building a site x species matrix (presence/absences)
2. Building a co-occurrence matrix for two groups of species such as symbionts and hosts (number of co-occurrences per pairwise association)

### Parameters
dataframe: data frame that stores the variables to be cast
variable1: variable that will conform the rows of the output matrix
variable2: variable that will conform the columns of the output matrix
pa: if pa = TRUE, the output matrix is a presences/absences matrix (0 and 1 data). If pa=FALSE (default option), it informs about the number of co-occurrences (0 and positive integers).

## czscores

### Description
R function to calculate the intermodular connectivity (c) and intramodular degree(z) for each node of a bipartite interaction network given a certain partition.

We apply the same formulas as described in bipartite package.
c = 1 - sum( (k.it/k.i)^2) # among-module connectivity
z = (k.is - ks.av) / SD.ks # within-module degree
k.is = number of links of i to other species in its own module s; ks.av = average k.is of all species in module s; SD.ks = standard deviation of k.is of all species in module s; k.it = number of links of species i to module t; k.i = degree of species i Note that for any species alone (in its level) in a module the z-value will be NaN, since then SD.ks is 0. This is a limitation of the way the z-value is defined (in multiples of degree/strength standard deviations).

Following Krasnov2012, we assign z = 0 to cases wherein sd(z) = 0 (reasoning well explained in Krasnov2012).

### Parameters
matrix: MxN incidence matrix where M (i.e. rows) represent species of level one (e.g. hosts), and N (i.e. cols) reperesent species of the level 2 (e.g. symbionts).
partition: vector of size M+N that stores the labels of module indentities of nodes, with species of level one first, followed by species of level two.

### Output
A list with two elements:
c is a vector with the c values
z is a vector with the z values

## getGRS_DynamicAlpha

### Description
R function to estimate geographic range size (GRS) of a species from occurrence data.

### Methods
Extent of occurrence (EOO) based on alpha-convex hull calculated with dynamic alpha (see getDynamicAlphaHull function from packagage rangeBuilder)
Area of occupancy (AOO) calculated with the conR package

### Input
"data" is a data frame with the occurrence data. It should have three columns (it is mandatory to respect column positions):
The first column is the coordx
The second column is the coordy
The third column is the species names
"sp"
 Character value indicating the name of the target species for GRS calculations
"Cell_size_AOO"
 Numeric value indicating the grid size (km) used for AOO estimations. By default, equal to 2

### Output
Object of class SpatialPolygons storing the alpha-convex hull
EOO in km² 
AOO in km²

## int.to.com.matrix

### Description
R function to convert an incidence matrix for bipartite networks (i.e. sp_level1 x sp_level2) into a community matrix (modules x sp) based on a partition (i.e. a vector storing module indentities of the species).

### Parameters
matrix: incidence matrix with MxN dimensions where M (rows) are species of level 1 (e.g. hosts), and M (cols) are species of level 2 (e.g. symbionts)
partition: vector that stores the identities of the interaction modules for the species belonging to the target level (either level 1 or 2). Elements should keep the same order than species labels in the matrix.
level: If level = 1, we will obtain the nri for species of level 1 (i.e rows of the incidence matrix); if level = 2, we obtain nri for species of level 2 (columns the incidence matrix).

### Output
Community matrix with MxN dimensions where M (rows) are community entities (i.e. interaction modules), and N (cols) are species of the target level.

## nfri.R

### Description
R function to calculate the geographically corrected NFRI (net functional relatedness index, also called nearest relative index) for each module of a modular interaction network (for a given partition).

### Methods
The NFRI is akin to the NRI index of phylogenetic relatedness.
In this case, the difference is that we use a functional distance matrix instead of a phylogenetic distance matrix.
The nfri() function is an adjustment of the nri() that directly uses a distance matrix as input instead of a phylo tree. Otherwise, is equivalent to nri().
The NRI is calculated using the function ses.mpd() of picante. This functions calculates the standardized effect size of MPD for each community (each module in our case).
Then, to obtain the NRI, we multiply the sesMPD by -1.
The geographically corrected version of NRI is based on null models that assume geographic restrictions on availability of taxa. For that, null models are built from species pools composed only of those species belonging to all the geographic modules represented in the target interaction module.
If we do not provide a geographic partition, then the standard nri is calculated.
For the geo-corrected version, we use the option 'taxa.labels' for the null models. In that case, the function creates a specific distance matrix for each interaction module, which is composed of the pool of geographically available species. For the standard version of nri, we use the option "taxa.labels".

### Parameters
matrix: MxN incidence matrix where M (i.e. rows) represent species of level one (e.g. hosts), and N (i.e. cols) reperesent species of the level 2 (e.g. symbionts).
# dist_matrix: matrix of pairwise functional distances ("as.matrix" format).
partition: vector that stores the identities of the interaction modules for the species belonging to the target level (either level 1 or 2)
geo_partition: vector that stores the identities of the geographic modules for species belonging to the target level (either level 1 or 2). If geo_partition = NA, then standard nri is provided instead of a geo-corrected version.
level: If level = 1, we will obtain the nri for species of level 1 (i.e rows of the incidence matrix); if level = 2, we obtain nri for species of level 2 (columns the incidence matrix). 
runs: number of randomizations

### Output
A data frame of results for each interaction module:
ntaxa Number of taxa in community
mpd.obs Observed mpd in community
mpd.rand.mean Mean mpd in null communities
mpd.rand.sd Standard deviation of mpd in null communities
mpd.obs.rank Rank of observed mpd vs. null communities
mpd.obs.z Standardized effect size of mpd vs. null communities (= (mpd.obs - mpd.rand.mean)/ mpd.rand.sd, equivalent to -NRI)
mpd.obs.p P-value (quantile) of observed mpd vs. null communities (= mpd.obs.rank / runs + 1)
runs Number of randomizations

## nfri.sp.R

### Description
R function to calculate the geographically corrected NFRI (net functional relatedness index, also called nearest relative index) for each set of hosts used by a symbiont species. It can be used to analyse the diversity of host body sizes used by a symbiont species.

### Methods
The NFRI is akin to the NRI index of phylogenetic relatedness.
In this case, the difference is that we use a functional distance matrix instead of a phylogenetic distance matrix.
The nfri.sp() function is an adjustment of the nri() that directly uses a distance matrix as input instead of a phylo tree. Otherwise, is equivalent to nri().
The NRI is calculated using the function ses.mpd() of picante. This functions calculates the standardized effect size of MPD for each community (each module in our case).
Then, to obtain the NRI, we multiply the sesMPD by -1.
The geographically corrected version of NRI is based on null models that assume geographic restrictions on availability of taxa. For that, null models are built from species pools composed only of those species belonging to all the geographic modules represented in the target interaction module.
If we do not provide a geographic partition, then the standard nri is calculated.
For the geo-corrected version, we use the option 'taxa.labels' for the null models. In that case, the function creates a specific distance matrix for each interaction module, which is composed of the pool of geographically available species. For the standard version of nri, we use the option "taxa.labels"

### Parameters
com.matrix: MxN incidence matrix where M (i.e. rows) represent species of symbionts, and N (i.e. cols) reperesent species of hosts.
dist_matrix: matrix of pairwise functional distances ("as.matrix" format).
geo_partition: vector that stores the identities of the geographic modules for species belonging to the target level (i.e. hosts). If geo_partition = NA, then standard nri is provided instead of a geo-corrected version.
runs: number of randomizations

### Output
A data frame of results for each interaction module:
ntaxa Number of taxa in community
mpd.obs Observed mpd in community
mpd.rand.mean Mean mpd in null communities
mpd.rand.sd Standard deviation of mpd in null communities
mpd.obs.rank Rank of observed mpd vs. null communities
mpd.obs.z Standardized effect size of mpd vs. null communities (= (mpd.obs - mpd.rand.mean)/ mpd.rand.sd, equivalent to -NRI)
mpd.obs.p P-value (quantile) of observed mpd vs. null communities (= mpd.obs.rank / runs + 1)
runs Number of randomizations

## nri.R

### Description
R function to calculate the geographically corrected NRI (net relatedness index, also called nearest relative index) for each module of a modular interaction network (for a given partition). 

### Methods
The NRI is calculated using the function ses.mpd() of picante. This functions calculates the standardized effect size of MPD for each community (each module in our case).
Then, to obtain the NRI, we multiply the sesMPD by -1.
The geographically corrected version of NRI is based on null models that assume geographic restrictions on availability of taxa. For that, null models are built from species pools composed only of those species belonging to all the geographic modules represented in the target interaction module.
If we do not provide a geographic partition, then the standard nri is calculated.

For the geo-corrected version, we use the option 'taxa.labels' for the null models. In that case, the function creates a specific distance matrix for each interaction module, which is composed of the pool of geographically available species. For the standard version of nri, we use the option "taxa.labels"

### Parameters
matrix: MxN incidence matrix where M (i.e. rows) represent species of level one (e.g. hosts), and N (i.e. cols) reperesent species of the level 2 (e.g. symbionts).
tree: phylogenetic tree of species belonging to the target level
partition: vector that stores the identities of the interaction modules for the species belonging to the target level (either level 1 or 2)
geo_partition: vector that stores the identities of the geographic modules for species belonging to the target level (either level 1 or 2). If geo_partition = NA, then standard nri is provided instead of a geo-corrected version.
level: If level = 1, we will obtain the nri for species of level 1 (i.e rows of the incidence matrix); if level = 2, we obtain nri for species of level 2 (columns the incidence matrix). 
runs: number of randomizations

### Output
A data frame of results for each interaction module:
ntaxa Number of taxa in community
mpd.obs Observed mpd in community
mpd.rand.mean Mean mpd in null communities
mpd.rand.sd Standard deviation of mpd in null communities
mpd.obs.rank Rank of observed mpd vs. null communities
mpd.obs.z Standardized effect size of mpd vs. null communities (= (mpd.obs - mpd.rand.mean)/ mpd.rand.sd, equivalent to -NRI)
mpd.obs.p P-value (quantile) of observed mpd vs. null communities (= mpd.obs.rank / runs + 1)
runs Number of randomizations

#################
## null.weights.geo.R

### Description 
R function to calculate a matrix of weighted probabilities of trait frequencies based on geographic availability of species in the target community.

### Methods
The function uses a community matrix and a vector of geographic modules. For each community (i.e. row of the community matrix), it obtains the set of species that belong to the geographic modules represented in the target community. Then, it calculates the frequencies of trait variants in the selected species and store them as a row of the output matrix.

### Parameters
com.matrix: MxN matrix where M (rows) are community entities (e.g. interaction modules), and N (columns) are species.
trait: vector with trait variants of species composing the global pool (same order than columns of com.matrix).
geo_partition: vector with the geographic module identities off all the species composing the global pool (same order than columns of com.matrix).

#### Output
Matrix with the weighted probabilities, where rows are community entities and columns are trait variants.

## phyloSpec.rpd.R

### Description
R function to calculate the RPD (relative phylogenetic distinctiveness) for each species of symbiont based on the hosts used by the symbiont.

### Methods
The RPD is defined in Calatayud (2016; they provide the equation):
First, we calculate the phylogenetic mean pairwise distances (MPD) between the taxon i and all other taxa belonging to its module (i.e., MPD i intracommunity ) minus the mean pairwise distances of taxon i and all other taxa from different modules (MPD i intercommunity ). Then, we computed relative phylogenetic distinctiveness (RPD) as the inverse of the average between all taxa in a module. Here, modules represent the set of hosts used by a symbiont.

### Parameters
dist_mat: Matrix of pairwise phylogenetic distances between species. The same function can be used to calculate relative functional distinctiveness if we use a distance matrix of interspecific differences in a functional trait.
int_mat: matrix storing the symbiont-host interactions, with host species in rows and symbiont species in columns

### Output
RPD: vector with the rpd values for the species (length equal to the number of symbiont species present in the interaction matrix)

## phyloSpec.rpd_geo.R

### Description
R function to calculate the geographically corrected RPD (relative phylogenetic distinctiveness) for each species of symbiont based on the hosts used by the symbiont.
This function uses the version of the geocorrected metric that excludes intramodular species for the intercommunity comparisons.

### Methods:
The RPD is defined in Calatayud (2016; they provide the equation):
First, we calculate the phylogenetic mean pairwise distances (MPD) between the taxon i and all other taxa belonging to its module (i.e., MPD i intracommunity ) minus the mean pairwise distances of taxon i and all other taxa from different modules (MPD i intercommunity ). Then, we computed relative phylogenetic distinctiveness (RPD) as the inverse of the average between all taxa in a module.
For the geographically corrected version, the source of species for the null models is constrained by geographic restrictions based on a geographic partition.  In this case, the mpd_inter for each species is calculate based on all the species that belong to the any of the geographic modules represented in the target community.

### Parameters
dist_mat: Matrix of pairwise phylogenetic distances between species. The same function can be used to calculate relative functional distinctiveness if we use a distance matrix of interspecific differences in a functional trait.
int_mat: matrix storing the symbiont-host interactions, with host species in rows and symbiont species in columns
geo_partition: vector storing the identity of the geographic module for each of the species (same species order than "partition" vector).

### Output
RPD_geo: vector with the geographically corrected rpd values for the modules (length equal to the number of modules of the parittion)

## rpd.R

### Description
R function to calculate the RPD (relative phylogenetic distinctiveness) for each module of a modular interaction network (for a given partition).

### Methods
The RPD is defined in Calatayud (2016; they provide the equation):
First, we calculate the phylogenetic mean pairwise distances (MPD) between the taxon i and all other taxa belonging to its module (i.e., MPD i intracommunity ) minus the mean pairwise distances of taxon i and all other taxa from different modules (MPD i intercommunity ). Then, we computed relative phylogenetic distinctiveness (RPD) as the inverse of the average between all taxa in a module.

### Parameters
dist_mat: Matrix of pairwise phylogenetic distances between species. The same function can be used to calculate relative functional distinctiveness if we use a distance matrix of interspecific differences in a functional trait.
partition: vector storing the module indentity for each of the species. The vector should be labeled with the species names, because labels are used to order distance matrix according to the partition vector.

### Output
RPD: vector with the rpd values for the modules (length equal to the number of modules of the parittion)

## rpd_geo.R

### Description
R function to calculate the geographically corrected RPD (relative phylogenetic distinctiveness) for each module of a modular interaction network (for a given partition).
This function calculates the version of the geocorrected metric that excludes intramodular species for the intercommunity comparisons. Moreover, for intercommunity comparisons, it includes all the species that belong to any of the geogrpahic modules represented in the target community (either an interaction module or the set of hosts used by a symbiont)

### Methods:
The RPD is defined in Calatayud (2016; they provide the equation):
First, we calculate the phylogenetic mean pairwise distances (MPD) between the taxon i and all other taxa belonging to its module (i.e., MPD i intracommunity ) minus the mean pairwise distances of taxon i and all other taxa from different modules (MPD i intercommunity ). Then, we computed relative phylogenetic distinctiveness (RPD) as the inverse of the average between all taxa in a module.
For the geographically corrected version, the source of species for the null models is constrained by geographic restrictions based on a geographic partition. In this case, the mpd_inter for each species is calculate based on all the species that belong to the any of the geographic modules represented in the target community.

### Parameters
dist_mat: Matrix of pairwise phylogenetic distances between species. The same function can be used to calculate relative functional distinctiveness if we use a distance matrix of interspecific differences in a functional trait.
partition: vector storing the module indentity for each of the species. The vector should be labeled with the species names, because labels are used to order distance matrix according to the partition vector.
geo_partition: vector storing the identity of the geographic module for each of the species (same species order than "partition" vector).

### Output
RPD_geo: vector with the geographically corrected rpd values for the modules (length equal to the number of modules of the parittion)

## ses.pd.fixed

### Description
Corrected function of ses.pd function of package picante.
Source:
https://www.biorxiv.org/content/10.1101/579300v1.supplementary-material

## ses.pd_geo.R

### Description
R function to calculate the geographically corrected version of the "standardized effect size of phylogenetic diversity" (sesPD).

### Methods:
The sesPD is calculated using the function ses.pd.fixed() of picante. This function calculates the sesPD for each community (each set of hosts used by a target symbiont species in our case).
The geographically corrected version of sesPD is based on null models that assume geographic restrictions on availability of taxa. For that, null models are built from species pools composed only of those species belonging to all the geographic modules represented in the set of hosts used by the target symbiont.

Community matrix as input:
Rows of the community matrix = communities (symbiont species in our case)
Cols of the community matrix = species (host species in our case)

We use the option 'taxa.labels' for the null models. The function creates a specific distance matrix for each symbiont, which is composed of the pool of geographically available hosts to the symbiont.

### Parameters
matrix: MxN incidence matrix where M (i.e. rows) represent species of level one (e.g. hosts), and N (i.e. cols) reperesent species of the level 2 (e.g. symbionts).
tree: phylogenetic tree of species belonging to the target level
geopartition: vector that stores the identities of the geographic modules for species belonging to the level 1 (hosts in our case).
runs: number of randomizations

### Output
A data frame of results for each community with the following content:
ntaxa Number of taxa in community
pd.obs Observed PD in community
pd.rand.mean Mean PD in null communities
pd.rand.sd Standard deviation of PD in null communities
pd.obs.rank Rank of observed PD vs. null communities
pd.obs.z Standardized effect size of PD vs. null communities (= (pd.obs - pd.rand.mean) / pd.rand.sd)
pd.obs.p P-value (quantile) of observed PD vs. null communities (= mpd.obs.rank / runs + 1)
#runs Number of randomizations

## ses.phyloSpec.rpd.R

### Description
R function to calculate the standardised effects of RPD (relative phylogenetic distinctiveness) for each symbiont species, based on the hosts used by the symbiont.

### Methods:
The RPD is defined in Calatayud (2016; they provide the equation):
First, we calculate the phylogenetic mean pairwise distances (MPD) between the taxon i and all other taxa belonging to its module (i.e., MPD i intracommunity ) minus the mean pairwise distances of taxon i and all other taxa from different modules (MPD i intercommunity ). Then, we computed relative phylogenetic distinctiveness (RPD) as the inverse of the average between all taxa in a module. Here the module represents the set of hosts used by a symbiont, and the RPD is a metric of symbiont specialization.
For the geographically corrected version, the source of species for the null models is constrained by geographic restrictions based on a geographic partition. In this case, the mpd_inter for each species is calculate based on all the species that belong to the any of the geographic modules represented in the target community.

The function provides random mean and sd obtained from null models based on randomisation of labels in the distance matrix. We also provide the z-score (z.rpd) and the p-value of the observed rpd.

The z.rpd is calculated as follows: z.rpd = obs_value - mean(null_values)/sd(null_values)
The p-value is the proportion of null values below (if obs_value - mean(null_values) < 0) or above (if obs_value - mean(null_values) < 0) the observed value.

### Parameters
dist_mat: Matrix of pairwise phylogenetic distances between species. The same function can be used to calculate relative functional distinctiveness if we use a distance matrix of interspecific differences in a functional trait.
int_mat: matrix storing the symbiont-host interactions, with host species in rows and symbiont species in columns
geo_partition: vector storing the identity of the geographic module for each of the species (same species order than "partition" vector). If geo_partition = NA, then output is based on the standard rpd (i.e. without geo-correction). If a geo_partition is provided, then the output is based on the rpd_geo2.
runs: number of null randomisations (i.e. number of null values that will be created to test for significance of the observed value)

### Output
A data frame of results for each interaction module. With the following content:
ntaxa Number of taxa in the interaction module
rpd.obs Observed mpd in the interaction module
rpd.rand.mean Mean mpd in null communities
rpd.rand.sd Standard deviation of mpd in null communities
rpd.obs.rank Rank of observed mpd vs. null communities
rpd.obs.z Standardized effect size of rpd vs. null modules (= (rpd.obs - rpd.rand.mean)/ rpd.rand.sd)
rpd.obs.p P-value (quantile) of observed rpd vs. null modules (= rpd.obs.rank / runs + 1)
runs Number of randomizations

## ses.rpd.R

### Description
R function to calculate the standardised effects of RPD (relative phylogenetic distinctiveness) for each module of a modular interaction network (for a given partition).

### Methods
The RPD is defined in Calatayud (2016; they provide the equation):
First, we calculate the phylogenetic mean pairwise distances (MPD) between the taxon i and all other taxa belonging to its module (i.e., MPD i intracommunity ) minus the mean pairwise distances of taxon i and all other taxa from different modules (MPD i intercommunity ). Then, we computed relative phylogenetic distinctiveness (RPD) as the inverse of the average between all taxa in a module.
The function provides random mean and sd obtained from null models based on randomisation of labels in the distance matrix. We also provide the z-score (z.rpd) and the p-value of the observed rpd.  Moreover, for intercommunity comparisons, it includes all the species that belong to any of the geographic modules represented in the target community (either an interaction module or the set of hosts used by a symbiont)
The z.rpd is calculated as follows: z.rpd = obs_value - mean(null_values)/sd(null_values)
The p-value is the proportion of null values below (if obs_value - mean(null_values) < 0) or above (if obs_value - mean(null_values) < 0) the observed value.

### Parameters
dist_mat: Matrix of pairwise phylogenetic distances between species. The same function can be used to calculate relative functional distinctiveness if we use a distance matrix of interspecific differences in a functional trait.
partition: vector storing the module indentity for each of the species (same order than rows and columns of the distance matrix).
geo_partition: vector storing the identity of the geographic module for each of the species (same species order than "partition" vector). If geo_partition = NA, then output is based on the standard rpd (i.e. without geo-correction). If a geo_partition is provided, then the output is based on the rpd_geo.
runs: number of null randomisations (i.e. number of null values that will be created to test for significance of the observed value)

### Output
A data frame of results for each interaction module, with the following content:
ntaxa Number of taxa in the interaction module
rpd.obs Observed mpd in the interaction module
rpd.rand.mean Mean mpd in null communities
rpd.rand.sd Standard deviation of mpd in null communities
rpd.obs.rank Rank of observed mpd vs. null communities
rpd.obs.z Standardized effect size of rpd vs. null modules (= (rpd.obs - rpd.rand.mean)/ rpd.rand.sd)
rpd.obs.p P-value (quantile) of observed rpd vs. null modules (= rpd.obs.rank / runs + 1)
runs Number of randomizations

## sort_dist_mat

### Description
R function to sort a distance matrix (or a regular matrix)
Source: https://rdrr.io/cran/dendextend/src/R/cor_cophenetic.R

## sp.to.tfreq.matrix

### Description
R function to convert a community matrix with species occurrences into a community matrix with frequencies of variants a categorigcal trait e.g. (geographic module, genus...).

### Parameters
com.matrix: community matrix with MxN dimensions where M (rows) are community entities (e.g. interaction modules), and N (cols) are species.
trait.values: vector storing trait variants of all the species of the global pool (same order than columns of com.matrix).

### Output
Community matrix with MxN dimensions where M (rows) are community entities (e.g. interaction modules), and N (cols) are trait variants. For each community, it stores the frequency of occurrence of a species trait variant in the set of species composing the community.
