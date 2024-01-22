############################
## ses.pd.fixed
############################

############################
### Description:
# Corrected function of ses.pd function of package picante.
# Source:
#https://www.biorxiv.org/content/10.1101/579300v1.supplementary-material

ses.pd.fixed<-function(samp, tree, null.model = c("taxa.labels", "richness",
"frequency", "sample.pool", "phylogeny.pool", "independentswap", "trialswap"),
runs = 999, iterations = 1000, include.root = TRUE)
{
if(include.root == TRUE) {
pd.obs <- as.vector(pd(samp, tree, include.root = TRUE)$PD)
null.model <- match.arg(null.model)
pd.rand <- switch(null.model, taxa.labels = t(replicate(runs,
as.vector(pd(samp, tipShuffle(tree), include.root = TRUE)$PD))), richness =
t(replicate(runs,as.vector(pd(randomizeMatrix(samp, null.model = "richness"),
tree, include.root = TRUE)$PD))), frequency = t(replicate(runs, as.vector
(pd(randomizeMatrix(samp,null.model = "frequency"), tree, include.root =
TRUE)$PD))), sample.pool = t(replicate(runs,as.vector(pd(randomizeMatrix(samp,
null.model = "richness"), tree, include.root = TRUE)$PD))), phylogeny.pool =
t(replicate(runs, as.vector(pd(randomizeMatrix(samp, null.model =
"richness"),tipShuffle(tree), include.root = TRUE)$PD))), independentswap =
t(replicate(runs,as.vector(pd(randomizeMatrix(samp, null.model =
"independentswap", iterations), tree, include.root = TRUE)$PD))), trialswap =
t(replicate(runs,as.vector(pd(randomizeMatrix(samp, null.model = "trialswap",
iterations), tree, include.root = TRUE)$PD))))
pd.rand.mean <- apply(X = pd.rand, MARGIN = 2, FUN = mean, na.rm = TRUE)
pd.rand.sd <- apply(X = pd.rand, MARGIN = 2, FUN = sd, na.rm = TRUE)
pd.obs.z <- (pd.obs - pd.rand.mean)/pd.rand.sd
pd.obs.rank <- apply(X = rbind(pd.obs, pd.rand), MARGIN = 2, FUN = rank)[1, ]
pd.obs.rank <- ifelse(is.na(pd.rand.mean), NA, pd.obs.rank)
Result<-data.frame(ntaxa = specnumber(samp), pd.obs, pd.rand.mean,
pd.rand.sd, pd.obs.rank, pd.obs.z, pd.obs.p = pd.obs.rank/(runs + 1), runs =
runs, row.names = row.names(samp))
return(Result)
}
if(include.root == FALSE) {
pd.obs <- as.vector(pd(samp, tree, include.root = FALSE)$PD)
null.model <- match.arg(null.model)
pd.rand <- switch(null.model, taxa.labels = t(replicate(runs,
as.vector(pd(samp, tipShuffle(tree), include.root = FALSE)$PD))), richness =
t(replicate(runs,as.vector(pd(randomizeMatrix(samp, null.model = "richness"),
tree, include.root = FALSE)$PD))), frequency = t(replicate(runs, as.vector 
(pd(randomizeMatrix(samp,null.model = "frequency"), tree, include.root =
FALSE)$PD))),sample.pool = t(replicate(runs,as.vector(pd(randomizeMatrix(samp,
null.model = "richness"), tree, include.root = FALSE)$PD))), phylogeny.pool =
t(replicate(runs, as.vector(pd(randomizeMatrix(samp, null.model =
"richness"),tipShuffle(tree), include.root = FALSE)$PD))), independentswap =
t(replicate(runs,as.vector(pd(randomizeMatrix(samp, null.model =
"independentswap", iterations), tree, include.root = FALSE)$PD))), trialswap =
t(replicate(runs,as.vector(pd(randomizeMatrix(samp, null.model = "trialswap",
iterations), tree, include.root = FALSE)$PD))))
pd.rand.mean <- apply(X = pd.rand, MARGIN = 2, FUN = mean, na.rm = FALSE)
pd.rand.sd <- apply(X = pd.rand, MARGIN = 2, FUN = sd, na.rm = FALSE)
pd.obs.z <- (pd.obs - pd.rand.mean)/pd.rand.sd
pd.obs.rank <- apply(X = rbind(pd.obs, pd.rand), MARGIN = 2, FUN = rank)[1, ]
pd.obs.rank <- ifelse(is.na(pd.rand.mean), NA, pd.obs.rank)
Result<-data.frame(ntaxa = specnumber(samp), pd.obs, pd.rand.mean,
pd.rand.sd, pd.obs.rank, pd.obs.z, pd.obs.p = pd.obs.rank/(runs + 1), runs =
runs, row.names = row.names(samp))
return(Result)
}
}
