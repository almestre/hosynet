#################
# sort_dist_mat
#################

#################
## DESCRIPTION: function to sort a distance matrix (or a regular matrix)
# Source: https://rdrr.io/cran/dendextend/src/R/cor_cophenetic.R

sort_dist_mat <- function(dist_mat, by_rows = TRUE, by_cols = TRUE, ...) {
  dist1 <- as.matrix(dist_mat)
  if (by_rows) dist1 <- dist1[order(rownames(dist1)), ]
  if (by_cols) dist1 <- dist1[, order(colnames(dist1))]
  dist1 <- as.dist(dist1)
  attributes(dist1)[c("Diag", "Upper")] <- attributes(dist_mat)[c("Diag", "Upper")]
  if (class(dist_mat) == "matrix") {dist1 <- as.matrix(dist1)}
  dist1
}
