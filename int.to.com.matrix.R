#################
# int.to.com.matrix
#################

#################
## Created by Alexandre Mestre
## Date of creation: 19/03/2020

#################
## DESCRIPTION: function to convert an incidence matrix for bipartite networks (i.e. sp_level1 x sp_level2) into a community matrix (modules x sp) based on a partition (i.e. a vector storing module indentities of the species).

#################
## Parameters
# matrix: incidence matrix with MxN dimensions where M (rows) are species of level 1 (e.g. hosts), and M (cols) are species of level 2 (e.g. symbionts)
# partition: vector that stores the identities of the interaction modules for the species belonging to the target level (either level 1 or 2). Elements should keep the same order than species labels in the matrix.
# level: If level = 1, we will obtain the nri for species of level 1 (i.e rows of the incidence matrix); if level = 2, we obtain nri for species of level 2 (columns the incidence matrix).

#################
## Output
# Community matrix with MxN dimensions where M (rows) are community entities (i.e. interaction modules), and N (cols) are species of the target level

int.to.com.matrix <- function( matrix , partition, level = 1)
   {
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
      return(com.matrix)
   }
