#################
# czscores
#################

#################
## DESCRIPTION: function to calculate the intermodular connectivity (c) and intramodular degree(z) for each node of a bipartite interaction network given a certain partition.

# We apply the same formulas as described in bipartite package.
# c = 1 - sum( (k.it/k.i)^2) # among-module connectivity = participation coefficient P in Guimerà & Amaral
# z = (k.is - ks.av) / SD.ks # within-module degree
# k.is = number of links of i to other species in its own module s; ks.av = average k.is of all species in module s; SD.ks = standard deviation of k.is of all species in module s; k.it = number of links of species i to module t; k.i = degree of species i Note that for any species alone (in its level) in a module the z-value will be NaN, since then SD.ks is 0. This is a limitation of the way the z-value is defined (in multiples of degree/strength standard deviations).

# Following Krasnov2012, we assigned z = 0 to cases wherein sd(z) = 0 (reasoning well explained in Krasnov2012).

#################
## Created by Alexandre Mestre
## Date of creation: 06/03/2020

#################
## Parameters
# matrix: MxN incidence matrix where M (i.e. rows) represent species of level one (e.g. hosts), and N (i.e. cols) reperesent species of the level 2 (e.g. symbionts).
# partition: vector of size M+N that stores the labels of module indentities of nodes, with species of level one first, followed by species of level two.

#################
## Output
# A list with two elements:
# c is a vector with the c values
# z is a vector with the z values

czscores <- function(matrix, partition, minModSize = NA)
{
   # Partition labels for rows and columns
   part <- as.factor(partition)
   partid <- levels(part)
   npart <- length(levels(part))
   part_row <- part[1:nrow(matrix)]
   ModSizes_row <- tapply(part_row, part_row, length)
   part_col <- part[-c(1:nrow(matrix))]
   ModSizes_col <- tapply(part_col, part_col, length)

   # filter data according to the minModSize
   if(!is.na(minModSize))
      {
         ModSizes <- tapply(partition, part, length)
         selModules <- levels(part)[ModSizes >= minModSize]
         matrix <- matrix[part_row %in% selModules, part_col %in% selModules]
         part <- partition[ part %in% selModules ]
         part <- as.factor(part)
         partid <- levels(part)
         npart <- length(levels(part))
         part_row <- part[1:nrow(matrix)]
         ModSizes_row <- tapply(part_row, part_row, length)
         part_col <- part[-c(1:nrow(matrix))]
         ModSizes_col <- tapply(part_col, part_col, length)
      }

   ############## cz values for species of first level (rows)
   #### c1
   ki1 <- rowSums (matrix)
   mod_counts1.ki <- matrix(rep(NA,nrow(matrix)*npart),nrow=nrow(matrix),ncol=npart)

   for (i in 1:npart)
      {
                    # Cases with only one species in the module
          if (ModSizes_col[i] > 1)
             {
                mod_counts1.ki[,i] <- (rowSums(matrix[,part_col == partid[i]])/ki1)^2
             }
          else
             {
                mod_counts1.ki[,i] <- (matrix[,part_col == partid[i]]/ki1)^2
             }
       }
   c1 <- 1-(rowSums(mod_counts1.ki))

   #### z1
   z1 <- rep (NA,nrow(matrix))
   mod_counts1 <- matrix(rep(NA,nrow(matrix)*npart),nrow=nrow(matrix),ncol=npart)

   for (i in 1:npart)
      {
          if (ModSizes_col[i] > 1)
             {
                mod_counts1[,i] <- rowSums(matrix[,part_col == partid[i]])
             }
          else
             {
                mod_counts1[,i] <- matrix[,part_col == partid[i]]
             }
       }

   for (i in 1:nrow(matrix))
      {
         si <- part_row[i] # identity of module s of species i
         kis <- mod_counts1[i, partid == si]
         ks_av <- mean(mod_counts1[part_row == si, partid == si])
         ks_sd <- sd(mod_counts1[part_row == si, partid == si])
         z1[i] <- (kis-ks_av)/ks_sd
      }

   ############## cz values for species of second level (cols)
   #### c2
   ki2 <- colSums (matrix)
   mod_counts2.ki <- matrix(rep(NA,ncol(matrix)*npart),nrow=ncol(matrix),ncol=npart)

   for (i in 1:npart)
      {
          if (ModSizes_row[i] > 1)
             {
                mod_counts2.ki[,i] = (colSums(matrix[part_row == partid[i],])/ki2)^2
             }
          else
             {
                mod_counts2.ki[,i] = (matrix[part_row == partid[i],]/ki2)^2
             }
       }
   c2 <- 1-(rowSums(mod_counts2.ki))

   #### z2
   z2 <- rep (NA,ncol(matrix))
   mod_counts2 <- matrix(rep(NA,ncol(matrix)*npart),nrow=ncol(matrix),ncol=npart)

   for (i in 1:npart)
      {
          if (ModSizes_row[i] > 1)
             {
                mod_counts2[,i] = colSums(matrix[part_row == partid[i],])
             }
          else
             {
                mod_counts2[,i] = matrix[part_row == partid[i],]
             }
       }

   for (i in 1:ncol(matrix))
      {
         si <- part_col[i] # identity of module s of species i
         kis <- mod_counts2[i, partid == si]
         ks_av <- mean(mod_counts2[part_col == si, partid == si])
         ks_sd <- sd(mod_counts2[part_col == si, partid == si])
         z2[i] <- (kis-ks_av)/ks_sd
      }
  
   c <- c(c1,c2)
   z <- c(z1,z2)

   # Replace NA and NaN by zero
   z[is.na(z)] <- 0
   z[z == NaN] <- 0

   return(list(c,z))
}
