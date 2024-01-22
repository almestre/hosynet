#################
# cast_to_matrix
#################

#################
## DESCRIPTION: function to build a matrix of pairwise associations of two "character" variables stored in a dataframe. Examples of applications:
# 1. Building a site x species matrix (presence/absences)
# 2. Building a co-occurrence matrix for two groups of species such as symbionts and hosts (number of co-occurrences per pairwise association)

### Parameters
# dataframe: dataframe that stores the variables to be casted
# variable1: variable that will conform the rows of the output matrix
# variable2: variable that will conform the columns of the output matrix
# pa: if pa = TRUE, the output matrix is a presences/absences matrix (0 and 1 data). If pa=FALSE (default option), it informs about the number of co-occurrences (0 and positive integers).

cast_to_matrix <- function(dataframe, variable1, variable2, pa=FALSE)
{
# 1. We build the matrix with empty data.
# 1.1 Obtain rownames and colnames (sorted), and numbers of rows

  # variables
  var1 <- dataframe[,variable1]
  var2 <- dataframe[,variable2]

  # merging variables:
  merge_var <- paste(var1,var2)

  # rownames and colnames of the output matrix
  mrow <- sort(unique(var1))
  mcol <- sort(unique(var2))

  # number of rows and cols of the matrix
  lmrow <- length(mrow)
  lmcol <- length(mcol)

  # Data on var1-var2 pairwise combinations, stored in a matrix the same structure than the outcome
  mrow_data <- matrix(rep(mrow,lmcol),nrow=lmrow,ncol=lmcol, byrow=FALSE)
  mcol_data <- matrix(rep(mcol,lmrow),nrow=lmrow,ncol=lmcol, byrow=TRUE)
  pw_comb <- paste(mrow_data,mcol_data)
  dim(pw_comb) <- c(lmrow,lmcol)

  # buiding the output matrix
  m_output <- sapply(pw_comb,function(x) sum(merge_var==x)) ## it returns a vector
  dim(m_output) <- c(lmrow,lmcol) # adding dimensions to get a matrix
  rownames(m_output) <- mrow
  colnames(m_output) <- mcol

  if (pa==TRUE)
    {
      m_output <- ifelse (m_output > 0, 1,0) # ifelse (test,true,false)
    }

  return(m_output)
}
