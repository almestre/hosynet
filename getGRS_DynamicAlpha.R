########################
## getGRS_DynamicAlpha
########################
# Function to estimate geographic range size (GRS) of a species from occurrence data.

# Methods:
#Extent of occurrence (EOO) based on alpha-convex hull calculated with dynamic alpha (see getDynamicAlphaHull function from packagage rangeBuilder)
#Area of occupancy (AOO) calculated with the conR package

# Input:
# "data" is a dataframe with the occurrence data. It should have three columns (it is mandatory to respect column positions):
   #The first column is the coordx
   #The second column is the coordy
   #The third column is the species names
# "sp"
   #Character value indicating the name of the target species for GRS calculations
# "Cell_size_AOO"
   #Numeric value indicating the grid size (km) used for AOO estimations. By default, equal to 2

# Output:
#Object of class SpatialPolygons storing the alpha-convex hull
#EOO in km² 
#AOO in km²

#### Required packages, functions and data:
#library("rangeBuilder")
#library("ConR")
#library("maptools")
#library("alphahull")
#library("sp")
#source("functions/GRS_functions/ahull2lines.R")
#library("rgeos") # per utilitzar gIntersection
#library("geosphere") # to calculate geographic area from a SpatialPolygons object
#source("functions/GRS_functions/ah2sp.R")

getGRS_DynamicAlpha <- function(data, sp, Cell_size_AOO = 2)
   {
      ## 1. List to store the output
      output <- list()
      ## 2. Data selection
      selec <- data[,3] == sp
      coords <- unique(data[selec,c(1,2)])
      x <- coords[,1]
      y <- coords[,2]

      if (length(x) >= 3)
         {
            # Calculation of the alpha-hull polygon
            xy <- data.frame("x" = x,"y" = y)
            GR <- getDynamicAlphaHull(
               x = xy,
               fraction = 0.95,
               partCount = 1,
               buff = 5000,
               initialAlpha = 3,
               clipToCoast = "terrestrial"
               )
            GR <- GR[[1]]
            # GRS calculations
            EOO <- areaPolygon(GR)/1000000 # default value is in square meters!
            data_AOO <- data.frame("ddlat" = y,"ddlon" = x,"tax" = rep(sp,length(x)))
            AOO <- AOO.computing(data_AOO, Cell_size_AOO = Cell_size_AOO)
            output[[1]] <- GR
            output[[2]] <- EOO
            output[[3]] <- AOO
         }
      else if (length(x) < 3)
         {
            data_AOO <- data.frame("ddlat" = y,"ddlon" = x,"tax" = rep(sp,length(x)))
            AOO <- AOO.computing(data_AOO, Cell_size_AOO = Cell_size_AOO)
            output[[1]] <- NA
            output[[2]] <- NA
            output[[3]] <- AOO
         }
      else
         {
            output[[1]] <-NA
            output[[2]] <-NA
            output[[3]] <-NA
         }

      return(output)
   }
