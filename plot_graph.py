#!python
#!/usr/bin/env python

import matlab.engine 
import numpy as np
import igraph as ig
from igraph import * # per importar tots els objectes de igraph i evitem emprar sempre el handle
from itertools import groupby

##############
## Function to convert nested list onto a flattened list
##########################################################
from collections import Iterable
def flatten(lis):
     for item in lis:
         if isinstance(item, Iterable) and not isinstance(item, basestring):
             for x in flatten(item):
                 yield x
         else:        
             yield item
#########################################################

# Load incidence matrix
eng = matlab.engine.start_matlab()

eng.evalc("s = load('modtest_g1.mat');");
A = np.asmatrix(eng.eval("s.A"));
partg = np.asarray(eng.eval("s.S_best"));

# Transmormem l'objecte que conté la info dels moduls en un format adecuat
# L'lobjecte importat era un array que contenia dos nested lists, i els tipus de dades eren 'float'
partg = partg.astype(int) # convertim en type int
partg = partg.tolist() # eliminem l'array extern
partg = list(flatten(partg)) # passem de nested lists a una flattened list



#######################################################################
#### First select the species that belong to the target modules
# The target modules are the modules with a minimum amount of nodes given the focal resolution
# In this case, for gamma = 1.3 the target modules are those with more than 19 nodes
#######################################################################

# Create vector with node counts per module
partg_sorted = np.sort(partg).tolist() # sort before using groupby!
Nnod = [len(list(group)) for key, group in groupby(partg_sorted)]

# Create vector with module identities for Nnod
idmod = range(1,len(Nnod)+1) # range elimina l'últim element del rang!

# Identify modules with more than 19 nodes
# First convert the lists into arrays to allow for arithmetic opretaions
Nnod = np.array(Nnod)
idmod = np.array(idmod)
selmod = idmod[Nnod>19]

# From partg (non-sorted version!) identify nodes belonging to the selected modules
selnod = np.isin (partg, selmod) # boleean vector of node selection

# Obtain filters to select host species (i.e. rows of the incidence matrix) and symbiont species (i.e. cols of the incidence matrix)
sel_hosts = selnod[:A.shape[0]]
sel_symbionts = selnod[A.shape[0]:]

# Select rows and cols (separately)
A_sel = A[sel_hosts]
A_sel = A_sel[:,sel_symbionts]

###############################
## Create the bipartite graph
###############################
# Represent the incidence matrix as a list
A_sel = A_sel.tolist()

g= ig.Graph()
g = g.Incidence(matrix=A_sel,directed=True,mode="OUT",multiple=False)

###############################
## Plot the bipartite graph
###############################

# Define vertex shapes based on organism type (host versus symbiont)
vshapes = np.repeat('square',sum(sel_hosts)).tolist()
vshapes2 = np.repeat("circle",sum(sel_symbionts)).tolist()
vshapes.extend(vshapes2)

# Define vertex sizes based on organism type (host versus symbiont)
vsizes = np.repeat(25,sum(sel_hosts)).tolist()
vsizes2 = np.repeat(20,sum(sel_symbionts)).tolist()
vsizes.extend(vsizes2)

# Define colors for vertex based on the number of modules
colors = ["green","grey","red","blue","orange","pink","yellow","black"]
colors = np.array(colors)

## RANDOM COLORS
#############################################
#from random import randint
#colors = []

#for i in range(len(selmod)):
#    colors.append('#%06X' % randint(0, 0xFFFFFF))

#colors = np.array(colors)
#############################################

# Create the array of colors for the nodes
partg_sel = np.array(partg)[selnod]

nod_colors = []

for i in range (len(partg_sel)):
  item = colors[selmod==partg_sel[i]]
  nod_colors.append(item[0].tolist())

out_fig_name = "graph_g1.eps"
visual_style = {}

# Set bbox and margin
visual_style["bbox"] = (2000,2000)
visual_style["margin"] = 50

# Set vertex colours
visual_style["vertex_color"] = nod_colors

# Set vertex size
visual_style["vertex_size"] = vsizes

# set vertex shape
visual_style["vertex_shape"] = vshapes

# Set vertex lable size
visual_style["vertex_label_size"] = 8

# Don't curve the edges
visual_style["edge_curved"] = False

# Set the layout
my_layout = g.layout_fruchterman_reingold()
visual_style["layout"] = my_layout

# Plot the graph
plot(g, out_fig_name, **visual_style)


