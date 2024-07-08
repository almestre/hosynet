#!python
#!/usr/bin/env python

import matlab.engine 
import numpy as np
import champ
import matplotlib.pyplot as plt
import igraph as ig

eng1 = matlab.engine.start_matlab()
eng2 = matlab.engine.start_matlab()


# Use evalc to load the .mat file so that it is kept in the MATLAB Engine workspace.
# This stops the API from attempting to convert the variables Python types on load.
eng1.evalc("s = load('modtest_g4.mat');")
eng2.evalc("s = load('best_gamma.mat');")

# Now we can grab the data, as needed, from the MATLAB Engine workspace and pull it into Python.
# At this point the API will need to convert the data to Python types.
A_mat = np.asmatrix(eng1.eval("s.A_mat"));
P_mat = np.asmatrix(eng1.eval("s.P_mat"));
partition_array = np.asarray(eng2.eval("s.partition_array"));
Nmod20 = np.asarray(eng2.eval("s.Nmod20_gtest"));
Q = np.asarray(eng2.eval("s.Q_gtest"));
gamma = np.arange(0,6.01,0.01) # np.arange excludes the last element of the sequence

## Create the array of coefficients for the partitions
coeff_array=champ.champ_functions.create_coefarray_from_partitions(A_mat=A_mat,
   P_mat=P_mat,
   partition_array=partition_array)

# Create the ensemble object
ensemble1 = champ.PartitionEnsemble()
ensemble1._partitions = partition_array
ensemble1.int_edges = coeff_array[:,0]
ensemble1.exp_edges = coeff_array[:,1]
ensemble1.numparts = 601
ensemble1.orig_mods = Q*1020
ensemble1.resolutions = gamma
ensemble1.numcoms = Nmod20
ensemble1.min_com_size = 20

# Apply the champ algorithm to the partitions
ensemble1.apply_CHAMP(maxpt=6)

# Obtain the indices of the unique partitions
ensemble1.get_unique_partition_indices(reindex=True)

# Obtain the indices of the twin partitions
ensemble1.twin_partitions

# Plot modularity mapping
#ensemble1.plot_modularity_mapping(champ_only=True)
#ensemble1.plot_modularity_mapping()
#plt.show()

# Plot modularity domains
#champ.plot_single_layer_modularity_domains(ensemble1.ind2doms)
#plt.show()

# Plot modularity similarity heatmap
#champ.plot_similarity_heatmap_single_layer(ensemble1.partitions,ensemble1.ind2doms)
#plt.show()
