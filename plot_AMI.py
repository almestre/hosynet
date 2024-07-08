#!python
#!/usr/bin/env python

import matlab.engine
import numpy as np
import champ
import matplotlib.pyplot as plt
import igraph as ig
from sklearn.metrics import adjusted_mutual_info_score

eng = matlab.engine.start_matlab()

# Use evalc to load the .mat file so that it is kept in the MATLAB Engine workspace.
# This stops the API from attempting to convert the variables Python types on load.
eng.evalc("s = load('modtest_g1.mat');")


# Now we can grab the data, as needed, from the MATLAB Engine workspace and pull it into Python.
# At this point the API will need to convert the data to Python types.
part = np.asarray(eng.eval("s.part_g1"));
idx_best = np.asarray(eng.eval("s.idx_best"));
idx_best = int(idx_best)

# Select the best partition
best_part = part[idx_best-1] # python indices start from 0!!

# Remove the best partition from the part array of 100 replicates
# obj species the index; axis indicates if we remove rows (0) or cols (1)
part_rep = np.delete(arr=part, obj=idx_best-1, axis=0)

# Calculate AMI index best_part each part_rep[k]
AMI = [0 for i in range(99)] # list comprehension

for i in range(99):
   AMI[i] = adjusted_mutual_info_score(best_part,part_rep[i])

meanAMI = np.mean(AMI)
print(meanAMI)

stdAMI = np.std(AMI)
print(stdAMI)
