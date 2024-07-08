#!python
#!/usr/bin/env python

import best_g1 as bg1
ensemble1 = bg1.ensemble1

#import matlab.engine 
#import numpy as np
#import champ
#import igraph as ig
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patheffects as path_effects


# Plot modularity mapping

plot1 = ensemble1.plot_modularity_mapping(champ_only=True)

# Recover the twinx handle and reset the Twin axis label
plot1_twinx = plot1.get_shared_x_axes().get_siblings(plot1)[0]
plot1_twinx.remove()

# Reset axes
plot1.axes.clear()

# Reset parameters for mk3:
a2 = plot1.twinx()
a2.grid('off')
rand_ind=range(len(ensemble1.partitions))
allgams = [ensemble1.resolutions[ind] for ind in rand_ind]
allcoms = [ensemble1.numcoms[ind] for ind in rand_ind]
sct2 = a2.scatter(allgams, allcoms, marker="^", color="#91AEC1",
   alpha=1, label=r'Number of modules with $\geq 20$ nodes', zorder=1)

# Reset parameters for mk5:
plot1.plot(color="#FA8072")
#Get lists for the champ subset
champ_inds = ensemble1.get_CHAMP_indices()
# take the x-coord of first point in each domain
gammas=[ ensemble1.ind2doms[ind][0][0] for ind in champ_inds  ]
# take the y-coord of first point in each domain
mods = [ensemble1.ind2doms[ind][0][1] for ind in champ_inds]
mk5 = plot1.plot(gammas, mods, ls='--', color="green", lw=3, zorder=3)

# Reset parameters for mk4:
champ_coms = [ensemble1.numcoms[ind] for ind in champ_inds]
stp = a2.step(gammas, champ_coms, color="black", where='post', lw=2)
a2.set_ylabel(r'Number of modules with $\geq 20$ nodes')

# Reset parameters for mk2:
mk2 = plot1.scatter(gammas, mods, marker="v", color='blue', s=100, zorder=4)

# Set Axis lables
plot1.set_ylabel(u'Modularity (${Q}$)')
plot1.set_xlabel(u'Resolution (${\gamma}$)')


plot1.set_zorder(a2.get_zorder() + 1)  # put ax in front of ax2
plot1.patch.set_visible(False)  # hide the 'canvas'
plot1.set_xlim(xmin=0, xmax=max(allgams))

#### Set legend
# Fake for legend symbol
mk2 = plot1.scatter([], [], marker="v", color='blue', s=80, zorder=4)

# Fake for legend symbol
mk3 = plot1.scatter([], [], marker="^", color="#91AEC1", alpha=1,
						 label=r'Number of modules with $\geq 20$ nodes',
						 zorder=1,
						 s=50)
# Fake for legend symbol
mk4 = mlines.Line2D([], [], color='black', lw=2)

# Fake for legend symbol
mk5 = mlines.Line2D([], [], color='green', ls='--', lw=3)

# Fake for legend symbol
mk6 = mlines.Line2D([], [], color='#B22222', lw=1)


leg_set=[mk2, mk3, mk4, mk5, mk6]
leg_text=[ r'Transitions,$\gamma$',
   r'Number of modules with $\geq 20$ nodes',
   r'Number of modules with $\geq 20$ nodes (optimal)',   
   r'Convex hull of $Q(\gamma)$',
   'Selected low resolution']

plot1.legend(leg_set,
   leg_text,
   bbox_to_anchor=[0.68, .91],
   loc='center',
   frameon=True, fontsize=10)

# Add vertical line indicating the selected resolution
plt.axvline(x=1.3, color = "#B22222")

plt.show()

