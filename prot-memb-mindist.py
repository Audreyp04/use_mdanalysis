#Written by Audrey D. Prendergast
#NOT FUNCTIONAL - WIP

import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as plt

#USER INPUTS
in_traj = '../traj/cat_pbc_nowat.xtc'
in_top = 'nowat.tpr'
outfile = 'hex1_mid_dist.png'
title = 'Protein to Membrane Midplane Distances (Hexamer, Rep1)'

#This must be a matplotlib set color OR a HEX Code formated this way: '#ABC123'
color = 'darkslateblue'

#prep for analysis
def prep():
    u = mda.Universe(in_top, in_traj)
    prot = u.select_atoms("protein and name CA")
    memb = u.select_atoms("not protein and (name P N C2) and not (resname CL K POPX)")

    return u, prot, memb

#calculate contacts
def calculate_distances(u, prot, memb):

    distances = []
    times=[]

    for ts in u.trajectory:
        z_mid = memb.positions[:, 2].mean() #membrane midplane
        z_prot = prot.center_of_mass()[2] #prot center of mass 
        distances.append(abs(z_prot - z_mid) / 10.0)  # nm
        times.append(ts.time / 1000) #convert to ns 

    return np.array(times), np.array(distances)

#plot contacts
def plot_moiety_contacts(time, dist):

    plt.figure(figsize = (8, 8))
    plt.plot(time, dist, c=color)
    plt.xlabel("Time (ns)")
    plt.ylabel("Distance to Membrane Midplane(nm)")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.close()

if __name__ == "__main__":
    u, prot, memb = prep()
    time, dist = calculate_distances(u, prot, memb)
    plot_moiety_contacts(time, dist)