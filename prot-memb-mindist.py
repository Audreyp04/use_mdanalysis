#Written by Audrey D. Prendergast
#NOT FUNCTIONAL - WIP

import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as plt
from collections import defaultdict
from MDAnalysis.lib.distances import distance_array


#USER INPUTS
in_traj = '../traj/cat_pbc_nowat.xtc'
in_top = 'nowat.tpr'
outfile = 'hex1_mindist.png'
title = 'Protein to Membrane Minimum Distances (Hexamer, Rep1)'

#This must be a matplotlib set color OR a HEX Code formated this way: '#ABC123'
color = 'darkslateblue'

#prep for analysis
def prep():
    u = mda.Universe(in_top, in_traj)
    prot = u.select_atoms("protein")
    memb = u.select_atoms("not protein and not resname CL and not resname K and not resname POPX")

    return u, prot, memb

#calculate contacts
def calculate_distances(u, prot, memb):

    nframes = 0
    mindist = []
    for ts in u.trajectory:
        nframes += 1

        dists = distance_array(
            prot.positions,
            memb.positions
        )
        mindist.append(np.min(dists))

    return mindist, nframes

#plot contacts
def plot_moiety_contacts(mindist, nframes):

    time = nframes / 10

    plt.plot(time, mindist, c=color)
    plt.figure(figsize=(8, 8))
    plt.xlabel("Time (ns)")
    plt.ylabel("Minimum Distance (nm)")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.close()

if __name__ == "__main__":
    u, prot, memb = prep()
    mindist = calculate_distances(u, prot, memb)
    plot_moiety_contacts(mindist)