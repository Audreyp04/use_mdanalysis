#Written by Audrey D. Prendergast
#Fully functional as of 4/24/2026

import numpy as np
import pandas as pd
import MDAnalysis as mda
import matplotlib.pyplot as plt


# =========================
# GLOBALS
# =========================
in_top = None
in_traj = None
outfile = None
name = None
window = 50
color = "darkslateblue"


def init(*, top, traj, out, title):
    global in_top, in_traj, outfile, name

    in_top = top
    in_traj = traj
    outfile = f"{out}_mid_dist.png"
    name = f"Protein to Membrane Midplane Distances - {title}"

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
    time_df = pd.DataFrame(time)
    dist_df = pd.DataFrame(dist)
    time_smooth = time_df.rolling(window).mean()
    dist_smooth = dist_df.rolling(window).mean()
    plt.figure(figsize = (8, 8))
    plt.plot(time_smooth, dist_smooth, c=color)
    plt.xlabel("Time (ns)")
    plt.ylabel("Distance to Membrane Midplane (nm)")
    plt.xlim(0,2000)
    plt.ylim(0,16)
    plt.title(name)
    plt.tight_layout()
    plt.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.close()

if __name__ == "__main__":
    u, prot, memb = prep()
    print("Prep complete")
    time, dist = calculate_distances(u, prot, memb)
    print("Distance calculations complete")
    plot_moiety_contacts(time, dist)
    print("Plotting complete")