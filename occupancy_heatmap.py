#Written by Audrey D. Prendergast 
#Fully functional as of 4/6/2026
import MDAnalysis as mda
import numpy as np
from numpy.lib.stride_tricks import sliding_window_view
from MDAnalysis.analysis.align import AlignTraj
import matplotlib.pyplot as plt

# Change only these items
replicates=[1,2,3]
in_top = "top.nowat.pdb" # topology file (path) that matches trajectory (.tpr, .pdb, .gro)
traj_template = "rep{}/cat.pbc.nowat.xtc"
title = "Decamer Eccentric Occupancy (w/o Free Lipids)" # Title to go on figure
out_filename = "wolip_dec_occupancy.png" # Name for figure .png file
window = 5 #window size for smoothed data, 50 is usually good

#CHANGE NOTHING BELOW THIS LINE
#------------------------------
def calculate_eccentricity_all_reps():
    for i in replicates:
        print(f"Processing replicate {i}")
        
        #load in trajectories
        traj_file = traj_template.format(i)
        u = mda.Universe(in_top, traj_file)
        ag = u.select_atoms("protein")

        # Align whole trajectory
        aligner = AlignTraj(u, u, select="protein")
        aligner.run()

        ecc = []

        for ts in u.trajectory:  
            p=ag.moment_of_inertia()
            e1,e2,e3 = np.linalg.eigvalsh(p)
            etop=e1+e2-e3
            ebot=-e1+e2+e3
            e = np.sqrt((1 - (etop / ebot)))
            ecc.append(e)
        ecc=np.array(ecc)
        np.save(f"eccentricity_rep{i}.npy",ecc)
        print(f"  Saved eccentricity_rep{i}.npy")


def plot_eccentricity():

    ecc_stack = []

    for i in replicates:
        data = np.load(f"eccentricity_rep{i}.npy")
        smooth = sliding_window_view(data, window).mean(axis=1)
        ecc_stack.append(smooth)

    ecc_stack = np.array(ecc_stack)
    mean_ecc = ecc_stack.mean(axis=0)[np.newaxis, :]  # shape: (1, time)

    dt = 0.1  # ns per frame
    tmax = mean_ecc.shape[1] * dt

    plt.figure(figsize=(7,2))

    im = plt.imshow(
        mean_ecc,
        aspect="auto",
        origin="lower",
        cmap="magma",
        vmin=0.0,
        vmax=1.0,
        extent=[0, tmax, 0, 1]
    )

    plt.colorbar(im, label="Eccentricity")
    plt.xlabel("Time (ns)")
    plt.yticks([])
    plt.title(title + " (Mean over replicates)")

    plt.xlim(0, 2000)

    plt.savefig(out_filename, dpi=300, bbox_inches="tight")
    plt.close()

if __name__ == "__main__":
    calculate_eccentricity_all_reps()
    plot_eccentricity()