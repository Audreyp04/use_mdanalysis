#Written by Audrey D. Prendergast 
#Fully functional as of 4/6/2026

import MDAnalysis as mda
import numpy as np
from numpy.lib.stride_tricks import sliding_window_view
from MDAnalysis.analysis import align
from MDAnalysis.analysis.align import AlignTraj
import matplotlib.pyplot as plt
from matplotlib import colormaps
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize


# =========================
# GLOBALS (initialized via init)
# =========================
u = None
ag = None
name = None
out_filename = None
window = 5


def init(*, top, traj, out, title):
    global u, ag, name, out_filename

    name = f"Eccentricity - {title}"
    out_filename = f"{out}_ecc.png"

    u = mda.Universe(top, traj)
    ag = u.select_atoms("protein")


def align_trajectory():
    align=AlignTraj(u, u, select="protein")
    align.run()

def calculate_eccentricity():
    eccentricity=[]
    for ts in u.trajectory:  
        p=ag.moment_of_inertia()
        e1,e2,e3 = np.linalg.eigvalsh(p)
        etop=e1+e2-e3
        ebot=-e1+e2+e3
        e = np.sqrt((1 - (etop / ebot)))
        eccentricity.append(e)
    eccentricity=np.array(eccentricity)
    print(eccentricity)
    np.save('eccentricity.npy',eccentricity)

def plot_eccentricity():
    data = np.load("eccentricity.npy")

    smoothdata = sliding_window_view(data,window).mean(axis=1)

    # Time axis correction to correspond to smoothed data
    time_ns = np.arange(len(data)) / 10.0
    time_smooth = time_ns[:len(smoothdata)]

    #color the line by eccentricity value
    points = np.array([time_smooth, smoothdata]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    norm = Normalize(vmin=0.0, vmax=1.0)
    lc = LineCollection(segments, cmap="magma", norm=norm)
    lc.set_array(smoothdata)
    lc.set_linewidth(1)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.add_collection(lc)

    ax.set_ylim(0, 1)
    ax.set_xlim(0,2000)
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Eccentricity")
    ax.set_title(name)

    cbar = fig.colorbar(lc, ax=ax)
    cbar.set_label("Eccentricity")

    plt.savefig(out_filename, dpi=300, bbox_inches="tight")
    plt.close()

if __name__ == "__main__":
    align_trajectory()
    calculate_eccentricity()
    plot_eccentricity()