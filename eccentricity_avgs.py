import MDAnalysis as mda
import numpy as np
from numpy.lib.stride_tricks import sliding_window_view
from MDAnalysis.analysis import align
from MDAnalysis.analysis.align import AlignTraj
import matplotlib.pyplot as plt
from matplotlib import colormaps
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize

u=mda.Universe("production_nowat.tpr", "cat_pbc_nowat_rep1.xtc")
ag=u.select_atoms("protein")

def align_trajectory():
    AlignTraj(u, u, select="protein", memoryview=True)

def calculate_eccentricity():
    eccentricity=[]
    for [i],ts in u.trajectory:  
        I=ag.moment_of_inertia()
        
        evals=np.linalg.eigvalsh(I)
        evals=np.sort(evals)

        Imin=evals[0]
        Iavg=np.mean(evals)

        e = 1-(Imin/Iavg)
        eccentricity.append(e)
    eccentricity = np.array(eccentricity)
    np.save('eccentricity_rep{i}.npy',eccentricity)

def compile_reps():
    ecc_files=['eccentricity_rep1.npy','eccentricity_rep2.npy','eccentricity_rep3.npy']
    ecc_stack=[]
    for fname in ecc_files:
        data=np.load(fname)
        smooth=sliding_window_view(data,window)
def plot_eccentricity():
    data = np.load("eccentricity.npy")

    smoothdata = sliding_window_view(data, window).mean(axis=1)

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

    ax.set_xlim(time_smooth.min(), time_smooth.max())
    ax.set_ylim(0, 1)
    ax.set_xlim(0,2000)
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Eccentricity")
    ax.set_title("Hexamer Eccentricity (Rep 1)")

    cbar = fig.colorbar(lc, ax=ax)
    cbar.set_label("Eccentricity")

    plt.savefig("eccentricity.png", dpi=300, bbox_inches="tight")
    plt.close()

if __name__ == "__main__":
    align_trajectory()
    calculate_eccentricity()
    plot_eccentricity()