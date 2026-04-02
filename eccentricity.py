import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis import align
import matplotlib.pyplot as plt
from matplotlib import colormaps

u=mda.Universe("nowat.tpr", "../traj/cat_pbc.xtc")
ag=u.select_atoms("protein")

def calculate_eccentricity():
    eccentricity=[]
    for ts in u.trajectory:
        I=ag.moment_of_inertia(weights='mass')
        
        evals=np.linalg.eigvals(I)
        evals=np.sort(evals)

        Imin=evals[0]
        Iavg=np.mean(evals)

        e = 1-(Imin/Iavg)
        eccentricity.append(e)
        np.save('eccentricity.npy',eccentricity)

def plot_eccentricity():
    data=np.load("eccentricity.npy")
    plt.imshow(data,cmap='magma')
    plt.colorbar('Oligomer Eccentricity')
    plt.xlabel('Time (µs)')
    plt.ylabel('Eccentricity')
    plt.yscale(min=0,max=1)
    plt.savefig('eccentricity.png', bbox_inches='tight')
