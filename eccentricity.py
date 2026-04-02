import MDAnalysis as mda
import numpy as np
from numpy.lib.stride_tricks import sliding_window_view
from MDAnalysis.analysis import align
import matplotlib.pyplot as plt
from matplotlib import colormaps

u=mda.Universe("nowat.tpr", "../traj/cat_pbc.xtc")
ag=u.select_atoms("protein")

def calculate_eccentricity():
    eccentricity=[]
    for ts in u.trajectory:
        I=ag.moment_of_inertia()
        
        evals=np.linalg.eigvalsh(I)
        evals=np.sort(evals)

        Imin=evals[0]
        Iavg=np.mean(evals)

        e = 1-(Imin/Iavg)
        eccentricity.append(e)
    eccentricity = np.array(eccentricity)
    np.save('eccentricity.npy',eccentricity)

def plot_eccentricity():
    data=np.load("eccentricity.npy")
    smoothdata=sliding_window_view(data, window_shape=5).mean()
    time_ns = np.arange(len(data))/10
    plt.figure(figsize=(6, 4))
    plt.plot(time_ns, smoothdata, color='black')
    plt.xlabel('Time (ns)')
    plt.ylabel('Eccentricity')
    plt.xlim(0,2000)
    plt.ylim(0, 1)
    plt.title('Oligomer Eccentricity vs Time')

    plt.savefig('eccentricity.png', bbox_inches='tight', dpi=300)
    plt.close()


if __name__ == "__main__":
    calculate_eccentricity()
    plot_eccentricity()