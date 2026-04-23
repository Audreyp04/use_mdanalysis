#Written by Audrey D. Prendergast
#Untested - intended for use with GROMACS mindist output (mindist.xvg)

import matplotlib.pyplot as plt
import numpy as np

#USER INPUTS
infile = 'mindist.xvg'
outfile = 'mindist.png'

from matplotlib.pyplot import minorticks_on
# Load data from the .xvg file, skipping lines starting with '#', '@', or '&'
try:
    x, y = np.loadtxt(infile, comments=("#", "@", "&"), unpack=True)
except FileNotFoundError:
    print("Error: File not found. Please check the file name and path.")
except ValueError:
    print("Error: Invalid data format in the file.")
else:

    # Plot the data
    time = x / 1000 #convert ps to ns
    plt.plot(time, y, ls = "-", color='C4')
    plt.xlim(left=0, right=2000)

    minorticks_on()
    plt.xlabel("Time (ns)", fontfamily='arial', fontsize=12)
    plt.ylabel("MinDist (nm)", fontfamily='arial', fontsize=12)
    plt.title("Minimum Distance of Protein to Membrane Atoms", fontfamily='arial', fontsize=12)
    plt.rcParams['font.family'] = 'arial'
    plt.rcParams['font.size'] = 12
    plt.savefig(outfile)
    plt.close()
