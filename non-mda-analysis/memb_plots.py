#For use with Marion Q. LoPresti's membrane analysis 
#Written by Audrey D. Prendergast

import numpy as np
import matplotlib.pyplot as plt

# =====================
# filenames (edit if needed)
# =====================
THICKNESS_FILE = "thickness.hydrophobic.mat.npy"
NEIGHBOR_FILE  = "neighbor_density_frames.npy"

thick_title = "hex1_thickness.png"
dens_title = "hex1_density.png"

SYSTEM_NAME = "Hexamer with Free Lipids"
REP = 1

# =====================
# load data
# =====================
# thickness: shape (time, x, y)
thickness_mat = np.load(THICKNESS_FILE)
thickness = np.nanmean(thickness_mat, axis=0)

# neighbor density: list[dict], dict["top"] -> (x, y)
neighbor_frames = np.load(NEIGHBOR_FILE, allow_pickle=True)
neighbor_stack = np.stack([f["top"] for f in neighbor_frames])
neighbor_density = np.nanmean(neighbor_stack, axis=0)

# =====================
# plotting helper
# =====================
def plot_map(data,filename, title, cmap, cbar_label, vmin, vmax):
    plt.figure(figsize=(6, 6))
    im = plt.imshow(
        data.T,
        origin="lower",
        cmap=cmap,
        aspect="equal",
        vmin=vmin,
        vmax=vmax
    )
    plt.xlabel("Grid X")
    plt.ylabel("Grid Y")
    plt.title(title)
    plt.colorbar(im, label=cbar_label)
    plt.tight_layout()
    plt.savefig(filename)

# =====================
# plots
# =====================
plot_map(
    thickness,thick_title,
    title=f"{SYSTEM_NAME} – Rep {REP}\nAverage Membrane Thickness",
    cmap="magma",
    cbar_label="Thickness (nm)",
    vmin=4.55,
    vmax=4.65
)

plot_map(
    neighbor_density,dens_title,
    title=f"{SYSTEM_NAME} – Rep {REP}\nAverage Neighbor Density (Top Leaflet)",
    cmap="magma",
    cbar_label="Neighbor count",
    vmin=4,
    vmax=5
)
