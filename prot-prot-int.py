import numpy as np
import seaborn as sns
import MDAnalysis as mda
import matplotlib.pyplot as plt
from collections import defaultdict 
from MDAnalysis.lib.distances import distance_array

#USER INPUTS
in_traj = '../traj/cat_pbc_nowat.xtc'
in_top = 'nowat.tpr'
outfile = 'hex1_ppi.png'
title = 'Protein-Protein Interactions (Hexamer, Rep1)'
cutoff = 4.0 #in angstroms

#prep for analysis
def prep():
    u = mda.Universe(in_top, in_traj)
    prot = u.select_atoms("protein")

    prot_resi = prot.residues
    prot_atoms = prot.atoms

    contact_matrix = np.zeros((len(prot_atoms), len(prot_atoms)))
    return u, prot_resi, prot_atoms, contact_matrix

#calculate contacts
def calculate_popx_popx_contacts(u, prot_resi, prot_atoms, contact_matrix):

    nframes = 0

    for ts in u.trajectory:
        nframes += 1

        dists = distance_array(
            prot_atoms.positions,
            prot_atoms.positions
        )

        # Exclude self-contacts
        same_chain = (
            prot_resi[:, None] == prot_resi[None, :]
        )
        atom_contacts = (dists < cutoff) & ~same_chain

        contact_matrix += atom_contacts.astype(int)

    contact_freq = contact_matrix / nframes
    return contact_freq

#plot contacts
def plot_moiety_contacts(prot_resi, contact_freq):

    plt.figure(figsize=(8, 8))
    sns.heatmap(
        contact_freq,
        xticklabels=prot_resi,
        yticklabels=prot_resi,
        cmap="magma",
        cbar_kws={"label": "Contact Frequency"},
        vmin=0, vmax=0.25
    )

    plt.xlabel("Protein Residue")
    plt.ylabel("Protein Residue")
    plt.title(title)

    plt.tight_layout()
    plt.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.close()

if __name__ == "__main__":
    u, prot_resi, prot_atoms, contact_matrix = prep()
    contact_freq = calculate_popx_popx_contacts(u, prot_resi, prot_atoms, contact_matrix)
    plot_moiety_contacts(prot_resi, contact_freq)