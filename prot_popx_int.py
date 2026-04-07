import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from MDAnalysis.lib.distances import distance_array

#USER INPUTS
#----------#
in_top = "nowat.tpr" # topology file (path) that matches trajectory (.tpr, .pdb, .gro)
in_traj = "../traj/cat_clus_pbc.xtc" # trajectory file (path) (.xtc, .trr)
title = "Hexamer-POPX Interactions (Rep 1)" # Title to go on figure
out_filename = "hex1_popx_int.png" # Name for figure .png file
cutoff = 4.0 # Distance cutoff for interactions in angstroms
#----------#

#DO NOT MODIFY BELOW THIS POINT

def prep():
    u=mda.Universe(in_top, in_traj)
    prot=u.select_atoms('protein')
    popx=u.select_atoms('resname POPX')

    prot_resi=prot.residues
    popx_atoms=popx.atoms

    n_res=len(prot_resi)
    n_popx=len(popx_atoms)

    print(n_res,n_popx)

    contact_matrix = np.zeros((n_res, n_popx))

    nframes=0

def calculate_contacts(u, prot_resi, popx_atoms, contact_matrix, nframes):
    for ts in u.trajectory:
        nframes += 1

        #Compute PROT x POPX distance matrix
        for i, res in enumerate(prot_resi):
            res_atoms=res.atoms.positions
            popx_pos = popx_atoms.positions

            dists=distance_array(res_atoms, popx_pos)

            # contacts per POPX atom
            contacts=np.any(dists < cutoff, axis=0)
            contact_matrix[i] += contacts.astype(int)

            # normalize contact frequencies
        contact_freq = contact_matrix / nframes
        np.save("protein_popx_contacts.npy", contact_freq)

def plot_contacts(prot_resi, popx_atoms, contact_freq):
    res_labels = [f'{res.resname}{res.resid}' for res in prot_resi]
    popx_labels = [f'{atom.name}{atom.index}' for atom in popx_atoms]

    plt.figure(figsize=(12,8))
    sns.heatmap(contact_freq,
                xticklabels=popx_labels,
                yticklabels=res_labels,
                cmap='magma',
                cbar_kws={'label': 'Contact Frequency'})
    
    plt.xlabel("Free POPC Atoms")
    plt.ylabel("Protein Residues")
    plt.title(title)

    plt.tight_layout()
    plt.savefig(out_filename, dpi=300, bbox_inches="tight")
    plt.close()


if __name__ == '__main__':
    prep()
    calculate_contacts(u, prot_resi, popx_atoms, contact_matrix, nframes)
    plot_contacts(prot_resi, popx_atoms, contact_freq)