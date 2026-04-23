#Written by Audrey D. Prendergast
#NOT FUNCTIONAL - WIP

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

    contact_matrix = np.zeros((len(prot_resi), len(prot_resi)))
    return u, prot_resi, prot_atoms, contact_matrix

#calculate contacts
def calculate_residue_contacts(u, prot_resi, prot_atoms, contact_matrix):
    nframes = 0

    atom_to_res = prot_atoms.resindices     # shape (N_atoms,)
    atom_segids = prot_atoms.segids         # chain IDs

    n_res = len(prot_resi)

    for ts in u.trajectory:
        nframes += 1

        # atom–atom distances
        dists = distance_array(prot_atoms.positions, prot_atoms.positions)

        # exclude self contacts
        np.fill_diagonal(dists, np.inf)

        # contact mask
        contact_mask = dists < cutoff

        # exclude intra‑chain contacts (protein–protein only)
        inter_chain = atom_segids[:, None] != atom_segids[None, :]
        contact_mask &= inter_chain

        # atom indices forming contacts
        ai, aj = np.where(contact_mask)

        # map atoms to residues
        ri = atom_to_res[ai]
        rj = atom_to_res[aj]

        # accumulate contacts into residue matrix
        np.add.at(contact_matrix, (ri, rj), 1)

    contact_freq = contact_matrix / nframes
    return contact_freq

#plot contacts
def plot_moiety_contacts(prot_resi, contact_freq):
    labels = [f"{res.resname}{res.resid}" for res in prot_resi]
    plt.figure(figsize=(8, 8))
    sns.heatmap(
        contact_freq,
        xticklabels=labels,
        yticklabels=labels,
        cmap="magma",
        cbar_kws={"label": "Contact Frequency"},
        vmin=0, vmax=1
    )
    plt.xlabel("Residue")
    plt.ylabel("Residue")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outfile, dpi=300)
    plt.close()

if __name__ == "__main__":
    u, prot_resi, prot_atoms, contact_matrix = prep()
    contact_freq = calculate_residue_contacts(u, prot_resi, prot_atoms, contact_matrix)
    plot_moiety_contacts(prot_resi, contact_freq)