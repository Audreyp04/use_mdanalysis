#Written by Audrey D. Prendergast
#work in progress

import numpy as np
import seaborn as sns
import MDAnalysis as mda
import matplotlib.pyplot as plt
from collections import defaultdict 
from MDAnalysis.lib.distances import distance_array

#USER INPUTS
in_traj = '../traj/cat_pbc_nowat.xtc'
in_top = 'nowat.tpr'
outfile = 'hex1_resi_ppi.png'
title = 'Residue Averaged Contact Frequency (Hexamer, Rep1)'
cutoff = 4.0 #in angstroms
residues_per_chain = 42

#prep for analysis
def prep():
    u = mda.Universe(in_top, in_traj)
    prot = u.select_atoms("protein")
    residues = prot.residues
    nchains = len(residues) // residues_per_chain
    chains = []

    for i in range(nchains):
        start = i * residues_per_chain
        end = (i + 1) * residues_per_chain
        chains.append(residues[start:end])

    residue_labels = [
        f"{res.resname}{res.resid}" for res in chains[0]
    ]

    return u, chains, residue_labels

#calculate contacts
def calculate_residue_contacts(u, chains):
    nchains = len(chains)
    nframes = 0

    contact_counts = np.zeros((residues_per_chain, residues_per_chain), dtype=np.float64)

    chain_atom_resindex = []
    for chain in chains:
        atom_resindex = chain.atoms.resindices % residues_per_chain
        chain_atom_resindex.append(atom_resindex)

    for ts in u.trajectory:
        nframes += 1
        for c1 in range(nchains):
            for c2 in range(c1 + 1, nchains):
                ag1 = chains[c1].atoms
                ag2 = chains[c2].atoms
                dmat = distance_array(ag1.positions, ag2.positions)
                contact_mat = dmat < cutoff
                if not np.any(contact_mat):
                    continue

                r1 = chain_atom_resindex[c1]
                r2 = chain_atom_resindex[c2]

                atom_i, atom_j = np.nonzero(contact_mat)
                res_i = r1[atom_i]
                res_j = r2[atom_j]
                unique_pairs = np.unique(
                    np.stack((res_i, res_j), axis=1),
                    axis=0
                )
                np.add.at(contact_counts,
                          (unique_pairs[:, 0], unique_pairs[:, 1]),
                          1)
                np.add.at(contact_counts,
                          (unique_pairs[:, 1], unique_pairs[:, 0]),
                          1)

    n_chain_pairs = nchains * (nchains - 1) / 2
    contact_freq = contact_counts / (nframes * n_chain_pairs)
    return contact_freq

#plot contacts
def plot_moiety_contacts(contact_freq, residue_labels):

    plt.figure(figsize=(8, 8))
    sns.heatmap(
        contact_freq,
        xticklabels=residue_labels,
        yticklabels=residue_labels,
        cmap="magma",
        cbar_kws={"label": "Contact Frequency"},
        vmin=0, vmax=0.5
    )
    plt.xlabel("Residue")
    plt.ylabel("Residue")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outfile, dpi=300)
    plt.close()

if __name__ == "__main__":
    u, chains, residue_labels = prep()
    print("Prep complete!")
    contact_freq = calculate_residue_contacts(u, chains)
    print("Contact calculations complete!")
    plot_moiety_contacts(contact_freq, residue_labels)
    print("Plotting complete!")