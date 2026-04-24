#Written by Audrey D. Prendergast
#Fully functional as of 4/24/2026

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
title = 'Chain-Chain Contact Frequency (Hexamer, Rep1)'
cutoff = 4.0 #in angstroms
residues_per_chain = 42

#prep for analysis
def prep():
    u = mda.Universe(in_top, in_traj)
    prot = u.select_atoms("protein")
    prot_atoms = prot.atoms
    residues = prot_atoms.residues
    nchains = len(residues) // residues_per_chain
    chains = {}

    for i in range(nchains):
        start = i * residues_per_chain
        end = (i + 1) * residues_per_chain

        chain_residues = residues[start:end]
        chains[f"Chain{i+1}"] = chain_residues.atoms

    return u, chains

#calculate contacts
def calculate_chain_contacts(u, chains):
    chain_ids = list(chains.keys())
    n_chains = len(chain_ids)
    nframes = 0

    contact_counts = np.zeros((n_chains, n_chains))

    for ts in u.trajectory:
        nframes += 1
        for i in range(n_chains):
            for j in range(i+1, n_chains):
                ai = chains[chain_ids[i]].positions
                aj = chains[chain_ids[j]].positions

                dists = distance_array(ai, aj)

                if np.any(dists < cutoff):
                    contact_counts[i,j] += 1
                    contact_counts[j,i] += 1

    contact_freq = contact_counts / nframes
    return chain_ids, contact_freq

#plot contacts
def plot_moiety_contacts(chain_ids, contact_freq):

    plt.figure(figsize=(8, 8))
    sns.heatmap(
        contact_freq,
        xticklabels=chain_ids,
        yticklabels=chain_ids,
        cmap="magma",
        cbar_kws={"label": "Contact Frequency"},
        vmin=0, vmax=1
    )
    plt.xlabel("Chain")
    plt.ylabel("Chain")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outfile, dpi=300)
    plt.close()

if __name__ == "__main__":
    u, chains = prep()
    chain_ids, contact_freq = calculate_chain_contacts(u, chains)
    plot_moiety_contacts(chain_ids, contact_freq)