import numpy as np
import seaborn as sns
import MDAnalysis as mda
import matplotlib.pyplot as plt
from collections import defaultdict 
from MDAnalysis.lib.distances import distance_array

#USER INPUTS
in_traj = ''
in_top = ''
outfile = ''
title = ''
cutoff = 4.0 #in angstroms

#prep for analysis
def prep():
    u = mda.Universe(in_top, in_traj)
    popx = u.select_atoms("resname POPX and not name H*")

    popx_resi = popx.residues
    popx_atoms = popx.atoms

    print(f"POPX molecules: {len(popx_resi)}")
    print(f"POPX atoms: {len(popx_atoms)}")

    contact_matrix = np.zeros((len(popx_atoms), len(popx_atoms)))
    return u, popx_resi, popx_atoms, contact_matrix

#calculate contacts
def calculate_popx_popx_contacts(u, popx_atoms, atom_to_lipid, contact_matrix):

    nframes = 0

    for ts in u.trajectory:
        nframes += 1

        dists = distance_array(
            popx_atoms.positions,
            popx_atoms.positions
        )

        # Exclude self-contacts and same-lipid contacts
        same_lipid = (
            atom_to_lipid[:, None] == atom_to_lipid[None, :]
        )
        atom_contacts = (dists < cutoff) & ~same_lipid

        contact_matrix += atom_contacts.astype(int)

    contact_freq = contact_matrix / nframes
    return contact_freq


#plot contacts