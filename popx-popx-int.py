#Written by Audrey D. Prendergast
#

import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from MDAnalysis.lib.distances import distance_array
from collections import defaultdict

# =========================
# USER INPUTS
# =========================
in_top = "../../nowat.tpr"
in_traj = "../traj/cat_pbc_nowat.xtc"
title = "Free Lipid Interactions (Rep 1)"
out_filename = "hex-popx-int.png"
cutoff = 4.0          # Å

# =========================
# PREP
# =========================
def prep():
    u = mda.Universe(in_top, in_traj)
    popx = u.select_atoms("resname POPX and not name H*")

    popx_resi = popx.residues
    popx_atoms = popx.atoms

    print(f"POPX molecules: {len(popx_resi)}")
    print(f"POPX atoms: {len(popx_atoms)}")

    contact_matrix = np.zeros((len(popx_atoms), len(popx_atoms)))
    return u, popx_resi, popx_atoms, contact_matrix


# =========================
# CONTACT CALCULATION
# =========================

def build_atom_to_lipid_map(popx_atoms):
    resid_indices = popx_atoms.resindices
    _, atom_to_lipid = np.unique(resid_indices, return_inverse=True)
    return atom_to_lipid

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

# =========================
# MOIETY COLLAPSING
# =========================
def build_moiety_map(popx_atoms):
    atom_to_moiety = []
    moiety_names = ["Headgroup", "Palmitoyl", "Oleoyl"]

    for atom in popx_atoms:
        name = atom.name

        # sn-1 chain (palmitoyl)
        if name in {"C3","C31","C32","C33","C34","C35" , 
                      "C36", "C37", "C38","C39","C310",
                      "C311","C312" , "C313", "C314", "C315",
                      "C316"}:
            atom_to_moiety.append("Palmitoyl")

        # sn-2 chain (oleoyl)
        elif name in {"C2","C21","C22","C23","C24","C25" , 
                      "C26", "C27", "C28","C29","C210",
                      "C211","C212" , "C213", "C214", "C215",
                      "C216", "C217", "C218"}:
            atom_to_moiety.append("Oleoyl")

        # Headgroup
        else:
            atom_to_moiety.append("Headgroup")

    return atom_to_moiety, moiety_names


def collapse_contacts_moiety_to_moiety(contact_freq_atom,atom_to_moiety,moiety_names):
    n_moiety = len(moiety_names)
    collapsed = np.zeros((n_moiety, n_moiety))

    for i, mi in enumerate(moiety_names):
        idx_i = [k for k, m in enumerate(atom_to_moiety) if m == mi]
        for j, mj in enumerate(moiety_names):
            idx_j = [k for k, m in enumerate(atom_to_moiety) if m == mj]

            if idx_i and idx_j:
                collapsed[i, j] = np.max(
                    contact_freq_atom[np.ix_(idx_i, idx_j)]
                )

    return collapsed

def write_moiety_table(popx_atoms, atom_to_moiety, outfile="moiety_assignments.dat"):
    with open(outfile, "w") as f:
        f.write("# resid resname atomname moiety\n")
        for atom, moiety in zip(popx_atoms, atom_to_moiety):
            f.write(
                f"{atom.resid:5d} {atom.resname:6s} "
                f"{atom.name:8s} {moiety}\n"
            )

# =========================
# PLOTTING
# =========================
def plot_moiety_contacts(moiety_names, collapsed_contacts):

    plt.figure(figsize=(8, 8))
    sns.heatmap(
        collapsed_contacts,
        xticklabels=moiety_names,
        yticklabels=moiety_names,
        cmap="magma",
        cbar_kws={"label": "Contact Frequency"},
        vmin=0, vmax=0.25
    )

    plt.xlabel("POPC Moiety")
    plt.ylabel("POPC Moiety")
    plt.title(title)

    plt.tight_layout()
    plt.savefig(out_filename, dpi=300, bbox_inches="tight")
    plt.close()

# =========================
# MAIN
# =========================
if __name__ == "__main__":

    u, popx_resi, popx_atoms, contact_matrix = prep()

    atom_to_lipid = build_atom_to_lipid_map(popx_atoms)
    
    contact_freq_atom = calculate_popx_popx_contacts(
        u, popx_atoms, atom_to_lipid, contact_matrix
    )

    atom_to_moiety, moiety_names = build_moiety_map(popx_atoms)

    write_moiety_table(popx_atoms, atom_to_moiety)

    moiety_contacts = collapse_contacts_moiety_to_moiety(
        contact_freq_atom,
        atom_to_moiety,
        moiety_names
    )

    plot_moiety_contacts(
        moiety_names, moiety_contacts
    )

    print("✅ POPX–POPX moiety interaction analysis complete.")
