#Written by Audrey D. Prendergast
#Fully functional as of 4/15/2026

import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from MDAnalysis.lib.distances import distance_array

# =========================
# USER INPUTS
# =========================
in_top = "nowat.tpr"
in_traj = "../traj/cat_pbc_nowat.xtc"
title = "Hexamer & Free Lipid Interactions (Rep 1)"
out_filename = "hex1_popx_moiety_contacts.png"
cutoff = 4.0          # Å
subunit_length = 42
threshold=0.80

# =========================
# PREP
# =========================
def prep():
    u = mda.Universe(in_top, in_traj)
    prot = u.select_atoms("protein")
    popx = u.select_atoms("resname POPX and not name H*")

    prot_resi = prot.residues
    prot_resids = prot.residues.resids 
    popx_atoms = popx.atoms

    print(f"Protein residues: {len(prot_resi)}")
    print(f"POPX atoms: {len(popx_atoms)}")

    contact_matrix = np.zeros((len(prot_resi), len(popx_atoms)))
    return u, prot_resi, prot_resids, popx_atoms, contact_matrix


# =========================
# CONTACT CALCULATION
# =========================
def build_atom_to_residue_map(prot_resi):
    
    #Returns an array mapping each protein atom index
    #to its residue index (0..n_res-1)
    
    atom_to_res = np.empty(len(prot_resi.atoms), dtype=int)
    idx = 0
    for i, res in enumerate(prot_resi):
        n_atoms = len(res.atoms)
        atom_to_res[idx:idx + n_atoms] = i
        idx += n_atoms
    return atom_to_res

def calculate_contacts_vectorized(u, prot_resi, popx_atoms, contact_matrix):

    protein_atoms = prot_resi.atoms

    # Build residue slices
    res_atom_counts = np.array([len(r.atoms) for r in prot_resi])
    res_atom_starts = np.r_[0, np.cumsum(res_atom_counts[:-1])]

    nframes = 0

    for ts in u.trajectory:
        nframes += 1

        dists = distance_array(
            protein_atoms.positions,
            popx_atoms.positions
        )

        atom_contacts = dists < cutoff

        # Collapse atom → residue
        residue_contacts = np.logical_or.reduceat(
            atom_contacts,
            res_atom_starts,
            axis=0
        )

        contact_matrix += residue_contacts.astype(int)

    contact_freq = contact_matrix / nframes
    return contact_freq

# =========================
# MOIETY COLLAPSING
# =========================
def build_moiety_map(popx_atoms):
    atom_to_moiety = []
    moiety_names = ["headgroup", "glycerol", "palmitoyl", "oleoyl"]

    for atom in popx_atoms:
        name = atom.name

        # Headgroup
        if name in {"P", "N", "O11", "O12", "O13", "O14"}:
            atom_to_moiety.append("headgroup")

        # Glycerol / linker
        elif name in {"C1", "C2", "C3", "O21", "O22", "O31", "O32"}:
            atom_to_moiety.append("glycerol")

        # sn-1 chain (palmitoyl)
        elif name.startswith("C1"):
            atom_to_moiety.append("palmitoyl")

        # sn-2 chain (oleoyl)
        elif name.startswith("C2"):
            atom_to_moiety.append("oleoyl")

        # Anything else (should now be empty)
        else:
            atom_to_moiety.append("other")

    return atom_to_moiety, moiety_names


def collapse_contacts_by_moiety(contact_freq, atom_to_moiety, moiety_names):
    n_res = contact_freq.shape[0]
    collapsed = np.zeros((n_res, len(moiety_names)))

    for j, moiety in enumerate(moiety_names):
        atom_idx = [i for i, m in enumerate(atom_to_moiety) if m == moiety]
        if atom_idx:
            # Use MAX: any atom in moiety contacts residue
            collapsed[:, j] = np.max(contact_freq[:, atom_idx], axis=1)

    return collapsed

def split_residue_index(prot_resids, subunit_length):
    prot_resids = np.asarray(prot_resids)
    subunit = (prot_resids - 1) // subunit_length + 1
    local_resi = ((prot_resids - 1) % subunit_length) + 1
    return subunit, local_resi

def collapse_across_subunits(contact_freq_moiety, local_resi):
    """
    Collapse residue contacts so residues 1–42 appear once.
    Uses max across subunits.
    """
    n_local = local_resi.max()
    n_moieties = contact_freq_moiety.shape[1]

    collapsed = np.zeros((n_local, n_moieties))

    for r in range(1, n_local + 1):
        mask = local_resi == r
        collapsed[r - 1] = np.max(contact_freq_moiety[mask], axis=0)

    return collapsed

def interacting_subunits(contact_freq_moiety, local_resi, subunit, threshold):
    """
    Returns a dict: {local_residue: sorted list of interacting subunits}
    """
    interactions = {}

    for r in np.unique(local_resi):
        mask = local_resi == r
        interacting = np.any(contact_freq_moiety[mask] > threshold, axis=1)
        interactions[r] = sorted(subunit[mask][interacting])

    return interactions

def compact_subunits(subunits):
    if not subunits:
        return ""

    ranges = []
    start = prev = subunits[0]

    for s in subunits[1:]:
        if s == prev + 1:
            prev = s
        else:
            ranges.append((start, prev))
            start = prev = s
    ranges.append((start, prev))

    formatted = []
    for a, b in ranges:
        if a == b:
            formatted.append(f"S{a}")
        else:
            formatted.append(f"S{a}–{b}")

    return ", ".join(formatted)

def build_residue_labels(prot_resi, local_resi, interacting_map):
    labels = []

    for r in range(1, int(local_resi.max()) + 1):
        idx = np.where(local_resi == r)[0][0]
        res = prot_resi[idx]

        base = f"{res.resname}{r}"
        subs = interacting_map[r]

        if subs:
            base += f" ({compact_subunits(subs)})"

        labels.append(base)


    return labels

# =========================
# PLOTTING
# =========================
def plot_moiety_contacts(res_labels, moiety_names, collapsed_contacts):

    plt.figure(figsize=(8, 8))
    sns.heatmap(
        collapsed_contacts,
        xticklabels=moiety_names,
        yticklabels=res_labels,
        cmap="magma",
        cbar_kws={"label": "Contact Frequency"},
        vmin=0, vmax=1
    )

    plt.xlabel("POPC Moiety")
    plt.ylabel("Protein Residue")
    plt.title(title)

    plt.tight_layout()
    plt.savefig(out_filename, dpi=300, bbox_inches="tight")
    plt.close()

# =========================
# MAIN
# =========================
if __name__ == "__main__":

    u, prot_resi, prot_resids, popx_atoms, contact_matrix = prep()

    contact_freq = calculate_contacts_vectorized(
        u, prot_resi, popx_atoms, contact_matrix
    )

    atom_to_moiety, moiety_names = build_moiety_map(popx_atoms)

    contact_freq_moiety = collapse_contacts_by_moiety(
        contact_freq, atom_to_moiety, moiety_names
    )

    subunit, local_resi = split_residue_index(
        prot_resids, subunit_length
    )

    collapsed_contacts = collapse_across_subunits(
        contact_freq_moiety, local_resi
    )

    interacting_map = interacting_subunits(
        contact_freq_moiety, local_resi, subunit
    )

    res_labels = build_residue_labels(
        prot_resi, local_resi, interacting_map
    )

    plot_moiety_contacts(
        res_labels, moiety_names, collapsed_contacts
    )

    print("✅ Analysis complete. Output saved to:", out_filename)