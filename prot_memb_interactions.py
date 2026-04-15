

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
title = "Hexamer & Membrane Lipid Interactions (Rep 1)"
out_filename = "hex_memb_moiety_contacts.png"
cutoff = 4.0          # Å
subunit_length = 42
threshold=0.80

# =========================
# PREP
# =========================
def prep():
    u = mda.Universe(in_top, in_traj)
    prot = u.select_atoms("protein")
    memb = u.select_atoms("not resname POPX and not name H* and not protein and not resname CL and not resname K")

    prot_resi = prot.residues
    prot_resids = prot.residues.resids 
    memb_resi = memb.residues
    memb_atoms = memb.atoms

    contact_matrix = np.zeros((len(prot_resi), len(memb_atoms)))
    return u, prot_resi, prot_resids, memb_resi, memb_atoms, contact_matrix


# =========================
# CONTACT CALCULATION
# =========================
def collapse_membrane_atoms_to_residues(contact_freq, memb_atoms, memb_resi):
    """
    Collapse membrane atom contacts into membrane residue contacts.
    """
    n_prot = contact_freq.shape[0]
    n_memb_res = len(memb_resi)

    collapsed = np.zeros((n_prot, n_memb_res))

    # atom → residue index for membrane atoms
    atom_to_res = memb_atoms.resindices  # 0..n_memb_res-1

    for atom_idx, res_idx in enumerate(atom_to_res):
        collapsed[:, res_idx] = np.maximum(
            collapsed[:, res_idx],
            contact_freq[:, atom_idx]
        )

    return collapsed

def calculate_contacts_vectorized(u, prot_resi, memb_atoms, contact_matrix):

    protein_atoms = prot_resi.atoms

    # Build residue slices
    res_atom_counts = np.array([len(r.atoms) for r in prot_resi])
    res_atom_starts = np.r_[0, np.cumsum(res_atom_counts[:-1])]

    nframes = 0

    for ts in u.trajectory:
        nframes += 1

        dists = distance_array(
            protein_atoms.positions,
            memb_atoms.positions
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

def average_contacts_by_lipid_type(collapsed_contacts, memb_resi):
    """
    Average membrane contacts by lipid type (resname).
    """
    # Get lipid type per membrane residue
    lipid_names = np.array([r.resname for r in memb_resi])
    lipid_types = np.unique(lipid_names)

    avg_contacts = np.zeros((collapsed_contacts.shape[0], len(lipid_types)))

    for i, lipid in enumerate(lipid_types):
        mask = lipid_names == lipid
        avg_contacts[:, i] = collapsed_contacts[:, mask].mean(axis=1)

    return avg_contacts, lipid_types


def split_residue_index(prot_resids, subunit_length):
    prot_resids = np.asarray(prot_resids)
    subunit = (prot_resids - 1) // subunit_length + 1
    local_resi = ((prot_resids - 1) % subunit_length) + 1
    return subunit, local_resi

def collapse_across_subunits(contact_freq, local_resi):
    """
    Collapse residue contacts so residues appear once.
    Uses max across subunits.
    """
    n_local = local_resi.max()
    n_contacts = contact_freq.shape[1]

    collapsed_contacts = np.zeros((n_local, n_contacts))

    for r in range(1, n_local + 1):
        mask = local_resi == r
        collapsed_contacts[r - 1] = np.max(contact_freq[mask], axis=0)

    return collapsed_contacts

def interacting_subunits(contact_freq, local_resi, subunit, threshold):
    """
    Returns a dict: {local_residue: sorted list of interacting subunits}
    """
    interactions = {}

    for r in np.unique(local_resi):
        mask = local_resi == r
        interacting = np.any(contact_freq[mask] > threshold, axis=1)
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
def plot_contacts(res_labels, lipid_types, collapsed_contacts):

    plt.figure(figsize=(8, 8))
    sns.heatmap(
        collapsed_contacts,
        xticklabels=lipid_types,
        yticklabels=res_labels,
        cmap="magma",
        cbar_kws={"label": "Mean Contact Frequency"},
        vmin=0, vmax=1
    )

    plt.xlabel("Membrane Component")
    plt.ylabel("Protein Residue")
    plt.title(title)

    plt.tight_layout()
    plt.savefig(out_filename, dpi=300, bbox_inches="tight")
    plt.close()

# =========================
# MAIN
# =========================
if __name__ == "__main__":

    u, prot_resi, prot_resids, memb_resi, memb_atoms, contact_matrix = prep()

    contact_freq = calculate_contacts_vectorized(
        u, prot_resi, memb_atoms, contact_matrix
    )

    subunit, local_resi = split_residue_index(
        prot_resids, subunit_length
    )

    collapsed_contacts = collapse_across_subunits(
        contact_freq, local_resi
    )
    
    collapsed_contacts = collapse_membrane_atoms_to_residues(
    collapsed_contacts,
    memb_atoms,
    memb_resi
    )

    
    collapsed_contacts, lipid_types = average_contacts_by_lipid_type(
        collapsed_contacts,
        memb_resi
    )


    interacting_map = interacting_subunits(
        contact_freq, local_resi, subunit, threshold
    )

    res_labels = build_residue_labels(
        prot_resi, local_resi, interacting_map
    )

    plot_contacts(
        res_labels, lipid_types, collapsed_contacts
    )

    print("✅ Analysis complete. Output saved to:", out_filename)