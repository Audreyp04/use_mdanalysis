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

# =========================
# PREP
# =========================
def prep():
    u = mda.Universe(in_top, in_traj)
    prot = u.select_atoms("protein")
    popx = u.select_atoms("resname POPX and not name H*")

    prot_resi = prot.residues
    popx_atoms = popx.atoms

    print(f"Protein residues: {len(prot_resi)}")
    print(f"POPX atoms: {len(popx_atoms)}")

    contact_matrix = np.zeros((len(prot_resi), len(popx_atoms)))
    return u, prot_resi, popx_atoms, contact_matrix


# =========================
# CONTACT CALCULATION
# =========================
def calculate_contacts(u, prot_resi, popx_atoms, contact_matrix):
    nframes = 0

    for ts in u.trajectory:
        nframes += 1
        popx_pos = popx_atoms.positions

        for i, res in enumerate(prot_resi):
            res_pos = res.atoms.positions
            dists = distance_array(res_pos, popx_pos)

            contacts = np.any(dists < cutoff, axis=0)
            contact_matrix[i] += contacts.astype(int)

    contact_freq = contact_matrix / nframes
    np.save("protein_popx_contacts.npy", contact_freq)

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


# =========================
# PLOTTING
# =========================
def plot_moiety_contacts(prot_resi, moiety_names, contact_freq):
    res_labels = [f"{res.resname}{res.resid}" for res in prot_resi]

    plt.figure(figsize=(8, 8))
    sns.heatmap(
        contact_freq,
        xticklabels=moiety_names,
        yticklabels=res_labels,
        cmap="magma",
        cbar_kws={"label": "Contact Frequency"}
    )

    plt.xlabel("POPX Moiety")
    plt.ylabel("Protein Residue")
    plt.title(title)

    plt.tight_layout()
    plt.savefig(out_filename, dpi=300, bbox_inches="tight")
    plt.close()


# =========================
# MAIN
# =========================
if __name__ == "__main__":
    u, prot_resi, popx_atoms, contact_matrix = prep()
    contact_freq = calculate_contacts(u, prot_resi, popx_atoms, contact_matrix)

    atom_to_moiety, moiety_names = build_moiety_map(popx_atoms)

    contact_freq_moiety = collapse_contacts_by_moiety(
        contact_freq, atom_to_moiety, moiety_names
    )

    plot_moiety_contacts(
        prot_resi, moiety_names, contact_freq_moiety
    )

    print("✅ Analysis complete. Output saved to:", out_filename)