import argparse
import MDAnalysis as mda
#analysis functions
import eccentricity_line
import popx_popx_int
import ppi_resi
import prot_memb_interactions
import prot_popx_int
import prot_memb_mindist
import prot_prot_int

# EVALUATE IF POPX ANALYSIS IS NEEDED

def system_contains_popx(top, traj):
    """
    Returns True if POPX residues exist in the system.
    """
    u = mda.Universe(top, traj)
    return len(u.select_atoms("resname POPX")) > 0


p = argparse.ArgumentParser(description="MD Analysis Driver")

# trajectory, structure, and system
p.add_argument("-f", "--traj", required=True, help="Trajectory file (path), preferrably concatenated, pbc corrected, and without waters.")
p.add_argument("-s", "--top", required=True, help="Topology file (path), should match the atoms in traj (no waters)")
p.add_argument("-o","--out", required=True, help="System name to append to outfile names (e.g. Hex1_lipids)")
p.add_argument("-t","--title", required=True, help="System name for title of plots, contained within quotes (i.e. 'Hexamer (Rep 1, W/ Free Lipids))")

# other required variables
p.add_argument("-l","--length",required=True, help="Chain length of protein. All subunits must be identical", type=int)
p.add_argument("-c","--cutoff", required=False, default=4.0, help="Cutoff distance (A) for interactions", type=float)

args = p.parse_args()


has_popx = system_contains_popx(args.top, args.traj)

if not has_popx:
    print("⚠️  POPX not detected — skipping POPX analyses.")

# =========================
# ECCENTRICITY
# =========================
eccentricity_line.init(
    top=args.top,
    traj=args.traj,
    out=args.out,
    title=args.title,
)

eccentricity_line.align_trajectory()
eccentricity_line.calculate_eccentricity()
eccentricity_line.plot_eccentricity()


# =========================
# POPX–POPX INTERACTIONS
# =========================
if has_popx:
    popx_popx_int.init(
        top=args.top,
        traj=args.traj,
        out=args.out,
        title=args.title,
        cutoff_val=args.cutoff,
    )

    atom_to_lipid = popx_popx_int.build_atom_to_lipid_map()
    contact_freq = popx_popx_int.calculate_popx_popx_contacts(atom_to_lipid)

    atom_to_moiety, moiety_names = popx_popx_int.build_moiety_map()
    collapsed = popx_popx_int.collapse_contacts(
        contact_freq, atom_to_moiety, moiety_names
    )

    popx_popx_int.plot_moiety_contacts(moiety_names, collapsed)


# =========================
# RESIDUE–RESIDUE PPI
# =========================
ppi_resi.init(
    top=args.top,
    traj=args.traj,
    out=args.out,
    title=args.title,
    cutoff_val=args.cutoff,
    length=args.length,
)

u, chains, residue_labels = ppi_resi.prep()
contact_freq = ppi_resi.calculate_residue_contacts(u, chains)
ppi_resi.plot_moiety_contacts(contact_freq, residue_labels)


# =========================
# PROTEIN–MEMBRANE
# =========================
prot_memb_interactions.init(
    top=args.top,
    traj=args.traj,
    out=args.out,
    title=args.title,
    cutoff_val=args.cutoff,
    length=args.length,
)

u, prot_resi, prot_resids, memb_resi, memb_atoms, cm = (
    prot_memb_interactions.prep()
)

contact_freq, nframes = prot_memb_interactions.calculate_contacts_capped(
    u, prot_resi, memb_atoms, cm
)

subunit, local_resi = prot_memb_interactions.split_residue_index(
    prot_resids, args.length
)

collapsed = prot_memb_interactions.collapse_across_subunits(
    contact_freq, local_resi
)

collapsed = prot_memb_interactions.collapse_membrane_atoms_to_residues(
    collapsed, memb_atoms, memb_resi
)

collapsed, lipid_types = prot_memb_interactions.average_contacts_by_lipid_type(
    collapsed, memb_resi
)

interacting_map = prot_memb_interactions.interacting_subunits(
    contact_freq, local_resi, subunit, threshold=0.80
)

res_labels = prot_memb_interactions.build_residue_labels(
    prot_resi, local_resi, interacting_map
)

prot_memb_interactions.plot_contacts(
    res_labels, lipid_types, collapsed, nframes
)


# =========================
# PROTEIN–POPX
# =========================
if has_popx:
    prot_popx_int.init(
        top=args.top,
        traj=args.traj,
        out=args.out,
        title=args.title,
        cutoff_val=args.cutoff,
        length=args.length,
    )

    u, prot_resi, prot_resids, popx_atoms, cm = prot_popx_int.prep()

    contact_freq = prot_popx_int.calculate_contacts_vectorized(
        u, prot_resi, popx_atoms, cm
    )

    atom_to_moiety, moiety_names = prot_popx_int.build_moiety_map(popx_atoms)
    contact_freq_moiety = prot_popx_int.collapse_contacts_by_moiety(
        contact_freq, atom_to_moiety, moiety_names
    )

    subunit, local_resi = prot_popx_int.split_residue_index(
        prot_resids, args.length
    )

    collapsed = prot_popx_int.collapse_across_subunits(
        contact_freq_moiety, local_resi
    )

    interacting_map = prot_popx_int.interacting_subunits(
        contact_freq_moiety, local_resi, subunit
    )

    res_labels = prot_popx_int.build_residue_labels(
        prot_resi, local_resi, interacting_map
    )

    prot_popx_int.plot_moiety_contacts(
        res_labels, moiety_names, collapsed
    )

# =========================
# MIDPLANE DISTANCE
# =========================
prot_memb_mindist.init(
    top=args.top,
    traj=args.traj,
    out=args.out,
    title=args.title,
)

u, prot, memb = prot_memb_mindist.prep()
time, dist = prot_memb_mindist.calculate_distances(u, prot, memb)
prot_memb_mindist.plot_moiety_contacts(time, dist)


# =========================
# CHAIN–CHAIN PPI
# =========================
prot_prot_int.init(
    top=args.top,
    traj=args.traj,
    out=args.out,
    title=args.title,
    cutoff_val=args.cutoff,
    length=args.length,
)

u, chains = prot_prot_int.prep()
chain_ids, contact_freq = prot_prot_int.calculate_chain_contacts(u, chains)
prot_prot_int.plot_moiety_contacts(chain_ids, contact_freq)


print("✅ All analyses complete.")