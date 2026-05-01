import argparse
#analysis functions
import eccentricity_line
import popx_popx_int
import ppi_resi
import prot_memb_interactions
import prot_popx_int
import prot_memb_mindist
import prot_prot_int

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

eccentricity_line(top=args.top, traj=args.traj, out=args.out, title=args.title, length=args.length, cutoff=args.cutoff)
popx_popx_int(top=args.top, traj=args.traj, out=args.out, title=args.title, length=args.length, cutoff=args.cutoff)
ppi_resi(top=args.top, traj=args.traj, out=args.out, title=args.title, length=args.length, cutoff=args.cutoff)
prot_memb_interactions(top=args.top, traj=args.traj, out=args.out, title=args.title, length=args.length, cutoff=args.cutoff)
prot_popx_int(top=args.top, traj=args.traj, out=args.out, title=args.title, length=args.length, cutoff=args.cutoff)
prot_memb_mindist(top=args.top, traj=args.traj, out=args.out, title=args.title, length=args.length, cutoff=args.cutoff)
prot_prot_int(top=args.top, traj=args.traj, out=args.out, title=args.title, length=args.length, cutoff=args.cutoff)