import argparse
#analysis functions
import eccentricity_avgs
import eccentricity_line
import occupancy_heatmap
import popx_popx_int
import ppi_resi
import prot_memb_avgs
import prot_memb_interactions
import prot_popx_avgs
import prot_popx_int
import prot_memb_mindist
import prot_prot_int

p = argparse.ArgumentParser(description="MD Analysis Driver")

# trajectory, structure, and system
p.add_argument("-f", "--traj", required=True, help="Trajectory file, preferrably concatenated, pbc corrected, and without waters.")
p.add_argument("-s", "--top", required=True, help="Topology file, should match the atoms in traj (no waters)")
p.add_argument("-o","--out", required=True, help="System name to append to outfile names (e.g. Hex1_lipids)")
p.add_argument("-t","--title", required=True, help="System name for title of plots, contained within quotes (i.e. 'Hexamer (Rep 1, W/ Free Lipids))")

# other required variables
p.add_argument("-l","--length",required=True, help="Chain length of protein. All subunits must be identical")
p.add_argument("-c","--cutoff", required=False, default=4.0, help="Cutoff distance (A) for interactions")

args = p.parse_args()

eccentricity_avgs(args)
eccentricity_line(args)
occupancy_heatmap(args)
popx_popx_int(args)
ppi_resi(args)
prot_memb_avgs(args)
prot_memb_interactions(args)
prot_popx_avgs(args)
prot_popx_int(args)
prot_memb_mindist(args)
prot_prot_int(args)