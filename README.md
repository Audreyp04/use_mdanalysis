# use_mdanalysis
This analysis set runs: 
eccentricity, eccentric occupancy, free popc interactions (**IF resname = POPX ONLY), protein-protein interactions (chain and residue), protein-membrane interactions by membrane component and protein residue, distance of protein-membrane midplane

Entire per replicate analysis set can be run using the run_analysis.py script

USAGE:

-f, --traj = Trajectory file (path), preferrably concatenated, pbc corrected, and without waters.

-s, --top = Topology file (path), should match the atoms in traj (no waters)

-o,--out = System name to append to outfile names (e.g. Hex1_lipids)

-t,--title = System name for title of plots, contained within quotes (i.e. 'Hexamer (Rep 1, W/ Free Lipids)')

-l,--length = Chain length of protein. All subunits must be identical

-c,--cutoff = Cutoff distance (A) for interactions (optional, default = 4.0)

average analyses must be run manually, modifying user inputs in the corresponding .py script

If there is additional analysis you would like added, submit an issue labeled "enhancement"

If there is an issue you encounter in any scripts, submit an issue labeled "bug"

If you uncover a problem with script logic/calculations/etc, submit an issue labeled "invalid"