import MDAnalysis as mda 
from MDAnalysis.analysis import distances
import numpy as np
import matplotlib.pyplot as plt

f=mda.Universe('nowat.tpr', '../traj/cat_pbc.xtc')
#Includes every membrane component I have used thusfar,
#may need to be modified in the future if addnl components are used
mem_ca = f.select_atoms('name C1 CA and resname *PG *PC POPE PSM POPS POPI CHOL GM1 GLPA CER1 GLPB') 
pro_ca = f.select_atoms('name CA and protein') #selects alpha carbons of protein residues

n_mem=len(mem_ca)
n_pro=len(pro_ca)

print(n_mem)
print(n_pro)

cutoff = 10.0 # 5 Angstrom Cutoff

# initialize contact accumulator
contact_sum = np.zeros((n_mem, n_pro))

n_frames = 0

for ts in f.trajectory:
    dist_arr = distances.distance_array(
        mem_ca.positions,
        pro_ca.positions,
        box=f.dimensions
    )
    
    contacts = (dist_arr <= cutoff).astype(int)
    contact_sum += contacts
    n_frames += 1

# average contact frequency
contact_avg = contact_sum / n_frames

contact_avg[contact_avg < 0.01] = np.nan

cutoff_display = 0.50  # e.g. at least 80% of frames

print(contact_avg)
# find rows/columns that have any interaction
mem_mask = np.any(contact_avg > cutoff_display, axis=1)
pro_mask = np.any(contact_avg > cutoff_display, axis=0)

# reduce the matrix
contact_reduced = contact_avg[mem_mask][:, pro_mask]

# also reduce labels
mem_resids_reduced = mem_ca.resnames[mem_mask]
pro_resids_reduced = pro_ca.resnames[pro_mask]

fig, ax = plt.subplots(figsize=(15,10))

im = ax.imshow(contact_reduced, cmap='viridis', vmin=0, vmax=1)

tick_interval = 1  # now you can show all since it's smaller
ax.set_yticks(np.arange(len(mem_resids_reduced)))
ax.set_xticks(np.arange(len(pro_resids_reduced)))

ax.set_yticklabels(mem_resids_reduced)
ax.set_xticklabels(pro_resids_reduced, rotation=45, ha='right')

plt.ylabel('Membrane residues')
plt.xlabel('Protein residues')
plt.title('Trajectory-averaged contacts')

cbar = fig.colorbar(im, fraction=0.025, pad=0.04)
cbar.ax.set_ylabel('Contact frequency')

plt.savefig('reduced_contact_map.png', bbox_inches='tight')
plt.close()