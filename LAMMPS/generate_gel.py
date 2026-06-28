import numpy as np
from scipy.spatial import cKDTree

# ---------------------------------------------------------------
# Step 1: Settings for this SHORT-CHAIN gel variant
# ---------------------------------------------------------------
# chain_len=2 means each "chain" is just a single bond (dumbbell),
# so there is almost NO slack -- every bond is under tension almost
# immediately when stretched. n_chains=300 keeps a reasonable total
# bead count. box_size=7.8 was found by testing -- it's the smallest
# box that still gives full connectivity while keeping the volume
# per bead safely above a bead's own physical volume (avoiding the
# severe LJ overlap that smaller boxes would cause).
np.random.seed(0)
n_chains = 300
beads_per_chain = 2
bond_length = 0.7
box_size = 7.8
crosslink_cutoff = 1.1   # safely below FENE max bond length (1.5)

# ---------------------------------------------------------------
# Step 2: Build chains as random walks
# ---------------------------------------------------------------
positions = []
bonds = []
chain_id_of_bead = []

bead_id = 0
for chain in range(n_chains):
    pos = np.random.uniform(0, box_size, size=3)
    positions.append(pos)
    chain_id_of_bead.append(chain)
    bead_id += 1

    for i in range(1, beads_per_chain):
        direction = np.random.normal(size=3)
        direction = direction / np.linalg.norm(direction)
        pos = pos + direction * bond_length
        positions.append(pos)
        chain_id_of_bead.append(chain)
        bonds.append((bead_id, bead_id + 1, 1))
        bead_id += 1

positions = np.array(positions)
positions = positions % box_size

# ---------------------------------------------------------------
# Step 3: Add ALL safe cross-chain crosslinks within cutoff
# (with only 2 beads per chain, every bead is an "end" -- no
# interior-bead spacing restriction is needed/possible here)
# ---------------------------------------------------------------
tree = cKDTree(positions, boxsize=box_size)
pairs = tree.query_pairs(r=crosslink_cutoff)

crosslinks_added = 0
for (i, j) in pairs:
    chain_i = chain_id_of_bead[i]
    chain_j = chain_id_of_bead[j]
    if chain_i != chain_j:
        bonds.append((i + 1, j + 1, 2))
        crosslinks_added += 1

n_beads = len(positions)
n_bonds = len(bonds)
print(f"Beads: {n_beads} | Bonds: {n_bonds} | Crosslinks: {crosslinks_added}")

# ---------------------------------------------------------------
# Step 3b: Verify the network is fully connected
# ---------------------------------------------------------------
parent = list(range(n_beads))
def find(x):
    while parent[x] != x:
        x = parent[x]
    return x
def union(x, y):
    px, py = find(x), find(y)
    if px != py:
        parent[px] = py

for a, b, _ in bonds:
    union(a - 1, b - 1)

cluster_sizes = {}
for i in range(n_beads):
    r = find(i)
    cluster_sizes[r] = cluster_sizes.get(r, 0) + 1

sizes = sorted(cluster_sizes.values(), reverse=True)
print(f"Number of separate clusters: {len(sizes)}")
print(f"Largest cluster size: {sizes[0]} / {n_beads} beads "
      f"({100*sizes[0]/n_beads:.1f}% of the gel)")
if len(sizes) > 1:
    print("WARNING: gel is NOT fully connected.")
else:
    print("Gel forms ONE single connected network. Good to proceed.")

# ---------------------------------------------------------------
# Step 3c: Verify INITIAL bond lengths are safe
# ---------------------------------------------------------------
box = box_size
max_len = 0.0
n_unsafe = 0
for a, b, btype in bonds:
    diff = positions[a-1] - positions[b-1]
    diff = diff - box * np.round(diff / box)
    length = np.linalg.norm(diff)
    max_len = max(max_len, length)
    if length > 1.4:
        n_unsafe += 1
print(f"Initial max bond length: {max_len:.4f} (FENE limit is 1.5)")
print(f"Bonds within 0.1 of the limit: {n_unsafe}")
if max_len > 1.4:
    print("WARNING: some bonds start dangerously close to the FENE limit!")
else:
    print("All initial bond lengths are safe.")

vol_per_bead = box_size**3 / n_beads
print(f"\nVolume per bead: {vol_per_bead:.4f} (a sigma=1.0 bead occupies ~0.52)")
if vol_per_bead < 0.6:
    print("WARNING: packing is dense -- expect a harder equilibration.")

# ---------------------------------------------------------------
# Step 4: Write the LAMMPS data file
# ---------------------------------------------------------------
with open("gel.data", "w") as f:
    f.write("LAMMPS data file -- short-chain coarse-grained gel network\n\n")
    f.write(f"{n_beads} atoms\n")
    f.write(f"{n_bonds} bonds\n\n")
    f.write("1 atom types\n")
    f.write("2 bond types\n\n")

    f.write(f"0.0 {box_size} xlo xhi\n")
    f.write(f"0.0 {box_size} ylo yhi\n")
    f.write(f"0.0 {box_size} zlo zhi\n\n")

    f.write("Masses\n\n")
    f.write("1 1.0\n\n")

    f.write("Atoms # id mol type x y z\n\n")
    for idx, pos in enumerate(positions):
        mol_id = chain_id_of_bead[idx] + 1
        f.write(f"{idx+1} {mol_id} 1 {pos[0]:.4f} {pos[1]:.4f} {pos[2]:.4f}\n")

    f.write("\nBonds # id type atom1 atom2\n\n")
    for idx, (a, b, btype) in enumerate(bonds):
        f.write(f"{idx+1} {btype} {a} {b}\n")

print("\nLAMMPS data file written as gel.data")
