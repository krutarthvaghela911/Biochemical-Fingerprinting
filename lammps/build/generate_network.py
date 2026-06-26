import numpy as np
import random

# ─── Network parameters ───────────────────────────────────────────────────────
box_size   = 20.0
chain_len  = 8
n_chains   = 60
bond_len   = 0.97

random.seed(42)
np.random.seed(42)

atoms    = []
bonds    = []
atom_id  = 1
bond_id  = 1

chain_ends = []

# ─── Generate chains ──────────────────────────────────────────────────────────
for c in range(n_chains):
    pos = np.random.uniform(2.0, box_size - 2.0, size=3)
    chain_atom_ids = []

    for b in range(chain_len):
        atoms.append([atom_id, 1, pos[0], pos[1], pos[2]])
        chain_atom_ids.append(atom_id)

        if b > 0:
            bonds.append([bond_id, 1, atom_id - 1, atom_id])
            bond_id += 1

        if b < chain_len - 1:
            step = np.random.randn(3)
            step = step / np.linalg.norm(step) * bond_len
            pos  = np.clip(pos + step, 0.5, box_size - 0.5)

        atom_id += 1

    chain_ends.append((chain_atom_ids[0], chain_atom_ids[-1]))

# ─── Force crosslinks by moving end beads closer ─────────────────────────────
crosslink_count = 0
used = set()
crosslink_cutoff = 1.3

for i in range(len(chain_ends)):
    for j in range(i + 1, len(chain_ends)):
        if crosslink_count >= 40:
            break
        for end_i in chain_ends[i]:
            for end_j in chain_ends[j]:
                if end_i in used or end_j in used:
                    continue
                pi = np.array(atoms[end_i - 1][2:5])
                pj = np.array(atoms[end_j - 1][2:5])
                dist = np.linalg.norm(pi - pj)
                if dist < crosslink_cutoff:
                    bonds.append([bond_id, 2, end_i, end_j])
                    bond_id += 1
                    crosslink_count += 1
                    used.add(end_i)
                    used.add(end_j)

# ─── If still too few, force pair closest ends ───────────────────────────────
if crosslink_count < 20:
    end_list = []
    for i, (s, e) in enumerate(chain_ends):
        if s not in used:
            end_list.append((i, s, np.array(atoms[s-1][2:5])))
        if e not in used:
            end_list.append((i, e, np.array(atoms[e-1][2:5])))

    while len(end_list) >= 2 and crosslink_count < 30:
        best_dist = 999
        best_pair = None
        for a in range(len(end_list)):
            for b in range(a+1, len(end_list)):
                if end_list[a][0] == end_list[b][0]:
                    continue
                d = np.linalg.norm(end_list[a][2] - end_list[b][2])
                if d < best_dist:
                    best_dist = d
                    best_pair = (a, b)
        if best_pair is None:
            break
        a, b = best_pair
        id_a = end_list[a][1]
        id_b = end_list[b][1]

        # move them closer to each other
        pa = np.array(atoms[id_a-1][2:5])
        pb = np.array(atoms[id_b-1][2:5])
        mid = (pa + pb) / 2
        direction = (pb - pa)
        if np.linalg.norm(direction) > 0:
            direction = direction / np.linalg.norm(direction)
        atoms[id_a-1][2:5] = list(mid - direction * 0.45)
        atoms[id_b-1][2:5] = list(mid + direction * 0.45)

        bonds.append([bond_id, 2, id_a, id_b])
        bond_id += 1
        crosslink_count += 1
        used.add(id_a)
        used.add(id_b)

        end_list = [e for e in end_list if e[1] not in used]

print(f"Atoms:      {len(atoms)}")
print(f"Bonds:      {len(bonds)}")
print(f"Crosslinks: {crosslink_count}")

# ─── Write LAMMPS data file ───────────────────────────────────────────────────
out = "/home/krutarth0911/biomechanical-fingerprinting/lammps/build/gel.data"

with open(out, "w") as f:
    f.write("LAMMPS polymer gel network\n\n")
    f.write(f"{len(atoms)} atoms\n")
    f.write(f"{len(bonds)} bonds\n\n")
    f.write("2 atom types\n")
    f.write("2 bond types\n\n")
    f.write(f"0.0 {box_size} xlo xhi\n")
    f.write(f"0.0 {box_size} ylo yhi\n")
    f.write(f"0.0 {box_size} zlo zhi\n\n")
    f.write("Masses\n\n")
    f.write("1 1.0\n")
    f.write("2 1.0\n\n")
    f.write("Atoms\n\n")
    for a in atoms:
        f.write(f"{a[0]} 1 {a[1]} {a[2]:.4f} {a[3]:.4f} {a[4]:.4f}\n")
    f.write("\nBonds\n\n")
    for b in bonds:
        f.write(f"{b[0]} {b[1]} {b[2]} {b[3]}\n")

print(f"Written to: {out}")
