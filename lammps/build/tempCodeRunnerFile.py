

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