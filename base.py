"""
base
"""


def indexring(atoms, ring, startatoms):
    """
    Index ring with first two atoms in "startatoms"
    """
    lr = len(ring)-1
    j0, j1 = startatoms
    idx = list(enumerate(startatoms))
    nk = 1
    while nk < lr:
        jk = [k for k in atoms[j1].neiset if k in ring and k != j0][0]
        nk += 1
        idx.append((nk, jk))
        j0, j1 = j1, jk
    return idx