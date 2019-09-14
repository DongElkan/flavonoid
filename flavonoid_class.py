"""
Flavonoid classes

TO DO:
Move flavonoids with complex structures to this module
"""

import itertools

from Skeleton import _Skeleton
from base import indexring
from flavonoid_class_check import (_flavoneclass,
                                   _flavanclass_chain,
                                   _flavanclass_ring,
                                   _stilbenoflavonoid
                                   )


def _getsk(molobj, ringAC, ringB):
    """
    Get the skeleton of the flavonoid, with the atom in ring C that
    connects to ring B (for flavonoids and its derivatives), or the
    atom C between ring C and ring B (for homoflavonoids).
    """
    flavsk = []
    for rgix, j, c, t in ringB:
        for rgixi, rgixk, _ in ringAC:
            skj = _Skeleton()
            # ring B should not have shared atoms with ring A and C
            if set(rgix) & set(rgixi) or set(rgix) & set(rgixk):
                continue
            skj.ringA = rgixi
            skj.ringC = rgixk
            skj.ringB = rgix
            if c in rgixk:
                skj.connects = {"atomc": c, "atomb": j, "topology": t}
                # this is flavonoid derivatives
                skj.label = "f"
                flavsk.append(skj)
            else:
                for nei in molobj.atoms[c].neighbors:
                    if nei.index in rgixk and nei.symbol == "c":
                        skj.connects = {"atomc": (c, nei.index),
                                        "atomb": j,
                                        "topology": (t, nei.bondtopology)}
                        # this is homoflavonoid derivatives
                        skj.label = "hf"
                        flavsk.append(skj)

    return tuple(flavsk)


def _getringABC(molobj):
    """
    Ring A, B and C for the flavonoid
    """
    ringac, ringb = [], []
    for rgix, rgobj in molobj.benzylrings.items():
        # potential ring A and C, which connected together
        if rgobj.isringAC is None:
            continue

        # potential ring B, which has a branch with the bond
        # destination being atom symbol C
        for j in rgix:
            ringb += [(rgix, j, nei.index, nei.bondtopology)
                      for nei in molobj.atoms[j].neighbors
                      if nei.symbol == "c" and nei.index not in rgix]

        for ck, (rgixk, rgobjk) in itertools.product(rgobj.ringACnns,
                                                     molobj.thps.items()):
            if rgobjk.index.issuperset(ck):
                ringac.append((rgix, rgixk, ck))

    return ringac, ringb


def _indexringAC(molobj, ringac):
    """
    Index ring A and C
    """
    rac_index = {}
    for rgixi, rgixk, connects in ringac:
        ack = {}

        # reindex ring C with starting atom at oxygen
        c9, c10 = connects[2], connects[0]
        oxy = connects[3]
        nextc = [k for k in molobj.atoms[oxy].neiset
                 if k != c9 and k in rgixk][0]
        reindicesC = indexring(molobj, rgixk, (oxy, nextc))

        # reindex ring A with starting atom at C(5), next to C(10)
        c5 = [k for k in molobj.atoms[c10].neighbors
              if k != c9 and k in rgixi][0]
        reindicesA = indexring(molobj, rgixi, (c10, c5))

        # reindex ring A and C together
        nj = 0
        for _, k in reindicesC[:4]:
            nj += 1
            ack[k] = nj
        for _, k in reindicesA[1:-1]:
            nj += 1
            ack[k] = nj

        # C(9) and C(10)
        ack[c9] = nj+1
        ack[c10] = nj+2
        try:
            rac_index[rgixi].append(ack)
        except KeyError:
            rac_index[rgixi] = [ack]

    return rac_index


def _flavClass(molobj, pskeletons, ringACindex):
    """
    Check whether it is the skeleton of a C6-C3-C6 structured flavonoid,
    flavonolignan or homoisoflavonoids, as well as other flavonoid
    derivatives, such as dehydropeltogynoid, rotenoid, pterocarpan
    etc.
    """
    # is no skeleton is identified, return None
    if len(pskeletons) == 0:
        return None

    # get skeleton
    for skj in pskeletons:
        name, params = None, {}
        # set up parameters for defining flavonoid classes
        for key in ["ringA", "ringB", "ringC"]:
            params[key] = getattr(skj, key)

        # connections between ring B and C
        if skj.label == "f":
            cj = skj.connects["atomc"]
        else:
            cj = skj.connects["atomc"][0]

        # return atoms around in ring C possibly connecting to ring B
        catoms = {}
        for ack in ringACindex[params["ringA"]]:
            try:
                # atom order in ring C connects to ring B
                j = ack[cj]
            except KeyError:
                continue
            # get o(1), c(3) and c(4) in ring C
            for ck, rk in ack.items():
                catoms[rk] = ck
        if not catoms:
            continue

        # parameter values
        #     index of atom in C connecting to ring B
        params["k"] = j
        #     atoms between ring B and C, including the atom in ring C
        params["atomsbc"] = skj.connects["atomc"]
        #     topology of bond between ring B and C
        params["ischain"] = False
        if (isinstance(skj.connects["atomc"], int)
                and skj.connects["topology"] == molobj.RING):
            params["ischain"] = True
        elif (isinstance(skj.connects["atomc"], tuple)
                and all(bd == molobj.CHAIN
                        for bd in skj.connects["topology"])):
            params["ischain"] = True

        #     the oxygen atoms around ring C
        for i in range(2, 5):
            bdoxyi = molobj.getoxygen(catoms[i])
            if bdoxyi is not None:
                params['o%d' % i], params['bdoc%d' % i] = bdoxyi[:2]
            else:
                params['o%d' % i], params['bdoc%d' % i] = None, None
        #     bond order of C(2) and C(3)
        params['bdc23'] = molobj.getbondorder(catoms[2], catoms[3])
        #     bond order of C(3) and C(4)
        params['bdc34'] = molobj.getbondorder(catoms[3], catoms[4])

        # assign names to the skeleton
        if params['bdc23'] == 2:
            name, sdatoms, sdrings = _flavoneclass(molobj, params)
        elif params['bdc23'] == 1:
            if params['bondbc'] == molobj.CHAIN:
                name, sdatoms, sdrings = _flavanclass_chain(molobj, params)
            else:
                name, sdatoms, sdrings = _flavanclass_ring(molobj, params)

        if name is None:
            skj.name = None
            continue

        # put into structure
        skj.name = name
        skj.siderings = sdrings
        skj.sdatoms = sdatoms
        skj.corders = catoms
        skj.derivatives = None

        if name == 'flavan-3-ol':
            derivative, sidegroups = _stilbenoflavonoid(molobj, skj)
            if derivative is not None:
                skj.name = derivative
                skj.derivatives = sidegroups

    return [skj for skj in pskeletons if skj.name is not None]


def _tetrahydroflav(molobj, skeleton, ringACindex):
    """
    Tetrahydroflavanones
    """
    # pseudo ring B
    for rgix, rgobj in molobj.benzylrings.items():
        if (rgobj.label is None
                and all(molobj.atoms[i].symbol == "c" for i in rgix)
                and len(rgix) == 6):
            if rgobj.smiles.count("=") == 1:
                rgobj = rgobj._replace(label="bpseudo")

    # get skeleton
    rgA = skeleton.ringA
    skix = ringACindex[rgA][0]
    if (molobj.benzylrings[rgA].label == "bpseudo"
            and skeleton.label == "f"
            and skix[skeleton.connects['atomc']] == 2):
        skeleton.name = 'tetrahydroflavanones'

    return skeleton


def Flavonoid(molobj):
    """
    Flavonoid class
    """
    # set up ring A, B, C and indices
    ringAC, ringB = _getringABC(molobj)
    skeleton = _getsk(molobj, ringAC, ringB)
    ringACindex = _indexringAC(molobj, ringAC)

    # check whether it is tetrahydroflavanones
    if len(skeleton) == 1:
        skeleton = _tetrahydroflav(molobj, skeleton, ringACindex)

        if skeleton.name is not None:
            return [skeleton]

    # get the skeletons with names
    skeleton = _flavClass(molobj, skeleton, ringACindex)

    return skeleton
