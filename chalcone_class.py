"""
Categorize chalcones

TO DO:
Move chalcone polymers and other chalcones with complex structure
to this module.
"""

import itertools
from Skeleton import _Skeleton


def _ringconnections(molobj):
    """
    Get connections between any two benzyl rings. Only the chain with
    all Cs is considered here, for, e.g., identifying chalcones
    """
    chainconnections = {}

    benzylrings = molobj.benzylrings
    if len(benzylrings) < 2:
        return None

    # rings with side chains
    ringsides = {}
    for rg in benzylrings.keys():
        sides = []
        for i in rg:
            sidecs = [nei for nei in molobj.atoms[i]['neighbors']
                      if nei['index'] not in rg]
            if len(sidecs) != 1:
                continue
            nei = list(sidecs)[0]
            if nei['symbol'] == 'c':
                sides.append(((i, nei['index']),
                              nei['bondorder'],
                              nei['bondtopology']))
        if sides:
            ringsides[rg] = tuple(sides)

    if len(ringsides) < 2:
        return None

    # check connections between two benzyl rings
    for rg1, rg2 in itertools.combinations(ringsides.keys()):
        _temp = {}
        for side1, side2 in itertools.product(ringsides[rg1], ringsides[rg2]):
            ((i1, k1), bd1, bt1) = side1
            ((i2, k2), bd2, bt2) = side2
            if k1 in molobj.atoms[k2]['neisetidx']:
                nei = [nj for nj in molobj.atoms[k2]['neighbors']
                       if nj['index'] == k1][0]
                _temp['atoms'] = (i1, k1, k2, i2)
                _temp['bondorders'] = (bd1, nei['bondorder'], bd2)
                _temp['topology'] = (bt1, nei['bondtopology'], bt2)
            else:
                jx = (molobj.atoms[k1]['neisetidx']
                      & molobj.atoms[k2]['neisetidx'])
                if len(jx) == 1:
                    j = list(jx)[0]
                    if molobj.atoms[j]['symbol'] == 'c':
                        nei1 = [nj for nj in molobj.atoms[k1]['neighbors']
                                if nj['index'] == j][0]
                        nei2 = [nj for nj in molobj.atoms[k2]['neighbors']
                                if nj['index'] == j][0]
                        _temp['atoms'] = (i1, k1, j, k2, i2)
                        _temp['bondorders'] = (bd1, nei1['bondorder'],
                                               nei2['bondorder'], bd2)
                        _temp['topology'] = (bt1, nei1['bondtopology'],
                                             nei2['bondtopology'], bt2)

        if _temp:
            try:
                chainconnections[rg1][rg2] = _temp
            except KeyError:
                chainconnections[rg1] = {}
                chainconnections[rg1][rg2] = _temp
    return chainconnections


def _chalconeclass(molobj, cchains):
    """
    Identify chalcone classes
    """
    if cchains is None:
        return None, None, None

    if not all(bt == molobj.CHAIN for bt in cchains["topology"]):
        return None, None, None

    k1, k2 = cchains["atoms"][1], cchains["atoms"][-2]
    oxy1 = molobj.getoxygen(k1)
    oxy2 = molobj.getoxygen(k2)
    bd2 = None if oxy2 is None else oxy2["bondorder"]
    bd1 = None if oxy1 is None else oxy1["bondorder"]

    # two atoms between
    if (len(cchains["atoms"]) == 4
            and all(bd == 1 for bd in cchains["bondorders"])):
        if not (bd1 == 2 or bd2 == 2):
            return None, None, None

        if bd1 == 2:
            # alpha-methyldeoxybenzoin
            for nei in molobj.atoms[k2]["neighbors"]:
                j = nei["index"]
                if (nei["symbol"] == "c"
                        and nei["bondorder"] == 1
                        and len(molobj.atoms[j]["neighbors"]) == 1):
                    return ("alpha-methyldeoxybenzoin",
                            True, (oxy1["index"], j))

        if bd2 == 2:
            for nei in molobj.atoms[k1]["neighbors"]:
                j = nei["index"]
                if (nei["symbol"] == "c"
                        and nei["bondorder"] == 1
                        and len(molobj.atoms[j]["neighbors"]) == 1):
                    return ("alpha-methyldeoxybenzoin",
                            False, (oxy2["index"], j))

    # three atoms between
    elif len(cchains["atoms"]) == 5:
        # oxodihydrochalcone
        if bd1 == 2 and bd2 == 2:
            return "oxodihydrochalcone", True, (oxy1["index"], oxy2["index"])

        # chalcone
        if bd1 == 2:
            if cchains["bondorders"][2] == 2:
                return "chalcone", True, oxy1["index"]
            # dihydrochalcone
            if bd2 != 2 and all(bd == 1 for bd in cchains["bondorders"]):
                return "dihydrochalcone", True, oxy1["index"]
        if bd2 == 2:
            if cchains["bondorders"][1] == 2:
                return "chalcone", False, oxy2["index"]
            # dihydrochalcone
            if bd1 != 2 and all(bd == 1 for bd in cchains["bondorders"]):
                return "dihydrochalcone", False, oxy2["index"]

        # chalcone cyclopentenedione
        if (bd1 == 1
                and cchains["bondorders"][2] == 2
                and cchains["bondorders"][0] == 2):
            return "chalcone cyclopentenedione", True, oxy1["index"]

        if (bd2 == 1
                and cchains["bondorders"][1] == 2
                and cchains["bondorders"][3] == 2):
            return "chalcone cyclopentenedione", False, oxy2["index"]

    return None, None, None


def _chalClass(molobj, connections):
    """
    check whether it is a chalcone or its derivatives, including
    chalcone, dihydrochalcone, oxodihydrochalcone, chalcone
    cyclopentenedione, chalcone-quinol, alpha-methyldeoxybenzoin
    and so on.
    """
    skeletons = []
    if not connections:
        return skeletons

    # assign names to chalcones
    for rg1 in connections.keys():
        for rg2 in connections[rg1].keys():
            connectchain = connections[rg1][rg2]
            n = len(connectchain['atoms'])
            skj = _Skeleton()
            name, ordertag, sideatoms = _chalconeclass(molobj, connectchain)

            if name is None:
                continue

            # match the name to the benzyl ring assigned, for name correction
            rga = rg1 if ordertag else rg2
            rgb = rg2 if ordertag else rg1
            if molobj.benzylrings[rgb]['label'] != 'b':
                continue

            if molobj.benzylrings[rga]['label'] == 'b':
                pass
            elif molobj.benzylrings[rga]['label'].startswith('q'):
                name = 'chalcone-quinol' if name == 'chalcone' else None
            elif name == 'chalcone cyclopentenedione':
                if not molobj.benzylrings[rga]['label'] == 'cc':
                    name = None
            if name is None:
                continue

            skj.ringA = rga
            skj.ringB = rgb
            skj.connects['atomc'] = connectchain['atoms'][0]
            skj.connects['atomb'] = connectchain['atoms'][-1]
            skj.connects['connectatoms'] = connectchain[1:n-1]
            skj.connects['topology'] = connectchain['bondtopology']
            skj.corders = connectchain['bondorders']
            if not ordertag:
                skj.connects['atomc'] = connectchain['atoms'][-1]
                skj.connects['atomb'] = connectchain['atoms'][0]
                skj.connects['connectatoms'] = connectchain[n-1:0:-1]
                skj.connects['topology'] = connectchain['bondtopology'][::-1]
                skj.corders = connectchain['bondorders'][::-1]
            skj.label = 'c'
            skj.name = name
            skj.sideatoms = sideatoms

            skeletons.append(skj)

    return skeletons    


def Chalcone(molobj):
    """
    Chalcone class
    """
    connections = _ringconnections(molobj)
    return _chalClass(molobj, connections)