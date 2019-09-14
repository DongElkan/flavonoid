"""
Define flavonoid classes according to the structure parameters
provided by the module "flavonoid_categorizer"
"""


def _checklignan(molobj, ringBobj):
    """
    check whether it is a lignan
    """

    for rgk, rgobjk in molobj.rings.items():
        ckk = rgobjk.index & ringBobj.index
        if not ckk:
            continue
        kkoxys = [molobj.getoxygen(ck) for ck in ckk]
        if (rgobjk.smiles.count("O") == 2
                and len(ckk) == 2 and all(kkoxys)
                and all(ck[0] in rgk for ck in kkoxys)):
            return "flavonolignan", rgk

    return None, None


def _flavoneclass(molobj, params):
    """
    flavone class: in this class, bonder order between C(2) and C(3) is 2
    """
    # ring objects
    rgbobj = molobj.benzylrings[params["ringB"]]
    rgcobj = molobj.thps[params["ringC"]]
    # in this kind of classification, only the order of the bond between
    # C(2) and C(3) equals to 2 is considered
    if (params["bdc23"] != 2
            or (isinstance(params["atomsbc"], tuple)
                and len(params["atomsbc"]) > 2)
            or rgbobj.index & rgcobj.index):
        return None, None, None

    # assign name to the skeleton
    if (len(params["ringB"]) == 5
            and isinstance(params["atomsbc"], int)
            and params["ischain"]):
        # "2-arylbenzofuran"
        return "2-arylbenzofuran", None, None

    if isinstance(params["atomsbc"], int):

        if not params["ischain"]:
            return None, None, None

        if params["k"] == 2:
            # flavone
            if params["bdoc4"] == 2 and params["bdoc3"] == 1:
                return "flavonol", (params["o3"], params["o4"]), None

            # flavonol
            if params["bdoc4"] == 2:
                return "flavone", params["o4"], None

        # isoflavone and dehydrorotenoid
        if params["k"] == 3:
            return "isoflavone", params["o4"], None

    else:

        # homoisoflavone
        if (params["k"] == 3 and params["ischain"]):
            return ("homoisoflavone",
                    (params["o4"], params["atomsbc"][1]),
                    None)

    # coumaronochromone
    if isinstance(params["atomsbc"], int) and params["ischain"]:
        for rg in molobj.rings:
            if not (len(rg.index & rgbobj.index) == 2
                    and len(rg.index & rgcobj.index) == 2):
                continue

            if (rgcobj["ringAc"]
                    and rg.index.issuperset(rgcobj["ringAcnns"])):
                if len(rg) == 5 and params["bdoc4"] == 2:
                    return "coumaronochromone", params["o4"], rg
                elif len(rg) == 6:
                    if params["k"] == 2:
                        return "dehydrorotenoid", None, rg
                    elif params["k"] == 3:
                        return "flavone C3-C2' Ether Linkage (R6)", None, rg
                elif len(rg) == 7:
                    return "flavone C3-C2' Ether Linkage (R7)", None, rg
            elif params["k"] == 2 and params["bdoc3"] == 1:
                return "dehydropeltogynoid", None, rg

    return None, None, None


def _flavanclass_chain(molobj, params):
    """
    Flavan class, in which order of C(3)-C(2) bond is 1
    For chain, the topology of the bond connecting ring B and C
    should be chain, that is, no ring exists between the two rings.
    """
    if not params["ischain"]:
        return None, None, None

    # ring objects
    rgbobj = molobj.benzylrings[params["ringB"]]
    rgcobj = molobj.thps[params["ringC"]]

    # homo-flavonoids
    if (params["k"] == 3
            and len(params["atomsbc"]) == 2
            and params["ischain"]):
        bdk = molobj.getbondorder(*params["atomsbc"])
        if params["o4"] is None:
            if bdk != 1:
                return None, None, None

            if params["bdc34"] == 1:
                return "homoisoflavan", params["atomsbc"][1], None
            elif params["bdc34"] == 2:
                return "homoisoflavene", params["atomsbc"][1], None

        if params["o4"] == 2 and params["bdc34"] == 1:
            if bdk == 1:
                return ("homoisoflavanone",
                        (params["o4"], params["atomsbc"][1]), None)
            elif bdk == 2:
                return ("3-benzylidene-4-chromanone",
                        (params["o4"], params["atomsbc"][1]), None)

    if isinstance(params["atomsbc"], tuple):
        return None, None, None

    # in this kind of classification, only the order of the bond between
    # C(2) and C(3) equals to 2 is considered
    if len(params["ringB"]) == 5 and isinstance(params["atomsbc"], int):
        # "2-arylbenzofuran"
        return "2-arylbenzofuran", None, None

    # aurone
    if len(rgcobj.index) == 5 and params["bdoc3"] == 2:
        bcbd = molobj.getbondorder(*params["atomsbc"])
        if bcbd == 2:
            return "aurone", params["o3"], None

        bco = molobj.getoxygen(params["atomsbc"][0])
        if bcbd == 1 and bco is not None:
            return "dihydroaurone", params["o3"], None

    # isoflavene, 3-arylcoumarin, neoflavonoid
    if params["bdc34"] == 2 and isinstance(params["atomsbc"], int):
        o1 = molobj.benzylrings[params["ringA"]]["ringAcnns"][3]
        if molobj.atoms[o1]["charge"] == 0:

            if params["bdoc2"] == 2:
                if params["k"] == 3:
                    return "3-arylcoumarin", params["o2"], None
                elif params["k"] == 4:
                    return "neoflavonoid", params["o2"], None
            return "isoflavene", None, None
        else:
            # anthocyanidin
            c2 = [nei["index"] for nei in molobj.atoms[o1]["neighbors"]
                  if nei["index"] not in params["ringA"]]
            if molobj.getbondorder(o1, c2) == 2:
                return "anthocyanidin", None, None

    # assign names to the skeleton
    if params["bdoc3"] is None:

        if params["bdoc4"] == 2:

            # flavanone and flavonolignan
            if params["k"] == 2:
                name, rk = _checklignan(molobj, rgbobj)
                if name is not None:
                    return "flavonolignan", params["o4"], rk
                return "flavanone", params["o4"], None

            # isoflavanone
            if params["k"] == 3:
                name, rk = _checklignan(molobj, rgbobj)
                if name is not None:
                    return "flavonolignan", params["o4"], rk
                return "isoflavanone", params["o4"], None

        if params["bdoc4"] == 1:

            # flavan-4-ol
            if params["k"] == 2:
                return "flavan-4-ol", params["o4"], None

            # isoflavanol
            if params["k"] == 3:
                return "isoflavanol", params["o4"], None

        if params["bdoc4"] is None:

            # flavan
            if params["k"] == 2:
                return "flavan", None, None

            if params["k"] == 3:
                # isoflavanquinone
                if rgbobj.label == "q":
                    return "isoflavanquinone", None, None
                # "isoflavan"
                return "isoflavan", None, None

    elif params["bdoc3"] == 1:

        # flavan-3,4-diol
        if params["k"] == 2 and params["bdoc4"] == 1:
            return "flavan-3,4-diol", (params["o4"], params["o3"]), None

        # flavan-3-ol
        if params["k"] == 2:
            return "flavan-3-ol", params["o3"], None

        # dihydroflavonol
        if params["k"] == 2 and params["bdoc4"] == 2:
            name, rk = _checklignan(molobj, rgbobj)
            if name is not None:
                return "flavonolignan", params["o4"], rk
            return "dihydroflavonol", (params["o4"], params["o3"]), None

    return None, None, None


def _flavanclass_ring(molobj, params):
    """
    Flavan class, in which order of C(3)-C(2) bond is 1
    For ring, the topology of the bond connecting ring B and C
    should be ring, that is, there is a ring between the two rings.
    """
    if params["ischain"] or isinstance(params["atomsbc"], tuple):
        return None, None, None

    # ring objects
    rgbobj = molobj.benzylrings[params["ringb"]]
    rgcobj = molobj.thps[params["ringc"]]

    candidaterings = dict(molobj.benzylrings, **molobj.thps)
    rgds = []
    for rgk, rgkobj in candidaterings.keys():
        if (len(rgkobj.index & rgbobj.index) == 2
                and len(rgkobj.index & rgcobj.index) == 2):
            rgds.append(rgk)

    # if there are more than 1 rings between ring B and C,
    # reject it
    if len(rgds) > 1:
        return None, None, None

    # find the ring between ring B and C
    # .. 4-arylchroman
    rgd = rgds[0]
    if rgbobj.ringAC is None or "O" not in candidaterings[rgd].smiles:
        if (len(rgcobj.index) == 5
                and params["bdoc3"] == 1
                and params["k"] == 4
                and all(molobj.atoms[i].symbol == "c" for i in rgd)):
            return "4-arylchroman", params["o3"], rgd

        if (len(rgcobj.index) == 6
                and params["bdoc3"] == 1
                and params["o3"] in rgd and params["k"] == 2):
            return "peltogynoid", None, rgd

    else:
        if not any(candidaterings[rgd].index.issuperset(atm)
                   for atm in rgbobj["ringAcnns"]):
            return None, None, None

        if params["k"] == 2 and len(rgd) == 5 and params["bdoc4"] == 2:
            return ("dihydroflavonol C3-C2' Ether Linkage (R5)",
                    params["o4"], rgd)

        if params["k"] == 3 and len(rgd) == 6:
            return "rotenoid", None, rgd

        if params["k"] == 3 and len(rgd) == 5:

            if params["bdoc4"] == 1 and params["o4"] in rgd:
                if params["bdc34"] == 2:
                    if params["bdoc2"] == 2 and params["o2"] not in rgd:
                        return "coumestan", params["o2"], rgd
                    return "pterocarpan", None, rgd

                if params["bdc34"] == 1:
                    return "dihydropterocarpan", None, rgd

            if (params["bdoc4"] == 2
                    and params["bdoc2"] == 1
                    and params["o2"] in rgd):
                return "dihydrocoumaronochromone", params["o4"], rgd

    return None, None, None


def _stilbenoflavonoid(molobj, pskeleton):
    """
    Stilbeno-flavonoid. This is the derivative of flavan-3-ol, in which a
    stilbenoid is connected to C(4) of the flavonoid
    """
    benzylrings, atoms = molobj.benzylrings, molobj.atoms
    if pskeleton.name != "flavan-3-ol":
        return None, None

    for rg1 in benzylrings.keys():
        if (benzylrings[rg1].label != "b"
                or rg1 == pskeleton.ringA
                or rg1 == pskeleton.ringB):
            continue
        # the ring connects to ring C
        j, chains1 = None, {}
        for i in rg1:
            if pskeleton.corders[4] in atoms[i].neiset:
                j = i
            else:
                nx = atoms[i].neiset-benzylrings[rg1].index
                for nk in nx:
                    if atoms[nk].symbol == "c":
                        chains1[nk] = i
        if j is None:
            continue

        # the other ring
        for rg2 in benzylrings.keys():
            if (rg2 == rg1 or benzylrings[rg2].label != "b"
                    or rg2 == pskeleton.ringA
                    or rg2 == pskeleton.ringC):
                continue
            for i in rg2:
                nx = atoms[i].neiset - benzylrings[rg2].index
                if not nx:
                    continue
                cnns = [nei for j in nx for nei in atoms[j].neighbors
                        if atoms[nk].symbol == "c"
                        and nei.index in chains1
                        and nei.bondorder == 2
                        ]
                if len(cnns) != 1:
                    continue
                cnei = cnns[0]
                return ("stilbeno-flavonoid",
                        (rg1, j, chains1[cnei.index], cnei.index, nk, i, rg2)
                        )

    return None, None
