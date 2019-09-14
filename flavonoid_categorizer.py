"""
Categorize flavonoids according to the input molecule.
"""
import itertools
from indigo import Indigo
from itertools import combinations, product

from chalcone_class import Chalcone


idg = Indigo()
CHAIN = idg.CHAIN
RING = idg.RING


def getskidx(ringidx, sk):
    """ get the indices of atoms in identified skeletons """
    skix = []
    for ski in sk:
        tix = []
        for i in ski[0]:
            tix += ringidx[i]
        tix += ski[1:-1]
        skix.append(tix)
    return skix


# TO DO
def _checkXanthSK(molobj):
    """
    Get skeleton of Xanthones.
    
    """
    
    sk = []
    nr, setrg = len(ringidx), [set(rg) for rg in ringidx]
    cringix = []
    for i in xrange(nr):
        if all(fatom[j][1]=='c' for j in ringidx[i]):
            cringix.append(i)

    for i,j in iters(bridx, cringix, True):
        for k in xrange(nr):
            if len(setrg[i]&setrg[k])==2 and len(setrg[j]&setrg[k])==2\
               and not setrg[j]&setrg[i]:
                rc = setrg[k].copy()
                su = setrg[i].union(setrg[j])
                rc.difference_update(su)
                t1, t2 = False, False
                if len(rc) == 2:
                    for c in rc:
                        ci = [ni[0] for ni in fneis[c] if ni[0] in ringidx[i]]
                        cj = [ni[0] for ni in fneis[c] if ni[0] in ringidx[j]]
                            
                        if not ci or not cj: break
                            
                        ca = [ni[2] for ni in fneis[ci[0]] if ni[0]==c][0]
                        if ca == 'o':
                            t1, o, ci1, cj1 = True, c, ci[0], cj[0]
                        elif ca == 'c':
                            co2 = [ni[0] for ni in fneis[c] if ni[1]==2 and ni[2]=='o']
                            co2 = co2[0] if co2 else None
                            if co2:
                                t2, cc, ci2, cj2 = True, c, ci[0], cj[0]

                    if t1 and t2:
                        if j in bridx:
                            sk.append([[i,k,j], ci1,o,cj1,ci2,cc,cj2,co2, 'x'])
                        else:
                            rjk = list(setrg[j]&setrg[k])
                            b = [ni[1] for ni in fneis[rjk[0]] if ni[0]==rjk[1]][0]
                            if any(len(setrg[j]&setrg[jk])==2 for jk in bridx) and b==2:
                                sk.append([[i,k,j], ci1,o,cj1,ci2,cc,cj2,co2, 'x'])
                                
    return sk


def _dalbergione(molobj):
    """
    check the structure of dalbergione
    """
    benzylrings = molobj.benzylrings
    if not any(rg["label"]=="q" for rg in benzylrings.values()):
        return None

    # get the structure
    skeletons = []
    for rgix1, rg1 in benzylrings.values():
        if rg1["label"]!="q":
            continue
        cc1 = set()
        for j in rgix1:
            cc1.update(nei["index"] for nei in molobj.atoms[j]["neighbors"]
                       if nei["index"] not in rgix1 and nei["symbol"]=="c"])

        skj = _Skeleton()

        for rgix2, rg2 in benzylrings.items():
            if rg["label"] != "b":
                continue
            cc2 = set()
            for j in rgix2:
                cc2.update(nei["index"]
                           for nei in molobj.atoms[j]["neighbors"]
                           if nei["index"] not in rgix2
                              and nei["symbol"]=="c"])
            if len(cc1&cc2) != 1:
                continue
            # the only atom connecting the quinol and benzyl ring
            j, _is_candidate = list(cc1&cc2)[0], False
            for nei in molobj.atoms[j]["neighbors"]:
                if nei["index"] in rgix1:
                    k1 = nei["index"]
                elif nei["index"] in rgix2:
                    k2 = nei["index"]
                elif nei["symbol"]=="c":
                    # side chain
                    kc = nei["index"]
                    if len(molobj.atoms[kc]["neighbors"]) == 2:
                        kc2 = [nei for nei in molobj.atoms[kc]["neighbors"]
                               if nei["index"] != j and nei["degree"] == 2
                               and nei["topology"] == molobj.CHAIN]
                        if kc2:
                            kc2 = kc2[0]["index"]
                            _is_candidate = True
            if _is_candidate:
                skj.ringA = rgix1
                skj.ringB = rgix2
                skj.connects["atomb"] = j
                skj.connects["connectatoms"] = [k1, j, k2]
                skj.connects["topology"] = [molobj.CHAIN, molobj.CHAIN]
                skj.corders = [1, 1]
                skj.label = "d"
                skj.name = "dalbergione"
                skj.sideatoms = [kc, kc2]
                skeletons.append(skj)

    return skeletons


class FlavonoidClass():
    """
    Get skeleton and class name of the input molecule described by string
    such as canonical SMILES and InchI or molecular file.
    """
    def __init__(self, molobj):
        for key, value in molobj.__dict__.items():
            setattr(self, key, value)
        self.molobj = molobj
        self._ringABC()
        self._ringACindexing()
        self._ringBindexing()
        self._skeleton()
        self._tetrahydroflav()

    def _ringABC(self):
        """
        identify potential ring A, B and C
        call function ringABC(molobj)
        """
        self.ringAC, self.ringB = _getringABC(self.molobj)

    def _ringACindexing(self):
        """
        Index the identified potential ring A, C
        call function indexringAC(ringac)
        """
        self.ringACindex = _indexringAC(self.molobj, self.ringAC)

    def _skeleton(self):
        """
        flavonoid skeleton with structure C6-C3-C6
        """
        self.flavsk = _getsk(self.molobj, self.ringAC, self.ringB)

    def _flav(self):
        """
        assign names to flavonoid skeleton identified
        """
        self.flavsk = _flavClass(self.molobj, self.ringACindex, self.flavsk)

    def _chalcone(self):
        """
        assign names to chalcone
        """
        self.chalcone = Chalcone(self.molobj)
