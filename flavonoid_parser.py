"""
Parse flavonoid molecule object read from Indigo toolkit. The molecular
identifiers accepted are canonical SMILES, InchI or file like format
..mol (file stored in ChemSpider database), and ..sdf (in PubChem database).
"""

import collections
from re import findall, sub
from indigo import Indigo, IndigoException
from indigo_inchi import IndigoInchi
from indigo_renderer import IndigoRenderer as renderer

idg = Indigo()
idgin = IndigoInchi(idg)
idgrender = renderer(idg)

Atom = collections.namedtuple("Atom",
                              ["symbol", "charge", "neighbors", "neiset"])
Nei = collections.namedtuple("Nei",
                             ["index", "symbol", "bondorder",
                              "degree", "bondtopology"])
Ring = collections.namedtuple("Ring",
                              ["obj", "index", "smiles", "label", "supports",
                               "isringAC", "ringACnns"])


class FlavonoidException(Exception):
    """ flavonoid exceptions, to identify whether the input molecule
    is valid flavonoid candidate, and other related exceptions."""
    def __init__(self, errid):
        self.id = errid

    def __str__(self):
        if self.id == 1:
            return " ".join(["Multiple components are unsupported in",
                             "identifying flavonoids."])
        elif self.id == 2:
            return " ".join(["Elements except CHNOS are not allowed as a",
                             "candidate of flavonoid."])
        elif self.id == 3:
            return "Group matching times out."
        elif self.id == 4:
            return "This is not a flavonoid."
        elif self.id == 5:
            return "Can not load molecule from the input identifier by Indigo."
        elif self.id == 6:
            return " ".join(["Charged atoms except oxygen or charges except",
                             "+1 are not allowed as a flavonoid."])
        elif self.id == 7:
            return " ".join(["Unexcepted side groups, elements or radical",
                             "atoms exist in flavonoids."])
        elif self.id == 8:
            return " ".join(["Unexpected side group including nitrogens in",
                             "a flavonoid."])
        elif self.id == 9:
            return " ".join(["Too condensed distribution of benzene rings",
                             "which unlikely appears as a flavonoid."])
        elif self.id == 10:
            return " ".join(["Not sufficient rings are found as a flavonoid"
                             "candidate."])
        elif self.id == 11:
            return " ".join(["Multiple oxygens bonding to the same atom"
                             "in a ring is unacceptable as a flavonoid."])


def loadmol(identifier):
    """
    Load molecule using the identifier to construct an Indigo object.
    """
    # It's a molecular formatted file, only ..mol and ..sdf are accepted.
    if identifier.endswith(".mol") or identifier.endswith(".sdf"):
        try:
            mol = idg.loadMoleculeFromFile(identifier)
        except IndigoException:
            raise FlavonoidException(5)

    # Otherwise it's an identifier, SMILES or InchI
    try:
        mol = idg.loadMolecule(identifier)
    except IndigoException:
        try:
            mol = idgin.loadMolecule(identifier)
        except IndigoException:
            raise FlavonoidException(5)

    return mol


class MoleculeParser(object):
    """
    Parse molecular object for subsequent analysis
    """
    def __init__(self, identifier):
        mol = loadmol(identifier)
        mol.dearomatize()
        self.molecule = mol
        self.CHAIN = idg.CHAIN
        self.RING = idg.RING
        self._atoms()
        self._Rings()
        self._assignbenzylabel()
        self._RingA()

    def _atoms(self):
        """
        Get information of atoms and corresponding neighbor atoms in
        current molecule and set the obtained variables as global.
        The information for neighbors (in variable 'neis') include:
            .. index of neighbor atom;
            .. bond order the this neighbor atom to current atom;
            .. atom symbol in lower case;
            .. degree of neighbor atom.
        The information for atoms (in )
        """
        atoms = {}
        for atom in self.molecule.iterateAtoms():
            k = atom.index()
            neis = tuple([Nei(nei.index(),
                              nei.symbol().lower(),
                              nei.bond().bondOrder(),
                              nei.degree(),
                              nei.bond().topology())
                          for nei in atom.iterateNeighbors()])
            atoms[k] = Atom(atom.symbol().lower(), atom.charge(), neis,
                            set(nei.index for nei in neis))
        self.atoms = atoms

    def _Rings(self):
        """
        Get ring objects and atom indices in each ring
        """
        coset = set("co")
        # rings containing atom objects
        rings, benzylrings, thps = {}, {}, {}
        for i, r in enumerate(self.molecule.iterateRings(3, 8)):
            rjx = tuple([atom.index() for atom in r.iterateAtoms()])
            s = r.clone().smiles()
            setjx = set(rjx)

            tempobj = Ring(r, setjx, s, None, None, None, None)

            rings[rjx] = tempobj

            if 5 <= len(rjx) <= 6:
                # to avoid internal ring
                if any(len(self.atoms[j].neiset & set(rjx)) >= 3 for j in rjx):
                    continue
                # ring smiles
                sc = ''.join(findall("[a-zA-Z]+", s)).lower()
                if all(ak == "c" for ak in set(sc)):
                    # all are C and benzenes
                    benzylrings[rjx] = tempobj
                elif sc.count("o") == 1 and coset.issuperset(set(sc)):
                    # tetrahydrofuran and tetrahydropyran
                    thps[rjx] = tempobj

        self.rings = rings
        self.benzylrings = benzylrings
        self.thps = thps

    def _assignbenzylabel(self):
        """
        Assign names to benzyl rings
        """
        for rgix, rgobj in self.benzylrings.items():
            rsmile = rgobj.smiles
            # For a candidate benzene ring, at least one neighbor not in
            # the same ring is C
            if not any(any(nei.symbol == "c" and nei.index not in rgix
                           for nei in self.atoms[j].neighbors)
                       for j in rgix):
                continue

            label, supports = None, None
            # benzene ring
            n, db = len(rgobj.index), rsmile.count("=")
            if n == 6 or (db == 3 and n == 6):
                label = "b"
            # check whether the ring is a benzene if it is surrounded by other
            # benzene rings
            elif n == 6 and db >= 1:

                cbs = []
                for rgixk, rgobjk in self.benzylrings.items():
                    ck = rgobj.index & rgobjk.index
                    if len(ck) == 2 and rgobjk.smiles.count("=") > 2:
                        cbs.append(ck)
                if (len(cbs) > 0 and db == 2) or (len(cbs) > 1 and db == 1):
                    label = "b"
                elif db == 2:
                    # candidate benzene ring
                    label = "bx"

            # other types of aromatic rings
            if label is None or label == "bx":
                c2o = []
                for j in rgix:
                    for nei in self.atoms[j].neighbors:
                        if nei.symbol == "o" and nei.bondorder == 2:
                            c2o.append(nei.index)
                if db >= 1:
                    if db == 2 and len(c2o) == 2:
                        # quinone
                        if n == 6:
                            label, supports = "q", tuple(c2o)
                    elif db == 2 and c2o:
                        # methyldeoxybenzoin if 6 membered ring
                        # else chalcone cyclopentenedione
                        label, supports = "m" if n == 6 else "cc", c2o[0]
                    elif db == 1 and len(c2o) == 2:
                        label, supports ="q2" if n == 6 else "cc", tuple(c2o)
                else:
                    if len(c2o) == 3 and not set(c2o) & rgobj.index:
                        # cyclohexane-1,3,5-trione
                        label, supports ="c", tuple(c2o)

            if label is not None:
                rgobj = rgobj._replace(label=label)
                rgobj = rgobj._replace(supports=supports)

                self.benzylrings[rgix] = rgobj

        # recheck those unidentified or candidate benzene rings
        while 1:
            t = True
            for rgix, rgobj in self.benzylrings.items():
                if rgobj.label is None or rgobj.label == "bx":
                    _n_adj, label = 0, rgobj.label
                    for rgixk, rgobjk in self.benzylrings.items():
                        ck = rgobj.index & rgobjk.index
                        if len(ck) == 2 and rgobjk.label == "b":
                            _n_adj += 1
                    if rgobj.label == 'bx' and _n_adj > 0:
                        label, t ="b", False
                    elif rgobj.label is None and _n_adj == 3:
                        label, t = "b", False

                    if label == "bx":
                        label = None

                    rgobj = rgobj._replace(label=label)
                    self.benzylrings[rgix] = rgobj
            if t:
                break

        return self

    def _RingA(self):
        """
        Identify candidate ring A, i.e., a benzene like ring with
        side chains -C- and -O- bonding same benzene ring bond.
        Return information around ring A:
            index of benzene ring that ring A belongs to,
            object of atom in benzene ring that bonds to C,
            object of C,
            object of atom in benzene ring that bonds to O,
            object of O.
        """
        for rgix, rgobj in self.benzylrings.items():
            for j in rgix:
                # C(10)-C(4)
                cc = [nei.index for nei in self.atoms[j].neighbors
                      if nei.index not in rgix and nei.symbol == "c"]
                if not cc:
                    continue
                for k in self.atoms[j].neiset&rgobj.index:
                    # C(9)-O(1)
                    bo = [nei.index for nei in self.atoms[k].neighbors
                          if nei.symbol == "o"]
                    if not bo:
                        continue
                    rgobj = rgobj._replace(isringAC=True)
                    rgobj = rgobj._replace(ringACnns=tuple([(j, jk, k, bo[0])
                                                            for jk in cc]))
            self.benzylrings[rgix] = rgobj

        return self

    def getoxygen(self, k):
        """
        Get oxygen atom bonded to atom k, return the atom index,
        ring order and ring topology
        """
        neioxys = []
        for nei in self.atoms[k].neighbors:
            if nei.symbol == 'o':
                neioxys.append((nei.index, nei.bondorder, nei.bondtopology))

        if len(neioxys) > 1 and any(bdt == self.RING for _, _, bdt in neioxys):
            raise FlavonoidException(11)

        return neioxys if neioxys else None

    def getbondorder(self, j, k):
        """
        Get order of the bond between atom j and k
        """
        for nei in self.atoms[j].neighbors:
            if nei.index == k:
                return nei.bondorder
        # if atom j and k are not connected, return None
        return None


class ValidFlavonoidCheck(object):
    """
    Check whether the input identifier is valid as a flavonoid.
    All attributes are inheritated from class MoleculeParser which
    creates an Indigo object for the molecule and parses the object
    to get necessary information for facilitating the subsequent
    analysis. However, prior to that, the check of the input molecule
    to be a valid flavonoid is necessary.
    """
    def __init__(self, molobj):
        for key, value in molobj.__dict__.items():
            setattr(self, key, value)
        self._isvalidelements()
        # self._isvalidcomponentcomposition()
        self._isradical()
        self._isvalidnitrogen()
        self._isvalidringdist()
        self._isvalidmolecule()

    def _isvalidelements(self):
        """
        raise exception if the element composition contains exceptional
        elements as a flavonoid.
        """
        smiles = self.molecule.canonicalSmiles()
        # get all letters with lower representing
        letts = "".join(findall("[a-zA-Z]+", smiles)).lower()
        # remove letters associated with elements composing a flavonoid
        # molecule and identify whether any letter still exists, if
        # yes, raise the exception, i.e., this is not a flavonoid.
        if sub("[chons]", "", letts):
            raise FlavonoidException(2)

    def _isvalidcomponentcomposition(self):
        """
        raise exception if too many components are found in the molecule.
        """
        if self.molecule.countComponents() > 1:
            raise FlavonoidException(1)

    def _isradical(self):
        """
        raise exception if it is a radical
        """
        # raise exception if it is a radical
        if any(atom.radicalElectrons()
               for atom in self.molecule.iterateAtoms()):
            raise FlavonoidException(7)

    def _isvalidnitrogen(self):
        """
        raise exception if nitrogen linkage is invalid as a flavonoid:
        -- N-N or N-O linkage is not allowed
        -- double bond with one destination as N is not allowed except
            that the bond locates in a ring
        """
        for atom in self.atoms.values():
            if atom.symbol != "n":
                continue
            for nei in atom.neighbors:
                if nei.symbol != "c":
                    raise FlavonoidException(8)
                if nei.bondorder >= 2 and nei.bondtopology == idg.CHAIN:
                    raise FlavonoidException(8)

    def _isvalidringdist(self):
        """
        check whether distribution of benzene rings is valid in the
        molecule as a flavonoid, if three benzene rings bonding
        together(i.e., having shared bonds between any two) or more
        then 2 bonded benzene groups exist, raise exception.
        """
        setrg = [rg.index for rg in self.benzylrings.values()
                 if rg.label is not None]
        phcluster = None
        for i, rgi in enumerate(setrg):
            bci = set([j for j, rgj in enumerate(setrg)
                       if rgi & rgj and j != i])
            if len(bci) >= 2:
                if phcluster is not None:
                    phcluster = bci
                elif phcluster != bci:
                    raise FlavonoidException(9)

    def _isvalidmolecule(self):
        """
        check whether the molecule is valid.
        """
        if len(self.rings) == 0 or len(self.benzylrings) == 0:
            raise FlavonoidException(10)
