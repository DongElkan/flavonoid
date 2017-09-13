"""
This module interprets canonical SMILES to get skeleton and substituent information
of flavonoid for automatic fragmentation of the molecule. The obtained fragments are
applied to annotating MS/MS spectrum. Indigo library
    lifescience.opensource.epam.com/indigo/index.html
with python wrapper is adopted without any modification.

LICENSE
This program is under the terms adopted by Indigo library, i.e. GNU General Public
License, version 3.

AUTHOR
Naiping Dong

EMAIL
naiping.dong@hotmail.com
"""

import os, json, sys
from re import findall
from indigo import Indigo
from indigo_inchi import IndigoInchi
from indigo_renderer import IndigoRenderer as renderer

idg         = Indigo()
idgin       = IndigoInchi(idg)
idgrender   = renderer(idg)
CHAIN = idg.CHAIN

# exception names that will not perform indexing
name_exceptions = ['theaflavin', 'biflavonoid', 'stilbeno-flavonoid']

with open('sidechain.json','r') as f:
    sideChains = json.load(f)
    

class FlavonoidException(Exception):
    """ flavonoid exceptions """
    def __init__(self, errid):
        self.id = errid
        
    def __str__(self):
        if self.id==1:
            return 'Multiple components are unsupported in identifying flavonoids.'
        elif self.id==2:
            return 'Elements except CHNOS are not allowed as a candidate of flavonoid.'
        elif self.id==3:
            return 'Group matching times out.'
        elif self.id==4:
            return 'This is not a flavonoid.'
        elif self.id==5:
            return 'Can not load molecule from the input identifier.'
        elif self.id==6:
            return 'Charged atoms except oxygen or charges except +1 are not allowed as a candidate of flavonoid.'
        elif self.id==7:
            return 'Unexcepted side groups, elements or radical atoms exist in flavonoids.'
        elif self.id==8:
            return 'Unexpected side group including nitrogens in a flavonoid.'
        elif self.id==9:
            return 'Too condensed distribution of benzene rings which unlikely appears as a flavonoid.'            



def mis(l):
    """
    Maximum independent set of a list.
    To solve this problem, a brute force procedure is adopted with some constraints
    which can be used to reduce the number of combinations.
    
    Command:
        mis_list,idx=MIS(l)
        
    Output:
        list contains largest independent list set and indices of the independent list.
    """
    n = len(l)
    rdidx = sorted(range(n), key=lambda k: l[k])
    l.sort()
    
    # Get the indices of lists independent from other lists
    # matrix containing the indices of lists in list "l" that have uncommon elements
    # comparing with current list, e.g. ith list
    M=[]
    sets = [set(li) for li in l]
    for i in xrange(n):
        M.append([j for j in xrange(i,n) if not sets[i]&sets[j]])

    # Find the largest independent list set
    C, t = [], []
    for i in xrange(n-1):
        t[:]=[[i]]
        ct=[] # current path
        while True:
            stop=True # identify whether any new index is found, if not, stop
            for j in xrange(len(t)):
                step=t[j][-1]
                if step<n-1:
                    for k in M[step]:
                        # this list must be independent from all lists have been
                        # selected, i.e., in t[j]
                        if all(k in M[v] for v in t[j]):
                            ct.append(t[j]+[k])
                            stop=False

            if stop:
                break
            else:
                t[:]=ct[:]
                ct[:]=[]

            if all(j[-1]==n for j in t):
                break

        C+=t

    # Get the largest independent list set
    ll=[len(v) for v in C]  # length of independent list set
    idx=ll.index(max(ll))   # index of largest independent list set
    return [(l[i], rdidx[i]) for i in C[idx]]


def numberprefix(n):
    """ Greek prefixes of number """
    greekprefix = ['mono-','di-','tri-','tetra-','penta-','hexa-','hepta-','octa-',
                   'nona-', 'deca-','hendeca-','dodeca-','triskaideca-',
                   'tetrakaideca-','pentakaideca-','hexadeca-','heptakaideca-',
                   'octakaideca-']
    return greekprefix[n]


def validmolcheck(smiles):
    """ check whether input canonical smile for flavonoid is valid """
    # get all letters with lower representing
    letts = ''.join(findall("[a-zA-Z]+", smiles)).lower()
    # remove letters associated with elements composing a flavonoid molecule
    # and identify whether any letter still exists, if yes, return FALSE, i.e.,
    # this is not a flavonoid
    if letts.translate(None, "chons"):
        return False
    return True


def validnitrogencheck(molobj):
    """
    check whether nitrogen linkage is valid in the molecule:
    -- N-N or N-O linkage is not allowed
    -- double bond with one destination as N is not allowed except that
       the bond locates in a ring
    """
    neic = []
    for atom in molobj.iterateAtoms():
        r = False
        if atom.symbol().lower()=='n':
            neic[:] = [(ni.bond().bondOrder(), ni.symbol(), ni.bond().topology(), ni)
                         for ni in atom.iterateNeighbors()]
            for d, s, p, a in neic:
                if s.lower()!='c':
                    return False
                if d>=2 and p==CHAIN:
                    return False
    
    return True


def validringdistcheck(ringidx, aromaticName):
    """
    check whether distribution of benzene rings is valid in the molecule
    as a flavonoid, if three benzene rings bonding together (i.e., having
    shared bonds between any two) or more then 2 bonded benzene groups
    exist, raise exception.
    """
    nr, setrg = len(ringidx), [set(rg) for rg in ringidx]
    bcs = []
    for i in xrange(nr):
        if not aromaticName[i]: continue
        bci = [j for j in xrange(nr) if setrg[i]&setrg[j] and aromaticName[j] and j!=i]
        if len(bci)>=2:
            if not bcs:
                bcs.append(set(bci))
            elif not any(bi.issuperset(bci) or bi.issubset(bci) for bi in bcs):
                bcs.append(set(bci))

    if len(bcs) > 1:
        raise FlavonoidException(9)
    

def loadmol(string):
    """ Load molecule from the input string to identify whether it is valid,
        this can be a canonical SMILES, InchI identifier, or a valid molecular
        file (e.g., ..mol or ..sdf).
        If the input is a file name, full path should be input if the file
        does not locate under same directory as this code file.
    """

    if string.endswith('.mol') or string.endswith('.sdf'):
        return idg.loadMoleculeFromFile(string)
        
    try:
        return idg.loadMolecule(string)
    except:
        try:
            return idgin.loadMolecule(string)
        except:
            raise FlavonoidException(5)


def generatemolimage(string, filename):
    """ generate image for molecule 'string' """
    mol = loadmol(string)
    idgrender.renderToFile(mol,filename)


def getinchi(string):
    """ get inchi identifier of the molecule """
    mol = loadmol(string)
    return idgin.getInchi(mol)


def getcanonicalsmiles(string):
    """ get canonical SMILES of the molecule """
    mol = loadmol(string)
    return mol.canonicalSmiles()


def getnames(names):
    """ get name of the molecule """
    if isinstance(names,str):
        return names
    
    elif all(x=='flavan-3-ol' for x in names) and len(names)==2:
        return 'proanthocyanidin'
    
    elif len(names) > 1:
        if any('adduct' in s for s in names):
            if all('chalcone' in s for s in names):
                return 'chalcone %smer'%(numberprefix(len(names)-1))[:-1]
            else:
                nc = [s for s in names if 'chalcone' in s]
                nf = [s for s in names if 'chalcone' not in s]
                if len(nc) == 1 and len(nf) == 1:
                    return 'chalcone-%s'%nf[0]
                elif len(nf) == 1:
                    return 'chalcone %smer-%s'%(numberprefix(len(nc)-1)[:-1], nf[0])
                else:
                    nn1, nn2 = numberprefix(len(nc)-1)[:-1], numberprefix(len(nf)-1)
                    return 'chalcone %s-%sflavonoid'%(nn1, nn2)
        return '%sflavonoid: '%(numberprefix(len(names)-1))+', '.join(names)
    
    elif len(names)==1:
        if names[0]=='chalcone_adduct':
            return 'chalcone adduct'
        elif names[0]=='chalcone_da_adduct':
            return 'chalcone Diels-Alder adduct'
        return names[0]
    
    
def miters(iter1,iter2,flag=False):
    """
    iterator for pair iterating of two iterables.
    """
    if flag:
        return iter([(x,y) for x in iter1 for y in iter2 if x!=y])
    else:
        return iter([(x,y) for x in iter1 for y in iter2])


def numrgass(assrings, aromaticName):
    """ Return number of rings assigned """
    t = []
    for info in assrings:
        t += info[0]
    return len(set([i for i in t if aromaticName[i]]))


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


def atominfo(mol):
    """
    Get information of atoms and corresponding neighbor atoms in current
    molecule and set the obtained variables as global.
    The information for neighbors (in variable 'neis') include:
        .. index of neighbor atom;
        .. bond order the this neighbor atom to current atom;
        .. atom symbol in lower case;
        .. degree of neighbor atom.
        
    The information for atoms (in )
    """
    global neis, atoms
    neis, atoms = [], []
    for atom in mol.iterateAtoms():
        neis.append([(nei.index(),
                      nei.bond().bondOrder(),
                      nei.symbol().lower(),
                      nei.degree(),
                      nei.bond().topology())
                     for nei in atom.iterateNeighbors()])
        atoms.append([atom.index(), atom.symbol().lower(), atom.charge()])
    neis = tuple(neis)
    atoms = tuple(atoms)


def benzeneRingID(ringObjs, ringIdx):
    """
    Identify whether a ring is a benzene ring or a ring with conjugated structure,
    i.e., a ring has a phenolic structure
    """
    # pre-declaration to faster the use of global variable
    fneis = neis

    nrg, setrg = len(ringIdx), [set(rg) for rg in ringIdx]
    rsmiles = [r.clone().smiles() for r in ringObjs]

    aromaticName = [None]*nrg

    brc = []
    for i in xrange(nrg):
        brc.append(rsmiles[i].lower().count('c')==ringObjs[i].countAtoms())
    
    for i in xrange(nrg):

        if not brc[i]: continue

        f = False
        for j in ringIdx[i]:
            if any(ni[2]=='c' and ni[0] not in ringIdx[i] for ni in fneis[j]):
                f = True
                break
        # For a candidate benzene ring, at least one neighbor not in the
        # same ring is C
        if not f: continue

        rsmile = rsmiles[i]
        
        if rsmile.count('c')==6:
            aromaticName[i] = 'b'
            continue

        nc, db, tp = rsmile.count('C'), rsmile.count("="), None
                
        if db==3 and nc==6:
            tp = 'b'

        # check whether the ring is a benzene if it is surrounded by other
        # benzene rings
        if not tp and nc == 6 and db >= 1:
            cbs = []
            for j in xrange(nrg):
                if len(setrg[i]&setrg[j]) == 2 and rsmiles[j].count("=") > 2:
                    cbs.append(setrg[i]&setrg[j])

            t = True
            for j in setrg[i]:
                if not any(j in cbi for cbi in cbs):
                    k = [ni[0] for ni in fneis[j] if ni[0] in setrg[i] and ni[1]==2]
                    if any(cbi.issuperset(k) for cbi in cbs):
                        t = False
                        break
            if t:
                if len(cbs)>0 and db==2:
                    tp = 'b'
                elif len(cbs)>1 and db==1:
                    tp = 'b'
                elif db==2:
                    tp = 'bx' # candidate benzene ring
                    
        if not tp or tp=='bx':
            if db>=1:
                b2 = [b for b in ringObjs[i].iterateBonds() if b.bondOrder() == 2]
                a2idx = [b.destination().index() for b in b2]+\
                        [b.source().index() for b in b2]
                neio, neic = [], []
                for j in ringIdx[i]:
                    neio += [v[1]==2 for v in fneis[j] if v[2] == 'o']
                    neic += [v[0] for v in fneis[j] if v[0] in a2idx]
                    
                if nc==6:
                    if db==2 and sum(neio)==2 and len(set(neic))==4:
                        tp = 'q'                              # quinone
                    elif db==2 and any(neio):
                        tp = 'm'
                    elif db==1 and sum(neio)==2:
                        tp = 'q2'
                else:
                    if neio and all(neio) and len(neio)==2 and db==1:
                        tp = 'cc'                             # chalcone cyclopentenedione
                    elif any(neio) and db == 2:
                        tp = 'cc'
                    
            else:
                c2o, c2nei = [], []
                for j in ringIdx[i]:
                    c2o += [v[0] for v in fneis[j] if v[2] == 'o' and v[1]==2]
                    c2nei += [v[0] for v in fneis[j] if v[0] in ringIdx[i]]
                if len(c2o) == 3 and len(set(c2o)&set(c2nei)) == 0:
                    tp = 'c'                                  # cyclohexane-1,3,5-trione

        aromaticName[i] = tp

    # recheck those unidentified or candidate benzene rings
    while True:
        t = True
        for i in xrange(nrg):
            if brc[i] and (not aromaticName[i] or aromaticName[i]=='bx'):
                nn = 0
                for j in xrange(nrg):
                    if len(setrg[i]&setrg[j])==2 and aromaticName[j]=='b':
                        nn += 1
                if aromaticName[i] == 'bx' and nn > 0:
                    aromaticName[i] = 'b'
                    t = False
                elif not aromaticName[i] and nn == 3:
                    aromaticName[i] = 'b'
                    t = False

                if aromaticName[i] == 'bx': aromaticName[i] = None

        if t: break

    return aromaticName


def getRings(mol, checkBenzene = True):
    """ Get ring objects and indices for benzenes or benzene like rings """
    
    ringidx =[]           # rings containing atom objects
    
    if checkBenzene:

        ringObjs = []
        for r in mol.iterateRings(5,6):
            s = r.clone().smiles().lower()
            s = ''.join(findall("[a-zA-Z]+", s))
            if all(k in 'co' for k in s):
                rix = [atom.index() for atom in r.iterateAtoms()]
                # to avoid internal ring
                t = True
                for atmi in r.iterateAtoms():
                    if sum(ni.index() in rix for ni in atmi.iterateNeighbors())>=3:
                        t = False
                        break
                if t:
                    ringidx.append(rix)
                    ringObjs.append(r)

        if len(ringidx)<2:  return None, None

        aromaticName = benzeneRingID(ringObjs, ringidx)

        return ringidx, aromaticName
    
    else:
        
        for r in mol.iterateRings(3,8):
            ringidx.append([atom.index() for atom in r.iterateAtoms()])
        return ringidx


def getRingA(ringidx, aromaticName):
    """ Get candidate ring A, i.e., a benzene like ring with side chains -C- and -O-
        bonding same benzene ring bond.
        Return information around ring A: index of benzene ring that ring A belongs
        to, object of atom in benzene ring that bonds to C, object of C, object of
        atom in benzene ring that bonds to O, object of O.
    """
    # pre-declaration to faster the use of global variable
    fneis = neis
    
    ringA = []
    for i in xrange(len(ringidx)):
        if not aromaticName[i]:
            continue
        
        idxs = ringidx[i]
        for j in idxs:
            # C(10)-C(4)
            cc = [v[0] for v in fneis[j] if v[0] not in idxs and v[2]=='c']
            if cc:
                for nei in fneis[j]:
                    # C(9)-O
                    bo = [v[0] for v in fneis[nei[0]] if v[2]=='o' and nei[0] in idxs]
                    if bo:
                        for k in cc:
                            ringA.append([i, j, k, nei[0], bo[0]])
                    
    return ringA


def classFlav(ringidx, sk, aromaticName):
    """
    Classify flavonoids into subclasses
    """
    if isinstance(sk,str):
        if sk=='c':
            return 'chalcone'
        elif sk=='dc':
            return 'dihydrochalcone'
        elif sk=='cq':
            return 'chalcone-quinol'
        elif sk=='oc':
            return 'oxodihydrochalcone'
        elif sk=='pc':
            return 'pentahydrochalcone'
        elif sk=='cc':
            return 'chalcone_cyclopentenedione'
        elif sk=='au':
            return 'aurone'
        elif sk=='da':
            return 'dihydroaurone'
        else:
            return None
    
    if len(sk[0])>4: return None
    
    # pre-declaration to faster the use of global variable
    fneis, fatoms = neis, atoms

    if sk[-1] == 'b':
        return '2-arylbenzofuran'
    elif sk[-1] == 'm':
        return 'alpha-methyldeoxybenzoin'
    elif sk[-1] == 'd':
        return 'dalbergione'
    elif sk[-1] == 'fl':
        return 'flavonolignan'
    elif sk[-1] == 'x':
        return 'xanthone'
    elif sk[-1] == 'arl':
        return '4-arylchroman'
    elif sk[-1] == 'fd':
        return 'flavonoid-derivative'

    tsk = sk[0]             # to avoid in-place replacement
    oxyset = set([atom[0] for atom in fatoms if atom[1]=='o'])
    
    idxA, idxC = ringidx[tsk[0]], ringidx[tsk[1]]
    brix = [i for i in xrange(len(ringidx)) if aromaticName[i]]
    
    # reindex of ring C, starting from O; C in ring A bonding to O is
    # indexed as 2, then other C in ring A is 3, ..., and C in ring C bonding
    # to O is 6
    noC = idxC[:]
    for i in xrange(6):
        if idxC[i] in oxyset:
            noC[i] = 1
        elif idxC[i] == sk[1]:
            noC[i] = 3
        elif idxC[i] in idxA:
            noC[i] = 2
        elif idxC[i] == sk[2]:
            noC[i] = 4
        elif any(nei[0]==sk[2] for nei in fneis[idxC[i]]):
            noC[i] = 5
        else:
            noC[i] = 6

    # bond order in ring C
    bnoC = idxC[:]
    allb1 = True
    for i in xrange(6):
        j = noC.index(i+1)
        ncs = [a for a in fneis[idxC[j]]]
        k = idxC[noC.index(i+2)] if i < 5 else idxC[noC.index(1)]
        b = [a[1] for a in ncs if a[0]==k]
        if any(idxC[j] in ringidx[ii] and k in ringidx[ii] for ii in brix):
            b[0] = 2
        if i >= 3:
            bt = [a[1] for a in ncs if a[2]=='o' and a[0] not in idxC and a[3]<2]
            if i == 3:
                # bond order of C(4)-O and whether bond C(4)-C(3) is
                # double (=2) or aromatic (=4)
                b4, bc4 = bt[0] if bt else None, b[0]==2 or b[0]==4
            if i == 4:
                # bond order of C(3)-O and C(3)-C(2)
                b5, bc5 = bt[0] if bt else None, b[0]==2 or b[0]==4
            if i == 5:
                # bond order of C(2)-O and C(2)-O
                b6, bc6 = bt[0] if bt else None, b[0]==2 or b[0]==4
            allb1 = allb1 and b[0]==1

    # assign name to flavonoids
    i, j, lj = idxC.index(sk[3]), tsk[-1], len(ringidx[tsk[-1]])
    if len(tsk) == 4:
        oc = list(oxyset&set(ringidx[j]))
        if oc:
            ot = any(ni[0] in idxC for ni in fneis[oc[0]])
            if noC[i] == 6:
                if ot:
                    if lj == 6:
                        return 'peltogynoid' if allb1 else 'dehydropeltogynoid'
                    elif lj == 5:
                        if allb1:
                            return 'dihydroflavonol C3-C2\' ether linkage (R5)'
                        elif bc5:
                            return 'flavone C3-C2\' ether linkage (R5)'
                else:
                    if lj == 6:
                        if allb1:
                            return 'dihydroflavonol C3-C2\' ether linkage (R6)'
                        elif bc5:
                               return 'flavone C3-C2\' ether linkage (R6)'
                    elif lj == 7:
                        if any(ni[0] in ringidx[tsk[2]] for ni in fneis[oc[0]]):
                            if allb1:
                                return 'dihydroflavonol C3-C2\' ether linkage (R7)'
                            elif bc5:
                               return 'flavone C3-C2\' ether linkage (R7)'
            elif noC[i] == 5 and oc:
                if lj == 6 and not ot:
                    return 'dehydrorotenoid' if bc5 else 'rotenoid'
                elif lj == 5 and ot:
                    if b6==2 and bc4:
                        return 'coumestan'
                    else:
                        return 'pterocarpan' if bc4 else 'dihydropterocarpan'
                elif lj == 5 and b4==2 and ot:
                    if bc5:
                        return 'coumaronochromone'
                    else:
                        return 'dihydrocoumaronochromone'
    if len(tsk) == 4: return None

    typeA, typeB = aromaticName[tsk[0]], aromaticName[tsk[2]]
    allbr = typeA=='b' and typeB=='b'
    
    if sk[-1] == 'hf' and allbr:
        bc = [nei[1] for nei in fneis[sk[3]] if nei[0]==sk[4]][0]
        if b4==2 and bc5:
            return 'homoisoflavone'
        elif allb1 and bc != 2:
            return 'homoisoflavanone' if b4==2 else 'homoisoflavan'
        elif bc4:
            return 'homoisoflavene'
        elif b4==2 and bc == 2:
            return '3-benzylidene-4-chromanone'
            
    if noC[i]==6:
        if typeB=='q' and typeA=='b' and allb1:
            return 'flavanquinone'
        elif allbr:
            if b4==2 and bc5:
                return 'flavonol' if b5==1 else 'flavone'
            elif b4==2 and allb1:
                return 'dihydroflavonol' if b5==1 else 'flavanone'
            elif b4==1 and allb1:
                return 'flavan-3,4-diol' if b5==1 else 'flavan-4-ol'
            elif bc4 and bc6:
                return 'anthocyanidin'
            elif bc4:
                return 'flav-3-ene'
            elif allb1:
                return 'flavan-3-ol' if b5==1 else 'flavan'
        
    elif noC[i]==5:
        
        if b4==2 and bc5 and allbr:
            return 'isoflavone'
        elif b4==2 and allb1 and allbr:
            return 'isoflavanone'
        elif b4==1 and allb1 and allbr:
            return 'isoflavanol/isoflavan-4-ol'
        elif allb1:
            if typeA=='q' or typeB=='q':
                return 'isoflavanquinone'
            return 'isoflavan'
        elif bc4 and allbr:
            return '3-arylcoumarin' if b6 and b6==2 else 'isoflav-3-ene'
        
    elif noC[i]==4:
        
        if b6==2 and bc4 and allbr:
            return 'neoflavonoid'

    return None


def checkStilbenoSK(ringidx, brix, sk, name):
    """ Check stilbeno skeleton to check whether it exists, which is a candidate
        stilbeno-flavonoid
        If Stilbeno group is found, the indices of rings are added to the ring
        and indices of
            .. C in ring C bonds to Stilbene
            .. C in benzene ring that bonds to ring C, denoted as C(1)
            .. C in benzene ring that bonds to C(1), denoted as C(2)
            .. C in chain that bonds to C(2), which as double bond
            .. C in chain that bonds to above C
            .. C in another benzene ring
        and
            .. indicator of the skeleton
    """
    # pre-declaration to faster the use of global variable and function 
    fneis = neis
    iters = miters

    a, assbr = -1, []
    for skj in sk:
        a += 1
        if not name[a] != 'flavan-3-ol' or len(skj[0]) != 3: continue
        
        for i,j in iters(brix, brix, flag=True):

            if i in skj[0] or j in skj[0]: continue
            if i in assbr or j in assbr: continue

            # the benzene ring that bonds to ring C of skj
            ki = skj[2]
            ix = [ni[0] for ni in fneis[ki] if ni[0] in ringidx[i] and ni[4]==CHAIN]
            
            if not ix: continue
            ix = ix[0]

            atomnei = [ni[0] for ni in fneis[ix] if ni[0] in ringidx[i]]
            for ix, iy in iters(atomnei, ringidx[j]):
                for n1, n2 in iters(fneis[ix], fneis[iy]):
                    c = [ni for ni in fneis[n1[0]] if ni[0]==n2[0]]
                    g = True
                    for rix in ringidx:
                        if n1[0] in rix or n2[0] in rix:
                            g = False
                            break
                    if c and c[0][1]==2 and g:
                        skj.remove(skj[-1])
                        skj[0] += [i,j]
                        skj += [idx[0], k, ix, nei1[0], c[0][0], iy, 'sf']
                        name[a] = 'stilbeno-flavonoid'
                        assbr += [i, j]
    
    return sk, name


def checkBiflavonoidSK(ringidx, brix, ringA):
    """
    Check for biflavonoids which have complex structure like Lophirone A.
    """
    # pre-declaration to faster the use of global variable and function 
    iters = miters
    fneis, fatoms = neis, atoms
    
    sk = []
    
    # get tertiary carbon
    d3 = []
    for i in xrange(len(fneis)):
        if len(fneis[i])>=3 and fatoms[i][1] == 'o' and\
           len([ni for ni in fneis[i] if ni[2]=='c' and ni[1]==1])>=3:
            d3.append(i)

    # get ring index
    cridx, rno, setrg = [], [], []
    for idx in ringidx:
        cridx += idx
        rno.append(len([i for i in idx if fatoms[i][1]=='o']))
        setrg.append(set(idx))

    d3 = [i for i in d3 if i not in cridx]
    if len(d3) < 2: return sk, None

    # get bonded tertiary carbons and two benzene rings bonding to one
    # of the tertiary carbon
    p3 = []
    for i, j in miters(d3, d3, True):
        if any(ni[0] == j for ni in fneis[i]) and any(ni[0] == i for ni in fneis[j]):
            bz3 = []
            for k in brix:
                if any(ni[0] in ringidx[k] for ni in fneis[i]):
                    bz3.append(k)

            if len(bz3) == 2:
                p3.append((i, j, bz3))
    if not p3: return sk, None

    # get chromane ring
    cr = []
    for i in brix:
        for j in xrange(len(ringidx)):
            if i != j and len(setrg[i]&setrg[j])==2 and rno[j]==1:
                c = setrg[i]&setrg[j]
                nx = 0
                for k in c:
                    if any(ni[0] in ringidx[j] and ni[2]=='o' for ni in fneis[k]):
                        nx += 1
                    else:
                        c2=[ni for ni in fneis[k] if ni[0] in ringidx[j] and ni[2]=='c']
                        if c2:
                            c2 = c2[0][0]
                            if any(ni[1]==2 and ni[2]=='o' for ni in fneis[c2]):
                                nx += 1
                if nx == 2:
                    cr.append((j,i))
    if not cr: return sk, None

    # get other two rings
    for i, j, bzs in p3:
        for c in cr:
            k = c[0]
            
            if not any(ni[0] in ringidx[k] for ni in fneis[j]): continue
            
            c2 = [ni[0] for ni in fneis[j] if ni[0] not in ringidx[k] and ni[2]=='c']
            for l in c2:
                for m in brix:
                    if m!=i and m!=j and m!=c[1] and\
                       any(ni[0] in ringidx[m] for ni in fneis[l]):
                        sk += [[c[1], k, m, bzs[0], bzs[1]], i,j]
                        sk += [ni[0] for ni in fneis[i] if ni[0] in ringidx[bzs[0]]]
                        sk += [ni[0] for ni in fneis[i] if ni[0] in ringidx[bzs[1]]]
                        sk.append(l)
                        sk += [ni[0] for ni in fneis[l] if ni[0] in ringidx[m]]
                        sk.append('bi')
                        return sk, 'biflavonoid'
        
    return sk, None


def checkTheaflavinSK(sk, ringidx, brid, ringA, ring7, names):
    """
    Check whether it is the skeleton of a theaflavin
    ringA.append([i, j, k, nei[0], bo[0]])
    """
    # pre-declaration
    iters = miters
    fneis, fatoms = neis, atoms

    setrg, nr = [set(r) for r in ringidx], len(ringidx)

    # get rings containing 7 atoms with aromatic structures
    br7 = []
    for r in ring7:
        
        if all(nei[1]==4 for i in r for nei in fneis[i] if nei[0] in r):
            br7.append(set(r))
        else:
            k = 0
            for i in r:
                k += sum([nei[1]==2 and nei[0] in r for nei in fneis[i]])
            if k == 6:
                br7.append(set(r))

    if not br7: return sk, names
    
    # get ring A of Chromone that bonds to benzyl ring 7
    tx = []
    for i in xrange(len(names)):
        if names[i] != 'flavan-3-ol': continue
        for rA in ringA:
            j, t = rA[0], rA[1:]
            
            if j in sk[i][0]: continue
            
            # get ring C of the other flavan-3-ol like structure
            tag = False
            for k in xrange(nr):
                if setrg[k].issuperset(t) and len(setrg[j]&setrg[k])==2:
                    tag = True
                    break
            if not tag: continue

            # get ring B of br7
            tag = False
            for a, r in iters(ringidx[k], br7):
                b = [nei[0] for nei in fneis[a] if nei[0] in r and nei[-1]==CHAIN]
                if b:
                    b, tag = b[0], True
                    break
            if not tag: continue

            # check the identified ring B, it must be adjacent to ring
            # B of the skeleton i
            if len(r&setrg[sk[i][0][-1]]) != 2: continue

            # insert the identified structure to the flavan-3-ol skeleton
            ij = br7.index(r)
            sk[i][0] += [j, k, ij+nr]
            sk[i][:-1] += [t[0],t[1],a,b,t[2],t[3]]
            sk[i][-1] = 'tf'
            tx.append(i)

    # if only one skeleton is found, return anyway, or no theaflavin
    # is found, return
    if not tx: return sk, names
    
    # return idenfied theaflavin and ignore other skeletons since
    # in theaflavin unlikely exist other flavonoid skeletons
    return [sk[i] for i in tx], 'theaflavin'


def checkFlavSK(ringidx, aromaticName, ringA, ring7):
    """
    Check whether it is the skeleton of a flavonoid,  flavonolignan or
    homoisoflavonoids, and return indices of benzene rings and ring C.
    
    Inputs:
        ringidx         ring list containing atom indices
        bridx           indices of benzen rings
        ringA           candidate ring A
        brid            identifier of aromatic ring
        
    Outputs:
        sk, ring indices that make up C6-C3-C6 skeletons
    Element in sk: a list contains ring indices making up the skeleton, indics of:
        .. C in ring A bonding to C in ring C
        .. this C in ring C
        .. C in ring C bonding to ring B
        .. C in ring B bonding to this C
        .. O in ring C bonding to ring A
        .. O of carboxyl group bonding to ring B if existed
        .. identifier for flavonoid (i.e., 'f')
    """
    sk = []

    # pre-declaration to faster the use of global variable and function 
    iters = miters
    fneis, fatoms = neis, atoms

    # set up structure for storing information of rings and ring C
    nr = len(ringidx)
    brix = [i for i in xrange(nr) if aromaticName[i]]
    setrg, lr, clr, rno, i = [], [], [], [], -1
    for ix in ringidx:
        i += 1
        setrg.append(set(ix))
        lr.append(len(ix))
        rno.append(len([j for j in ix if fatoms[j][1]=='o']))
        if len(ix)==6 and rno[-1]==1:
            clr.append(i)

    # if no ring C is found, this is not a flavonoid, return empty list
    if not clr:
        return [], []

    # set up ring Bs to include the identification of homoisoflavonoids
    # i.e., a (phenyl)methyl structure, which will be used as a benzyl as for
    # identification of flavonoids or isoflavonoids
    rgb = brix[:]
    for i in brix:
        for j in ringidx[i]:
            c = [v[0] for v in fneis[j]
                 if v[2]=='c' and v[0] not in ringidx[i] and v[4]==CHAIN]
            for k in c:
                rgb.append((i,j,k))

    # get skeleton
    for rA, j in iters(ringA, rgb):
        i = rA[0]               # candidate ring A
        t = rA[1:]
        jx = j if isinstance(j, int) else j[0]      # candidate ring B
        for k in clr:
            if not setrg[k].issuperset(t) or len(setrg[i]&setrg[k])!=2 or\
               setrg[jx]&setrg[k]:
                continue

            cs = [ix for ix in setrg[k] if ix not in setrg[i] and ix != rA[-1]]
            # remove C in C=O in ring C
            o = [ni[0] for ni in fneis[t[1]] if ni[1]==2 and ni[2]=='o']
            if o:
                cs.remove(t[1])

            # flavonoids or isoflavonoids
            if isinstance(j, int) and j != i:
                for ic in cs:
                    ix = [v[0] for v in fneis[ic] if v[0] in setrg[j] and v[2]=='c']
                    if ix:
                        ix = ix[0]
                        skt = [[i, k, j], t[0],t[1],ic,ix,t[2],t[3]]+o
                        skt.append('f' if aromaticName[j] in 'bq' else 'fd')
                        sk.append(skt)
                        
                        # check whether this is a flavonolignan, which a lignan
                        # group locates adjacent ring B (i.e., j)
                        g = False
                        for p in xrange(nr):
                            if lr[p]!=6 or p==k or len(setrg[p]&setrg[j])!=2:
                                continue

                            u1, u2 = tuple(setrg[p]&setrg[j])
                                
                            o1 = [v for v in fneis[u1] if v[0] in setrg[p] and v[2]=='o']
                            o2 = [v for v in fneis[u2] if v[0] in setrg[p] and v[2]=='o']
                            if o1 and o2:
                                g = True
                                break

                        if g and not any(ic in rg and ix in rg for rg in setrg):
                            sk[-1][0].append(p)
                            sk[-1][-1] = 'fl'

            # homoisoflavonoids         
            elif isinstance(j, tuple) and jx != i:
                jx, j0, j1 = j
                for ic in cs:
                    ix = [v for v in fneis[ic] if v[0]==j1 and v[0] not in setrg[k]]
                    if ix:
                        ix = ix[0][0]
                        sk.append([[i,k,jx], t[0],t[1],ic,ix,j0,t[2],t[3]]+o+['hf'])

    # chech whether multiple skeletons share same benzene rings; these skeletons
    # are possibly coumestans, pterocarpans or peltogynoid;
    # and check whether there exists a ring among ring C and ring B;
    # remove structures that can not be assigned by known class name;
    # if the chain in a flavonoid (i.e., the bond between ring B and C) is
    # one of the edges of another ring but the flavonoid is also not a
    # known class having four rings, remove it.
    if ring7:
        setrg2 = setrg+[set(rg) for rg in ring7]
        ridx2 = ringidx+ring7
        rid = aromaticName+[None]*len(ring7)
        nr2 = len(ridx2)
        for ix in ring7:
            rno.append(len([j for j in ix if fatoms[j][1]=='o']))
            lr.append(len(ix))
    else:
        setrg2, ridx2, nr2, rid = setrg, ringidx, nr, aromaticName

    taidx, n, names = [ra[0] for ra in ringA], len(sk), [None]*len(sk)
    if n>0:
        sk4 = []
        for i in xrange(n):
            if sk[i][-1] != 'f': continue
            skr = sk[i][0]
            
            for j in xrange(nr2):
                if len(setrg2[j]&setrg2[skr[1]])==2 and len(setrg2[j]&setrg2[skr[2]])==2:
                    if rno[j]==1 and not setrg2[j]&setrg2[skr[0]] and\
                       not any(setrg2[j]&setrg2[k] and not k in skr and k!=j
                               for k in xrange(nr2)):
                        skt = sk[i][:]
                        skt[0].append(j)
                        name = classFlav(ridx2, skt, rid)
                        if name:
                            sk[i][:] = skt[:]
                            sk4.append(i)
                            names[i] = name
                    elif rno[j]==0 and lr[j]==5:
                        kt, kc = sk[i][2], setrg[j]&setrg[skr[1]]
                        for k in kc:
                            o = [ni[0] for ni in fneis[k] if ni[2] == 'o']
                            t = any(ni[0] in ringidx[i] for ni in fneis[k])
                            if len(o) == 1 and not t and (kt!=k and kt in kc):
                                sk[i][0].append(j)
                                sk[i] = sk[i][:-1]+[o[0], 'arl']     # 4-arylchroman
                                names[i] = classFlav(ringidx, sk[i], aromaticName)
                                sk4.append(i)

        setsk = [set(ski[0]) for ski in sk]
        for i in xrange(n):
            lr = len(sk[i][0])
            if lr>4:
                if names[i]: names[i] = None
                continue

            if lr<4 and sk4:
                # if other skeletons share any ring with skeleton composed
                # by 4 rings, ignore this skeleton
                sk4c = [j for j in sk4 if setsk[i]&setsk[j]]
                if sk4c:
                    continue

            if lr == 4:
                if sk[i][-1] == 'fl':
                    names[i] = classFlav(ringidx, sk[i], aromaticName)
                    if n > 1:
                        # if this is flavonolignan, part of its rings is also
                        # a skeleton
                        skt = [sk[i][0][:-1]]+sk[i][1:]
                        skt[-1] = 'f'
                        name = classFlav(ringidx, skt, aromaticName)
                        if name:
                            names.append(name)
                            sk.append(skt)
                
                # check whether a third ring exists between ring A and C or
                # ring B and ring between ring B and C, if so, remove this
                # skeleton
                r1, r2, r3, r4 = tuple(sk[i][0])
                for k in xrange(nr):
                    if k not in sk[i][0]:
                        if len(setrg[r2]&setrg[k])>0 and len(setrg[r3]&setrg[k])>0:
                            names[i] = None
                            break
                    
            else:
                # bond topology between ring B and C
                tp = [ni[4] for ni in fneis[sk[i][3]] if ni[0]==sk[i][4]][0]
                if tp!=CHAIN: continue
                
                if sk[i][-1] != 'f':
                    # if this is not standard flavonoid (e.g., flavone, with
                    # name notation of 'f') and share same ring with standard
                    # flavonoid, ignore it.
                    if all(not setsk[i]&setsk[j] for j in xrange(n) if sk[j][-1]=='f'):
                        names[i] = classFlav(ringidx, sk[i], aromaticName)
                else:
                    names[i] = classFlav(ringidx, sk[i], aromaticName)

        n = len(sk)
        sk[:], names[:]= [sk[i] for i in xrange(n) if names[i]], [m for m in names if m]

    return sk, names


def checkChalSK(ringidx, aromaticName):
    """
    Check the skeleton of chain flavonoid, in which the ring C is a chain, not a
    ring. This type of flavonoids is mainly composed by chalcone, dihydrochalcone
    and alpha-methyldeoxybenzoin.

    Return sk and the class of chalcone, i.e., chalcone, dihydrochalcone or
    aurone.
    Elements in sk: a list contains ring indices making up the skeleton, indices of:
        .. C in ring A making up the chain between ring A and ring C
        .. C of the chain bonding to ring A and C(alpha), denoted as C(1)
        .. C(alpha)
        .. C(beta)
        .. C in ring C bonding to C(beta), i.e., C in ring C makeing up the chain
        .. O in C(1)=O
        .. O in C(beta)=O, for oxodihydrochalcone only
        .. identifier, 'c' for chalcone and dihydrochalcone, and 'a' for aurone
    """
    # pre-declaration to faster the use of global variable and function 
    fneis = neis
    iters = miters
    
    sk, name, nameid, skid =[], [], '', []
    setrg,nr,lr =[set(rg) for rg in ringidx],len(ringidx),[len(rg) for rg in ringidx]

    # get neighbor Cs of a benzene ring
    brneis = []
    for i in xrange(nr):
        tnei=[]
        if aromaticName[i]:
            for j in ringidx[i]:
                tnei+=[(n,j) for n in fneis[j] if n[2]=='c' and n[0] not in ringidx[i]]
        brneis.append(tnei)

    sbr = []
    for i, j in iters(xrange(nr), xrange(nr), flag=True):
        
        if not brneis[i] or not brneis[j] or any(si.issuperset([i, j]) for si in sbr):
            continue

        t = False
        for iinfo, jinfo in iters(brneis[i], brneis[j]):
            ni, i0 = iinfo      # neighbors and corresponding atom in the ring
            nj, j0 = jinfo
            i1, j1 = ni[0], nj[0]   # indices of the neighbor atoms
            
            if i1 == j1 or setrg[i]&setrg[j]: continue

            for ni2, nj2 in iters(fneis[i1], fneis[j1]):
                if ni2[0]==nj2[0] and ni2[2]=='c' and ni2[0]!=i0 and nj2[0]!=j0:
                    t = True
                    break
            if t: break

        if t:

            # ketone at chalcone's chain or five-membered ring of aurone
            noxy1 = [v[0] for v in fneis[i1] if v[2]=='o' and v[1]==2]
            noxy2 = [v[0] for v in fneis[j1] if v[2]=='o' and v[1]==2]

            if noxy1 or noxy2:
                # aurone
                nr5, o = [], []
                for k in xrange(nr):
                    if lr[k]==5 and ni2[0] in setrg[k]:
                        cks = []
                        if len(setrg[k]&setrg[i])==2 and noxy1:
                            nr5.append([i, k, j])
                            cks, bj = setrg[k]&setrg[i], nj[-1]
                            tpj = nj2[4]==CHAIN and nj[4]==CHAIN
                        elif len(setrg[k]&setrg[j])==2 and noxy2:
                            nr5.append([j, k, i])
                            cks, bj = setrg[k]&setrg[j], ni[-1]
                            tpj = ni2[4]==CHAIN and ni[4]==CHAIN
                        for ik in cks:
                            oc=[v for v in fneis[ik] if v[2]=='o' and v[0] in setrg[k]]
                            if oc: o.append((oc[0][0], ik))
                        
                if o and len(nr5)==1 and tpj:
                    o, ik = o[0]
                    nameid = 'au' if bj == 2 else 'da'
                    name.append(classFlav([],nameid,[]))
                    nr5, i2 = nr5[0], ni2[0]
                    if nr5[0]==i:
                        sk.append([nr5, i0, i1, i2, j1, j0, ik, o, noxy1[0], 'au'])
                    else:
                        sk.append([nr5, j0, j1, i2, i1, i0, ik, o, noxy2[0], 'au'])
                    skid.append(False)
                    sbr.append(set(nr5))

            if nameid: continue

            # other types of chalcone
            tname = []
            if noxy1 and not noxy2:
                if aromaticName[i] == 'b':
                    nameid = 'c' if nj2[1] >= 2 else 'dc'
                elif aromaticName[i] in 'mq2':
                    nameid = 'cq'
                if nameid:
                    sk.append([[i,j], i0, i1, ni2[0], j1, j0, noxy1[0], 'c'])
        
            elif noxy2 and not noxy1:
                if aromaticName[j] == 'b':
                    namid = 'c' if ni2[1] >= 2 else 'dc'
                elif aromaticName[j] in 'mq2':
                    nameid = 'cq'
                if nameid:
                    sk.append([[j,i], j0, j1, nj2[0], i1, i0, noxy2[0], 'c'])
                    
            elif noxy1 and noxy2:
                sk.append([[i,j], i0, i1, ni2[0], j1, j0, noxy1[0], noxy2[0], 'o'])
                nameid = 'oc'
                
            elif ni2[1]==1 and nj2[1]==1:
                noxy1 = [v[0] for v in fneis[i1] if v[2]=='o' and v[1]==1]
                noxy2 = [v[0] for v in fneis[j1] if v[2]=='o' and v[1]==1]
                if noxy1:
                    sk.append([[i,j], i0, i1, ni2[0], j1, j0, noxy1[0], 'c'])
                    nameid = 'pc'
                elif noxy2:
                    sk.append([[j,i], j0, j1, nj2[0], i1, i0, noxy2[0], 'c'])
                    nameid = 'pc'
            
            elif lr[i] == 6 and aromaticName[i]=='q2':
                if ni[1]==2 and aromaticName[j]=='b' and nj2[1]==2:
                    noxy = [v[0] for v in fneis[i1] if v[2]=='o' and v[1]==1]
                    if noxy:
                        sk.append([[i,j], i0, i1, ni2[0], j1, j0, noxy[0], 'c'])
                        nameid = 'cq'
            
            elif lr[i] == 5:
                noxy = [v[0] for v in fneis[i1] if v[2]=='o' and v[1]==1]
                if noxy:
                    o2 = []
                    for k in ringidx[i]:
                        xt = [ni[0] for ni in fneis[k] if ni[2]=='o' and ni[1]==2]
                        if xt:  o2.append(xt[0])
                    sk.append([[i,j], i0, i1, ni2[0], j1, j0, noxy[0]]+o2+['c'])
                    nameid = 'cc'

            if nameid:
                name.append(classFlav([], nameid, []))
                sbr.append(set([i, j]))
                nameid = ''

                # for these types of flavonoids except aurone and,
                # dihydroaurone, no any atom in the chain between the
                # two benzenes exists in other rings.
                tps, t  = [i1, ni2[0], j1], True
                for k in tps:
                    if any(ni[4]!=CHAIN for ni in fneis[k]):
                        t = False
                        break
                if t:
                    # for alpha and beta carbons, only -OH or -OCH3 is allowed
                    ji, js = sk[-1][3:5], sk[-1][1:-1]
                    for j in ji:
                        nj = [v[0] for v in fneis[j] if v[2]!='o' and v[0] not in js]
                        njo = [v[0] for v in fneis[j] if v[2]=='o' and v[0] not in js]
                        if nj:
                            t = False
                            break
                        elif njo:
                            oj = njo[0]
                            if len(fneis[oj])>1:
                                nj2 = [v[0] for v in fneis[oj] if v[2]=='c' and v[0]!=j]
                                if len(nj2) != 1 or len(fneis[nj2[0]]) != 1:
                                    t = False
                                    break
                
                skid.append(t)

    return sk, name, skid


def checkChalpolymer(csk, cskid, name, sk, ringidx):
    """
    check chalcone dimers and oligomers and validate chalcone
    """
    # preallocation and parameter setting
    fatom, fnei = atoms, neis
    namex = ['chalcone', 'dihydrochalcone', 'pentahydrochalcone']
    ncsk, nr = len(csk), len(ringidx)
    cix = [i for i in xrange(ncsk) if name[i] in namex]
    setrg, lenrg = [set(rg) for rg in ringidx], [len(rg) for rg in ringidx]
    
    # get number of double bonds in each ring
    b2s = []
    for i in xrange(nr):
        b = set()
        for k in ringidx[i]:
            b.update(ni[0] for ni in fnei[k] if ni[0] in ringidx[i] and ni[1]==2)
        b2s.append(b)

    # get assigned benzene rings
    brs, brixs = set(), set()
    for ski in csk: brs.update(ski[0])
    for ski in sk: brs.update(ski[0])

    cr, cc = [], []
    for ski in csk:
        cr.append(set(ski[0]))
        cc.append(set(ski[1:-1]))

    dix, sbr = set(), set()
    # since very loose criteria are applied in identifying chalcone,
    # it is probably some flavonoids are also identified as chalcone,
    # if this situation occurs, remove the chalcone
    for i in xrange(ncsk):
        if any(cr[i].issubset(ski[0]) for ski in sk):
            dix.add(i)

    # identify chalcone adducts obtained by biological synthesis reaction
    # such as Diels-Alder reaction.
    for i in xrange(ncsk):
        js = set(csk[i][3:5])
        if i in cix and not cskid[i] and i not in dix:
            
            if cr[i]&sbr:
                if '_a' not in csk[i][-1]: dix.add(i)
                continue

            for k in xrange(nr):
                if k not in brs and js.issubset(setrg[k]):
                    tj = True
                    for ij in cc[i]:
                        if not all(ni[4]==CHAIN for ni in fnei[ij]):
                            tj = False
                            break
                        
                    if not tj:
                        dix.add(i)
                        csk[i][-1] = 'c'
                        break

                    # form six-membered ring
                    # .. connect to other skeleton or form single ring
                    if lenrg[k] == 6 and len(b2s[k]) >= 2:
                        if not all(fatom[j][1]=='c' for j in setrg[k]):
                            dix.add(i)
                            break
                        
                        for j in cix:
                            if j != i and not cr[i]&cr[j] and not cc[i]&cc[j] and\
                               j not in dix:
                                if len(setrg[k]&setrg[csk[j][0][-1]])==2 and \
                                   name[i] in namex[:2] and name[j] in namex[:2]:
                                    sbr.update(csk[i][0]+csk[j][0])
                                    csk[i][-1] += '_a'
                                    csk[j][-1] += '_a'
                                    break
                
                    # form five-membered ring tetrohydro-furan   
                    elif lenrg[k]==5 and sum(fatom[j][1]=='o' for j in ringidx[k])==1:
                        t = False
                        oi = [j for j in ringidx[k] if fatom[j][1]=='o'][0]
                        for j in brs:
                            if len(setrg[j]&setrg[k]) == 2 and j not in csk[i][0]:
                                if any(ni[0] in setrg[j] for ni in fnei[oi]):
                                    sbr.update(cr[i])
                                    sbr.add(j)
                                    csk[i][-1] += '_a'
                                else:
                                    dix.add(i)
                                t = True
                                break
                        if not t:
                            for j in cix:
                                if j != i and len(setrg[k]&cc[j])>=2 and not cc[j]&cc[i]:
                                    t = True
                                    sbr.update(csk[i][0]+csk[j][0])
                                    csk[i][-1] += '_a'
                                    csk[j][-1] += '_a'
                                    break
                        if not t:
                            dix.add(i)

                    # form four-membered ring cyclobutane
                    elif lenrg[k] == 4:
                        if not all(fatom[j][1]=='c' for j in ringidx[k]):
                            dix.add(i)
                            break
                        
                        t = False
                        for j in cix:
                            if j != i and setrg[k].issuperset(csk[j][3:5]):
                                t = True
                                sbr.update(csk[i][0]+csk[j][0])
                                csk[i][-1] += '_a'
                                csk[j][-1] += '_a'
                                break
                        if not t:
                            dix.add(i)
                            break

                elif lenrg[k]==5 and sum(fatom[j][1]=='o' for j in ringidx[k])==1:
                    if len(setrg[k]&cc[i])>=2:
                        t = False
                        for j in brs:
                            if len(setrg[k]&setrg[j])==2 and j not in cr[i]:
                                t = True
                                sbr.update(cr[i])
                                sbr.add(j)
                                csk[i][-1] += '_a'
                                break
                        if not t:
                            dix.add(i)
            
        elif not 'aurone' in name[i] and not cskid[i]:
            t = False
            for j in js:
                nj = [v[0] for v in fnei[j] if v[2]=='c' and v[0] not in cc[i]]
                if len(nj)==1:
                    for ski in sk:
                        if any(j in ringidx[k] for k in ski[0]):
                            t = True
                            break
                if t: break
            if not t:
                dix.add(i)

    for i in dix:
        if '_a' in csk[i][-1]:
            dix.remove(i)

    csk[:] = [csk[i] for i in xrange(ncsk) if i not in dix]
    name[:] = [name[i] for i in xrange(ncsk) if i not in dix]
    for i in xrange(len(csk)):
        if '_a' in csk[i][-1]:
            name[i] = 'chalcone_adduct'

    return csk, name


def checkDARflav(ringidx, aromaticName, ringA, sk, names):
    """
    check skeleton with Diels-Alder reactions, this function mainly
    focuses to chalcone-isoprene adducts
    """
    fneis, fatoms = neis, atoms
    setrg, nr, nsk = [set(rg) for rg in ringidx], len(ringidx), len(sk)
    
    # find two six-member rings having 3 common atoms for identifying
    # chalcone Diels-Alder adducts
    das = []
    for i, j in miters(xrange(nr), xrange(nr), flag=True):
        if len(setrg[i]&setrg[j]) == 3:
            for rA in ringA:
                if setrg[j].issuperset(rA[1:]) and aromaticName[rA[0]]=='b'\
                   and all(fatoms[k][1]=='c' for k in ringidx[i]):
                    das.append([i, j, rA[0], (rA[1:3])])

    for rs in das:
        t = False
        if any(rs[2] in ski[0] for ski in sk):
            for i in xrange(nsk):
                if sk[i][-1]=='c' and setrg[rs[0]].issuperset(sk[i][3:5]):
                    sk[i][-1] += '_da'
                    names[i] = 'chalcone_da_adduct'
                    t = True

        if not t:    
            for rA in ringA:
                i = rA[0]
                for j in xrange(nr):
                    if i not in rs and j not in rs and setrg[j].issuperset(rA[1:])\
                       and len(setrg[j]&setrg[rs[0]])==2:
                        j0, j1 = rs[-1]
                        i2 = [ni[0] for ni in fneis[rA[2]]
                              if ni[0] in setrg[j] and ni[0] not in setrg[i]][0]
                        ox = [ni[0] for ni in fneis[rA[2]] if ni[1]==2 and ni[2]=='o']
                        if ox:
                            sk.append([[i,j],rA[1],rA[2],i2,j1,j0,ox[0],'c_da'])
                            names.append('chalcone_da_adduct')

    # form six-membered ring
    for i in xrange(nsk):
        if sk[i][-1]=='c':
            js, t = set(sk[i][3:5]), True
            for j in js:
                if any(ni[1]==2 for ni in fneis[j]):
                    t = False
                    break
            if not t:
                continue
            
            for k in xrange(nr):
                if js.issubset(setrg[k]) and len(setrg[k])==6 and aromaticName[k]!='b':
                    b2s, t = set(), True
                    for j in ringidx[k]:
                        b2s.update(ni[0] for ni in fneis[j]
                                   if ni[0] in setrg[k] and ni[1]==2)
                        if fatoms[j][1]!='c':
                            t = False
                            break
                        
                    if not any(setrg[k]&setrg[j] for j in xrange(nr) if j != k)\
                       and len(b2s)==2 and t:
                        sk[i][-1] += '_da'
                        names[i] = 'chalcone_da_adduct'

    return sk, names


def checkAnthoSK(ringidx, brix, oix):
    """
    Check the skeleton of anthocyanidins
    """
    fneis, fatoms = neis, atoms

    sk, names, bs = [], [], []
    nrg, nnei, setrg = len(ringidx), len(fneis), [set(rg) for rg in ringidx]
    setrg, rno = [], []
    for ix in ringidx:
        setrg.append(set(ix))
        rno.append(len([i for i in ix if fatoms[i][1]=='o']))

    for o, k in oix:
        
        if rno[k] != 1: continue
        
        ij = [None, None]
        for i in brix:
            t = False
            for ni in fneis[o]:
                
                if ni[1] == 2: t=True
                
                if ni[0] in ringidx[i] and len(setrg[i]&setrg[k]) == 2:
                    # ring A, use string to avoid false identification when i = 0
                    ij[0] = str(i)
                elif t:
                    for nj in fneis[ni[0]]:
                        jx = nj[0]
                        if jx in ringidx[i] and nj[4] == CHAIN:
                            ij[1] = (i, ni[0], jx)      # ring B
                            break

        if all(ij):
            j, j1, j0 = ij[1]
            c1 = [ni[0] for ni in fneis[j1] if ni[0]!=o and ni[0] in ringidx[k]][0]
            c2 = [ni[0] for ni in fneis[c1] if ni[0] in ringidx[k] and ni[1]==2]
            if c2:
                i = int(ij[0])
                ik = list(setrg[i]&setrg[k])
                ix = 0 if [ni[0] for ni in fneis[ik[0]] if ni[2]=='o'] else 1
                i0, i1 = ik[1-ix], ik[ix]                   # C(C)-C(O) of ring A
                c1 = [ni[0] for ni in fneis[i1] if ni[0] in ringidx[k] and ni[0]!=i0][0]
                sk.append([[i, k, j], i0, c1, j1, j0, i1, o, 'a'])
                names.append('anthocyanidin')
                bs += [i, j]

    return sk, names


def checkXanthSK(ringidx, bridx):
    """
    Get skeleton of Xanthones.
    Elements in sk: a list of indices of rings making up the skeleton, indices of:
        .. C in benzene A bonding to oxygen O
        .. O
        .. C in benzene B bonding to oxygen O
        .. C in benzene A bonding to carboxy C
        .. carboxy C
        .. C in benzene B bonding to carboxy C
    """
    fneis, fatom = neis, atoms
    iters = miters

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


def checkSpecialSK(ringidx, aromaticName, ringA):
    """ Get skeleton for special flavonoids
    Elements in sk: a list of indices of rings making up the skeleton, indices of:
        .. C in ring A bonding to ring C (for 'b') or chain (for 'm')
        .. C in ring C ring C (for 'b') or chain (for 'm') bonding to ring A,
           denoted as C(1)
        .. C in ring C (for 'b') or chain (for 'm') bonding to ring B, denoted as C(2)
        .. C in ring B bonding to ring C (for 'b') or chain (for 'm')
        .. O in C(1)=O (for 'm' only)
        .. C in C(2)-CH3
        .. identifier, 'b' for 2-arylbenzofuran and 'm' for alpha-methyldeoxybenzoin
        ringA.append([i, j, k, nei[0], bo[0]])
    """
    # pre-declaration to faster the use of global variable and function 
    fneis, fatoms = neis, atoms
    iters = miters

    sk = []

    nr = len(ringidx)
    if nr >= 3:
        brs = [i for i in xrange(nr) if aromaticName[i]=='b']
        if len(brs)>=2:
            setrg, rno, k = [], [], []
            for i in xrange(nr):
                if i not in brs and len(ringidx[i]) == 5:
                    k.append(i)
                    setrg.append(set(ringidx[i]))
                    rno.append(len([j for j in ringidx[i] if fatoms[j][1]=='o'])==1)
            nr = len(k)
            
            for rA in ringA:
                r, r0, r1, r2, r3 = rA[0], rA[1], rA[2], rA[3], rA[4]
                rgc = set(ringidx[r])
                for i in xrange(nr):
                    if len(rgc&setrg[i])==2 and r3 in setrg[i] and rno[i]:
                        # C(alpha)
                        c = [n[0] for n in fneis[r1] if n[0] in setrg[i] and n[1]==2]
                        if not c: continue
                        c = c[0]
                        for j in brs:
                            if j==r: continue
                            c1 = [ni[0] for ni in fneis[c]
                                  if ni[0] in ringidx[j] and ni[4]==CHAIN]
                            if c1:
                                c1 = c1[0]
                                sk.append([[r,k[i],j], r0,r1,c,c1,r2,r3, 'b'])
                                break

    for i,j in iters(xrange(nr), xrange(nr), flag=True):
        
        if not aromaticName[i] or not aromaticName[j]:
            continue
        
        for i0, j0 in iters(ringidx[i],ringidx[j]):
            for nei1, nei2 in iters(fneis[i0], fneis[j0]):
                if nei1[4]==CHAIN and nei2[4]==CHAIN:
                    i1, j1 = nei1[0], nei2[0]
                    if i1 != j1:
                        o = [ni[0] for ni in fneis[i1]
                             if ni[2]=='o' and ni[1]==2 and ni[4]==CHAIN]
                        if not o: continue
                        c = [ni[0] for ni in fneis[i1]
                             if ni[0]==j1 and ni[2]=='c' and ni[4]==CHAIN]
                        if not c: continue
                        c1 = [ni[0] for ni in fneis[j1] if ni[2]=='c' and ni[1]==1
                              and ni[0]!=j0 and ni[4]==CHAIN and ni[0]!=i1]
                        for ci in c1:
                            if len(fneis[ci]) == 1:
                                sk.append([[i, j], i0, i1, j1, j0, o[0], ci, 'm'])
                                break
                    elif aromaticName[i]=='q':
                        c = [ni[0] for ni in fneis[i1]
                             if ni[2]=='c' and ni[0]!=i0 and ni[0]!=j0 and ni[4]==CHAIN]
                        if c:
                            c = c[0]
                            c1 = [ni[0] for ni in fneis[c] if ni[0]!=c
                                   and ni[1]==2 and ni[4]==CHAIN and ni[2]=='c']
                            if c1:
                                sk.append([[i, j], i0, i1, j0, c, c1[0], 'd'])

    return sk


def checkValidSKs(ringidx, sk, names):
    """
    Check skeletons to identify whether those skeletons share same rings or
    bonds, if so, maximum independent set algorithm is used to get skeleton set
    with most number of independent skeletons (i.e., no sharing bonds between
    any two skeletons).
    Besides, a priori check is performed for chalcones to ensure valid
    groups around three-membered chains.
    Output: a list with indices retained
    """
    fnei = neis

    delix = set()
    skrs, n = [set(ski[0]) for ski in sk], len(sk)
    sets = [set(rg) for rg in ringidx]
    
    # remove chalcones in which atoms in the three-membered
    # chain is part of other rings, or any other skeleton sharing same
    # ring or atom with chalcone_adduct
    # for pentahydrochalcones, very strict criterion is set, at which
    # not any other group links to -OH group
    for i in xrange(n):
        if sk[i][-1]=='c' or sk[i][-1]=='o':
            if names[i]=='pentahydrochalcone':
                if len(fnei[sk[i][-2]])>1:
                    delix.add(i)
                else:
                    for k in sk[i][2:5]:
                        if any(ni[0] not in sk[i][1:] or ni[4]!=CHAIN for ni in fnei[k]):
                            delix.add(i)
                            break
                continue
            
            for k in sk[i][2:5]:
                if any(ni[4]!=CHAIN for ni in fnei[k]):
                    delix.add(i)
                    break
        elif '_a' in sk[i][-1] or '_da' in sk[i][-1]:
            delix.update(j for j in xrange(n) if skrs[i]&skrs[j] and j!=i)

    # if they share same atoms or rings, select one
    c = lambda sk1, sk2, s: set(sk1[0])&set(sk2[0]) or\
        any(set(sk1[1:])&s[k] for k in sk2[0]) or\
        any(set(sk2[1:])&s[k] for k in sk1[0])

    for i in xrange(n):
        if i in delix: continue
        t = [c(sk[i],sk[j],sets) for j in xrange(n) if j!=i and j not in delix]
        if t and all(t):
            delix.add(i)
    
    retains = [i for i in xrange(n) if i not in delix]

    # get maximum independent set
    if len(retains)>1:
        s = mis([skrs[i] for i in retains])
        retains = [retains[v[1]] for v in s]

    return retains


def getSkeleton(mol):
    """
    Get skeleton and class name of the input molecule described by string such as
    canonical SMILES and InchI or molecular file.
    """
    na = numrgass
    atominfo(mol)
    fatom, fneis = atoms, neis
    
    sk, skix, names =[], [], []
    
    ringidx, aromaticName = getRings(mol)
    
    if not ringidx or not aromaticName or not any(aromaticName):
        return None, None, None, None
    validringdistcheck(ringidx, aromaticName)
    
    nr = len(ringidx)
    brix = [i for i in xrange(nr) if aromaticName[i]]
    n = len(brix)

    # tetrahydroflavanones
    pseudbr, j = [None]*nr, -1
    for r in ringidx:
        j += 1
        if not aromaticName[j] and all(fatom[i][1]=='c' for i in r) and len(r)==6:
            b = []
            for i in r:
                b += [ni[1]==2 for ni in fneis[i] if ni[0] in r]
            if sum(b)==2:
                pseudbr[j] = 'b'

    if any(pseudbr):
        pseudbr2 = [pseudbr[i] or aromaticName[i] for i in xrange(nr)]
        ringAx = getRingA(ringidx, pseudbr)
        fsk, name = checkFlavSK(ringidx, pseudbr2, ringAx, [])
        if len(fsk)==1 and aromaticName[fsk[0][0][2]] and pseudbr[fsk[0][0][0]]:
            if name[0] == 'flavanone':
                skix = getskidx(ringidx, fsk)
                return fsk,skix,'tetrahydroflavanones',(ringidx,aromaticName)
    
    if n<2: return None, None, None, None

    # anthocyanidin
    oix, t = [], False
    for atom in fatom:
        if atom[1] == 'o' and atom[2] == 1:
            t = True
            ri = [i for i in xrange(nr) if atom[0] in ringidx[i]]
            if ri:
                oix.append((atom[0], ri[0]))
        elif atom[2] != 0:
            raise FlavonoidException(6)

    if oix:
        ask, name = checkAnthoSK(ringidx, brix, oix)
        if ask:
            sk[:], names[:] = ask, name
            n -= na(sk, aromaticName)
            if n < 2:
                skix = getskidx(ringidx, sk)
                return sk, skix, names[0], (ringidx, aromaticName)
            # remove assigned benzene rings
            rmix, brid = [], aromaticName[:]
            for ski in sk:  rmix += ski[0]
            for i in brix:
                if i not in rmix:
                    brid[i] = None
        else:
            if t: return None, None, None, None

    if 'brid' not in locals():
        brid = aromaticName

    ringA = getRingA(ringidx,brid)
    benzenes = [i for i in xrange(nr) if brid[i]=='b']

    # flavonoids
    # .. get indices of seven-membered rings
    ring7 = []
    for r in mol.iterateRings(7,7):
        r7ix = [atm.index() for atm in r.iterateAtoms()]
        if not any(set(r7ix).issuperset(rg) for rg in ringidx):
            ring7.append(r7ix)
    brid7 = [None]*len(ring7)

    # combine rings with 7 elements
    ringidx2, brid2 = ringidx+ring7, brid+brid7

    # flavonoid and its derivatives
    flavsk, name = checkFlavSK(ringidx, brid, ringA, ring7)
    if flavsk:
        # stilbeno-flavonoid
        flavsk, name = checkStilbenoSK(ringidx, benzenes, flavsk, name)
        sk[:] = [flavsk[i] for i in xrange(len(flavsk)) if name[i]]
        names[:] = [name[i] for i in xrange(len(flavsk)) if name[i]]
        # if this is stilbeno-flavonoid, return the identified skeletons directly
        nf = len(flavsk)
        if nf>1:
            dix = []
            j = [sk[i][-1]=='sf' for i in xrange(nf)]
            if j:
                j = j[0]
                for i in xrange(nf):
                    if i != j and set(sk[i][0])&set(sk[j][0]):
                        dix.append(i)
            if dix:
                sk[:] = [sk[i] for i in xrange(nf) if i not in dix]
                names[:] = [names[i] for i in xrange(nf) if i not in dix]
                skix = getskidx(ringidx, sk)
                return sk, skix, names, (ringidx, aromaticName)

        # Theaflavin
        if ring7 and 'flavan-3-ol' in names:
            sk, name = checkTheaflavinSK(sk, ringidx, brid, ringA, ring7, names)
            if isinstance(name, str):
                skix = getskidx(ringidx2, sk)
                return sk, skix, name, (ringidx2, aromaticName+brid7)

    # chalcone
    csk, cname, cskid = checkChalSK(ringidx, brid)
    if csk:
        ncsk = len(csk)
        if ncsk>1 or (sk and ncsk):
            csk, cname = checkChalpolymer(csk, cskid, cname, sk, ringidx2)
        sk = sk+csk
        names += cname
    sk, names = checkDARflav(ringidx2, brid2, ringA, sk, names)
    if sk:
        retains = checkValidSKs(ringidx, sk, names)
        sk[:] = [sk[i] for i in retains]
        names[:] = [names[i] for i in retains]

    # special flavonoids
    for ski in sk:
        for i in ski[0]:
            if i in benzenes:
                benzenes.remove(i)
                brid2[i] = None
    if len(benzenes)>=1:
        xsk = checkXanthSK(ringidx, benzenes)
        osk = checkSpecialSK(ringidx, brid, ringA)
        osk += xsk
        for ski in osk:
            if not any(set(ski[0])&set(skj[0]) for skj in sk):
                name = classFlav(ringidx, ski, brid)
                names.append(name)
                sk.append(ski)

    if not sk:
        sk, names = checkBiflavonoidSK(ringidx, benzenes, ringA)

    skix = getskidx(ringidx2, sk)
    
    return sk, skix, names, (ringidx2, aromaticName+brid7)


"""
Indexing the skeleton to get the site of side chains. This procedure follows
the IUPAC provisional recommendation for Nomenclature of Flavonoids
(Rauter A. P., et al. 2009). Note that this function only considers the side
groups that bond to skeleton and saccharides.
"""

def indexglc(glcrings, ringidx):
    """ index glycosides of flavonoids"""
    fnei, fatom = neis, atoms
    glcix = []

    for glcring in glcrings:
        glcixi = []
        rix = [rg for rg in ringidx if set(glcring).issuperset(rg)][0]
        cix = set([i for i in glcring if fatom[i][1]=='c'])
        epc = cix.difference(rix)
        
        # candidate anomeric carbon
        for i in cix:
            if any(ni[2]=='o' and ni[1]==2 for ni in fnei[i]):
                c = [ni[0] for ni in fnei[i] if ni[2]=='c' and ni[0] in cix][0]
                if sum(ni[2]=='o' for ni in fnei[c])==2 and\
                   any(ni[2]=='o' and ni[0] in rix for ni in fnei[c]):
                    glcixi += [i, c]
                    break
                
            if any(ni[2]=='o' and ni[0] in rix for ni in fnei[i]) and not\
               any(ni[0] in epc for ni in fnei[i]):
                if sum(ni[2]=='o' for ni in fnei[i])>=2:
                    glcixi.append(i)
                    break
                elif len(fnei[i])>2 and not glcixi:
                    # in this situation, terminal O is bounded to the first C
                    # and OH at this C carbon is replaced by other group
                    glcixi.append(i)
                    break

        while True:
            flag = True
            for i in cix:
                if i not in glcixi and any(ni[0]==glcixi[-1] for ni in fnei[i]):
                    glcixi.append(i)
                    flag = False
                    break

            if flag:
                break
        glcix.append(glcixi)

    return glcix


def indexring(ring, stinfo, sgx, flag):
    """
    index ring according to start atom information (stinfo) and
    flag. "flag" is an indicator to indicate the type of ring, the
    valid input is A, B or C.
    "sgx" is the indices of all atoms in side chains
    """
    fneis = neis
    
    sti, l = stinfo[-1], len(ring)-1
    if flag != 'b':
        
        ix1, ix2 = ring.index(stinfo[0]), ring.index(stinfo[1])
        if flag=='a':
            if ((ix2==0 and ix1==l) or (ix2>ix1 and not (ix2==l and ix1==0))):
                ring.reverse()
        elif flag=='c':
            if ((ix1==0 and ix2==l) or (ix2<ix1 and not (ix1==l and ix2==0))):
                ring.reverse()
        ix = ring.index(sti)
        ring[:] = ring[ix:] + ring[:ix]
        
    elif flag == 'b':

        ix = ring.index(sti)
        ring[:] = ring[ix:] + ring[:ix]
        nc1, nc2 = [], []
        for i in ring:
            if any(ni[0] in sgx and ni[0] not in ring for ni in fneis[i]) and i!=sti:
                ix2 = ring.index(i)+1
                nc1.append(ix2)
                nc2.append(len(ring)-ix2+2)

        if sum(nc2)<sum(nc1):
            ring.reverse()
            ring[:] = [ring[-1]]+ring[:-1]

    return ring


def indexsg(sk, skix, sg, ringidx, names):
    """ indexing side chains """
    fneis = neis

    skgix, sgix, skatoms, skatomorder = [], [], [], []
    
    setrg, lr, lsk = [set(rg) for rg in ringidx], [len(rg) for rg in ringidx], len(sk)

    sgs = []
    for key in sg.keys():
        for g in sg[key]['atomIndex']:
            sgs += g
    if not sgs and len(sk)==1:
        return skgix, sgix
    
    exs = ['coumestan', 'pterocarpan', 'dihydropterocarpan']
    for ii in xrange(lsk):
        skt, ixt, sgsi, ixx = [], [], [], []
        ski = sk[ii]
        
        sgsi[:] = sgs[:]
        for i in xrange(lsk):
            if i != ii:
                sgsi += skix[i]
                
        l, skgixi, sgixi = len(ski[0]), [], []
        if (l == 3 and ski[-1] != 'x') or (l == 4 and names[ii] not in exs):
            
            i, k, j = tuple(ski[0][:3])
            
            # get indices of atoms between rings A and C
            if len(ski[1:])==7 or (len(ski[1:])==8 and ski[-1] in 'farlflfd'):
                i0, i1, o, j0 = ski[1], ski[5], ski[6], ski[4]
            else:
                i0, i1, o, j0 = ski[1], ski[6], ski[7], ski[5]
            
            # return the rank of the atoms
            # ring A
            c = [x[0] for x in fneis[i0] if x[0] not in setrg[k] and x[0] in setrg[i]]
            rgt = indexring(ringidx[i][:], [i0,i1,c[0]], [], 'a')
            skt += rgt
            ixt += [str(i) for i in xrange(lr[k]-1,lr[k]+3)]+[None, None]
            ixx += [i for i in xrange(lr[k]-1, lr[k]+3)]+[9,10]

            # ring C
            rgt = indexring(ringidx[k][:], [i0,i1,o], [], 'c')
            skt += rgt
            ixt += [str(i+1) for i in xrange(lr[k]-2)]+[None, None]
            ixx += [i+1 for i in xrange(lr[k]-2)]+[5,6]

            # ring B
            rgt = indexring(ringidx[j][:], [None,None,j0], sgsi, 'b')
            skt += rgt
            ixt += [str(i+1)+"'" for i in xrange(lr[j])]
            ixx += [i+1 for i in xrange(lr[j])]
            
        elif l == 2:

            i, j = ski[0][0], ski[0][1]

            # ring A
            rgt = indexring(ringidx[i][:], [None,None,ski[1]], sgsi, 'b')
            skt += rgt
            ixt += [str(i+1)+"'" for i in xrange(lr[i])]
            ixx += [i+1 for i in xrange(lr[i])]

            # ring C
            j0 = [k for k in ski[1:] if k in ringidx[j]][0]
            rgt = indexring(ringidx[j][:], [None,None,j0], sgsi, 'b')
            skt += rgt
            ixt += [str(i+1) for i in xrange(lr[j])]
            ixx += [i+1 for i in xrange(lr[j])]

        elif l == 4:            # pterocarpan and coumestan

            ks = ski[0]
            # get start atom of each ring to find the atom that participates
            # in the indexing of atoms
            ipr = [(ski[1],ski[5]), (ski[5],ski[1])]
            ci = [i for i in setrg[ks[2]]&setrg[ks[3]] if i!=ski[4]][0]
            ipr += [(ski[4],ci), (ci,ski[4])]

            s = []
            for i in ks:
                s.append([x for x in setrg[i]
                          if not any(x in setrg[j] for j in ks if j!=i)])

            ixs = []
            for i in xrange(l):
                i0, i1 = ipr[i]
                while True:
                    xi = [ni[0] for ni in fneis[i0] if ni[0] in s[i] and ni[0]!=i1]
                    if xi:
                        ixs += xi
                        i0, i1 = xi[0], i0
                    else:
                        break
            skt += ixs
            ixt += [str(i+1) for i in xrange(len(ixs))]
            ixx += [i+1 for i in xrange(len(ixs))]

        else:                   # xanthones

            i, j, k = ski[0][0], ski[0][2], ski[0][1]
            i1, i2, j1, j2 = ski[1], ski[4], ski[3], ski[6]
            c1 = [x[0] for x in fneis[i1] if x[0] not in setrg[k]][0]
            c2 = [x[0] for x in fneis[j2] if x[0] not in setrg[k]][0]
            
            rg1 = indexring(ringidx[i][:], [i1,i2,c1], [], 'a')
            rg2 = indexring(ringidx[j][:], [j1,j2,c2], [], 'c')
            # combine atoms in ring i and j not in ring k as a large ring
            rg = rg2[:4]+rg1[:4]
            rg = indexring(rg, [None, None, c2], sgsi, 'b')
            if not any(ni[0]==rg[1] for ni in fneis[rg[0]]):
                rg[:] = rg[1:]+[rg[0]]

            skt += rg
            ixt += [str(i+1) for i in xrange(len(rg))]
            ixx += [i+1 for i in xrange(len(rg))]

        skatoms.append(skt)
        skatomorder.append(ixx)

        for i in xrange(len(skt)):
            if ixt[i]:
                for key in sg.keys():
                    j = -1
                    for g in sg[key]['atomIndex']:
                        j += 1
                        if skt[i] in g or any(ni[0] in g and ni[0] not in skt for ni in fneis[skt[i]]):
                            skgixi.append(ixt[i])
                            sgixi.append((key, j, skt[i]))

                for ij in xrange(lsk):
                    if ij != ii:
                        if any(ni[0] in skix[ij] for ni in fneis[skt[i]]):
                            skgixi.append(ixt[i])
                            sgixi.append(('sk', ij, skt[i]))
        
        skgix.append(skgixi)
        sgix.append(sgixi)

    return skgix, sgix, (skatoms, skatomorder)
 

"""
Get the side chains of current molecule. List of side chains can be found in
side_chain module. If new side chains are found, it can be added to the end
of the list.

Note:
    It should be noted that, the list is orderred as ring contained groups
    followed by non-ring groups which is in order of most complexity (i.e.,
    contain more atoms) to simple groups (i.e., have less atoms). Thus, if
    new side chain is added, it should be added according this order.

Modify:
    To avoid complexity in adding new groups according to order of group list,
    side chain list is stored in a dictionary with keys are 'chainGroup' and
    'ringGroup':
            >>> import side_chains as sd
            >>> sd.sideChains.keys()
            ['chainGroup', 'ringGroup']
    Each element in 'ringGroup' list is a dictionary with keys: 'smiles',
    'num', 'name', 'sinfo', 'ckr':
            >>> sd.sideChains['ringGroup'][0].keys()
            ['smiles', 'num', 'name', 'sinfo', 'ckr']
    which store canonical SMILES, number of atoms, name, structure information
    (number of atoms in chain, number of rings) of the group and a check tag
    to indicate whether matching the group of which one atom is removed.
    Each element in 'chainGroup' list is also a dictionary with keys 'smiles',
    'num', 'name':
            >>> sd.sideChains['chainGroup'][0].keys()
            ['smiles', 'num', 'name']
    which store canonical SMILES, number of atoms and name of the group.

"""

def getsidenets(sk, skix, sginfo, glcnumbering):
    """ get networks of side groups """
    sgps, skgix, sgix = sginfo
    fneis = neis
    if not sgix:
        if not sgps:
            return 'There is not any side group around the skeleton.\n'
        # if the class is exceptional, statistics of side groups are printed only.
        else:
            sx = ''
            for key in sgps.keys():
                gname = sgps[key]['groupName']
                if not sgps[key]['groupName']: continue
                sx += 'Chain groups include: ' if 'chain' in key \
                      else 'Ring groups include: '
                unique_gname = set(sgps[key]['groupName'])
                kx = []
                for name in unique_gname:
                    nn = sum(name==x for x in sgps[key]['groupName'])
                    kx.append('%d %s groups'%(nn, name))
                sx += ','.join(kx)
                sx += '.\n'
            if not sx:
                return 'There is not any side group around the skeleton.\n'
            return sx.rstrip()

    lsk, skixset = len(sk), [set(ix) for ix in skix]
    glcix = sginfo[0]['ringGroup']['glcIndex']
    
    sdneis, sdsets, sdkey, sdnum = {}, {}, {}, {}
    for key in sgps.keys():
        if not sgps[key]['atomIndex']: continue
        kx = 'ring' if 'ring' in key else 'chain'

        cns = []            # neighbors of all side groups
        for gi in sgps[key]['atomIndex']:
            t = set()
            for j in gi:
                t.update(ni[0] for ni in fneis[j] if ni[0] not in gi)
            cns.append(t)
        sdneis[key] = cns
        sdsets[key] = [set(gi) for gi in sgps[key]['atomIndex']]
        sdkey[key] = kx
        sdnum[key] = len(sgps[key]['groupName'])
    
    s = ''
    for key in sgps.keys():
        if not sgps[key]['atomIndex']: continue

        s += key+':\n'
        n = len(sgps[key]['groupName'])
        for i in xrange(n):
            s += 'side group %s (%d) connects to ' %(sgps[key]['groupName'][i],i)
            t = ''
            for j in xrange(lsk):
                if skixset[j]&sdneis[key][i]:
                    t += ' %d,'%j
            if t:
                s += 'skeletons%s ' %t if len(t)>3 else 'skeleton%s ' %t

            for key2 in sgps.keys():
                kx, t, sc = sdkey[key2], '', ''
                for j in xrange(sdnum[key2]):
                    if not (j==i and key2==key) and sdsets[key][i]&sdneis[key2][j]:
                        # identify the positions of glycosides linking to other glycosides
                        if 'ring' in key and kx=='ring' and i in glcix and j in glcix:
                            il, jl = glcix.index(i), glcix.index(j)
                            c = list(sdsets[key][i]&sdneis[key2][j])[0]
                            o = sdsets[key][i]&sdsets[key2][j]
                            if o:
                                xi = glcnumbering[il].index(c)+1
                                o=list(o)[0]
                                oi=[ni[0] for ni in fneis[o] if ni[0] in sdsets[key][j]]
                                oi = oi[0]
                                xj = glcnumbering[jl].index(oi)+1
                                t += ' [(%d%s%d)-O-] ring group %d,'%(xi,u'\u2192',xj,j)
                            else:
                                e = 'C'
                                ci=[ni[0] for ni in fneis[c] if ni[0] in sdsets[key][i]]
                                z = [ni[2]=='c' for ni in fneis[ci[0]] if ni[0]==c][0]
                                if z:
                                    xi = glcnumbering[il].index(c)+1
                                else:
                                    e = 'O'
                                    xi = glcnumbering[il].index(ci[0])+1
                                z=[(ni[2]=='c', ni[0]) for ni in fneis[c]
                                   if ni[0] in sdsets[key][j]][0]
                                if z[0]:
                                    ci = z[1]
                                else:
                                    e = 'O'
                                    ci = [ni[0] for ni in fneis[z[1]] if ni[0]!=c][0]
                                xj = glcnumbering[jl].index(ci)+1
                                t+=' [(%d%s%d)-%s-] ring group %d,'%(xi,u'\u2192',xj,e,j)
                        else:
                            t += ' %d,'%j
                if t:
                    s += '%s groups%s' %(kx,t) if len(t)>3 else '%s group%s' %(kx,t)
            
            s = s[:-1] + '\n'
    return s


def getsidegroupdist(sk, names, sginfo):
    """ get distributions of side groups around skeleton """
    sgps, skgix, sgix = sginfo
    if not sgix:
        if not sgps:
            return 'There is not any side group around the skeleton.\n'
        # if the class is exceptional, statistics of side groups are printed.
        else:
            sx = ''
            for key in sgps.keys():
                gname = sgps[key]['groupName']
                if not sgps[key]['groupName']: continue
                sx += 'Chain groups are: ' if 'chain' in key else 'Ring groups are: '
                unique_gname = set(sgps[key]['groupName'])
                kx = []
                for name in unique_gname:
                    nn = sum(name==x for x in sgps[key]['groupName'])
                    kx.append('%d %s groups'%(nn, name))
                sx += ','.join(kx)
                sx += '.\n'
            if not sx:
                return 'There is not any side group around the skeleton.\n'
            return sx.rstrip()

    fatoms, fneis = atoms, neis
    
    lsg, glcix = [len(g) for g in sgix], sgps['ringGroup']['glcIndex']
    s, i  = '', -1
    for ski in sk:
        i += 1
        sc, kix, j = [], [], -1
        for gi in sgix[i]:
            j += 1
            gk, ik, gmi = gi
            if gk == 'sk':
                kx = [(x,y,z) for x,y,z in sgix[ik] if x=='sk' and y==i][0]
                k = sgix[ik].index(kx)
                si = "(%s%s%s) skeleton %d"%(skgix[i][j],u'\u2192',skgix[ik][k],ik+1)
                if not kix:
                    sc.insert(0, si)
                    kix.append(ik)
                else:
                    if ik > max(kix):
                        sc.insert(len(kix), si)
                        kix.append(ik)
                    for xi in xrange(len(kix)):
                        if ik < kix[xi]:
                            sc.insert(xi, si)
                            kix.insert(xi, ik)
            else:
                gpname = sgps[gk]['groupName'][ik]
                if 'ring' in gk:
                    for ix in sgps[gk]['atomIndex'][ik]:
                        if any(ni[0] == gmi for ni in fneis[ix]):
                            if ik in glcix:
                                namepre = 'O-' if fatoms[ix][1]=='o' else 'C-'
                            elif fatoms[ix][1]=='o':
                                namepre = 'O-'
                            gpname = '%s%s'%(namepre,gpname)
                            break
                si = "%s-%s (%d)"%(skgix[i][j], gpname, ik)
                sc.append(si)
        if sc:
            s += 'skeleton %d (%s): '%(i+1, names[i]) + '; '.join(sc)
    return s


def isaromatic(brs):
    """
    check whether there exist large aromatic systems, which can dramatically
    slow down the speed of finding sub groups containing aromatic system
    """
    fneis = neis
    n = len(brs)
    brneis = []
    for r in brs:
        nx = []
        for i in r:
            nx += [ni[0] for ni in fneis[i] if ni[0] not in r and ni[2]=='c']
        brneis.append(set(nx))
        
    for i in xrange(n):
        a = [i]
        while True:
            t = 0
            for j in xrange(n):
                if j not in a:
                    for k in a:
                        if len(brs[j]&brs[k])>0 or len(brneis[j]&brs[k])>0:
                            a.append(j)
                            t += 1
                            break
                if len(a) >= 3:
                    return True
            if t == 0:
                break
    return False


def isvalidringmatch(mhi, rgset, qix, mol, qinfo, checkgly=False):
    """
    check whether the match of ring group is valid
    """
    # check whether any of ring in rgset is subset of matched group
    # if not, return False
    if not any(rg.issubset(mhi) for rg in rgset):
        return False

    # check whether the chain is part of a ring that is not assigned,
    # if so, return False, thus the match is rejected
    fneis = neis
    _, chainix, _, qsml, isbr, num = qinfo
    chainset = set([mhi[i] for i in xrange(len(qix)) if qix[i] in chainix])
    
    for rg in rgset:
        cmr = chainset&rg
        for j in cmr:
            if any(ni[0] in rg and ni[0] in mhi for ni in fneis[j]):
                return False
    
    if checkgly: return True
    
    # check the exact match of the ring group to avoid false assignment
    # due to the structure match resulted from tautomer match, if not passed
    # the exact match, return False
    q = idg.loadMolecule(qsml)
    m = mol.createSubmolecule(mhi)
    #.. if contains benzene ring, check the ring to identify whether
    #.. the extracted sub molecule still contains enough double bonds
    #.. if not and this is a benzene ring, set the lost double bond
    #.. to double
    if isbr:
        for r in m.iterateRings(5,6):
            sml = r.clone().smiles()
            if sml.count('c')==6 or sml.count('=')==3: continue
            ix = None
            aix = [a.index() for a in r.iterateAtoms()]
            for a in r.iterateAtoms():
                if all(ni.bond().bondOrder()==1 and ni.index() in aix
                       for ni in a.iterateNeighbors()) and a.countHydrogens()>0:
                    if not ix:
                        ix = a.index()
                    else:
                        for ni in a.iterateNeighbors():
                            if ni.index()==ix and ni.countHydrogens()>0:
                                ni.bond().setBondOrder(2)
                                ix = None
                    
    if not idg.exactMatch(m, q):
        if not idg.exactMatch(m, q, 'TAU'):
            return False
    
    return True
    

def getsmallchains(ix, igix):
    """
    get small chain groups containing only single or two atoms that bond
    to the input atom index 'ix'
    """
    fneis = neis

    gps = []
    
    nei1 = [nei[0] for nei in fneis[ix] if nei[0] not in igix]
    for j in nei1:
        
        if len(fneis[j])==1:              # this means that the only neighbor is ix
            gps.append(nei1)
            continue
        
        nei2 = [nei[0] for nei in fneis[j] if nei[0] != ix]
        if len(nei2)==1 and len(fneis[j])==2 and len(fneis[nei2[0]])==1:
            gps.append(nei1+nei2)

    return gps


def separategroups(mol):
    """ remove atom indices in 'rmix' and get separated side groups,
        note that atoms in molecule object input (i.e., mol) here have
        been removed by function "removeatoms"
    """
    fnei = neis
    
    sdix = [atom.index() for atom in mol.iterateAtoms()]
    # group retained indices according to side chains
    sdix.sort()
    sidegroups  = [[i] for i in sdix]                            # initial groups
    nsd         = len(sdix)
    
    neiset, sdset = [], []
    for sd in sidegroups:
        sdset.append(set(sd))
        neiset.append(set([x[0] for i in sd for x in fnei[i]]))

    while True:
        t = True
        for i in xrange(nsd-1):
            for j in xrange(i+1,nsd):
                if len(sdset[i]&neiset[j])>0 and not sdset[i].issuperset(sdset[j]):
                    sdset[i].update(sdset[j])
                    neiset[i].update(neiset[j])
                    t = False
        if t:
            break

    # remove combined side groups
    for i in xrange(nsd-1):
        for j in xrange(i+1,nsd):
            if sdset[i].issuperset(sdset[j]):
                sdset[j].clear()
    sidegroups[:] = [list(sd) for sd in sdset if sd]
                
    return sidegroups


def removeatoms(mol, setig, rgset, rmrginfo, qr):
    """ remove atoms that have been assigned """
    fnei = neis
    getAtom = mol.getAtom
    allatomset = set([atom.index() for atom in mol.iterateAtoms()])
    retainset = allatomset.difference(setig)
    setig.intersection_update(allatomset)
    
    rmrgset, brset = rmrginfo
    rmix = []
    for i in setig:
        
        if any(i in rg for rg in rgset): continue
        
        if not any(nei[0] in retainset for nei in fnei[i]) or i in rmrgset:
            getAtom(i).remove()
            rmix.append(i)
        else:
            g = getsmallchains(i, setig)
            cnei = [nei for nei in fnei[i] if nei[0] not in setig]
            if len(g)==len(cnei):
                getAtom(i).remove()
                rmix.append(i)

    # if an edge of a ring bonds to a benzene ring that has been assigned and
    # bond order of this edge is 1, set it to 2
    if qr:
        mol.dearomatize()
        
        if not mol.dearomatize():
            idg.setOption('dearomatize-verification',False)
            mol.dearomatize()
            idg.setOption('dearomatize-verification',True)

        for rg in rgset:
            for br in brset:
                b = list(rg&br)
                if len(b)==2:
                    bo = [ni[1] for ni in fnei[b[0]] if ni[0]==b[1]][0]
                    t = any(ni[1]==2 or ni[1]==4 for k in b for ni in fnei[k]
                            if ni[0] not in b and ni[0] not in rmix)
                    if all(getAtom(i).countHydrogens()>0 for i in b) and not t and\
                       (bo==1 or bo==4):
                        for nei in getAtom(b[0]).iterateNeighbors():
                            if nei.index() == b[1]:
                                nei.bond().setBondOrder(2)
##        idgrender.renderToFile(mol,'m2.png')


def getQueryGroup(ringidx, name, ckr=True):
    """
    Check side group
    ringidx     indices of rings out of skeleton rings
    asname      assigned names
    ckr         identifier to indicate the type of side chain, True for ring group
    """
    # declare local variables from global variables and functions
    loadQueryMolecule   = idg.loadQueryMolecule
    key                 = 'ringGroup' if ckr else 'chainGroup'
    sidechains          = sideChains[key]
    qnames              = [g['name'] for g in sidechains]
    
    querysds, qidxs     = [], []
    if ckr:
        
        lr = [len(rg) for rg in ringidx]
        cx, i, nameix = [], -1, [None]*len(name)
        for g in sidechains:
            i += 1
            gp = g['property']
            n = gp['number']

            if name and gp['checkaromatic']:
                for j in xrange(len(name)):
                    if g['name'] == name[j]:
                        nameix[j] = i
                continue
            
            q = loadQueryMolecule(g['smiles'])
            crix = getRings(q, checkBenzene = False)
            if gp['checkglc']:
                cx[:] = [j for j in xrange(n) if j not in crix[0]]
                delatm = []
                for j in cx:
                    catm = q.getAtom(j)
                    if catm.symbol() == 'C':
                        delatm.append(j)
                    elif catm.symbol() == 'O' and any(nei.bond().bondOrder()==2 for
                                                      nei in catm.iterateNeighbors()):
                        delatm.append(j)
                for j in delatm:
                    cx.remove(j)
                qidxs.append([cx[:],i])
            else:
                if len(crix) == 1:
                    cx[:] = [j for j in xrange(n) if j not in crix[0]]
                else:
                    rx = []
                    for r in crix: rx += r
                    cx[:] = [j for j in xrange(n) if j not in rx]
            
            querysds.append((i, set(cx), q, g['smiles'], gp['haveBenzene'], n))

        if nameix:
            nameix[:] = nameix if any(nameix) else []
            return querysds, qidxs, qnames, nameix

        return querysds, qidxs, qnames
        
    else:
        
        for g in sidechains:
            q = loadQueryMolecule(g['smiles'])
            querysds.append((q, g['property']['number'], g['smiles']))
    
        return querysds, qnames


def aromaticRingMatch(benzenerings, rmix):
    """
    match side chains having aromatic ring system
    """
    fneis = neis

    gp, name, gpc, namec = [], [], [], []
    
    for rg in benzenerings:
        gpc[:], namec[:] = [], []
        for i in rg:
            for nei1 in fneis[i]:
                if nei1[2] == 'c' and nei1[4]==CHAIN:
                    c1 = nei1[0]
                    c2 = [ni for ni in fneis[c1] if ni[2]=='c'
                          and ni[0] not in rmix and ni[4]==CHAIN]
                    if not c2: continue

                    for k0 in c2:
                        if k0[1]!=1 or\
                           sum(ni[1]==1 and ni[2]=='c' for ni in fneis[k0[0]])<3:
                            continue
                        kk0 = k0[0]
                        c3 = [ni[0] for ni in fneis[kk0] if ni[2]=='c' and ni[0]!=c1
                              and ni[1]==1 and ni[4]==CHAIN]
                        k1, k2, k3 = None, None, []
                        for k in c3:
                            if sum(ni[1]==1 and ni[2]=='c' and ni[4]==CHAIN
                                   for ni in fneis[k])<3:
                                k1 = k
                            elif len(k3)<2:
                                k3[:] = [ni[0] for ni in fneis[k] if ni[0]!=kk0
                                      and ni[2]=='c' and ni[1]==1]
                                k2 = k
                                
                            if k1 and len(k3)>=2:
                                ix = (rg.index(i)+3)%6
                                o3 = [ni[0] for ni in fneis[rg[ix]] if ni[2]=='o'
                                      and ni[4]==CHAIN and ni[1]==1]
                                if o3:
                                    gpc.append([c1, kk0, k1, k2, k3[0], k3[1], o3[0]])
                                    namec.append('2,3-dimethylbutylphenol')
     
                    c2 = [x[0] for x in c2 if x[1]==2]
                    if not c2: continue
                    c2 = c2[0]
                    c3 = [ni[0] for ni in fneis[c2] if ni[2]=='c' and ni[0]!=c1
                          and ni[0] not in rmix and ni[4]==CHAIN]
                    if not c3:
                        gpc.append([c1, c2])
                        namec.append('phenylethylene')
                        continue
                    
                    for k in c3:
                        o1 = [ni[0] for ni in fneis[k] if ni[2]=='o' and ni[1]==1
                             and ni[4]==CHAIN]
                        o2 = [ni[0] for ni in fneis[k] if ni[2]=='o' and ni[1]==2
                             and ni[4]==CHAIN]
                        if o1 and o2:
                            o1, o2 = o1[0], o2[0]
                            oc = [(rg.index(i)+len(rg)+2+j)%6 for j in xrange(3)]
                            o3 = [ni[0] for ni in fneis[oc[1]] if ni[2]=='o'
                                  and ni[4]==CHAIN and ni[1]==1]
                            if o3:
                                oh = [ni[0] for ni in fneis[oc[0]] if ni[2]=='o'
                                      and ni[4]==CHAIN and ni[1]==1]
                                if oh:
                                    gpc.append([c1,c2,k,o1,o2,o3[0],oh[0]])
                                    namec.append('caffeic acid')
                                else:
                                    oh = [ni[0] for ni in fneis[oc[2]] if ni[2]=='o'
                                          and ni[4]==CHAIN and ni[1]==1]
                                    if oh:
                                        gpc.append([c1,c2,k,o1,o2,o3[0],oh[0]])
                                        namec.append('caffeic acid')
                                    else:
                                        gpc.append([c1,c2,k,o1,o2,o3[0]])
                                        namec.append('p-Coumaric acid')
                            else:
                                gpc.append([c1, c2, k, o1, o2])
                                namec.append('cinnamic acid')
                            
                        elif o1:
                            gpc.append([c1, c2, k, o1[0]])
                            namec.append('cinnamyl alcohol')
                        else:
                            gpc.append([c1, c2, k])
                            namec.append('hydroxystyryl')
        if gpc:
            l = [len(g) for g in gpc]
            i = l.index(max(l))
            gp.append(rg+gpc[i])
            name.append(namec[i])

    return gp, name


def glycosideMatch(mol, querys, qidxs, rgset, sdset):
    """ match glycoside groups """
    fnei, rmatoms       = neis, removeatoms
    g, gix = [], []
    # prior check whether candiate glycosyl rings exist
    begly = []
    for rg in rgset:
        # get number of O in and beside the ring
        oinr, ober, b2 = set(), set(), set()
        for i in rg:
            for ni in fnei[i]:
                if ni[0] in rg:
                    if ni[1] >= 2:
                        b2.add(i)
                    if ni[2] == 'o':
                        oinr.add(ni[0])
                else:
                    if ni[1]==1 and ni[2]=='o':
                        ober.add(i)
                    elif ni[1]>=2:
                        b2.add(i)
        if len(oinr)==1 and len(ober)>=2 and not b2:
            begly.append(rg)

    if not begly:
        return g, gix
            
    substructureMatcher = idg.substructureMatcher
    loadQueryMolecule   = idg.loadQueryMolecule
    exactMatch          = idg.exactMatch
    loadMolecule        = idg.loadMolecule
    iis = [qi[1] for qi in qidxs]
    hasChiral = mol.isChiral()
    # indices of last query group in all groups that have same number of atoms
    terminix = []
    for i in set([qi[-1] for qi in querys]):
        terminix.append([qi[0] for qi in querys if qi[-1]==i][-1])
    #idgrender.renderToFile(mol,'m2.png')
    
    # match groups
    rmrg, delatms, cg, cgix, lp, chiralTag = set(), set(), [], [], 0, []
    for qc in querys:
        i, _, q, qsml, _, num = qc
        ij = iis.index(i)

        qams  =  [atom for atom in q.iterateAtoms()]
        qamix = [atm.index() for atm in qams]
        for j in xrange(len(qidxs[ij][0])+1):
            atmc, atmic = list(qams), list(qamix)
            # further checking by delete an atom around the ring
            q2 = loadMolecule(qsml)
            if j:
                q = loadQueryMolecule(qsml)
                atmic.pop(qidxs[ij][0][j-1])
                atmc.pop(qidxs[ij][0][j-1])
                for ni in qams[qidxs[ij][0][j-1]].iterateNeighbors():
                    ix = atmic.index(ni.index())
                q.getAtom(qidxs[ij][0][j-1]).remove()
                q2.getAtom(qidxs[ij][0][j-1]).remove()
            nqstereo = q2.countStereocenters()

            if substructureMatcher(mol).countMatches(q):
                
                matcher = substructureMatcher(mol)
                t = False
            
                for mhi in matcher.iterateMatches(q):
                    mhc = [mhi.mapAtom(atm).index() for atm in atmc]
                    if isvalidringmatch(mhc, rgset, atmic, mol, qc, True):
                        # further check, they should be matched exactly
                        # by stereostructure
                        submol = mol.createSubmolecule(mhc)
                        nmstereo = submol.countStereocenters()
                        if (nqstereo and nmstereo and exactMatch(submol,q2,'STE')) or\
                            nqstereo<nmstereo or (nqstereo==0 and nmstereo==0):
                            # the deleted atom should be replaced by other groups
                            if (j and not all(ni[0] in mhc for ni in fnei[mhc[ix]]))\
                               or not j:
                                cg.append(mhc)
                                cgix.append(i)
                                chiralTag.append(q.isChiral())
                                lp, t = num, True

                if t and j: break

        # in the above iterations, all possible matches are obtained
        # to select possibly correct matches, following criteria are
        # applied:
        # 1. match obtaining maximum atoms is retained;
        # 2. if several maximums are found:
        # .... matches having chiral strucutre are preferential,
        # .... then, lowest index in 'gix' is retained.
        if cg and i in terminix:
            lg = [len(gi) for gi in cg]
            selix = []
            for rg in rgset:
                ixs, ncg = [], len(cg)
                for j in xrange(ncg):
                    if rg.issubset(cg[j]):
                        ixs.append(j)

                nc = len(ixs)
                if nc == 0: continue

                # add current ring into the set which contains all
                # atoms should be removed in next iteration of query
                # group match.
                rmrg.update(rg)
                
                if nc==1:
                    selix.append(ixs[0])
                    continue

                # check whether any query group has chiral structure,
                # if yes, retain the structure
                ix2 = [j for j in ixs if chiralTag[j]]
                if ix2: ixs = ix2

                # if multiple structure matched with same length
                # minimum query group index is adopted
                ml = max(lg[j] for j in ixs)
                ix2 = [j for j in ixs if lg[j]==ml]
                if len(ix2)>1:
                    gix2 = [cgix[j] for j in ix2]
                    selix.append(ix2[gix2.index(min(gix2))])
                else:
                    selix.append(ix2[0])

            g += [cg[j] for j in selix]
            gix += [cgix[j] for j in selix]
            for gi in g:
                delatms.update(gi)
                
                # update side groups
                for k in xrange(len(sdset)):
                    if sdset[k].issuperset(gi):
                        sdset[k].difference_update(gi)
                        for j in gi:
                            if any(ni[0] in sdset[k] for ni in fnei[j]):
                                sdset[k].add(j)

            # remove atoms from the molecule and stored ring groups
            rgset[:] = [rg for rg in rgset if not rg.issubset(delatms)]
            rmatoms(mol, delatms, rgset, (rmrg, []), False)

            if not rgset:
                return g, gix
            else:
                t = False
                for rg in rgset:
                    if any(rg==bi for bi in begly):
                        t = True
                        break
                if not t:
                    return g, gix
            
            cg, cgix, chiralTag = [], [], []

    return g, gix


def ringGroupMatch(mol, rgsets, querys, qidxs):
    """
    Find ring groups from current molecule
    """
    matchedix, mgix  = [], []
    if mol.countSSSR()==0:
        return matchedix, mgix
    # pre-declare local functions and variables
    loadQueryMolecule   = idg.loadQueryMolecule
    substructureMatcher = idg.substructureMatcher
    rmatoms             = removeatoms
    fnei                = neis

    # preallocation and initialization
    rgset, brset= rgsets
    # identifier for identifying whether the ring is a subset of another object
    brsubgroupid = [rg.issubset for rg in brset]

    # get separated side groups
    sdgroups = separategroups(mol)
    sdset = [set(g) for g in sdgroups]

    # get ring groups and properties such as number of atoms in a subgroup
    # and whether there exist benzene rings
    havebr, rix = False, []
    for g in sdset:
        if any(rg.issubset(g) for rg in rgset):
            rix.append(g)
        if not havebr and any(subi(g) for subi in brsubgroupid):
            havebr = True
    sdset[:] = rix
    msd, nl = max(len(g) for g in sdset), len(sdset)

    # match glycosides
    iis = [qi[1] for qi in qidxs]
    glyqs = [qi for qi in querys if qi[0] in iis]
    if any(qi[-1]<=msd+1 for qi in glyqs):
        matchedix, mgix = glycosideMatch(mol, glyqs, qidxs, rgset, sdset)
    if not rgset:
        return matchedix, mgix
    msd, nl = max(len(sdi) for sdi in sdset), len(sdset)

    # match groups
    rmrg, delatms, delrg, ismatched = set(), set(), [], False
    while True:
        nm = 0
        for qc in querys:
            i, _, q, qsml, bi, num = qc
            # if length of ring group smaller than, or number of rings
            # larger than that of the query, or no benzene ring exists in
            # side chain while the query group has, continue
            if havebr-bi or num>msd or 'N' in qsml or i in iis: continue
            
            qatoms  =  [atom for atom in q.iterateAtoms()]
            qatomix = [atm.index() for atm in qatoms]

            if substructureMatcher(mol).countMatches(q):
                delatms.clear()
                matcher = substructureMatcher(mol)
                
                for mhi in matcher.iterateMatches(q):
                    mhc = [mhi.mapAtom(atm).index() for atm in qatoms]
                    
                    if isvalidringmatch(mhc, rgset, qatomix, mol, qc):
                        delatms.update(mhc)
                        matchedix.append(mhc)
                        mgix.append(i)
                        ismatched = True
                        nm += 1

                        for k in xrange(nl):
                            if sdset[k].issuperset(mhc):
                                sdset[k].difference_update(mhc)
                                for j in mhc:
                                    if any(ni[0] in sdset[k] for ni in fnei[j]):
                                        sdset[k].add(j)

                # remove atoms from the molecule and stored ring groups
                if delatms:
                    rmrg.clear()
                    delrg[:] = []
                    for rg in rgset:
                        if rg.issubset(delatms):
                            rmrg.update(rg)
                            delrg.append(rg)
                    rgset[:] = [rg for rg in rgset if rg not in delrg]
                    rmatoms(mol, delatms, rgset, (rmrg, brset), True)
                    
                    if not rgset:
                        if mol.countSSSR()>0:
                            delrg[:] = []
                            for r in mol.iterateRings(3,8):
                                sr = set(a.index() for a in r.iterateAtoms())
                                if any(sr.issubset(g) for g in matchedix):
                                    delrg += list(sr)
                            rmatoms(mol, set(delrg), [], ([], []), False)
                        return matchedix, mgix
                    
            # check tautomer matching
            while True:
                mh = substructureMatcher(mol,'TAU').match(q)
                if not mh: break
                c = [mh.mapAtom(atm).index() for atm in qatoms]

                if not isvalidringmatch(c, rgset, qatomix, mol, qc):
                    break

                rmrg.clear()
                delrg[:] = []
                for rg in rgset:
                    if rg.issubset(c):
                        rmrg.update(rg)
                        delrg.append(rg)
                if not delrg: break

                matchedix.append(c)
                mgix.append(i)
                ismatched = True
                nm += 1

                # remove atoms from the molecule and stored ring groups
                rgset[:] = [rg for rg in rgset if rg not in delrg]
                rmatoms(mol, set(c), rgset, (rmrg, brset), True)
                if not rgset:
                    if mol.countSSSR()>0:
                        delrg[:] = []
                        for r in mol.iterateRings(3,8):
                            sr = set(a.index() for a in r.iterateAtoms())
                            if any(sr.issubset(g) for g in matchedix):
                                delrg += list(sr)
                        rmatoms(mol, set(delrg), [], ([], []), False)
                    return matchedix, mgix

                for k in xrange(nl):
                    if sdset[k].issuperset(c):
                        sdset[k].difference_update(c)
                        for j in c:
                            if any(ni[0] in sdset[k] for ni in fnei[j]):
                                sdset[k].add(j)

            if ismatched:
                if havebr:
                    t = False
                    for g in rgset:
                        if any(subi(g) for subi in brsubgroupid):
                            t = True
                            break
                    havebr = t

                # get the max length of side groups containing rings
##                idgrender.renderToFile(mol,'m%d.png'%i)
                rix[:] = []
                for rg in rgset:
                    t = True
                    for i in xrange(nl):
                        if sdset[i].intersection(rg):
                            t = False
                            sdset[i].update(rg)
                            break
                    if t:
                        rix.append(rg)

                sdset += rix
                sdset[:] = [g for g in sdset if any(g.issuperset(rg) for rg in rgset)]
                msd, nl = max(len(g) for g in sdset), len(sdset)
                ismatched = False
        
        if nm==0:
            break

    return matchedix, mgix


def alkaloidMatch(mol, querys, rgset):
    """ match alkaloid group """
    fatom, fnei = atoms, neis
    substructureMatcher = idg.substructureMatcher
    rmatoms             = removeatoms

    g, gix = [], []

    if rgset:
        t = False
        for rg in rgset:
            if any(fatom[i][1]=='n' for i in rg):
                t = True
                break
        if not t:
            return g, gix

        rmrg, delrg = set(), []
        for qc in querys:
            i, chainix, q, qsml, bi, num = qc

            qatoms  =  [atom for atom in q.iterateAtoms()]
            while True:
                mh = substructureMatcher(mol,'TAU').match(q)
                if not mh: break
                c = [mh.mapAtom(atm).index() for atm in qatoms]

                rmrg.clear()
                delrg[:] = []
                for rg in rgset:
                    if rg.issubset(c):
                        rmrg.update(rg)
                        delrg.append(rg)
                if not delrg: break

                g.append(c)
                gix.append(i)

                # remove atoms from the molecule and stored ring groups
                rgset[:] = [rg for rg in rgset if rg not in delrg]
                rmatoms(mol, set(c), rgset, (rmrg, []), False)
                if not rgset:
                    return g, gix
        
        for rg in rgset:
            if any(fatom[i][1]=='n' for i in rg):
                raise FlavonoidException(8)
            
    else:
        if not any(fatom[a.index()][1]=='n' for a in mol.iterateAtoms()):
            return g, gix

        for qc in querys:
            i, q, num, qsml = qc

            qatoms  =  [atom for atom in q.iterateAtoms()]
            while True:
                mh = substructureMatcher(mol,'TAU').match(q)
                if not mh: break
                c = [mh.mapAtom(atm).index() for atm in qatoms]
                
                # check whether C=O is near N, otherwise, exceptions
                # will be probably raised
                t = False
                for j in c:
                    if fatom[j][1]=='c' and any(x[1]==2 and x[2]=='o' for x in fnei[j]):
                        t = True
                        break
                    elif fatom[j][1]=='n':
                        for ni in fnei[j]:
                            if ni[2]=='c' and\
                               any(x[1]==2 and x[2]=='o' for x in fnei[ni[0]]):
                                t = True
                                break
                        if t: break
                        
                if not t:
                    if qsml == 'CNC':
                        # this should be dimethylamino, thus C should be
                        # the terminal methyl group
                        for j in c:
                            if fatom[j][1]=='c' and len(fnei[j]) > 1:
                                raise FlavonoidException(8)
                    elif qsml == 'N':
                        # this should be the terminal NH2
                        if any(len(fnei[j])>1 for j in c):
                            raise FlavonoidException(8)
                
                g.append(c)
                gix.append(i)

                # remove atoms from the molecule
                for j in c: mol.getAtom(j).remove()
                # if no atom is retained, return the matched indices
                if mol.countAtoms() == 0:
                    return g, gix
    
    return g, gix
    

def chainGroupMatch(mol, querys):
    """ Find side chains from molecule """
    matchedix, mgix = [], []
    if mol.countAtoms()==0:
        return matchedix, mgix
    # pre-declare local functions and variables
    substructureMatcher = idg.substructureMatcher
    iterateMatches      = substructureMatcher(mol).iterateMatches
    
    # preallocation and initialization
    n, delatms  = len(querys), set()

    # get separated side groups
    sdgroups = separategroups(mol)
    sdset = [set(g) for g in sdgroups]
    maxm, nl = max(len(g) for g in sdset), len(sdset)

    # match groups
    ismatched = False
    for i in xrange(n):
        q, num, qsml = querys[i]
        if num>maxm or 'N' in qsml: continue
        
        qatoms  =  [atom for atom in q.iterateAtoms()]

        if substructureMatcher(mol).countMatches(q):
            delatms.clear()
            for mhi in substructureMatcher(mol).iterateMatches(q):
                mhc = [mhi.mapAtom(atm).index() for atm in qatoms]

                delatms.update(mhc)
                matchedix.append(mhc)
                mgix.append(i)
                ismatched = True

                for k in xrange(nl):
                    if sdset[k].issuperset(mhc):
                        sdset[k].difference_update(mhc)

            # remove atoms from the molecule
            if delatms:
                for j in delatms: mol.getAtom(j).remove()
                # if no atom is retained, return the matched indices
                if mol.countAtoms() == 0: return matchedix, mgix
        
        # check tautomer matching
        while True:
            mh = substructureMatcher(mol,'TAU').match(q)
            if not mh: break
            c = [mh.mapAtom(atm).index() for atm in qatoms]
            
            matchedix.append(c)
            mgix.append(i)
            ismatched = True

            # remove atoms from the molecule
            for j in c: mol.getAtom(j).remove()
            # if no atom is retained, return the matched indices
            if mol.countAtoms() == 0: return matchedix, mgix

            for k in xrange(nl):
                if sdset[k].issuperset(c):
                    sdset[k].difference_update(c)

        if ismatched:
            sdset[:] = [g for g in sdset if g]
            maxm, nl = max(len(g) for g in sdset), len(sdset)
                                        
    return matchedix, mgix


def getSidechains(mol, skeleton, skix, ringinfo):
    """
    Get side chains of the molecule after specifying the skeletons
    """
    rmatoms = removeatoms
    n = mol.countAtoms()
    rgix = getRings(mol, checkBenzene = False)
    skidx, rmrg, rmrgix = set(), set(), set()
    nr = len(ringinfo[0])
    for sk in skeleton:
        rmrgix.update(sk[0])
        for i in sk[0]:
            rmrg.update(ringinfo[0][i])
    for skixi in skix: skidx.update(skixi)
    rgset = [set(ringinfo[0][i]) for i in xrange(nr) if i not in rmrgix]
    rgset += [set(rg) for rg in rgix if rg not in ringinfo[0]]
    
    # obtain benzene rings
    brs = [ringinfo[0][i] for i in xrange(nr) if ringinfo[1][i] == 'b']

    if len(set(skidx)) == mol.countAtoms():
        return {'ringGroup':{'atomIndex':[], 'groupName':[], 'glcIndex': []},
                'chainGroup':{'atomIndex':[], 'groupName':[]}}, True

    if len(rgix) > len(rmrgix):
        # delete the ring composed by other two or more rings
        delix, nr = [], len(rgset)
        for i in xrange(nr):
            if i in delix: continue
            k = 0
            for j in xrange(nr):
                if j != i and rgset[i].issuperset(rgset[j]):
                    k += 1
            if k > 0: delix.append(i)
        if delix:
            rgset[:] = [rgset[i] for i in xrange(nr) if i not in delix]
            nr -= len(delix)
            
        # if there exists large aromatic system, extract the side chains
        # directly
        brset = [set(rg) for rg in brs]
        sdbrset = [rg for rg in brset if rg in rgset]
        g = []
        if len(sdbrset) >= 3 and isaromatic(sdbrset):
            g, name = aromaticRingMatch(brs, skidx)
            name[:] = name if name else [None]
            querysds, qidxs, qnames, nix = getQueryGroup(rgset, name)
        else:
            querysds, qidxs, qnames = getQueryGroup(rgset, [])

        rgset2 = list(rgset)

        rmatoms(mol, skidx, rgset, (rmrg, brset), True)
        # match alkaloid groups
        querysdsN = []
        for qi in querysds:
            if 'N' in qi[3]:
                querysdsN.append(qi)
        ng, ngix = alkaloidMatch(mol, querysdsN, rgset)
##        idgrender.renderToFile(mol,'m2.png')
        mhidx, mhgs = ringGroupMatch(mol, (rgset, brset), querysds, qidxs)
        
        if g:
            # combine the pre-assigned groups to current identified groups
            # if the ring assigned previously also is identified, compare the
            # length of two groups, retain the group having more atoms
            setgs = [set(gi) for gi in g]
            setmhg = [set(gi) for gi in mhidx]
            delix1, delix2 = [], []
            for rg in rgset2:
                t1, t2 = False, False
                for gi in setgs:
                    if len(gi&rg) == len(rg):
                        t1 = True
                        break
                for gj in setmhg:
                    if len(gj&rg) == len(rg):
                        t2 = True
                        break
                if t1 and t2:
                    if len(gi) <= len(gj):
                        delix1.append(gi)
                    else:
                        delix2.append(gj)

            mhgs[:] = [mhgs[i] for i in xrange(len(mhgs)) if setmhg[i] not in delix2]
            mhidx[:] = [mhidx[i] for i in xrange(len(mhidx)) if setmhg[i] not in delix2]
            nix[:] = [nix[i] for i in xrange(len(g)) if setgs[i] not in delix1]
            g[:] = [g[i] for i in xrange(len(g)) if setgs[i] not in delix1]
                
            mhidx += g
            mhgs += nix

        mhidx += ng
        mhgs += ngix

        # get matched saccharides
        glcix = [qx[1] for qx in qidxs]
        mhglc = [i for i in xrange(len(mhgs)) if mhgs[i] in glcix]
        
        mhrnames = [qnames[i] for i in mhgs]
        mhs = set()
        for mhi in mhidx: mhs.update(mhi)
        setig = set([atom.index() for atom in mol.iterateAtoms()])&mhs
        for rg in rgset: rmrg.update(rg)
    else:
        setig = set(skidx)
        mhidx, mhrnames, mhglc = [], [], []

    if mol.countAtoms()>0:
        rmatoms(mol, setig, [], (rmrg, []), False)
        
        if mol.countAtoms() > 0:
            querysds, qnames = getQueryGroup([], [], ckr=False)
            querysdsN = []
            for i in xrange(len(querysds)):
                qi1, qi2, qi3 = querysds[i]
                if 'N' in qi3:
                    querysdsN.append((i, qi1, qi2, qi3))
            ng, ngix = alkaloidMatch(mol, querysdsN, [])

            mhidxc, mhgsc = chainGroupMatch(mol, querysds)
            # delete groups that have been assigned in other groups
            if mhidx:
                delix = []
                mhrset = [set(g) for g in mhidx]
                for gi in mhidxc:
                    if any(len(set(gi)&g)==len(gi) for g in mhrset):
                        delix.append(gi)
                nr = len(mhgsc)
                mhgsc[:]=[mhgsc[i] for i in xrange(nr) if mhidxc[i] not in delix]
                mhidxc[:]=[g for g in mhidxc if g not in delix]
            mhgsc += ngix
            mhidxc += ng
            mhcnames = [qnames[i] for i in mhgsc]
        else:
            mhidxc, mhcnames = [], []
            
    else:
        mhidxc, mhcnames = [], []

    mg = {'ringGroup':  {'atomIndex': mhidx, 'groupName':mhrnames, 'glcIndex': mhglc},
          'chainGroup': {'atomIndex': mhidxc,   'groupName':mhcnames}}

    if mol.countAtoms()>0:
        return mg, False
    
    return mg, True


def getmoleculeInfo(string):
    """ Get information of molecule """
    
    mol = loadmol(string)

    if any(atom.radicalElectrons() for atom in mol.iterateAtoms()):
        raise FlavonoidException(7)
    
    mol.dearomatize()               # convert aromatic structure to Kekule form
    if not validnitrogencheck(mol):
        raise FlavonoidException(8)

    # get main skeleton
    sk, skix, names, ringinfo =  getSkeleton(mol)

    if not sk or not names:
        raise FlavonoidException(4)

    # get side chains
    mg, validsd = getSidechains(mol, sk, skix, ringinfo)
    if not validsd:
        raise FlavonoidException(7)

    # indexing skeletons and glycon rings, however, if the skeleton
    # belongs to some exceptions, skip this
    if (isinstance(names, str) and names in name_exceptions) or\
       any(name in name_exceptions for name in names):
        return sk, skix, names, (mg, [], []), []
    
    skgix, sgix, skorderinfo = indexsg(sk, skix, mg, ringinfo[0], names)
    if mg['ringGroup']['glcIndex']:
        glcs = []
        for i in mg['ringGroup']['glcIndex']:
            glcs.append(mg['ringGroup']['atomIndex'][i])

        glcnumering = indexglc(glcs, ringinfo[0])
    else:
        glcnumering = []

    return sk, skix, names, (mg, skgix, sgix, glcnumering), skorderinfo


if __name__=='__main__':
    s = 'C1C(OC2=C(C1=O)C=CC(=C2)OC3C(C(C(C(O3)CO)O)O)O)C4=CC=CC=C4'
    sk, name, mg, sg = getmoleculeInfo(s)
    print "skeleton:", sk, '\nname:', name, '\nside chains:', mg
