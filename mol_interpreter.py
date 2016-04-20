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
np.dong572@gmail.com
"""

import os, json
from re import findall
from itertools import combinations
from indigo import Indigo, IndigoException
from indigo_inchi import IndigoInchi
from indigo_renderer import IndigoRenderer as renderer

idg         = Indigo()
idgin       = IndigoInchi(idg)
idgrender   = renderer(idg)

CHAIN = idg.CHAIN

with open('sidechain.json','r') as f:
    sideChains = json.load(f)
    

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
            raise ValueError('Can not load molecule from identifier %s'%string)


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
        return '%sflavonoid: '%(numberprefix(len(names)-1))+', '.join(names)
    elif len(names)==1:
        return names[0]
    
    
def miters(iter1,iter2,flag=False):
    """
    iterator for pair iterating of two iterables.
    """
    if flag:
        return iter([(x,y) for x in iter1 for y in iter2 if x!=y])
    else:
        return iter([(x,y) for x in iter1 for y in iter2])


def numrgass(assrings, idx):
    """ Return number of rings assigned """
    t = []
    for info in assrings:
        t += info[0]
    return len(set([i for i in t if i in idx]))


def getskidx(ringidx, sk):
    """ get the indices of identified skeletons """
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
    
    benzeneidx, brid = [], []      # indices and types of benzene ring in all rings
    nrg, setrg = len(ringIdx), [set(rg) for rg in ringIdx]
    rsmiles = [r.clone().smiles() for r in ringObjs]

    for i in xrange(nrg):

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
            benzeneidx.append(i)
            brid.append('b')
            continue

        if rsmile.count('C') != ringObjs[i].countAtoms(): continue
        
        nc, db, tp = rsmile.count('C'), rsmile.count("="), None
        
        if db==3 and nc==6:
            tp = 'b'

        if not tp and nc == 6 and db == 2:
            for j in xrange(nrg):
                if len(setrg[i]&setrg[j]) == 2 and rsmiles[j].count("=") >= 2:
                    tp = 'b'
                    break
                    
        if not tp and nc >= 5:
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

        if tp:
            benzeneidx.append(i)
            brid.append(tp)

    return benzeneidx, brid


def getRings(mol, checkBenzene = True):
    """ Get ring objects and indices for benzenes or benzene like rings """
    
    ringidx =[]           # rings containing atom objects
    
    if checkBenzene:

        ringObjs = []
        for r in mol.iterateRings(5,6):
            s = r.clone().smiles().lower()
            s = ''.join(findall("[a-zA-Z]+", s))
            if all(k in 'co' for k in s):
                ringidx.append([atom.index() for atom in r.iterateAtoms()])
                ringObjs.append(r)

        if len(ringidx)<2:  return None, None, None

        benzeneidx, brid = benzeneRingID(ringObjs, ringidx)

        return ringidx, benzeneidx, brid
    
    else:
        
        for r in mol.iterateRings(3,8):
            ringidx.append([atom.index() for atom in r.iterateAtoms()])
        return ringidx


def getRingA(ringidx, benzeneidx):
    """ Get candidate ring A, i.e., a benzene like ring with side chains -C- and -O-
        bonding same benzene ring bond.
        Return information around ring A: index of benzene ring that ring A belongs
        to, object of atom in benzene ring that bonds to C, object of C, object of
        atom in benzene ring that bonds to O, object of O.
    """
    # pre-declaration to faster the use of global variable
    fneis = neis
    
    ringA = []
    for i in benzeneidx:
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


def classFlav(ringidx, sk, brix, brid):
    """
    Classify flavonoids into subclasses
    """
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
            if noC[i] == 6 and oc:
                if ot and allb1:
                    if lj == 6:
                        return 'peltogynoid'
                    elif lj == 5:
                        return 'dihydroflavonol C3-C2\' ether linkage'
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

    typeA, typeB = brid[brix.index(tsk[0])], brid[brix.index(tsk[2])]
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
            
    if noC[i]==6 and allbr:
        
        if b4==2 and bc5:
            return 'flavonol' if b5==1 else 'flavone'
        elif b4==2 and allb1:
            return 'dihydroflavonol' if b5==1 else 'flavanone'
        elif b4==1 and allb1:
            return 'flavan-3,4-diol' if b5==1 else 'flavan-4-ol'
        elif bc4 and bc6:
            return 'anthocyanidin'
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


def checkStilbenoSK(ringidx, benzeneidx, sk, name):
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

    a = -1
    for skj in sk:
        a += 1
        if not name[a] != 'flavan-3-ol': continue
        
        for i,j in iters(benzeneidx, benzeneidx, flag=True):

            if i in skj[0] or j in skj[0]: continue

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
    
    return sk, name


def checkBiflavonoidSK(ringidx, benzeneidx, ringA):
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
            for k in benzeneidx:
                if any(ni[0] in ringidx[k] for ni in fneis[i]):
                    bz3.append(k)

            if len(bz3) == 2:
                p3.append((i, j, bz3))
    if not p3: return sk, None

    # get chromane ring
    cr = []
    for i in benzeneidx:
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
                for m in benzeneidx:
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


def checkFlavSK(ringidx, bridx, ringA, brid):
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

    nr = len(ringidx)
    setrg, lr, clr, rno, i = [], [], [], [], -1
    for ix in ringidx:
        i += 1
        setrg.append(set(ix))
        lr.append(len(ix))
        rno.append(len([j for j in ix if fatoms[j][1]=='o']))
        if lr[-1]==6 and rno[-1]==1:
            clr.append(i)

    if not clr:
        return [], []

    # set up ring Bs to include the identification of homoisoflavonoids
    # i.e., a (phenyl)methyl structure, which will be used as a benzyl as for
    # identification of flavonoids or isoflavonoids
    rgb = bridx[:]
    for i in bridx:
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
            if t[-1] not in setrg[k] or t[1] not in setrg[k]: continue
            if len(setrg[i]&setrg[k]) != 2: continue
            if len(setrg[jx]&setrg[k]) > 0: continue

            cs = [ix for ix in setrg[k] if ix not in setrg[i] and ix != rA[-1]]
            # remove C in C=O in ring C
            o = [ni[0] for ni in fneis[t[1]] if ni[1]==2 and ni[2]=='o']
            if o:
                o = o[0]
                cs.remove(t[1])

            # flavonoids or isoflavonoids
            if isinstance(j, int) and j != i:
                for ic in cs:
                    ix = [v[0] for v in fneis[ic] if v[0] in setrg[j] and v[2]=='c']
                    if ix:
                        ix = ix[0]
                        skt = [[i, k, j], t[0],t[1],ic,ix,t[2],t[3]]
                        if brid[rgb.index(j)] in 'bq':
                            skt += [o, 'f'] if o else ['f']
                        else:
                            skt += [o, 'fd'] if o else ['fd']
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
                        if o:
                            sk.append([[i,k,jx], t[0],t[1],ic,ix,j0,t[2],t[3],o,'hf'])
                        else:
                            sk.append([[i,k,jx], t[0],t[1],ic,ix,j0,t[2],t[3],'hf'])

    # chech whether multiple skeletons share same benzene rings; these skeletons
    # are possibly coumestans, pterocarpans or peltogynoid;
    # and check whether there exists a ring among ring C and ring B;
    # remove structures that can not be assigned by known class name;
    # if the chain in flavonoids (i.e., bond between ring B and C) is one of
    # the edges of another ring but the flavonoid is also not known class
    # having four rings, remove it.
    taidx, n, names = [ra[0] for ra in ringA], len(sk), [None]*len(sk)
    if n>0:
        sk4 = []
        for i in xrange(n):
            if sk[i][-1] != 'f': continue
            skr = sk[i][0]
            
            for j in xrange(nr):
                if len(setrg[j]&setrg[skr[1]])==2 and len(setrg[j]&setrg[skr[2]])==2:
                    if rno[j]==1 and not setrg[j]&setrg[skr[0]]:
                        skt = sk[i][:]
                        skt[0].append(j)
                        name = classFlav(ringidx, skt, bridx, brid)
                        if name:
                            sk[i][:] = skt[:]
                            sk4.append(i)
                            names[i] = name
                    elif rno[j]==0 and lr[j]==5:
                        for k in setrg[j]&setrg[skr[1]]:
                            o = [ni[0] for ni in fneis[k] if ni[2] == 'o']
                            t = any(ni[0] in ringidx[i] for ni in fneis[k])
                            if len(o) == 1 and not t:
                                sk[i][0].append(j)
                                sk[i] = sk[i][:-1]+[o[0], 'arl']     # 4-arylchroman
                                names[i] = classFlav(ringidx, sk[i], bridx, brid)
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
                    names[i] = classFlav(ringidx, sk[i], bridx, brid)
                    if n > 1:
                        # if this is flavonolignan, part of its rings is also
                        # a skeleton
                        skt = [sk[i][0][:-1]]+sk[i][1:]
                        skt[-1] = 'f'
                        name = classFlav(ringidx, skt, bridx, brid)
                        if name:
                            names.append(name)
                            sk.append(skt)
                
                # check whether a third ring exists between ring A and C or
                # ring B and ring between ring B and C, if so, remove this
                # skeleton
                g = False
                r1, r2, r3, r4 = tuple(sk[i][0])
                for k in xrange(nr):
                    if k not in sk[i][0]:
                        if len(setrg[r1]&setrg[k])>0 and len(setrg[r4]&setrg[k])>0:
                            g = True
                            break
                        if len(setrg[r2]&setrg[k])>0 and len(setrg[r3]&setrg[k])>0:
                            g = True
                            break
                        
                if g:
                    names[i] = None
            else:
                # bond topology between ring B and C
                tp = [ni[4] for ni in fneis[sk[i][3]] if ni[0]==sk[i][4]][0]
                if tp!=CHAIN: continue
                
                if sk[i][-1] != 'f':
                    # if this is not standard flavonoid (e.g., flavone, with
                    # name notation of 'f') and share same ring with standard
                    # flavonoid, ignore it.
                    if all(not setsk[i]&setsk[j] for j in xrange(n) if sk[j][-1]=='f'):
                        names[i] = classFlav(ringidx, sk[i], bridx, brid)
                else:
                    names[i] = classFlav(ringidx, sk[i], bridx, brid)

        n = len(sk)
        sk[:], names[:]= [sk[i] for i in xrange(n) if names[i]], [m for m in names if m]

    return sk, names


def checkChalSK(ringidx, bridx, brid):
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
    
    sk, name =[], []
    setrg,nr,lr =[set(rg) for rg in ringidx],len(ringidx),[len(rg) for rg in ringidx]

    # get neighbor Cs of a benzene ring
    brneis, delidx = [], []
    for i in bridx:
        tnei=[]
        for j in ringidx[i]:
            tnei+=[(n,j) for n in fneis[j] if n[2]=='c' and n[0] not in ringidx[i]]
        brneis.append(tnei)

    assignedbr = set([])
    for i, j in iters(bridx, bridx, flag=True):
        
        ia, ib = bridx.index(i), bridx.index(j)
        
        if not brneis[ia] or not brneis[ib] or i in assignedbr or j in assignedbr:
            continue

        t = False
        for iinfo, jinfo in iters(brneis[ia], brneis[ib]):
            ni, i0 = iinfo      # neighbors and corresponding atom in the ring
            nj, j0 = jinfo
            i1, j1 = ni[0], nj[0]   # indices of the neighbor atoms
            
            if i1 == j1: continue

            for ni2, nj2 in iters(fneis[i1], fneis[j1]):
                if ni2[0] == nj2[0] and ni2[2] == 'c':
                    t = True
                    break
            if t: break

        if t:

            noxy1 = [v[0] for v in fneis[i1] if v[2]=='o' and v[1]==2]
            noxy2 = [v[0] for v in fneis[j1] if v[2]=='o' and v[1]==2]

            t1 = False
            if noxy1 or noxy2:
                # aurone
                ix, jx, nxi, nxj = (i, j, ni2[0], nj2) if noxy1 else (j, i, nj2[0], ni2)
                ix0, jx0, oi, ix1, jx1 = (i0, j0, noxy1[0], i1, j1) if noxy1 else\
                                          (j0, i0, noxy2[0], j1, i1)
                
                for k in xrange(nr):
                    o = []
                    if lr[k]==5 and nxi in setrg[k] and len(setrg[k]&setrg[ix])==2:
                        for ik in setrg[k]&setrg[ix]:
                            o+=[v for v in fneis[ik] if v[2]=='o' and v[0] in setrg[k]]
                            break
                            
                        if o:
                            o = o[0][0]
                            t1 = True
                            name.append('aurone' if nxj[1] == 2 else 'dihydroaurone')
                            sk.append([[ix, k, jx],ix0,ix1,nxi,jx1,jx0,ik,o,oi,'au'])
                            break

            if t1:
                assignedbr.update(sk[-1][0])
                continue
            
            # for these types of flavonoids except aurone and dihydroaurone,
            # no any atom in the chain between two benzenes exists in other
            # rings.
            tps  = [i1, ni2[0], j1]
            t = False
            for k in tps:
                if any(ni[4]!=CHAIN for ni in fneis[k]):
                    t = True
                    break
            if t: continue

            tname = []
            if noxy1 and not noxy2:
                if brid[ia] == 'b':
                    tname.append('chalcone' if nj2[1] >= 2 else 'dihydrochalcone')
                elif brid[ia] in 'mq2':
                    tname.append('chalcone-quinol')
                if tname:
                    sk.append([[i,j], i0, i1, ni2[0], j1, j0, noxy1[0], 'c'])
                    name.append(tname)
        
            elif noxy2 and not noxy1:
                if brid[ib] == 'b':
                    tname.append('chalcone' if ni2[1] >= 2 else 'dihydrochalcone')
                elif brid[ib] in 'mq2':
                    tname.append('chalcone-quinol')
                if tname:
                    sk.append([[j,i], j0, j1, nj2[0], i1, i0, noxy2[0], 'c'])
                    name.append(tname)
                    
            elif noxy1 and noxy2:
                sk.append([[i,j], i0, i1, ni2[0], j1, j0, noxy1[0], noxy2[0], 'o'])
                name.append('oxodihydrochalcone')
                
            elif ni2[1]==1 and nj2[1]==1:
                noxy1 = [v[0] for v in fneis[i1] if v[2]=='o' and v[1]==1]
                noxy2 = [v[0] for v in fneis[j1] if v[2]=='o' and v[1]==1]
                if noxy1:
                    sk.append([[i,j], i0, i1, ni2[0], j1, j0, noxy1[0], 'c'])
                    name.append('pentahydrochalcone')
                elif noxy2:
                    sk.append([[j,i], j0, j1, nj2[0], i1, i0, noxy2[0], 'c'])
                    name.append('pentahydrochalcone')
            
            elif lr[i] == 6 and brid[ia]=='q2':
                if ni[1]==2 and brid[ib]=='b' and nj2[1]==2:
                    noxy = [v[0] for v in fneis[i1] if v[2]=='o' and v[1]==1]
                    if noxy:
                        sk.append([[i,j], i0, i1, ni2[0], j1, j0, noxy[0], 'c'])
                        name.append('chalcone-quinol')
            
            elif lr[i] == 5:
                noxy = [v[0] for v in fneis[i1] if v[2]=='o' and v[1]==1]
                if noxy:
                    o2 = []
                    for k in ringidx[i]:
                        xt = [ni[0] for ni in fneis[k] if ni[2]=='o' and ni[1]==2]
                        if xt:  o2.append(xt[0])
                    sk.append([[i,j], i0, i1, ni2[0], j1, j0, noxy[0]]+o2+['c'])
                    name.append('chalcone_cyclopentenedione')

            if sk: assignedbr.update(sk[-1][0])

    return sk, name


def checkAnthoSK(ringidx, benzeneidx, oix):
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
        for i in benzeneidx:
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


def checkSpecialSK(ringidx, benzeneidx, brid, ringA):
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
        brs = [benzeneidx[i] for i in xrange(len(brid)) if brid[i]=='b']
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

    for i,j in iters(benzeneidx, benzeneidx, flag=True):
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
                              and ni[0]!=i0 and ni[4]==CHAIN and ni[0]!=i1]
                        if c1:
                            sk.append([[i, j], i0, i1, j1, j0, o[0], c1[0], 'm'])
                    elif brid[benzeneidx.index(i)]=='q':
                        c = [ni[0] for ni in fneis[i1]
                             if ni[2]=='c' and ni[0]!=i0 and ni[0]!=j0 and ni[4]==CHAIN]
                        if c:
                            c = c[0]
                            c1 = [ni[0] for ni in fneis[c] if ni[0]!=c
                                   and ni[1]==2 and ni[4]==CHAIN and ni[2]=='c']
                            if c1:
                                sk.append([[i, j], i0, i1, j0, c, c1[0], 'd'])

    return sk


def checkMultipleSKs(ringidx, sk):
    """
    Check skeletons to identify whether those skeletons share same rings or
    bonds, if so, maximum independent set algorithm is used to get skeleton set
    with most number of independent skeletons (i.e., no sharing bonds between
    any two skeletons). Also, a priori check is performed if two skeletons with
    different number of rings share same bonds. In this case, skeleton with
    more rings is retained.
    Output: a list with indices retained
    """
    n = len(sk)
    skrs, lenskr = [ski[0] for ski in sk], [len(ski[0]) for ski in sk]
    sets = [set(ring) for ring in ringidx]

    c = lambda sk1, sk2, s: len(set(sk1[0])&set(sk2[0]))>0 or\
        any(len(set(sk1[1:])&s[k])>0 for k in sk2[0]) or\
        any(len(set(sk2[1:])&s[k])>0 for k in sk1[0])

    delix = []
    for i in xrange(n):
        if all(c(sk[i], sk[j], sets) and j!=i for j in xrange(n)):
            delix.append(i)
    
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
    
    ringidx, benzeneidx, brid = getRings(mol)
    if not ringidx or not benzeneidx:
        return None, None, None, None

    n, nr = len(benzeneidx), len(ringidx)

    # tetrahydroflavanones
    pseudbr, nx = [], -1
    for r in ringidx:
        nx += 1
        if nx not in benzeneidx and all(fatom[i][1]=='c' for i in r) and len(r)==6:
            b = []
            for i in r:
                b += [ni[1]==2 for ni in fneis[i] if ni[0] in r]
            if sum(b)==2:
                pseudbr.append(nx)

    if pseudbr:
        bridx = ['b']*(len(pseudbr)+n)+brid
        ringAx = getRingA(ringidx, pseudbr)
        fsk, name = checkFlavSK(ringidx, pseudbr+benzeneidx, ringAx, bridx)
        if len(fsk)==1 and fsk[0][0][2] in benzeneidx and fsk[0][0][0] in pseudbr:
            if name[0] == 'flavanone':
                skix = getskidx(ringidx, fsk)
                return fsk,skix,'tetrahydroflavanones',(ringidx,benzeneidx,brid)
    
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
            raise ValueError('Charged atoms except oxygen or charges except +1 are not allowed as a candidate of flavonoid.')

    if oix:
        ask, name = checkAnthoSK(ringidx, benzeneidx, oix)
        if ask:
            sk[:], names[:] = ask, name
            n -= na(sk, benzeneidx)
            if n < 2:
                skix = getskidx(ringidx, sk)
                return sk, skix, names[0], (ringidx, benzeneidx, brid)
            # remove assigned benzene rings
            rmix, bzix, brid2 = [], [], []
            for ski in sk:  rmix += ski[0]
            for i in benzeneidx:
                if i not in rmix:
                    bzix.append(i)
                    brid2.append(brid[benzeneidx.index(i)])
        else:
            if t: return None, None, None, None

    if 'bzix' not in locals():
        bzix, brid2 = benzeneidx, brid

    ringA = getRingA(ringidx,bzix)
    benzenes = [bzix[i] for i in xrange(n) if brid2[i]=='b']

    # flavonoids
    flavsk, name = checkFlavSK(ringidx, bzix, ringA, brid2)
    if flavsk:
        flavsk, name = checkStilbenoSK(ringidx, benzenes, flavsk, name)
        nf = len(flavsk)
        sk[:] = [flavsk[i] for i in xrange(nf) if name[i]]
        names[:] = [name[i] for i in xrange(nf) if name[i]]

    # chalcone
    csk, name = checkChalSK(ringidx, bzix, brid2)
    if csk:
        sk += csk
        names += name

    n -= na(sk, benzeneidx)
    # special flavonoids
    if n > 1:
        xsk = checkXanthSK(ringidx, benzenes)
        osk = checkSpecialSK(ringidx, bzix, brid2, ringA)
        osk += xsk
        for ski in osk:
            name = classFlav(ringidx, ski, bzix, brid2)
            names.append(name)
        sk += osk

    # specify name to flavonoids
    nf = len(sk)
    if nf > 1:
        sfsk, sfname = [], []
        for ski in sk:
            if ski[-1] == 'sf':
                sfsk[:] = ski[:]
                sfname[:] = 'stilbeno-flavonoid'
                break
        if sfsk:
            rix = []
            for i in xrange(nf):
                if not set(sk[i][0])&set(sfsk[0]):
                    rix.append(i)
        else:
            rix = list(xrange(nf))

        if len(rix)>1:
            retains = checkMultipleSKs(ringidx, [sk[i] for i in rix])
            rix[:] = [rix[i] for i in retains]
            sk[:] = sfsk+[sk[i] for i in rix]
            names[:] = sfname+[names[i] for i in rix]

    if nf==0:
        sk, name = checkBiflavonoidSK(ringidx, benzenes, ringA)
        skix = getskidx(ringidx, sk)
        return sk, skix, name, (ringidx, benzeneidx, brid)

    skix = getskidx(ringidx, sk)

    return sk, skix, names, (ringidx, benzeneidx, brid)


"""
Indexing the skeleton to get the site of side chains. This procedure follows
the IUPAC provisional recommendation for Nomenclature of Flavonoids
(Rauter A. P., et al. 2009). Note that this function only considers the side
groups that bond to skeleton and saccharides.
"""

def indexring(ring, stinfo, sgx, flag):
    """
    Index ring according to start atom information (stinfo) and
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

    skgix, sgix = [], []
    
    setrg, lr, lsk = [set(rg) for rg in ringidx], [len(rg) for rg in ringidx], len(sk)

    sgs = []
    for key in sg.keys():
        for g in sg[key]['atomIndex']:
            sgs += g
    if not sgs and len(sk)==1:
        return skgix, sgix
    
    exs = ['coumestan', 'pterocarpan', 'dihydropterocarpan']
    skt, ixt, sgsi = [], [], []
    for ii in xrange(lsk):
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

            # ring C
            rgt = indexring(ringidx[k][:], [i0,i1,o], [], 'c')
            skt += rgt
            ixt += [str(i+1) for i in xrange(lr[k]-2)]+[None, None]

            # ring B
            rgt = indexring(ringidx[j][:], [None,None,j0], sgsi, 'b')
            skt += rgt
            ixt += [str(i+1)+"'" for i in xrange(lr[j])]
            
        elif l == 2:

            i, j = ski[0][0], ski[0][1]

            # ring A
            rgt = indexring(ringidx[i][:], [None,None,ski[1]], sgsi, 'b')
            skt += rgt
            ixt += [str(i+1)+"'" for i in xrange(lr[i])]

            # ring C
            j0 = [k for k in ski[1:] if k in ringidx[j]][0]
            rgt = indexring(ringidx[j][:], [None,None,j0], sgsi, 'b')
            skt += rgt
            ixt += [str(i+1) for i in xrange(lr[j])]

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

        for i in xrange(len(skt)):
            if ixt[i]:
                for key in sg.keys():
                    j = -1
                    for g in sg[key]['atomIndex']:
                        j += 1
                        if skt[i] in g or any(ni[0] in g for ni in fneis[skt[i]]):
                            skgixi.append(ixt[i])
                            sgixi.append((key, j, skt[i]))

                for ij in xrange(lsk):
                    if ij != ii:
                        if any(ni[0] in skix[ij] for ni in fneis[skt[i]]):
                            skgixi.append(ixt[i])
                            sgixi.append(('sk', ij, skt[i]))
        
        skgix.append(skgixi)
        sgix.append(sgixi)

        skt[:], ixt[:], sgsi[:] = [], [], []

    return skgix, sgix
 

"""
Get the side chains of current molecule. List of side chains can be found in
side_chain module. So if new side chains are found, it can be added to the end
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

def getsidenets(sk, skix, sginfo):
    """ get networks of side groups """
    sgps, skgix, sgix = sginfo
    fneis = neis
    if not sgix: return 'There is not any side group around the skeleton.\n'

    lsk, skixset = len(sk), [set(ix) for ix in skix]
    
    sdneis, sdsets, sdinfo = {}, {}, {}
    for key in sgps.keys():
        if not sgps[key]['atomIndex']: continue
        kx = 'ring' if 'ring' in key else 'chain'

        cns = []
        for gi in sgps[key]['atomIndex']:
            t = set()
            for j in gi:
                t.update([ni[0] for ni in fneis[j]])
            cns.append(t)
        sdneis[key] = cns
        sdsets[key] = [set(gi) for gi in sgps[key]['atomIndex']]
        sdinfo[key] = kx
    
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
                kx, t = sdinfo[key2], ''
                for j in xrange(len(sgps[key2]['groupName'])):
                    if j != i and key2!=key and sdsets[key][i]&sdneis[key2][j]:
                        t += ' %d,'%j
                if t:
                    s += '%s groups%s' %(kx,t) if len(t)>3 else '%s group%s' %(kx,t)
            
            s = s[:-1] + '\n'
    return s


def getsidegroupdist(sk, names, sginfo):
    """ get distributions of side groups around skeleton """
    sgps, skgix, sgix = sginfo
    if not sgix: return 'There is not any side group around the skeleton.\n'

    fatoms, fneis = atoms, neis
    glcnames = ['glucuronide','alloside','galactoside','glucoside',
                'rhamnoside','apioside','arabinofuranoside','xyloside',
                'xylopyranoside','arabinoside']
    
    lsg = [len(g) for g in sgix]
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
                if gk=='ringGroup' and any(sn in gpname for sn in glcnames):
                    for ix in sgps[gk]['atomIndex'][ik]:
                        if any(ni[0] == gmi for ni in fneis[ix]):
                            namepre = 'O-' if fatoms[ix][1]=='o' else 'C-'
                            if 'or' in gpname:
                                gpname = '%s(%s)'%(namepre,gpname)
                            else:
                                gpname = '%s%s'%(namepre,gpname)
                            break
                si = "%s-%s (%d)"%(skgix[i][j], gpname, ik)
                sc.append(si)
        if sc:
            s += 'skeleton %d (%s): '%(i+1, names[i]) + '; '.join(sc) + '\n'
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


def isvalidringmatch(mhi, rgset, qix, mol, qinfo):
    """
    check whether the match of ring group is valid
    """
    # check whether any of ring in rgset is subset of matched group
    # if not, return False
    if not any(rg.issubset(mhi) for rg in rgset):
        return False

    # check whether the chain is part of a ring, if so, return False thus
    # the match is rejected.
    fneis = neis
    i, chainix, q, qsml, isbr, num = qinfo
    chainset = set([mhi[i] for i in xrange(num) if qix[i] in chainix])
    
    for rg in rgset:
        cmr = chainset&rg
        if len(cmr)>=2:
            for j in cmr:
                if any(nei[0] in cmr for nei in fneis[j]):
                    return False
    
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


def checkQueryGroup(ringidx, name, ckr=True):
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
                        namec.append('Phenylethylene')
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
                                    namec.append('Caffeic acid')
                                else:
                                    oh = [ni[0] for ni in fneis[oc[2]] if ni[2]=='o'
                                          and ni[4]==CHAIN and ni[1]==1]
                                    if oh:
                                        gpc.append([c1,c2,k,o1,o2,o3[0],oh[0]])
                                        namec.append('Caffeic acid')
                                    else:
                                        gpc.append([c1,c2,k,o1,o2,o3[0]])
                                        namec.append('p-Coumaric acid')
                            else:
                                gpc.append([c1, c2, k, o1, o2])
                                namec.append('Cinnamic acid')
                            
                        elif o1:
                            gpc.append([c1, c2, k, o1[0]])
                            namec.append('Cinnamyl alcohol')
                        else:
                            gpc.append([c1, c2, k])
                            namec.append('Hydroxystyryl')
        if gpc:
            l = [len(g) for g in gpc]
            i = l.index(max(l))
            gp.append(rg+gpc[i])
            name.append(namec[i])

    return gp, name


def ringGroupMatch(mol, rgsets, querys, qidxs):
    """
    Find ring groups from current molecule
    """
    # pre-declare local functions and variables
    loadQueryMolecule   = idg.loadQueryMolecule
    substructureMatcher = idg.substructureMatcher
    rmatoms             = removeatoms

    # preallocation and initialization
    matchedix, mgix  = [], []
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

    # match groups
    rmrg, delatms, delrg, ismatched = set(), set(), [], False
    while True:
        nm = 0
        for qc in querys:
            i, chainix, q, qsml, bi, num = qc
            # if length of ring group smaller than, or number of rings
            # larger than that of the query, or no benzene ring exists in
            # side chain while the query group has, continue
            if not havebr and bi and num>msd: continue
            
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

            if ismatched:
                if havebr:
                    t = True
                    for g in rgset:
                        if any(subi(g) for subi in brsubgroupid):
                            t = False
                            break
                    if t:
                        havebr = False

                # get the max length of side groups containing rings
                rix[:] = []
                for g in sdset:
                    apdrg = []
                    for rg in rgset:
                        if g.intersection(rg):
                            apdrg += rg
                    if apdrg:
                        rix.append(g.union(apdrg))
                sdset[:] = rix
                msd, nl = max(len(g) for g in sdset), len(sdset)
                ismatched = False
        
        if nm==0:
            break

    # further checking by delete an atom around the ring, this is restricted to
    # glucose and similar groups
    iis = [qi[0] for qi in querys]
    for qx in qidxs:
        i = qx[1]
        ix = iis.index(i)
        for j in qx[0]:
            q = loadQueryMolecule(querys[ix][3])
            q.getAtom(j).remove()
            qatoms  =  [atom for atom in q.iterateAtoms()]
            delatms.clear()
            for mhi in substructureMatcher(mol).iterateMatches(q):
                c = [mhi.mapAtom(atm).index() for atm in qatoms]
                matchedix.append(c)
                mgix.append(i)
                delatms.update(c)
            if delatms:
                rmatoms(mol, delatms, rgset, (rmrg, brset), False)

##    idgrender.renderToFile(mol,'m2.png')
    return matchedix, mgix


def chainGroupMatch(mol, querys):
    """ Find side chains from molecule """
    # pre-declare local functions and variables
    substructureMatcher = idg.substructureMatcher
    iterateMatches      = substructureMatcher(mol).iterateMatches
    
    # preallocation and initialization
    n, matchedix, mgix, delatms  = len(querys), [], [], set([])

    # get separated side groups
    sdgroups = separategroups(mol)
    sdset = [set(g) for g in sdgroups]
    maxm, nl = max(len(g) for g in sdset), len(sdset)

    # match groups
    ismatched = False
    for i in xrange(n):
        q, num, qsml = querys[i]
        if num>maxm: continue
        
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
    skidx, rmrg, brs, rmrgix = set(), set(), [], set()
    for sk in skeleton:
        rmrgix.update(sk[0])
        for i in sk[0]:
            rmrg.update(ringinfo[0][i])
    for skixi in skix: skidx.update(skixi)
    rgset = [set(ringinfo[0][i]) for i in xrange(len(ringinfo[0])) if i not in rmrgix]
    rgset += [set(rg) for rg in rgix if rg not in ringinfo[0]]
    
    # obtain benzene rings
    for i in ringinfo[1]:
        if ringinfo[2][ringinfo[1].index(i)] == 'b':
            brs.append(ringinfo[0][i])

    if len(set(skidx)) == mol.countAtoms():
        return {'ringGroup':{'atomIndex':[], 'groupName':[]},
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
            querysds, qidxs, qnames, nix = checkQueryGroup(rgset, name)
        else:
            querysds, qidxs, qnames = checkQueryGroup(rgset, [])

        rgset2 = rgset[:]

        rmatoms(mol, skidx, rgset, (rmrg, brset), True)
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

        mhrnames = [qnames[i] for i in mhgs]
        mhs = set()
        for mhi in mhidx: mhs.update(mhi)
        setig = set([atom.index() for atom in mol.iterateAtoms()])&mhs
        for rg in rgset: rmrg.update(rg)
    else:
        setig = set(skidx)
        mhidx, mhrnames = [], []

    if mol.countAtoms()>len(setig):
        rmatoms(mol, setig, [], (rmrg, []), False)
        querysds, qnames = checkQueryGroup([], [], ckr=False)
        mhidxc, mhgsc = chainGroupMatch(mol, querysds)
        # delete groups that have been assigned in other groups
        if mhidx:
            delix = []
            mhrset = [set(g) for g in mhidx]
            for gi in mhidxc:
                if any(len(set(gi)&g)==len(gi) for g in mhrset):
                    delix.append(gi)
            mhgsc[:] = [mhgsc[i] for i in xrange(len(mhgsc)) if mhidxc[i] not in delix]
            mhidxc[:] = [g for g in mhidxc if g not in delix]
        mhcnames = [qnames[i] for i in mhgsc]
            
    else:
        mhidxc, mhcnames = [], []

    mg = {'ringGroup':  {'atomIndex': mhidx,    'groupName':mhrnames},
          'chainGroup': {'atomIndex': mhidxc,   'groupName':mhcnames}}

    if mol.countAtoms()>0:
        return mg, False
    
    return mg, True    


def getmoleculeInfo(string):
    """ Get information of molecule """
    
    mol = loadmol(string)
    
    mol.dearomatize()               # convert aromatic structure to Kekule form
    sk, skix, names, ringinfo =  getSkeleton(mol)
    
    if not sk or not names:
        return None, None, None, None

    mg, validsd = getSidechains(mol, sk, skix, ringinfo)
    if not validsd:
        raise ValueError('Unexcepted side groups or elements exist in flavonoids.')

    skgix, sgix = indexsg(sk, skix, mg, ringinfo[0], names)

    return sk, skix, names, (mg, skgix, sgix)


if __name__=='__main__':
    s = 'C1C(OC2=C(C1=O)C=CC(=C2)OC3C(C(C(C(O3)CO)O)O)O)C4=CC=CC=C4'
    sk, name, mg, sg = getmoleculeInfo(s)
    print "skeleton:", sk, '\nname:', name, '\nside chains:', mg
