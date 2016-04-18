"""
Read and write side chain properties to a json file (sidechain.json) used to
store the information of side chains. The write function can also add
properties and then store them into the json file.
"""

import json
from indigo import Indigo

idg = Indigo()


def addproperty(name, value, ringclass='ring'):
    """
    add property to the side chains
    """
    with open('sidechain.json','r') as f:
        sd = json.load(f)

    key = 'ringGroup' if ringclass == 'ring' else 'chainGroup'

    for i in xrange(len(sd[key])):
        sd[key][i]['property'][name] = value[i]
    

    with open('sidechain.json','w') as f:
        json.dump(sd, f, indent=4, separators=(',',': '), sort_keys=True)



def addgroup(name, smiles):
    """
    add groups as candiates of side groups of flavonoids
    "name": "Oxirane",
            "property": {
                "checkaromatic": false,
                "checkglc": false,
                "haveBenzene": false,
                "number": 3,
                "ringnumber": 1
            },
            "smiles": "C1CO1"
    """
    loadmol = idg.loadMolecule
    exactMatch = idg.exactMatch
    
    m = loadmol(smiles)
    
    with open('sidechain.json','r') as f:
        sd = json.load(f)

    ks = sd.keys()
    sid = 'ring' if m.countSSSR()>0 else 'chain'
    key = [s for s in ks if sid in s][0]

    n = m.countAtoms()

    # if the group has been stored, ignore it
    for g in sd[key]:
        if g['property']['number'] == n:
            q = loadmol(g['smiles'])
            if exactMatch(m, q, 'TAU'):
                return

    rgset, havebr = [], False
    for r in m.iterateRings(3,6):
        rs = r.clone().smiles()
        if rs.count('c')==6 or (rs.count('=')==3 and rs.count('C')==6):
            havebr = True
        rgset.append(set(a.index() for a in r.iterateAtoms()))

    # delete rings that is super set of other two or more rings
    nr = len(rgset)
    if nr>1:
        delix = []
        for i in xrange(nr):
            if i in delix: continue
            k = 0
            for j in xrange(nr):
                if j != i and rgset[i].issuperset(rgset[j]):
                    k += 1
            if k>0: delix.append(i)
        if delix:
            rgset[:] = [rgset[i] for i in xrange(nr) if i not in delix]
            nr -= len(delix)
        
    # assign property to the molecule
    p = {}
    p['property'] = {}
    p['name'] = name
    p['smiles'] = smiles
    p['property']['number'] = n
    p['property']['ringnumber'] = nr
    p['property']['haveBenzene'] = havebr
    p['property']['checkaromatic'] = False
    p['property']['checkglc'] = False

    sd[key].append(p)

    # sort side groups
    nsd = []
    ns = [g['property']['number'] for g in sd[key]]
    uns = list(set(ns))
    uns.sort(reverse=True)
    for i in uns:
        sdcurr = [sd[key][j] for j, ni in enumerate(ns) if ni==i]
        sdcurr[:] = sorted(sdcurr, key=lambda k: k['name'])
        nsd += sdcurr
    sd[key] = nsd

    with open('sidechain.json','w') as f:
        json.dump(sd, f, indent=4, separators=(',',': '), sort_keys=True)
        
