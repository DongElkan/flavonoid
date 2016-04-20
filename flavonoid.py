"""
    This module collects all the information of the input flavonoid via the molecular
identifier like canonical SMILES, InchI or a file. Currently ..mol (file stored in
ChemSpider database) and ..sdf (in PubChem database) are accepted.
"""

import os, sys
from multiprocessing import Process
import mol_interpreter as mi
from sidechain_props import addgroup

if sys.version_info[0]==3:
    try:
        from importlib import reload
    except ImportError:
        from imp import reload


def validflavonoid(string):
    """ check whether the input string is valid for identifying a flavonoid """
    mol = mi.loadmol(string)

    if mol.countComponents()>1:
        raise ValueError('Multiple components are unsupported in identifying flavonoids.')

    if not mi.validmolcheck(mol.canonicalSmiles()):
        raise ValueError('Elements except CHNOS are not allowed as a candidate of flavonoid.')

    if mol.countSSSR() >= 20:
            p = Process(target=mi.getmoleculeInfo, args=(string, ))
            p.start()

            # Wait 8 seconds for matching groups
            p.join(6)

            if p.is_alive():
                # Terminate
                p.terminate()
                p.join()
                raise RuntimeError('Group matching times out.')
            
    return mol


class Flavonoids(object):
    """ flavonoid class to get the information of the molecule """
    def __init__(self, string):

        mol = validflavonoid(string)

        sk, skix, names, sginfo = mi.getmoleculeInfo(string)
        if not sk:
            raise ValueError('This is not a flavonoid.')

        self.molecule = mol
        self.identifier = string
        self.skinfo = (sk, skix, names)
        self.sginfo = sginfo

    def getcanonicalSmiles(self):
        """ canonical SMILES """
        return self.molecule.canonicalSmiles()

    def className(self):
        """ class name of the flavonoid """
        return mi.getnames(self.skinfo[2])
    
    def exactMass(self):
        """ exact mass of input molecule """
        return self.molecule.monoisotopicMass()

    def mass(self):
        """ average mass """
        return self.molecule.molecularWeight()

    def formula(self):
        """ formula of the molecule """
        return ''.join(self.molecule.grossFormula().split())

    def getInchi(self):
        """ get InChI """
        return mi.getinchi(self.identifier)

    def getallSidegroups(self):
        """ get all side groups """
        print mi.getsidenets(self.skinfo[0], self.skinfo[1], self.sginfo)

    def getSidegroups(self):
        """ get side groups around the skeleton """
        print mi.getsidegroupdist(self.skinfo[0], self.skinfo[2], self.sginfo)

    def generateImage(self, filename):
        """ Generate image of current molecule """
        mi.generatemolimage(self.identifier, filename)

    def generateSidegroups(self, path=None):
        """ generate all side groups """
        cwd = os.getcwd()
        if not path:
            if not os.path.isdir('sidegroups'):
                os.makedirs('sidegroups')
            os.chdir(os.path.join(cwd, 'sidegroups'))
        else:
            if not os.path.isabs(path):
                path = os.path.join(cwd, path)

            if not os.path.isdir(path):
                os.makedirs(path)
            os.chdir(path)
                
        for key in sideChains.keys():
            for g in sideChains[key]:
                gm = loadmol(g['smiles'])
                mi.idgrender.renderToFile(gm,'%s.png'%g['name'])

        os.chdir(cwd)

    def addsidegroup(self, name, smiles):
        """ add side groups to current side group library """
        addgroup(name, smiles)
        print "The group %s has been added.\n"%name
        reload(mi)
        
