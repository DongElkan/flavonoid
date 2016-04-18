"""
    This module collects all the information of the input flavonoid via the molecular
identifier like canonical SMILES, InchI or a file. Currently ..mol (file stored in
ChemSpider database) and ..sdf (in PubChem database) are accepted.
"""

import os
from mol_interpreter import getmoleculeInfo, loadmol, getnames, getinchi, getsidegroupdist, generatemolimage, sideChains, idgrender


class Flavonoids(object):
    """ flavonoid class to get the information of the molecule """
    def __init__(self, string):
        sk, names, sginfo = getmoleculeInfo(string)
        if not sk:
            raise ValueError('This is not a flavonoid.')
        
        mol = loadmol(string)

        self.molecule = mol
        self.identifier = string
        self.skinfo = (sk, names)
        self.sginfo = sginfo

    def getcanonicalSmiles(self):
        """ canonical SMILES """
        return self.molecule.canonicalSmiles()

    def className(self):
        """ class name of the flavonoid """
        return getnames(self.skinfo[1])
    
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
        return getinchi(self.identifier)

    def getSidegroups(self):
        """ get side groups around the skeleton """
        print getsidegroupdist(self.skinfo[0], self.skinfo[1], self.sginfo)

    def generateImage(self, filename):
        """ Generate image of current molecule """
        generatemolimage(self.identifier, filename)

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
                idgrender.renderToFile(gm,'%s.png'%g['name'])

        os.chdir(cwd)
                
