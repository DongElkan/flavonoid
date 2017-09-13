"""
    This module collects all the information of the input flavonoid via the molecular
identifier like canonical SMILES, InchI or a file. Currently ..mol (file stored in
ChemSpider database) and ..sdf (in PubChem database) are accepted.
"""

import os, sys
from multiprocessing import Process
import mol_interpreter as mi
from mol_interpreter import FlavonoidException, loadmol, validmolcheck, getnames, getinchi, getsidenets, getsidegroupdist, generatemolimage, sideChains
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
        raise FlavonoidException(1)

    if not validmolcheck(mol.canonicalSmiles()):
        raise FlavonoidException(2)

    if mol.countSSSR() >= 15:
            p = Process(target=mi.getmoleculeInfo, args=(string, ))
            p.start()

            # Wait 6 seconds for matching groups
            p.join(6)

            if p.is_alive():
                # Terminate
                p.terminate()
                p.join()
                raise FlavonoidException(3)
            
    return mol


class Flavonoids(object):
    """ flavonoid class to get the information of the molecule """
    def __init__(self, string):

        mol = validflavonoid(string)

        sk, skix, names, sginfo, skorderinfo = mi.getmoleculeInfo(string)
        if not sk:
            raise FlavonoidException(4)

        self.__molecule = mol
        self.__identifier = string
        self.__skinfo = (sk, skix)
        self.__names = names
        self.__sginfo = sginfo[:-1]
        self.__glcnumbering = sginfo[-1]

    def getcanonicalSmiles(self):
        """ canonical SMILES """
        return self.__molecule.canonicalSmiles()

    def className(self):
        """ class name of the flavonoid """
        return getnames(self.__names)
    
    def exactMass(self):
        """ exact mass of input molecule """
        return self.__molecule.monoisotopicMass()

    def mass(self):
        """ average mass """
        return self.__molecule.molecularWeight()

    def formula(self):
        """ formula of the molecule """
        return ''.join(self.__molecule.grossFormula().split())

    def getInchi(self):
        """ get InChI """
        return getinchi(self.__identifier)

    def getallSidegroups(self):
        """ get all side groups """
        s = 'Due to complex structure of the skeleton, the indexing ' +\
            'of side groups is ignored, only all side group nanes ' +\
            'are printed. For the classes that ignore this, please ' +\
            'see attribute "idxNameExceptions".'
        names = [self.__names] if isinstance(self.__names, str) else self.__names
        if any(name in mi.name_exceptions for name in names):
            print getsidenets(self.__skinfo[0], self.__skinfo[1],
                              self.__sginfo, self.__glcnumbering)
            sys.stderr.write('Warning: %s'%s)
            return None
            
        print getsidenets(self.__skinfo[0], self.__skinfo[1],
                          self.__sginfo, self.__glcnumbering)

    def getSidegroups(self):
        """ get side groups around the skeleton """
        s = 'Due to complex structure of the skeleton, the indexing ' +\
            'of side groups is ignored, only all side group nanes ' +\
            'are printed. For the classes that ignore this, please ' +\
            'see attribute "idxNameExceptions".'
        names = [self.__names] if isinstance(self.__names, str) else self.__names
        if any(name in mi.name_exceptions for name in names):
            print getsidegroupdist(self.__skinfo[0], names, self.__sginfo)
            sys.stderr.write('Warning: %s'%s)
            return None
        
        print getsidegroupdist(self.__skinfo[0], self.__names, self.__sginfo)

    def generateImage(self, filename):
        """ Generate image of current molecule """
        generatemolimage(self.__identifier, filename)

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
                generatemolimage(g['smiles'],'%s.png'%g['name'])

        os.chdir(cwd)

    def idxNameExceptions(self):
        """ Name exceptions in performing indexing of side groups """
        s = 'Currently, indexing of side groups for the following ' +\
            'classes are not considered: \n'
        print s + '\t' + ';\n\t'.join(mi.name_exceptions)

    def addsidegroup(self, name, smiles):
        """ add side groups to current side group library """
        addgroup(name, smiles)
        print "The group %s has been added.\n"%name
        reload(mi)
        
