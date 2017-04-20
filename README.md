"flavonoid" is a pure Python module used to identify the class of flavonoids according to the input molecular identifiers and facilitate the interpretation of mass spectra derived from the molecules, with the help of indigo toolkit.

**Usage**
    
    from flavonoid import Flavonoids
    fv = Flavonoids(s)

Where *s* is an string of molecular identifier. Currently, canonical SMILES, InChI, and chemical file format ..mol and ..sdf are supported. Once Flavonoids object "fv" is constructed, properties of the molecule can be accessed by the methods built.
For example, the class of the input molecule can be obtained via

    fv.className()

And groups except the skeletons (which is called as side groups) are

    fv.getSidegroups()
    
Note that the side groups obtained here are the groups around the skeleton, thus the site of each side group is also provided according to the recommended nomenclature of flavonoids [1]. All identified side groups with their connection information can be printed out via

    fv.getallSidegroups()
    
According to [2] and [3], a side group library with total 129 groups is constructed, which can cover most frequently occurred side groups in flavonoids. The structures of all side groups can be generated via

    fv.generateSidegroups(path)

Thus all groups can be found in *path*. If *path* is empty, the groups can be found under *sidegroups* folder locating under current working directory. **_If there exists any side group unidentified, a ValueError is raised._** To solve this, new side group can be added to the library by

    fv.addsidegroup(name,smiles)
    
Once the side group is added to the library, classifying flavonoid and matching side groups should be re-implemented, using `fv = Flavonoids(s)`.
	
In addition, basic properties such as molecular weight, exact molecular weight, formula can also be obtained, and the image can be generated.

**Classes of flavonoids**
The classes of flavonoids are summarized in **flavonoids_skeleton.pdf** provided in current package. The name of each class is collected from [2], [3] and references in the ..pdf file. If the input molecular identifier is not matched with the skeleton of defined class, `None` will be output. Therefore, some special flavonoids like *Mbamichalcone* can not be identified.

**Prerequisites**
Python 2.7

* This module is not tested on Python 3. If Python 3 is installed to run the module, *2to3* can be adopted to transform the modules *mol_interpreter*, *sidechain_props* and *flavonoid* into valid Python 3 program.

**Contact**
Nai-ping Dong

naiping.dong@hotmail.com

**References**  
[1] Rauter AP. [Nomenclature of Flavonoids (IUPAC Provisional Recommendations)](http://iupac.org/recommendation/nomenclature-of-flavonoids/).

[2] Harborne JB, Baxter H. The Handbook of Natural Flavonoids (2 Volume Set). Wiley, 1999.

[3] Buckingham J, Munasinghe VRN. Dictionary of Flavonoids. CRC Press, 2015.
