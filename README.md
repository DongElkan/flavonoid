Module "flavonoid" is used to identify the class of flavonoids according to the input molecular identifiers and facilitate the interpretation of mass spectra derived from the molecules, with the help of indigo toolkit.

#**Usage**
    
    from flavonoid import Flavonoids
    fv = Flavonoids(s)

Where "s" is an string of molecular identifier. Currently, canonical SMILES, InChI, and chemical file format ..mol and ..sdf are supported. Once Flavonoids object "fv" is constructed, properties of the molecule can be accessed by the methods built.
For example, the skeleton information of the input molecule can be obtained via

    fv.className()

And groups except the skeletons (which is called as side groups) are

    fv.getSidegroups()
    
Note that the side groups obtained here are the groups around the skeleton, thus the site of each side group is also provided according to the recommended nomenclature of flavonoids [1]. All side groups with their connecting network can be extracted via

    fv.getallSidegroups()
    
Need to continue......
