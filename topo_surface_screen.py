'''
@author: wenhao sun
'''

from ScreenSurf import ScreenSurf
from pymatgen.io.vasp import Poscar
from pymatgen.io.cif import CifParser
from pymatgen.matproj.rest import MPRester
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

# Structures can be directly imported from the Materials Project, using the Materials Project API

# Here are some examples from the paper
# TiO2-Anatase:   'mp-390'
# Cu2O-Cuprite:   'mp-361'
# Al2O3-Corundum: 'mp-1143'

with MPRester("YOUR_API_KEY") as m:
        bulkstruct = m.get_structure_by_material_id('mp-1143')

# Initiate the surface screening algorithm
SS=ScreenSurf(bulkstruct)

index=SS.ScreenSurfaces(natomsinsphere=20,keep=1,samespeconly=True,ignore=[])
# ScreenSurfaces uses the following parameters:
#    natomsinsphere = Increases radius of sphere until there are #natomsinsphere in the sphere. For simple structures this can be as low as 20. 
#    keep = What percentage of the histogram to keep up to. Using keep=0.8 means we only keep the most common 80% of planes. 
#    samespeconly = This is whether or not to search for only atoms of the same element in the sphere. 
#    ignore = these are atoms to ignore. These often include small atoms on complex anions (e.g., O on CO_3^2-), or other moieties, like hydrogen. 


#Structures can also be read from CIF files in pymatgen
parser = CifParser("AnataseTiO2.cif")
Anatase=parser.get_structures()[0]
SS=ScreenSurf(Anatase)
index=SS.ScreenSurfaces(natomsinsphere=20,keep=1,samespeconly=True,ignore=[])


#Structures can also be read from VASP POSCAR files in pymatgen
LiFePO4=Poscar.from_file("Poscar_LiFePO4").structure
SS=ScreenSurf(LiFePO4)
index=SS.ScreenSurfaces(natomsinsphere=100,keep=1,samespeconly=True,ignore=['O'])
