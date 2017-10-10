# Surface-Screening
A Topological Screening Heuristic for Low-Energy, High-Index Surfaces

Authors: Wenhao Sun1,2, Gerbrand Ceder1,2,3

1 Department of Materials Science and Engineering, Massachusetts Institute of Technology, Cambridge, MA 02139 USA.

2 Materials Sciences Division, Lawrence Berkeley National Laboratory, Berkeley, California 94720, USA 

3 Department of Materials Science and Engineering, UC Berkeley, Berkeley, California 94720, USA

Supplementary Information

In this repository are four files: 

•	ScreenSurf.py

•	topo_surface_screen.py

•	AnataseTiO2.cif

•	Poscar_LiFePO4


To run the surface screening algorithm, execute and topo_surface_screen.py in your favorite Python development environment. Note that executing this code requires the PyMatGen (Python Materials Genomics) package, which can be attained at http://pymatgen.org/#getting-pymatgen


The algorithm intakes pymatgen structure objects. We have provided 3 examples that create PyMatGen structure objects from the 1) Materials Project, 2) a CIF file, 3) a VASP POSCAR file. If you wish to pull structures from the Materials Project, attain a (free) API key at https://materialsproject.org/open
