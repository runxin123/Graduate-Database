from ase import Atoms
from ase.build import fcc111,add_adsorbate,bulk
from ase.io import read,write
import numpy as np
def test():
	d = 1.1
	adsorbate = Atoms('NO')
	print(adsorbate)
	# co.set_cell(2*np.identity(3))
	# co.set_positions([(0,0,0),(0,0,1.1)])
	slab = fcc111('Pb',(2,2,3), a=3.7,vacuum=15.0)
	adsorbate[1].z = 1.1
	add_adsorbate(slab,adsorbate,1.9,'ontop')
	write('POSCAR',slab*(3,3,1))
	pass