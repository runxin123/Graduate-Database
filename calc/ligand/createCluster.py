from __future__ import division, unicode_literals, print_function

"""
Wulff construction to create the nanoparticle
"""
import os
import sys
import itertools
from fractions import gcd
from functools import reduce
from collections import defaultdict
import logging
import copy
import numpy as np
import pymatgen as mp
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure, Molecule
from pymatgen.core.surface import get_symmetrically_distinct_miller_indices
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.coord import in_coord_list
from pymatgen.ext.matproj import MPRester
# from mpinterfaces import get_struct_from_mp

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(levelname)s:%(name)s:%(message)s')
sh = logging.StreamHandler(stream=sys.stdout)
sh.setFormatter(formatter)
logger.addHandler(sh)


class Nanoparticle(Molecule):
    """
    Construct nanoparticle using wulff construction
    """
    def __init__(self, structure, rmax=15, hkl_family=[(1,0,0), (1,1,1)],
                  surface_energies=[28,25], scell=None):
        self.structure = structure
        self.rmax = rmax
        self.hkl_family = hkl_family
        self.surface_energies = surface_energies
        if scell is None:
            ncella = int(np.ceil(2*rmax/structure.lattice.a))
            ncellb = int(np.ceil(2*rmax/structure.lattice.b))
            ncellc = int(np.ceil(2*rmax/structure.lattice.c))
            self.scell = [ncella, ncellb, ncellc]
        else:
            self.scell = scell
        self.structure.make_supercell(self.scell)
        recp_lattice = self.structure.lattice.reciprocal_lattice_crystallographic
        self.recp_lattice = recp_lattice.scale(1)
        self.set_miller_family()        
        Molecule.__init__(self, [site.species_and_occu 
                                 for site in self.structure.sites],
                           self.structure.cart_coords, charge=0 )
        #self._sites = list(self._sites)

    def set_miller_family(self):
        """
        get all miller indices for the given maximum index
        get the list of indices that correspond to the given family
        of indices
        """
        recp_structure = Structure(self.recp_lattice, ["H"], [[0, 0, 0]])
        analyzer = SpacegroupAnalyzer(recp_structure, symprec=0.001)
        symm_ops = analyzer.get_symmetry_operations()

        max_index = max( max(m) for m in self.hkl_family)
        r = list(range(-max_index, max_index + 1))
        r.reverse()
        miller_indices = []
        self.all_equiv_millers = []
        self.all_surface_energies = []
        for miller in itertools.product(r, r, r):
            if any([i != 0 for i in miller]):
                d = abs(reduce(gcd, miller))
                miller_index = tuple([int(i / d) for i in miller]) # get the simplex miller index
                for op in symm_ops:
                    for i, u_miller in enumerate(self.hkl_family):
                        if in_coord_list(u_miller, op.operate(miller_index)):
                            self.all_equiv_millers.append(miller_index)
                            self.all_surface_energies.append(self.surface_energies[i])
    def get_normals(self):
        """
        get the normal to the plane (h,k,l)
        """
        normals = []
        for hkl in self.all_equiv_millers:
            normal = self.recp_lattice.matrix[0,:]*hkl[0] + \
              self.recp_lattice.matrix[1,:]*hkl[1] + \
              self.recp_lattice.matrix[2,:]*hkl[2] 
            normals.append(normal/np.linalg.norm(normal))
        return normals

    def get_centered_molecule(self):
        center = self.center_of_mass
        new_coords = np.array(self.cart_coords) - center
        return Molecule(self.species_and_occu, new_coords,
                        charge=self._charge, 
                        spin_multiplicity=self._spin_multiplicity,
                        site_properties=self.site_properties)
    def create(self):
        """
        creates the nanoparticle by chopping of the corners normal to the
        specified surfaces.
        the distance to the surface from the center of the particel =
        normalized surface energy * max radius
        """
        mol = self.get_centered_molecule()    
        normalized_surface_energies = \
          np.array(self.all_surface_energies)/ float(max(self.all_surface_energies))
        surface_normals = self.get_normals()
        remove_sites = []
        for i, site in enumerate(mol):
            for j, normal in enumerate(surface_normals):
                n = np.array(normal)
                n = n/np.linalg.norm(n)
                if np.dot(site.coords, n) + self.rmax*normalized_surface_energies[j] <= 0 : # 核心判断依据
                    remove_sites.append(i)
                    break
        self.remove_sites(remove_sites)
        #new_sites = [site for k, site in enumerate(mol) if k not in remove_sites]
        #return Molecule.from_sites(new_sites)
    def getCoords(self):
        mol =self.get_centered_molecule()
        mol1 = copy.deepcopy(mol)
        return mol1


def mkdirs(path):
    isExists = os.path.exists(path)
    if not isExists:
        os.makedirs(path)
    else:
        print('dir was exists')
    pass
def setCell(filename,cell=[[20,0,0],[0,20,0],[0,0,20]]):
    lines = []
    with open(filename,'r') as f:
        for line in f:
            lines.append(line)
        lines[2] = str(cell[0][0])+'  '+str(cell[0][1])+'  '+str(cell[0][2])+'\n'
        lines[3] = str(cell[1][0])+'  '+str(cell[1][1])+'  '+str(cell[1][2])+'\n'
        lines[4] = str(cell[2][0])+'  '+str(cell[2][1])+'  '+str(cell[2][2])+'\n'
    s = ''.join(lines)
    with open(filename,'w') as f:
        f.write(s)
    pass
def createCluster(dirpath = './Database/dimercluster/',MPI_key='ccxfUSVCIvKrMF9W',symbols=['PbS'],surface_families=[(1,0,0),(1,1,1)],surface_energies=[28,25],rmax=15,cell=[[20,0,0],[0,20,0],[0,0,20]]):
    structures = []
    formulas = []
    with MPRester(MPI_key) as m:
        for i in symbols:
            data = m.get_data(i)
            x = {}
            for d in data:
                mpid = str(d['material_id'])
                x[mpid] = d['e_above_hull']
                structure = m.get_structure_by_material_id(mpid)
                structures.append(structure)
                formulas.append(i)
    # print(formulas)
    print(structures)
    for j,struct in enumerate(structures):
        sa = SpacegroupAnalyzer(structures[j])
        structure_conventional = sa.get_conventional_standard_structure()

        nanoparticle = Nanoparticle(structure_conventional, rmax=rmax,
                                hkl_family=surface_families,
                                surface_energies=surface_energies)
        surface_families_str = ''
        surface_energies_str = ''
        for m in surface_energies:
            surface_energies_str = surface_energies_str+str(m)+'_'
        for n in surface_families:
            surface_families_str = surface_families_str+str(n[0])+str(n[1])+str(n[2])+'_'
        fullpath = dirpath+str(formulas[j])+'/'+surface_energies_str+'/'+surface_families_str+'/'+str(j)+'/'
        print(fullpath)
        mkdirs(fullpath)
        #create nanoparticle 
        nanoparticle.create()

        print(nanoparticle.getCoords().cart_coords)
        ## get the coords

        nanoparticle.to(fmt='VASP', filename=fullpath+str(formulas[j])+str(j)+'.VASP')
        setCell(filename=fullpath+str(formulas[j])+str(j)+'.VASP',cell=cell)
    pass

            
if __name__ == '__main__':
    #nanopartcle settings
    #max radius in angstroms
    MPI_key = "ccxfUSVCIvKrMF9W"

    rmax = 15 
    # symbols=['PdPt','AuPd','AuPt','AuRh','AgPt','AgPd','AuPd']
    # metal = ["Pt","Pd","Au","Ag","Rh","Ru","Ir","Os"]
    metal = ["Zn","S"]
    # symbols=['PdPt']
    dimer_metal = []
    for i in metal:
        for x in metal:
            if x!=i:
                if metal.index(x)>metal.index(i):
                    dimer_metal.append(i+x)

    #surface families to be chopped off
    for i in range(1,10):
        surface_families = [(i,0,0), (i,1,1),(i,1,0)]
        for j in range(10,30):
            surface_energies = [j-5,j-1,j-5]
            createCluster(symbols=dimer_metal,surface_families=surface_families,surface_energies=surface_energies,rmax=6)

    # surface_families = [(1,0,0), (1,1,1)]
    #could be in any units, will be normalized 
    # surface_energies = [28,25]  

    # createCluster(symbols=dimer_metal,surface_families=surface_families,surface_energies=surface_energies,rmax=6)
    #caution: set the structure wrt which the the miller indices are specified
    #use your own API key



    #primitve ---> conventional cell
