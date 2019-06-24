from ase import Atoms
from ase.build import fcc111,add_adsorbate,bulk,bcc100,bcc110,bcc111,fcc110,fcc100,fcc211,hcp0001
from ase.io import read,write
# import ase.build as bd
from ase.constraints import FixAtoms,FixedPlane
import numpy as np
from ase.collections import g2
from ase.build import molecule
import os
import copy
class autogen(object):
	"""docstring for autogen"""
	'''
	auto genarator the slab and adsorbate
	'''
	# slab_name = ['Pt','Pd','Rh','Ru','Ag','Au','Pb','Cu','Fe','Ni','Co','Ir','Mo']
	slab_name = ['Rh','Ru','Ir']
	slabs = [] # build slabs array()
	s_types = {} # slab types 
	# adsorbate_name = ['CO','CO2','NO','H2','H2O','NO2','O','H','OH','N2O','NH2','CH3','NH']
	adsorbate_name = ['CO','NO','OH']
	adsorbates = {}
	vacuum = 15.0
	lattice_indice = {'fcc':{'111':['ontop','fcc','hcp','bridge'],'110':['ontop','longbridge','shortbridge','hollow'],'100':['ontop','bridge','hollow']},
	'bcc':{'111':['ontop','hollow'],'110':['ontop','longbridge','shortbridge','hollow'],'100':['ontop','bridge','hollow']},'hcp':{'0001':['ontop','bridge','fcc','hcp']}}
	# {fcc:{111:ontop}} # lattice indice adsorption site
	slab_lattice_constant = {'Pt':3.92,'Pd':3.89,'Rh':3.803,'Ru':2.706,'Ag':4.086,'Au':4.0781,'Pb':4.95,'Fe':2.866,'Ni':3.524,'Cu':3.615,'Co':2.505,'Zn':2.6648,'Ir':3.840,'Mo':3.147,'Nb':3.301,'Ti':2.951,'V':3.024,'Y':3.647,'Zr':3.231,'Hf':3.196,'Os':2.734} #{element:lattice_constant}

	def b_adsorbates(self,nameList=[]):
		# default adsorbates was CO
		# setting the building the adsorbates from the ase collections
		adsorbates_t = {}
		if len(nameList) !=0:
			for i in nameList:
				adsorbates_t[i] = molecule(i)
		else:
			pass
		return adsorbates_t

	def b_slab(self,name='Pt',nameList=[],size_t=3,types={},vacuum=15):
		# size = (4,4,4)
		# types =fcc111
		if len(types)!=0:
			for t in types:
				for j in types[t]:
					self.s_types[str(t)+str(j)] = []
					for m in types[t][j]:
						self.s_types[str(t)+str(j)].append(m)
		else:
			pass

		if (len(self.s_types)!=0) and (len(nameList)!=0):
			for i in nameList:
				element_nm = i
				for j in self.s_types:
					nm = j
					# print(element_nm)
					slabs_t = {}
					slabs_t['slab'] = build_slab(element_nm,(size_t,size_t,4),self.slab_lattice_constant[element_nm],vacuum,nm)
					slabs_t['element'] = element_nm
					slabs_t['types'] = nm
					print(slabs_t)
					self.slabs.append(slabs_t)
			pass
		# print(self.slabs)
		# print(self.s_types)
	

	def add_ligand_list(self,name):
		self.adsorbate_name.append(name)

	def bind_Ligand_slab(self,adsorbate_t):
		# todo bind ligand and slab
		for ads in adsorbate_t:
			for i in self.slabs:
				slabdir =i['element']+'/slab/'+i['types']+'/'
				slabdir_t = './Database/'+slabdir
				mkdirs(slabdir_t)
				# print(slabdir_t)
				slab_b = copy.deepcopy(i['slab'])
				t_slab = fixlayers(slab_b,2)
				write(str(slabdir_t+'POSCAR'),t_slab)
				modifyPOSCAR(slabdir_t,'POSCAR')
				for j in self.s_types[i['types']]:
					ads_dir = i['element']+'/'+ads+'/'+i['types']+'/'+j
					filedir_t = './Database'+'/'+ads_dir
					mkdirs(filedir_t)
					print('adsorbate_t[ads] = '+str(adsorbate_t[ads]))
					slab_t = copy.deepcopy(i['slab'])
					d_slab = fixlayers(slab_t,2)
					adsorption(filedir_t,'POSCAR',d_slab,adsorbate_t[ads],j)
					modifyPOSCAR(filedir_t,'POSCAR')
	# get slab_name
	def get_slab_name(self):
		return self.slab_name

	# get lattice_indice
	def get_lattice_indice(self):
		return self.lattice_indice

	# get the adsorbates	
	def get_adsorbates(self,nameList=[]):
		if len(nameList)==0:
			self.adsorbates = self.b_adsorbates(self.adsorbate_name)
		else:
			self.adsorbates = self.b_adsorbates(nameList)
		return self.adsorbates
		
def test():
	# d = 1.1
	# adsorbate = Atoms('NO')
	# print(adsorbate)
	# # co.set_cell(2*np.identity(3))
	# # co.set_positions([(0,0,0),(0,0,1.1)])

	# # slab = fcc111('Al', size=(2, 2, 3))
	# # add_adsorbate(slab, 'Au', 1.7, 'ontop')
	# # slab.center(axis=2, vacuum=4.0)

	# slab = fcc111('Pb',(2,2,3), a=3.7,vacuum=15.0)

	# mask = [atom.tag > 1 for atom in slab]
	# #print(mask)
	# fixlayers = FixAtoms(mask=mask)
	# # plane = FixedPlane(-1, (1, 0, 0))
	# slab.set_constraint(fixlayers)

	# # adsorbate[1].z = 1.1
	# # add_adsorbate(slab,adsorbate,1.9,'ontop')
	# write('POSCAR',slab)
	# pass
	c = autogen()
	c.b_slab('Pt',c.get_slab_name(),3,c.get_lattice_indice(),15)
	c.bind_Ligand_slab(c.get_adsorbates()) # add the dir and get POSCAR 


def build_slab(symbol = 'Pt',size=(3,3,4),a=3.5,vacuum=15,types='fcc111',c=0):
	c= (8/3)*a
	# default the c was 8/3 a
	print('a = '+str(a))
	if types == 'fcc111':
		slab = fcc111(symbol,size,a,vacuum)
		return slab
	elif types == 'fcc110':
		slab = fcc110(symbol,size,a,vacuum)
		return slab
	elif types == 'fcc100':
		slab = fcc100(symbol,size,a,vacuum)
		return slab
	elif types == 'bcc100':
		slab = bcc100(symbol,size,a,vacuum)
		return slab
	elif types == 'hcp10m10':
		slab = hcp10m10(symbol,size,a,vacuum)
		return slab
	elif types == 'diamond100':
		slab = diamond100(symbol,size,a,vacuum)
		return slab
	elif types == 'fcc211':
		slab = fcc211(symbol,size,a,vacuum)
		return slab
	elif types == 'bcc110':
		slab = bcc110(symbol,size,a,vacuum)
		return slab
	elif types == 'bcc111':
		slab = bcc111(symbol,size,a,vacuum)
		return slab
	elif types == 'hcp0001':
		slab = hcp0001(symbol,size,a,c,vacuum)
		return slab
	elif types == 'diamond111':
		slab = diamond111(symbol,size,a,vacuum)
		return slab
	else:
		pass


def adsorption(path,filename,slab,adsorbate,position = 'ontop',height=1.5,size=(1,1,1)):
	# size supercell size
	'''
	@ path was file that create
	@ filename was file that create
	'''
	# add_adsorbate will change the slab
	slab_t = slab
	# print(slab)
	# print(adsorbate)
	add_adsorbate(slab_t,adsorbate,height,position)
	write(str(path+'/'+filename),slab_t*size)
	pass

def mkdirs(path):
	isExists = os.path.exists(path)
	if not isExists:
		os.makedirs(path)
	else:
		print('dir was exists')
	pass

def modifyPOSCAR(filedir,filename):
	# modify POSCAR with element 
	lines = []
	with open(filedir+'/'+filename,'r') as f:
		for line in f:
			lines.append(line)
	line_t =lines[0]
	line_str = "  "+lines[0]
	lines.insert(5,line_str)
	s = ''.join(lines)
	with open(filedir+'/'+filename,'w') as f:
		f.write(s)

def fixlayers(slab,layer):
	mask = [atom.tag > layer for atom in slab]
	#print(mask)
	fixlayers = FixAtoms(mask=mask)
	# plane = FixedPlane(-1, (1, 0, 0))
	slab.set_constraint(fixlayers)
	return slab



test()