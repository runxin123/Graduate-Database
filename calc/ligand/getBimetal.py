from ase import Atoms
from ase.build import fcc111,add_adsorbate,bulk,bcc100,bcc110,bcc111,fcc110,fcc100,fcc211,hcp0001
from ase.io import read,write,iread
# import ase.build as bd
from ase.cluster.cubic import FaceCenteredCubic
# from catkit.gen.adsorption import AdsorptionSites
# import os
from scipy.spatial import ConvexHull,Delaunay
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 
from ase.visualize import view
from ase.constraints import FixAtoms,FixedPlane
import numpy as np
from ase.collections import g2
from ase.build import molecule,surface,make_supercell
import os
import copy
from scipy.spatial import ConvexHull,Delaunay
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

def writePOSCAR(filename,atoms,cell,title,selective=False,fix_index=[]):
	with open(filename,'w') as f:
		f.write(str(title)+'\n')
		f.write(str(1.0)+'\n')
		for i in cell:
			f.write(' '+str(i[0])+'  '+str(i[1])+'  '+str(i[2])+'\n')
		for j in atoms:
			f.write('  '+str(j))
			pass
		f.write('\n')
		for m in atoms:
			f.write('  '+str(len(atoms[m])))
			pass

		f.write('\n')
		# f.write(str('Selective')+'\n')
		if selective==False:
			f.write(str('Cartesian')+'\n')
		for k in atoms:
			for n in atoms[k]:
				f.write(' '+str(n[0])+'  '+str(n[1])+'  '+str(n[2])+'\n')
			pass
def cutPlane(atoms,indice=(1,1,1),layer=4,fix_layer=0,vacuum=10,size=(1,1,1)):
	tmp_surface = surface(atoms,indice,layer)
	tmp_surface = tmp_surface*size
	tmp_surface.center(vacuum=vacuum,axis=2)
	return tmp_surface

def getLattice(filename,fmt='cif'):
	# print(filename)
	tmp_lattice = read(filename,format=fmt)
	return tmp_lattice
def getALattice(atoms):
	return atoms

def func4(one_list):
    '''
    使用排序的方法
    '''
    result_list=[]
    temp_list=sorted(one_list)
    i=0
    while i<len(temp_list):
        if temp_list[i] not in result_list:
            result_list.append(temp_list[i])
        else:
            i+=1
    return result_list

def getInformation(slab):
	atomss = {}
	atoms = [] 
	title = ''
	tit = ''
	for i in slab:
		# print(i)
		# print(i.get_raw('symbol'))
		k = i.get_raw('symbol')
		if k in atoms:
			pass
		else:
			atoms.append(k)

	for k in atoms:
		atomss[k] = []
	for i in slab:
		if i.get_raw('symbol') in atoms:
			atomss[i.get_raw('symbol')].append(i.get_raw('position').tolist())
	for k in atoms:
		title  = title+str(k)+str(len(atomss[k]))
		tit = tit +str(k)
		pass
	atoms_t = copy.deepcopy(atoms)
	atomss_t = copy.deepcopy(atomss)
	title_t = copy.deepcopy(title)
	tit_t = copy.deepcopy(tit)
	# print(str(atoms_t))
	# print(str(atomss_t))
	# print(str(title_t))
	return atoms_t,atomss_t,title_t,tit_t
def getDimerMetalPath(filepath = './Database/'):
	list_dir = os.listdir(filepath)
	list_dir1 = []
	# print(list_dir)
	for i in list_dir:
		filepath1 = filepath+i+'/'
		for j in os.listdir(filepath1):
			list_dir1.append(filepath1+j)
	return list_dir1
def add_adsorbate(surf,atom='CO',positions=[],d=1.8,adsorb_vector=[0,0,1],bond=1.1459):
	positions_t = copy.deepcopy(positions)
	positions[2] = positions[2]+d
	positions_t = np.array(positions)+np.array(adsorb_vector)*1.1459
	positions_t = positions_t.tolist()
	pos = [positions,positions_t]
	# pos = [positions]
	# print(positions)

	surf1 = copy.deepcopy(surf)
	surf1.extend(Atoms(atom,positions=pos))
	return surf1
	pass

def add_adsorbate1(surf,atom='CO',positions=[],d=1.8,adsorb_vector=[0,0,1],bond=1.1459):
	positions_t = copy.deepcopy(positions)
	positions[2] = positions[2]+d
	positions_t = np.array(positions)+np.array(adsorb_vector)*bond
	positions_t = positions_t.tolist()
	pos = [positions,positions_t]
	# pos = [positions]
	# print(positions)

	surf1 = copy.deepcopy(surf)
	surf.extend(Atoms(atom,positions=pos))
	return surf
	pass

def mkdirs(path):
	isExists = os.path.exists(path)
	if not isExists:
		os.makedirs(path)
	else:
		print('dir was exists')
	pass

def getXYcoord(pts):
	Zcoord = []
	XYcount = []
	XYcoord = {}
	for i in pts:
		if i[2] in Zcoord:
			pass
		else:
			Zcoord.append(i[2])
	for n in Zcoord:
		XYcoord[Zcoord.index(n)]=[]

	for j in pts:
		if j[2] in Zcoord:
			XYcoord[Zcoord.index(j[2])].append([j[0],j[1]])
			pass
		pass
	for m in XYcoord:
		XYcount.append(len(XYcoord[m]))
	return Zcoord,XYcoord,XYcount

def triMesh(pots=[]):
	tri = Delaunay(pots,qhull_options='Qx')
	bri = []
	trid = []
	top2 = []
	tri_result = []
	facet = []
	top = []
	pots = np.array(pots)
	# print('coplanar = '+str(tri.coplanar))
	for simple in tri.simplices:
		print(simple)
		for i in range(0,3):
			tmp_top1 = []
			for j in range(0,2):
				tmp_top1.append(pots[simple,j][i])
			top.append(tmp_top1)
		bri.append((np.array(top[-1])+np.array(top[-2]))/2)
		bri.append((np.array(top[-1])+np.array(top[-3]))/2)
		bri.append((np.array(top[-2])+np.array(top[-3]))/2)
		trid.append((np.array(top[-1])+np.array(top[-2])+np.array(top[-3]))/3)
		tri_result.append(top)
	return np.array(top),bri,trid

def setBriSites(surf,atom='CO',d=2.0,layer=2,bond=1.1459,indice=(1,1,1),index=0):
	pots = surf.get_positions().tolist()
	XYcoord = getXYcoord(pots)[1]
	Zcoord = getXYcoord(pots)[0]
	i1 = 0
	pts = []
	database_dir = './Database/dimermetal/adsorption/'
	atomt = getInformation(surf)
	dir_t = database_dir+'bri'+'/'+str(atomt[2])+'/'+str(indice[0])+str(indice[1])+str(indice[2])+'/'+str(d)+'/'
	s5 = add_adsorbate(surf=surf,atom=atom,positions=[0,0,0],d=d)
	mkdirs(dir_t)
	for j in range(1,layer+1):
		pts = pts+XYcoord[Zcoord.index(Zcoord[-j])]
	pts = func4(pts)
	sitess = triMesh(pots)[1]
	for i in sitess:

		i = i.tolist()
		i.append(Zcoord[-1])
		s4 = add_adsorbate(surf=surf,atom=atom,positions=i,d=d)
		atomt4 = getInformation(s4)
		atomt = getInformation(surf)	
		# print(atomt)
		mkdirs(dir_t+'/'+str(i1)+'/')
		writePOSCAR(dir_t+'/'+str(i1)+'/'+str(atomt4[2])+'.vasp',atoms=atomt4[1],cell=s4.get_cell(),title=atomt4[2])
		i1 = i1+1
	surf_static = copy.deepcopy(surf)	
	for j in sitess:
		j = j.tolist()
		# ax.scatter(j[0],j[1],c='b')
		j.append(Zcoord[-1])
		s5 = add_adsorbate1(surf=surf_static,atom=atom,positions=j,d=d)
	atomt5 = getInformation(s5)
	# writePOSCAR(dir_t+str(atomt5[2])+str(index)+'/'+str(atomt5[2])+'.vasp',atoms=atomt5[1],cell=s5.get_cell(),title=atomt5[2])

def setTopSites(surf,atom='CO',d=2.0,layer=2,bond=1.1459,indice=(1,1,1)):
	pots = surf.get_positions().tolist()
	XYcoord = getXYcoord(pots)[1]
	Zcoord = getXYcoord(pots)[0]
	pts = []

	database_dir = './Database/dimermetal/adsorption/'
	atomt = getInformation(surf)
	dir_t = database_dir+'top'+'/'+str(atomt[2])+'/'+str(indice[0])+str(indice[1])+str(indice[2])+'/'+str(d)+'/'
	s5 = add_adsorbate(surf=surf,atom=atom,positions=[0,0,0],d=d)
	mkdirs(dir_t)
	for j in XYcoord:
		pts = pts + XYcoord[j]
	pts = func4(pts)
	sitess = triMesh(pts)[0]
	for i in sitess:
		i = i.tolist()
		i.append(Zcoord[-1])
		s4 = add_adsorbate(surf=surf,atom=atom,positions=i,d=d)
		atomt4 = getInformation(s4)
		print(atomt)
		writePOSCAR(dir_t+'/'+str(atomt4[2])+str(i)+'.vasp',atoms=atomt4[1],cell=s4.get_cell(),title=atomt4[2])
	surf_static = copy.deepcopy(surf)
	for j in sitess:
		j = j.tolist()
		j.append(Zcoord[-1])
		s5 = add_adsorbate1(surf=surf_static,atom=atom,positions=j,d=d)
	atomt5 = getInformation(s5)
	writePOSCAR(dir_t+'/'+str(atomt5[2])+str(i)+'.vasp',atoms=atomt5[1],cell=s5.get_cell(),title=atomt5[2])

def setHollowSites(surf,atom='CO',d=2.0,layer=2,bond=1.1459,indice=(1,1,1)):
	pots = surf.get_positions().tolist()
	XYcoord = getXYcoord(pots)[1]
	Zcoord = getXYcoord(pots)[0]
	pts = []
	for i in pots:
		i.pop()
	database_dir = './Database/dimermetal/adsorption/'
	atomt = getInformation(surf)
	dir_t = database_dir+'hollow'+'/'+str(atomt[2])+'/'+str(indice[0])+str(indice[1])+str(indice[2])+'/'+str(d)+'/'
	s5 = add_adsorbate(surf=surf,atom=atom,positions=[0,0,0],d=d)
	mkdirs(dir_t)
	for j in XYcoord:
		pts = pts + XYcoord[j]
	pts = func4(pts)
	print('pts = '+str(len(pts)))
	sitess = triMesh(pots)[2]
	for i in sitess:
		# ax.scatter(i[0],i[1],c='g')
		i = i.tolist()
		i.append(Zcoord[-1])
		s4 = add_adsorbate(surf=surf,atom=atom,positions=i,d=d)
		atomt4 = getInformation(s4)
		# print(atomt)
		writePOSCAR(dir_t+'/'+str(atomt4[2])+str(i)+'.vasp',atoms=atomt4[1],cell=s4.get_cell(),title=atomt4[2])
	surf_static = copy.deepcopy(surf)
	for j in sitess:
		j = j.tolist()
		ax.scatter(j[0],j[1],c='g')
		j.append(Zcoord[-1])
		s5 = add_adsorbate1(surf=surf_static,atom=atom,positions=j,d=d)
	atomt5 = getInformation(s5)
	writePOSCAR(dir_t+'/'+str(atomt5[2])+str(i)+'.vasp',atoms=atomt5[1],cell=s5.get_cell(),title=atomt5[2])
# AuPt = Atoms('AuPt',cell=[[3.1763689518,0.0000000000,0.0000000000],[0.0000000000,3.1763689518,0.0000000000],[0.0000000000,0.0000000000,3.1763689518]],positions=[[0.000000000,0.000000000,0.000000000],
# [1.588181257,1.588181257,1.588181257]],pbc=True)


# file_dir = './Database/dimermetal/VASP/'
# filepath  = './Database/dimermetal/cif/'
# filedir  = getDimerMetalPath(filepath=filepath)
# indice = (4,1,1)
# size = (2,2,1)
# for i in filedir:
# 	slab = getLattice(i,fmt='cif')
# 	# print(tmp_AuPt)
# 	# tmp_AuPt = tmp_AuPt*(4,4,1)
# 	# s3 = surface(tmp_AuPt,(2,1,1),4)
# 	# s3.center(vacuum=10,axis=2)
# 	s3 = cutPlane(atoms=slab,indice=indice,layer=4,size=size)
# 	# atomss = {}
# 	# atoms = [] 
# 	atomt = getInformation(s3)
# 	mkdirs(file_dir+str(atomt[3])+'/'+str(indice[0])+str(indice[1])+str(indice[2])+'/'+str(size[0])+str(size[1])+str(size[2])+'/'+str(str(atomt[2]))+'/')
# 	writePOSCAR(file_dir+str(atomt[3])+'/'+str(indice[0])+str(indice[1])+str(indice[2])+'/'+str(size[0])+str(size[1])+str(size[2])+'/'+str(atomt[2])+'/'+str(atomt[2])+'.vasp',atoms=atomt[1],cell=s3.get_cell(),title=atomt[2])

# AuPt.set_chemical_symbols('PtAu')
# s4 = surface(AuPt,(1,1,1),4)

##########test the sites###############
# indice  = (1,1,1)
# slab = getLattice('./Database/dimermetal/cif/PtAu/Au1Pt1-183337.cif',fmt='cif')
# # print(tmp_AuPt)
# # tmp_AuPt = tmp_AuPt*(4,4,1)
# # s3 = surface(tmp_AuPt,(2,1,1),4)
# # s3.center(vacuum=10,axis=2)
# s3 = cutPlane(atoms=slab,indice=indice,layer=4,size=(2,2,1))
# # atomss = {}
# # atoms = [] 
# atomt = getInformation(s3)
# # mkdirs(file_dir+str(indice)+'/')
# # writePOSCAR(file_dir+str(indice)+'/'+str(atomt[2])+'.vasp',atoms=atomt[1],cell=s3.get_cell(),title=atomt[2])






#######test the sites###############
# points = s3.get_positions().tolist()
# for i in points:
# 	i.pop()
# bri = triMesh(points)[1]
# top = triMesh(points)[0]
# hollow = triMesh(points)[2]
# print(s3.get_positions())


# s4 = make_supercell(AuPt,[[3,0,0],[0,4,0],[0,0,4]])
# s4 = surface(s4,(1,1,0),4)
# s4.center(vacuum=10,axis=2)

# print(s4.get_positions())
# surfaces = [(1, 0, 0), (1, 1, 0), (2, 1, 0)]
# layers = [4, 4, 4]
# lc = 3.61000
# atoms = FaceCenteredCubic('PtAu', surfaces, layers, latticeconstant=lc)
# p = atoms.get_positions()

# write('./PtAu_slab.vasp',s3*(1,1,1))
# modifyPOSCAR('.','PtAu_slab.vasp')

# view(s4)
# print(tmp_AuPt)


# fig = plt.figure()
# ax= Axes3D(fig)



# print(pts)


file_dir = './Database/dimermetal/VASP/'
filepath  = './Database/dimermetal/cif/'
filedir  = getDimerMetalPath(filepath=filepath)
indice = (1,1,1)
i1 = 0
for i in filedir:
	slab_t = getLattice(i,fmt='cif')
	print(slab_t)
	surf = copy.deepcopy(cutPlane(slab_t,indice,4,size=(2,2,1)))
	# print(surf)
	# atom_t = (getInformation(surf))
	setBriSites(surf=surf,atom='CO',d=1.0,layer=4,indice=indice,index=i1)
	i1 = i1+1


# setBriSites(s3,atom='CO',d=2.0,layer=2,bond=1.1459,indice=(1,1,1))

# for i in top:
# 	ax.scatter(i[0],i[1],c='r')
# for j in bri:
# 	ax.scatter(j[0],j[1],c='g')


# for k in hollow:
# 	ax.scatter(k[0],k[1],c='b')


# plt.show()
# view(s3)
