# from catkit.pawprint import Fingerprinter
# from catkit.build import surface
# from ase.build import bulk
# import toppy as np

# bulk = bulk('Pd', cubic=True)
# bulk[3].symbol = 'Pt'

# slab = surface('Al', size=(2, 2, 3), a=3.8, vacuum=10)

# images = [bulk, slab]

# parameters = [
#     'atomic_topber',
#     'atomic_radius',
#     'atomic_volume',
# ]

# operations = [
#     'periodic_convolution',
#     ['periodic_convolution', {'d': 1}]
# ]

# fp = Fingerprinter(images)
# fingerprints = fp.get_fp(parameters, operations)

# print('#+attr_latex: :mode math :environment matrix')
# for fp in np.array(fingerprints):
#     fp_table = np.round(fp, 3)
#     print('|' + '|'.join(fp_table.astype(str)))


# from catkit import Gratoms
# from catkit.build import molecule, surface
# import toppy as np

# atoms0 = surface('Al', size=(2, 2, 3), a=3.8, vacuum=10)

# print(type(atoms0))
# print('\nAtomic topbers')

# topbers = atoms0.topbers

# print('#+attr_latex: :mode math :environment matrix')
# print('|' + '|'.join(topbers.astype(str)))

# print('\nConnectivity matrix:'.format())
# matrix = atoms0.connectivity
# table_matrix = np.array(matrix, dtype=str)

# print('#+attr_latex: :mode math :environment matrix')
# for row in table_matrix:
#     print('|' + '|'.join(row.astype(str)))

# atoms1 = atoms0.copy()
# del atoms1[0]

# print('\nAtomic topbers')
# topbers = atoms1.topbers

# print('#+attr_latex: :mode math :environment matrix')
# print('|' + '|'.join(topbers.astype(str)))

# print('\nConnectivity matrix:'.format())
# matrix = atoms1.connectivity
# table_matrix = np.array(matrix, dtype=str)

# print('#+attr_latex: :mode math :environment matrix')
# for row in table_matrix:
#     print('|' + '|'.join(row.astype(str)))

# print('\nAre isomorphs: {}'.format(atoms0.is_isomorph(atoms1)))

from catkit.gen.surface import SlabGenerator
from ase.build import bulk
from ase.visualize import view
import numpy as np
# Make a test slab
# atoms = bulk('Pt', 'fcc', cubic=True)
# atoms[2].symbol = 'Cu'

# gen = SlabGenerator(
#     atoms,
#     miller_index=(1, 1, 1),
#     layers=9,
#     fixed=5,
#     vacuum=4)

# terminations = gen.get_unique_terminations()

# images = []
# for i, t in etoperate(terminations):
#     images += [gen.get_slab(iterm=i)]

from ase.cluster.cubic import FaceCenteredCubic
from catkit.gen.adsorption import AdsorptionSites
from scipy.spatial import ConvexHull,Delaunay
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 
# from matplotlib.collections import PolyCollection
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
surfaces = [(1, 0, 0), (1, 1, 0), (1, 1, 1)]
layers = [4, 4, 4]
lc = 3.61000
atoms = FaceCenteredCubic('Pd', surfaces, layers, latticeconstant=lc)
p = atoms.get_positions()
facet = []
def minusList(list1,list2):
	lisindex = min(len(list1),len(list2))
	c = [list1[i] - list2[i] for i in range(lisindex)]
	return c

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
def getNormalVertax(face_array):
	'''
	face_array [[x1,y1,z1],[x2,y2,z2],[x3,y3,z3]]
	'''
	face = []
	face.append((face_array[0][1]-face_array[1][1])*(face_array[1][2]-face_array[2][2])-(face_array[0][2]-face_array[2][2])*(face_array[1][1]-face_array[2][1]))
	face.append((face_array[1][0]-face_array[2][0])*(face_array[0][2]-face_array[2][2])-(face_array[0][0]-face_array[1][0])*(face_array[1][2]-face_array[2][2]))
	face.append((face_array[0][0]-face_array[1][0])*(face_array[1][1]-face_array[2][1])-(face_array[0][1]-face_array[1][1])*(face_array[1][0]-face_array[2][0]))
	return face
def getTopCenter(top):
	x=0
	y=0
	z=0
	TopCenter = []
	for i in top:
		x = x+i[0]
		y = y+i[1]
		z = z+i[2]
	avg_x = x/len(top)
	avg_y = y/len(top)
	avg_z = z/len(top)
	print(TopCenter)
	TopCenter = [avg_x,avg_y,avg_z]
	return TopCenter

def getTopNormal(top,TopCenter=[0,0,0]):
	TopNormal = []
	for i in top:
		TopNormal.append(minusList(i,TopCenter))
		# print(TopNormal)
		pass
	return TopNormal

# view(atoms)
# coordinates, connectivity = gen.adsorption_sites(atoms)
# sites = AdsorptionSites(atoms)
# sites.plot('./adsorption-sites.png')
# for i in p:
print(len(p))
hull = ConvexHull(p,False,'Qx')
# using coplaner 
fig = plt.figure()
ax = Axes3D(fig)
top  = []
bri = []
top1 = []
trid = []
x = []
y = []
z = []
t_normal = []
# for simple2 in hull.coplanar:
# 	# tmp_top =[]
# 	# print(p[simple2,0])
# 	# print(p[simple2,1])
# 	# print(p[simple2,2])
# 	# x.append(p[simple2,0])
# 	# y.append(p[simple2,1])
# 	# z.append(p[simple2,2])

# 	# ax.scatter(p[simple2,0],p[simple2,1],p[simple2,2],c='g')
# 	for i in range(0,3):
# 		tmp_top = []
# 		for j in range(0,3):
# 			tmp_top.append(p[simple2,j][i])
# 		top.append(tmp_top)
# 		x.append(p[simple2,0][i])
# 		y.append(p[simple2,1][i])
# 		z.append(p[simple2,2][i])
# 	# bri.append((np.array(top[-1])+np.array(top[-2]))/2)
# 	# bri.append((np.array(top[-1])+np.array(top[-3]))/2)
# 	# bri.append((np.array(top[-3])+np.array(top[-2]))/2)
# 	# trid.append((np.array(top[-1])+np.array(top[-2])+np.array(top[-3]))/3)
# 	# top.append(p[simple2,0].tolist())
# 	# top.append(p[simple2,1].tolist())
# 	# top.append(p[simple2,2].tolist())


# print(top)
pointNumber = len(top)
points = np.zeros((pointNumber, 3))
points[:, 0] = np.array(x)
points[:, 1] = np.array(y)
points[:, 2] = np.array(z)
topn = np.array(top)
tri = Delaunay(p,qhull_options='Qx')
x1 = []
y1 = []
z1 = []
vertice = []
# s_top = set()
for simple in tri.convex_hull:
	top2 = []
	for i in range(0,3):

		tmp_top1 = []
		for j in range(0,3):
			tmp_top1.append(p[simple,j][i])
		top1.append(tmp_top1)
		# s_top.add(tmp_top1)
		top2.append(tmp_top1)
		# print(getNormalVertax(tmp_top1))
		x1.append(p[simple,0][i])
		y1.append(p[simple,1][i])
		z1.append(p[simple,2][i])
	bri.append((np.array(top1[-1])+np.array(top1[-2]))/2)
	bri.append((np.array(top1[-1])+np.array(top1[-3]))/2)
	bri.append((np.array(top1[-3])+np.array(top1[-2]))/2)
	trid.append((np.array(top1[-1])+np.array(top1[-2])+np.array(top1[-3]))/3)
	# print(getNormalVertax(top2))
	facet.append(getNormalVertax(top2))

	# print(simple)
for i in tri.convex_hulls:
	print(i)
# get cluster center get the average
center = getTopCenter(top1)

# print(center)
# get top adsorption sites and minus the center

t_normal = getTopNormal(top1,center)

for i in range(0,len(x1)):
	vertice.append([x1[i],y1[i],z1[i]])	

vertice = [list(zip(x1,y1,z1))]
# s_top = set(top1)
# print(len(s_top))
s_top = func4(top1)
print(len(s_top))
# ax.add_collection3d(Poly3DCollection(vertice))
# ax.view_init(20,20)

for i,k in zip(top1,t_normal):
	ax.scatter(i[0],i[1],i[2],c='r')
	ax.quiver(i[0],i[1],i[2],k[0],k[1],k[2],length=0.1)
# for j in bri:
# 	ax.scatter(j[0],j[1],j[2],c='g')

# for m,n in zip(trid,facet):
# 	ax.scatter(m[0],m[1],m[2],c='b')
	# ax.quiver(m[0],m[1],m[2],n[0],n[1],n[2],length=0.1)
# for i in range(0,len(top)-2):
# 	ax.scatter(top[i],top[i+1],top[i+2])
	# ax.scatter(topn[simple[0]-1][0],topn[simple[1]-1][1],topn[simple[2]-1][2])	


plt.show()

# view(atoms)

def v_denaunary (points):

	pass


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
def test():
	pass