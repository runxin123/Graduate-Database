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

# from catkit.gen.surface import SlabGenerator
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
# from catkit.gen.adsorption import AdsorptionSites
from scipy.spatial import ConvexHull,Delaunay
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 
# from matplotlib.collections import PolyCollection
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import copy
from modeltest import *
import pandas as pd
surfaces = [(1, 0, 0), (1, 1, 0), (1, 1, 1)]
layers = [4, 4, 4]
lc = 3.61000
atoms = FaceCenteredCubic('Pt', surfaces, layers, latticeconstant=lc)
p = atoms.get_positions()
print(p)
facet = []
def minusList(list1,list2):
	lisindex = min(len(list1),len(list2))
	c = [list1[i] - list2[i] for i in range(lisindex)]
	return c
# def func4(one_list):
#     '''
#     使用排序的方法
#     '''
#     result_list=[]
#     temp_list=sorted(one_list)
#     i=0
#     while i<len(temp_list):
#         if temp_list[i] not in result_list:
#             result_list.append(temp_list[i])
#         else:
#             i+=1
#     return result_list
def getNormalVertax(face_array):
	'''
	face_array [[x1,y1,z1],[x2,y2,z2],[x3,y3,z3]]
	'''
	face = []
	x_1 = face_array[0][0]-face_array[1][0]
	y_1 = face_array[0][1]-face_array[1][1]
	z_1 = face_array[0][2]-face_array[1][2]
	x_2 = face_array[0][0]-face_array[2][0]
	y_2 = face_array[0][1]-face_array[2][1]
	z_2 = face_array[0][2]-face_array[2][2]
	vector_1 = [x_1,y_1,z_1]
	vector_2 = [x_2,y_2,z_2]
	face = np.cross(vector_1,vector_2).tolist()
	# print(face)
	# face.append((face_array[0][1]-face_array[1][1])*(face_array[1][2]-face_array[2][2])-(face_array[0][2]-face_array[2][2])*(face_array[1][1]-face_array[2][1]))
	# face.append((face_array[1][0]-face_array[2][0])*(face_array[0][2]-face_array[2][2])-(face_array[0][0]-face_array[1][0])*(face_array[1][2]-face_array[2][2]))
	# face.append((face_array[0][0]-face_array[1][0])*(face_array[1][1]-face_array[2][1])-(face_array[0][1]-face_array[1][1])*(face_array[1][0]-face_array[2][0]))
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


def writePoints(ligand,points,filedir,filename,num):
	with open(filedir+'/'+filename,'w') as f:
		f.write(str(num)+'\n')
		for i in points:
			f.write(str(i[0])+' '+str(i[1])+' '+str(i[2])+'\n')
		for j in ligand:
			f.write(str(j[0])+' '+str(j[1])+' '+str(j[2])+'\n')
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


def add_abs(vertex,d,point_array):
	abs_array = []
	if d > 0 :
		for i,j in zip(vertex,point_array):
			abs_array_t = []
			abs_array_t.append(j[0]+i[0]*d)
			abs_array_t.append(j[1]+i[1]*d)
			abs_array_t.append(j[2]+i[2]*d)
			abs_array.append(abs_array_t)
		pass
	return abs_array

def getBriCenter(point1,point2):
	Bri = []
	for i,j in zip(point1,point2):
		Bri.append((i+j)/2)
	return Bri

def getAngleBisector(vector1,vector2):
	vector1 = np.array(vector1)
	vector2 = np.array(vector2)
	v_t1 = np.linalg.norm(vector1, axis=1, keepdims=True)
	v_t2 = np.linalg.norm(vector2, axis=1, keepdims=True)
	vector1 = vector1/v_t1
	vector2 = vector2/v_t2
	c = (vector1+vector2).tolist
	return c

def getFacetCenter(vector1,vector2,vector3):
	v1 = np.array(vector1)
	v2 = np.array(vector2)
	v3 = np.array(vector3)
	result = (v1+v2+v3)/3
	result = result.tolist()
	return result


def equalVector(vector1,vector2):
	if (vector1[0]==vector2[0]) and (vector1[1]==vector2[1]) and (vector1[2]==vector2[2]):
		return True
	else:
		return False
def not_equal_v(vector1,vector2):
	'''
	[[[x,y,z],[x1,y1,z1]],[]]
	'''
	v1 = np.array(vector1)
	v2 = np.array(vector2)
	if (v1 == v2).all():
		return False
	else:
		return True
		pass


def getBriFacet(facet,cluster_center):
	result = []
	bri = []
	biosector = []
	bri_sector = ()
	a_sum = 0
	facet_t = copy.deepcopy(facet)
	for i in facet:
		for j in facet_t:
			# print('i = '+str(i)+'\n'+'j = '+str(j)+'\n')
			if not_equal_v(i,j):
				a = []
				a_sum = 0
				for m in i:
					for n in j:
						# print('m = '+str(m))
						if equalVector(m,n) == True:
							a.append(m)
				a_sum = len(a)
				# print('a = '+str(a))
				if a_sum == 2:
					# print(a)
					bri = getBriCenter(a[0],a[1])
					face_normal_t1 = getNormalVertax(i)
					face_normal_t2 = getNormalVertax(j)
					face_center_1 = getFacetCenter(i[0],i[1],i[2])
					face_center_2 = getFacetCenter(j[0],j[1],j[2])
					face_normal1 = tri_normal_single(cluster_center,face_center_1,face_normal_t1)
					face_normal2 = tri_normal_single(cluster_center,face_center_2,face_normal_t2)
					biosector = getAngleBisector(face_normal1,face_normal2)
					bri_sector = (bri,biosector)
					result.append(bri_sector)
	# print('result = '+str(result))
	return result


def bri_test(arg,bri,pointNumber,top_array):
	print(bri[0])
	# my_res=[[0,0,0,0,0,0,0]]
	my_res=[[0,0,0,0,0]]
	# my_columns=['G1','G2','G3','G4','G5','GCN_mean','GCN_cjl']
	my_columns = ['G1','G2','G3','G4','G5']
	my_res=pd.DataFrame(my_res,columns=my_columns)
	tmp_array = []
	for i in bri:
		st = copy.deepcopy(top_array)
		st.append(i)
		tmp_array.append(calc_G_function(st,pointNumber))
	my_res = np.array(func4(tmp_array))
	# print(len(my_res))
	my_res = pd.DataFrame(my_res,columns=my_columns)
	my_res.to_csv(arg+'.csv')



def tri_test(arg,tri,pointNumber,top_array):
	# my_res=[[0,0,0,0,0,0,0]]
	my_res = [[0,0,0,0,0]]
	# my_columns=['G1','G2','G3','G4','G5','GCN_mean','GCN_cjl']
	my_columns = ['G1','G2','G3','G4','G5']
	my_res=pd.DataFrame(my_res,columns=my_columns)
	tmp_array = []
	for i in tri:
		st = copy.deepcopy(top_array)
		st.append(i)
		tmp_array.append(calc_G_function(st,pointNumber))
	my_res = np.array(func4(tmp_array))
	# print(len(my_res))
	my_res = pd.DataFrame(my_res,columns=my_columns)
	my_res.to_csv(arg+'.csv')


def top_test(arg,top,pointNumber,top_array):
	tmp_array = []
	# my_res=[[0,0,0,0,0,0,0]]
	my_res=[[0,0,0]]
	# my_columns=['G1','G2','G3','G4','G5','GCN_mean','GCN_cjl']
	my_columns=['G1','G2','G3']
	my_res=pd.DataFrame(my_res,columns=my_columns)
	for i in top:
		st = copy.deepcopy(top_array)
		st.append(i)
		tmp_array.append(calc_G_function(st,pointNumber))
	my_res = np.array(func4(tmp_array))
	# print(len(my_res))
	my_res = pd.DataFrame(my_res,columns=my_columns)
	my_res.to_csv(arg+'.csv')

# view(atoms)
# coordinates, connectivity = gen.adsorption_sites(atoms)
# sites = AdsorptionSites(atoms)
# sites.plot('./adsorption-sites.png')
# for i in p:

print(len(p))
hull = ConvexHull(p,False,'Qc')

# using coplaner 
# print('hull.points'+str(len(hull.points)))

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


points = np.zeros((0, 3))
points[:, 0] = np.array(x)
points[:, 1] = np.array(y)
points[:, 2] = np.array(z)
topn = np.array(top)
tri = Delaunay(p,qhull_options='Qc')
# print(len(p))
x1 = []
y1 = []
z1 = []
vertice = []
top2 = []
tri_result = []
# s_top = set()
# for x in tri.vertex_neighbor_vertices:
	# for i in x:
	# print(x[1])


# print(tri.vertex_to_simplex)
	# print('h '+str(i))
	# pass
for simple in tri.convex_hull:
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
	tri_result.append(top2)
	facet.append(getNormalVertax(top2))

	# print(simple)
# for i in tri.convex_hull:
# 	print(i)
# get cluster center get the average

print('top2'+str(top2)+'\n')

for i in range(0,len(x1)):
	vertice.append([x1[i],y1[i],z1[i]])	

vertice = [list(zip(x1,y1,z1))]
# s_top = set(top1)
# print(len(s_top))
# global s_top
s_top = func4(top1)
print('stop = '+str(len(s_top)))

pointNumber = len(s_top)
center = getTopCenter(s_top)

# print(center)
# get top adsorption sites and minus the center

t_normal = getTopNormal(s_top,center)
# ax.add_collection3d(Poly3DCollection(vertice))
# ax.view_init(20,20)

for i,k in zip(s_top,t_normal):
	ax.scatter(i[0],i[1],i[2],c='r')
	ax.quiver(i[0],i[1],i[2],k[0],k[1],k[2],length=0.1)


# for j in bri:
# 	ax.scatter(j[0],j[1],j[2],c='g')





pointNumber = pointNumber+1




# tri_test('trid',trid,pointNumber,s_top)

arr = []
# arr = []



def tri_normal(trid_center,trid):
	modify_trid = []
	mis = np.array([-1,-1,-1])
	for i,j in zip(trid_center,trid):
		i_dot = np.array(i)
		j_dot = np.array(j)
		result = np.dot(i_dot,j_dot)
		if result < 0:
			j_dot = j_dot*mis
			modify_trid.append(j_dot.tolist())
		else:
			modify_trid.append(j_dot.tolist())
	return modify_trid
	pass

def tri_normal_single(center,face_center,trid):
	modify_trid = []
	mis = np.array([-1,-1,-1])
	i_dot = np.array(face_center)-np.array(center)
	j_dot = np.array(trid)
	result = np.dot(i_dot,j_dot)
	if result > 0:
		j_dot = j_dot*mis
		modify_trid = j_dot.tolist()
	else:
		modify_trid = j_dot.tolist()

	return modify_trid
	pass



# top_test('hh',arr,pointNumber,s_top)

trid_center = getTopNormal(trid,center)

print(trid_center)
facet_normal = tri_normal(trid_center,facet)

# print(facet)

# print('tri = '+str(tri_result))

# getBriFacet(tri_result,center)



for m,n in zip(trid,trid_center):
	# ax.scatter(m[0],m[1],m[2],c='b')
	ax.quiver(m[0],m[1],m[2],n[0],n[1],n[2],length=0.5)
	
arr = add_abs(trid_center,0.4,trid)

for i,j in zip(arr,trid):
	ax.scatter(j[0],j[1],j[2],c='b')
	ax.scatter(i[0],i[1],i[2],c='g')


plt.show()

# view(atoms)

def v_denaunary (points):

	pass



