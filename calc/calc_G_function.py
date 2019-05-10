# coding:utf-8
import numpy as np
import re
import pandas as pd
import math


in_file_dir='E:/cluster/'
# in_file='./data source/AuCu(single).txt'
out_file_name='res_ASCF_GCN_cjl_muti.csv'
input_filename = ['AgRh','PdAg','PdAu','PdCu','PdRh','PtAg','PtAu','PdRu','PtCu','PtPd','PtRh','PtRu']
# out_file_name = ['AgRh.csv',]
# atom_total_num=30
# atom_total_num=30	
distanz_err=0.1
Rc=8.0
yita=5.0
Rs=4.0
k=0.5
namuda=-1
yipselong=1
cn_max=12


def my_cos(x):
	return 1-0.5*(x**2)+(1/24.0)*(x**4)
	pass

def cut_function(x):
	out=0
	if x<=Rc:
		out=0.5*(math.cos(math.pi*x/Rc)+1)
		pass
	return out
	pass

def make_Rij(point_array,atom_total_num=0):
	point_array=point_array.astype(np.float64)

	ones=np.matrix(np.ones(atom_total_num)).T
	H_matrix=ones*np.matrix(point_array[atom_total_num-1])
	point_array_matrix=np.matrix(point_array)
	H_delta_matrix=point_array_matrix-H_matrix
	Rij=np.array(H_delta_matrix*(H_delta_matrix.T)).diagonal()
	Rij=np.sqrt(Rij)

	return Rij

	pass

def make_CUT_Rij(Rij,atom_total_num=0):
	r=map(cut_function,Rij)
	CUT_Rij=np.array(list(r))
	CUT_Rij[atom_total_num-1]=0

	return CUT_Rij

	pass

def make_G1(point_array,atom_total_num=0):
	Rij=make_Rij(point_array,atom_total_num)
	#print(Rij)
	CUT_Rij=make_CUT_Rij(Rij,atom_total_num)
	#print(CUT_Rij)
	return CUT_Rij.sum()
	pass

def make_G2(point_array,atom_total_num=0):
	Rij=make_Rij(point_array,atom_total_num)
	CUT_Rij=make_CUT_Rij(Rij,atom_total_num)

	f=lambda x:math.exp(-1*yita*(x-Rs))
	r=map(f,Rij)
	exp_yita_Rij=np.array(list(r))
	#print(Rij)
	#print(exp_yita_Rij)
	#print(CUT_Rij)
	add_all=np.array(np.matrix(exp_yita_Rij)*np.matrix(CUT_Rij).T)
	#print(add_all)
	return add_all.sum()

	pass

def make_G3(point_array,atom_total_num=0):
	Rij=make_Rij(point_array,atom_total_num)
	CUT_Rij=make_CUT_Rij(Rij,atom_total_num)

	f=lambda x:math.cos(k*x)
	r=map(f,Rij)
	exp_yita_Rij=np.array(list(r))
	add_all=np.array(np.matrix(exp_yita_Rij)*np.matrix(CUT_Rij).T)
	return add_all.sum()

	pass

def make_cos_ijk(point_array,atom_total_num=1):
	point_array=point_array.astype(np.float64)

	ones=np.matrix(np.ones(atom_total_num)).T
	H_matrix=ones*np.matrix(point_array[atom_total_num-1])
	point_array_matrix=np.matrix(point_array)
	H_delta_matrix=point_array_matrix-H_matrix
	dot_ij=np.array(H_delta_matrix*(H_delta_matrix.T))

	Rij=make_Rij(point_array,atom_total_num)
	Rij[atom_total_num-1]=17362
	#print(Rij)

	Rij_Rik=np.array(np.matrix(Rij).T*np.matrix(Rij))
	#print(Rij_Rik)

	cos_ijk=dot_ij/Rij_Rik

	return cos_ijk


	pass

# def f_exp(x):
# 	cut=0.0
# 	if x<=Rc:
# 		cut=0.5*(math.cos(math.pi*x/Rc)+1)
# 		pass
# 	return (math.e**(-1*yita*x*x))*cut
# 	pass

def make_G4(point_array,distanz_narray,atom_total_num=0):
	Rij=make_Rij(point_array,atom_total_num)
	CUT_Rij=make_CUT_Rij(Rij,atom_total_num)

	#print(Rij)

	f_exp=lambda x:cut_function(x)*(math.e**(-1*yita*x*x))
	#f_exp=lambda x:cut_function(x)

	r=map(f_exp,Rij)
	F_Rij=np.array(list(r))
	F_Rij[atom_total_num-1]=0
	#print(F_Rij)

	F_Rij_Rik=np.array(np.matrix(F_Rij).T*np.matrix(F_Rij))
	#print(F_Rij_Rik)

	distanz_narray[distanz_narray>Rc]=Rc
	#1-0.5*((math.pi*x/Rc)**2)+(1/24.0)*((math.pi*x/Rc)**4)
	f_exp_jk=lambda x:(math.e**(-1*yita*x*x))*0.5*((1-0.5*((math.pi*x/Rc)**2)+(1/24.0)*((math.pi*x/Rc)**4))+1)
	r=map(f_exp_jk,distanz_narray)
	F_Rjk=np.array(list(r))
	#print(F_Rjk)

	F_Rij_Rik_Rjk=F_Rij_Rik*F_Rjk
	#print(F_Rij_Rik_Rjk)

	cos_ijk=make_cos_ijk(point_array)
	#print(cos_ijk)

	F_exp_namuda=(1+namuda*cos_ijk)**yipselong

	G4_matrix=(2**(1-yipselong))*np.array(np.matrix(F_exp_namuda)*np.matrix(F_Rij_Rik_Rjk).T)
	#print(G4_matrix)

	#print(G4_matrix.sum())

	return G4_matrix.sum()

	pass

def make_G5(point_array,atom_total_num=0):
	Rij=make_Rij(point_array,atom_total_num)
	CUT_Rij=make_CUT_Rij(Rij,atom_total_num)

	#print(Rij)

	f_exp=lambda x:cut_function(x)*(math.e**(-1*yita*x*x))
	#f_exp=lambda x:cut_function(x)

	r=map(f_exp,Rij)
	F_Rij=np.array(list(r))
	F_Rij[atom_total_num-1]=0
	#print(F_Rij)

	F_Rij_Rik=np.array(np.matrix(F_Rij).T*np.matrix(F_Rij))
	#print(F_Rij_Rik)

	cos_ijk=make_cos_ijk(point_array)
	#print(cos_ijk)

	F_exp_namuda=(1+namuda*cos_ijk)**yipselong

	G5_matrix=(2**(1-yipselong))*np.array(np.matrix(F_exp_namuda)*np.matrix(F_Rij_Rik).T)
	#print(G4_matrix)

	#print(G4_matrix.sum())

	return G5_matrix.sum()

	pass

def make_angle_array(min_index,coor_array_min,point_array,atom_total_num=0):
	point_array=point_array.astype(np.float64)
	#angle_array=[-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2]
	angle_array= np.linspace(-2,-2,atom_total_num).astype(np.int16)
	#print(point_array[13])
	min_H=point_array[atom_total_num-1]-point_array[min_index]
	L_min_H=np.sqrt(min_H.dot(min_H))
	i=0
	while i<atom_total_num:
		if coor_array_min[i]!=0:
			min_dort=point_array[i]-point_array[min_index]
			L_min_dort=np.sqrt(min_dort.dot(min_dort))
			angle_array[i]=min_dort.dot(min_H)/(L_min_dort*L_min_H)
			pass
		i+=1
		pass
	return angle_array
	pass

def make_norm_matrix(point_array,atom_total_num=0):

	#ein_array=np.array([1,1,1,1,1,1,1,1,1,1,1,1,1,1]).astype(np.float64)
	ein_array=np.linspace(1,1,atom_total_num).astype(np.float64)
	minus_ein_array=np.linspace(-1,-1,atom_total_num).astype(np.float64)
	x_array=np.matrix(point_array.T[0]).astype(np.float64)
	y_array=np.matrix(point_array.T[1]).astype(np.float64)
	z_array=np.matrix(point_array.T[2]).astype(np.float64)

	#print(ein_array+ein_array)

	x_2d=np.array(((x_array.T)*np.matrix(ein_array))+((np.matrix(minus_ein_array).T)*x_array))
	y_2d=np.array(((y_array.T)*np.matrix(ein_array))+((np.matrix(minus_ein_array).T)*y_array))
	z_2d=np.array(((z_array.T)*np.matrix(ein_array))+((np.matrix(minus_ein_array).T)*z_array))

	distanz=np.sqrt(x_2d*x_2d+y_2d*y_2d+z_2d*z_2d)

	return distanz
	pass

def find_nearst_atom_index(distanz_narray,atom_total_num=0):
	atom_array=distanz_narray[atom_total_num-1]
	atom_arr=np.vstack((atom_array,np.arange( 0, atom_total_num, 1 ).astype(np.int16))).T
	atom_pd=pd.DataFrame(atom_arr,columns=['distanz','index'])
	atom_pd=atom_pd.sort_values(by='distanz')
	return int(np.array(atom_pd)[1][1])
	pass

def make_coor_matrix(distanz_narray,atom_total_num=0):
	coor_matrix=[]
	i=0
	while i<atom_total_num:
		#print("******************"+str(i)+"**********************")
		add_array=np.linspace(0,0,atom_total_num).astype(np.int16)
		atom_array=distanz_narray[i]
		atom_arr=np.vstack((atom_array,np.arange( 0, atom_total_num, 1 ).astype(np.int16))).T
		atom_pd=pd.DataFrame(atom_arr,columns=['distanz','index'])
		#print(atom_pd)
		atom_pd=atom_pd.sort_values(by='distanz')
		#print(atom_pd)
		up_distanz=np.array(atom_pd)[1][0]
		#判断最近的是不是H
		if int(np.array(atom_pd)[1][1])==atom_total_num-1:
			up_distanz=np.array(atom_pd)[2][0]
			pass
		up_distanz=up_distanz*(1+distanz_err)
		#print(up_distanz)
		#print('index')
		for x in atom_pd.index:
			at_index=int(atom_pd['index'][x])
			if atom_pd['distanz'][x]<=up_distanz:
				add_array[at_index]=1
				#print(at_index)
			else:
				add_array[at_index]=0
				pass
			pass
		add_array[atom_total_num-1]=0
		add_array[i]=0
		coor_matrix.append(add_array)
		i+=1
		pass
	coor_matrix=np.array(coor_matrix)
	return coor_matrix
	pass

def make_GCN(coor_matrix):
	CN=np.sum(coor_matrix,axis=1)
	CN_matrix=np.matrix(CN)
	coor_matrix_np=np.matrix(coor_matrix)
	GCN=(np.array(CN_matrix*coor_matrix_np.T)/cn_max)[0]
	#print(GCN)
	return GCN
	pass

def make_GCN_mean_ft(GCN,coor_matrix,atom_total_num=0):
	H_matrix=coor_matrix[atom_total_num-1]
	H_matrix=H_matrix/(H_matrix.sum())
	GCN_mean=np.array(np.matrix(H_matrix)*np.matrix(GCN).T).sum()
	#print(GCN_mean)
	return GCN_mean
	pass

def make_cjl_GCN(GCN,coor_matrix,atom_total_num=0):
	H_matrix=coor_matrix[atom_total_num-1]
	# print(H_matrix)
	# print(coor_matrix)
	cjl_array=np.array(np.matrix(H_matrix)*np.matrix(coor_matrix))
	cjl_array[cjl_array>0]=1
	cjl_array=cjl_array/(cjl_array.sum())
	GCN_cjl=np.array(np.matrix(cjl_array)*np.matrix(GCN).T).sum()
	# print(GCN_cjl)
	return GCN_cjl
	pass

def make_result(point_array):
	pass

def make_feature(in_file,out_file):
	data_file=open(in_file)
	str_line=data_file.readline()

	# my_res=[[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]]
	# my_columns=['dE','min','min_CE0','min_CE1','min_CE2','min_CE3','min_CE4','min_CE5','min_CE6','min_CE7','min_CE8','min_CE9','min_CE10','min_CE11','min_CE12','min_CE13','CE0','CE1','CE2','CE3','CE4','CE5','CE6','CE7','CE8','CE9','CE10','CE11','CE12','CE13','ang0','ang1','ang2','ang3','ang4','ang5','ang6','ang7','ang8','ang9','ang10','ang11','ang12','ang13']
	# my_res=pd.DataFrame(my_res,columns=my_columns)

	# my_res=[[0,0,0,0,0,0,0]]
	# my_columns=['G1','G2','G3','G4','G5','GCN_mean','GCN_cjl']
	# my_res=pd.DataFrame(my_res,columns=my_columns)


	atom_cnt=0

	while len(str_line)!=0:
		# dE=re.findall(r'-?\d+\.?\d*e?-?\d*?', data_file.readline())
		atom_total_num=0
		# dE=data_file.readline().split('=')[1]

		# print(dE)
		# dE=float(dE)
		j = 0
		for j in range(0,6):
			str_line = data_file.readline()
			if j == 4:
				str_line = data_file.readline().strip().split()
				atom_total_num = int(str_line[0])+int(str_line[1])
			pass
		# print(atom_total_num)
		point_array=[]
		i=0
		while i<atom_total_num:
			str_line=data_file.readline()
			str_split=str_line.strip().split()
			# print(str_split)
			point=np.array([str_split[0],str_split[1],str_split[2]])
			point_array.append(point)
			i+=1
			pass
		point_array=np.array(point_array)
		distanz_narray=make_norm_matrix(point_array,atom_total_num)
		coor_matrix=make_coor_matrix(distanz_narray,atom_total_num)
		min_index=find_nearst_atom_index(distanz_narray,atom_total_num)
		angle_array=make_angle_array(min_index,coor_matrix[min_index],point_array,atom_total_num)
		#print(coor_matrix)
		GCN=make_GCN(coor_matrix)
		GCN_mean=make_GCN_mean_ft(GCN,coor_matrix,atom_total_num)
		GCN_cjl=make_cjl_GCN(GCN,coor_matrix,atom_total_num)
		# print(make_G1(point_array))
		# print(make_G2(point_array))
		# print(make_G3(point_array))
		# print(make_G4(point_array,distanz_narray))
		# print(make_G5(point_array))
		G1=make_G1(point_array,atom_total_num)
		G2=make_G2(point_array,atom_total_num)
		G3=make_G3(point_array,atom_total_num)
		G4=make_G4(point_array,distanz_narray,atom_total_num)
		G5=make_G5(point_array,atom_total_num)
		result=[G1,G2,G3,G4,G5,GCN_mean,GCN_cjl]
		# result=np.hstack((result,coor_matrix[min_index]))
		# result=np.hstack((result,np.sum(coor_matrix,axis=1)))
		# result=np.hstack((result,angle_array))
		res=[]
		res.append(result)
		res=np.array(res)
		# res=pd.DataFrame(res,columns=my_columns)
		#print(res)
		# my_res=my_res.append(res,ignore_index=True)
		#print(my_res)
		atom_cnt+=1
		print("data count:"+str(atom_cnt))
		str_line=data_file.readline()
		# str_line=''
		pass

	# my_res.to_csv(out_file)

	data_file.close()
	return result
	pass



my_res=[[0,0,0,0,0,0,0]]
my_columns=['G1','G2','G3','G4','G5','GCN_mean','GCN_cjl']
my_res=pd.DataFrame(my_res,columns=my_columns)
tmp_res = []
for i in input_filename:

	filename = in_file_dir+i+'.vasp'
	# tmp_res = np.array(make_feature(filename,out_file_name))
	tmp_res.append(make_feature(filename,out_file_name))
	print(make_feature(filename,out_file_name))
	print(tmp_res)
my_res=np.array(tmp_res)
my_res=pd.DataFrame(my_res,columns=my_columns)
my_res.to_csv(out_file_name)


######################################################################################################################################

def makeHAboPoint(point_array):


	pass





#一组数据
# dE=re.findall(r'-?\d+\.?\d*e?-?\d*?', data_file.readline())
# print(dE)
# point_array=[]
# i=0
# while i<14:
# 	str_line=data_file.readline()
# 	str_split=str_line.split()
# 	point=np.array([str_split[1],str_split[2],str_split[3]])
# 	point_array.append(point)
# 	i+=1
# 	pass





# print(np.array(point_array))

# point_array=np.array(point_array)

# distanz_narray=make_norm_matrix(point_array)

# print(distanz_narray)

# coor_matrix=make_coor_matrix(distanz_narray)

# print(coor_matrix)

# print(np.sum(coor_matrix,axis=1))

# np.savetxt("coor_matrix.txt",coor_matrix)

# min_index=find_nearst_atom_index(distanz_narray)

# print(min_index)

# angle_array=make_angle_array(min_index,coor_matrix[min_index],point_array)

# print(angle_array)


