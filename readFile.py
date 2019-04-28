# from ase import *
# import sys
# import json
# import os
import numpy as np
import re
from scipy.integrate import simps
#
import plotly.plotly as pl
import plotly.graph_objs as go
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
class POSCAR(object):
    first_line = ''
    element_types = []
    element = {}  # {'Pt':1,'Pd':1}
    element_num = []
    element_nums = 0  # count elements numbers
    lattice_constant = {'x': 1, 'y': 1, 'z': 1, 'xy': 0, 'xz': 0, 'yz': 0}
    lattice_plane = []
    matrix_context = 1.0
    cart_or_direct = 1  # 1 was cart 0 was direct
    selective = 0  # 1 was selective 0 was no-selective
    coord = {}  # {'1':['Pt',1,1,1,'T','T','T']} # LINES NUMS ,ELEMENTS TYPES, COORDINATIONS SELECTIVE FIX OR NOT
    filedir = ''
    filename = ''
    content = ''
    suffix_file = '.txt'

    def read(self, filedir, filename):
        # filedir the poscar file dir
        # filename the poscar contcar file name
        try:
            f = open(filedir + '/' + filename, 'r')
            i = 0
            for line in f:
                line_t = line.strip()
                print(line)
                i = i + 1
                if i == 1:
                    self.first_line = line_t
                elif i == 2:
                    self.matrix_context = float(line_t)
                elif i == 3:
                    self.lattice_constant['x'] = float(
                        line_t.split()[0])  # using the default split to split the many space
                    self.lattice_constant['xy'] = float(line_t.split()[1])
                    self.lattice_constant['xz'] = float(line_t.split()[2])
                elif i == 4:
                    self.lattice_constant['y'] = float(line_t.split()[1])
                    self.lattice_constant['yz'] = float(line_t.split()[2])
                elif i == 5:
                    self.lattice_constant['z'] = float(line_t.split()[2])
                elif i == 6:
                    eletypes = line_t.split()
                    for x in eletypes:
                        self.element_types.append(x)  # get the elements_types
                elif i == 7:
                    elenum = line_t.split()
                    num = 0
                    for t in elenum:
                        num = int(t)
                        self.element_num.append(num)  # get the elemt nums
                        self.element_nums += num
                    for m in range(len(self.element_types)):
                        self.element[self.element_types[m]] = self.element_num[m]  # get element json
                elif i == 8:
                    if line_t[0] == 's' or line_t[0] == 'S':
                        # print(line_t[1])
                        self.selective = 1
                    else:
                        # self.selective = 0
                        if line_t[0] == 'C' or line_t[0] == 'c':
                            self.cart_or_direct = 1

                        else:
                            self.cart_or_direct = 0
                        self.coord[i - 8] = line_t
                elif i == 9:
                    # self.selective = 0
                    if line_t[0] == 'C' or line_t[0] == 'c':
                        self.cart_or_direct = 1
                        self.coord[i - 9] = line_t
                    elif line_t[0] == 'D' or line_t[0] == 'd':
                        self.cart_or_direct = 0
                        self.coord[i - 9] = line_t
                    else:
                        self.coord[i - 8] = line_t.split()
                else:
                    self.coord[i - 9] = line_t.split()
        except IOError:
            print("Error")
        else:
            print("content read file success")
            f.close()

    def read_content(self,content):
        # with open('E:/CONTCAR','r') as f:
        #     self.content = f.read()
            # self.content.split('\n')
        print(self.content.split('\n'))
    def write(self, filedir, filename):
        try:
            f = open(filedir + '/' + filename, 'w')
            i = 0
            space = '  '
            f.write(self.first_line + '\n')

            f.write(space + str(self.matrix_context) + '\n')
            f.write(space + str(self.lattice_constant['x']) + space + str(self.lattice_constant['xy']) + space + str(
                self.lattice_constant['xz']) + '\n')
            f.write(space + str(self.lattice_constant['xy']) + space + str(self.lattice_constant['y']) + space + str(
                self.lattice_constant['yz']) + '\n')
            f.write(space + str(self.lattice_constant['xz']) + space + str(self.lattice_constant['yz']) + space + str(
                self.lattice_constant['z']) + '\n')


            for i in self.element_types:
                f.write(i + space)
            f.write('\n')
            for t in self.element_num:
                f.write(str(t) + space)
            f.write('\n')


            if self.selective == 1:
                f.write("Selective\n")
            else:
                pass
            if self.cart_or_direct == 1:
                f.write('Cartisian\n')
            else:
                f.write('Direct\n')

            # write the coordination
            if self.selective == 1:
                for i in range(1, self.element_nums+1):
                    f.write(space + str(self.coord[i][0]) + space + str(self.coord[i][1]) + space + str(
                        self.coord[i][2]) + space + str(self.coord[i][3]) + space + str(self.coord[i][4]) + space + str(
                        self.coord[i][5]) + '\n')
                    print(space + str(self.coord[i][0]) + space + str(self.coord[i][1]) + space + str(
                        self.coord[i][2]) + space + str(self.coord[i][3]) + space + str(self.coord[i][4]) + space + str(
                        self.coord[i][5]) + '\n')  # test for the seq
            else:
                for i in range(1, self.element_nums+1):
                    f.write(space + str(self.coord[i][0]) + space + str(self.coord[i][1]) + space + str(
                        self.coord[i][2]) + '\n')
                    print(space + str(self.coord[i][0]) + space + str(self.coord[i][1]) + space + str(
                        self.coord[i][2]) + '\n')
        except IOError:
            print("Error")
        else:
            print("content read file success")
            f.close()

    def dir2cart(self):
        try:
            if self.cart_or_direct == 0:
                for i in range(1, self.element_nums):
                    x_t = float(self.coord[i][0])
                    # print(x_t)
                    y_t = float(self.coord[i][1])
                    # print(y_t)
                    z_t = float(self.coord[i][2])
                    # print(z_t)
                    self.coord[i][0] = x_t * self.lattice_constant['x'] + y_t * self.lattice_constant['xy'] + z_t * \
                                                                                                              self.lattice_constant[
                                                                                                                  'xz']
                    self.coord[i][1] = x_t * self.lattice_constant['xy'] + y_t * self.lattice_constant['y'] + z_t * \
                                                                                                              self.lattice_constant[
                                                                                                                  'yz']
                    self.coord[i][2] = x_t * self.lattice_constant['xz'] + y_t * self.lattice_constant['yz'] + z_t * \
                                                                                                               self.lattice_constant[
                                                                                                                   'z']
                self.cart_or_direct = 1 # convert the dir to cart successfully
            else:
                print('was cart not to dir')
        except IOError:
            print('Error')
        else:
            print("convert successfully")

    def cart2dir(self):
        try:
            if self.cart_or_direct == 1:
                a = np.array([[self.lattice_constant['x'], self.lattice_constant['xy'], self.lattice_constant['xz']],
                              [self.lattice_constant['xy'], self.lattice_constant['y'], self.lattice_constant['yz']],
                              [self.lattice_constant['xz'], self.lattice_constant['yz'], self.lattice_constant['z']]])
                print(np.linalg.inv(a))
                verse_a = np.linalg.inv(a)  # verse_a = a^-1
                for i in range(1, self.element_nums):
                    x_t = float(self.coord[i][0])
                    y_t = float(self.coord[i][1])
                    z_t = float(self.coord[i][2])
                    self.coord[i][0] = x_t * verse_a[0][0] + y_t * verse_a[0][1] + z_t * verse_a[0][2]
                    self.coord[i][1] = x_t * verse_a[1][0] + y_t * verse_a[1][1] + z_t * verse_a[1][2]
                    self.coord[i][2] = x_t * verse_a[2][0] + y_t * verse_a[2][1] + z_t * verse_a[2][
                        2]  # get reverse of the matrix
                self.cart_or_direct = 0
            else:
                print('was dir not to cart')

        except IOError:
            print('was dir not to cart')
        else:
            print('convert successfully')
    # fix slab layer param was layer = 2 calc the functions for the poscar
    # change poscar to
    def fixlayer(self,layer):
        pass

class INCAR(object):
    content = ''
    incar_prop = {}
    IALGO = 48
    ENCUT = 400.00 # DEFAULT PROPERTY OF THE INCAR
    SYSTEM = 'unknown'
    PREC = 'Normal'
    NELM = 100
    # IONIC Relaxation
    NSW = 400
    ISIF = 0
    IBRION = 2
    POTIM = 0.1
    EDIFFG = -0.05
    # DOS RELATED VALUES
    ISMEAR = 1
    SIGMA = 0.2
    IDIPOL = 3
    # SGI ORIGIN OPTIMIZED VALUES:
    LWAVE = ''
    LCHARG = ''
    LVTOT = ''
    LREAL = ''
    NPAR = 2
    EDIFF = 0.1
    def read(self,filedir,filename):
        with open(filedir+'/'+filename,'r') as f:
            for line in f:
                line_t = line.strip()
                if line_t.find('=')>1:
                    array_t = line_t.strip().split('=')
                    # print(array_t[1].strip().split()[0])
                    self.incar_prop[array_t[0].strip().upper()] = array_t[1].strip().split()[0]
                else:
                    pass
            for i in self.incar_prop:
                if hasattr(self,i):
                    pass
                else:
                    setattr(self,i,self.incar_prop[i])
        self.IALGO = int(self.incar_prop['IALGO'])
        self.SYSTEM = self.incar_prop['SYSTEM']
        self.ENCUT = float(self.incar_prop['ENCUT'])
        self.PREC = self.incar_prop['PREC']
        self.EDIFFG = float(self.incar_prop['EDIFFG'])
        self.EDIFF = float(self.incar_prop['EDIFF'])
        self.NSW = int(self.incar_prop['NSW'])
        self.IBRION = int(self.incar_prop['IBRION'])
        self.POTIM = float(self.incar_prop['POTIM'])
        self.ISMEAR = int(self.incar_prop['ISMEAR'])
        self.SIGMA = float(self.incar_prop['SIGMA'])
        self.IDIPOL = int(self.incar_prop['IDIPOL'])
        self.LWAVE = self.incar_prop['LWAVE'].replace('.','')
        self.LVTOT = self.incar_prop['LVTOT'].replace('.','')
        self.LREAL = self.incar_prop['LREAL'].replace('.','')
        self.LCHARG = self.incar_prop['LCHARG'].replace('.','')
        # self.NELM = self.incar_prop['NELM']
        self.NPAR = int(self.incar_prop['NPAR'])

    # read from the database
    def read_content(self,content):
        self.content = content
        print(self.content.split('\n'))


    def write(self,filedir,filename):
        with open(filedir+'/'+filename,'w') as f:
            f.write('IALGO'+' = '+str(self.IALGO)+'\n')
            f.write('SYSTEM'+' = '+self.SYSTEM+'\n')
            f.write('ENCUT'+' = '+str(self.ENCUT)+'\n')
            f.write('PREC'+' = '+self.PREC+'\n')
            f.write('EDIFF'+' = '+str(self.EDIFF)+'\n')
            f.write('EDIFFG'+' = '+str(self.EDIFFG)+'\n')
            f.write('POTIM'+' = '+str(self.POTIM)+'\n')
            f.write('SGIMA'+' = '+str(self.SIGMA)+'\n')
            f.write('IDIOL'+' = '+str(self.IDIPOL)+'\n')
            f.write('LWAVE'+' = '+'.'+self.LWAVE+'.'+'\n')
            f.write('LREAL'+'='+'.'+self.LREAL+'.'+'\n')
            f.write('LCHARG'+' = '+'.'+self.LCHARG+'.'+'\n')
            f.write('LVTOT'+' = '+'.'+self.LVTOT+'.'+'\n')
            f.write('NPAR'+'='+str(self.NPAR)+'\n')
            f.write('NELM'+'='+str(self.NELM)+'\n')
            # f.write('SYSTEM'+'='+self.SYSTEM+'\n')
    # write from the incar_prop file
    def write2(self,filedir,filename):
        with open(filedir+'/'+filename,'w') as f:
            for i in self.incar_prop:
                f.write(i+'='+self.incar_prop[i]+'\n')

class KPOINTS(object):
    content = ''
    first_line = 'unknown' # KPOINTS LINE 1
    K_2 = 0 # KPOINTS LINE 2
    TYPE = 0 # TYPE 0 WAS GAMMA TYPE 1 WAS MACK
    TYPE_A = ['G','M'] # grid type GAMMA AND Mack
    GRID = [5,5,5] # grid some things
    TEMP = [0,0,0] # tmp
    def read(self,filedir,filename):
        with open(filedir+'/'+filename,'r') as f:
            self.first_line = f.readline().strip()
            self.K_2 = int(f.readline().strip())
            if f.readline().strip()[0].upper() == 'G':
                TYPE = 0
            elif 'M' == f.readline().strip()[0].upper():
                TYPE = 1
            else:
                pass
            grid_t = f.readline().strip().split()
            for i in range(0,len(grid_t)-1):
                self.GRID[i] = int(grid_t[i])
            temp_t = f.readline().strip().split()
            for t in range(0,len(temp_t)-1):
                self.TEMP[t] = int(temp_t[t])

    # read from the database
    def read_content(self,content):
        self.content = content
        print(content.split().strip())


    # write from the database
    def write(self,filedir,filename):
        space = ' '
        # WRITE THE KPOINTS
        with open(filedir+'/'+filename,'w') as f:
            f.write(self.first_line+'\n')
            f.write(str(self.K_2)+'\n')
            f.write(self.TYPE_A[self.TYPE]+'\n')
            f.write(str(self.GRID[0])+space+str(self.GRID[1])+space+str(self.GRID[2])+'\n')
            f.write(str(self.TEMP[0])+space+str(self.TEMP[1])+space+str(self.TEMP[2])+'\n')

    # change the Grid and TYPE of Grid
    def alterGrid(self,grid):
        self.GRID = grid
    def alterType(self,type):
        self.TYPE = type

class OSZICAR(object):
    content = ''
    N_step = [] #default N_step
    e0 = []# default N_step
    energy = []
    force = [] #
    eK =[]
    sp = []
    sk = []
    T = []
    energy_min = 0
    def read(self,filedir,filename):
        with open(filedir+'/'+filename,'r') as f:
            for line in f:
                line_t = line.strip()
                if re.match(r'\d',line_t,re.I):
                    # print(line_t.split())
                    line_trim = line_t.split()
                    self.N_step.append(int(line_trim[0]))
                    self.T.append(int(line_trim[2].replace('.','')))
                    self.energy.append(float(line_trim[4]))
                    self.e0.append(float(line_trim[8]))
                    self.force.append(float(line_trim[6]))
                    self.eK.append(float(line_trim[10]))
                    self.sp.append(float(line_trim[12]))
                    self.sk.append(float(line_trim[14]))
    # WRITE FROM THE DATABASE

    def read_content(self,content):
        self.content = content
        content_t = self.content.split()
        print(self.content)
    # get the Min energy from the slab or cluster
    def getMinEnergy(self):
        self.energy_min = self.energy[-1]
        return self.energy_min


class OUTCAR(object):
    e_fermi = []
    volume = 0
    is_check = 0 # default convergence from ionic was the is check was 1
    energy = []
    e0 = []
    energy_s = [] # free energy
    freq = [] # frequency of OUTCAR
    # read from the OUTCAR
    def read(self,filedir,filename):
        with open(filedir+'/'+filename,'r') as f:
            for line in f:
                line_t = line.strip()
                if re.match(r'E-fermi',line_t,re.I):
                    line_trimmed = line_t.split()
                    # line_ret = line_trimmed[2].split(':')
                    # print(line_trimmed)
                    self.e_fermi.append(float(line_trimmed[2])) # get E-fermi
                elif re.match(r'volume of cell',line_t,re.I):
                    line_trimmed = line_t.split(':')
                    # print(line_trimmed)
                    self.volume = float(line_trimmed[1].strip())
                elif re.match(r'reached required accuracy - stopping structural energy minimisation',line_t,re.I):
                    self.is_check = 1
                    print(line_t)
                elif re.match(r'free energy    TOTEN',line_t,re.I):
                    line_trimmed = line_t.split('=')
                    # print(line_trimmed)
                    self.energy.append(float(line_trimmed[1].strip().replace('eV','').strip()))
                elif re.match(r'energy  without entropy', line_t, re.I):
                    line_trimmed = line_t.split()
                    self.energy_s.append(float(line_trimmed[3]))
                    self.e0.append(float(line_trimmed[-1]))
                    # print(line_trimmed)
                else:
                    pass

    # get the e-fermi
    def getEfermi(self):
        return self.e_fermi[-1]

    # get the Energy
    def getEnergy(self):
        return self.energy[-1]

    # get Energy without entropy
    def getEnergy_s(self):
        return self.energy_s[-1]

    # get Energy sigma -> 0
    def getE0(self):
        return self.e0[-1]


class LOCPOT(object):
    rawdata = []
    interval = 1000 # linenumber <1000
    def read(self,filedir,filename):
        with open(filedir+'/'+filename,'r') as f:
            i =1
            for line in f:
                i = i+1
                if line == '':
                    self.interval = i
                if i>self.interval:
                    line_t = line.strip().split()
                    self.rawdata.append(line_t)


class DOSCAR(object):
    rawdata = {}
    data_x = {}
    data_y = {}
    interval = 7
    nedos_step =301 #default step of NEDOS get from the NEDOS
    signal = 0
    d_center = 0
    def read(self,filedir,filename):
        with open(filedir+'/'+filename,'r') as f:
            i = 1
            line_a = []
            for line in f:
                i = i+1
                if i>self.interval:
                    line_t = line.strip().split()
                    # print(line_t)
                    line_a.append(line_t)
                    # print(line_a)
                    if (i-self.interval)%self.nedos_step == 0:
                        self.rawdata[self.signal] = line_a
                        # print(line_a)
                        # print(self.rawdata[1])
                        self.signal = self.signal+1
                        line_a = []


        # print(self.rawdata[1])
    def parse(self):
        data_x_t = []
        data_y_t = []
        for j in range(1,self.signal):
            for i in self.rawdata[j]:
                data_x_t.append(float(i[0]))
                data_y_t.append(float(i[1]))
            self.data_x[j] = data_x_t
            self.data_y[j] = data_y_t
            data_x_t = []
            data_y_t = []
    def test(self):
        # N = 1000
        # random_x = np.random.randn(N)
        # random_y = np.random.randn(N)
        x_t  ={}
        y_t = {}
        trace_t = []
        for i in range(1,self.signal-1):
            x_t[i] = self.data_x[i]
            y_t[i] = self.data_y[i]
        # Create a trace
        for i in range(1,self.signal-1):
            trace = go.Scatter(
                x=x_t[i],
                y=y_t[i],
                mode='lines+markers'
            )
            trace_t.append(trace)

        data = trace_t

        # Plot and embed in ipython notebook!
        plot(data, filename='basic-scatter')

        # or plot with: plot_url = py.plot(data, filename='basic-line')

    def dcenter(self,start,end):
        # intergrate from start to end [start ,end]
        # calculate the dos intergrate and the dcenter
        intergrate_x = np.arange(start,end)
        dcenter_datax = np.array(self.data_x[1])
        dcenter_datay = np.array(self.data_y[1])
        # dcenter_dataxy = dcenter_datax*dcenter_datay
        # dcenter_dataxy = intergrate_x*dcenter_datay


        dcenter_devided = simps(self.data_y[1],intergrate_x)
        # dcenter_devide = simps(dcenter_dataxy,self.data_x[1])
        # dcenter_result = dcenter_devide/dcenter_devided

        # x = np.arange(-25,7)
        print(dcenter_devided)# intergrate from the begin to end
