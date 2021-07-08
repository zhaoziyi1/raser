"""

Author: Tao Yang <yangtao@ihep.ac.cn> 
Created [2021-07-07 Wed 08:47] 
Based on Raser C++ https://github.com/dt-np/raser


This program solves 2D Poisson's equation for SiC PIN & LGAD doping profile

    - div grad u(x,y) = f(x,y)

on the unit interval with source f given by

    f(x,y) = e0*Q(x,y)/perm0/perm

    PIN:
    f(x,y) = e0*1e13*1e6/perm0/perm_sic

    LGAD:
    f(x,y) = (6.1942*pow(10,16)/sqrt(2*3.1415926)/0.13821*exp(-pow(x-0.67,2))/2/0.13821/0.13821 + pow(10,13))*(-(1.60217733e-19)*(1e6)/(8.854187817e-12)/(9.76))

and boundary conditions given by
    u(x,y) = 0 for x = 0
    u(x,y) = Bias_voltage for x = thin

"""


import fenics
import matplotlib.pyplot as plt
import numpy as np
import ROOT



""" Global Constant """
E0 = 1.60217733e-19 # C
PERM0 = 8.854187817e-12 #F/m



""" Define Material """

class Material:

    def __init__(self,mat_name):
        self.mat_name = mat_name

    def mat_database(self):

        # material par
        self.si_par_dict = {'Permittivity' : 11.5,\
                             'Avalanche': 'unknown',\
                             'Mobility' : 'unknown'\
                            }

        self.sic_par_dict = {'Permittivity' : 9.76,\
                             'Avalanche': 'unknown',\
                             'Mobility' : 'unknown'\
                            }

        # global data base
        self.mat_db_dict = {'SiC' : self.sic_par_dict,\
                            'Si' : self.si_par_dict\
                            }

        return self.mat_db_dict[self.mat_name]



""" Define Detector Geometry """

class R2dDetector:

    def __init__(self,det_width,det_thin):
        self.det_width = det_width #um X-axis
        self.det_thin = det_width #um Y-axis
    
    def mesh(self,x_step,y_step):
        self.nx = int(self.det_width/x_step)
        self.ny = int(self.det_thin/y_step)

    def set_mat(self,mat):
        self.mat =mat
        self.par_dict = self.mat.mat_database()

    def set_doping(self,doping_epr):
        self.doping_epr = doping_epr

    def set_bias_voltage(self,bias_voltage):
        self.bias_voltage = bias_voltage


class FenicsPossion:

    def __init__(self,det):
        
        self.det = det

        # poential & field
        self.potential_value_2d = []
    
        self.electric_field_x_value = [ [] for n in range(self.det.ny+1) ]
        self.electric_field_y_value = [ [] for n in range(self.det.ny) ]

        self.electric_field_x_position = [ [] for n in range(self.det.ny+1) ]
        self.electric_field_y_position = [ [] for n in range(self.det.ny) ]

        # weighting poential & weighting field
        self.weighting_potential_value_2d = []
    
        self.weighting_electric_field_x_value = [ [] for n in range(self.det.ny+1) ]
        self.weighting_electric_field_y_value = [ [] for n in range(self.det.ny) ]

        self.weighting_electric_field_x_position = [ [] for n in range(self.det.ny+1) ]
        self.weighting_electric_field_y_position = [ [] for n in range(self.det.ny) ]

    def cal_possion(self):
        
        width = self.det.det_width
        thin = self.det.det_thin
        
        nx = self.det.nx
        ny = self.det.ny

        # Create mesh and function space
        mesh = fenics.RectangleMesh(fenics.Point(0, 0), fenics.Point(width, thin), nx, ny)
        V = fenics.FunctionSpace(mesh, "P", 1)

        # Define boundary condition
        u_D = fenics.Expression('x[1] < tol? det_voltage : 0', degree = 2,tol = 1E-14,det_voltage = self.det.bias_voltage)

        def boundary(x, on_boundary):
            return abs(x[1])<1E-14 or abs(x[1]-thin)<1E-14

        bc = fenics.DirichletBC(V, u_D, boundary)

        # Define variational problem
        u = fenics.TrialFunction(V)
        v = fenics.TestFunction(V)

        f = fenics.Expression(self.det.doping_epr,degree=2)
        a = fenics.dot(fenics.grad(u), fenics.grad(v))*fenics.dx
        L = f*v*fenics.dx #+ g*v*ds

        # Compute solution
        u = fenics.Function(V)
        fenics.solve(a == L, u, bc)

        potential_value_1d = u.compute_vertex_values()
        potential_value_2d = np.array(potential_value_1d).reshape(ny+1,nx+1)

        self.potential_value_2d = potential_value_2d

        # print(potential_value_2d)


    def cal_weighting_possion(self):

        width = self.det.det_width
        thin = self.det.det_thin
        
        nx = self.det.nx
        ny = self.det.ny

        # Create mesh and function space
        mesh = fenics.RectangleMesh(fenics.Point(0, 0), fenics.Point(width, thin), nx, ny)
        V = fenics.FunctionSpace(mesh, "P", 1)

        # Define boundary condition
        u_D = fenics.Expression('x[1] < tol? det_voltage : 0', degree = 2,tol = 1E-14,det_voltage = -1)

        def boundary(x, on_boundary):
            return abs(x[1])<1E-14 or abs(x[1]-thin)<1E-14

        bc = fenics.DirichletBC(V, u_D, boundary)

        # Define variational problem
        u = fenics.TrialFunction(V)
        v = fenics.TestFunction(V)

        f = fenics.Constant(0)
        a = fenics.dot(fenics.grad(u), fenics.grad(v))*fenics.dx
        L = f*v*fenics.dx #+ g*v*ds

        # Compute solution
        u = fenics.Function(V)
        fenics.solve(a == L, u, bc)

        weighting_potential_value_1d = u.compute_vertex_values()
        weighting_potential_value_2d = np.array(weighting_potential_value_1d).reshape(ny+1,nx+1)

        self.weighting_potential_value_2d = weighting_potential_value_2d

        return weighting_potential_value_2d

    def cal_electric_field(self):

        width = self.det.det_width
        thin = self.det.det_thin
        
        nx = self.det.nx
        ny = self.det.ny

        x_step = width/nx
        y_step = thin/ny

        # x direction
        for j in range(ny+1):
            for i in range(nx):

                # electric field
                tmp_xpos = 0.5*x_step*(2*i+1)
                tmp_ef = (self.potential_value_2d[j][i] - self.potential_value_2d[j][i+1])/x_step
                self.electric_field_x_position[j].append(tmp_xpos)
                self.electric_field_x_value[j].append(tmp_ef)

                # weighting field
                tmp_wxpos = 0.5*x_step*(2*i+1)
                tmp_wef = (self.weighting_potential_value_2d[j][i] - self.weighting_potential_value_2d[j][i+1])/x_step
                self.weighting_electric_field_x_position[j].append(tmp_wxpos)
                self.weighting_electric_field_x_value[j].append(tmp_wef)


        # y direction 
        for j in range(ny):
            for i in range(nx+1):

                # electric field
                tmp_ypos = 0.5*y_step*(2*j+1)
                tmp_ef = (self.potential_value_2d[j][i]- self.potential_value_2d[j+1][i])/y_step
                self.electric_field_y_position[j].append(tmp_ypos)
                self.electric_field_y_value[j].append(tmp_ef)

                # weighting field
                tmp_wypos = 0.5*y_step*(2*j+1)
                tmp_wef = (self.weighting_potential_value_2d[j][i] - self.weighting_potential_value_2d[j+1][i])/y_step
                self.weighting_electric_field_y_position[j].append(tmp_wypos)
                self.weighting_electric_field_y_value[j].append(tmp_wef)     


    def solve(self):

        self.cal_possion()
        self.cal_weighting_possion()
        self.cal_electric_field()

    def draw(self):

        # plot electric field at x = middle
        # ep = [depth[250] for depth in self.electric_field_y_position]
        # ev = [ef[250] for ef in self.electric_field_y_value]
        # g_e = ROOT.TGraph(len(ep),ep,ev)
        
        # plot weighting electric field at x = middle
        # wep = [depth[250] for depth in self.weighting_electric_field_y_position]
        # wev = [wef[250] for wef in self.weighting_electric_field_y_value]
        # g_we = ROOT.TGraph(len(wep),wep,wev)

        # c = ROOT.TCanvas( 'c', 'c', 500, 500, 500, 500 )

        # c.Divide(1,2)

        # c.cd(1)
        # g_e.Draw()

        plt.figure(figsize=(4,4), dpi=80)
        plt.figure(1)
                
        plt.plot( [depth[250] for depth in self.electric_field_y_position], [ef[250] for ef in self.electric_field_y_value])
        plt.xlabel('depth [um]')
        plt.ylabel('electric field [V/um]')
        plt.tight_layout()
        plt.show()
        

    
if __name__ == '__main__':

    my_lgad = R2dDetector(50,50)

    my_sic_mat = Material('SiC')
    my_lgad.set_mat(my_sic_mat)

    my_lgad.set_doping("(6.1942*(1e14)/sqrt(2*3.1415926)/0.13821*exp(-pow(x[1]-0.67,2))/2/0.13821/0.13821 + (1e13))*((1.60217733e-19)*(1e6)/(8.854187817e-12)/(9.76))*1e-12")

    my_lgad.set_bias_voltage(-200)

    my_lgad.mesh(0.1, 0.1)

    my_possion_solver = FenicsPossion(my_lgad)
    my_possion_solver.solve()
    my_possion_solver.draw()



# perm_sic = 9.76  #Permittivity
# e0 = 1.60217733e-19
# perm0 = 8.854187817e-12   #F/m
# 
# width = 50.0 # um X-axis
# thin = 50.0 # um Y-axis
# nx = 500
# ny = 500
# 
# bias_voltage = -200 # V
# 
# # Create mesh and function space
# 
# mesh = fenics.RectangleMesh(fenics.Point(0, 0), fenics.Point(width, thin), nx, ny)
# V = fenics.FunctionSpace(mesh, "P", 1)
# 
# # Define variational problem
# u = fenics.TrialFunction(V)
# v = fenics.TestFunction(V)
# 
# # LGAD
# # f = Expression("(6.1942*(1e14)/sqrt(2*3.1415926)/0.13821*exp(-pow(x[1]-0.67,2))/2/0.13821/0.13821 + # (1e13))*((1.60217733e-19)*(1e6)/(8.854187817e-12)/(9.76))*1e-12", degree=2)
# 
# # PIN
# value = e0*10*1e6/perm0/perm_sic
# f= fenics.Constant(value)
# 
# a = fenics.dot(fenics.grad(u), fenics.grad(v))*dx
# L = f*v*fenics.dx #+ g*v*ds
# 
# 
# # Define boundary condition
# u_D = fenics.Expression('x[1] < tol? det_voltage : 0', degree=2,tol=1E-14,det_voltage=bias_voltage)
# 
# def boundary(x, on_boundary):
#     return abs(x[1])<1E-14 or abs(x[1]-thin)<1E-14
# 
# bc = fenics.DirichletBC(V, u_D, boundary)
# 
# # Compute solution
# u = fenics.Function(V)
# fenics.solve(a == L, u, bc)
# 
# # Save solution to file
# #file = File("poisson.pvd")
# #file << u
# 
# # Calculate electric field
# potential_value_1d = u.compute_vertex_values()
# potential_value_2d = np.array(potential_value_1d).reshape(nx+1,ny+1)
# # print(potential_value_2d)
# 
# 
# potential_value = []
# potential_position = []
# 
# for i in range(ny+1):
#     potential_value.append(potential_value_2d[i][250])
#     potential_position.append(i*thin/ny)
# 
# 
# electric_field_value = []
# electric_field_position = []
# 
# for i in (range(ny)):
#     tmp_e = (potential_value_2d[i][250]-potential_value_2d[i+1][250])/(50.0/ny)
#     electric_field_value.append(tmp_e)
#     electric_field_position.append(thin/ny*(i+1))
# 
# #print(len(electric_field_value))
# 
# 
# 
# # Plot solution
# depth = np.linspace(0, 50, 500)
# doping = (6.1942*(1e14)/sqrt(2*3.1415926)/0.13821*np.exp(-pow(depth-0.67,2))/2/0.13821/0.13821 + (1e13))
# 
# #plot(u)
# # plt.figure(figsize=(4,4), dpi=80)
# # plt.plot(depth,doping)
# # plt.yscale('symlog')
# # plt.xlabel('depth [um]')
# # plt.ylabel('doping [cm-3]')
# 
# plt.figure(figsize=(12,4), dpi=80)
# plt.figure(1)
# 
# ax1 = plt.subplot(131)
# ax1.invert_yaxis()
# plt.xlabel('width [um]')
# plt.ylabel('depth [um]')
# plot(u)
# 
# ax2 = plt.subplot(132)
# plt.plot(potential_position,potential_value)
# plt.xlabel('depth [um]')
# plt.ylabel('potential [V]')
# 
# ax2 = plt.subplot(133)
# plt.plot(electric_field_position,electric_field_value)
# plt.xlabel('depth [um]')
# plt.ylabel('electric field [V/um]')
# 
# 
# plt.tight_layout()
# 
# plt.show()