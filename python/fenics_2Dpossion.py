from fenics import *
import matplotlib.pyplot as plt
import ROOT
from array import array
import sys
#link raser and python fenics
args = sys.argv[1:]
SetUpVolume=float(args[0])  #um
detector_y=float(args[1])
det_voltage=float(args[2])
Neff=-float(args[3])#/cm3   //n type is negetive and p type is positive
out_root=args[4]#/cm3 
print(args)
nx=1
perm_sic=9.76
detector_x=100.0
ny=int(detector_y/SetUpVolume)
#####Poisson equationc

#### constant value
e0=1.60217733e-19
perm0=8.854187817e-12   #F/m
# Create mesh and define function space
# divide n parts

mesh = RectangleMesh(Point(0, 0), Point(detector_x, detector_y), nx, ny)
V = FunctionSpace(mesh, 'P', 1)
# Define boundary condition
tol=1E-14
u_D = Expression('x[1] < tol? det_voltage : 0', degree=2,tol=1E-14,det_voltage=det_voltage)
def boundary(x, on_boundary):
    return abs(x[1])<tol or abs(x[1]-detector_y)<tol

bc = DirichletBC(V, u_D, boundary)
# # Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
#possion
f_value=-e0*Neff*1e6/perm0/perm_sic
# print("f_value:",f_value)
f = Constant(f_value)
a = dot(grad(u), grad(v))*dx
L = f*v*dx

# # Compute solution
u = Function(V)
solve(a == L, u, bc)

#####Laplace's equation
u_w_D = Expression('x[1] < tol? 0 : 1', degree=2,tol=1E-14)

def boundary_w(x, on_boundary):
    return abs(x[1])<tol or abs(x[1]-detector_y)<tol
bc_w = DirichletBC(V, u_w_D, boundary_w)

# # Define variational problem
u_w = TrialFunction(V)
v_w = TestFunction(V)

f_w = Constant(0)
a_w = dot(grad(u_w), grad(v_w))*dx
L_w = f_w*v_w*dx

# # Compute solution
u_w = Function(V)
solve(a_w == L_w, u_w, bc_w)

# write data to the root file
Events=array('i',[0])
sen_depth=array('d',[0])
E_ele=array('d',[0])
P_w_ele=array('d',[0])
P_ele=array('d',[0])
Events[0]=0
sen_depth[0]=0.0
E_ele[0]=0.0
P_w_ele[0]=0.0
P_ele[0]=0.0

outfile="python/"+out_root
out_file=ROOT.TFile(outfile,"RECREATE")
tree_out=ROOT.TTree('tree','tree')
tree_out.Branch('Events',Events,'Events/I')
tree_out.Branch('sen_depth',sen_depth,'sen_depth/D')
tree_out.Branch('E_ele',E_ele,'E_ele/D')
tree_out.Branch('P_w_ele',P_w_ele,'P_w_ele/D')
tree_out.Branch('P_ele',P_ele,'P_ele/D')

number_events=0
#plot potential
vertex_values_u = u.compute_vertex_values()
vertex_values_u_w = u_w.compute_vertex_values()
x_value=[]
y_value=[]
y_w_value=[]
Ex_value=[]
E_value=[]
E_w_value=[]
l_step=detector_y/(len(vertex_values_u)/2-1)
l_start=0.0
#calculate the electric field and weighting potential
for i in range(0,int(len(vertex_values_u)/2)-1):
    number_events+=1

    l_x=l_start+i*l_step
    x_value.append(l_x)
    y_value.append(vertex_values_u[2*i])
    y_w_value.append(vertex_values_u_w[2*i])  #weighting potenital
    x_E_f=(vertex_values_u[2*i+2]-vertex_values_u[2*i])/l_step
    x_w_E=(vertex_values_u_w[2*i+2]-vertex_values_u_w[2*i])/l_step 
    E_value.append(x_E_f)  #electric field
    E_w_value.append(x_w_E) # weighting potential

    #fill to root
    Events[0]=number_events
    sen_depth[0]=x_value[i]
    E_ele[0]=x_E_f
    P_w_ele[0]=y_w_value[i]
    P_ele[0]=y_value[i]
    tree_out.Fill()
tree_out.Write()
out_file.Close()
