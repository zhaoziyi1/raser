import matplotlib.pyplot as plt
from array import array
import fenics
import mshr
import numpy as np
import math
import random
import ROOT
import time
import sys
import os
FACTOR = 1000.0 #um->mm
#define the detector parameter 
class R3dDetector:
    def __init__(self,l_x,l_y,l_z):
        self.l_x = l_x / FACTOR #size
        self.l_y = l_y / FACTOR
        self.l_z = l_z / FACTOR
    def set_para(self,doping,voltage,temperature):
        self.d_neff = doping #dopingX1e12 cm^-3
        self.v_voltage = voltage #Voltage
        self.temperature = temperature #Voltage
        self.n_bin = 1000
        self.t_end = 3e-9
        self.positive_cu = ROOT.TH1F("charge+","Positive Current",self.n_bin,0,self.t_end)
        self.negtive_cu = ROOT.TH1F("charge-","Negative Current",self.n_bin,0,self.t_end)
        self.sum_cu = ROOT.TH1F("charge","Total Current",self.n_bin,0,self.t_end)
    def set_electrode(self):
        e_r = 15.0 / FACTOR
        e_int = 120.0 / FACTOR
        e_t_y = R3dDetector.infor_ele(self,e_r,e_int)
        self.e_tr=[]
        self.e_t_1 = [self.l_x*0.5          ,self.l_y*0.5      ,e_r,0,self.l_z,"p"]
        self.e_t_2 = [self.l_x*0.5-e_int    ,self.l_y*0.5      ,e_r,0,self.l_z,"n"]
        self.e_t_3 = [self.l_x*0.5+e_int    ,self.l_y*0.5      ,e_r,0,self.l_z,"n"]
        self.e_t_4 = [self.l_x*0.5-e_int*0.5,self.l_y*0.5+e_t_y,e_r,0,self.l_z,"n"]
        self.e_t_5 = [self.l_x*0.5+e_int*0.5,self.l_y*0.5+e_t_y,e_r,0,self.l_z,"n"]
        self.e_t_6 = [self.l_x*0.5-e_int*0.5,self.l_y*0.5-e_t_y,e_r,0,self.l_z,"n"]
        self.e_t_7 = [self.l_x*0.5+e_int*0.5,self.l_y*0.5-e_t_y,e_r,0,self.l_z,"n"]
        for i in range(7):
           n_e = eval('self.e_t_' + str(i+1))
           self.e_tr.append(n_e)
    def infor_ele(self,e_r,e_int):
        e_x_gap = self.l_x - 2*e_r - 2*e_int
        if e_x_gap < 0:
            print("the electrode at x position is large than sensor length")
            sys.exit(0)
        e_t_y = math.sqrt(e_int*e_int*0.75)
        if 2*e_t_y > self.l_y:
            print("the electrode at y position is large than sensor length")
            sys.exit(0)            
        return e_t_y

#Calculate the weighting potential and electric field
class Fenics_cal:
    #parameter of SiC
    def __init__(self,my_d,mesh_v):
        self.p_electric = []
        self.w_p_electric = []
        perm_sic = 9.76  #Permittivity
        e0 = 1.60217733e-19
        perm0 = 8.854187817e-12   #F/m
        self.f_value = -e0*my_d.d_neff*1e6/perm0/perm_sic*FACTOR*FACTOR
        self.tol = 1e-14
        #fenics space        
        m_sensor =  mshr.Box(fenics.Point(0, 0, 0), fenics.Point(my_d.l_x, my_d.l_y, my_d.l_z))
        for i in range(len(my_d.e_tr)):
            e_t_i = my_d.e_tr[i]
            elec_n=mshr.Cylinder(fenics.Point(e_t_i[0], e_t_i[1], e_t_i[3]), fenics.Point(e_t_i[0], e_t_i[1], e_t_i[4]),e_t_i[2],e_t_i[2])
            m_sensor =m_sensor - elec_n
        self.mesh3D = mshr.generate_mesh(m_sensor,mesh_v)      
        self.V = fenics.FunctionSpace(self.mesh3D, 'P', 1)

    def fenics_p_electric(self,my_d):    #get the electric potential
        # Define boundary condition
        bc_u=[]
        for i in range (len(my_d.e_tr)):
            e_t_i = my_d.e_tr[i]
            str_e =  "x[0]>={elec_1_0}-{elec_1_2} && x[0]<={elec_1_0}+"+\
                "{elec_1_2} && x[1]>={elec_1_1}-{elec_1_2} && "+\
                "x[1]<={elec_1_1}+{elec_1_2} && x[2]>={elec_1_3} && x[2]<={elec_1_4} && on_boundary"
            elec_p = str_e.format(elec_1_0=e_t_i[0],elec_1_1=e_t_i[1],elec_1_2=e_t_i[2],elec_1_3=e_t_i[3],elec_1_4=e_t_i[4])
            if e_t_i[5] == "p":
                bc = fenics.DirichletBC(self.V, my_d.v_voltage, elec_p)
            else:
                bc = fenics.DirichletBC(self.V, 0.0, elec_p)
            bc_u.append(bc)
        # # Define variational problem
        u = fenics.TrialFunction(self.V)
        v = fenics.TestFunction(self.V)
        f = fenics.Constant(self.f_value)
        a = fenics.dot(fenics.grad(u), fenics.grad(v))*fenics.dx
        L = f*v*fenics.dx
        # # Compute solution
        self.u = fenics.Function(self.V)
        fenics.solve(a == L, self.u, bc_u,solver_parameters=dict(linear_solver='gmres', preconditioner='ilu'))
        W = fenics.VectorFunctionSpace(self.mesh3D, 'P', 1)
        self.E_field = fenics.project(fenics.as_vector((self.u.dx(0),self.u.dx(1),self.u.dx(2))),W)

    def fenics_p_w_electric(self,my_d):  #get the electric weighting potential
        #####Laplace's equation
        bc_w=[]
        for i in range (len(my_d.e_tr)):
            e_t_i = my_d.e_tr[i]
            str_e =  "x[0]>={elec_1_0}-{elec_1_2} && x[0]<={elec_1_0}+"+\
                "{elec_1_2} && x[1]>={elec_1_1}-{elec_1_2} && "+\
                "x[1]<={elec_1_1}+{elec_1_2} && x[2]>={elec_1_3} && x[2]<={elec_1_4} && on_boundary"
            elec_p = str_e.format(elec_1_0=e_t_i[0],elec_1_1=e_t_i[1],elec_1_2=e_t_i[2],elec_1_3=e_t_i[3],elec_1_4=e_t_i[4])
            if e_t_i[5] == "p":
                bc = fenics.DirichletBC(self.V, 0.0, elec_p)
            else:
                bc = fenics.DirichletBC(self.V, 1.0, elec_p)
            bc_w.append(bc)
        # # Define variational problem
        u_w = fenics.TrialFunction(self.V)
        v_w = fenics.TestFunction(self.V)
        f_w = fenics.Constant(0)
        a_w = fenics.dot(fenics.grad(u_w), fenics.grad(v_w))*fenics.dx
        L_w = f_w*v_w*fenics.dx
        # # Compute solution
        self.u_w = fenics.Function(self.V)
        fenics.solve(a_w == L_w, self.u_w, bc_w)

    def get_e_field(self,px,py,pz):

        try:
            x_value,y_value,z_value = self.E_field(px,py,pz)
        except RuntimeError:
            x_value,y_value,z_value = 0,0,0
        return x_value,y_value,z_value

    def get_w_p(self,px,py,pz):
        try:
            f_w_p = self.u_w(px,py,pz)
        except RuntimeError:
            f_w_p = 0.0
        return f_w_p

#####The tracks of the different incident paticle (Mip or laser or others)
class Tracks:
    def t_mip(self,p_entry,p_exit,n_div):
        self.p_entry = p_entry
        self.p_exit = p_exit
        self.p_tracks = [ [] for n in range(n_div-1) ]
        self.d_tracks = []
        p_x=p_entry[0]
        p_y=p_entry[1]
        p_z=p_entry[2]
        # self.py_entry=[]
        for i in range(n_div-1):
            x_div_point = (p_exit[0]-p_entry[0])/n_div*i+p_entry[0]+(p_exit[0]-p_entry[0])/(2*n_div)
            y_div_point = (p_exit[1]-p_entry[1])/n_div*i+p_entry[1]+(p_exit[1]-p_entry[1])/(2*n_div)
            z_div_point = (p_exit[2]-p_entry[2])/n_div*i+p_entry[2]+(p_exit[2]-p_entry[2])/(2*n_div)
            self.p_tracks[i].append(x_div_point)
            self.p_tracks[i].append(y_div_point)
            self.p_tracks[i].append(z_div_point)
            self.d_tracks.append(math.sqrt(math.pow(x_div_point-p_x,2)+math.pow(y_div_point-p_y,2))+math.pow(z_div_point-p_z,2))
            p_x = x_div_point
            p_y = y_div_point
            p_z = z_div_point
            
# mobility model
def sic_mobility(charge,aver_e,my_detector):
    T=my_detector.temperature
    E=aver_e
    Neff=abs(my_detector.d_neff)
    if(charge>0):
        alpha = 0.34
        ulp = 124.0 * math.pow(T / 300.0, -2.0)
        uminp = 15.9
        Crefp = 1.76e19
        betap = 1.213 * math.pow(T / 300.0, 0.17)
        vsatp = 2e7 * math.pow(T / 300.0, 0.52)
        lfm = uminp + ulp/(1.0 + math.pow(Neff*1e12 / Crefp, alpha))
        hfm = lfm / (math.pow(1.0 + math.pow(lfm * E / vsatp, betap), 1.0 / betap))                        
    if(charge<0):
        alpha = 0.61
        ulp = 947.0 * math.pow(T / 300.0, -2)
        Crefp = 1.94e19
        betap = 1.0 * math.pow(T / 300.0, 0.66)
        vsatp = 2e7 * math.pow(T / 300.0, 0.87)
        lfm = ulp/ (1.0 + math.pow(Neff*1e12 / Crefp, alpha))
        hfm = lfm / (math.pow(1.0 + math.pow(lfm * E / vsatp, betap), 1.0/betap))                      
    return hfm
    
#The drift of generated particles
class Drifts:
    def __init__(self,my_track):
        self.muhh=1650.0   #mobility related with the magnetic field (now silicon useless)
        self.muhe=310.0
        self.BB=np.array([0,0,0])
        self.sstep=0.1/FACTOR #drift step
        self.kboltz=8.617385e-5 #eV/K
        self.max_drift_len=1e9/FACTOR #maximum diftlenght [um]
        self.d_dic_n = {}
        self.d_dic_p = {}
        for n in range(len(my_track.p_tracks)):
            self.d_dic_n["tk_"+str(n+1)] = [ [] for n in range(5) ]
            self.d_dic_p["tk_"+str(n+1)] = [ [] for n in range(5) ]
    def initial_parameter(self):
        self.end_cond=0
        self.d_time=0
        self.path_len=0
        self.n_step=0
        self.charge=0

    def delta_p(self):
        # magnetic field effect
        if(self.charg)>0:
            FF=self.e_field+self.muhh*np.cross(self.e_field,self.BB)
        else:
            FF=self.e_field-self.muhe*np.cross(self.e_field,self.BB)
        #the delta x with electric field
        if(np.linalg.norm(FF)!=0):
            self.delta_x=-self.sstep*self.charg*FF[0]/np.linalg.norm(FF)
            self.delta_y=-self.sstep*self.charg*FF[1]/np.linalg.norm(FF)
            self.delta_z=-self.sstep*self.charg*FF[2]/np.linalg.norm(FF)
        else:
            self.delta_x=0
            self.delta_y=0
            self.delta_z=0

    def drift_v(self,my_d,my_f):
        e_delta_f = np.array(my_f.get_e_field(self.d_x+self.delta_x,self.d_y+self.delta_y,self.d_z+self.delta_z))
        aver_e=(np.linalg.norm(self.e_field)+np.linalg.norm(e_delta_f))/2.0*1e4/FACTOR  #V/cm
        self.v_drift=sic_mobility(self.charg,aver_e,my_d)*aver_e  # mobility cm2/(V s) v : cm/s
        #drift part
        if(self.v_drift==0):
            self.delta_x=0.0
            self.delta_y=0.0
            self.delta_z=0.0
            self.dif_x=0.0
            self.dif_y=0.0
            self.dif_z=0.0
            self.end_cond=9
        else:
            #off when the field gets large enough
            DiffOffField=100.0*FACTOR  # the silicon value ???              
            if(np.linalg.norm(e_delta_f)<DiffOffField):
                self.s_time=self.sstep*1e-4*FACTOR/self.v_drift
                s_sigma=math.sqrt(2.0*self.kboltz*sic_mobility(self.charg,aver_e,my_d)*my_d.temperature*self.s_time)
                self.dif_x=random.gauss(0.0,s_sigma)*1e4/FACTOR
                self.dif_y=random.gauss(0.0,s_sigma)*1e4/FACTOR
                self.dif_z=random.gauss(0.0,s_sigma)*1e4/FACTOR          
            else:
                self.dif_x=0.0
                self.dif_y=0.0
                self.dif_z=0.0

    def drift_s_step(self,my_d):
        # x axis   
        if((self.d_x+self.delta_x+self.dif_x)>=my_d.l_x): 
            self.d_cx = my_d.l_x
        elif((self.d_x+self.delta_x+self.dif_x)<0):
            self.d_cx = 0
        else:
            self.d_cx = self.d_x+self.delta_x+self.dif_x
        # y axis
        if((self.d_y+self.delta_y+self.dif_y)>=my_d.l_y): 
            self.d_cy = my_d.l_y
        elif((self.d_y+self.delta_y+self.dif_y)<0):
            self.d_cy = 0
        else:
            self.d_cy = self.d_y+self.delta_y+self.dif_y
        # z axis
        if((self.d_z+self.delta_z+self.dif_z)>=my_d.l_z): 
            self.d_cz = my_d.l_z
        elif((self.d_z+self.delta_y+self.dif_z)<0):
            self.d_cz = 0
        else:
            self.d_cz = self.d_z+self.delta_z+self.dif_z

    def drift_end_condition(self):    
        if(self.wpot>(1-1e-5)):
            self.end_cond=1
        if(self.d_x<=0):
            self.end_cond=2
        if(self.d_y<=0):
            self.end_cond=4
        if(self.d_z<=0):
            self.end_cond=5
        if(self.path_len>self.max_drift_len):
            self.end_cond=6
        if(self.n_step>10000):
            self.end_cond=7
    def save_inf_track(self):
        e0 = 1.60217733e-19
        if(self.charge<0):
            if(self.charg>0):
                self.d_dic_p["tk_"+str(self.n_track)][0].append(self.d_x)
                self.d_dic_p["tk_"+str(self.n_track)][1].append(self.d_y)
                self.d_dic_p["tk_"+str(self.n_track)][2].append(self.d_z)
                self.d_dic_p["tk_"+str(self.n_track)][3].append(self.charge)
                self.d_dic_p["tk_"+str(self.n_track)][4].append(self.d_time)
            else:
                self.d_dic_n["tk_"+str(self.n_track)][0].append(self.d_x)
                self.d_dic_n["tk_"+str(self.n_track)][1].append(self.d_y)
                self.d_dic_n["tk_"+str(self.n_track)][2].append(self.d_z)
                self.d_dic_n["tk_"+str(self.n_track)][3].append(self.charge)
                self.d_dic_n["tk_"+str(self.n_track)][4].append(self.d_time)

    def cal_current(self,my_d,my_t):
        my_d.positive_cu.Reset()
        my_d.negtive_cu.Reset()
        my_d.sum_cu.Reset()
        self.sum_p_current = []
        test_p = ROOT.TH1F("test+","test+",my_d.n_bin,0,my_d.t_end)
        test_n = ROOT.TH1F("test-","test-",my_d.n_bin,0,my_d.t_end)
        test_sum = ROOT.TH1F("test sum","test sum",my_d.n_bin,0,my_d.t_end)
        total_pairs=0
        for j in range(len(my_t.p_tracks)):
            for i in range(len(self.d_dic_p["tk_"+str(j+1)][2])):
                test_p.Fill(self.d_dic_p["tk_"+str(j+1)][4][i],self.d_dic_p["tk_"+str(j+1)][3][i])
            test_p = Drifts.get_current_his(self,test_p)           
            for i in range(len(self.d_dic_n["tk_"+str(j+1)][2])):
                test_n.Fill(self.d_dic_n["tk_"+str(j+1)][4][i],self.d_dic_n["tk_"+str(j+1)][3][i])
            test_n = Drifts.get_current_his(self,test_n)
            n_pairs=Drifts.get_energy_loss(self,my_t.d_tracks[j])
            total_pairs+=n_pairs
            test_p.Scale(n_pairs)
            test_n.Scale(n_pairs)            
            my_d.positive_cu.Add(test_p)
            my_d.negtive_cu.Add(test_n)
            test_p.Reset()
            test_n.Reset()
        laudau_t_pairs = Drifts.get_laudau_dis(self,my_t)
        n_scale = laudau_t_pairs/total_pairs
        my_d.positive_cu.Scale(n_scale)
        my_d.negtive_cu.Scale(n_scale)
        my_d.sum_cu.Add(my_d.positive_cu)
        my_d.sum_cu.Add(my_d.negtive_cu)
    def get_current_his(self,histo):
        e0 = 1.60217733e-19
        hist = ROOT.TH1F()
        histo.Copy(hist)
        for i in range(hist.GetNbinsX()):
            histo.SetBinContent(i, hist.GetBinContent(i) \
                /((hist.GetXaxis().GetXmax() - hist.GetXaxis().GetXmin()) \
                / hist.GetNbinsX())*e0)
        return histo

    def get_energy_loss(self,d_drift):
        #Drift distance: d_drift
        loss_energy = 8.4 #silicon carbide ev /one pairs
        Myf = ROOT.TFile("SiC_5um.root")
        hist = Myf.Get("Edep_device")
        gRandom = ROOT.TRandom3(0)
        ran_energy = hist.GetRandom(gRandom)
        n_pairs = d_drift*ran_energy*1e6/(5*loss_energy)*FACTOR
        return n_pairs

    def get_laudau_dis(self,my_t):
        loss_energy = 8.4
        d_x = math.pow(my_t.p_exit[0]-my_t.p_entry[0],2)
        d_y = math.pow(my_t.p_exit[1]-my_t.p_entry[1],2)
        d_z = math.pow(my_t.p_exit[2]-my_t.p_entry[2],2)
        d_dis = math.sqrt(d_x+d_y+d_z)*FACTOR
        d_mpv = 0.04 * math.log(d_dis) + 0.27
        d_FWHM = 0.31 * math.pow(d_dis, -0.17)
        d_da = d_FWHM / 4.
        gRandom = ROOT.TRandom3(0)
        LanConst = gRandom.Landau(d_mpv, d_da)
        if (LanConst > 5.0 * d_mpv):
            LanConst = gRandom.Landau(d_mpv, d_da)
        Lan_pairs = LanConst*1000.0/loss_energy*d_dis
        return Lan_pairs

    def ionized_drift(self,my_t,my_f,my_d):       
        for i in range(len(my_t.p_tracks)):
            self.n_track=i+1
            for j in range(2):
                if (j==0):
                    self.charg=1 #hole
                if (j==1):
                    self.charg=-1 #electron
                Drifts.initial_parameter(self)
                #generated particles positions
                self.d_x=my_t.p_tracks[i][0]#initial position
                self.d_y=my_t.p_tracks[i][1]
                self.d_z=my_t.p_tracks[i][2]
                while (self.end_cond==0):
                    if (self.d_y>=(my_d.l_y-1.0/FACTOR) or self.d_x>=(my_d.l_x-1.0/FACTOR) or self.d_z>=(my_d.l_z-1.0/FACTOR)):
                        self.end_cond=3  
                    else:                                       
                        # field of the position
                        self.e_field = np.array(my_f.get_e_field(self.d_x,self.d_y,self.d_z))
                        # print(self.field)
                        #delta_poisiton
                        Drifts.delta_p(self)
                        #drift_position
                        Drifts.drift_v(self,my_d,my_f)
                        #drift_next_posiiton
                        Drifts.drift_s_step(self,my_d)
                        #charge_collection
                        delta_Ew=my_f.get_w_p(self.d_cx,self.d_cy,self.d_cz)-my_f.get_w_p(self.d_x,self.d_y,self.d_cz)
                        self.charge=self.charg*delta_Ew
                        if(self.v_drift!=0):
                            self.d_time=self.d_time+self.sstep*1e-4*FACTOR/self.v_drift
                            self.path_len+=self.sstep
                        self.d_x=self.d_cx
                        self.d_y=self.d_cy
                        self.d_z=self.d_cz
                        self.wpot=my_f.get_w_p(self.d_x,self.d_y,self.d_z)
                        Drifts.save_inf_track(self)  
                        Drifts.drift_end_condition(self)                                                                         
                    self.n_step+=1
        Drifts.cal_current(self,my_d,my_t)

    def draw_drift_path(self,my_d,my_f,n_e_v):
        ROOT.gStyle.SetOptStat(0)
        # # ROOT.gROOT.SetBatch(1)
        c1 = ROOT.TCanvas("c", "canvas1", 200,10,1000, 1000)
        nx_e = n_e_v[0]
        ny_e = n_e_v[1]
        nz_e = n_e_v[2]
        structrue = ROOT.TH3D("","",nx_e,0,my_d.l_x,ny_e,0,my_d.l_y,nz_e,0,my_d.l_z)
        for k in range(nz_e):
            for j in range (ny_e):
                for i in range(nx_e):
                    x_v = (i+1)*(my_d.l_x/nx_e)
                    y_v = (j+1)*(my_d.l_y/ny_e)
                    z_v = (k+1)*(my_d.l_z/nz_e)
                    try:
                        x_value,y_value,z_value = my_f.E_field(x_v,y_v,z_v)
                        structrue.SetBinContent(i+1,j+1,k+1,0)
                    except RuntimeError:
                        structrue.SetBinContent(i+1,j+1,k+1,2)
        structrue.Draw("ISO")

        x_array=array('f')
        y_array=array('f')
        z_array=array('f')
        for i in range(len(self.d_dic_p)):
            n=len(self.d_dic_p["tk_"+str(i+1)][0])
            if(n>0):
                x_array.extend(self.d_dic_p["tk_"+str(i+1)][0])
                y_array.extend(self.d_dic_p["tk_"+str(i+1)][1]) 
                z_array.extend(self.d_dic_p["tk_"+str(i+1)][2])              
                gr_p = ROOT.TPolyLine3D(n,x_array,y_array,z_array)
                gr_p.SetLineColor(4)
                gr_p.SetLineStyle(1)
                gr_p.Draw("SAME")
                del x_array[:]
                del y_array[:]
                del z_array[:]
        for j in range(len(self.d_dic_n)):
            m=len(self.d_dic_n["tk_"+str(j+1)][0])
            if(m>0):
                x_array.extend(self.d_dic_n["tk_"+str(j+1)][0])
                y_array.extend(self.d_dic_n["tk_"+str(j+1)][1])
                z_array.extend(self.d_dic_n["tk_"+str(j+1)][2])                
                gr_n = ROOT.TPolyLine3D(m,x_array,y_array,z_array)
                gr_n.SetLineColor(2)
                gr_n.SetLineStyle(1)
                gr_n.Draw("SAME")
                del x_array[:]
                del y_array[:]
                del z_array[:]
        c1.SaveAs("fig/drift_path.root")
        del c1

class Amplifier:
    def CSA_amp(self,my_d,t_rise,t_fall,trans_imp):
        
        hist = ROOT.TH1F()
        my_d.sum_cu.Copy(hist)
        max_num = hist.GetNbinsX()
        preamp_Q = [0.0]*max_num
        itot = [0.0]*max_num
        shaper_out_Q = [0.0]*max_num
        shaper_out_V = [0.0]*max_num
        qtot = 0.0
        sh_max = 0.0

        tau_rise = t_rise / 2.2*1e-9
        tau_fall = t_fall / 2.2*1e-9   
        if (tau_rise == tau_fall):
            tau_rise *= 0.9
        time_unit = (hist.GetXaxis().GetXmax()- hist.GetXaxis().GetXmin()) / hist.GetNbinsX()
        for i in range(max_num-1):
            if (i>0):
                preamp_Q[i] = 0.0
                itot[i] = hist.GetBinContent(i)
                preamp_Q[i] = itot[i] * time_unit
                qtot += itot[i] *time_unit
            elif (i==0):
                preamp_Q[i]==0.0
        for i in range(max_num-1):
            dif_shaper_Q = preamp_Q[i]
            if(dif_shaper_Q != 0):
                for j in range(max_num-i):
                    shaper_out_Q[i+j] += tau_fall/(tau_fall+tau_rise)*dif_shaper_Q \
                        *(math.exp(-j*time_unit/tau_fall)-math.exp(-j*time_unit/tau_rise))
            if (abs(shaper_out_Q[i]) > abs(sh_max)):
                sh_max = shaper_out_Q[i]
        for i in range(max_num):
            shaper_out_V[i] = shaper_out_Q[i] * trans_imp * 1e15 *qtot
            hist.SetBinContent(i,shaper_out_V[i])

        charge_t=my_d.sum_cu.Integral() \
            * ((my_d.sum_cu.GetXaxis().GetXmax() \
            - my_d.sum_cu.GetXaxis().GetXmin()) \
            / my_d.sum_cu.GetNbinsX()) * 1e15
        print("total_charge:%s fc"%(qtot*1e15))
        print("total_charge:%s fc"%charge_t)
        return charge_t,qtot,hist

def draw_plot(my_detector,ele_current,qtot,my_drift,my_field,drift_path):

    ROOT.gStyle.SetOptStat(0)
    c = ROOT.TCanvas("c", "canvas", 200,10,1000, 1000)
    c.SetLeftMargin(0.12)
    # c.SetTopMargin(0.12)
    c.SetBottomMargin(0.14)
    my_detector.sum_cu.GetXaxis().SetTitle("Time [s]")
    my_detector.sum_cu.GetYaxis().SetTitle("Current [A]")
    my_detector.sum_cu.Draw("HIST")
    my_detector.positive_cu.Draw("SAME HIST")
    my_detector.negtive_cu.Draw("SAME HIST")
    c.Update()
    rightmax = 1.1*ele_current.GetMinimum()
    n_scale = ROOT.gPad.GetUymin() / rightmax
    ele_current.Scale(n_scale)
    ele_current.Draw("SAME HIST")
    my_detector.sum_cu.SetLineColor(3)
    my_detector.positive_cu.SetLineColor(2)
    my_detector.negtive_cu.SetLineColor(4)
    ele_current.SetLineColor(6)
    
    axis = ROOT.TGaxis(ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymin(), ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymax(), rightmax, 0, 510, "+L")
    axis.SetLineColor(6)
    axis.SetTextColor(6)
    axis.SetTextSize(0.02)
    axis.SetLabelColor(6)
    axis.SetLabelSize(0.02)
    axis.SetTitle("Ampl [mV]")
    axis.CenterTitle()
    axis.Draw("same")

    legend = ROOT.TLegend(0.5, 0.3, 0.9, 0.6)
    legend.AddEntry(my_detector.sum_cu, "sum", "l")
    legend.AddEntry(my_detector.negtive_cu, "electron", "l")
    legend.AddEntry(my_detector.positive_cu, "hole", "l")
    legend.AddEntry(ele_current, "current after electric", "l")
    legend.SetBorderSize(0)
    legend.SetTextFont(43)
    legend.SetTextSize(30)
    legend.Draw("same")

    c.Update()
    c.SaveAs("fig/basic_infor.root")
    del c

    if drift_path == 1:
        my_drift.draw_drift_path(my_detector,my_field,[100,100,100])

def draw_ele_field(my_d,my_f,plane,depth):

    ROOT.gStyle.SetOptStat(0)
    c1 = ROOT.TCanvas("c", "canvas",1000, 1000)
    c1.SetLeftMargin(0.12)
    c1.SetRightMargin(0.2)
    c1.SetBottomMargin(0.14)
    c1.Divide(2,2)
    c1.GetPad(1).SetRightMargin(0.2)
    c1.GetPad(2).SetRightMargin(0.2)
    c1.GetPad(3).SetRightMargin(0.2)
    c1.cd(1)
    e_field=fill_his("E",depth,my_d,my_f,plane)
    e_field.Draw("COLZ")
    
    c1.Update()
    c1.cd(2)
    p_field=fill_his("P",depth,my_d,my_f,plane)
    p_field.Draw("COLZ")
    c1.SetRightMargin(0.12)
    c1.Update()
    c1.cd(3)
    w_p_field=fill_his("WP",depth,my_d,my_f,plane)
    w_p_field.Draw("COLZ")
    c1.SetRightMargin(0.12)
    c1.Update()
    c1.SaveAs("fig/ele_field"+plane+str(depth)+".root")
    del c1

def fill_his(model,depth,my_d,my_f,plane):

    if plane == "xy":
        l_x = my_d.l_x
        l_y = my_d.l_y
        t_name = plane + " at z = " + str(depth)
    elif plane == "yz":
        l_x = my_d.l_y
        l_y = my_d.l_z
        t_name = plane + " at x = " + str(depth)
    elif plane == "xz":
        l_x = my_d.l_x
        l_y = my_d.l_z
        t_name = plane + " at y = " + str(depth)
    else:
        print("the draw plane is not existing")
    nx_e = 500
    ny_e = 500
    e_v = ROOT.TH2F("","",nx_e,0,l_x,ny_e,0,l_y)
    if model == "E":
        v_sample = my_f.E_field
        e_v.SetTitle("electric field "+t_name)
    elif model == "P":
        v_sample = my_f.u
        e_v.SetTitle("potential "+t_name)
    elif model == "WP":
        v_sample = my_f.u_w  
        e_v.SetTitle("weigthing potential "+t_name)      

    for j in range (ny_e):
        for i in range(nx_e):
            x_v = (i+1)*(l_x/nx_e)
            y_v = (j+1)*(l_y/ny_e)
            try:
                if plane == "xy":
                    f_v = v_sample(x_v,y_v,depth)
                elif plane == "yz":
                    f_v = v_sample(depth,x_v,y_v)
                elif plane == "xz":
                    f_v = v_sample(x_v,depth,y_v)
                if model == "E":
                    f_v = math.sqrt(math.pow(f_v[0],2)+math.pow(f_v[1],2)+math.pow(f_v[2],2))              
            except RuntimeError:
                f_v = 0.0
            e_v.SetBinContent(i+1,j+1,f_v)
    if plane == "xy":
        e_v.GetXaxis().SetTitle("x")
        e_v.GetYaxis().SetTitle("y")
    elif plane == "yz":
        e_v.GetXaxis().SetTitle("y")
        e_v.GetYaxis().SetTitle("z")
    elif plane == "xz":
        e_v.GetXaxis().SetTitle("x")
        e_v.GetYaxis().SetTitle("z") 
    return e_v
def save_charge(charge_t,qtot,x_v,y_v):
    now = int(time.time())     # 1533952277
    timeArray = time.localtime(now)
    otherStyleTime = time.strftime("%Y--%m--%d %H:%M:%S", timeArray)
    with open("test.txt",'w') as f:
        f.write(otherStyleTime)
        f.write(str(charge_t)+','+str(qtot)+'\n')
def threeD_time():
    ### define the structure of the detector
    my_detector = R3dDetector(300.0,300.0,300.0)
    my_detector.set_para(doping=-10.0,voltage=-500.0,temperature=300.0)
    my_detector.set_electrode()
    ### get the electric field and weighting potential
    my_field = Fenics_cal(my_detector,mesh_v=55)
    my_field.fenics_p_electric(my_detector) 
    my_field.fenics_p_w_electric(my_detector)
    # ### define the tracks and type of incident particles
    my_track = Tracks()
    my_track.t_mip([my_detector.e_t_4[0],my_detector.e_t_4[1],0.0/FACTOR],[my_detector.e_t_6[0],my_detector.e_t_6[1],300.0/FACTOR],300)
    #my_track.t_mip([115,150,0],[115,150,300],300)
    # # ### drift of ionized particles
    my_drift = Drifts(my_track)
    my_drift.ionized_drift(my_track,my_field,my_detector)
    # # ### after the electronics
    my_electronics = Amplifier()
    charge_t,qtot,ele_current=my_electronics.CSA_amp(my_detector,t_rise=0.4,t_fall=0.2,trans_imp=10)
    # # ### electric plot
    draw_ele_field(my_detector,my_field,"xy",my_detector.l_z*0.5)
    # # ###  current plot
    draw_plot(my_detector,ele_current,qtot,my_drift,my_field,drift_path=1)

### get the 2D time resolution
def threeD_time_scan(output,numbers,t_numbers,n_step):
    ### define the structure of the detector
    my_detector = R3dDetector(300,300,300)
    my_detector.set_para(doping=-10,voltage=-2000,temperature=300)
    my_detector.set_electrode()
    ### get the electric field and weighting potential
    my_field = Fenics_cal(my_detector,mesh_v=100)
    my_field.fenics_p_electric(my_detector) 
    my_field.fenics_p_w_electric(my_detector)
    # ### define the tracks and type of incident particles
    my_track = Tracks()
    x_min = my_detector.e_t_2[0]-my_detector.e_t_2[2]
    x_max = my_detector.e_t_3[0]+my_detector.e_t_3[2]
    y_min = my_detector.e_t_7[1]-my_detector.e_t_7[2]
    y_max = my_detector.e_t_3[1]+my_detector.e_t_3[2]
    nx = ny = int(math.sqrt(t_numbers))
    x_step = (x_max-x_min)/(nx-1)
    y_step = (y_max-y_min)/(ny-1)
    for i in range (numbers-n_step,numbers):
        n_y = int(i/nx)
        n_x = i%nx
        x_v = x_min + x_step*n_x
        y_v = y_min + y_step*n_y
        my_track.t_mip([x_v,y_v,0],[x_v,y_v,300],300)
    # ### drift of ionized particles
        my_drift = Drifts(my_track)
        my_drift.ionized_drift(my_track,my_field,my_detector)
    #     ### after the electronics
        my_electronics = Amplifier()
        charge_t,qtot,ele_current=my_electronics.CSA_amp(my_detector,t_rise=0.4,t_fall=0.2,trans_imp=10)
        ROOT.gROOT.SetBatch(1)
        c = ROOT.TCanvas("Plots", "Plots", 1000, 1000)
        ele_current.Draw("HIST")
        c.Update()
        output_path = output
        os.system("mkdir %s -p"%(output_path))
        c.SaveAs(output+"/t_"+str(i)+"_events.C")
        del c
        save_charge(charge_t,qtot,x_v,y_v,)

def main():
    args = sys.argv[1:]
    model = args[0]
    if model == "3D":
        threeD_time()
    if model == "3D_scan":
        output = args[1]
        numbers = int(args[2])
        t_numbers = int(args[3])
        n_step = int(args[4])
        threeD_time_scan(output,numbers,t_numbers,n_step)
    print("run end")
    
if __name__ == '__main__':
    #record run time
    starttime = time.time()
    main() 
    endtime=time.time()
    dtime = endtime-starttime
    print ("the process run time %.8s s" %dtime) 