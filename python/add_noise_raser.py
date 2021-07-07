'''
author: tanyuhang
time: 2021.3.8
Use: Read the data of KDetsim induced current
'''
import os
import sys
import re
from ROOT import TCanvas,TGraph,TMultiGraph,TLegend,gRandom,TF1,TLatex
# from tools import set_root_style
from array import array
import ROOT
from ROOT import gStyle
import math

Events=array('i',[0])
h_pulse_time=ROOT.std.vector(float)()
h_pulse_height=ROOT.std.vector(float)()
h_noise_height=ROOT.std.vector(float)()
h_bc_noise_height=ROOT.std.vector(float)()
h_bc_max_pulse_height=array('f',[0])
h_max_pulse_height=array('f',[0])
h_max_pulse_time=array('f',[0])
h_time_resolution=array('f',[0])
h_noise_height_jitter=array('f',[0])

h_pulse_time.clear()
h_pulse_height.clear()
h_noise_height.clear()
h_bc_noise_height.clear()
h_bc_max_pulse_height[0]=0.0
h_max_pulse_height[0]=0.0
h_max_pulse_time[0]=0.0
h_time_resolution[0]=0.0
h_noise_height_jitter[0]=0.0

def main(): 
    args = sys.argv[1:]
    input_file=args[0]
    out_p=args[1]
    thre_vth=10 #mv
    CFD=0.5
    list = os.listdir(input_file) 
    mg = TMultiGraph("mg","")
    leg=TLegend(0.55,0.65,0.75,0.85)
    color = []
    marker = []

    out_file=ROOT.TFile(out_p,"RECREATE")
    tree_out=ROOT.TTree('tree','tree')

    tree_out.Branch('Events',Events,'Events/I')
    tree_out.Branch('h_pulse_time',h_pulse_time)
    tree_out.Branch('h_pulse_height',h_pulse_height)
    tree_out.Branch('h_noise_height',h_noise_height)
    tree_out.Branch('h_bc_noise_height',h_bc_noise_height)
    tree_out.Branch('h_bc_max_pulse_height',h_bc_max_pulse_height,'h_bc_max_pulse_height/F')    
    tree_out.Branch('h_max_pulse_height',h_max_pulse_height,'h_max_pulse_height/F')    
    tree_out.Branch('h_max_pulse_time',h_max_pulse_time,'h_max_pulse_time/F')
    tree_out.Branch('h_time_resolution',h_time_resolution,'h_time_resolution/F')
    tree_out.Branch('h_noise_height_jitter',h_noise_height_jitter,'h_noise_height_jitter/F')
    #start_write_to_txt(out_path) 
    c_number=0
    marker,color=defind_color_marker(marker,color)
    i=1
    jk=0
    CFD_time=[]
    file_n = []
    file_m = []
    bc_CFD_times = []
    for root,dirs,files in os.walk(input_file):
        for file in files:
            if i<10000:
                print("................Events:%s..............."%(Events[0])) 
                Events[0]+=1
                i+=1
                list_c=[]
                path = os.path.join(input_file, file)
                out_path=out_p+file #txt
                for line in open(path):
                    list_c.append(line)
                time_list=[]
                ampl_CSA_list=[]
                max_pulse_height=0.0
                max_pulse_time=0.0
                noise_height=[]
                time_r=0.0
                noise_height_jitter=0.0
                if len(list_c)>20:
                    bc_CFD_time=toa_bc(list_c,CFD)
                    bc_CFD_times.append(bc_CFD_time)
                    file_n.append(file)
                    time_list,ampl_CSA_list,max_pulse_height,max_pulse_time,noise_height=add_noise(list_c)
                    for j in range(0,len(time_list)):
                        h_bc_noise_height.push_back(noise_height[j])
                    h_bc_max_pulse_height[0]=max_pulse_height
                    if (max_pulse_height>thre_vth and  max_pulse_height<50):
                        file_m.append(file)
                        jk=jk+1
                        h_max_pulse_height[0]=max_pulse_height
                        h_max_pulse_time[0]=max_pulse_time
                        
                        time_r,noise_height_jitter=get_CFD_time(time_list,ampl_CSA_list,max_pulse_height,max_pulse_time,CFD)
                        CFD_time.append(time_r)                     
                        h_time_resolution[0]=time_r
                        h_noise_height_jitter[0]=noise_height_jitter
                        for j in range(0,len(time_list)):
                            h_pulse_time.push_back(time_list[j])
                            h_pulse_height.push_back(ampl_CSA_list[j])
                            h_noise_height.push_back(noise_height[j])
                        tree_out.Fill()
                        h_max_pulse_height[0]=0.0
                        h_bc_max_pulse_height[0]=0.0
                        h_max_pulse_time[0]=0.0  
                        h_time_resolution[0]=0.0
                        h_noise_height_jitter[0]=0.0
                        h_bc_noise_height.clear()
                        h_pulse_time.clear()
                        h_pulse_height.clear()	
                        h_noise_height.clear()
                #write_end_to_txt(out_path,time_list)
                    #gr=save_graph(time_list,ampl_CSA_list)
        #         #gr=set_color_marker(color,marker,c_number,gr)              
        #         #leg=fill_legend(leg,gr,out_path)   
        #         c_number=c_number+1  
                    # mg.Add(gr)
                    # mg1.Add(gr1)
                    # mg2.Add(gr2)
            else:
                break
    print(jk)
    tree_out.Write()
    out_file.Close()
    draw_toa(bc_CFD_times,file_n,input_file,out_p)
    draw_CFD_time(CFD_time,input_file,out_p,file_m)
            # i+=1
    # draw_mg(mg,leg,out_p,"BB")
    # draw_mg(mg1,leg,out_p,"CSA")
    # draw_mg(mg2,leg,out_p,"Itotal")
def draw_toa(CFD_time,input_n,input_file,out_p):
    output_file = input_file.split("/")[2]
    out_name = out_p.split("/")[1]
    c1 = TCanvas("c1", "canvas1", 800, 600)
    #gStyle.SetPalette(57)
    gStyle.SetOptStat(0)
    ROOT.gROOT.SetBatch(1)
    #set_root_style(stat=0, grid=0)
    c1.SetGrid()
    c1.SetLeftMargin(0.12)
    # c.SetTopMargin(0.12)
    c1.SetBottomMargin(0.14)
    histo1=ROOT.TH2F("","",100,0,100,100,0,100)
    histo1.GetXaxis().SetTitleOffset(1.2)
    histo1.GetXaxis().SetTitleSize(0.05)
    histo1.GetXaxis().SetLabelSize(0.05)
    histo1.GetXaxis().SetNdivisions(510)
    histo1.GetYaxis().SetTitleOffset(1.1)
    histo1.GetYaxis().SetTitleSize(0.05)
    histo1.GetYaxis().SetLabelSize(0.05)
    histo1.GetYaxis().SetNdivisions(505)
    histo1.GetXaxis().SetTitle("x [#mum]");
    histo1.GetYaxis().SetTitle("y [#mum]")
    histo1.GetXaxis().CenterTitle()
    histo1.GetYaxis().CenterTitle()
    gap_v=float(output_file.split('_')[1])
    voltage_v=float(output_file.split('_')[3])
    # e_r=5.0
    # l_x = 100
    # l_y = 100
    # e_t_y = math.sqrt(gap_v*gap_v*0.75)
    # x_min = l_x*0.5 - gap_v*0.5
    # x_max = l_x*0.5 + gap_v*0.5
    # y_min = l_y*0.5 - e_t_y
    # y_max = l_y*0.5 + e_t_y
    for i in range(len(CFD_time)):
        n_x = float(input_n[i].split("_")[3])
        n_y = float(input_n[i].split("_")[5])
        r_v = math.sqrt((n_x-50)*(n_x-50) + (n_y-50)*(n_y-50))
        # if CFD_time[i]*1e9>=0 and CFD_time[i]*1e9<=2:
        #     time = CFD_time[i]*1e9
        # else:
        #     time = 0.0*1e9
        if CFD_time[i]!=0 and r_v <= gap_v: 
            histo1.Fill(n_x,n_y,CFD_time[i]*1e9)

    ROOT.gStyle.SetPalette(107)
    histo1.Draw("COLZ")
    # histo1.SetMaximum(1.5)
    # histo1.SetMinimum(0)
    c1.Update()
    c1.SaveAs("result/"+out_name+"/toa_"+output_file+".pdf")
    c1.SaveAs("result/"+out_name+"/toa_"+output_file+".C")
def toa_bc(list_c,CFD):
    bc_time_list = []
    bc_pulse_height = []
    for j in range (0,len(list_c)):
        time= float(list(filter(None,list_c[j].split(",")))[0])
        ampl_CSA=float(list(filter(None,list_c[j].split(",")))[1])
        bc_time_list.append(time)
        bc_pulse_height.append(ampl_CSA)
    max_pulse_height=max(bc_pulse_height)
    max_index=bc_pulse_height.index(max(bc_pulse_height))
    max_pulse_time=bc_time_list[max_index]
    bc_CFD_time = 0
    for i in range (0,len(bc_time_list)):
        if bc_pulse_height[i]>=max_pulse_height*CFD and bc_time_list[i]<max_pulse_time:
            bc_CFD_time = bc_time_list[i]
            break
    return bc_CFD_time

def add_noise(list_c):
    time_list=[]
    ampl_CSA_list=[]
    noise_height_list=[]
    time=0.0
    ampl_CSA=0.0
    #noise
    gRandom.SetSeed(0)
    random_gauss = gRandom.Gaus
    for j in range (0,len(list_c)):
        time= float(list(filter(None,list_c[j].split(",")))[0])
        noise_height=random_gauss(0.293,2.606)
        ampl_CSA=float(list(filter(None,list_c[j].split(",")))[1])+noise_height
        noise_height_list.append(noise_height)
        time_list.append(time)
        ampl_CSA_list.append(ampl_CSA)
    max_pulse_height=max(ampl_CSA_list)
    max_index=ampl_CSA_list.index(max(ampl_CSA_list))
    max_pulse_time=time_list[max_index]
    return time_list,ampl_CSA_list,max_pulse_height,max_pulse_time,noise_height_list

def draw_CFD_time(CFD_time,out_put,out_p,file_m):
    out_put_1=out_put.split("/")[2]
    voltage_1 = out_put.split("/")[2]
    voltage_2 = float(voltage_1.split("_")[3])
    gap_v = float(voltage_1.split("_")[1])
    out_name = out_p.split("/")[1]
    c = TCanvas("c", "canvas", 800, 600)
    #gStyle.SetPalette(57)
    # e_r=5.0
    # l_x = 100
    # l_y = 100
    # e_t_y = math.sqrt(gap_v*gap_v*0.75)
    # x_min = l_x*0.5 - gap_v*0.5
    # x_max = l_x*0.5 + gap_v*0.5
    # y_min = l_y*0.5 - e_t_y
    # y_max = l_y*0.5 + e_t_y
    gStyle.SetOptStat(0)
    ROOT.gROOT.SetBatch(1)
    #set_root_style(stat=0, grid=0)
    c.SetGrid()
    c.SetLeftMargin(0.2)
    # c.SetTopMargin(0.12)
    c.SetBottomMargin(0.2)
    leg = ROOT.TLegend(0.25, 0.6, 0.40, 0.8)
    
    #FormatLegend(leg)

    histo=ROOT.TH1F("","",100,-0.2,0.6)
    gStyle.SetOptFit()
    #gStyle.SetStatY(0.6);
    for i in range(0,len(CFD_time)):
        n_x = float(file_m[i].split("_")[3])
        n_y = float(file_m[i].split("_")[5])
        r_v = math.sqrt((n_x-50)*(n_x-50) + (n_y-50)*(n_y-50))
        if CFD_time[i]!=0 and r_v <= gap_v:
            histo.Fill(CFD_time[i])
    histo.GetXaxis().SetTitle("ToA [ns]");
    histo.GetYaxis().SetTitle("Events")
    histo.GetXaxis().CenterTitle()
    histo.GetYaxis().CenterTitle()
    fit_func_1 = TF1('fit_func_1','gaus',-0.2,0.5)
    histo.Fit("fit_func_1","ROQ+","",-0.2,0.5)
    print("constant:%s"%fit_func_1.GetParameter(0))
    print("constant_error:%s"%fit_func_1.GetParError(0))
    print("mean:%s"%fit_func_1.GetParameter(1))
    print("mean_error:%s"%fit_func_1.GetParError(1))
    print("sigma:%s"%fit_func_1.GetParameter(2))
    print("sigma_error:%s"%fit_func_1.GetParError(2))
    sigma=fit_func_1.GetParameter(2)*1000
    error=fit_func_1.GetParError(2)*1000
    leg.AddEntry(fit_func_1,"Fit","L")
    leg.AddEntry(histo,"Sim","L")
    histo.GetXaxis().SetTitleOffset(1.2)
    histo.GetXaxis().SetTitleSize(0.07)
    histo.GetXaxis().SetLabelSize(0.05)
    histo.GetXaxis().SetNdivisions(510)
    histo.GetYaxis().SetTitleOffset(1.1)
    histo.GetYaxis().SetTitleSize(0.07)
    histo.GetYaxis().SetLabelSize(0.05)
    histo.GetYaxis().SetNdivisions(505)
    histo.GetXaxis().CenterTitle()
    histo.GetYaxis().CenterTitle()
    histo.SetLineWidth(2)
    fit_func_1.SetLineWidth(2)
    histo.Draw()
    fit_func_1.Draw("same")
    leg.Draw("same")
    tex = ROOT.TLatex()
    tex.SetNDC(1)
    tex.SetTextFont(43)
    tex.SetTextSize(25)
    tex.DrawLatexNDC(0.6, 0.8, "Electron Spacing:%.0f#mum"%(gap_v))
    tex.DrawLatexNDC(0.6, 0.7, "Volatge=%.0fV,CFD=0.5"%(voltage_2))
    tex.DrawLatexNDC(0.6, 0.6, "#sigma = %.0f #pm %.0f ps"%(sigma,error))
    c.Update()
    c.SaveAs("result/"+out_name+"/"+out_put_1+".pdf")
    c.SaveAs("result/"+out_name+"/"+out_put_1+".c")
    c.SaveAs("result/"+out_name+"/"+out_put_1+".root")

def get_CFD_time(time_list,ampli_list,max_pulse_height,max_pulse_time,CFD):

    random_gauss = gRandom.Gaus
    noise_height_jitter=random_gauss(7.82,1.01)
    CFD_time=0.0
    #print ("noise_height_1:%s" %noise_height)
    # noise_height=random_gauss(0.66,2.36)
    # print ("noise_height_2:%s" %noise_height)
    for i in range (0,len(time_list)):
        if ampli_list[i]>=max_pulse_height*CFD and time_list[i]<max_pulse_time:
            dVdt=(ampli_list[i+1]-ampli_list[i-1])/(time_list[i+1]-time_list[i-1])
            jit=0
            if (dVdt!=0):
                jitter=noise_height_jitter/dVdt*1e9
                CFD_time = time_list[i]*1e9+random_gauss(0,jitter)
            else:
            #CFD_time=time_list[i]
                CFD_time=time_list[i]*1e9
            break
    return CFD_time,noise_height_jitter


def FormatLegend(leg):
    
    leg.SetBorderSize(1)
    leg.SetTextFont(43)
    leg.SetTextSize(40)
    leg.SetFillStyle(1)
    leg.SetFillColor(1)
    leg.SetLineColor(2) 

def set_color_marker(color,marker,i,gr):
    f=marker[i]
    gr.SetMarkerStyle(f)
    gr.SetMarkerSize(1)
    k=color[i]
    gr.SetLineColor(k)
    gr.SetLineWidth(2)
    gr.SetMarkerColor(k)
    return gr

def fill_legend(leg,gr,name):
    # leg name define
    #leg_name="position:x"+str(x_number)+"_y"+str(y_number)
    # print c_number
    # if c_number == 0:
    #     leg_name="first_point"
    # if c_number == 1:
    #     leg_name="second_point"
    # if c_number == 2:
    #     leg_name="third_point"
    leg.AddEntry(gr,name,"LP")
    return leg

def defind_color_marker(marker,color):
    
    color.append(int(2))
    color.append(int(45))
    color.append(int(3))
    color.append(int(7))
    color.append(int(46))
    color.append(int(38))
    color.append(int(40))
    color.append(int(4))
    color.append(int(8))
    color.append(int(11))

    marker.append(int(20))
    marker.append(int(24))
    marker.append(int(23))
    marker.append(int(25))
    marker.append(int(33)) 
    marker.append(int(26))
    marker.append(int(30))
    marker.append(int(27))
    marker.append(int(28))
    marker.append(int(13))
    return marker,color
def FormatLegend(leg):
    
    leg.SetBorderSize(0)
    leg.SetTextFont(43)
    leg.SetTextSize(30)
    leg.SetFillStyle(0)
    leg.SetFillColor(0)
    leg.SetLineColor(0)    
if __name__ == '__main__':
    main()   