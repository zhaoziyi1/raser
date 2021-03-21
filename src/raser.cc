// 
// Author: SHI Xin <shixin@ihep.ac.cn> 
// Created [2021-03-31 Sun 15:05] 
// Based on KDetSim : http://kdetsim.org 
// 


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <TH2D.h> 
#include <TFile.h> 
#include <TCanvas.h> 
#include <TStyle.h> 
#include <Riostream.h>
#include <TTree.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TApplication.h> 
#include <TMath.h> 

// using namespace std; 

void set_root_style(int stat=1110, int grid=0){
  gROOT->Reset();

  gStyle->SetTitleFillColor(0) ; 
  gStyle->SetTitleBorderSize(0); 
    
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetCanvasDefX(0); 
  gStyle->SetCanvasDefY(0); 
  gStyle->SetFrameBorderMode(0); 
  gStyle->SetFrameBorderSize(1); 
  gStyle->SetFrameFillColor(0); 
  gStyle->SetFrameFillStyle(0); 
  gStyle->SetFrameLineColor(1); 
  gStyle->SetFrameLineStyle(1); 
  gStyle->SetFrameLineWidth(1); 

  // gStyle->SetPadTopMargin(PadTopMargin);  
  gStyle->SetPadLeftMargin(0.10);  
  gStyle->SetPadRightMargin(0.05);  

  gStyle->SetLabelSize(0.03, "XYZ");  
  gStyle->SetTitleSize(0.04, "XYZ");  
  gStyle->SetTitleOffset(1.2, "Y");  

  gStyle->SetPadBorderMode(0);  
  gStyle->SetPadColor(0);  
  gStyle->SetPadTickX(1); 
  gStyle->SetPadTickY(1); 
  gStyle->SetPadGridX(grid); 
  gStyle->SetPadGridY(grid); 

  gStyle->SetOptStat(stat); 
  gStyle->SetStatColor(0); 
  gStyle->SetStatBorderSize(1); 
}


TGraph * get_graph_from_log(TString inputFile, TString& err_msg) {
  std::ifstream in;
  in.open(inputFile); 
  Float_t x, y; 
  Float_t factor_I(1.0);
  if (inputFile.Contains("_uA_") ) {
    factor_I = 1e-6;
    std::cout << "Using current factor: " << factor_I << " for " << inputFile << std::endl;  
  } 
  std::vector<Float_t> voltages; 
  std::vector<Float_t> currents; 

  std::string line;
  int nlines = 0;

  // float i_v150, i_v100;
  double i_v150, i_v100;
  bool pass1(false), pass2(false); 

  while (getline(in, line)) {
    std::istringstream iss(line);
    if ( line.find("#") == 0 ) continue; 
    if (!(iss >> x >> y )) break; 
    if (!in.good()) break;
    voltages.push_back(fabs(x));
    currents.push_back(fabs(y)*factor_I);
    // Pick up values like: -100.043	
    if ( fabs(fabs(x)-150) < 1) i_v150 = fabs(y); 
    if ( fabs(fabs(x)-100) < 1) i_v100 = fabs(y); 
    nlines ++; 
  }

  if (nlines < 1) {
    std::cerr << "No valid data found in : " << inputFile << std::endl;
    return NULL; 
  }
  
  if ( i_v150 < 2E-6) pass1 = true;
  if ( i_v150/i_v100 < 2 ) pass2 = true;

  if (!pass1) err_msg = Form("I(150V) >= 2uA (%.1e)",i_v150) ;
  if (!pass2) err_msg += Form("I(150V)/I(100V) >= 2 (%.1f)", i_v150/i_v100) ;
  
  in.close();
  TGraph *gr = new TGraph(nlines, &voltages[0], &currents[0]);
  return gr; 
}

TCanvas* drawIV(std::vector<TString> inputFiles){
  set_root_style();

  TCanvas *c = new TCanvas("c", "IV scan", 800, 800);
  c->SetGrid();

  TMultiGraph *mg = new TMultiGraph();
  TLegend *leg = new TLegend(0.2, 0.6, 0.5, 0.8);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetNColumns(1);
  leg->SetTextSize(0.02);
  leg->SetTextSizePixels(25);

  for (std::vector<int>:: size_type i = 0; i != inputFiles.size(); i++) {
    TString err_msg = "";  
    TGraph *gr = get_graph_from_log(inputFiles[i], err_msg);
    if (gr == NULL) continue; 
    gr->SetMarkerStyle(20+i);
    gr->SetMarkerSize(0.9);
    int color = i+1;
    if (color >= 5) color ++; // bypass the yellow  
    if (color >= 10) color = color % 10 + 1 ; // reuse the first 9 colors
    gr->SetMarkerColor(color);
    leg->AddEntry(gr, Form("%s %s", inputFiles[i].Data(),
			   err_msg.Data()), "p"); 
    mg->Add(gr); 
  }
  
  mg->Draw("APL"); 
  mg->GetXaxis()->SetTitle("Bias Voltage [V]");
  mg->GetYaxis()->SetTitle("Leakage Current [A]");
  mg->GetYaxis()->SetRangeUser(1e-10, 1e-4); 
  leg->Draw(); 

  c->SetLogy();
  c->Update(); 
  return c;
}


#ifndef __CINT__ 

void print_usage(){
  printf("NAME\n\tdrawIV - draw IV curves from log files\n");
  printf("\nSYNOPSIS\n\tdrawIV [-b] [-h] file1 file2 ...\n");
  printf("\nOPTIONS\n");
  printf("\t%-5s  %-40s\n", "-h", "Print this message");
  printf("\t%-5s  %-40s\n", "-b", "Batch mode, save to pdf file directly");
  printf("\nAUTHOR\n\tXin Shi <Xin.Shi@cern.ch>\n");
}

int main(int argc, char** argv) {
  
  if (argc < 2) {
    print_usage() ;  
    return -1; 
  }

//    K3D *det = new K3D(7, 80, 80, 300);

  bool doBatch(false);
  std::vector<TString> inputFiles(argv+1, argv+argc);
  TString outFile = "test.pdf";
  
  for (int i = 0; i < argc; i++){
    if (!strcmp(argv[i], "-h")) {
      print_usage();
      break; 
    }

    if (!strcmp(argv[i], "-b")) {
      doBatch = true;
      inputFiles.erase(inputFiles.begin()+i-1);
    }
  }

  if (doBatch) { 
    std::cout << "Run in batch mode ... " << std::endl;
    TCanvas *c = drawIV(inputFiles);
    c->SaveAs(outFile);
    delete c;
    gSystem->Exit(0);
  }
  
  TApplication theApp("App", 0, 0);
  theApp.SetReturnFromRun(true);
  drawIV(inputFiles);
  theApp.Run();
}

#endif

