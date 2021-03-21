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


// KDetSim 


// KGeometry 

#include "TH3I.h"
#include "TH2F.h"

// #include "TMath.h"
#include "math.h"

class KGeometry
{
private:
public:

  TH3I    *EG;        //electrode geometry
  TH3I    *DM;        //detector material
  Int_t nx;           //x-divisions
  Int_t ny;           //y-divisions
  Int_t nz;           //z-divisions 
  
  // Constructors of the class
  KGeometry();
  ~KGeometry();
  KGeometry(TH3I *x){GetGrid(x,0); };
  KGeometry(TH3I *x, TH3I *y){GetGrid(x,0); GetGrid(x,1);};
  void GetGrid(TH3I *,Short_t =0); 
//   ClassDef(KGeometry,1) 
};



// ClassImp(KGeometry)
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// KGeometry                                                            //
//                                                                      //
// Class for description of detector geometry                           //
// The class defines the geometry of the detector.                      //
// It is based upon two TH3S histograms                                 // 
// It contains some functions defined to design the electrodes          //  
//                                                                      //
//////////////////////////////////////////////////////////////////////////

KGeometry::KGeometry()
{
  EG=NULL;
  DM=NULL;
  nx=1;
  ny=1;
  nz=1;
}

KGeometry::~KGeometry()
{
  if(EG!=NULL) delete EG;
  if(DM!=NULL) delete DM;
}

void KGeometry::GetGrid(TH3I *x, Short_t which)
{
  // EG grid
  //bit 1 = 1  -> 1st electrode
  //bit 2 = 2  -> 2nd electrode
 
  switch(which)
    {
    case 0: 
      if(EG!=NULL) delete EG;
	EG=new TH3I(); x->Copy(*EG);
	nx=EG->GetNbinsX();
	ny=EG->GetNbinsY();
	nz=EG->GetNbinsZ();
      break;
    case 1:
      if(DM!=NULL) delete DM;
	DM=new TH3I(); x->Copy(*DM); 
      break;
    }

   if(DM!=NULL) if(DM->GetNbinsX()!=nx) printf("Warning: dimenssions mismatch - X !\n");
   if(DM!=NULL) if(DM->GetNbinsY()!=ny) printf("Warning: dimenssions mismatch - Y !\n");
   if(DM!=NULL) if(DM->GetNbinsZ()!=nz) printf("Warning: dimenssions mismatch - Z !\n");

}

// KMaterial 


class KMaterial
{
private:
public:
    static Int_t Mat;              // Material index
    static Float_t Temperature;    // Temperature
    static Int_t Mobility;         // mobility model for each material
    static Int_t ImpactIonization; // impact ionization model

    //////////////////////////////////////////////////////

    KMaterial() { Mat = 1; } // MobMod=1;}
    ~KMaterial(){};
    // ClassDef(KMaterial, 1)
};

// ClassImp(KMaterial)
Int_t KMaterial::Mat = 1;
Float_t KMaterial::Temperature = 293;

// KField 

class KField
{
private:
  Int_t Method;   // Method to calculate the intermediate points
  Int_t dim;
public:
  TH3F *U;
  TH3F *Ex;
  TH3F *Ey;
  TH3F *Ez;
  TH3F *E;

  KField() {U=NULL; Ex=NULL; Ey=NULL; Ez=NULL;};
 ~KField();
//   ClassDef(KField,1) 
   };

KField::~KField()
{
  if(U!=NULL)  delete U; 
  if(Ex!=NULL) delete Ex; 
  if(Ey!=NULL) delete Ey; 
  if(Ez!=NULL) delete Ez;
}



// KDetector 

#include "TRandom.h"
#include "TF3.h"

class KDetector : public KGeometry, public KMaterial { 


private:
  Double_t Deps;
  TRandom *ran;               //random number generator
  Double_t CalErr;               //Error of the solver
  Int_t MaxIter;              //Maximum number of iterations in eq solver
  Short_t Debug;              //Print information of drift calculation etc.

public:
  Float_t Voltage;  //Voltage
  Float_t Voltage2; //Voltage2 
  TArrayF Voltages; //Array of voltages

  // Definition of space charge 
  TF3     *NeffF;     //effective dopping concentration function
  TH3F    *NeffH;     //effective dopping concentration histogram

  // Weigthing, electric and magnetic field
  KField *Ramo;       // ramo field 
  KField *Real;       // electric field
  Float_t B[3];      // magnetic field

  // Trapping and variables used for multiplication studies
  Float_t taue;      // effective trapping time constants 
  Float_t tauh;      // effective trapping time constants 
  TF3 *TauE;         // Function of TauE(x,y,z);
  TF3 *TauH;         // Function of TauH(x,y,z);

  Int_t BreakDown;     // if break down occurs it goes to 1 otherwise is 0
  Float_t MTresh;      // treshold for taking multiplication into account
  Float_t BDTresh;     // hole multiplication - break down treshold 
  Float_t DiffOffField;// electric field where diffusion is switched off [V/um] !!
  // Drift parameters

  Float_t enp[3];      //entry point for the charge drift
  Float_t exp[3];      //exit point for the cahrge drift
  Int_t diff;          // Diffusion simulation (yes=1, no=0)
  Int_t average;       // Average (over how many events)
  Float_t SStep;       // Simulation step size;
  Float_t MaxDriftLen; // Maximum drift lenght before stopping the drift

  // Output histograms
  TH1F *pos;           // contribution of the holes to the total drift current
  TH1F *neg;           // contribution of the electrons  to the total drift current
  TH1F *sum;	       // total drift current

  // Constructors and destructor
  KDetector();
  ~KDetector();
  //_______________________________________________________________________________
  void SetDriftHisto(Float_t x,Int_t=200);

  //Configuration functions 
    // ClassDef(KDetector,1) 
};

#include "TPolyLine3D.h"
// #include "KDetector.h"
#include "TFile.h"

#define ABS(x) x>0?x:-x
#define PREDZNAK(x) x>0?1:-1

#define C1(x,y,z) y3[n]=x+y+z;

#define L1(x) y2[n]=x;
#define R1(x) y4[n]=x;
#define U1(x) y5[n]=x;
#define D1(x) y6[n]=x;
#define I1(x) y7[n]=x;
#define O1(x) y8[n]=x;

#define L2(x) y2[n]=2*x;
#define R2(x) y4[n]=2*x;
#define U2(x) y5[n]=2*x;
#define D2(x) y6[n]=2*x;
#define I2(x) y7[n]=2*x;
#define O2(x) y8[n]=2*x;

#define C0 y3[n]=1.;
#define U0 y5[n]=0.;
#define D0 y6[n]=0.;
#define R0 y4[n]=0.;
#define L0 y2[n]=0.;
#define I0 y7[n]=0.;
#define O0 y8[n]=0.;

#define  PI  3.1415927
#define EPS 1.0e-14
#define STEP_DET_CUR 25e-9

// ClassImp(KDetector)


double **a,*b,*y6,*y2,*y3,*y4,*y5,*y7,*y8;


void KDetector::SetDriftHisto(Float_t x, Int_t numbins)
{
  if(x<0 || x>10000e-9) 
    printf("Selected time range out of scope !\n"); 
  else
    {
    if(pos!=NULL) delete pos;
	pos  = new TH1F("charge+","Positive Charge",numbins,0,x);
    if(neg!=NULL) delete neg;
        neg  = new TH1F("charge-","Negative Charge",numbins,0,x); 	
    if(sum!=NULL) delete sum;
	sum  = new TH1F("charge","Total Charge",numbins,0,x); 
    
    sum->SetXTitle("t [s]");  neg->SetXTitle("t [s]"); pos->SetXTitle("t [s]");
    sum->SetYTitle("I [arb.]");neg->SetYTitle("I [arb.]");pos->SetYTitle("I [arb.]");
    sum->GetYaxis()->SetTitleOffset(1.4); neg->GetYaxis()->SetTitleOffset(1.4); pos->GetYaxis()->SetTitleOffset(1.4);
    pos->SetLineColor(2);
    neg->SetLineColor(4);
    sum->SetLineColor(1);	

    }
}



KDetector::KDetector()
{
//////////////////////////////////////////////////////////////////////
// Author: Gregor Kramberger                                        //
// Default constructor for KDetector class                           //
// Default = no space charge -> default space charge is function     //
//////////////////////////////////////////////////////////////////////

NeffF=new TF3("Profile","x[0]*x[1]*x[2]*0+[0]",0,1000,0,1000,0,1000);
NeffF->SetParameter(0,0);
NeffH=NULL; //Neff histogram se to NULL

//set magnetic field to zero
for(Int_t i=0;i<3;i++) B[i]=0;

// setting up default random generator - diffusion
ran=new TRandom(33);  

// Calculation parameters
 CalErr=1e-6;
 MaxIter=2000;

// histograms for storing the drift
pos=NULL; neg=NULL; sum=NULL;
SetDriftHisto(25e-9);

// setting up general variables
// multiplication 
TauE=NULL;
TauH=NULL;
taue=-1;   //no electron trapping 
tauh=-1;   //no hole trapping 
MTresh=-1; //no multiplication 
BDTresh=-1;
DiffOffField=8;  // critical field for diffusion to be switched off

// drift

Deps=1e-5;       //precision of tracking
MaxDriftLen=1e9; // maximum driftlenght in [um]

//MobMod=1;  //Mobility parametrization
average=1; //average over waveforms
diff=0;    //diffusion
SStep=1;   //Step size of simulation
Temperature=263; //temperature
BreakDown=0; // no breakdown
Debug=0;     // bo printing of debug information
Voltage2=0;
 Ramo=new KField();
 Real=new KField();

}


KDetector::~KDetector()
{
  //Destructor of the detector class

  if(NeffF!=NULL) delete NeffF;
  if(NeffH!=NULL) delete NeffH;
  if(ran!=NULL) delete ran;
  if(pos!=NULL) delete pos;
  if(neg!=NULL) delete neg;
  if(sum!=NULL) delete sum;
  
}


// K3D 


#include "TArray.h"
#include "TArrayI.h"
#include "TArrayF.h"
#include "TGraph.h"
#include <stdio.h>
#include <stdlib.h>
//#include "KStruct.h"
#include "TMinuit.h"
// #include "KDetector.h"


class K3D : public KDetector
{
private:
public:
    Int_t Col;
    Float_t CellZ;
    Float_t CellX;
    Float_t CellY;

    Float_t *PosD; //[Col]
    Float_t *PosX; //[Col]
    Float_t *PosY; //[Col]
    Float_t *PosR; //[Col]
    Short_t *PosW; //[Col]
    Short_t *PosM; //[Col]

    K3D(Int_t, Float_t = 100, Float_t = 100, Float_t = 105);
    ~K3D();

    // ClassDef(K3D, 1)
};

// ClassImp(K3D)
    //////////////////////////////////////////////////////////////////////////
    //                                                                      //
    // K3D                                                                  //
    //                                                                      //
    // Description of the 3D detector.                                      //
    // The geometry of the                                                  //
    // is defined by                                                        //
    // dasdasd                                                              //
    //////////////////////////////////////////////////////////////////////////

    K3D::~K3D()
{
    delete PosD;
    delete PosX;
    delete PosY;
    delete PosR;
    delete PosW;
}

K3D::K3D(Int_t x1, Float_t x, Float_t y, Float_t z)
{
    Col = x1;
    PosD = new Float_t[Col];
    PosX = new Float_t[Col];
    PosY = new Float_t[Col];
    PosR = new Float_t[Col];
    PosW = new Short_t[Col];
    PosM = new Short_t[Col];

    for (Int_t i = 0; i < Col; i++) {
        PosD[i] = 0;
        PosX[i] = 0;
        PosY[i] = 0;
        PosR[i] = 0;
        PosW[i] = 0;
        PosM[i] = 0;
    }

    CellZ = z;
    CellX = x;
    CellY = y;
}






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

    // Start the Test3D_SiC_One
    gStyle->SetCanvasPreferGL(kTRUE);
    // define a 3D detector with 5 electrodes
    // x=100 , y is 50 and thickness 120
    K3D *det = new K3D(7, 80, 80, 300);
    
    theApp.Run();
    }

    #endif

