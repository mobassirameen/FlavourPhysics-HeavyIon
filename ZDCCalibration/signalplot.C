
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TF1.h>
#include <TH1D.h>
#include "TGraphAsymmErrors.h"
#include "TVirtualFFT.h"
#include "TBinomialEfficiencyFitter.h"
#include "TVectorF.h"
#include "TPaveText.h"
#include <vector>
#include "TTree.h"
#include "TFile.h"
#include "TEventList.h"
#include "Riostream.h"
#include "string.h"
#include "TList.h"
#include "TDirectory.h"
#include "TCut.h"
#include "TChain.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TEfficiency.h"
#include "TLorentzVector.h"
#include "TProfile.h"
#include "TLegend.h"

using namespace std;





//______________________________________________________________________________
void signalplot()
{
   //Double_t x1[10], x2[10], x3[10], z[10], errorz[10], TS[10] ;
   Double_t  errorz[10];
   Double_t a1=2195.17, b1= 17984.4, c1= 2165.44, d1=245.34;
   Double_t a2=4131.0825, b2= 32094.1250, c2=2335.6399, d2=-145.4329;
   const Int_t  n=10; 
   //Time slices 
   Float_t TS[n] = {25, 50,75, 100, 125, 150, 175, 200, 225, 250};
   Float_t chifit[n], anafit[n];

// The z values data points, 
   //Float_t z[n]={ 139.783,160.45, 200, 288.2, 34793.6, 6356 , 1571.73, 2434.73, 1045, 1391.97};
  // Float_t z[n]={120.8, 40.9667, 52.6, 61.9, 8173.15, 456.833, 599.233, 320.433, 206.567, 172.467};
  // Float_t z[n]={113.567,94.9667, 70.1667, 116.15, 12232.7, 2180.95, 471.033, 456.85, 402.6, 464.067};

 //Float_t z[n]={278.123, 23067.7, 1139.6, 581.871, 26731.9, 1771.06, 966.376, 2375.01, 591.682, 655.143};
 //Float_t z[n]={302.53, 9597.55, 679.543, 450.421, 28027.2, 1471.81, 827.59, 1206.05, 436.99, 459.506};
     //Float_t z[n]={275.002, 3569.48, 273.291,242.707,7502.73,409.144,273.855,190.283,156.368,187.645};
    //Float_t z[n]={408.786, 5485.92, 1116.4, 550.394, 41180.7, 7343.13, 2069, 3249.8, 1172.92, 1657.21};
    //Float_t z[n]={171.64, 3050.79, 404.998, 237.17, 14554.6, 1615.07,617.18, 477.718,268.29,429.551};
Float_t z[n]={259.527, 3660.31, 366.701, 318.493, 32096.7, 1804.68, 1017.79, 2106.04, 570.553, 689.849};
// The errors on z values 5% of the tprofile values of the data points
   
   errorz[0]=5.137;
    errorz[1]=113.0;
   errorz[2]=13.249;
   errorz[3]=10.5;
   errorz[4]=1394.99;
   errorz[5]=153.85;
   errorz[6]=45.69;
   errorz[7]=45.69;
   errorz[8]=24.06;
   errorz[9]=32.177;
// the x1 values
   Float_t x1[n]={ 0, 1, 0.1, 0.05, 0, 0, 0, 0, 0,0};
   
// the x2 values
   Float_t x2[n] = { 0, 0, 0, 0, 1.0, 0.1, 0.05,0,0,0};
   
   
// the x3 values
   Float_t x3[n] = {0, 0, 0, 0, 0, 0, 0,1,0.1,0.05};
    for(Int_t j=0; j<10; j++){
     chifit[j] = (a1 * x1[j]) + (b1 * x2[j]) + (c1 * x3[j]) + d1;
     anafit[j] = (a2 * x1[j]) + (b2 * x2[j]) + (c2 * x3[j]) + d2; 
  }
  // TCanvas *c11 = new TCanvas("c11","A Simple Graph Example",200,10,700,700);
   TGraph *gr1 = new TGraph(n,TS,z);
   gr1->SetMarkerColor(kRed);
   gr1->SetMarkerStyle(21);
   gr1->SetLineWidth(2);
   gr1->SetLineColor(kRed);
   //gr1->Draw("AC");

  // TCanvas *c22 = new TCanvas("c22","A Simple Graph Example",200,10,700,700);
   TGraph *gr2 = new TGraph(n,TS,chifit);
   gr2->SetMarkerColor(kBlue);
   gr2->SetMarkerStyle(20);
   gr2->SetLineWidth(2);
   gr2->SetLineColor(kBlue);
  // gr2->Draw("AC");
  // TCanvas *c3 = new TCanvas("c3","A Simple Graph Example",200,10,700,700);
   TGraph *gr3 = new TGraph(n,TS,anafit);
   gr3->SetMarkerColor(kGreen);
   gr3->SetMarkerStyle(22);
   gr3->SetLineWidth(2);
   gr3->SetLineColor(kGreen);
   //gr3->Draw("cloz");
   TCanvas *cm = new TCanvas("cm","Multigraph",200,10,700,700);
   gStyle->SetOptStat(0);
   TMultiGraph *mg = new TMultiGraph();
   mg->Add(gr1);
   mg->Add(gr2);
   mg->Add(gr3);
   mg->Draw("AC*");
   TLegend *legend = new TLegend(0.1,0.7,0.48,0.9);
   //legend->SetHeader("Header","C");             // option "C" allows to center the header
   legend->AddEntry(gr1,"Data");
   legend->AddEntry(gr2,"TMinuiteClass");
   legend->AddEntry(gr3,"Analytical method");
   legend->Draw();
   cm->Update();

}

