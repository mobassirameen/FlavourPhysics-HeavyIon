#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLeaf.h"
#include "TEfficiency.h"
#include "TROOT.h"
#include <sstream>
#include <fstream>
#include <iostream>
#include <cstdio>
#include "TStyle.h"
#include "TLegend.h"
#include "CMS_lumi.h"

using namespace std;

void initRootStyle();
void BinLogX(TH1*);
double fit_noise(double*, double*);
Double_t sumOfGauss(Double_t*, Double_t*);

TH1F* spectrum_noise;

const double C = 0.1;

void trigger_eff(int runnumber=286178){
  initRootStyle();

  writeExtraText = true;       // if extra text
  extraText  = "Preliminary";  // default extra text is "Preliminary"
  lumi_sqrtS = "PbPb 5.02 TeV";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
  int iPeriod = 0;

  TCanvas *c1 = new TCanvas();
  TH1::SetDefaultSumw2();

  // Name of directory to plot
  //TFile *f = new TFile(Form("digitree_%d.root",runnumber)); // opening the root file
  TFile *f = new TFile("forward_326483.root"); // opening the root file
  //TTree *ZDCRecTree = (TTree*)f->Get("ZDCRecTree"); // reading ZDC rec tree
  TTree *ZDCDigiTree = (TTree*)f->Get("analyzer/zdcdigi"); // reading ZDC digi tree
  const int NSIDE=2; const char* stit[NSIDE] = {"#minus","#plus"};  const char* stit2[NSIDE] = {"neg","pos"};
  const int NTYPE=2; const char* ttit[NTYPE] = {"EM","HAD"};
  const int NCH=5; const char* ctit[NTYPE][NCH] = {
                                                    {"1","2","3","4","5"}, //HD sections run only 1-4
                                                    {"1","2","3","4","5"} //EM sections run 1-5
                                                  };
 

  TEfficiency* rate = new TEfficiency("eff",";ADC;Rate [ZB rate]",250,0,250); 
  TH2F* trigger = new TH2F("","",2,0,2,2,0,2);

  const int NTS=10;            // number of timeslices
  TLeaf* bxidLeaf = (TLeaf*) ZDCDigiTree->GetLeaf("bxid");
  TLeaf* zsideLeaf = (TLeaf*) ZDCDigiTree->GetLeaf("zside");
  TLeaf* sectionLeaf = (TLeaf*) ZDCDigiTree->GetLeaf("section");
  TLeaf* channelLeaf = (TLeaf*) ZDCDigiTree->GetLeaf("channel");
  TLeaf* zerobias = (TLeaf*) ZDCDigiTree->GetLeaf("HLT_HIZeroBias_v1");
  TLeaf* ntrk = (TLeaf*) ZDCDigiTree->GetLeaf("nTrack");
  TLeaf* nHFneg = (TLeaf*) ZDCDigiTree->GetLeaf("nHFneg");
  TLeaf* nHFpos = (TLeaf*) ZDCDigiTree->GetLeaf("nHFpos");
  
  TLeaf* adcLeaf[NTS];
  TLeaf* fCleaf[NTS];

  TLeaf* minbiasLeaf[20];

  for(int i = 0; i < 20; i++)
    minbiasLeaf[i] = (TLeaf*) ZDCDigiTree->GetLeaf(Form("HLT_HIMinimumBias_part%d_v1",i));


  double w[2][2][5] = {
    //{{0.0198368,0.00362876,0.139074,0.165967,0.0825955},{1,0.924266,0.982671,0.894771}},
    //{{1.0,1.0,1.0,1.0,1.0},{1.0,0.617/0.4,0.315/0.23,0.259/0.17,0.0}},  // negative side
    {{1.0,1.0,1.0,1.0,1.0},{1.0,0.587/0.95,0.307/0.45,0.250/0.3,0.0}},  // negative side
    //{{0.1,0.1,0.1,0.1,0.1},{1.0,0.249669,0.706763,0.828292,0.0}},
    
    //{{1.0,1.0,1.0,1.0,1.0},{1.0,0.618,0.315/0.5,0.259/0.33,0.0}}}; // positive side
    {{1.0,1.0,1.0,1.0,1.0},{1.0,0.618/0.75,0.315/0.4,0.259/0.3,0.0}}}; // positive side

  for(int iTS = 0; iTS < NTS; iTS++){
    adcLeaf[iTS] = (TLeaf*) ZDCDigiTree->GetLeaf(Form("adc%d",iTS));
    fCleaf[iTS] = (TLeaf*) ZDCDigiTree->GetLeaf(Form("nfC%d",iTS));
  }

  std::ofstream outfile("singleNeutron.dat");

  for(int i = 0; i < 1000000/*ZDCDigiTree->GetEntries()*/; i++){
    ZDCDigiTree->GetEntry(i);

    bool isMinBias = false;
    bool zdc_pos = false, zdc_neg = false;

    for(int j = 0; j < 20; j++)
      if(minbiasLeaf[j]->GetValue() == 1) isMinBias = true;

    if(zerobias->GetValue()){
      for(int thr = 0; thr < 250; thr++){
        bool passed = false;
        for(int n = 0; n < 50; n++){
          int side = (int)((zsideLeaf->GetValue(n)+1)/2.0);
          int type = (int)(sectionLeaf->GetValue(n))-1;
          int channel = (int)(channelLeaf->GetValue(n))-1;
          int adc = (int)(adcLeaf[4]->GetValue(n));

          if(type == 1 && adc > thr) passed = true;
          if(type == 1 && side == 1 && (channel == 0 || channel == 1) && adc > 100) zdc_pos = true;
          if(type == 1 && side == 0 && channel == 0 && adc > 119) zdc_neg = true;
          if(type == 1 && side == 0 && channel == 1 && adc > 100) zdc_neg = true;

        }
        rate->Fill(passed,thr);
      }
    
      if((zdc_pos && zdc_neg) && isMinBias)
        trigger->Fill(1.5,1.5);
      if(!(zdc_pos && zdc_neg) && !isMinBias)
        trigger->Fill(0.5,0.5);
      if(!(zdc_pos && zdc_neg) && isMinBias)
        trigger->Fill(0.5,1.5);
      if((zdc_pos && zdc_neg) && !isMinBias)
        trigger->Fill(1.5,0.5);
    }



    if(i % 100000 == 0) std::cout << i << " events are processed." << std::endl;
  }

  ///////////////////////////////////////
  // Formatting and drawing histograms //
  ///////////////////////////////////////


  //c1->SetLogy();

  TFile f2("eff.root","RECREATE"); 

  rate->SetMarkerColor(1);
  rate->SetMarkerStyle(20);
  rate->SetMarkerSize(0.75);

  rate->Draw("p a e");

  c1->SaveAs("ZDC_figures/forward_326483/eff/rate.png");
  c1->SaveAs("ZDC_figures/forward_326483/eff/rate.pdf");

  rate->Write();

  trigger->SetMarkerColor(1);
  trigger->SetMarkerSize(1.8);
  trigger->Draw("col");
  trigger->Draw("text same");

  c1->SaveAs("trigger.png");


  return;
}

double fit_noise(double* x, double* p){
  int N = 10;
  double f = 0.0;

  for(int i = 1; i < N; i++)
    f += p[i+2] * TMath::Gaus(x[0],i*p[1],sqrt(i)*p[2],kTRUE);
  
  return p[0]*f + spectrum_noise->Interpolate(x[0]);
}

Double_t sumOfGauss(Double_t* x, Double_t* p){
  int N = 100;
  double f = 0.0;

  for(int i = 1; i < N; i++)
    f += p[i+1] * TMath::Gaus(x[0],i*p[0],sqrt(i)*p[1],kTRUE);
  
  return f;
}

void initRootStyle(){
  //  using namespace RooFit ;

  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  gStyle->SetOptFit(0);

  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetStatColor(0);

  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetFrameLineColor(kBlack);

  gStyle->SetCanvasColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasBorderSize(0);

  gStyle->SetPadColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadBorderSize(0);

  gStyle->SetLegendBorderSize(0);
  //gStyle->SetTextSize(0.04);
  gStyle->SetTextFont(42);
  gStyle->SetLabelFont(42,"xyz");
  gStyle->SetTitleFont(42,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadTopMargin(0.10);
  gStyle->SetPadRightMargin(0.12); // 0.10
  gStyle->SetPadLeftMargin(0.15); // 0.12

  gStyle->SetTitleXOffset(1.1);
  gStyle->SetTitleYOffset(1.5); // 1.2

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  gStyle->SetCanvasDefH(600);
  gStyle->SetCanvasDefW(600);
/*
  gStyle->SetStatX(0.92); // 0.36
  gStyle->SetStatY(0.92);
*/
  //gStyle->SetHistMinimumZero(kTRUE);

  //gStyle->SetErrorX(0); //disable if you want to draw horizontal error bars, e.g. when having variable bin size
  gStyle->SetEndErrorSize(1);

  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerColor(2);
  gStyle->SetMarkerSize(0.2);

  //gROOT->ForceStyle();

  std::cout << "ROOT style loaded." << std::endl;
}

void BinLogX(TH1* h){
  TAxis *axis = h->GetXaxis();
  int bins = axis->GetNbins();

  Axis_t from = axis->GetXmin();
  Axis_t to = axis->GetXmax();
  Axis_t width = (to - from) / (double) bins;
  Axis_t *new_bins = new Axis_t[bins + 1];

  for (int i = 0; i <= bins; i++){
    new_bins[i] = TMath::Power(10, from + i * width);
  }
  axis->Set(bins, new_bins);
  delete new_bins;
}

