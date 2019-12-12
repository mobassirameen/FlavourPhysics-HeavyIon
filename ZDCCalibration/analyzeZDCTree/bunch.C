#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLeaf.h"
#include "TROOT.h"
#include <sstream>
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

const double MU = 59.92;

void bunch(int runnumber=286178){
  initRootStyle();

  writeExtraText = true;       // if extra text
  extraText  = "Preliminary";  // default extra text is "Preliminary"
  lumi_sqrtS = "PbPb 5.02 TeV";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
  int iPeriod = 0;

  TCanvas *c1 = new TCanvas();
  TH1::SetDefaultSumw2();

  // Name of directory to plot
  //TFile *f = new TFile(Form("digitree_%d.root",runnumber)); // opening the root file
  TFile *f = new TFile("express_326176.root"); // opening the root file
  //TTree *ZDCRecTree = (TTree*)f->Get("ZDCRecTree"); // reading ZDC rec tree
  TTree *ZDCDigiTree = (TTree*)f->Get("analyzer/zdcdigi"); // reading ZDC digi tree
  const int NSIDE=2; const char* stit[NSIDE] = {"#minus","#plus"};  const char* stit2[NSIDE] = {"neg","pos"};
  const int NTYPE=2; const char* ttit[NTYPE] = {"EM","HAD"};
  const int NCH=5; const char* ctit[NTYPE][NCH] = {
                                                    {"1","2","3","4","5"}, //HD sections run only 1-4
                                                    {"1","2","3","4","5"} //EM sections run 1-5
                                                  };
  TH1F* bx_ZDC[2][2][5];
  /*TH1F* bx2_ZDC[2][2][5];
  TH1F* bx3_ZDC[2][2][5];
  TH1F* bx4_ZDC[2][2][5];*/

  TH1F* bx_RPD[2][16];
  /*TH1F* bx2_RPD[2][16];
  TH1F* bx3_RPD[2][16];
  TH1F* bx4_RPD[2][16];*/

  for(int i = 0; i < 2; i++){
    for(int j = 0; j < 2; j++)
      for(int k = 0; k < 5; k++){
        bx_ZDC[i][j][k] = new TH1F(Form("bx1 ZDC %s%s channel %d",ttit[j],stit[i],k+1),Form("ZDC %s%s channel %d;BX ID;sum Q_{TS4} [fC]",ttit[j],stit[i],k+1),3600,0,3600);
/*        bx2_ZDC[i][j][k] = new TH1F(Form("bx2 ZDC %s%s channel %d",ttit[j],stit[i],k+1),Form("ZDC %s%s channel %d;BX ID;sum Q_{TS4} [fC]",ttit[j],stit[i],k+1),1000,1000,2000);
        bx3_ZDC[i][j][k] = new TH1F(Form("bx3 ZDC %s%s channel %d",ttit[j],stit[i],k+1),Form("ZDC %s%s channel %d;BX ID;sum Q_{TS4} [fC]",ttit[j],stit[i],k+1),1000,2000,3000);
        bx4_ZDC[i][j][k] = new TH1F(Form("bx4 ZDC %s%s channel %d",ttit[j],stit[i],k+1),Form("ZDC %s%s channel %d;BX ID;sum Q_{TS4} [fC]",ttit[j],stit[i],k+1),1000,3000,4000);*/
      }
    for(int j = 0; j < 16; j++){
      bx_RPD[i][j] = new TH1F(Form("bx1 RPD%s channel %d",stit[i],j+1),Form("RPD%s channel %d;BX ID;sum Q_{TS4} [fC]",stit[i],j+1),3600,0,3600);
/*      bx2_RPD[i][j] = new TH1F(Form("bx2 RPD%s channel %d",stit[i],j+1),Form("RPD%s channel %d;BX ID;sum Q_{TS4} [fC]",stit[i],j+1),1000,1000,2000);
      bx3_RPD[i][j] = new TH1F(Form("bx3 RPD%s channel %d",stit[i],j+1),Form("RPD%s channel %d;BX ID;sum Q_{TS4} [fC]",stit[i],j+1),1000,2000,3000);
      bx4_RPD[i][j] = new TH1F(Form("bx4 RPD%s channel %d",stit[i],j+1),Form("RPD%s channel %d;BX ID;sum Q_{TS4} [fC]",stit[i],j+1),1000,3000,4000);*/
    }
  }

  TH1F* event = new TH1F("",";TS [25 ns]; Q [fC]",10,0,10);

  const int NTS=10;            // number of timeslices
  TLeaf* bxidLeaf = (TLeaf*) ZDCDigiTree->GetLeaf("bxid");
  TLeaf* zsideLeaf = (TLeaf*) ZDCDigiTree->GetLeaf("zside");
  TLeaf* sectionLeaf = (TLeaf*) ZDCDigiTree->GetLeaf("section");
  TLeaf* channelLeaf = (TLeaf*) ZDCDigiTree->GetLeaf("channel");
  TLeaf* random = (TLeaf*) ZDCDigiTree->GetLeaf("HLT_HIRandom_v1");
  TLeaf* ntrk = (TLeaf*) ZDCDigiTree->GetLeaf("nTrack");
  TLeaf* nHFneg = (TLeaf*) ZDCDigiTree->GetLeaf("nHFneg");
  TLeaf* nHFpos = (TLeaf*) ZDCDigiTree->GetLeaf("nHFpos");
  
  TLeaf* adcLeaf[NTS];
  TLeaf* fCleaf[NTS];

  for(int iTS = 0; iTS < NTS; iTS++){
    adcLeaf[iTS] = (TLeaf*) ZDCDigiTree->GetLeaf(Form("adc%d",iTS));
    fCleaf[iTS] = (TLeaf*) ZDCDigiTree->GetLeaf(Form("nfC%d",iTS));
  }
  
  for(int i = 0; i < ZDCDigiTree->GetEntries(); i++){
    ZDCDigiTree->GetEntry(i);

    if(random->GetValue()){
      for(int n = 0; n < 50; n++){
        int side = (int)((zsideLeaf->GetValue(n)+1)/2.0);
        int type = (int)(sectionLeaf->GetValue(n))-1;
        int channel = (int)(channelLeaf->GetValue(n))-1;

        if(type != 3){ // EM or HAD section
          for(int iTS = 0; iTS < 10; iTS++)
            bx_ZDC[side][type][channel]->Fill(bxidLeaf->GetValue()-4+iTS,fCleaf[iTS]->GetValue(n));
/*          bx2_ZDC[side][type][channel]->Fill(bxidLeaf->GetValue(),fCleaf[4]->GetValue(n));
          bx3_ZDC[side][type][channel]->Fill(bxidLeaf->GetValue(),fCleaf[4]->GetValue(n));
          bx4_ZDC[side][type][channel]->Fill(bxidLeaf->GetValue(),fCleaf[4]->GetValue(n));*/
       }
        else{ // RPD section
          for(int iTS = 0; iTS < 10; iTS++)
            bx_RPD[side][channel]->Fill(bxidLeaf->GetValue()-4+iTS,fCleaf[iTS]->GetValue(n));
/*          bx2_RPD[side][channel]->Fill(bxidLeaf->GetValue(),fCleaf[4]->GetValue(n));
          bx3_RPD[side][channel]->Fill(bxidLeaf->GetValue(),fCleaf[4]->GetValue(n));      
          bx4_RPD[side][channel]->Fill(bxidLeaf->GetValue(),fCleaf[4]->GetValue(n));*/
        }
      }
    }

    if(i % 100000 == 0) std::cout << i << " events are processed." << std::endl;
  }

  ///////////////////////////////////////
  // Formatting and drawing histograms //
  ///////////////////////////////////////


  //c1->SetLogy();

  TFile f2("bx.root","RECREATE"); 

  // ZDC em
  for(int i = 0; i < 2; i++)
    for(int k = 0; k < 5; k++){
      bx_ZDC[i][0][k]->Draw("hist");
      c1->SaveAs(Form("ZDC_figures/bx/bx_EM_%s_channel_%d.png",stit2[i],k+1));
      bx_ZDC[i][0][k]->Write();
  /*    bx2_ZDC[i][0][k]->Draw("hist");
      c1->SaveAs(Form("ZDC_figures/bx/bx2_EM_%s_channel_%d.png",stit2[i],k+1));
      bx3_ZDC[i][0][k]->Draw("hist");
      c1->SaveAs(Form("ZDC_figures/bx/bx3_EM_%s_channel_%d.png",stit2[i],k+1));
      bx4_ZDC[i][0][k]->Draw("hist");
      c1->SaveAs(Form("ZDC_figures/bx/bx4_EM_%s_channel_%d.png",stit2[i],k+1));*/
    }

  // ZDC em
  for(int i = 0; i < 2; i++)
    for(int k = 0; k < 4; k++){
      bx_ZDC[i][1][k]->Draw("hist");
      c1->SaveAs(Form("ZDC_figures/bx/bx_HAD_%s_channel_%d.png",stit2[i],k+1));
      bx_ZDC[i][1][k]->Write();
      /*bx2_ZDC[i][1][k]->Draw("hist");
      c1->SaveAs(Form("ZDC_figures/bx/bx2_HAD_%s_channel_%d.png",stit2[i],k+1));
      bx3_ZDC[i][1][k]->Draw("hist");
      c1->SaveAs(Form("ZDC_figures/bx/bx3_HAD_%s_channel_%d.png",stit2[i],k+1));
      bx4_ZDC[i][1][k]->Draw("hist");
      c1->SaveAs(Form("ZDC_figures/bx/bx4_HAD_%s_channel_%d.png",stit2[i],k+1));*/
    }

  // RPD
  for(int i = 0; i < 2; i++)
    for(int k = 0; k < 16; k++){
      bx_RPD[i][k]->Draw("hist");
      c1->SaveAs(Form("ZDC_figures/bx/bx_RPD_%s_channel_%d.png",stit2[i],k+1));
      bx_RPD[i][k]->Write();
      /*bx2_RPD[i][k]->Draw("hist");
      c1->SaveAs(Form("ZDC_figures/bx/bx2_RPD_%s_channel_%d.png",stit2[i],k+1));
      bx3_RPD[i][k]->Draw("hist");
      c1->SaveAs(Form("ZDC_figures/bx/bx3_RPD_%s_channel_%d.png",stit2[i],k+1));
      bx4_RPD[i][k]->Draw("hist");
      c1->SaveAs(Form("ZDC_figures/bx/bx4_RPD_%s_channel_%d.png",stit2[i],k+1));*/
    }


  /*TLegend* leg2 = new TLegend(0.50,0.65,0.70,0.85);
  leg2->SetTextFont(42);
  leg2->SetTextSize(0.033);
  leg2->AddEntry(spectrum,Form("Data, run: %d",runnumber));
  leg2->AddEntry(fit,"Sum of Gaussians:");
  leg2->AddEntry((TObject*)0,"#mu_{n} = n #mu_{0}","0");
  leg2->AddEntry((TObject*)0,"#sigma_{n} = #sqrt{n} #sigma_{0}","0");

  leg2->Draw("same");

  CMS_lumi(c1,0);*/

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
  gStyle->SetCanvasDefW(1200);
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

