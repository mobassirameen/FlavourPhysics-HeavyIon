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

const double C = 1.0;

void neutron_peaks_minbias(int runnumber=286178){
  initRootStyle();

  writeExtraText = true;       // if extra text
  extraText  = "Preliminary";  // default extra text is "Preliminary"
  lumi_sqrtS = "PbPb 5.02 TeV";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
  int iPeriod = 0;

  TCanvas *c1 = new TCanvas();
  TH1::SetDefaultSumw2();

  // Name of directory to plot
  //TFile *f = new TFile(Form("digitree_%d.root",runnumber)); // opening the root file
  TFile *f = new TFile("minbias_326483.root"); // opening the root file
  //TTree *ZDCRecTree = (TTree*)f->Get("ZDCRecTree"); // reading ZDC rec tree
  TTree *ZDCDigiTree = (TTree*)f->Get("analyzer/zdcdigi"); // reading ZDC digi tree
  const int NSIDE=2; const char* stit[NSIDE] = {"#minus","#plus"};  const char* stit2[NSIDE] = {"neg","pos"};
  const int NTYPE=2; const char* ttit[NTYPE] = {"EM","HAD"};
  const int NCH=5; const char* ctit[NTYPE][NCH] = {
                                                    {"1","2","3","4","5"}, //HD sections run only 1-4
                                                    {"1","2","3","4","5"} //EM sections run 1-5
                                                  };
 
  TH1F* totEnergy_pos = new TH1F("total energy pos",";Energy [TeV];Entries",100,0,200);
  TH1F* totEnergy_neg = new TH1F("total energy neg",";Energy [TeV];Entries",100,0,200);

  TH2F* totEnergy_2D = new TH2F("total energy 2d pos",";ZDC#plus Energy [TeV];ZDC#minus Energy [TeV]",100,0,20,100,0,20);

  TH2F* tune_2D_pos = new TH2F("tune 2d pos",";HAD + c #cdot EM [a.u.];EM [a.u.]",100,0,25000,100,0,2500);
  TH2F* tune_2D_neg = new TH2F("tune 2d neg",";HAD + c #cdot EM [a.u.];EM [a.u.]",100,0,25000,100,0,2500);


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

  double w[2][2][5] = {
    //{{0.0198368,0.00362876,0.139074,0.165967,0.0825955},{1,0.924266,0.982671,0.894771}},
    //{{1.0,1.0,1.0,1.0,1.0},{1.0,0.617/0.4,0.315/0.23,0.259/0.17,0.0}},  // negative side
    {{0.1,0.1,0.1,0.1,0.1},{1.0,0.587/0.8,0.307/0.65,0.250/0.45,0.0}},  // negative side
    //{{0.1,0.1,0.1,0.1,0.1},{1.0,0.249669,0.706763,0.828292,0.0}},
    
    //{{1.0,1.0,1.0,1.0,1.0},{1.0,0.618,0.315/0.5,0.259/0.33,0.0}}}; // positive side
    {{0.1,0.1,0.1,0.1,0.1},{1.0,0.618/1.0,0.315/0.65,0.259/0.5,0.0}}}; // positive side

  for(int iTS = 0; iTS < NTS; iTS++){
    adcLeaf[iTS] = (TLeaf*) ZDCDigiTree->GetLeaf(Form("adc%d",iTS));
    fCleaf[iTS] = (TLeaf*) ZDCDigiTree->GetLeaf(Form("nfC%d",iTS));
  }

  std::ofstream outfile_pos("singleNeutron_pos.dat");
  std::ofstream outfile_neg("singleNeutron_neg.dat");

  for(int i = 0; i < ZDCDigiTree->GetEntries(); i++){
    ZDCDigiTree->GetEntry(i);

    double had1_pos, had2_pos, had3_pos, had4_pos;
    double had1_neg, had2_neg, had3_neg, had4_neg;
    double totHAD[2] = {0,0};
    double totEM[2] = {0,0};

    bool isMaxHAD1_pos = false, isMaxHAD2_pos = false/*, isMaxHAD3 = false, isMaxHAD4 = false*/;
    bool isMaxHAD1_neg = false, isMaxHAD2_neg = false/*, isMaxHAD3 = false, isMaxHAD4 = false*/;

    int iMax;
    double qMax;
    int adcMaxHAD1_pos, adcMaxHAD2_pos;
    int adcMaxHAD1_neg, adcMaxHAD2_neg;

    double signal_array_pos[9] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
    double signal_array_neg[9] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

    bool pileup = false;

    for(int n = 0; n < 50; n++){
      int side = (int)((zsideLeaf->GetValue(n)+1)/2.0);
      int type = (int)(sectionLeaf->GetValue(n))-1;
      int channel = (int)(channelLeaf->GetValue(n))-1;
      double signal = fCleaf[4]->GetValue(n);
  
      qMax = 0;

      for(int iTS = 3; iTS < 8; iTS++){
        if(fCleaf[iTS]->GetValue(n) > qMax){
          qMax = fCleaf[iTS]->GetValue(n);
          iMax = iTS;
        }
      }

      if(side == 1 && type == 1 && channel == 0 && iMax == 4 && adcLeaf[iMax]->GetValue(n) < 254 && adcLeaf[iMax]->GetValue(n) > 120){
        isMaxHAD1_pos = true;
        adcMaxHAD1_pos = adcLeaf[iMax]->GetValue(n);
      }
      if(side == 1 && type == 1 && channel == 1 && iMax == 4 && adcLeaf[iMax]->GetValue(n) < 254 && adcLeaf[iMax]->GetValue(n) > 120){
        isMaxHAD2_pos = true;
        adcMaxHAD2_pos = adcLeaf[iMax]->GetValue(n);
      }

      if(side == 0 && type == 1 && channel == 0 && iMax == 4 && adcLeaf[iMax]->GetValue(n) < 254 && adcLeaf[iMax]->GetValue(n) > 120){
        isMaxHAD1_neg = true;
        adcMaxHAD1_neg = adcLeaf[iMax]->GetValue(n);
      }
      if(side == 0 && type == 1 && channel == 1 && iMax == 4 && adcLeaf[iMax]->GetValue(n) < 254 && adcLeaf[iMax]->GetValue(n) > 120){
        isMaxHAD2_neg = true;
        adcMaxHAD2_neg = adcLeaf[iMax]->GetValue(n);
      }

/*      if(side == 0 && type == 1 && channel == 2 && (iMax == 4 || iMax == 5) && adcLeaf[iMax]->GetValue(n) < 254 && adcLeaf[iMax]->GetValue(n) > 120)
        isMaxHAD3 = true;
      if(side == 0 && type == 1 && channel == 3 && (iMax == 4 || iMax == 5) && adcLeaf[iMax]->GetValue(n) < 254 && adcLeaf[iMax]->GetValue(n) > 120)
        isMaxHAD4 = true;*/

      if(type == 1/* && iMax == 4*/){
        totHAD[side] += w[side][type][channel] * signal;
        if(side == 0)
          signal_array_neg[5+channel] = signal;
        else
          signal_array_pos[5+channel] = signal;
      }
      else if(type == 0/* && iMax == 4*/){
        totEM[side] +=  w[side][type][channel] * signal;
        if(side == 0)
          signal_array_neg[channel] = signal;
        else
          signal_array_pos[5+channel] = signal;
      }
      /*else if(type == 3)
        totRPD[iside] += signal;*/
      if(fCleaf[0]->GetValue(n) > 1000) pileup = true;
    }

    if(isMaxHAD1_pos || isMaxHAD2_pos /*&& adcMaxHAD1_pos < 254 && adcMaxHAD2_pos < 254 && adcMaxHAD1_pos > 50 && adcMaxHAD2_pos > 50*/){
      tune_2D_pos->Fill(totHAD[1]+C*totEM[1],totEM[1]);
      totEnergy_pos->Fill((totHAD[1]+C*totEM[1])/7500.0 * 2.51);
      //totEnergy_2D->Fill(totHAD[1]+C*totEM[1]/*+R*totRPD[1]*/,totHAD[0]+C*totEM[0]/*+R*totRPD[0]*/);
      if(totHAD[1]+C*totEM[1] < 10000){
        for(int ich = 0; ich < 9; ich++)
          outfile_pos << signal_array_pos[ich] << "\t";
        outfile_pos << std::endl;
      }  
    }

    if(isMaxHAD1_neg || isMaxHAD2_neg /*&& adcMaxHAD1_neg < 254 && adcMaxHAD2_neg < 254 && adcMaxHAD1_neg > 50 && adcMaxHAD2_neg > 50 && !pileup*/){
      tune_2D_neg->Fill(totHAD[0]+C*totEM[0],totEM[0]);
      totEnergy_neg->Fill((totHAD[1]+C*totEM[1])/11000.0 * 2.51/*+R*totRPD[1]*/);
      //
      if(totHAD[0]+C*totEM[0] < 15000){
        for(int ich = 0; ich < 9; ich++)
          outfile_neg << signal_array_neg[ich] << "\t";
        outfile_neg << std::endl;
      }
    }

    //if((isMaxHAD1_pos || isMaxHAD2_pos) || (isMaxHAD1_neg || isMaxHAD2_neg))
      totEnergy_2D->Fill((totHAD[1]+C*totEM[1])/7500.0 * 2.51 /*+R*totRPD[1]*/,(totHAD[0]+C*totEM[0])/11000.0 * 2.51 /*+R*totRPD[0]*/);     

    if(i % 100000 == 0) std::cout << i << " events are processed." << std::endl;
  }

  ///////////////////////////////////////
  // Formatting and drawing histograms //
  ///////////////////////////////////////



  TFile f2("ratio.root","RECREATE"); 

  c1->SetLogz();
  tune_2D_pos->Draw("colz");
  c1->SaveAs("ZDC_figures/minbias_326483/neutron_peaks/tune_2D_pos.pdf");
  c1->SaveAs("ZDC_figures/minbias_326483/neutron_peaks/tune_2D_pos.png");
  tune_2D_pos->Write();

  tune_2D_neg->Draw("colz");
  c1->SaveAs("ZDC_figures/minbias_326483/neutron_peaks/tune_2D_neg.pdf");
  c1->SaveAs("ZDC_figures/minbias_326483/neutron_peaks/tune_2D_neg.png");
  tune_2D_neg->Write();

  //c1->SetLogy();

  totEnergy_neg->SetLineColor(2);
  totEnergy_neg->SetMarkerColor(2);
  totEnergy_neg->SetMarkerSize(0.75);
  totEnergy_neg->SetMarkerStyle(20);
  totEnergy_neg->Draw("p e");

  totEnergy_pos->SetMarkerColor(1);
  totEnergy_pos->SetMarkerSize(0.75);
  totEnergy_pos->SetMarkerStyle(20);
  totEnergy_pos->Draw("p e same");

  totEnergy_pos->Write();

  c1->SaveAs("ZDC_figures/minbias_326483/neutron_peaks/totEnergy.pdf");
  c1->SaveAs("ZDC_figures/minbias_326483/neutron_peaks/totEnergy.png");
  totEnergy_neg->Write();

  totEnergy_2D->Draw("colz");
  c1->SaveAs("ZDC_figures/minbias_326483/neutron_peaks/totEnergy_2D.pdf");
  c1->SaveAs("ZDC_figures/minbias_326483/neutron_peaks/totEnergy_2D.png");
  totEnergy_2D->Write();

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

