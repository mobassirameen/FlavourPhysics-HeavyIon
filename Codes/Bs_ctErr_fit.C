#include "TStyle.h"
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooAbsPdf.h"
#include "RooFitResult.h"
#include "RooMCStudy.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooPolynomial.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TTree.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include  "RooGamma.h"
#include "RooDataHist.h"
#include "TH1D.h"
using namespace RooFit ;

void Bs_ctErr_fit(){

             TChain* chain_data = new TChain("fit");
             chain_data->Add("/Users/ab/Documents/BsToJpsiPhi/Angle_Eff/BsMC17_JpsiMu.root");
             Int_t nevt = (int)chain_data->GetEntries();
    std::cout<<"Events"<<nevt<<"\n";

              RooRealVar *BsCt2DMC = new RooRealVar("BsCt2DMC", "ct(cm)", 0.007, 0.3);
              RooRealVar *BsCt2DMCErr = new  RooRealVar("BsCt2DMCErr", " #Deltact_{B0}",0.0002, 0.005,"cm");
              RooRealVar *svmass = new RooRealVar("svmass", "mass", 5.1, 5.5);
              //Creating a data set which we are going to fit with the variables defined above
              RooDataSet* data = new RooDataSet("data", "raw data1", RooArgSet(*BsCt2DMC,*BsCt2DMCErr,*svmass),Import(*chain_data));
    
    
    RooRealVar *ml = new RooRealVar("ml","mean landau",0.001,0.0009,0.005) ;
    RooRealVar *sl = new RooRealVar("sl","sigma landau",0.0003,0.00001,0.003) ;
    //RooLandau *Di_Gamma = new RooLandau("lx","lx",*BsCt2DMCErr,*ml,*sl) ;
    
    
              RooRealVar *g1 =new RooRealVar ("g1","g1",13.638, 1.0, 20.0) ;
              RooRealVar *g2 = new RooRealVar("g2","g2",6.082,1.0,15.00) ;
              RooRealVar *b1 = new RooRealVar("b1","b1",0.00009512,0, 0.001) ;
              RooRealVar *b2 = new RooRealVar("b2","b2",0.0002825,0,0.001) ;
              RooRealVar *m1 = new RooRealVar("m1","m1", 0.00019) ;
              RooRealVar *m2 = new RooRealVar("m2","m2", 0.00019) ;
              
              RooGamma  *G1 = new RooGamma("G1", "G1", *BsCt2DMCErr, *g1, *b1, *m1);
              RooGamma  *G2 = new RooGamma("G2", "G2", *BsCt2DMCErr, *g2, *b2, *m2);
              
              RooRealVar *frac = new RooRealVar("frac","frac",0.5,0,1.0) ;
              
              RooAddPdf *Di_Gamma = new RooAddPdf("Di_Gamma", "Di_Gamma", RooArgList(*G1,*G2), RooArgList(*frac));
              
               RooFitResult* fitRes = Di_Gamma->fitTo(*data,Save());
               //RooFitResult* fitRes = modelEff.fitTo(*data,Save());
               //fitRes->Print("v");
               RooPlot* decayuncer = BsCt2DMCErr->frame(Title("Decay uncertainty"),Bins(100));
               data->plotOn(decayuncer,DataError(RooAbsData::SumW2));
               Di_Gamma->plotOn(decayuncer) ;
               Di_Gamma->paramOn(decayuncer);
               Di_Gamma->plotOn(decayuncer, LineColor(kBlack), LineWidth(1));
               Double_t chisquare_time = decayuncer->chiSquare();
               cout<<"Chi square of ctErr fit is :"<<chisquare_time<<"\n";
               RooPlot* pullframe2 = BsCt2DMCErr->frame(100);
               RooHist* hpull12 = decayuncer->pullHist();
               //TCanvas* c = new TCanvas("Lifetime Fit","Lifetime Fit",800, 800);
              Di_Gamma->plotOn(decayuncer, Components(*G1), LineColor(3), LineWidth(2), LineStyle(4));
              Di_Gamma->plotOn(decayuncer, Components(*G2), LineColor(2), LineWidth(2), LineStyle(2));
              TCanvas *cc2 = new TCanvas("cc2", "cc2",0,0,800,600);
               TPad *pad12 = new TPad("pad12","pad12",0,0.33,1,1);
               TPad *pad22 = new TPad("pad22","pad22",0,0,1,0.33);
               pad12->SetBottomMargin(0.00001);
               pad12->SetBorderMode(0);
               pad22->SetTopMargin(0.00001);
               pad22->SetBottomMargin(0.1);
               pad22->SetBorderMode(0);
               pad12->Draw();
               pad22->Draw();
               pad12->cd();
               gStyle->SetOptTitle(0);
               // cc->Range(-0.3786885,-19.31166,0.1745902,168.5381);
               cc2->SetFillColor(0);
               cc2->SetBorderSize(2);
               cc2->SetLeftMargin(0.1422222);
               cc2->SetRightMargin(0.04444445);
               decayuncer->SetStats(0);
               decayuncer->Draw();
               auto cms12 = new TLatex(0.00001, 3700, "#bf{CMS} #it{Simulations}, #sqrt{#bf{s}} = #bf{13 TeV}");
               cms12->SetNDC(false);
               cms12->SetTextColor(12);
               cms12->SetTextFont(42);
               cms12->SetTextSize(0.055);
               cms12-> Draw();        
               pad22->cd();
              pullframe2->addPlotable(hpull12,"P");
              pullframe2->Draw();
              cc2->cd();
              
  
}
