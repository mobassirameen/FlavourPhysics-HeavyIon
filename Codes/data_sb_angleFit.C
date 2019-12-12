#include "RooRealVar.h"
#include "RooDataSet.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TMath.h"
#include "RooPlot.h"
using namespace RooFit;
using namespace std;


void data_sb_angleFit(string filename, int maxOrder=7, int xbins=70, int ybins = 70, int zbins = 30)
{
  Int_t nbins =20;
  if ( maxOrder < 0 ) return;
  if ( ybins<1 ) ybins = xbins;
  if ( zbins<1 ) zbins = xbins;
  if ( xbins<1 ) return;

  RooRealVar *BscosthetaMC=new RooRealVar("BscosthetaMC","cos(#theta_{T})",-1,1);
  RooRealVar *BscospsiMC=new RooRealVar("BscospsiMC","cos(#psi_{T})",-1,1);
  RooRealVar *BsphiMC = new RooRealVar("BsphiMC","#phi",-TMath::Pi(),TMath::Pi());
  RooArgSet vars(*BscosthetaMC, *BscospsiMC, *BsphiMC);

   TChain* chain_data = new TChain("fit");
   chain_data->Add(filename.c_str());
   Int_t nevt = (int)chain_data->GetEntries();
    
   RooDataSet* data = new RooDataSet("data", "raw data1", RooArgSet(*BscosthetaMC, *BscospsiMC, *BsphiMC), Import(*chain_data));

  
    
    
    RooRealVar *psipar1 = new RooRealVar("psipar1", "psipar1", 0.5);
    RooRealVar *psipar11 = new RooRealVar("psipar11", "psipar11", 0.375);
    RooRealVar *psipar2 = new RooRealVar("psipar2", "psipar2", -1, 1);
    RooRealVar *psipar3 = new RooRealVar("psipar3", "psipar3", -1, 1);
    RooLegendre *psiL2 = new RooLegendre("psiL2", "psiL2",*BscospsiMC, 2,0);
    RooLegendre *psiL4 = new RooLegendre("psiL4", "psiL4",*BscospsiMC, 4,0);
    RooLegendre *thetaL2 = new RooLegendre("thetaL2", "thetaL2",*BscosthetaMC, 2,0);
    RooLegendre *thetaL4 = new RooLegendre("thetaL4", "thetaL4",*BscosthetaMC, 4,0);
    RooLegendre *thetaL6 = new RooLegendre("thetaL6", "thetaL6",*BscosthetaMC, 6,0);
    
    RooRealVar *thephipar1 = new RooRealVar("thephipar1", "thephipar1", -1, 1);
    RooRealVar *thephipar2 = new RooRealVar("thephipar2", "thephipar2", -1, 1);
    RooRealVar *thephipar3 = new RooRealVar("thephipar3", "thephipar3", -1, 1);
    RooRealVar *thephipar4 = new RooRealVar("thephipar4", "thephipar4", -1, 1);
    RooRealVar *thephipar5 = new RooRealVar("thephipar5", "thephipar5", -1, 1);
    RooRealVar *thephipar6 = new RooRealVar("thephipar6", "thephipar6", -1, 1);
    RooRealVar *thephipar7 = new RooRealVar("thephipar7", "thephipar7", -1, 1);
    RooRealVar *thephipar8 = new RooRealVar("thephipar8", "thephipar8", 0.375);
    
    RooFormulaVar phinew("phinew","BsphiMC/TMath::Pi()",*BsphiMC);
    RooLegendre *phizero =new RooLegendre("phizero","phizero",phinew,0,0);
    RooGenericPdf *psibgpdf = new RooGenericPdf("psibgpdf", "psibgpdf", "@0+@1*@3+@2*@4", RooArgList(*psipar1,*psipar2,*psipar3,*psiL2,*psiL4));
    RooGenericPdf *thetphiEq1 = new RooGenericPdf("thetphiEq1", "thetphiEq1", "@0+@1*@3+@2*@4", RooArgList(*psipar1,*thephipar1,*thephipar2,*thetaL2,*thetaL6));
    RooGenericPdf *thetphiEq2 = new RooGenericPdf("thetphiEq2", "thetphiEq2", "@0*cos(2*@3)+@1*cos(4*@3)+@2*cos(6*@3)", RooArgList(*thephipar3,*thephipar4,*thephipar5,*BsphiMC));
    RooGenericPdf *thetphiEq3 = new RooGenericPdf("thetphiEq3", "thetphiEq3", "@0*@1*(sin(@2)*sin(@2)-@3)", RooArgList(*thephipar6,*thetaL2,*BsphiMC,*psipar1));
    RooGenericPdf *thetphiEq4 = new RooGenericPdf("thetphiEq4", "thetphiEq4", "@0*@1*(sin(@2)*sin(@2)-@3)", RooArgList(*thephipar7,*thetaL4,*BsphiMC,*psipar1));
    RooGenericPdf *thetphiEq5 = new RooGenericPdf("thetphiEq5", "thetphiEq5", "@0*@1+(sin(@2)*sin(@2)*sin(@2)*sin(@2)-@3)", RooArgList(*thephipar8,*thetaL2,*BsphiMC,*psipar11));
    
    //RooRealSumPdf *TwoDimThetPhi = new RooRealSumPdf("TwoDimThetPhi", "ti", RooArgList(*thetphiEq1,*thetphiEq2, *thetphiEq3, *thetphiEq4, *thetphiEq5));
    
    RooArgList bgfuntot(*thetphiEq1,*thetphiEq2);
    bgfuntot.add(*thetphiEq3);
    bgfuntot.add(*thetphiEq4);
    bgfuntot.add(*thetphiEq5);
    RooArgList bgcoeftot(*thephipar1,*thephipar4);
    bgcoeftot.add(*thephipar6);
    bgcoeftot.add(*thephipar7);
    bgcoeftot.add(*thephipar8);
    
    RooRealSumPdf *TwoDimThetPhi=new RooRealSumPdf("TwoDimThetPhi","TwoDimThetPhi",bgfuntot,bgcoeftot);
    RooProdPdf *bg_AngularPdf = new RooProdPdf("bg_AngularPdf", "bg_AngularPdf", RooArgList(*psibgpdf,*TwoDimThetPhi));
    
   /* RooFitResult* bgfitresults = bg_AngularPdf->fitTo(*data, NumCPU(8), Save(kTRUE));
    bgfitresults->Print("v");
    
    gStyle->SetOptStat(0) ;
    gStyle->SetPalette(1) ;
    TH2* hcorr = bgfitresults->correlationHist() ;
    TCanvas* c = new TCanvas("Correlation Matrix","Correlation Matrix",1000,600) ;
    hcorr->GetYaxis()->SetTitleOffset(1.4) ; hcorr->Draw("colz");
*/
    
    
    RooPlot* costheta = BscosthetaMC->frame(Title("cos#theta_{T} (a.u.)"),Bins(nbins));
    data->plotOn(costheta,DataError(RooAbsData::SumW2));
    RooPlot* cospsi = BscospsiMC->frame(Title("cos#psi_{T} (a.u.)"),Bins(nbins));
    data->plotOn(cospsi,DataError(RooAbsData::SumW2));
    RooPlot* phi = BsphiMC->frame(Title("#phi_{T} (a.u.)"),Bins(nbins));
    data->plotOn(phi,DataError(RooAbsData::SumW2));
    RooPlot* costheta_func = BscosthetaMC->frame(Title("cos#theta_{T} (a.u.)"),Bins(nbins));
    bg_AngularPdf->plotOn(costheta_func,DataError(RooAbsData::SumW2));
    RooPlot* cospsi_func = BscospsiMC->frame(Title("cos#psi_{T} (a.u.)"),Bins(nbins));
    bg_AngularPdf->plotOn(cospsi_func,DataError(RooAbsData::SumW2));
    RooPlot* phi_func = BsphiMC->frame(Title("#phi_{T} (a.u.)"),Bins(nbins));
    bg_AngularPdf->plotOn(phi_func,DataError(RooAbsData::SumW2));
    TCanvas *lz = new TCanvas();
    lz->Divide(3,2);
    lz->cd(1);costheta->Draw();lz->cd(2);cospsi->Draw();lz->cd(3);phi->Draw();
    lz->cd(4);costheta_func->Draw();lz->cd(5);cospsi_func->Draw();lz->cd(6);phi_func->Draw();
    
    
    
    

}
