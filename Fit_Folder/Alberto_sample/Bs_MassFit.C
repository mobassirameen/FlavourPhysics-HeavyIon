/// \author 07/2018 - Muhammad Alibordi 


#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TH1D.h"
#include "TRandom.h"
using namespace RooFit ;

TH1* makeTH1() ;
TTree* makeTTree() ;
void Bs_MassFit()
{
     TChain* chain_data = new TChain("treeBs");
     chain_data->Add("/Users/ab/Documents/BsToJpsiPhi/Fit_Folder/Alberto_sample/fittree_reco_v2.root");
     //chanin_data->Add("/home/cms/BMeson/FileRepository/Generalplots/BsJpsiPhi_dG0.root");
     Int_t nevt = (int)chain_data->GetEntries();
   
    //Creating a data set which we are going to fit with the variables defined above
   
    RooRealVar *svmass= new RooRealVar("svmass", "M_{B_{s}} GeV/c^{2}",5.28,5.45);
   
   RooDataSet* data = new RooDataSet("data", "raw data1", RooArgSet(*svmass),Import(*chain_data));
    RooPlot* bsmass1 = svmass->frame(Title("M_{B_{s}} (GeV/c^{2})"),Bins(100));
    data->plotOn(bsmass1,DataError(RooAbsData::SumW2));
    TCanvas *c1 = new TCanvas("c1", "c1",0,0,1000,600);
    bsmass1->Draw();
   
 
   

    RooRealVar *mean = new RooRealVar("mean","mean",5.36,5.28,5.4) ;
    RooRealVar *sigma1 = new RooRealVar("sigma1","sigma1",0.03,0.,0.5) ;
    RooRealVar *sigma2 = new RooRealVar("sigma2","sigma2",0.018,0.,0.5) ;
    RooRealVar *sigma3 = new RooRealVar("sigma3","sigma3",0.022,0.,0.5) ;
   
//===============================================================================================
    RooGaussian *gauss1 = new RooGaussian("gauss1","gauss1",*svmass,*mean,*sigma1) ;
    RooGaussian *gauss2 = new RooGaussian("gauss2","gauss2",*svmass,*mean,*sigma2) ;
    RooGaussian *gauss3 = new RooGaussian("gauss3","gauss3",*svmass,*mean,*sigma3) ;
   
    
//========================================================================================================== Construct decay(t) X Gaaussian resoltution fn
  
  


    RooRealVar *nsmass = new RooRealVar("nsmass","number signal events",14,0.,1000000);
    RooRealVar *nstime = new RooRealVar("nstime","number signal events",14,0.,1000000);

    RooRealVar *frac1 = new RooRealVar("frac1","frac1",0.1681);
    RooRealVar *frac2 = new RooRealVar("frac2","frac2",0.45,0.35,.60);
 

     
    RooAddPdf *T_gaus = new RooAddPdf("T_gaus","T_gaus",RooArgList(*gauss1,*gauss2,*gauss3),RooArgList(*frac1, *frac2));
    



    RooRealVar *nSig = new RooRealVar("nSig", "Number of Signal Events in SIGNAL MC",1400,0.,(int)chain_data->GetEntries());
   

    RooAddPdf *MTpdf = new RooAddPdf("MTpdf", "Total pdf",RooArgList(*T_gaus), RooArgList(*nSig));
    RooFitResult* fitRes = MTpdf->fitTo(*data,Save(),Extended(1));
    fitRes->Print("v");
    RooPlot* bsmass = svmass->frame(Title("M_{B_{s}} (GeV/c^{2})"),Bins(100));
   
     

    data->plotOn(bsmass,DataError(RooAbsData::SumW2));
    MTpdf->plotOn(bsmass) ;
    MTpdf->paramOn(bsmass);
    MTpdf->plotOn(bsmass, LineColor(kBlue), LineWidth(2));
    Double_t chisquare_mass = bsmass->chiSquare();
    cout<<"Chi square of mass fit is :"<< chisquare_mass<< endl;
    RooPlot* pullframe = svmass->frame(RooFit::Title("Mass pull"));
    RooHist* hpull1 = bsmass->pullHist();
    //pullframe->addPlotable(hpull1,"P");
    pullframe->addPlotable(hpull1,"P0") ;
    pullframe->SetMinimum(-3) ;
    pullframe->SetMaximum(+3) ;
    pullframe->SetYTitle("pull");
    pullframe->SetMarkerStyle(20);
    pullframe->SetNdivisions(10);
    

    MTpdf->plotOn(bsmass, Components(*gauss1), LineColor(3), LineWidth(4), LineStyle(4));
    MTpdf->plotOn(bsmass, Components(*gauss2), LineColor(2), LineWidth(4), LineStyle(2));
    MTpdf->plotOn(bsmass, Components(*gauss3), LineColor(6), LineWidth(4), LineStyle(2));
    

    TCanvas *c = new TCanvas("c", "c",0,0,1000,600);
    TPad *pad1 = new TPad("pad1","pad1",0,0.33,1,1);
    TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.33);
    pad1->SetBottomMargin(0.00001);
    pad1->SetBorderMode(0);
    pad2->SetTopMargin(0.00001);
    pad2->SetBottomMargin(0.1);
    pad2->SetBorderMode(0);
    pad1->Draw();
    pad2->Draw();
    pad1->cd();
    gStyle->SetOptTitle(0);
    // cc->Range(-0.3786885,-19.31166,0.1745902,168.5381);
    c->SetFillColor(0);
    c->SetBorderSize(2);
    c->SetLeftMargin(0.1422222);
    c->SetRightMargin(0.04444445);
    bsmass->SetStats(0);
    bsmass->Draw();
    auto cms1 = new TLatex(5.28, 23000, "#bf{CMS} #it{Simulations}, #sqrt{#bf{s}} = #bf{13 TeV}");
    cms1->SetNDC(false);
    cms1->SetTextColor(12);
    cms1->SetTextFont(42);
    cms1->SetTextSize(0.055);
    cms1-> Draw();
    pad2->cd();
    pullframe->SetStats(0);
    pullframe->Draw();
    c->cd();
    
    
    
   
    
   
    
    
    
    
   

}
  










    
