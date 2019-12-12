/// \author  - Muhammad Alibordi
// Test of RooKeyPDF has ability to discriminate the background and

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TH1D.h"
#include "TRandom.h"
#include "RooJohnsonLocal.cxx"
using namespace RooFit ;
using namespace std;

TH1* makeTH1()
{
    TFile *fIn1 = new TFile("/Users/ab/Documents/BsToJpsiPhi/FourTree/fittree_ntuBsMC2017.root");
    TTree* myFit = (TTree*)fIn1->Get("treeFit");
    Int_t n_entries = myFit->GetEntries();
    Float_t        mistaG;
    myFit->SetBranchAddress("mistag",&mistaG);
    TH1D *FlavMTag = new TH1D("FlavMTag","Bs tagged as BsBar", 100 ,0.04, 0.70);
    for (Int_t i=0;i<n_entries;i++) {
        myFit->GetEntry(i);
            FlavMTag->Fill(mistaG);
        }
    return FlavMTag;
}

void mass_mistag_bg()
{
     TChain* chain_data = new TChain("treeFit");
     chain_data->Add("/Users/ab/Documents/BsToJpsiPhi/FourTree/fittree_ntuBsMC2017.root");
     Int_t nevt = (int)chain_data->GetEntries();
    std::cout<<"Number of total events"<<nevt<<"\n";
    //Creating a data set which we are going to fit with the variables defined above
   
    RooRealVar *svmass= new RooRealVar("svmass", "M_{B_{s}} GeV/c^{2}",5.25,5.49);
    RooRealVar *mistag = new RooRealVar("mistag","Mistag fraction of original B and Bbar",0.04,0.70);
    
    TH1* hh = makeTH1() ;
    
     RooDataSet *data = new RooDataSet("data", "data", RooArgSet(*svmass, *mistag), Import(*chain_data));
     RooDataSet *dataPBG = new RooDataSet("dataPBG", "dataPBG", RooArgSet(*svmass, *mistag), Import(*chain_data));
    
    RooRealVar *mul = new RooRealVar("mul","mul",0.3);
    RooRealVar *sigma = new RooRealVar("sigma","sigma",0.4);
    RooRealVar *excon = new RooRealVar("excon","excon",0.3);
    RooExponential *p1 = new RooExponential("p1","p1",*svmass,*excon);
    RooGaussian *p2 = new RooGaussian("p2","p2",*mistag,*mul,*sigma);
    TF1 * massfunc =p1->asTF(RooArgList(*svmass));
    auto m =(TH1D*)massfunc->GetHistogram();
    TF1 * mtfunc =p2->asTF(RooArgList(*mistag));
    auto mt =(TH1D*)mtfunc->GetHistogram();
    TF1 *myfunc = mt->GetFunction("myfunc");
    TCanvas *ll = new TCanvas();
    ll->Divide(2,1);
    ll->cd(1);m->Draw();ll->cd(2);mt->Draw();
    
    RooProdPdf * bgpdfsample = new RooProdPdf("bgpdfsample", "bgpdfsample", RooArgList(*p1, *p2));
    RooDataSet* data_bkg = bgpdfsample->generate(RooArgSet(*svmass,*mistag),98895);
        
    dataPBG->append(*data_bkg);
    dataPBG->Print("v");
    
    
   
    /*RooPlot* bsmass = svmass->frame(Title("M_{B_{s}} (GeV/c^{2})"),Bins(40));
    data2->plotOn(bsmass,DataError(RooAbsData::SumW2));
    RooPlot* misTagpl = mistag->frame(Title("A-kernel-estimation for mistag w/o mirroring"),Bins(40));
    data2->plotOn(misTagpl);
    bkghistpdf->plotOn(misTagpl);
    TCanvas *lg = new TCanvas();
    lg->Divide(2,1);
    lg->cd(1);bsmass->Draw();lg->cd(2);misTagpl->Draw();
    */
   
   
    
    
   
   

   
   //==============RooHistPdf for mistag
    
     RooDataHist * bkgdatahist = new RooDataHist("bkgdatahist", "bkgdatahist", *mistag, *data_bkg,1);
    RooHistPdf * bkghistpdf = new RooHistPdf("bkghistpdf", "bkghistpdf", *mistag, *bkgdatahist,0);
    RooDataHist * sigdatahist = new RooDataHist("sigdatahist", "sigdatahist", *mistag, *data,1);
    RooHistPdf * sighistpdf = new RooHistPdf("sighistpdf", "sighistpdf", *mistag, *sigdatahist,0);
    
    
     TCanvas *lg = new TCanvas();
    lg->Divide(2,2);
     RooPlot* misfr1 = mistag->frame(Title("mistag and its histpdf"),Bins(40));
    data_bkg->plotOn(misfr1);   
    RooPlot* misfr2 = mistag->frame(Title("mistag and its histpdf"),Bins(40));
    data->plotOn(misfr2);
      RooPlot* misfr3 = mistag->frame(Title("mistag and its histpdf"),Bins(40));
    data_bkg->plotOn(misfr3);
    bkghistpdf->plotOn(misfr3);
      RooPlot* misfr4 = mistag->frame(Title("mistag and its histpdf"),Bins(40));
    data->plotOn(misfr4);
    sighistpdf->plotOn(misfr4);
    lg->cd(1);misfr1->Draw();lg->cd(2);misfr2->Draw(); lg->cd(3);misfr3->Draw();lg->cd(4);misfr4->Draw();
    
    
    RooDataHist *mistagSigData = new RooDataHist("mistagSigData","mistagSigData",*mistag,Import(*hh)) ;
    RooHistPdf *mtSigPdf = new RooHistPdf("mtSigPdf","mtSigPdf",*mistag, *mistagSigData, 0) ;
    RooDataHist *mistagBkgData = new RooDataHist("mistagBkgData","mistagBkgData",*mistag,Import(*mt)) ;
    RooHistPdf *mtBkgPdf = new RooHistPdf("mtBkgPdf","mtBkgPdf",*mistag, *mistagBkgData, 0) ;
    
    RooRealVar *conts = new RooRealVar("conts","conts", 0.1, 0.4);
    RooExponential *expoBg = new RooExponential("expoBg","expoBG",*svmass,*conts);

    
   



    RooRealVar *mean = new RooRealVar("mean","mean",5.36679, 5.35, 5.37) ;
    RooRealVar *sigma1 = new RooRealVar("sigma1","sigma1",0.03,0.,0.5) ;
    RooRealVar *sigma2 = new RooRealVar("sigma2","sigma2",0.018,0.,0.5) ;
    RooRealVar *sigma3 = new RooRealVar("sigma3","sigma3",0.022,0.,0.5) ;
   
     RooRealVar *mu= new RooRealVar("mu", "mu", 5.36679, 5.35, 5.37);
    RooRealVar *lambda = new RooRealVar("lambda", "lambda", 0.5, 0, 1);
    RooRealVar *gamma = new RooRealVar("gamma", "gamma", 0, 0.0, 0.01);
    RooRealVar *delta = new RooRealVar("delta", "delta", 1., 0, 10);
    
    
    
     RooJohnsonLocal *john = new RooJohnsonLocal("john", "john", *svmass, *mu, *lambda, *gamma, *delta);
//===============================================================================================
    RooGaussian *gauss1 = new RooGaussian("gauss1","gauss1",*svmass,*mean,*sigma1) ;
    RooGaussian *gauss2 = new RooGaussian("gauss2","gauss2",*svmass,*mean,*sigma2) ;
    RooGaussian *gauss3 = new RooGaussian("gauss3","gauss3",*svmass,*mean,*sigma3) ;
   

    RooRealVar *frac1 = new RooRealVar("frac1","frac1",0.1681,0.1, 0.4);
    RooRealVar *frac2 = new RooRealVar("frac2","frac2",0.3,0.1,0.60);
    RooRealVar *frac3 = new RooRealVar("frac3","frac3",0.2,0.1,1);
   
  
     // Create adaptive kernel estimation pdf. In this configuration the input data, is mirrored over the boundaries to minimize edge effects in distribution, that do not fall to zero towards the edges,An adaptive kernel estimation pdf on the same data without mirroring option
    //===========================================================================================================kernel estimation pdf
    //RooKeysPdf *kerestisig= new RooKeysPdf("kerestisig","kerestisig",*mistag,*data,RooKeysPdf::NoMirror);
    //RooKeysPdf *kerestiMTBG= new RooKeysPdf("kerestiMTBG","kerestiMBG",*mistag,*data_bkg,RooKeysPdf::NoMirror);
    
    //RooKeysPdf *kerestiMBG= new RooKeysPdf("kerestiMBG","kerestiMBG",*svmass,*data2,RooKeysPdf::NoMirror);
    
    
    RooAddPdf *T_gaus = new RooAddPdf("T_gaus","T_gaus",RooArgList(*gauss1,*gauss2,*gauss3),RooArgList(*frac1,*frac2));

    RooRealVar *nSig = new RooRealVar("nSig", "Number of Signal Events in SIGNAL MC",1400,0.,(int)chain_data->GetEntries());
    RooRealVar *nBkg = new RooRealVar("nBkg", "Number of Backgound Events in produced  MC",800,0.,(int)chain_data->GetEntries());
   
    RooProdPdf *sigpdf = new RooProdPdf("sigpdf", "mass*mistag",RooArgList(*john,*sighistpdf));
   //RooProdPdf *sigpdf = new RooProdPdf("sigpdf", "mass*mistag",RooArgList(*john,*kerestisig));
    //RooAddPdf *SIG_pdf = new RooAddPdf("SIG_pdf", "Total pdf",RooArgList(*sigpdf), RooArgList(*nSig));
    
    RooProdPdf *bkgpdf = new RooProdPdf("bkgpdf", "mass*mistag",RooArgList(*expoBg,*bkghistpdf));
  // RooProdPdf *bkgpdf = new RooProdPdf("bkgpdf", "mass*mistag",RooArgList(*expoBg,*kerestiMTBG));
    //RooAddPdf *BKG_pdf = new RooAddPdf("BKG_pdf", "Total pdf",RooArgList(*bkgpdf), RooArgList(*nBkg));
    
    

    RooAddPdf *MTpdf = new RooAddPdf("MTpdf","MTpdf",RooArgList(*sigpdf,*bkgpdf),RooArgList(*frac3));
    
    
    RooFitResult* fitRes = MTpdf->fitTo(*dataPBG,Save());
    fitRes->Print("v");
    
    RooPlot* bsmass = svmass->frame(Title("M_{B_{s}} (GeV/c^{2})"),Bins(100));
    dataPBG->plotOn(bsmass,DataError(RooAbsData::SumW2));
    MTpdf->plotOn(bsmass) ;
    MTpdf->paramOn(bsmass);
    //MTpdf->plotOn(bsmass, LineColor(kBlue), LineWidth(1));
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
    
    
    MTpdf->plotOn(bsmass, Components(*john), LineColor(3), LineWidth(1), LineStyle(4));
    //MTpdf->plotOn(bsmass, Components(*gauss2), LineColor(2), LineWidth(1), LineStyle(3));
    //MTpdf->plotOn(bsmass, Components(*gauss3), LineColor(6), LineWidth(1), LineStyle(5));
    MTpdf->plotOn(bsmass, Components(*expoBg), LineColor(46), LineWidth(1), LineStyle(6));
    
    RooPlot* misTagpl = mistag->frame(Title("A-kernel-estimation for mistag w/o mirroring"),Bins(100));
    dataPBG->plotOn(misTagpl);
    MTpdf->paramOn(misTagpl);
    //MTpdf->plotOn(misTagpl, Components(*kerestisig), LineColor(3), LineWidth(1), LineStyle(4));
    //MTpdf->plotOn(misTagpl, Components(*kerestiMTBG), LineColor(2), LineWidth(1), LineStyle(5));
    
   // MTpdf->plotOn(misTagpl, Components(*mtSigPdf), LineColor(3), LineWidth(1), LineStyle(4));
   // MTpdf->plotOn(misTagpl, Components(*mtBkgPdf), LineColor(2), LineWidth(1), LineStyle(5));
    MTpdf->plotOn(misTagpl, Components(*sighistpdf), LineColor(3), LineWidth(2), LineStyle(4));
    MTpdf->plotOn(misTagpl, Components(*bkghistpdf), LineColor(2), LineWidth(2), LineStyle(5));
    
    TCanvas *c = new TCanvas("c", "c",0,0,600,600);
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
    auto cms1 = new TLatex(5.25, 32000, "#bf{CMS} #it{Simulations}, #sqrt{#bf{s}} = #bf{13 TeV}");
    cms1->SetNDC(false);
    cms1->SetTextColor(12);
    cms1->SetTextFont(42);
    cms1->SetTextSize(0.055);
    cms1-> Draw();
    pad2->cd();
    pullframe->SetStats(0);
    pullframe->Draw();
    c->cd();
    
    TCanvas* mis = new TCanvas("mis","kernelestimation",600,600) ;
    gPad->SetLeftMargin(0.15) ; misTagpl->GetYaxis()->SetTitleOffset(1.4) ;misTagpl->Draw() ;
 
 
 /*
    Float_t MisTag, MassBs;
     chain_data->SetBranchAddress("svmass",&MassBs);
     chain_data->SetBranchAddress("mistag",&MisTag);
     TRandom * mom =new TRandom();
     Double_t p0=-0.083, p00=0.0765;
     Double_t p0m=-0.03, p00m=0.0065;
     TH1D *m=new TH1D("m","m",100, 5.28,5.45); TH1D* mt= new TH1D("mt","mt",100,0.0,1.0);
     TH1D *m1=new TH1D("m1","m1",100, 5.28,5.45); TH1D* mt1= new TH1D("mt1","mt1",100,0.0,1.0);
     for(int l =0 ; l<(nevt/4);++l)
     {   chain_data->GetEntry(l);
     Double_t ms = mom->Uniform(5.28, 5.45);m->Fill(ms);
     Double_t tg = mom->Uniform(0.0,1.0);mt->Fill(MisTag);
     }
     TCanvas * mplt = new TCanvas();
     mplt->Divide(2,1);mplt->cd(1);m->Draw();mplt->cd(2);mt->Draw();
     for(int i =0 ; i<nevt;++i){
     chain_data->GetEntry(i);
     //svmass->setVal(MassBs+p0*(mom->Binomial(nevt,0.7))+p00*(mom->Binomial(nevt,0.5))*(m->GetRandom()));
     //mistag->setVal(MisTag+p0*(mom->Binomial(nevt,0.6))+p00*(mom->Binomial(nevt,0.7))*(mt->GetRandom()));
     svmass->setVal(MassBs+p0m+p00m*(m->GetRandom()));
     mistag->setVal(MisTag+p0+p00*(mt->GetRandom()));
     m1->Fill(MassBs+p0m+p00m*(m->GetRandom()));
     mt1->Fill(MisTag+p0+p00*(mt->GetRandom()));
     dataPBG->add(RooArgSet(*svmass, *mistag));
     }
     TCanvas * mplt1 = new TCanvas();
     mplt1->Divide(2,1);
     mplt1->cd(1);m1->Draw();mplt1->cd(2);mt1->Draw();
     dataPBG->Print("v");
  */
    
    
    
    
}
  










    
