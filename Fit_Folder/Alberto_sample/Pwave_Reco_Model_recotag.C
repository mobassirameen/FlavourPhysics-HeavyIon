//Author Md Alibordi, Giacomo Fedi

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooWorkspace.h"
#include "RooCategory.h"
#include "TKey.h"
#include "RooExponential.h"
#include <map>
#include "TCut.h"
#include "RooHistPdf.h"
#include "RooHist.h"
#include "TVirtualPad.h"
#include "RooDataHist.h"
#include <string>
#include "TEventList.h"
#include "RooFit.h"
#include "RooRealVar.h"
#include "TFile.h"
#include "RooDataSet.h"
#include "TTree.h"
#include "TH2D.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooGaussian.h"
#include "RooGenericPdf.h"
#include "RooFitResult.h"
#include "RooAddPdf.h"
#include "RooAddition.h"
#include "RooAbsPdf.h"
#include "RooPolynomial.h"
#include "RooGaussModel.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooDecay.h"
#include "RooDataHist.h"
#include "RooAddModel.h"
#include "RooProdPdf.h"
#include "RooProduct.h"
#include "RooPlot.h"
#include "TH1D.h"
#include "TRandom.h"
#include "RooMinuit.h"
#include "RooExtendPdf.h"
#include "RooChi2Var.h"
#include "Math/Functor.h"
#include "TRandom3.h"
#include "Math/DistFunc.h"
#include "RooClassFactory.h"
#include "RooFitResult.h"
#include "RooDataSet.h"
#include "RooRealConstant.h"
#include "RooConstVar.h"
#include "Roo1DTable.h"
#include "RooBDecay.h"
#include "RooFormulaVar.h"
//#include "RoogmModel.h"
#include "RooRealSumPdf.h"
#include "Math/SpecFunc.h"
#include "RooBMixDecay.h"
#include "RooBCPEffDecay.h"
#include "Riostream.h"
#include "RooRandom.h"
#include "TMath.h"
/*#include "RooFun1TRPdf.h"
#include "RooFun2TRPdf.h"
#include "RooFun3TRPdf.h"
#include "RooFun4TRPdf.h"
#include "RooFun5TRPdf.h"
#include "RooFun6TRPdf.h"
#include "RooEffcthPdf.h"*/
#include "RooMCStudy.h"
#include "RooArgSet.h"
#include "RooLegendre.h"
#include "RooSpHarmonic.h"
#include "RooBifurGauss.h"

using namespace RooFit; 
using namespace std;




TH1* makeTH1() 
{
                                         auto hctaucut_eff = new TH1D("hctaucut_eff","ct(proper time) Efficiency Plot; ct (cm);#epsilon (ct) (a.u.) ", 100,0.007, 0.3);
                                         auto hctaunocut_eff = new TH1D("hctaunocut_eff","ct(proper time) Efficiency Plot  ;ct[cm];#epsilon [a.u.]", 100 ,0.007, 0.3);
                                         auto resobulk_eff = new TH1D("resobulk_eff", "Resobulk ; resolution; Events ",  100, -0.004, 0.004);
                                         TFile *fIn1 = new TFile("/afs/cern.ch/work/m/mumuhamm/private/CMSSW_9_4_0/src/pwave_jpsimu/RECOFITRECOTAG/fittree_reco_v2.root");
                                         TTree* mumu = (TTree*)fIn1->Get("treeBs");
                                         Int_t n_entries = mumu->GetEntries();
                                         std::cout<<n_entries<<"\n";
                                              Double_t   recolife, gencutlife;
                               
                                               mumu->SetBranchAddress("BsCt2DMC",&recolife);              
                                               mumu->SetBranchAddress("BsCt2DMC_GEN",&gencutlife);       
                                                     
                                               for (Int_t i=0;i<n_entries;i++) { 
                                                                                                      mumu->GetEntry(i);
                                                                                                      hctaucut_eff->Fill(recolife);
                                                                                                      resobulk_eff->Fill(gencutlife-recolife);                                                                                                                 
                                                                                               }
                                                
     TFile *fIn2 = new TFile("/afs/cern.ch/work/m/mumuhamm/private/CMSSW_9_4_0/src/BsGenV13.root");
     TTree* ffgen = (TTree*)fIn2->Get("OutTreeGEN");
     Int_t ff_entries = ffgen->GetEntries();
     Float_t   ctaugen_eff;
     ffgen->SetBranchAddress("ctau_GEN",&ctaugen_eff);
     for (Int_t j=0;j<ff_entries;j++) 
                                                        {
                                                                  ffgen->GetEntry(j);      
                                                                  if(ctaugen_eff != ctaugen_eff) continue;           
                                                                  Double_t actualctgen_eff =  ctaugen_eff +  resobulk_eff->GetRandom();                     
                                                                  hctaunocut_eff->Fill(actualctgen_eff);
                                                        }
                                                        
                             
                             auto effgraph = (TH1D*)hctaucut_eff->Clone("effgraph");
                             effgraph->Sumw2();
                             effgraph->Divide(hctaunocut_eff);               
                             return effgraph;                 
                                              
}



Double_t ctFunction(Double_t *x,Double_t *par)
   {
      Double_t func;
        Double_t expoinv;
        Double_t lin, bolic;
	Double_t xx = x[0];
        lin         = par[0] + par[3]*xx;
        expoinv     = par[1]*(1./(1.+exp(-xx/par[2])));
        bolic            = par[3]*xx*xx;
        Double_t  poly2          = par[0] + par[1]*xx + par[2]*xx*xx + par[3]*xx*xx*xx + par[4]*xx*xx*xx*xx ;
        Double_t  sigmoid     = par[5]*(1./(1.+exp(-xx/par[6])));
        func            = poly2+sigmoid; //-bolic;
        return func;
    } 





void Pwave_Reco_Model_recotag(string filename, int maxOrder=7, int xbins=25, int ybins = 0, int zbins = 0){
//==========================================================================================================Efficiencies

 
    TH1* hh = makeTH1() ;
    RooRealVar *effvar=new RooRealVar("effvar","ct_{reco}/ct_{gen-smeared}; ct (cm); #epsilon (a.u.)",0.007,0.3) ;
    RooDataHist *dh = new RooDataHist("dh","dh",*effvar,Import(*hh)) ;
    RooHistPdf *effvarpdf = new RooHistPdf("effvarpdf","effhistpdf",*effvar, *dh, 0) ;
    dh->Print("v") ;
    RooPolynomial *p =new RooPolynomial("p","p",*effvar,RooArgList(RooConst(10),RooConst(3.0),RooConst(0.0003),RooConst(-9.0) )) ;
    effvar->setBins(50) ;
    //RooDataSet *dataeff1 = p->generate(effvar , 1000000) ;
    //RooDataHist *dh1 = dataeff1->binnedClone() ;
    //RooHistPdf *effvarpolypdf = new RooHistPdf("effvarpdf","effhistpdf",*effvar, *dh1, 0) ;
    RooPlot *efficiency = effvar->frame(Title("Efficiency"),Bins(100)) ;
    dh->plotOn(efficiency) ;
    effvarpdf->plotOn(efficiency) ;
    auto ceff = new TCanvas("ceff", "ceff", 0, 0, 800, 600);
    efficiency->Draw();


 
 	
 	 auto ctform = new TF1("ctform",ctFunction,0.007,0.3,7);
         ctform->SetParameter(0,2.87);  // par[0] : offset
         ctform->SetParLimits(0,0.,5.0);
  	 ctform->SetParName(0,"p0");
 	 ctform->SetParameter(1,1.5); // par[1] : width at f(x) minimum
 	 ctform->SetParLimits(1,0.0,5.);
 	 ctform->SetParName(1,"p1");
 	 ctform->SetParameter(2,-136.48);  // par[2] : slope (negative side)
 	 ctform->SetParLimits(2,-200.0,0.0);
 	 ctform->SetParName(2,"p2");
 	 ctform->SetParameter(3,755.0);   // par[3] : inflection point abscissa
 	 ctform->SetParLimits(3,0.0,1000.0);
 	 ctform->SetParName(3,"p3");
 	 ctform->SetParameter(4,-1258.6);   // par[3] : inflection point abscissa
 	 ctform->SetParLimits(4,-2500.0,0.0);
 	 ctform->SetParName(4,"p4");	  
 	 ctform->SetParameter(5,3.); // par[1] : width at f(x) minimum
 	 ctform->SetParLimits(5,0.0,5.);
 	 ctform->SetParName(5,"p5");
 	 ctform->SetParameter(6,1);  // par[2] : slope (negative side)
 	 ctform->SetParLimits(6,-2.0,2.);
 	 ctform->SetParName(6,"p6");
 	 
  Int_t nbins;
  std::cout<<" Give the value of bin number "<<"\n";
  cin>>nbins;
  
   Double_t ctpar0,ctpar1,ctpar2,ctpar3,ctpar4, ctpar5, ctpar6;
 	  Double_t kappa_val;
 	  
 	Double_t bins[] = {0.007, 0.008, 0.009, 0.01,0.011,0.012, 0.013,0.014,0.015,0.016,0.017,0.018,0.019, 0.02, 0.021, 0.022, 0.023, 0.024, 0.025, 0.026, 0.027, 0.028, 0.03, 0.034, 0.038, 0.042, 0.046, 0.05, 0.054, 0.058,0.062, 0.066, 0.07, 0.082, 0.086,0.09,0.12, 0.15,0.18,0.21,0.24,0.27,0.3};
          Int_t  binnum = sizeof(bins)/sizeof(Double_t) - 1; 
          
          
   
      auto hctaucut = new TH1D("hctaucut","ct(cm);ct(cm);#epsilon (ct) (a.u.)", binnum, bins);
      auto hctaunocut = new TH1D("hctaunocut","ct(cm);ct(cm);#epsilon(ct) (a.u.)", binnum, bins);
      auto resobulk = new TH1D("resobulk", "Resobulk ;Resolution ; Events ",  nbins, -0.004, 0.004);
      auto pullbulk = new TH1D("pullbulk", "pullbulk ; ; Events ",  nbins, -5.0, 5.0); 
     
   TTree          *fChain, *copy;   //!pointer to the analyzed TTree or TChain
   TFile* f = new TFile(filename.c_str());
   fChain = new TTree;
   fChain=(TTree*)f->Get("treeBs");   
   if (!fChain){
    cout << "No TTree found in input file, returning" << endl;
    return;
  }



  Long64_t n_entries = fChain->GetEntries();
        cout << "Start Processing " << n_entries << " events" <<"\n";     
        Double_t ctreco, ctgencut, cterr;
       
        fChain->SetBranchAddress("BsCt2DMC",&ctreco);
        fChain->SetBranchAddress("BsCt2DMC_GEN",&ctgencut);
        fChain->SetBranchAddress("BsCt2DMCErr",&cterr);
        for (Int_t i=0;i<n_entries;i++) { 
                                                           fChain->GetEntry(i);
                                                           hctaucut->Fill(ctreco);
                                                           resobulk->Fill(ctgencut-ctreco);
                                                           pullbulk->Fill((ctgencut-ctreco)/cterr);
                                                      }
        
  
    TFile *fIn2 = new TFile("/afs/cern.ch/work/m/mumuhamm/private/CMSSW_9_4_0/src/BsGenV13.root");
    if (!fIn2){return;}
    TTree* ffgen = (TTree*)fIn2->Get("OutTreeGEN");
    Int_t ff_entries = ffgen->GetEntries();
    Float_t  ctaugen;
 
    ffgen->SetBranchAddress("ctau_GEN",&ctaugen);
    for (Int_t j=0;j<ff_entries;j++) {
        
        ffgen->GetEntry(j);
   
        Double_t actualctgen =  ctaugen +  resobulk->GetRandom();
        hctaunocut->Fill(actualctgen);
        
    }
    
    TCanvas *ckap= new TCanvas("ckap", "ckap",0,0,800,600);
    pullbulk->Draw();
    kappa_val = pullbulk->GetRMS();
    std::cout<<"Kappa_value"<<kappa_val<<"\n";
    
    TCanvas *creso= new TCanvas("creso", "creso",0,0,800,600);
    resobulk->Draw();
  
    
    
    
    auto c4 = new TCanvas("c4", "c4", 0, 0, 1000, 800);
    auto hdivideIV = (TH1D*)hctaucut->Clone("hdivideIV");
    hdivideIV->Sumw2();
    hdivideIV->Divide(hctaunocut);
    hdivideIV->Draw("colz");    //colz ep
    hdivideIV->Fit(ctform);
    
   
    ctpar0 = ctform->GetParameter(0);
    ctpar1 = ctform->GetParameter(1);
    ctpar2 = ctform->GetParameter(2);
    ctpar3 = ctform->GetParameter(3);
    ctpar4 = ctform->GetParameter(4);
    ctpar5 = ctform->GetParameter(5);
    ctpar6 = ctform->GetParameter(6);
      
    
   RooRealVar *ctp0 = new RooRealVar("ctp0","ctp0", ctpar0);
   RooRealVar *ctp1 = new RooRealVar("ctp1","ctp1", ctpar1);
   RooRealVar *ctp2 = new RooRealVar("ctp2","ctp2", ctpar2);
   RooRealVar *ctp3 = new RooRealVar("ctp3","ctp3", ctpar3);
   RooRealVar *ctp4 = new RooRealVar("ctp4","ctp4", ctpar4);  
   RooRealVar *ctp5 = new RooRealVar("ctp5","ctp5", ctpar5);  
   RooRealVar *ctp6 = new RooRealVar("ctp6","ctp6", ctpar6);  
   
   
   
     
  if ( maxOrder < 0 ) return;
  if ( ybins<1 ) ybins = xbins;
  if ( zbins<1 ) zbins = xbins;
  if ( xbins<1 ) return;


             TChain* chain_data = new TChain("treeBs");
             chain_data->Add(filename.c_str());
             Int_t nevt = (int)chain_data->GetEntries();



    RooRealVar *BsCt2DMC = new RooRealVar ("BsCt2DMC","Bs ct", 0.007,0.3,"cm");//
    RooRealVar *BsCt2DMCErr = new RooRealVar ("BsCt2DMCErr","Bs ct Err", 0.0002,0.005,"cm");//
    RooRealVar *BscosthetaMC= new RooRealVar ("BscosthetaMC","cos(#theta_{T})", -1,1);
    RooRealVar *BscospsiMC= new RooRealVar ("BscospsiMC","cos(#psi_{T})", -1,1);
    RooRealVar *BsphiMC= new RooRealVar ("BsphiMC","#phi_{T}", -TMath::Pi(),TMath::Pi(),"rad");//cosdelta2
    RooRealVar *svmass = new RooRealVar("svmass", "M_{B_{s}} GeV/c^{2}", 5.28, 5.45);
    RooArgSet vars(*BscosthetaMC, *BscospsiMC, *BsphiMC);
    double tagmis = 0;
    RooRealVar *mistag = new RooRealVar("mistag","Mistag fraction of original B and Bbar",0,1.0);
    //RooRealVar *tag = new RooRealVar("tag", "tag for original B and Bbar", -1,+1);  
    RooCategory *tag = new RooCategory("tag","Flavour tag of the B meson");
    tag->defineType("Bsbar",-1);
    tag->defineType("Bs",+1);
    tag->defineType("untag",0);
    //RooRealVar *tag = new RooRealVar("tag","tag",-2,2);
    RooDataSet* data = new RooDataSet("data", "raw data1", RooArgSet(*BscosthetaMC, *BscospsiMC, *BsphiMC, *BsCt2DMC, *BsCt2DMCErr, *svmass, *tag, *mistag), Import(*chain_data));
    Roo1DTable* btable = data->table(*tag) ;
   btable->Print() ;
   btable->Print("v") ;
   
    
    
  TChain* t_den = new TChain();
  TChain* t_num = new TChain();
  t_den->Add("/afs/cern.ch/work/m/mumuhamm/private/CMSSW_9_4_0/src/BsGenV13.root/OutTreeGEN");
  t_num->Add("/afs/cern.ch/work/m/mumuhamm/private/CMSSW_9_4_0/src/pwave_jpsimu/RECOFITRECOTAG/fittree_reco_v2.root/treeBs");
  Int_t denEntries = t_den->GetEntries();
  Int_t numEntries = t_num->GetEntries();
  Int_t counter;

  Float_t costhetagen, cospsigen, phigen;
  Double_t recoCosTheta, recoCosPsi, recoPhi;
  t_den->SetBranchAddress( "angle_costheta_GEN" , &costhetagen );
  t_den->SetBranchAddress( "angle_cospsi_GEN" , &cospsigen);
  t_den->SetBranchAddress( "angle_phi_GEN", &phigen);
  t_num->SetBranchAddress( "BscosthetaMC" , &recoCosTheta );
  t_num->SetBranchAddress( "BscospsiMC" , &recoCosPsi );
  t_num->SetBranchAddress( "BsphiMC", &recoPhi );

  RooDataSet* recodata    = new RooDataSet( "recodata"   , "GEN distribution", vars );
  RooDataSet* numData = new RooDataSet( "numData", "RECO distribution after selections", vars ); 
    counter=0;
  for (int iCand=0; iCand<denEntries; ++iCand) {
    t_den->GetEntry(iCand);
    BscosthetaMC->setVal(costhetagen);
    BscospsiMC->setVal(cospsigen);
    BsphiMC->setVal(phigen);
    recodata->add(vars);
  }
  counter =0;
  for (int iCand=0; iCand<numEntries; ++iCand) {
    t_num->GetEntry(iCand);
     if ( iCand > 1.0*counter*numEntries/100 ) { cout<<counter<<"%"<<endl; counter += 10;     }
   BscosthetaMC->setVal(recoCosTheta);
    BscospsiMC->setVal(recoCosPsi);
    BsphiMC->setVal(recoPhi);
    numData->add(vars);
  }
  
   
  // compute and print average efficiency
  double avgEff = numData->sumEntries()/recodata->sumEntries();
  cout<<"Average efficiency = "<<avgEff<<endl;
  
  
   TH3D* denHist = (TH3D*)recodata->createHistogram( "denHist", *BscosthetaMC,   Binning(xbins,-1,1) , YVar(*BscospsiMC,Binning(ybins,-1,1)), ZVar(*BsphiMC,Binning(zbins,-TMath::Pi(),TMath::Pi())) );
   TH3D* numHist = (TH3D*)numData->createHistogram("numHist", *BscosthetaMC,  Binning(xbins,-1,1) ,   YVar(*BscospsiMC,Binning(ybins,-1,1)),   ZVar(*BsphiMC,Binning(zbins,-TMath::Pi(),TMath::Pi())) );
 
   TEfficiency* Effthetapsi; TEfficiency* Effthetaphi ; TEfficiency* Effpsiphi;
   auto num_xy = numHist->Project3D("xy");
   auto num_xz = numHist->Project3D("xz");
   auto num_yz = numHist->Project3D("yz");
 
  auto den_xy = denHist->Project3D("xy");
  auto den_xz = denHist->Project3D("xz");
  auto den_yz = denHist->Project3D("yz");
 
 auto c1 = new TCanvas("c1", "c1", 0, 0, 1200, 800);
                     TH3D *hdivideIII = (TH3D*)numHist->Clone("hdivideIII");
                     hdivideIII->Sumw2();
                     hdivideIII->Divide(denHist);
                     auto proj2d1 =  hdivideIII->Project3D("xy");
                     auto proj2d2 =  hdivideIII->Project3D("xz");
                     auto proj2d3 =  hdivideIII->Project3D("yz");
                     proj2d1->SetTitle("Efficiency cos(#theta_{T})/cos(#psi_{T}) projection ");
                     proj2d2->SetTitle("Efficiency cos(#theta_{T})/#phi projection");
                     proj2d3->SetTitle("Efficiency cos(#psi_{T})/#phi projection");
                     c1->Divide(3,1);
                     c1->cd(1); proj2d1->GetXaxis()->SetTitleOffset(1.4);proj2d1->GetYaxis()->SetTitleOffset(2);proj2d1->SetMinimum(0.0);proj2d1->Draw("SURF3");
                     c1->cd(2);proj2d2->GetXaxis()->SetTitleOffset(1.4);proj2d2->GetYaxis()->SetTitleOffset(2);proj2d2->SetMinimum(0.0);proj2d2->Draw("SURF3");
                     c1->cd(3);proj2d3->GetXaxis()->SetTitleOffset(1.4);proj2d3->GetYaxis()->SetTitleOffset(2);proj2d3->SetMinimum(0.0);proj2d3->Draw("SURF3");

 
 
   vector < RooRealVar* > factors;
  vector < double > proj;
  vector < RooLegendre* > vectFuncLegCosThetaT;
  vector < RooLegendre* > vectFuncLegCosPsiT;
  vector < RooFormulaVar* > vectFuncPoly;
  vector < RooProduct* > vectFunc;

  RooArgList facList;
  RooArgList funList;
  
  for (int xOrder=0; xOrder<=maxOrder; ++xOrder)
    for (int yOrder=0; yOrder<=maxOrder; ++yOrder)
      for (int zOrder=-1*TMath::Min(xOrder,yOrder); zOrder<=TMath::Min(xOrder,yOrder); ++zOrder) {

	// vector of coefficients for the function basis
	factors.push_back( new RooRealVar( Form("l%i_k%i_m%i",xOrder,yOrder,zOrder),
					   Form("l%i_k%i_m%i",xOrder,yOrder,zOrder), 0 ) );

	RooArgList prodList;

	// phi terms by trigonometric polynomials (degree zOrder)
	if (zOrder>0) {
	  vectFuncPoly.push_back( new RooFormulaVar( Form("funcPoly%i_%i_%i",xOrder,yOrder,zOrder),
						     Form("funcPoly%i_%i_%i",xOrder,yOrder,zOrder),
						     Form("cos(%i*BsphiMC)",zOrder), *BsphiMC ) );
	  prodList.add( *vectFuncPoly.back() );
	}
	if (zOrder<0) {
	  vectFuncPoly.push_back( new RooFormulaVar( Form("funcPoly%i_%i_%i",xOrder,yOrder,zOrder),
						     Form("funcPoly%i_%i_%i",xOrder,yOrder,zOrder),
						     Form("sin(%i*BsphiMC)",-1*zOrder), *BsphiMC ) );
	  prodList.add( *vectFuncPoly.back() );
	}

	// costhetaT terms by associated Legendre polynomials (degree l=xOrder m=zOrder)
	vectFuncLegCosThetaT.push_back( new RooLegendre ( Form("funcLegctK%i_%i_%i",xOrder,yOrder,zOrder),
						     Form("funcLegctK%i_%i_%i",xOrder,yOrder,zOrder),
						     *BscosthetaMC, xOrder, abs(zOrder) ) );
	prodList.add( *vectFuncLegCosThetaT.back() );

	// cospsiT terms by associated Legendre polynomials (degree l=yOrder m=zOrder)
	vectFuncLegCosPsiT.push_back( new RooLegendre ( Form("funcLegctL%i_%i_%i",xOrder,yOrder,zOrder),
						     Form("funcLegctL%i_%i_%i",xOrder,yOrder,zOrder),
						     *BscospsiMC, yOrder, abs(zOrder) ) );
	prodList.add( *vectFuncLegCosPsiT.back() );

	// build member of the basis of 3D functions
	vectFunc.push_back( new RooProduct ( Form("func%i_%i_%i",xOrder,yOrder,zOrder),
					     Form("func%i_%i_%i",xOrder,yOrder,zOrder),
					     prodList ) );

	// coefficients values to be filled later
	proj.push_back(0);
	
	// preparation of RooArgList objects
	funList.add( *vectFunc.back() );
	facList.add( *factors .back() );
//std::cout<<"orders of ploy, legendre"<<xOrder<<yOrder<<zOrder<<"\n";
      } 
  
  cout<<"Number of parameters used: "<<factors.size()<<endl;

  // Sum function
  RooAddition*projectedFunc = new RooAddition("projectedFunc", "projectedFunc", funList, facList );
  
    //Compute and set the coefficients
  double fact;
  int iOrder=-1;

  TStopwatch t;
  t.Start();

  // loop over the coefficients
  for (int xOrder=0; xOrder<=maxOrder; ++xOrder)
    for (int yOrder=0; yOrder<=maxOrder; ++yOrder)
      for (int zOrder=-1*TMath::Min(xOrder,yOrder); zOrder<=TMath::Min(xOrder,yOrder); ++zOrder) {
	
	++iOrder;

	// project the binned efficiency on the [iOrder] function
	for (int xBin=1; xBin<=denHist->GetNbinsX(); ++xBin)
	  for (int yBin=1; yBin<=denHist->GetNbinsY(); ++yBin)
	    for (int zBin=1; zBin<=denHist->GetNbinsZ(); ++zBin) {
	      
	      BscosthetaMC->setVal( denHist->GetXaxis()->GetBinCenter(xBin) );
	      BscospsiMC->setVal( denHist->GetYaxis()->GetBinCenter(yBin) );
	      BsphiMC->setVal( denHist->GetZaxis()->GetBinCenter(zBin) );

	      // contribution of one bin
 	      if ( denHist->GetBinContent(xBin,yBin,zBin)>0 )
		proj[iOrder] += ( numHist->GetBinContent(xBin,yBin,zBin) / denHist->GetBinContent(xBin,yBin,zBin) *
				  denHist->GetXaxis()->GetBinWidth(xBin) *
				  denHist->GetYaxis()->GetBinWidth(yBin) *
				  denHist->GetZaxis()->GetBinWidth(zBin) *
				  vectFunc[iOrder]->getVal( vars ) );

	    }

	// normalization of 0-degree trigonometric polynomial differs by a factor 2
	if (zOrder==0) proj[iOrder] = proj[iOrder]/2.0;
	
	// set coefficient value, normalised
	factors[iOrder]->setVal( proj[iOrder]
				 * (2*xOrder+1)*TMath::Factorial(xOrder-abs(zOrder))/2/TMath::Factorial(xOrder+abs(zOrder)) // associated legendre poly
				 * (2*yOrder+1)*TMath::Factorial(yOrder-abs(zOrder))/2/TMath::Factorial(yOrder+abs(zOrder))
				 / TMath::Pi() // trigonometric polynomial
				 );

	cout<<xOrder<<" "<<yOrder<<" "<<zOrder<<"\t"<<iOrder<<"\t"<<factors[iOrder]->getValV()<<endl;

      }
      

  t.Stop();
  t.Print();    
  
  
  
   // Plot efficiency projections
  gStyle->SetOptStat(0);
  auto h3_xyz = (TH3D*)projectedFunc->createHistogram("BscosthetaMC,BscospsiMC,BsphiMC",50,50,50);
  auto h3_x  = h3_xyz->ProjectionX();
  auto h3_y  = h3_xyz->ProjectionY();
  auto h3_z  = h3_xyz->ProjectionZ();
  auto h3_xy = h3_xyz->Project3D("xy");
  auto h3_xz = h3_xyz->Project3D("xz");
  auto h3_yz = h3_xyz->Project3D("yz");
  
  
    // 2D projections
  auto c3new= new TCanvas("c3new", "c3new",1200,800);
  h3_xy->SetTitle("Efficiency cos(#theta_{T})/cos(#psi_{T}) projected function");
  h3_xz->SetTitle("Efficiency cos(#theta_{T})/#phi projected function");
  h3_yz->SetTitle("Efficiency cos(#psi_{T})/#phi projected function");
  h3_xy->GetXaxis()->SetTitleOffset(1.4);
  h3_xy->GetYaxis()->SetTitleOffset(2);
  h3_xz->GetXaxis()->SetTitleOffset(1.4);
  h3_xz->GetYaxis()->SetTitleOffset(2);
  h3_yz->GetXaxis()->SetTitleOffset(1.4);
  h3_yz->GetYaxis()->SetTitleOffset(2);
  h3_xy->SetMinimum(0.0);
  h3_xz->SetMinimum(0.0);
  h3_yz->SetMinimum(0.0);
  c3new->Divide(3,1);
  c3new->cd(1);
  h3_xy->Draw("SURF3");
  c3new->cd(2);
  h3_xz->Draw("SURF3");
  c3new->cd(3);
  h3_yz->Draw("SURF3");
  c3new->SaveAs(Form("EffProjectionFunction_%i_%i_%i_2DProj_SpH%iOrder.pdf",xbins,ybins,zbins,maxOrder));
  
  
       
    
    
  //=================================================================================Model
       
       
       
 

    RooRealVar *A_0=new RooRealVar ("A_0","|A_0|^2|",0.6,0.30, 0.70);//0.51,0.30,0.7);
    RooRealVar *A_S=new RooRealVar ("A_S","|A_0|^2|",0.05,0.,0.1);
    RooRealVar *A_pe=new RooRealVar ("A_pe","|A_pe|^2",0.16,0.1,0.4);
    RooFormulaVar *A_pa=new RooFormulaVar ("A_pa","|A_pa|*2","1-@0-@1",RooArgList(*A_0,*A_pe));
    //RooRealVar *A_pa=new RooRealVar ("A_pa","|A_pa|^2",0.24,0.2,0.3);
    //RooFormulaVar *A_pe=new RooFormulaVar ("A_pe","|A_pe|*2","1-@0-@1",RooArgList(*A_0,*A_pa));
    RooFormulaVar *A_4=new RooFormulaVar ("A_4","A_4","sqrt(@0)*sqrt(@1)",RooArgList(*A_pa,*A_pe));
    RooFormulaVar *A_5=new RooFormulaVar ("A_5","|A_0||A_pa|","sqrt(@0)*sqrt(@1)",RooArgList(*A_0,*A_pa));
    RooFormulaVar *A_6=new RooFormulaVar ("A_6","A_6","sqrt(@0)*sqrt(@1)",RooArgList(*A_0,*A_pe));
    RooFormulaVar *A_8=new RooFormulaVar ("A_8","A_8","sqrt(@0)*sqrt(@1)",RooArgList(*A_S,*A_pa));
    RooFormulaVar *A_9=new RooFormulaVar ("A_9","|A_9","sqrt(@0)*sqrt(@1)",RooArgList(*A_S,*A_pe));
    RooFormulaVar *A_10=new RooFormulaVar ("A_10","A_10","sqrt(@0)*sqrt(@1)",RooArgList(*A_S,*A_0)); 
    //physical parameters, Phi_s, delta_1, delta_2, dm=Mass_L-Mass_H, DGamma=(Gamma_L-Gamma_H)/2
    
    RooRealVar *phi_s=new RooRealVar ("phi_s","phi_s",-0.04,-1,1);//-0.04,-1.,0.);//-0.04
    RooRealVar *deltaPa=new RooRealVar("deltaPa","#deltaPa",2.5,0.,4.0);//,2.5);   
    RooRealVar *deltaPe=new RooRealVar("deltaPe","#deltaPe",0.17,-0.2,0.2);//,-1,4); // =0.17  
    RooRealVar *deltaSPe=new RooRealVar("deltaSPe","#deltaPe",-0.028,-3.,3.);//,-1,4); // =0.17  
    RooRealVar *dm=new RooRealVar ("dm","dm",589.7,500,700);//700); //590.1// in um*c, PDG value 17.69 in 1/ps
    RooRealVar *dGam=new RooRealVar ("dGam","dGam",0,-3,6);//2.3//original value in MC
    RooRealVar *tau=new RooRealVar ("tau","tau",0.0441,0.03,0.055,"ps");//,0.043//origunal value in MC


   
    RooRealVar *kappa = new RooRealVar("kappa", "kappa", kappa_val);
    RooFormulaVar *kapreso = new RooFormulaVar("kapreso", "kappa*BsCt2DMCErr", RooArgList(*kappa, *BsCt2DMCErr));
    RooRealVar *bias = new RooRealVar("bias","bias",0);
    RooRealVar *sigma = new RooRealVar("sigma","per-event error scale factor",1);
    RooGaussModel gm("gm","gauss model scaled bt per-event error", *BsCt2DMC, *bias, *sigma, *kapreso);
    
    
     RooGenericPdf *cteffFunc = new RooGenericPdf("cteffFunc", "ctp0+ctp1*@7+ctp2*@7*@7+ctp3*@7*@7*@7+ctp4*@7*@7*@7*@7+ctp5*(1/(1+exp(-@7/ctp6)))",RooArgList(*ctp0,*ctp1,*ctp2,*ctp3,*ctp4,*ctp5,*ctp6,*BsCt2DMC));
    

    RooFormulaVar coshGBasis("coshGBasis","exp(-@0/@1)*cosh(@0*@2/2)", RooArgList(*BsCt2DMC ,*tau,*dGam));
    RooFormulaVar sinhGBasis("sinhGBasis","exp(-@0/@1)*sinh(@0*@2/2)", RooArgList(*BsCt2DMC ,*tau,*dGam));
    RooFormulaVar cosGBasis("cosGBasis","exp(-@0/@1)*cos(@0*@2)",RooArgList(*BsCt2DMC ,*tau,*dm));
    RooFormulaVar sinGBasis("sinGBasis","exp(-@0/@1)*sin(@0*@2)",RooArgList(*BsCt2DMC ,*tau,*dm));
    
    
      RooAbsReal* coshGConv0 = gm.convolution(&coshGBasis,BsCt2DMC);
      RooAbsReal* sinhGConv0 = gm.convolution(&sinhGBasis,BsCt2DMC);
      RooAbsReal* cosGConv0 = gm.convolution(&cosGBasis,BsCt2DMC);
      RooAbsReal* sinGConv0 = gm.convolution(&sinGBasis,BsCt2DMC);
      RooProduct *coshGConv = new RooProduct("coshGConv","coshGConv",RooArgSet(*coshGConv0,*cteffFunc));
      RooProduct *sinhGConv = new RooProduct("sinhGConv","sinhGConv",RooArgSet(*sinhGConv0,*cteffFunc));
      RooProduct *cosGConv = new RooProduct("cosGConv","cosGConv",RooArgSet(*cosGConv0,*cteffFunc));//effvarpdf));//cteffFunc));
      RooProduct *sinGConv = new RooProduct("sinGConv","sinGConv",RooArgSet(*sinGConv0,*cteffFunc));

   

    RooRealVar *mean =new RooRealVar("mean","mean",5.36682);
    RooRealVar *sigma1=new RooRealVar("sigma1","sigma1",0.01798,0.0,0.2);
    RooRealVar *sigma2= new RooRealVar("sigma2","sigma2",0.00913,0.01);
    RooRealVar *sigma3=new RooRealVar("sigma3","sigma3",0.022,0.,0.5);

    RooGaussian *gauss1 =new RooGaussian("gauss1","gauss1",*svmass,*mean,*sigma1);
    RooGaussian *gauss2=new RooGaussian("gauss2","gauss2",*svmass,*mean,*sigma2);
    RooGaussian *gauss3=new RooGaussian("gauss3","gauss3",*svmass,*mean,*sigma3);
 
    RooRealVar *frac1=new RooRealVar("frac1","frac1",0.2,0.0,1.0); 
    RooRealVar *frac2=new RooRealVar("frac2","frac2",0.44674);
   
    RooAddPdf *T_gaus = new RooAddPdf("T_gaus","T_gaus",RooArgList(*gauss1,*gauss2,*gauss3),RooArgList(*frac1, *frac2));   
    RooRealVar *nSig=new RooRealVar("nSig", "Number of Signal Events in SIGNAL MC",1400,0.,(int)chain_data->GetEntries());
    RooAddPdf *MTpdf=new RooAddPdf("MTpdf", "Total pdf",RooArgList(*T_gaus), RooArgList(*nSig));
   
              RooRealVar *g1 =new RooRealVar ("g1","g1",8.0440) ;
              RooRealVar *g2 = new RooRealVar("g2","g2",3.798) ;
              RooRealVar *b1 = new RooRealVar("b1","b1",0.00014473) ;
              RooRealVar *b2 = new RooRealVar("b2","b2",0.000541) ;
              RooRealVar *m1 = new RooRealVar("m1","m1", 0.00019) ;
              RooRealVar *m2 = new RooRealVar("m2","m2", 0.00019) ;
              
              RooGamma  *G1 = new RooGamma("G1", "G1", *BsCt2DMCErr, *g1, *b1, *m1);
              RooGamma  *G2 = new RooGamma("G2", "G2", *BsCt2DMCErr, *g2, *b2, *m2);
              
              RooRealVar *errfrac = new RooRealVar("errfrac","errfrac",0.9187) ;              
              RooAddPdf *Di_Gamma = new RooAddPdf("Di_Gamma", "Di_Gamma", RooArgList(*G1,*G2), RooArgList(*errfrac)); 

   


    RooFormulaVar phinew("phinew","BsphiMC/TMath::Pi()",*BsphiMC);
    RooLegendre phizero("phizero","phizero",phinew,0,0);
    RooFormulaVar* fun1T = new RooFormulaVar("fun1T","2*BscospsiMC*BscospsiMC*(1-(1-BscosthetaMC*BscosthetaMC)*cos(BsphiMC)*cos(BsphiMC))",RooArgSet(*BscosthetaMC,*BsphiMC,*BscospsiMC));
    RooFormulaVar* fun2T = new RooFormulaVar("fun2T","(1-BscospsiMC*BscospsiMC)*(1-(1-BscosthetaMC*BscosthetaMC)*sin(BsphiMC)*sin(BsphiMC))",RooArgSet(*BscosthetaMC,*BsphiMC,*BscospsiMC));
    RooFormulaVar* fun3T = new RooFormulaVar("fun3T","(1-BscospsiMC*BscospsiMC)*(1-BscosthetaMC*BscosthetaMC)*phizero",RooArgSet(*BscosthetaMC,*BscospsiMC,phizero));
    RooFormulaVar* fun4T = new RooFormulaVar("fun4T","-(1-BscospsiMC*BscospsiMC)*2*BscosthetaMC*sqrt(1-BscosthetaMC*BscosthetaMC)*sin(BsphiMC)",RooArgSet(*BscosthetaMC,*BsphiMC,*BscospsiMC));
    RooFormulaVar* fun5T = new RooFormulaVar("fun5T","2/sqrt(2.)*BscospsiMC*sqrt(1-BscospsiMC*BscospsiMC)*(1-BscosthetaMC*BscosthetaMC)*sin(2*BsphiMC)",RooArgSet(*BscosthetaMC,*BsphiMC,*BscospsiMC));
    RooFormulaVar* fun6T = new RooFormulaVar("fun6T","2/sqrt(2.)*BscospsiMC*sqrt(1-BscospsiMC*BscospsiMC)*2*BscosthetaMC*sqrt(1-BscosthetaMC*BscosthetaMC)*cos(BsphiMC)",RooArgSet(*BscosthetaMC,*BsphiMC,*BscospsiMC));
    RooFormulaVar* fun7T = new RooFormulaVar("fun7T","2/3*(1-(1-@0*@0)*cos(@1)*cos(@1))",RooArgSet(*BscosthetaMC,*BsphiMC));
    RooFormulaVar* fun8T = new RooFormulaVar("fun8T","sqrt(6.)/3*sqrt(1-BscospsiMC*BscospsiMC)*(1-BscosthetaMC*BscosthetaMC)*sin(2*BsphiMC)",RooArgSet(*BscosthetaMC,*BsphiMC,*BscospsiMC));
    RooFormulaVar* fun9T = new RooFormulaVar("fun9T","sqrt(6.)/3*sqrt(1-BscospsiMC*BscospsiMC)*2*BscosthetaMC*sqrt(1-BscosthetaMC*BscosthetaMC)*cos(BsphiMC)",RooArgSet(*BscosthetaMC,*BsphiMC,*BscospsiMC));
    RooFormulaVar* fun10T = new RooFormulaVar("fun10T","sqrt(3.)*4/3*BscospsiMC*(1-(1-BscosthetaMC*BscosthetaMC)*cos(BsphiMC)*cos(BsphiMC))",RooArgSet(*BscosthetaMC,*BsphiMC,*BscospsiMC));
 
    
    
    
   /*//RooFormulaVar *coef11=new RooFormulaVar("coef11","coef11","@0",RooArgList(*A_0));
    RooFormulaVar *coef12=new RooFormulaVar("coef12","coef12","-cos(@0)",RooArgList(*phi_s));
    RooFormulaVar *coef14=new RooFormulaVar("coef14","coef14","-sin(@0)*@1*(1-2*@2)",RooArgList(*phi_s,*tag,*mistag));
    //RooFormulaVar *coef21=new RooFormulaVar("coef21","coef21","@0",RooArgList(*A_pa));
    RooFormulaVar *coef22=new RooFormulaVar("coef22","coef22","-cos(@0)",RooArgList(*phi_s));
    RooFormulaVar *coef24=new RooFormulaVar("coef24","coef24","-sin(@0)*@1*(1-2*@2)",RooArgList(*phi_s,*tag,*mistag));
    //RooFormulaVar *coef31=new RooFormulaVar("coef31","coef31","@0",RooArgList(*A_pe));
    RooFormulaVar *coef32=new RooFormulaVar("coef32","coef32","cos(@0)",RooArgList(*phi_s));
    RooFormulaVar *coef34=new RooFormulaVar("coef34","coef34","sin(@0)*@1*(1-2*@2)",RooArgList(*phi_s,*tag,*mistag));
    RooFormulaVar *coef42=new RooFormulaVar("coef42","coef42","-sin(@0)*cos(@1-@2)",RooArgList(*phi_s,*deltaPe,*deltaPa));
    RooFormulaVar *coef43=new RooFormulaVar("coef43","coef43","-sin(@0-@1)*@2*(1-2*@3)",RooArgList(*deltaPe,*deltaPa,*tag,*mistag));
    RooFormulaVar *coef44=new RooFormulaVar("coef44","coef44","cos(@0)*cos(@1-@2)*@3*(1-2*@4)",RooArgList(*phi_s,*deltaPe,*deltaPa,*tag,*mistag));
    RooFormulaVar *coef51=new RooFormulaVar("coef51","coef51","cos(@0)",RooArgList(*deltaPa));
    RooFormulaVar *coef52=new RooFormulaVar("coef52","coef52","-cos(@0)*cos(@1)",RooArgList(*deltaPa,*phi_s));
    RooFormulaVar *coef54=new RooFormulaVar("coef54","coef54","-cos(@0)*sin(@1)*@2*(1-2*@3)",RooArgList(*deltaPa,*phi_s,*tag,*mistag));
    RooFormulaVar *coef62=new RooFormulaVar("coef62","coef62","-sin(@0)*cos(@1)",RooArgList(*phi_s,*deltaPe));
    RooFormulaVar *coef63=new RooFormulaVar("coef63","coef63","-sin(@0)*@1*(1-2*@2)",RooArgList(*deltaPe,*tag,*mistag));
    RooFormulaVar *coef64=new RooFormulaVar("coef64","coef64","cos(@0)*cos(@1)*@2*(1-2*@3)",RooArgList(*phi_s,*deltaPe,*tag,*mistag));
    RooFormulaVar *coef72=new RooFormulaVar("coef72","coef72","cos(@0)",RooArgList(*phi_s));
    RooFormulaVar *coef74=new RooFormulaVar("coef74","coef74","sin(@0)*@1*(1-2*@2)",RooArgList(*phi_s,*tag,*mistag));
    RooFormulaVar *coef82=new RooFormulaVar("coef82","coef82","-sin(@0)*sin(@1-@2-@3)",RooArgList(*phi_s,*deltaPa,*deltaSPe,*deltaPe));
    RooFormulaVar *coef83=new RooFormulaVar("coef83","coef83","-cos(@0-@1-@2)*@3*(1-2*@4)",RooArgList(*deltaPa,*deltaSPe,*deltaPe,*tag,*mistag));
    RooFormulaVar *coef84=new RooFormulaVar("coef84","coef84","sin(@0-@1-@2)*@3*(1-2*@4)*cos(@5)",RooArgList(*deltaPa,*deltaSPe,*deltaPe,*tag,*mistag,*phi_s));
    RooFormulaVar *coef91=new RooFormulaVar("coef91","coef91","sin(-@0)",RooArgList(*deltaSPe));
    RooFormulaVar *coef92=new RooFormulaVar("coef92","coef92","sin(-@0)*cos(@1)",RooArgList(*deltaSPe,*phi_s));
    RooFormulaVar *coef94=new RooFormulaVar("coef94","coef94","sin(-@0)*sin(@1)*@2*(1-2*@3)",RooArgList(*deltaSPe,*phi_s,*tag,*mistag));
    RooFormulaVar *coef102=new RooFormulaVar("coef102","coef102","-sin(@0)*sin(-@1-@2)",RooArgList(*phi_s,*deltaSPe,*deltaPe));
    RooFormulaVar *coef103=new RooFormulaVar("coef103","coef103","-cos(-@0-@1)*@2*(1-2*@3)",RooArgList(*deltaSPe,*deltaPe,*tag,*mistag));
    RooFormulaVar *coef104=new RooFormulaVar("coef104","coef104","sin(-@0-@1)*@2*(1-2*@3)*cos(@4)",RooArgList(*deltaSPe,*deltaPe,*tag,*mistag,*phi_s));
*/

 
    RooFormulaVar *coef12=new RooFormulaVar("coef12","coef12","-cos(@0)",RooArgList(*phi_s));
    RooFormulaVar *coef14=new RooFormulaVar("coef14","coef14","sin(@0)*@1*(1-2*@2)",RooArgList(*phi_s,*tag,*mistag));  
    RooFormulaVar *coef22=new RooFormulaVar("coef22","coef22","-cos(@0)",RooArgList(*phi_s));
    RooFormulaVar *coef24=new RooFormulaVar("coef24","coef24","sin(@0)*@1*(1-2*@2)",RooArgList(*phi_s,*tag,*mistag));  
    RooFormulaVar *coef32=new RooFormulaVar("coef32","coef32","cos(@0)",RooArgList(*phi_s));
    RooFormulaVar *coef34=new RooFormulaVar("coef34","coef34","-sin(@0)*@1*(1-2*@2)",RooArgList(*phi_s,*tag,*mistag));
    RooFormulaVar *coef42=new RooFormulaVar("coef42","coef42","-sin(@0)*cos(@1-@2)",RooArgList(*phi_s,*deltaPe,*deltaPa));
    RooFormulaVar *coef43=new RooFormulaVar("coef43","coef43","sin(@0-@1)*@2*(1-2*@3)",RooArgList(*deltaPe,*deltaPa,*tag,*mistag));
    RooFormulaVar *coef44=new RooFormulaVar("coef44","coef44","-cos(@0)*cos(@1-@2)*@3*(1-2*@4)",RooArgList(*phi_s,*deltaPe,*deltaPa,*tag,*mistag));
    RooFormulaVar *coef51=new RooFormulaVar("coef51","coef51","cos(@0)",RooArgList(*deltaPa));
    RooFormulaVar *coef52=new RooFormulaVar("coef52","coef52","-cos(@0)*cos(@1)",RooArgList(*deltaPa,*phi_s));
    RooFormulaVar *coef54=new RooFormulaVar("coef54","coef54","cos(@0)*sin(@1)*@2*(1-2*@3)",RooArgList(*deltaPa,*phi_s,*tag,*mistag));
    RooFormulaVar *coef62=new RooFormulaVar("coef62","coef62","-sin(@0)*cos(@1)",RooArgList(*phi_s,*deltaPe));
    RooFormulaVar *coef63=new RooFormulaVar("coef63","coef63","sin(@0)*@1*(1-2*@2)",RooArgList(*deltaPe,*tag,*mistag));
    RooFormulaVar *coef64=new RooFormulaVar("coef64","coef64","-cos(@0)*cos(@1)*@2*(1-2*@3)",RooArgList(*phi_s,*deltaPe,*tag,*mistag));
    RooFormulaVar *coef72=new RooFormulaVar("coef72","coef72","cos(@0)",RooArgList(*phi_s));
    RooFormulaVar *coef74=new RooFormulaVar("coef74","coef74","-sin(@0)*@1*(1-2*@2)",RooArgList(*phi_s,*tag,*mistag));
    RooFormulaVar *coef82=new RooFormulaVar("coef82","coef82","-sin(@0)*sin(@1-@2-@3)",RooArgList(*phi_s,*deltaPa,*deltaSPe,*deltaPe));
    RooFormulaVar *coef83=new RooFormulaVar("coef83","coef83","cos(@0-@1-@2)*@3*(1-2*@4)",RooArgList(*deltaPa,*deltaSPe,*deltaPe,*tag,*mistag));
    RooFormulaVar *coef84=new RooFormulaVar("coef84","coef84","-sin(@0-@1-@2)*@3*(1-2*@4)*cos(@5)",RooArgList(*deltaPa,*deltaSPe,*deltaPe,*tag,*mistag,*phi_s));
    RooFormulaVar *coef91=new RooFormulaVar("coef91","coef91","sin(-@0)",RooArgList(*deltaSPe));
    RooFormulaVar *coef92=new RooFormulaVar("coef92","coef92","sin(-@0)*cos(@1)",RooArgList(*deltaSPe,*phi_s));
    RooFormulaVar *coef94=new RooFormulaVar("coef94","coef94","sin(-@0)*sin(@1)*@2*(1-2*@3)",RooArgList(*deltaSPe,*phi_s,*tag,*mistag));
    RooFormulaVar *coef102=new RooFormulaVar("coef102","coef102","-sin(@0)*sin(-@1-@2)",RooArgList(*phi_s,*deltaSPe,*deltaPe));
    RooFormulaVar *coef103=new RooFormulaVar("coef103","coef103","cos(-@0-@1)*@2*(1-2*@3)",RooArgList(*deltaSPe,*deltaPe,*tag,*mistag));
    RooFormulaVar *coef104=new RooFormulaVar("coef104","coef104","-sin(-@0-@1)*@2*(1-2*@3)*cos(@4)",RooArgList(*deltaSPe,*deltaPe,*tag,*mistag,*phi_s));
    
    



    
    
    
    RooProduct *eq11=new RooProduct("eq11","amp0*g1",RooArgSet(*coshGConv,*fun1T,*projectedFunc));
    RooProduct *eq12=new RooProduct("eq12","amp0*g1",RooArgSet(*coef12,*sinhGConv,*fun1T,*projectedFunc));
    RooProduct *eq14=new RooProduct("eq14","amp0*g1",RooArgSet(*coef14,*sinGConv,*fun1T,*projectedFunc));
    RooProduct *eq21=new RooProduct("eq21","amp0*g1",RooArgSet(*coshGConv,*fun2T,*projectedFunc));
    RooProduct *eq22=new RooProduct("eq22","amp0*g1",RooArgSet(*coef22,*sinhGConv,*fun2T,*projectedFunc));
    RooProduct *eq24=new RooProduct("eq24","amp0*g1",RooArgSet(*coef24,*sinGConv,*fun2T,*projectedFunc));
    RooProduct *eq31=new RooProduct("eq31","amp0*g1",RooArgSet(*coshGConv,*fun3T,*projectedFunc));
    RooProduct *eq32=new RooProduct("eq32","amp0*g1",RooArgSet(*coef32,*sinhGConv,*fun3T,*projectedFunc));
    RooProduct *eq34=new RooProduct("eq34","amp0*g1",RooArgSet(*coef34,*sinGConv,*fun3T,*projectedFunc));
    RooProduct *eq42=new RooProduct("eq42","amp0*g1",RooArgSet(*coef42,*sinhGConv,*fun4T,*projectedFunc));
    RooProduct *eq43=new RooProduct("eq43","amp0*g1",RooArgSet(*coef43,*cosGConv,*fun4T,*projectedFunc));
    RooProduct *eq44=new RooProduct("eq44","amp0*g1",RooArgSet(*coef44,*sinGConv,*fun4T,*projectedFunc));
    RooProduct *eq51=new RooProduct("eq51","amp0*g1",RooArgSet(*coef51,*coshGConv,*fun5T,*projectedFunc));
    RooProduct *eq52=new RooProduct("eq52","amp0*g1",RooArgSet(*coef52,*sinhGConv,*fun5T,*projectedFunc));
    RooProduct *eq54=new RooProduct("eq54","amp0*g1",RooArgSet(*coef54,*sinGConv,*fun5T,*projectedFunc));
    RooProduct *eq62=new RooProduct("eq62","amp0*g1",RooArgSet(*coef62,*sinhGConv,*fun6T,*projectedFunc));
    RooProduct *eq63=new RooProduct("eq63","amp0*g1",RooArgSet(*coef63,*cosGConv,*fun6T,*projectedFunc));
    RooProduct *eq64=new RooProduct("eq64","amp0*g1",RooArgSet(*coef64,*sinGConv,*fun6T,*projectedFunc));
    RooProduct *eq71=new RooProduct("eq71","amp0*g1",RooArgSet(*coshGConv,*fun7T,*projectedFunc));
    RooProduct *eq72=new RooProduct("eq72","amp0*g1",RooArgSet(*coef72,*sinhGConv,*fun7T,*projectedFunc));
    RooProduct *eq74=new RooProduct("eq74","amp0*g1",RooArgSet(*coef74,*sinGConv,*fun7T,*projectedFunc));
    RooProduct *eq82=new RooProduct("eq82","amp0*g1",RooArgSet(*coef82,*sinhGConv,*fun8T,*projectedFunc));
    RooProduct *eq83=new RooProduct("eq83","amp0*g1",RooArgSet(*coef83,*cosGConv,*fun8T,*projectedFunc));
    RooProduct *eq84=new RooProduct("eq84","amp0*g1",RooArgSet(*coef84,*sinGConv,*fun8T,*projectedFunc));
    RooProduct *eq91=new RooProduct("eq91","amp0*g1",RooArgSet(*coef91,*coshGConv,*fun9T,*projectedFunc));
    RooProduct *eq92=new RooProduct("eq92","amp0*g1",RooArgSet(*coef92,*sinhGConv,*fun9T,*projectedFunc));
    RooProduct *eq94=new RooProduct("eq94","amp0*g1",RooArgSet(*coef94,*sinGConv,*fun9T,*projectedFunc));
    RooProduct *eq102=new RooProduct("eq102","amp0*g1",RooArgSet(*coef102,*sinhGConv,*fun10T,*projectedFunc));
    RooProduct *eq103=new RooProduct("eq103","amp0*g1",RooArgSet(*coef103,*cosGConv,*fun10T,*projectedFunc));
    RooProduct *eq104=new RooProduct("eq104","amp0*g1",RooArgSet(*coef104,*sinGConv,*fun10T,*projectedFunc));
    
    RooArgList funtot(*eq11,*eq12,*eq14,*eq21,*eq22,*eq24,*eq31,*eq32,*eq34);
    funtot.add(*eq42);
    funtot.add(*eq43);
    funtot.add(*eq44);
    funtot.add(*eq51);
    funtot.add(*eq52);
    funtot.add(*eq54);
    funtot.add(*eq62);
    funtot.add(*eq63);
    funtot.add(*eq64);
    RooArgList coeftot(*A_0,*A_0,*A_0,*A_pa,*A_pa,*A_pa,*A_pe,*A_pe,*A_pe);
    coeftot.add(*A_4);
    coeftot.add(*A_4);
    coeftot.add(*A_4);
    coeftot.add(*A_5);
    coeftot.add(*A_5);
    coeftot.add(*A_5);
    coeftot.add(*A_6);
    coeftot.add(*A_6);
    coeftot.add(*A_6);
    
    RooArgList funtotP(*eq11,*eq12,*eq14,*eq21,*eq22,*eq24,*eq31,*eq32,*eq34);
    funtotP.add(*eq42);
    funtotP.add(*eq43);
    funtotP.add(*eq44);
    funtotP.add(*eq51);
    funtotP.add(*eq52);
    funtotP.add(*eq54);
    funtotP.add(*eq62);
    funtotP.add(*eq63);
    funtotP.add(*eq64);
    funtotP.add(*eq71);
    funtotP.add(*eq72);
    funtotP.add(*eq74);
    funtotP.add(*eq82);
    funtotP.add(*eq83);
    funtotP.add(*eq84);
    funtotP.add(*eq91);
    funtotP.add(*eq92);
    funtotP.add(*eq94);
    funtotP.add(*eq102);
    funtotP.add(*eq103);
    funtotP.add(*eq104);
    RooArgList coeftotP(*A_0,*A_0,*A_0,*A_pa,*A_pa,*A_pa,*A_pe,*A_pe,*A_pe);
    coeftotP.add(*A_4);
    coeftotP.add(*A_4);
    coeftotP.add(*A_4);
    coeftotP.add(*A_5);
    coeftotP.add(*A_5);
    coeftotP.add(*A_5);
    coeftotP.add(*A_6);
    coeftotP.add(*A_6);
    coeftotP.add(*A_6);
    coeftotP.add(*A_S);
    coeftotP.add(*A_S);
    coeftotP.add(*A_S);
    coeftotP.add(*A_8);
    coeftotP.add(*A_8);
    coeftotP.add(*A_8);
    coeftotP.add(*A_9);
    coeftotP.add(*A_9);
    coeftotP.add(*A_9);
    coeftotP.add(*A_10);
    coeftotP.add(*A_10);
    coeftotP.add(*A_10);
    
      RooRealSumPdf *PDFdefP=new RooRealSumPdf("PDFdefP","Signal PDF Amp_i*funiT",funtotP,coeftotP);
    RooRealSumPdf *PDFdef=new RooRealSumPdf("PDFdef","Signal PDF Amp_i*funiT",funtot,coeftot);
  //RooRealSumPdf *Angular_Model=new RooRealSumPdf("Angular_Model","Signal PDF Amp_i*funiT",funtot,coeftot);
    RooProdPdf *Angular_Model=new RooProdPdf("Angular_Model","PDF_Signaltag",RooArgList(*MTpdf,*PDFdef, *Di_Gamma),Conditional(*PDFdef,RooArgList(*BsCt2DMC,*tag,*BscosthetaMC,*BscospsiMC,*BsphiMC,*mistag)));
//RooProdPdf *Angular_Model=new RooProdPdf("Angular_Model","PDF_Signaltag",RooArgList(*MTpdf,*PDFdefP),Conditional(*PDFdefP,RooArgList(*BsCt2DMC,*BscosthetaMC,*BscospsiMC,*BsphiMC,*mistag)));

    //RooRealSumPdf *PDFtag=new RooRealSumPdf("PDFtag","Signal PDF Amp_i*funiT",RooArgList(*pdf1tag,*pdf2tag,*pdf3tag,*pdf4tag,*pdf5tag,*pdf6tag),RooArgList(*A_0,*A_pa,*A_pe,*A_4,*A_5,*A_6));//,Conditional(gm,*BdCt2DBS));
   // RooRealSumPdf *Angular_Model=new RooRealSumPdf("Angular_Model","Signal PDF Amp_i*funiT",RooArgList(*pdf1tag,*pdf2tag,*pdf3tag,*pdf4tag,*pdf5tag,*pdf6tag),RooArgList(*A_0,*A_pa,*A_pe,*A_4,*A_5,*A_6));//,Conditional(gm,*BdCt2DBS));
               //RooProdPdf *Angular_Model = new RooProdPdf("Angular_Model", "Angular_Model", *PDFtag, *T_gaus);
               RooFitResult* PwaveMCResults = Angular_Model->fitTo(*data, NumCPU(8), Save(kTRUE));//, ConditionalObservables(*BsCt2DMCErr));
               PwaveMCResults->Print("v");
    
    //Construct 2D color plot of correlation matrix
    gStyle->SetOptStat(0) ;
    gStyle->SetPalette(1) ;
    TH2* hcorr = PwaveMCResults->correlationHist() ;
    TCanvas* c = new TCanvas("Correlation Matrix","Correlation Matrix",1000,600) ;
    hcorr->GetYaxis()->SetTitleOffset(1.4) ; hcorr->Draw("colz"); c->Print("corrMatrix.png","png");
              TString title1 = "cos#theta_{T}";
               TString title2 = "cos#psi_{T}";
               TString title3 = "#phi_{T}";
               TString title4 = "B_{s} ct";               
               RooPlot* CosTheta = BscosthetaMC->frame(Title(title1),Bins(nbins));
               data->plotOn(CosTheta,DataError(RooAbsData::SumW2));
               Angular_Model->plotOn(CosTheta) ;
               Angular_Model->paramOn(CosTheta);
               Angular_Model->plotOn(CosTheta, LineColor(kBlack), LineWidth(1));
               Double_t chisquare_costheta = CosTheta->chiSquare();
    
               
  /*             RooPlot* CosPsi = BscospsiMC->frame(Title(title2),Bins(nbins));
               data->plotOn(CosPsi,DataError(RooAbsData::SumW2));
               Angular_Model->plotOn(CosPsi);
               Angular_Model->paramOn(CosPsi);
               Angular_Model->plotOn(CosPsi, LineColor(kBlack), LineWidth(1));
               Double_t chisquare_cospsi = CosPsi->chiSquare();
    
               
               RooPlot* PHI = BsphiMC->frame(Title(title3),Bins(nbins));
               data->plotOn(PHI, DataError(RooAbsData::SumW2));
               Angular_Model->plotOn(PHI);
               Angular_Model->paramOn(PHI);
               Angular_Model->plotOn(PHI , LineColor(kBlack), LineWidth(1));
               Double_t chisquare_phi = PHI->chiSquare();
             
 
                RooPlot* ct = BsCt2DMC->frame(Title(title4),Bins(nbins));
                data->plotOn(ct, DataError(RooAbsData::SumW2));
                Angular_Model->plotOn(ct);
                Angular_Model->paramOn(ct);
                Angular_Model->plotOn(ct , LineColor(kBlack), LineWidth(1));
                Double_t chisquare_ct = ct->chiSquare();

                 RooPlot* bsmass = svmass->frame(Title("M_{B_{s}} (GeV/c^{2})"),Bins(nbins));     
                 data->plotOn(bsmass,DataError(RooAbsData::SumW2));
                 Angular_Model->plotOn(bsmass) ;
                 Angular_Model->paramOn(bsmass);
                 Angular_Model->plotOn(bsmass, LineColor(kBlue), LineWidth(1));
                 Double_t chisquare_mass = bsmass->chiSquare();
   
    

    Angular_Model->plotOn(bsmass, Components(*gauss1), LineColor(3), LineWidth(2), LineStyle(4));
    Angular_Model->plotOn(bsmass, Components(*gauss2), LineColor(2), LineWidth(2), LineStyle(2));
    Angular_Model->plotOn(bsmass, Components(*gauss3), LineColor(6), LineWidth(2), LineStyle(2));            
*/

               /*BsCt2DMC->setRange("Asymme",0.02,0.08) ;
               RooPlot* asymplot=BsCt2DMC->frame(Title("BsCt Asymetry"),Range("Asymme"),Binning(300));
               data->plotOn(asymplot,LineStyle(kDashed),Binning(300),RooFit::Name("Bs ct Asym"),Asymmetry(*tag),Range("Asymme"));
               Angular_Model->plotOn(asymplot,RooFit::Name("Fitted ct Asym"),Asymmetry(*tag),Range("Asymme"));
 
               RooPlot* framepull = BsCt2DMC->frame(RooFit::Title("Ct pull"));
               RooHist* hpull = ct->pullHist() ;
               framepull->addPlotable(hpull,"P0") ;
               framepull->SetMinimum(-5) ;
               framepull->SetMaximum(+5) ;
               framepull->SetYTitle("pull");
               framepull->SetMarkerStyle(20);               
               framepull->SetNdivisions(10);

               std::cout<<"Num of fit param P wave reco results:"<<PwaveMCResults->floatParsFinal().getSize()<<endl;
               std::cout<<"Chi2 ct :"<< chisquare_ct<<"\n";
               std::cout<<"Chi2  costheta :"<< chisquare_costheta<<"\n";
               std::cout<<"Chi2  cospsi :"<< chisquare_cospsi<<"\n";
               std::cout<<"Chi2  phi :"<< chisquare_phi<<"\n";
               std::cout<<"Chi2 Mass :"<<chisquare_mass<<"\n";
*/
               TCanvas *cc = new TCanvas("cc", "cc",0,0,800,600);
               CosTheta->Draw();
               cc->Print("ctheta.root", "root");
               /*TCanvas *cc11 = new TCanvas("cc11", "cc11",0,0,800,600);
               CosPsi->Draw();
               //cc11->Print("cpsi.root", "root");
               TCanvas *cc22 = new TCanvas("cc22", "cc22",0,0,800,600);
               PHI->Draw();
               //cc22->Print("phit.root", "root");
               TCanvas *cc33 = new TCanvas("cc33", "cc33",0,0,800,600);
               ct->Draw();
               //cc33->Print("ct.root", "root");
               TCanvas *cc44 = new TCanvas("cc44", "cc44",0,0,800,600);
               asymplot->Draw();
               //cc44->Print("Asymmetry.root", "root");
               TCanvas *cc55 = new TCanvas("cc55", "cc55",0,0,800,600);
               framepull->Draw();
               //cc55->Print("framepull.root", "root");
               TCanvas *cc66 = new TCanvas("cc66", "cc66",0,0,800,600);
               bsmass->Draw();
               TFile *output = new TFile("FitResults_Reco_GenTag.root","recreate");
               PwaveMCResults->Write();
               asymplot->Write("AsymmetryPlot");
               framepull->Write("LifetimePull") ;
               CosTheta->Write("Costhetafit");
               CosPsi->Write("Cospsifit");
               PHI->Write("Phifit");
               ct->Write("Lifetimefit");
               bsmass->Write("B_{s} mass"); 
               hcorr->Write("CorrelationMatrix");
               output->Close();

               Angular_Model->graphVizTree ("FitResults_Reco_GenTag.dot") ;



*/
  /*TCanvas c1 ;
  RooPlot * frame = BscosthetaMC->frame () ; 
 Angular_Model->plotOn (frame, LineColor (kRed + 2)) ;  
 Angular_Model->plotOn (frame, LineColor (kBlue + 2)) ;
  frame->Draw () ;	
  c1.Print ("modelPlot.gif","gif") ;
    */
    }
