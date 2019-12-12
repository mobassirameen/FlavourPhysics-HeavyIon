#include "RooRealVar.h"
#include "RooDataSet.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TMath.h"
#include "RooPlot.h"

using namespace RooFit;
using namespace std;

static const int nBins = 9;
float binBorders [nBins+1] = { 1, 2, 4.3, 6, 8.68, 10.09, 12.86, 14.18, 16, 19};

void plotEffBin(bool tagFlag, int maxOrder, int xbins, int ybins, int zbins);

TCanvas* c  [2*nBins];
TCanvas* c1 [2*nBins];
TCanvas* c3 [2*nBins];

void plotEff(int tagFlag=1, int maxOrder=5, int xbins=25, int ybins = 0, int zbins = 0)
{
  

  if ( maxOrder < 0 ) return;

  if ( ybins<1 ) ybins = xbins;
  if ( zbins<1 ) zbins = xbins;
  if ( xbins<1 ) return;

 
    if (tagFlag < 2 && tagFlag > -1) {
      
      plotEffBin( (tagFlag==1), maxOrder, xbins, ybins, zbins);
    }
    if (tagFlag == 2) {
     
      plotEffBin(true,  maxOrder, xbins, ybins, zbins);
      plotEffBin(false, maxOrder, xbins, ybins, zbins);
    }
  

  
      if (tagFlag < 2 && tagFlag > -1) {
	
	plotEffBin((tagFlag==1), maxOrder, xbins, ybins, zbins);
      }
      if (tagFlag == 2) {
	
	plotEffBin(true,  maxOrder, xbins, ybins, zbins);
	plotEffBin(false, maxOrder, xbins, ybins, zbins);
      }
    
  
  
}

void plotEffBin(bool tagFlag, int maxOrder, int xbins, int ybins, int zbins)
{

  //string shortString = Form(tagFlag?"b%ict":"b%iwt",q2Bin);
  //string longString  = Form(tagFlag?"q2 bin %i correct-tag":"q2 bin %i wrong-tag",q2Bin);
  //int confIndex = (tagFlag?q2Bin:q2Bin+nBins);

  // open file with efficiency and import efficiency function and variables
  TFile* fin = new TFile("effProjections_25_25_25.root", "READ" ); 
  RooWorkspace* ws = (RooWorkspace*)fin->Get("ws"); 
  RooAbsReal* eff = ws->function("projectedFunc");
  
  RooRealVar* costhetaMC = ws->var("costhetaMC");
  RooRealVar* cospsiMC = ws->var("cospsiMC");
  RooRealVar* phiMC= ws->var("phiMC");
  RooArgSet vars (*costhetaMC,*cospsiMC,*phiMC);
  // set a variable for weights, to correct GEN distributions in closure test
  RooRealVar wei ("wei","weight",-1,1);

  // Plot efficiency projections
  gStyle->SetOptStat(0);
  auto h3_xyz = (TH3D*)eff->createHistogram("costhetaMC,cospsiMC,phiMC",50,50,50);
  auto h3_x  = h3_xyz->ProjectionX();
  auto h3_y  = h3_xyz->ProjectionY();
  auto h3_z  = h3_xyz->ProjectionZ();
  auto h3_xy = h3_xyz->Project3D("xy");
  auto h3_xz = h3_xyz->Project3D("xz");
  auto h3_yz = h3_xyz->Project3D("yz");

  // 1D projections
 auto c1new = new TCanvas("c1new", "c1new",1200,800);
  h3_x->SetTitle("Efficiency cos(#theta_{T}) projection");
  h3_y->SetTitle("Efficiency cos(#psi_{T}) projection");
  h3_z->SetTitle("Efficiency #phi projection");
  h3_x->SetLineColor(2);
  h3_y->SetLineColor(2);
  h3_z->SetLineColor(2);
  h3_x->SetMinimum(0.0);
  h3_y->SetMinimum(0.0);
  h3_z->SetMinimum(0.0);
  c1new->Divide(3,1);
  c1new->cd(1);
  h3_x->Draw();
  c1new->cd(2);
  h3_y->Draw();
  c1new->cd(3);
  h3_z->Draw();
  c1new->SaveAs(Form("EffProjections_%i_%i_%i_1DProj_sh%io.pdf",xbins,ybins,zbins,maxOrder) );

  // 2D projections
  auto c3new= new TCanvas("c3new", "c3new",1200,800);
  h3_xy->SetTitle("Efficiency cos(#theta_{T})/cos(#psi_{T}) projection ");
  h3_xz->SetTitle("Efficiency cos(#theta_{T})/#phi projection");
  h3_yz->SetTitle("Efficiency cos(#psi_{T})/#phi projection");
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
  c3new->SaveAs(Form("EffProj_%i_%i_%i_2DProj_sh%io.pdf",xbins,ybins,zbins,maxOrder));

  // Load ntuples
  TChain* t_den = new TChain();
  TChain* t_num = new TChain();
  t_den->Add("/home/cms/BMeson/ffrepo/Gen_Study/BsGenV13.root/OutTreeGEN");
  t_num->Add("/home/cms/BMeson/ffrepo/ntuple/BsMC17_JpsiMu.root/fit");
  int denEntries = t_den->GetEntries();
  int numEntries = t_num->GetEntries();
  int counter;

 Float_t costhetagen, cospsigen, phigen, gentag;
  Double_t recoCosTheta, recoCosPsi, recoPhi;
  Int_t recotag;
  t_den->SetBranchAddress( "angle_costheta_GEN" , &costhetagen );
  t_den->SetBranchAddress( "angle_cospsi_GEN" , &cospsigen);
  t_den->SetBranchAddress( "angle_phi_GEN", &phigen);
  t_den->SetBranchAddress( "MC_Flavour", &gentag);
  t_num->SetBranchAddress( "BscosthetaMC" , &recoCosTheta );
  t_num->SetBranchAddress( "BscospsiMC" , &recoCosPsi );
  t_num->SetBranchAddress( "BsphiMC", &recoPhi );
  t_num->SetBranchAddress( "tag", &recotag);

  RooDataSet* recodata    = new RooDataSet( "recodata"   , "GEN distribution",RooArgSet(vars,wei) );
  RooDataSet* numData = new RooDataSet( "numData", "RECO distribution after selections", vars ); 

  
  // Set counter to evaluate size of negatine-efficiency regions
  int badCounter = 0;
  int totalCounter = 0;

  // Prepare denominator datasets
  cout<<"Starting denominator dataset filling..."<<endl;
  counter=0;
  for (int iCand=0; iCand<denEntries; ++iCand) {
    t_den->GetEntry(iCand);
   
    // status display
    if ( iCand > 1.0*counter*denEntries/100 ) {
      cout<<counter<<"%"<<endl;
      counter += 10;
    }
    // fill
   costhetaMC->setVal(costhetagen);
    cospsiMC->setVal(cospsigen);
    phiMC->setVal(phigen);
    // set event weight based on local efficiency value
    wei.setVal( eff->getVal( vars ) );
    recodata->add( RooArgSet(vars,wei) );
    // count negative efficiency values
    ++totalCounter;
    if ( wei.getValV() < 0 ) ++badCounter;
  }
  // create the weighted dataset for GEN events
  RooDataSet* wdata = new RooDataSet(recodata->GetName(),recodata->GetTitle(),recodata,*recodata->get(),0,wei.GetName());

  // Prepare numerator datasets
  cout<<"Starting numerator dataset filling..."<<endl;
  counter=0;
  for (int iCand=0; iCand<numEntries; ++iCand) {
    t_num->GetEntry(iCand);
       // status display
    if ( iCand > 1.0*counter*numEntries/100 ) {
      cout<<counter<<"%"<<endl;
      counter += 10;
    }
    // fill
    costhetaMC->setVal(recoCosTheta);
    cospsiMC->setVal(recoCosPsi);
    phiMC->setVal(recoPhi);
    numData->add(vars);
  }
  cout<<"Dataset prepared"<<endl;

  // Print size of negative-efficiency regions
  cout<<"Negative efficiency phase-space fraction: "<<badCounter<<"/"<<totalCounter<<" -> "<<1.0*badCounter/totalCounter<<endl;

  // Plot projections for closure test (RECO vs. eff*GEN)
  auto cpnew = new TCanvas ("cnew","Closure test ",1200,800);
  TLegend* leg = new TLegend (0.25,0.8,0.9,0.9);
  RooPlot* xframe = costhetaMC->frame(Title("cos(#theta_{T}) distributions"));
  RooPlot* yframe = cospsiMC->frame(Title("cos(#psi_{T}) distributions"));
  RooPlot* zframe = phiMC->frame(Title("#phi distributions"));
  wdata->plotOn(xframe,MarkerColor(kRed+1),LineColor(kRed+1),Binning(30),DataError(RooAbsData::SumW2),Name("plDenDist"));
  wdata->plotOn(yframe,MarkerColor(kRed+1),LineColor(kRed+1),Binning(30),DataError(RooAbsData::SumW2));
  wdata->plotOn(zframe,MarkerColor(kRed+1),LineColor(kRed+1),Binning(30),DataError(RooAbsData::SumW2));
  numData->plotOn(xframe,Binning(30),Name("plNumDist"));
  numData->plotOn(yframe,Binning(30));
  numData->plotOn(zframe,Binning(30));
  xframe->GetYaxis()->SetTitleOffset(1.6);
  yframe->GetYaxis()->SetTitleOffset(1.6);
  zframe->GetYaxis()->SetTitleOffset(1.6);
  xframe->SetMaximum(xframe->GetMaximum()*1.15);
  yframe->SetMaximum(yframe->GetMaximum()*1.15);
  zframe->SetMaximum(zframe->GetMaximum()*1.15);
  leg->SetTextSize(0.03);
  leg->AddEntry(xframe->findObject("plNumDist"),"Post-selection RECO distribution" ,"lep");
  leg->AddEntry(xframe->findObject("plDenDist"),"Efficiency-corrected GEN distribution" ,"lep");

  cpnew->Divide(3,1);
  cpnew->cd(1);
  gPad->SetLeftMargin(0.17); 
  xframe->Draw();
  leg->Draw("same");
  cpnew->cd(2);
  gPad->SetLeftMargin(0.17); 
  yframe->Draw();
  leg->Draw("same");
  cpnew->cd(3);
  gPad->SetLeftMargin(0.17); 
  zframe->Draw();
  leg->Draw("same");

  cpnew->SaveAs(Form("Closure_%i_%i_%i_sh%io.png",xbins,ybins,zbins,maxOrder));

}
