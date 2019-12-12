#include "RooRealVar.h"
#include "RooDataSet.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TMath.h"
#include "RooPlot.h"
using namespace RooFit;
using namespace std;


void angular_eff_3dim(string filename, int maxOrder=7, int xbins=70, int ybins = 70, int zbins = 30)
{
  if ( maxOrder < 0 ) return;
  if ( ybins<1 ) ybins = xbins;
  if ( zbins<1 ) zbins = xbins;
  if ( xbins<1 ) return;

RooRealVar *BscosthetaMC=new RooRealVar("BscosthetaMC","cos(#theta_{T})",-1,1);
RooRealVar *BscospsiMC=new RooRealVar("BscospsiMC","cos(#psi_{T})",-1,1);
RooRealVar *BsphiMC = new RooRealVar("BsphiMC","#phi",-TMath::Pi(),TMath::Pi());
RooArgSet vars(*BscosthetaMC, *BscospsiMC, *BsphiMC);

    TChain* chain_data = new TChain("treeFit");
    chain_data->Add(filename.c_str());
    Int_t nevt = (int)chain_data->GetEntries();
    
    RooDataSet* data = new RooDataSet("data", "raw data1", RooArgSet(*BscosthetaMC, *BscospsiMC, *BsphiMC), Import(*chain_data));
    
    
  TChain* t_den = new TChain();
  TChain* t_num = new TChain();
    t_den->Add("/Users/ab/Documents/Data_MC_Sample/ntuBsGEN.root/OutTree");
    t_num->Add("/Users/ab/Documents/Data_MC_Sample/fittree_ntuBsMC2017.root/treeFit");
  Int_t denEntries = t_den->GetEntries();
  Int_t numEntries = t_num->GetEntries();
  Int_t counter;

  Float_t costhetagen, cospsigen, phigen;
  Float_t recoCosTheta, recoCosPsi, recoPhi;
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
                     c1->SaveAs(Form("EffProjectionHist_%i_%i_%i_2DProj_SpH%iOrder.pdf",xbins,ybins,zbins,maxOrder));
                     
  vector<RooRealVar*> factors;
  vector <double> proj;
  vector <RooLegendre*> vectFuncLegCosThetaT;
  vector <RooLegendre*> vectFuncLegCosPsiT;
  vector <RooFormulaVar*> vectFuncPoly;
  vector <RooProduct*> vectFunc;

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
  RooAddition *projectedFunc = new RooAddition("projectedFunc", "projectedFunc", funList, facList );
//Compute and set the coefficients=================
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
      
    RooWorkspace *ws = new RooWorkspace("ws","Workspace with efficiency parameterisation");
    ws->import( *projectedFunc, Silence() );
    double costheta, cospsi, phi;
    RooArgList* EffCoeff = new RooArgList("EffCoeff");
    int k_ord, l_ord, m_ord;
    for (k_ord=0; k_ord<=maxOrder; ++k_ord)
        for (l_ord=0; l_ord<=maxOrder; ++l_ord)
            for (m_ord=-1*TMath::Min(k_ord,l_ord); m_ord<=TMath::Min(k_ord,l_ord); ++m_ord)
                EffCoeff->add(*ws->var(Form("l%i_k%i_m%i",k_ord,l_ord,m_ord)));
    
    int degree = 7;
    std::vector< double > LegCTheta;
    std::vector< double > LegCPsi;
    for (int l=0; l<=TMath::Max(2,degree); ++l) for (int m=0; m<=l; ++m) {
        costheta = BscosthetaMC->getValV();
        cospsi =BscospsiMC->getValV();
        phi = BsphiMC->getValV();
        LegCTheta.push_back(ROOT::Math::assoc_legendre(l,m,costheta));
        LegCPsi.push_back(ROOT::Math::assoc_legendre(l,m,cospsi));
    }
    double eff = 0;
    int i_ord = 0;
    int ctheta_ord, cpsi_ord, phi_ord;
    for (ctheta_ord=0; ctheta_ord<=degree; ++ctheta_ord)
        for (cpsi_ord=0; cpsi_ord<=degree; ++cpsi_ord)
            for (phi_ord=-1*TMath::Min(ctheta_ord,cpsi_ord); phi_ord<=TMath::Min(ctheta_ord,cpsi_ord); ++phi_ord) {
                double ff = ((RooAbsReal*)EffCoeff->at(i_ord))->getVal();
                // cout<<"++++++++++++++++++++++"<<ff<<"\n";
                if (phi_ord<0)      eff += ((RooAbsReal*)EffCoeff->at(i_ord))->getVal()*LegCTheta[cpsi_ord*(cpsi_ord+1)/2.0-phi_ord]*LegCPsi[ctheta_ord*(ctheta_ord+1)/2.0-phi_ord]*sin(-1.0*phi_ord*phi);
                
                else if (phi_ord>0) eff += ((RooAbsReal*)EffCoeff->at(i_ord))->getVal()*LegCPsi[cpsi_ord*(cpsi_ord+1)/2.0+phi_ord]*LegCTheta[ctheta_ord*(ctheta_ord+1)/2.0+phi_ord]*cos(phi_ord*phi);
                else eff += ((RooAbsReal*)EffCoeff->at(i_ord))->getVal()*LegCPsi[cpsi_ord*(cpsi_ord+1)/2.0] *LegCTheta[ctheta_ord*(ctheta_ord+1)/2.0];
                ++i_ord;
            }
    std::cout<<"Efficiency: "<<"===="<<eff<<"\n";

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
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

}
