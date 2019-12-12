#include "RooRealVar.h"
#include "RooDataSet.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TMath.h"
#include "RooPlot.h"
#include "TRandom3.h"
#include "TH3.h"
#include "TF3.h"
#include "TMath.h"
using namespace RooFit;
using namespace std;


Double_t g1(Double_t *x, Double_t *par) {
  Double_t thetf = x[0], psif = x[1], phif =x[2];//Double_t thetf = 1 - x[0]*x[0]; Double_t psif = 1- x[1]*x[1]; Double_t phif = x[2];
 
  Double_t func1 = par[0] + par[1]*(1-psif*psif)*(1-(1-thetf*thetf)*sin(phif)*sin(phif));
  Double_t func2 = par[2]*(1-psif*psif)*(1-thetf*thetf)+par[3]*pow(sin(phif),4)*pow(cos(thetf),4);
  Double_t func3 = par[4]*cos(psif)*cos(psif)*sqrt(1-cos(psif)*cos(psif))*(1-(1-cos(thetf)*cos(thetf))*sin(2*phif));
  Double_t func4 =par[5]*pow(sin(phif),4)+par[6]*(1-(1-cos(thetf)*cos(thetf))*cos(phif)*cos(phif));
  Double_t func = func1+func2+func3+func4;
  return func;
  
  
}



void angefffit (int maxOrder=7, int degree =7, int xbins=25, int ybins = 0, int zbins = 0)
{


    TF3 * angform= new TF3("angform",g1, -1,1,-1,1,-TMath::Pi(),TMath::Pi(),7);
    angform->SetParameter(0,0.06);  // par[0] : offset
    angform->SetParLimits(0,0.,2.0);
  	 angform->SetParName(0,"kBG0");
 	 angform->SetParameter(1,-0.01); // par[1] : width at f(x) minimum
 	 angform->SetParLimits(1,-1,1.);
 	 angform->SetParName(1,"kBG1");
 	 angform->SetParameter(2,0.007);  // par[2] : slope (negative side)
 	 angform->SetParLimits(2,-1,1.0);
 	 angform->SetParName(2,"kBG2");
 	 angform->SetParameter(3,-0.003);   // par[3] : inflection point abscissa
 	 angform->SetParLimits(3,-1,1.0);
 	 angform->SetParName(3,"kBG3");
 	 angform->SetParameter(4,-0.008);   // par[3] : inflection point abscissa
 	 angform->SetParLimits(4,-1.0,1.0);
 	 angform->SetParName(4,"kBG4");	  
 	 angform->SetParameter(5,-0.03); // par[1] : width at f(x) minimum
 	 angform->SetParLimits(5,-1.0,1.);
 	 angform->SetParName(5,"kBG5");
 	 angform->SetParameter(6,0.1);  // par[2] : slope (negative side)
 	 angform->SetParLimits(6,-1.0,1.);
 	 angform->SetParName(6,"kBG6");





if ( maxOrder < 0 ) return;
  if ( ybins<1 ) ybins = xbins;
  if ( zbins<1 ) zbins = xbins;
  if ( xbins<1 ) return;
RooRealVar *BscosthetaMC=new RooRealVar("BscosthetaMC","cos(#theta_{T})",-1,1);
RooRealVar *BscospsiMC=new RooRealVar("BscospsiMC","cos(#psi_{T})",-1,1);
RooRealVar *BsphiMC = new RooRealVar("BsphiMC","#phi",-TMath::Pi(),TMath::Pi());
RooArgSet vars(*BscosthetaMC, *BscospsiMC, *BsphiMC);
//RooRealVar wei ("wei","weight",-1,1);

  TChain* t_den = new TChain();
  TChain* t_num = new TChain();
    t_den->Add("/Users/ab/Documents/BsToJpsiPhi/Angle_Eff/BsGenV13.root/OutTreeGEN");
    t_num->Add("/Users/ab/Documents/BsToJpsiPhi/Angle_Eff/BsMC17_JpsiMu.root/fit");
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


//Prepare the date and assign to roorealVar
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
  
   TH3D* denHist = (TH3D*)recodata->createHistogram( "denHist", *BscosthetaMC,   Binning(xbins,-1,1) , YVar(*BscospsiMC,Binning(ybins,-1,1)), ZVar(*BsphiMC,Binning(zbins,-TMath::Pi(),TMath::Pi())) );
   TH3D* numHist = (TH3D*)numData->createHistogram("numHist", *BscosthetaMC,  Binning(xbins,-1,1) ,   YVar(*BscospsiMC,Binning(ybins,-1,1)),   ZVar(*BsphiMC,Binning(zbins,-TMath::Pi(),TMath::Pi())) );
 
 
 //Project and plot the angular efficiencies
// auto c1 = new TCanvas("c1", "c1", 0, 0, 1200, 800);
                     TH3D *hdivideIII = (TH3D*)numHist->Clone("hdivideIII");
                     hdivideIII->Sumw2();
                     hdivideIII->Divide(denHist);
                     hdivideIII->Fit(angform);
                     TFitResultPtr r =    hdivideIII->Fit(angform, "S");
                     TMatrixDSym cov = r->GetCovarianceMatrix();   //  to access the covariance matrix
                     Double_t chi2   = r->Chi2();                  // to retrieve the fit chi2
                     Double_t par0   = r->Parameter(0);            // retrieve the value for the parameter 0
                     Double_t err0   = r->ParError(0);             // retrieve the error for the parameter 0
                     r->Print("V");                                // print full information of fit including covariance matrix
                     r->Write();                                   // store the result in a file
                      Double_t angpar0,angpar1,angpar2, angpar3,angpar4, angpar5, angpar6;
                     auto proj2d1 =  hdivideIII->Project3D("xy");
                     auto proj2d2 =  hdivideIII->Project3D("xz");
                     auto proj2d3 =  hdivideIII->Project3D("yz");
                     proj2d1->SetTitle("Efficiency cos(#theta_{T})/cos(#psi_{T}) projection ");
                     proj2d2->SetTitle("Efficiency cos(#theta_{T})/#phi projection");
                     proj2d3->SetTitle("Efficiency cos(#psi_{T})/#phi projection");
                     auto cproj= new TCanvas("cproj", "cproj",1200,800);           
                     cproj->Divide(3,2);
                     cproj->cd(1); proj2d1->GetXaxis()->SetTitleOffset(1.4);proj2d1->GetYaxis()->SetTitleOffset(2);proj2d1->SetMinimum(0.0);proj2d1->Draw("SURF3");
                     cproj->cd(2);proj2d2->GetXaxis()->SetTitleOffset(1.4);proj2d2->GetYaxis()->SetTitleOffset(2);proj2d2->SetMinimum(0.0);proj2d2->Draw("SURF3");
                     cproj->cd(3);proj2d3->GetXaxis()->SetTitleOffset(1.4);proj2d3->GetYaxis()->SetTitleOffset(2);proj2d3->SetMinimum(0.0);proj2d3->Draw("SURF3");
                    // cproj->SaveAs(Form("EffProjectionHist_%i_%i_%i_2DProj_SpH%iOrder.pdf",xbins,ybins,zbins,maxOrder));
    angpar0 = angform->GetParameter(0);
    angpar1 = angform->GetParameter(1);
    angpar2 = angform->GetParameter(2);
    angpar3 = angform->GetParameter(3);
    angpar4 = angform->GetParameter(4);
    angpar5 = angform->GetParameter(5);
    angpar6 = angform->GetParameter(6);
    
    
    RooRealVar *angp0 = new RooRealVar("angp0","angp0", angpar0);
    RooRealVar *angp1 = new RooRealVar("angp1","angp1", angpar1);
    RooRealVar *angp2 = new RooRealVar("angp2","angp2", angpar2);
    RooRealVar *angp3 = new RooRealVar("angp3","angp3", angpar3);
    RooRealVar *angp4 = new RooRealVar("angp4","angp4", angpar4);
    RooRealVar *angp5 = new RooRealVar("angp5","angp5", angpar5);
    RooRealVar *angp6 = new RooRealVar("angp6","angp6", angpar6);
    
                     
    RooFormulaVar *funang1 = new RooFormulaVar("funang1", "angp0+angp1*(1-@2*@2)*(1-(1-@3*@3)*sin(@4)*sin(@4))", RooArgSet(*angp0,*angp1,*BscospsiMC,*BscosthetaMC,*BsphiMC));
    RooFormulaVar *funang2 = new RooFormulaVar("funang2", "angp2*(1-@1*@1)*(1-@2*@2)+@3*pow(sin(@4),4)*pow(cos(@2),4)",RooArgSet(*angp2,*BscospsiMC,*BscosthetaMC,*angp3,*BsphiMC));
    RooFormulaVar *funang3 = new RooFormulaVar("funang3", "angp4*cos(@1)*cos(@1)*sqrt(1-cos(@1)*cos(@1))*(1-(1-cos(@2)*cos(@2))*sin(2*@3))", RooArgSet(*angp4,*BscospsiMC,*BscosthetaMC,*BsphiMC));
    RooFormulaVar *funang4 = new RooFormulaVar("funang4", "angp5*pow(sin(@3),4)+angp6*(1-(1-cos(@2)*cos(@2))*cos(@3)*cos(@3))",RooArgSet(*angp5,*angp6,*BscosthetaMC,*BsphiMC));
  
    RooAddition *angTotFun = new RooAddition("angTotFun","3DAngEffFunction", RooArgSet(*funang1,*funang2,*funang3,*funang4));
  
    // Plot efficiency projections
    gStyle->SetOptStat(0);
    auto h3_xyz = (TH3D*)angTotFun->createHistogram("BscosthetaMC,BscospsiMC,BsphiMC",50,50,50);
    auto h3_x  = h3_xyz->ProjectionX();
    auto h3_y  = h3_xyz->ProjectionY();
    auto h3_z  = h3_xyz->ProjectionZ();
    auto h3_xy = h3_xyz->Project3D("xy");
    auto h3_xz = h3_xyz->Project3D("xz");
    auto h3_yz = h3_xyz->Project3D("yz");
    
    
    // 2D projections
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
    cproj->cd(4);
    h3_xy->Draw("SURF3");
    cproj->cd(5);
    h3_xz->Draw("SURF3");
    cproj->cd(6);
    h3_yz->Draw("SURF3");
    
    
    
  /*                   
                     
                     
                      TH1* projcthet =  hdivideIII->Project3D("x");
                  RooRealVar *ctheteff = new RooRealVar("ctheteff", "ctheteff", -1,1);
                  RooDataHist *cthethist = new RooDataHist("cthethist","cthethist",*ctheteff,Import(*projcthet)) ;
                  RooPlot* cthetfr = ctheteff->frame(Title("Binned cos#theta efficiency"),Bins(25)) ;
                  cthethist->plotOn(cthetfr) ; 
 
//=======================================================================================================================================Creat the function
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
  RooRealSumPdf *projectedFunc = new RooRealSumPdf("projectedFunc", "projectedFunc", funList, facList );
  
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
       if( factors[iOrder]->getValV() < 0.0009)continue;
	cout<<xOrder<<" "<<yOrder<<" "<<zOrder<<"\t"<<iOrder<<"\t"<<factors[iOrder]->getValV()<<endl;

      }
      
  t.Stop();
  t.Print();    
        RooWorkspace *ws = new RooWorkspace("ws","Workspace with efficiency parameterisation");
        ws->import( *projectedFunc, Silence() );                          
             gStyle->SetOptStat(0);
  auto h3_xyz = (TH3D*)projectedFunc->createHistogram("BscosthetaMC,BscospsiMC,BsphiMC",50,50,50);
  auto h3_x  = h3_xyz->ProjectionX();
  auto h3_y  = h3_xyz->ProjectionY();
  auto h3_z  = h3_xyz->ProjectionZ();
  auto h3_xy = h3_xyz->Project3D("xy");
  auto h3_xz = h3_xyz->Project3D("xz");
  auto h3_yz = h3_xyz->Project3D("yz");
  
 
  
 
    RooRealVar *cthetefffunc=new RooRealVar("cthetefffunc","CosTheta;#cos#theta_{T} (a.u.); #epsilon(cos#theta)", -1,1) ;
    RooDataHist *cthetfunc = new RooDataHist("dh","dh",*cthetefffunc,Import(*h3_x)) ;
    RooHistPdf *ctheteffpdf = new RooHistPdf("ctheteffpdf","effhistpdf",*cthetefffunc, *cthetfunc, 0) ;
           
   //ctheteffpdf->fitTo(*cthethist, Save(kTRUE));    
   ctheteffpdf->plotOn(cthetfr) ;
   ctheteffpdf->paramOn(cthetfr);
   ctheteffpdf->plotOn(cthetfr, LineColor(kBlue), LineWidth(2));
   //Double_t chisquare = cthetfr->chiSquare();
   //cout<<"Chi square of fit is :"<< chisquare<< endl;      
   cthetfr->Draw();    
    TCanvas *cc33 = new TCanvas("cc33", "cc33",0,0,800,600);
               h3_x->Draw();  
               */
               
} 


