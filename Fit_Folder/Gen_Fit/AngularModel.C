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
#include "RooTruthModel.h"
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

void AngularModel(string filename) {
  
             TChain* chain_data = new TChain("gen");
             chain_data->Add(filename.c_str());
             Int_t nevt = (int)chain_data->GetEntries();

    RooRealVar *BsCt2DMC = new RooRealVar ("BsCt2DMC","Bs ct", 0.,0.3,"cm");//
    RooRealVar *BscosthetaMC= new RooRealVar ("BscosthetaMC","cos(#theta_{T})", -1,1);
    RooRealVar *BscospsiMC= new RooRealVar ("BscospsiMC","cos(#psi_{T})", -1,1);
    RooRealVar *BsphiMC= new RooRealVar ("BsphiMC","#phi_{T}", -TMath::Pi(),TMath::Pi(),"rad");//cosdelta2
    double tagmis = 1;
    RooRealVar *mistag = new RooRealVar("mistag","Mistag fraction of original B and Bbar",tagmis);
    RooCategory *tag = new RooCategory("tag","Flavour tag of the B meson");
    tag->defineType("Bsbar",-1);
    tag->defineType("Bs",+1);
    tag->defineType("untag",0);

    RooDataSet* data = new RooDataSet("data", "raw data1", RooArgSet(*BscosthetaMC, *BscospsiMC, *BsphiMC, *BsCt2DMC,*tag), Import(*chain_data));

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
    RooRealVar *deltaPe=new RooRealVar("deltaPe","#deltaPe",0.17,-0.5,.5);//,-1,4); // =0.17  
    RooRealVar *deltaSPe=new RooRealVar("deltaSPe","#deltaPe",-0.028,-3.,3.);//,-1,4); // =0.17  
    RooRealVar *dm=new RooRealVar ("dm","dm",589.7,500,700); //590.1// in um*c, PDG value 17.69 in 1/ps
    RooRealVar *dGam=new RooRealVar ("dGam","dGam",0,-3.,6.);//2.3//original value in MC
    RooRealVar *tau=new RooRealVar ("tau","tau",0.0441,0.03,0.055,"ps");//,0.043//origunal value in MC

  

    RooTruthModel truth("truth","truth",*BsCt2DMC ); 

    RooFormulaVar fsinh1("fsinh1","fsinh1","-cos(@0)",RooArgList(*phi_s));
    RooFormulaVar fsin1bar("fsin1bar","fsin1bar","-sin(@0)",RooArgList(*phi_s));
    RooFormulaVar fsin1("fsin1","fsin1","sin(@0)",RooArgList(*phi_s));
    RooFormulaVar fsin1tag("fsin1tag","fsin1tag","-(1-2*@2)*@1*sin(@0)",RooArgList(*phi_s,*tag,*mistag));
    
    RooFormulaVar fsinh2("fsinh2","fsinh2","cos(@0)",RooArgList(*phi_s));
    RooFormulaVar fsin2bar("fsin2bar","fsin2bar","sin(@0)",RooArgList(*phi_s));
    RooFormulaVar fsin2("fsin2","fsin2","-sin(@0)",RooArgList(*phi_s));
    RooFormulaVar fsin2tag("fsin2tag","fsin2tag","(1-2*@2)*@1*sin(@0)",RooArgList(*phi_s,*tag,*mistag));

    RooFormulaVar fsinh4("fsinh4","fsinh4","-sin(@0)*cos(@1-@2)/sin(@1-@2)",RooArgList(*phi_s,*deltaPe,*deltaPa));
    RooFormulaVar fsin4bar("fsin4bar","fsin4bar","cos(@0)*cos(@1-@2)",RooArgList(*phi_s,*deltaPe,*deltaPa));
    RooFormulaVar fcos4bar("fcos4bar","fcos4bar","-sin(@0-@1)",RooArgList(*deltaPe,*deltaPa));
    RooFormulaVar fsin4("fsin4","fsin4","-cos(@0)*cos(@1-@2)",RooArgList(*phi_s,*deltaPe,*deltaPa));
    RooFormulaVar fcos4("fcos4","fcos4","sin(@0-@1)",RooArgList(*deltaPe,*deltaPa));
    RooFormulaVar fsin4tag("fsin4tag","fsin4tag","(1-2*@4)*@3*cos(@0)*cos(@1-@2)/sin(@1-@2)",RooArgList(*phi_s,*deltaPe,*deltaPa,*tag,*mistag));
    RooFormulaVar fcos4tag("fcos4tag","fcos4tag","-(1-2*@1)*@0",RooArgList(*tag,*mistag));
    
    RooFormulaVar fsinh5("fsinh5","fsinh5","-cos(@0)",RooArgList(*phi_s));
    RooFormulaVar fcosh5("fcosh5","fcosh5","cos(@0)",RooArgList(*deltaPa));
    RooFormulaVar fsin5bar("fsin5bar","fsin5bar","-sin(@0)*cos(@1)",RooArgList(*phi_s,*deltaPa));
    RooFormulaVar fsin5("fsin5","fsin5","sin(@0)*cos(@1)",RooArgList(*phi_s,*deltaPa));
    RooFormulaVar fsin5tag("fsin5tag","fsin5tag","-(1-2*@2)*@1*sin(@0)",RooArgList(*phi_s,*tag,*mistag));

    RooFormulaVar fsinh6("fsinh6","fsinh6","-sin(@0)*cos(@1)/sin(@1)",RooArgList(*phi_s,*deltaPe));
    RooFormulaVar fsin6bar("fsin6bar","fsin6bar","cos(@0)*cos(@1)",RooArgList(*phi_s,*deltaPe));
    RooFormulaVar fcos6bar("fcos6bar","fcos6bar","-sin(@0)",RooArgList(*deltaPe));
    RooFormulaVar fsin6("fsin6","fsin6","-cos(@0)*cos(@1)",RooArgList(*phi_s,*deltaPe));
    RooFormulaVar fcos6("fcos6","fcos6","sin(@0)",RooArgList(*deltaPe));
    RooFormulaVar fsin6tag("fsin6tag","fsin6tag","(1-2*@3)*@2*cos(@0)*cos(@1)/sin(@1)",RooArgList(*phi_s,*deltaPe,*tag,*mistag));
    RooFormulaVar fcos6tag("fcos6tag","fcos6tag","-(1-2*@1)*@0",RooArgList(*tag,*mistag));



    RooFormulaVar coshGBasis("coshGBasis","exp(-@0/@1)*cosh(@0*@2/2)", RooArgList(*BsCt2DMC ,*tau,*dGam));
    RooFormulaVar sinhGBasis("sinhGBasis","exp(-@0/@1)*sinh(@0*@2/2)", RooArgList(*BsCt2DMC ,*tau,*dGam));
    RooFormulaVar cosGBasis("cosGBasis","exp(-@0/@1)*cos(@0*@2)",RooArgList(*BsCt2DMC ,*tau,*dm));
    RooFormulaVar sinGBasis("sinGBasis","exp(-@0/@1)*sin(@0*@2)",RooArgList(*BsCt2DMC ,*tau,*dm));

    RooAbsReal* coshGConv = truth.convolution(&coshGBasis,BsCt2DMC );
    RooAbsReal* sinhGConv = truth.convolution(&sinhGBasis,BsCt2DMC );
    RooAbsReal* cosGConv = truth.convolution(&cosGBasis,BsCt2DMC );
    RooAbsReal* sinGConv = truth.convolution(&sinGBasis,BsCt2DMC );




    RooAddition *myAmp0tag=new RooAddition("myAmp0tag","myAmp0tag",RooArgList(*sinhGConv,*coshGConv,*sinGConv),RooArgList(fsinh1,RooConst(1),fsin1tag));
    RooAddition *myAmpetag=new RooAddition("myAmpetag","myAmpetag",RooArgList(*sinhGConv,*coshGConv,*sinGConv),RooArgList(fsinh2,RooConst(1),fsin2tag));
    RooAddition *myAmpa4tag=new RooAddition("myAmpa4tag","myAmpa4tag",RooArgList(*sinhGConv,*cosGConv,*sinGConv),RooArgList(fsinh4,fcos4tag,fsin4tag));
    RooAddition *myAmpa5tag=new RooAddition("myAmpa5tag","myAmpa5tag",RooArgList(*sinhGConv,*coshGConv,*sinGConv),RooArgList(fsinh5,RooConst(1),fsin5tag));
    RooAddition *myAmpa6tag=new RooAddition("myAmpa6tag","myAmpa6tag",RooArgList(*sinhGConv,*cosGConv,*sinGConv),RooArgList(fsinh6,fcos6tag,fsin6tag));



    RooFormulaVar *cosdpa=new RooFormulaVar("cosdpa","cosdpa","cos(@0)",RooArgList(*deltaPa));
    RooFormulaVar *sindpe=new RooFormulaVar("sindpe","sindpe","sin(@0)",RooArgList(*deltaPe));
    RooFormulaVar *sindpadpe=new RooFormulaVar("sindpadpe","sindpadpe","sin(@0-@1)",RooArgList(*deltaPe,*deltaPa));



    RooFormulaVar phinew("phinew","BsphiMC/TMath::Pi()",*BsphiMC);
    RooLegendre phizero("phizero","phizero",phinew,0,0);
    RooFormulaVar* fun1T = new RooFormulaVar("fun1T","2*BscospsiMC*BscospsiMC*(1-(1-BscosthetaMC*BscosthetaMC)*cos(BsphiMC)*cos(BsphiMC))",RooArgSet(*BscosthetaMC,*BsphiMC,*BscospsiMC));
    RooFormulaVar* fun2T = new RooFormulaVar("fun2T","(1-BscospsiMC*BscospsiMC)*(1-(1-BscosthetaMC*BscosthetaMC)*sin(BsphiMC)*sin(BsphiMC))",RooArgSet(*BscosthetaMC,*BsphiMC,*BscospsiMC));
    RooFormulaVar* fun3T = new RooFormulaVar("fun3T","(1-BscospsiMC*BscospsiMC)*(1-BscosthetaMC*BscosthetaMC)*phizero",RooArgSet(*BscosthetaMC,*BscospsiMC,phizero));
    RooFormulaVar* fun4T = new RooFormulaVar("fun4T","-(1-BscospsiMC*BscospsiMC)*2*BscosthetaMC*sqrt(1-BscosthetaMC*BscosthetaMC)*sin(BsphiMC)",RooArgSet(*BscosthetaMC,*BsphiMC,*BscospsiMC));
    RooFormulaVar* fun5T = new RooFormulaVar("fun5T","2/sqrt(2.)*BscospsiMC*sqrt(1-BscospsiMC*BscospsiMC)*(1-BscosthetaMC*BscosthetaMC)*sin(2*BsphiMC)",RooArgSet(*BscosthetaMC,*BsphiMC,*BscospsiMC));
    RooFormulaVar* fun6T = new RooFormulaVar("fun6T","2/sqrt(2.)*BscospsiMC*sqrt(1-BscospsiMC*BscospsiMC)*2*BscosthetaMC*sqrt(1-BscosthetaMC*BscosthetaMC)*cos(BsphiMC)",RooArgSet(*BscosthetaMC,*BsphiMC,*BscospsiMC));



    RooFormulaVar tagfun("tagfun","tagfum","(1-@0*(1-2*@1))/2",RooArgList(*tag,*mistag));
    RooFormulaVar tagfunbar("tagfunbar","tagfumbar","(1+@0*(1-2*@1))/2",RooArgList(*tag,*mistag));

    RooProduct *pdf1tag=new RooProduct("pdf1tag","amp0*g1",RooArgSet(*myAmp0tag,*fun1T));//,*eff_ct));// eff_ct + *poly3D
    RooProduct *pdf2tag=new RooProduct("pdf2tag","amp0*g2",RooArgSet(*myAmp0tag,*fun2T));//,*eff_ct));//
    RooProduct *pdf3tag=new RooProduct("pdf3tag","ampe*g3",RooArgSet(*myAmpetag,*fun3T));//,*eff_ct));//
    RooProduct *pdf4tag=new RooProduct("pdf4tag","Ampa4*fun4T",RooArgSet(*myAmpa4tag,*fun4T,*sindpadpe));//,*eff_ct));   
    RooProduct *pdf5tag=new RooProduct("pdf5tag","ampa0*g5",RooArgSet(*myAmpa5tag,*fun5T,*cosdpa));//,*eff_ct));//
    RooProduct *pdf6tag=new RooProduct("pdf6tag","Ampa6*fun6T",RooArgSet(*myAmpa6tag,*fun6T,*sindpe));//,*eff_ct)); 


   // RooRealSumPdf *PDFtag=new RooRealSumPdf("PDFtag","Signal PDF Amp_i*funiT",RooArgList(*pdf1tag,*pdf2tag,*pdf3tag,*pdf4tag,*pdf5tag,*pdf6tag),RooArgList(*A_0,*A_pa,*A_pe,*A_4,*A_5,*A_6));//,Conditional(gm,*BdCt2DBS)); 
               RooRealSumPdf *Angular_Model=new RooRealSumPdf("Angular_Model","Signal PDF Amp_i*funiT",RooArgList(*pdf1tag,*pdf2tag,*pdf3tag,*pdf4tag,*pdf5tag,*pdf6tag),RooArgList(*A_0,*A_pa,*A_pe,*A_4,*A_5,*A_6));//,Conditional(gm,*BdCt2DBS));
               RooFitResult* fitRes = Angular_Model->fitTo(*data, NumCPU(8), Save(kTRUE));
               fitRes->Print("v");
               TString title1 = "cos#theta_{T}";
               TString title2 = "cos#psi_{T}";
               TString title3 = "#phi_{T}";
               TString title4 = "B_{s} ct";               
               RooPlot* CosTheta = BscosthetaMC->frame(Title(title1),Bins(80));
               data->plotOn(CosTheta,DataError(RooAbsData::SumW2));
               Angular_Model->plotOn(CosTheta) ;
               Angular_Model->paramOn(CosTheta);
               Angular_Model->plotOn(CosTheta, LineColor(kBlack), LineWidth(1));
               Double_t chisquare_costheta = CosTheta->chiSquare();
               cout<<"Chi square of costheta  fit is :"<< chisquare_costheta<<"\n";
               
               RooPlot* CosPsi = BscospsiMC->frame(Title(title2),Bins(80));
               data->plotOn(CosPsi,DataError(RooAbsData::SumW2));
               Angular_Model->plotOn(CosPsi);
               Angular_Model->paramOn(CosPsi);
               Angular_Model->plotOn(CosPsi, LineColor(kBlack), LineWidth(1));
               Double_t chisquare_cospsi = CosPsi->chiSquare();
               cout<<"Chi square of cospsi  fit is :"<< chisquare_cospsi<<"\n"; 
               
               RooPlot* PHI = BsphiMC->frame(Title(title3),Bins(80));
               data->plotOn(PHI, DataError(RooAbsData::SumW2));
               Angular_Model->plotOn(PHI);
               Angular_Model->paramOn(PHI);
               Angular_Model->plotOn(PHI , LineColor(kBlack), LineWidth(1));
               Double_t chisquare_phi = PHI->chiSquare();
               cout<<"Chi square of phi  fit is :"<< chisquare_phi<<"\n";
                
                RooPlot* ct = BsCt2DMC->frame(Title(title4),Bins(80));
                data->plotOn(ct, DataError(RooAbsData::SumW2));
                Angular_Model->plotOn(ct);
                Angular_Model->paramOn(ct);
                Angular_Model->plotOn(ct , LineColor(kBlack), LineWidth(1));
               Double_t chisquare_ct = ct->chiSquare();
               cout<<"Chi square of lifetime  fit is :"<< chisquare_ct<<"\n";


                 TCanvas *cc = new TCanvas("cc", "cc",0,0,800,600);
                CosTheta->Draw();
                cc->Print("ctheta.root", "root");
               TCanvas *cc1 = new TCanvas("cc1", "cc1",0,0,800,600);
               CosPsi->Draw();
               cc1->Print("cpsi.root", "root");
               TCanvas *cc2 = new TCanvas("cc2", "cc2",0,0,800,600);
               PHI->Draw();
               cc2->Print("phit.root", "root");
               TCanvas *cc3 = new TCanvas("cc3", "cc3",0,0,800,600);
               ct->Draw();
               cc3->Print("ct.root", "root");
               Angular_Model->graphVizTree ("modelScheme.dot") ;


  /*TCanvas c1 ;
  RooPlot * frame = BscosthetaMC->frame () ; 
 Angular_Model->plotOn (frame, LineColor (kRed + 2)) ;  
 Angular_Model->plotOn (frame, LineColor (kBlue + 2)) ;
  frame->Draw () ;	
  c1.Print ("modelPlot.gif","gif") ;
    */
    }
