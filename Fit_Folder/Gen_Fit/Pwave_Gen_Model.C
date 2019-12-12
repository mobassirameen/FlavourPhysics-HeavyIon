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





void Pwave_Gen_Model(string filename){
//==========================================================================================================Efficiencies



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
    tag->defineType("Bs",1);
    tag->defineType("untag",0);

    RooDataSet* data = new RooDataSet("data", "raw data1", RooArgSet(*BscosthetaMC, *BscospsiMC, *BsphiMC, *BsCt2DMC,*tag), Import(*chain_data));

    RooRealVar *A_0=new RooRealVar ("A_0","|A_0|^2|",0.6,0.30, 0.70);//0.51,0.30,0.7);
    RooRealVar *A_S=new RooRealVar ("A_S","|A_0|^2|",0.05,0.,0.1);
    RooRealVar *A_pe=new RooRealVar ("A_pe","|A_pe|^2",0.16,0.0,0.3);
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
    RooRealVar *deltaPe=new RooRealVar("deltaPe","#deltaPe",0.17,-0.5,.50);//,-1,4); // =0.17  
    RooRealVar *deltaSPe=new RooRealVar("deltaSPe","#deltaPe",-0.028,-3.,3.);//,-1,4); // =0.17  
    RooRealVar *dm=new RooRealVar ("dm","dm",589.7,500,700); //590.1// in um*c, PDG value 17.69 in 1/ps
    RooRealVar *dGam=new RooRealVar ("dGam","dGam",0,-3.,6.);//2.3//original value in MC
    RooRealVar *tau=new RooRealVar ("tau","tau",0.0441,0.03,0.055,"ps");//,0.043//origunal value in MC


    RooTruthModel tm("tm","truth",*BsCt2DMC );
 
    



    RooFormulaVar coshGBasis("coshGBasis","exp(-@0/@1)*cosh(@0*@2/2)", RooArgList(*BsCt2DMC ,*tau,*dGam));
    RooFormulaVar sinhGBasis("sinhGBasis","exp(-@0/@1)*sinh(@0*@2/2)", RooArgList(*BsCt2DMC ,*tau,*dGam));
    RooFormulaVar cosGBasis("cosGBasis","exp(-@0/@1)*cos(@0*@2)",RooArgList(*BsCt2DMC ,*tau,*dm));
    RooFormulaVar sinGBasis("sinGBasis","exp(-@0/@1)*sin(@0*@2)",RooArgList(*BsCt2DMC ,*tau,*dm));
    
    
      RooAbsReal* coshGConv = tm.convolution(&coshGBasis,BsCt2DMC);
      RooAbsReal* sinhGConv = tm.convolution(&sinhGBasis,BsCt2DMC);
      RooAbsReal* cosGConv = tm.convolution(&cosGBasis,BsCt2DMC);
      RooAbsReal* sinGConv = tm.convolution(&sinGBasis,BsCt2DMC);
    
   
    

   


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
 
    
    
    
    //RooFormulaVar *coef11=new RooFormulaVar("coef11","coef11","@0",RooArgList(*A_0));
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

    
    
    
    
    RooProduct *eq11=new RooProduct("eq11","amp0*g1",RooArgSet(*coshGConv,*fun1T));
    RooProduct *eq12=new RooProduct("eq12","amp0*g1",RooArgSet(*coef12,*sinhGConv,*fun1T));
    RooProduct *eq14=new RooProduct("eq14","amp0*g1",RooArgSet(*coef14,*sinGConv,*fun1T));
    RooProduct *eq21=new RooProduct("eq21","amp0*g1",RooArgSet(*coshGConv,*fun2T));
    RooProduct *eq22=new RooProduct("eq22","amp0*g1",RooArgSet(*coef22,*sinhGConv,*fun2T));
    RooProduct *eq24=new RooProduct("eq24","amp0*g1",RooArgSet(*coef24,*sinGConv,*fun2T));
    RooProduct *eq31=new RooProduct("eq31","amp0*g1",RooArgSet(*coshGConv,*fun3T));
    RooProduct *eq32=new RooProduct("eq32","amp0*g1",RooArgSet(*coef32,*sinhGConv,*fun3T));
    RooProduct *eq34=new RooProduct("eq34","amp0*g1",RooArgSet(*coef34,*sinGConv,*fun3T));
    RooProduct *eq42=new RooProduct("eq42","amp0*g1",RooArgSet(*coef42,*sinhGConv,*fun4T));
    RooProduct *eq43=new RooProduct("eq43","amp0*g1",RooArgSet(*coef43,*cosGConv,*fun4T));
    RooProduct *eq44=new RooProduct("eq44","amp0*g1",RooArgSet(*coef44,*sinGConv,*fun4T));
    RooProduct *eq51=new RooProduct("eq51","amp0*g1",RooArgSet(*coef51,*coshGConv,*fun5T));
    RooProduct *eq52=new RooProduct("eq52","amp0*g1",RooArgSet(*coef52,*sinhGConv,*fun5T));
    RooProduct *eq54=new RooProduct("eq54","amp0*g1",RooArgSet(*coef54,*sinGConv,*fun5T));
    RooProduct *eq62=new RooProduct("eq62","amp0*g1",RooArgSet(*coef62,*sinhGConv,*fun6T));
    RooProduct *eq63=new RooProduct("eq63","amp0*g1",RooArgSet(*coef63,*cosGConv,*fun6T));
    RooProduct *eq64=new RooProduct("eq64","amp0*g1",RooArgSet(*coef64,*sinGConv,*fun6T));
    RooProduct *eq71=new RooProduct("eq71","amp0*g1",RooArgSet(*coshGConv,*fun7T));
    RooProduct *eq72=new RooProduct("eq72","amp0*g1",RooArgSet(*coef72,*sinhGConv,*fun7T));
    RooProduct *eq74=new RooProduct("eq74","amp0*g1",RooArgSet(*coef74,*sinGConv,*fun7T));
    RooProduct *eq82=new RooProduct("eq82","amp0*g1",RooArgSet(*coef82,*sinhGConv,*fun8T));
    RooProduct *eq83=new RooProduct("eq83","amp0*g1",RooArgSet(*coef83,*cosGConv,*fun8T));
    RooProduct *eq84=new RooProduct("eq84","amp0*g1",RooArgSet(*coef84,*sinGConv,*fun8T));
    RooProduct *eq91=new RooProduct("eq91","amp0*g1",RooArgSet(*coef91,*coshGConv,*fun9T));
    RooProduct *eq92=new RooProduct("eq92","amp0*g1",RooArgSet(*coef92,*sinhGConv,*fun9T));
    RooProduct *eq94=new RooProduct("eq94","amp0*g1",RooArgSet(*coef94,*sinGConv,*fun9T));
    RooProduct *eq102=new RooProduct("eq102","amp0*g1",RooArgSet(*coef102,*sinhGConv,*fun10T));
    RooProduct *eq103=new RooProduct("eq103","amp0*g1",RooArgSet(*coef103,*cosGConv,*fun10T));
    RooProduct *eq104=new RooProduct("eq104","amp0*g1",RooArgSet(*coef104,*sinGConv,*fun10T));
    
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
    //RooRealSumPdf *PDFdef=new RooRealSumPdf("PDFdef","Signal PDF Amp_i*funiT",funtot,coeftot);
    RooRealSumPdf *Angular_Model=new RooRealSumPdf("Angular_Model","Signal PDF Amp_i*funiT",funtot,coeftot);
   

    //RooProdPdf *Angular_Model=new RooProdPdf("Angular_Model","PDF_Signaltag",RooArgList(*MTpdf,*PDFdef),Conditional(*PDFdef,RooArgList(*BsCt2DMC,*BscosthetaMC,*BscospsiMC,*BsphiMC,*tag,*mistag)));

    //RooRealSumPdf *PDFtag=new RooRealSumPdf("PDFtag","Signal PDF Amp_i*funiT",RooArgList(*pdf1tag,*pdf2tag,*pdf3tag,*pdf4tag,*pdf5tag,*pdf6tag),RooArgList(*A_0,*A_pa,*A_pe,*A_4,*A_5,*A_6));//,Conditional(gm,*BdCt2DBS));
   // RooRealSumPdf *Angular_Model=new RooRealSumPdf("Angular_Model","Signal PDF Amp_i*funiT",RooArgList(*pdf1tag,*pdf2tag,*pdf3tag,*pdf4tag,*pdf5tag,*pdf6tag),RooArgList(*A_0,*A_pa,*A_pe,*A_4,*A_5,*A_6));//,Conditional(gm,*BdCt2DBS));
               //RooProdPdf *Angular_Model = new RooProdPdf("Angular_Model", "Angular_Model", *PDFtag, *T_gaus);
               RooFitResult* PwaveMCResults = Angular_Model->fitTo(*data, NumCPU(8), Save(kTRUE));//, ConditionalObservables(*BsCt2DMCErr));
               PwaveMCResults->Print("v");
    
    // Construct 2D color plot of correlation matrix
    gStyle->SetOptStat(0) ;
    gStyle->SetPalette(1) ;
    TH2* hcorr = PwaveMCResults->correlationHist() ;
    TCanvas* c = new TCanvas("Correlation Matrix","Correlation Matrix",1000,600) ;
    hcorr->GetYaxis()->SetTitleOffset(1.4) ; hcorr->Draw("colz");
               /*TString title1 = "cos#theta_{T}";
               TString title2 = "cos#psi_{T}";
               TString title3 = "#phi_{T}";    
               TString title4 = "ct_{B_{s}} (cm)";
               RooPlot* CosTheta = BscosthetaMC->frame(Title(title1),Bins(100));
               data->plotOn(CosTheta,DataError(RooAbsData::SumW2));
               Angular_Model->plotOn(CosTheta) ;
               Angular_Model->paramOn(CosTheta);
               Angular_Model->plotOn(CosTheta, LineColor(kBlack), LineWidth(1));
               Double_t chisquare_costheta = CosTheta->chiSquare();
    
               
               RooPlot* CosPsi = BscospsiMC->frame(Title(title2),Bins(100));
               data->plotOn(CosPsi,DataError(RooAbsData::SumW2));
               Angular_Model->plotOn(CosPsi);
               Angular_Model->paramOn(CosPsi);
               Angular_Model->plotOn(CosPsi, LineColor(kBlack), LineWidth(1));
               Double_t chisquare_cospsi = CosPsi->chiSquare();
    
               
               RooPlot* PHI = BsphiMC->frame(Title(title3),Bins(100));
               data->plotOn(PHI, DataError(RooAbsData::SumW2));
               Angular_Model->plotOn(PHI);
               Angular_Model->paramOn(PHI);
               Angular_Model->plotOn(PHI , LineColor(kBlack), LineWidth(1));
               Double_t chisquare_phi = PHI->chiSquare();
              
                
                RooPlot* ct = BsCt2DMC->frame(Title(title4),Bins(100));
                data->plotOn(ct, DataError(RooAbsData::SumW2));
                Angular_Model->plotOn(ct);
                Angular_Model->paramOn(ct);
                Angular_Model->plotOn(ct , LineColor(kBlack), LineWidth(1));
                Double_t chisquare_ct = ct->chiSquare();

               
   

               BsCt2DMC->setRange("Asymme",0.02,0.08) ;
               RooPlot* asymplot=BsCt2DMC->frame(Title("BsCt Asymetry"),Range("Asymme"),Binning(100));
               data->plotOn(asymplot,LineStyle(kDashed),Binning(100),RooFit::Name("Bs ct Asym"),Asymmetry(*tag),Range("Asymme"));
               Angular_Model->plotOn(asymplot,RooFit::Name("Fitted ct Asym"),Asymmetry(*tag),Range("Asymme"));
 
               RooPlot* framepull = BsCt2DMC->frame(RooFit::Title("Ct pull"));
               RooHist* hpull = ct->pullHist() ;
               framepull->addPlotable(hpull,"P0") ;
               framepull->SetMinimum(-5) ;
               framepull->SetMaximum(+5) ;
               framepull->SetYTitle("pull");
               framepull->SetMarkerStyle(20);               
               framepull->SetNdivisions(10);

               cout<<"Num of fit param P wave reco results:"<<PwaveMCResults->floatParsFinal().getSize()<<endl;
               cout<<"Chi2 ct :"<< chisquare_ct<<"\n";
               cout<<"Chi2  costheta :"<< chisquare_costheta<<"\n";
               cout<<"Chi2  cospsi :"<< chisquare_cospsi<<"\n";
               cout<<"Chi2  phi :"<< chisquare_phi<<"\n";

               TCanvas *cc = new TCanvas("cc", "cc",0,0,800,600);
               CosTheta->Draw();
               //cc->Print("Gen_ctheta.root", "root");
               TCanvas *cc11 = new TCanvas("cc11", "cc11",0,0,800,600);
               CosPsi->Draw();
               //cc11->Print("Gen_cospsi.root", "root");
               TCanvas *cc22 = new TCanvas("cc22", "cc22",0,0,800,600);
               PHI->Draw();
               //cc22->Print("Gen_phit.root", "root");
               TCanvas *cc33 = new TCanvas("cc33", "cc33",0,0,800,600);
               ct->Draw();
               //cc33->Print("Gen_ct.root", "root");
               TCanvas *cc44 = new TCanvas("cc44", "cc44",0,0,800,600);
               asymplot->Draw();
               //cc44->Print("Gen_Asymmetry.root", "root");
               TCanvas *cc55 = new TCanvas("cc55", "cc55",0,0,800,600);
               framepull->Draw();
              // cc55->Print("Gen_framepull.root", "root");
             
               TFile *output = new TFile("fitresults_gen.root","recreate");
               PwaveMCResults->Write();
               asymplot->Write("AsymmetryPlot");
               framepull->Write("LifetimePull") ;
               CosTheta->Write("Costhetafit");
               CosPsi->Write("Cospsifit");
               PHI->Write("Phifit");
               ct->Write("Lifetimefit");
               hcorr->Write("CorrelationMatrix");
              
               output->Close();

               Angular_Model->graphVizTree ("genfit_results.dot") ;
*/

  /*TCanvas c1 ;
  RooPlot * frame = BscosthetaMC->frame () ; 
 Angular_Model->plotOn (frame, LineColor (kRed + 2)) ;  
 Angular_Model->plotOn (frame, LineColor (kBlue + 2)) ;
  frame->Draw () ;	
  c1.Print ("modelPlot.gif","gif") ;
    */
    }
