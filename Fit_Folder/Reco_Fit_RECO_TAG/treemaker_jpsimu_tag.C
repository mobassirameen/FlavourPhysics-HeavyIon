#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <iostream>
#include <fstream>
#include <TH2.h>
#include <TH1.h>
#include <string>
#include <TMath.h>
#include <vector>

#include <TStyle.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TBuffer.h>

using namespace std;

void treemaker_jpsimu_tag(){
   // TFile *fIn1 = new TFile("/home/cms/BMeson/FileRepository/Generalplots/BsJpsiPhi_dG0.root");
    //TFile *fIn1 = new TFile("/home/cms/BMeson/ffrepo/ntuple/Tot_2017MC_withTag.root");
    TFile *fIn1 = new TFile("BsMC17V13_ter.root");
  //  TFile *fIn1 = new TFile("B0sample.root");
    if (!fIn1){return;}
    TTree* bs = (TTree*)fIn1->Get("OutTree");
    Int_t n_entries = bs->GetEntries();
    std::cout<<n_entries<<"\n";
    Float_t  isBs, HLT_JpsiMu, HLT_JpsiTkTk, angle_cospsi, angle_costheta, angle_phi, ctau, ctauErr, B_VProb, B_Mass;
    Float_t B_Pt, B_MassFromSV, Mum_Pt, Mum_Eta, Mum_Phi, Mup_Pt, Mup_Eta, Mup_Phi;
    Float_t KmK_Pt, KmK_Eta, KmK_Phi, KpPi_Pt, KpPi_Eta, KpPi_Phi, N_PV; 
    Float_t  angle_cospsi_GEN, angle_costheta_GEN, angle_phi_GEN, ctau_GEN, Jpsi_Mass;
    Float_t Kp_Hits, Km_Hits;
    Float_t MC_Flavour, PhiKstar_Mass, MisTag;
    Int_t Tag;
    
    bs->SetBranchAddress("isBs",&isBs);
    bs->SetBranchAddress("HLT_JpsiMu",&HLT_JpsiMu);
    bs->SetBranchAddress("HLT_JpsiTkTk",&HLT_JpsiTkTk);
    bs->SetBranchAddress("angle_cospsi",&angle_cospsi);
    bs->SetBranchAddress("angle_costheta",& angle_costheta);
    bs->SetBranchAddress("angle_phi",&angle_phi);
    bs->SetBranchAddress("ctau",&ctau);
    bs->SetBranchAddress("ctauErr",&ctauErr);
    bs->SetBranchAddress("B_VProb",&B_VProb);
    bs->SetBranchAddress("B_Mass",&B_Mass);
    bs->SetBranchAddress("B_Pt",&B_Pt);
    bs->SetBranchAddress("B_MassFromSV",&B_MassFromSV);
    bs->SetBranchAddress("Mum_Pt",&Mum_Pt);
    bs->SetBranchAddress("Mum_Eta",&Mum_Eta);
    bs->SetBranchAddress("Mum_Phi",&Mum_Phi);
    bs->SetBranchAddress("Mup_Pt",&Mup_Pt);
    bs->SetBranchAddress("Mup_Eta",&Mup_Eta);
    bs->SetBranchAddress("Mup_Phi",&Mup_Phi);
    bs->SetBranchAddress("KmK_Pt",&KmK_Pt);
    bs->SetBranchAddress("KmK_Eta",&KmK_Eta);
    bs->SetBranchAddress("KmK_Phi",& KmK_Phi);
    bs->SetBranchAddress("KpPi_Pt",&KpPi_Pt);
    bs->SetBranchAddress("KpPi_Eta",&KpPi_Eta);
    bs->SetBranchAddress("KpPi_Phi",& KpPi_Phi);
    bs->SetBranchAddress("N_PV",&N_PV);
    bs->SetBranchAddress("angle_cospsi_GEN",&angle_cospsi_GEN);
    bs->SetBranchAddress("angle_costheta_GEN",&angle_costheta_GEN);
    bs->SetBranchAddress("angle_phi_GEN",&angle_phi_GEN);
    bs->SetBranchAddress("ctau_GEN",&ctau_GEN);
    bs->SetBranchAddress("Jpsi_Mass",&Jpsi_Mass);
    bs->SetBranchAddress("KpPi_Hits",&Kp_Hits);
    bs->SetBranchAddress("KmK_Hits",&Km_Hits);
    bs->SetBranchAddress("MC_Flavour",&MC_Flavour);
    bs->SetBranchAddress("PhiKstar_Mass",&PhiKstar_Mass);
    bs->SetBranchAddress("Tag",&Tag);
    bs->SetBranchAddress("MisTag",&MisTag);
    
       
    Double_t svmass,BsCt2DMC, BscosthetaMC, BscospsiMC, BsphiMC, BsCt2DMCErr,  BsCt2DMC_GEN, BscosthetaMC_GEN, BscospsiMC_GEN, BsphiMC_GEN, mistag;
    Int_t Bs_NPV, tag;   
       
    //TFile *f = new TFile("BsMC17WithCuts_JpsiTrkTrk_NoJpsiMu_Hits.root","recreate");
     TFile *f = new TFile("fittree_reco_recotag.root","recreate");
    TTree *fit = new TTree("fit","selected ntcltestle");
    fit->Branch("svmass",&svmass,"svmass/D");
    fit->Branch("BsCt2DMC",&BsCt2DMC,"BsCt2DMC/D");
    fit->Branch("BsCt2DMCErr",&BsCt2DMCErr,"BsCt2DMCErr/D");
    fit->Branch("BscosthetaMC",&BscosthetaMC,"BscosthetaMC/D");
    fit->Branch("BscospsiMC",&BscospsiMC,"BscospsiMC/D");
    fit->Branch("BsphiMC",&BsphiMC,"BsphiMC/D");
    fit->Branch("BsCt2DMC_GEN",&BsCt2DMC_GEN,"BsCt2DMC_GEN/D");
    fit->Branch("BscosthetaMC_GEN",&BscosthetaMC_GEN,"BscosthetaMC_GEN/D");
    fit->Branch("BscospsiMC_GEN",&BscospsiMC_GEN,"BscospsiMC_GEN/D");
    fit->Branch("BsphiMC_GEN",&BsphiMC_GEN,"BsphiMC_GEN/D");
    fit->Branch("tag",&tag,"tag/I");
    fit->Branch("mistag", &mistag, "mistag/D");
   /* fit->Branch("ct_gen",&ct_gen,"ct_gen/F");
    fit->Branch("costheta_gen",&costheta_gen,"costheta_gen/F");
    fit->Branch("cospsi_gen",&cospsi_gen,"cos_psi/F");
    fit->Branch("phi_gen",&phi_gen,"phi_gen/F");*/
    fit->Branch("Bs_NPV",&Bs_NPV,"Bs_NPV/I");
         svmass=-999999,  BsCt2DMC=-999999,  BsCt2DMCErr=-999999,  BscosthetaMC = -999999,  BscospsiMC=-999999,  BsphiMC=-999999,  BsCt2DMC_GEN=-999999,  BscosthetaMC_GEN=-999999,  BscospsiMC_GEN=-999999,  BsphiMC_GEN=-999999, mistag=-999999;
         tag=-999999, Bs_NPV=-999999; 
   for (Int_t i=0;i<n_entries;i++) {
    // for (Int_t i=0;i<50000;i++) {
          
        
                                               bs->GetEntry(i); 
                                              Int_t kpnpst = ( int(Kp_Hits) / 100 ) % 10000; 
                                              Int_t kpntrk = kpnpst % 100;
                                              Int_t kmnpst = ( int(Km_Hits) / 100 ) % 10000; 
                                              Int_t kmntrk = kmnpst % 100;
                                               if( isBs == 1 && B_MassFromSV>5.24 && B_MassFromSV<5.45 && Mum_Pt>3.5 && Mup_Pt>3.5 && B_Pt > 11 && KmK_Pt>1.2 && KpPi_Pt>1.2 && fabs(Mum_Eta) < 2.5 && fabs(Mup_Eta) < 2.5&& fabs(KmK_Eta) < 2.5 && fabs(KpPi_Eta) < 2.5&& ctau>0.007 && abs(Jpsi_Mass-3.0969)<0.150 && abs(PhiKstar_Mass-1.01946) < 0.010 && B_VProb>0.001 && HLT_JpsiMu==1  && kpntrk>=4 && kmntrk>=4 ){
 
 
                                                                                               if(B_MassFromSV != B_MassFromSV)continue;
                                                                                               if(ctau != ctau)continue;
                                                                                               if(ctauErr != ctauErr)continue;
                                                                                               if (MisTag !=MisTag)continue;
                                               
                                                                                              
                                                                                               svmass = B_MassFromSV;                                                                                               
                                                                                               BsCt2DMC          = ctau;
                                                                                               BsCt2DMCErr    = ctauErr;
                                                                                               BscosthetaMC = angle_costheta;
                                                                                               BscospsiMC    = angle_cospsi;
                                                                                               BsphiMC         = angle_phi;
                                                                                               Bs_NPV       =(int) N_PV;
                                                                                               BscosthetaMC_GEN = angle_costheta_GEN;
                                                                                               BscospsiMC_GEN    = angle_cospsi_GEN;
                                                                                               BsphiMC_GEN         = angle_phi_GEN;
                                                                                               BsCt2DMC_GEN           = ctau_GEN;
                                                                                               mistag                           = MisTag;
                                                                                               tag                                 =(int)Tag;
                                                                                               
                                                                                                
                                                                                                fit->Fill(); 
    
   }
                                                                                               
                                                                                               
                                                                                           /* costheta_gen = angle_costheta_GEN;
                                                                                               cospsi_gen    = angle_cospsi_GEN;
                                                                                               phi_gen         = angle_phi_GEN;
                                                                                               ct_gen           = ctau_GEN; */
                                                                                               
                                                                                               
                                                                                               
   
    }
    fit->Print();
    f->Write();
    }
