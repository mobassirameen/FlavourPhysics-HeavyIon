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

void resolution_study(){
   // TFile *fIn1 = new TFile("/home/cms/BMeson/FileRepository/Generalplots/BsJpsiPhi_dG0.root");
    TFile *fIn1 = new TFile("/eos/cms/store/user/fiori/Bs_Ntuples_Run2/Tagged_Ntuples/Tot_2017MC_withTag.root");
    if (!fIn1){return;}
    TTree* bs = (TTree*)fIn1->Get("OutTree");
    Int_t n_entries = bs->GetEntries();
    std::cout<<n_entries<<"\n";
    Float_t  isBs, HLT_JpsiMu, angle_cospsi, angle_costheta, angle_phi, ctau, ctauErr, B_VProb, B_Mass;
    Float_t B_Pt, B_MassFromSV, Mum_Pt, Mum_Eta, Mum_Phi, Mup_Pt, Mup_Eta, Mup_Phi;
    Float_t KmK_Pt, KmK_Eta, KmK_Phi, KpPi_Pt, KpPi_Eta, KpPi_Phi, N_PV; 
    Float_t  angle_cospsi_GEN, angle_costheta_GEN, angle_phi_GEN, ctau_GEN, Jpsi_Mass;
    
    bs->SetBranchAddress("isBs",&isBs);
    bs->SetBranchAddress("HLT_JpsiMu",&HLT_JpsiMu);
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
    TH1D *resolution = new TH1D("resolution", "Resolution; ct_{GEN}- ct_{RECO} (cm); Events ", 100, -0.01, 0.01);
    TH1D *cterror = new TH1D("cterror", "CtError; CtErrorEvents ", 100, 0.00, 0.005); 
    TH1D *Recoct = new TH1D("Recoct", "ct_{reco};ct_{reco} (cm); Events",100, 0.007, .40);
    TH1D *Gencutct = new TH1D("Gencutct","ct_{gencut};ct_{gencut};Events",100, 0.007, .40);
    TH1D *Genct = new TH1D("Genct", "ct_{gen};ct_{gen};Events", 100, 0.0, .40);
    TH1D *ctpull = new TH1D("ctpull", "ctpull;Pull_{ct};Events", 100, -4.0, 4.0);
    
    Float_t svmass, ct, ctErr, costheta, cospsi, phi, npv;   
       
     for (Int_t i=0;i<n_entries;i++) {
      //for (Int_t i=0;i<50000;i++) {
            svmass=-999999,  ct=-999999,  ctErr=-999999,  costheta = -999999,  cospsi=-999999,  phi=-999999,   npv=-999999; 
                                               bs->GetEntry(i); 
                                                                                               
                                             //if(isBs == 1){                                           
                                             /*if (B_MassFromSV < 5.65 && B_MassFromSV > 5.20) {
                                             if (fabs(Mum_Eta) < 2.5 && fabs(Mup_Eta) < 2.5&& fabs(KmK_Eta) < 2.5 && fabs(KpPi_Eta) < 2.5&& KmK_Pt > 1.2)  && KpPi_Pt >1.2 ){         
                                             if( B_VProb >  0.001){
                                             if (B_Pt > 11.0){
                                             if (Mup_Pt > 3.5 && Mum_Pt > 3.5){
                                             if (HLT_JpsiMu == 1){
                                             if(N_PV < 30){*/
 if(isBs==1 && B_MassFromSV>5.24 && B_MassFromSV<5.45 && Mum_Pt>3.5 && Mup_Pt>3.5 && B_Pt > 11 && KmK_Pt>1.2 && KpPi_Pt>1.2 && fabs(Mum_Eta) < 2.5 && fabs(Mup_Eta) < 2.5&& fabs(KmK_Eta) < 2.5 && fabs(KpPi_Eta) < 2.5&& ctau>0.007 && abs(Jpsi_Mass-3.0969)<0.150 && B_VProb>0.01 && HLT_JpsiMu==1 ){
                                                                                           
                                                                                               Recoct->Fill(ctau); Gencutct->Fill(ctau_GEN);cterror->Fill(ctauErr);

    
    
                                                                                           }
                                                                                               
                                                                                            Genct->Fill(ctau_GEN);   
                                                                                                
                                                                                               
                                                                                               
   
    }
   
TCanvas *cc = new TCanvas("cc", "cc",0,0,800,600);
cc->Divide(2,4);
cc->cd(1);
Recoct->Draw();
cc->cd(2);
Gencutct->Draw();
cc->cd(3);
cterror->Draw();
cc->cd(4);
Genct->Draw();

for(int j=0; j<=100; j++){
//resolution->Fill((Genct->GetBinContent(j))-(Recoct->GetBinContent(j)));
ctpull->Fill(((Genct->GetBinContent(j))-(Recoct->GetBinContent(j)))/(cterror->GetBinContent(j)));
}
cc->cd(5);
resolution->Draw();
cc->cd(6);
ctpull->Draw();


    }
