#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TF1.h>
#include <TH1D.h>
#include "TGraphAsymmErrors.h"
#include "TVirtualFFT.h"
#include "TBinomialEfficiencyFitter.h"
#include "TVectorF.h"
#include "TPaveText.h"
#include <vector>
#include "TTree.h"
#include "TFile.h"
#include "TEventList.h"
#include "Riostream.h"
#include "string.h"
#include "TList.h"
#include "TDirectory.h"
#include "TCut.h"
#include "TChain.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TEfficiency.h"
#include "TLorentzVector.h"
#include "TProfile.h"

using namespace std;

void zdcdigi()
{

        
  //TCanvas *c = new TCanvas("c", "c", 0, 0, 1000, 600);
  TFile *fIn1 = new TFile("/Users/ab/Documents/ZDCCalibration/zdc_digi_minbias.root");
  //if (!fIn1){return;}
  TTree* mumu = (TTree*)fIn1->Get("analyzer/zdcdigi");
  Int_t nentries = mumu->GetEntries();
  std::cout<<"No of entries"<<nentries<<"\n";

  const int NSIDE=2; const char* stit[NSIDE] = {"#minus","#plus"};  const char* stit2[NSIDE] = {"neg","pos"};
  const int NTYPE=2; const char* ttit[NTYPE] = {"EM","HAD"};
  const int NCH=5; const char* ctit[NTYPE][NCH] = {
    {"1","2","3","4","5"}, //HD sections run only 1-4                                          
    {"1","2","3","4","5"} //EM sections run 1-5                                                
  };

  TProfile* chargetime_ZDC[2][2][5];
  TProfile* chargetime_ZDC_Tot[2][2][5];
  TProfile* chargetime_ZDC_PU[2][2][5];
  TH2F* adc_ZDC[2][2][5];
  TH2F* adc_RPD[2][16];

  TH2F* fC_ZDC[2][2][5];
  TH2F* fC_RPD[2][16];
  TH1D* EFC_ZDC[10];
   TH1D* PU[10];
  
 /* TH1D *pull[50];
 TH1D *resolution[50];
       for (int l = 0; l<5; l++){
                                                       resolution[l] = new TH1D(Form("Resolution%dth_Bin", l),Form("Resolution %dth_Bin; ct_{GEN}-ct_{RECO}; Events",l),50,-0.004,0.004);
                                                       pull[l] = new TH1D(Form("Pull%dth_Bin", l),Form("Pull%dth_Bin; ct_{GEN}-ct_{RECO}/#Delta ct; Events",l),  50, -5., 5.);
                                               }
 */ 
  for(int l =0 ; l < 10; l++){
  EFC_ZDC[l] = new TH1D(Form("Eventbyevent_Charge_for_%d_slice",l), Form("ZDCChargeValues_%d_slice",l), 10, 0, 250);
 // PU[l] = new TH1D(Form("Eventbyevent_Charge_for_%d_slice",l), Form("ZDCChargeValues_%d_slice",l), 10, 0, 10);
  } 

  for(int i = 0; i < 2; i++){
    for(int j = 0; j < 2; j++)
      for(int k = 0; k < 5; k++){
     adc_ZDC[i][j][k] = new TH2F(Form("adc ZDC %s%s channel %d",ttit[j],stit[i],k+1),Form("ZDC %s%s channel %d;TS [25 ns];ADC",ttit[j],stit\
											     [i],k+1),10,-0.5,9.5,256,0,256);
     fC_ZDC[i][j][k] = new TH2F(Form("fC ZDC %s%s channel %d",ttit[j],stit[i],k+1),Form("ZDC %s%s channel %d;TS [25 ns];Q [fC]",
ttit[j],stit[i],k+1),10,-0.5,9.5,2500,1,350000);
     //   chargetime_ZDC[i][j][k] = new TProfile(Form("fC ZDC %s%s channel %d",ttit[j],stit[i],k+1),Form("ZDC %s%s channel %d;TS [25 ns];Q [fC]",
     //											    ttit[j],stit[i],k+1),10,-0.5,9.5,50,1,350000);

     chargetime_ZDC[i][j][k] = new TProfile(Form("fCZDC%s%schannel%d",ttit[j],stit[i],k+1),Form("ZDC %s%s channel %d;TS [25 ns];Q [fC]",ttit[j],stit[i],k+1),10,-0.5,9.5,1,350000);
     chargetime_ZDC_Tot[i][j][k] = new TProfile(Form("fCZDC%s%schannel_SB%d",ttit[j],stit[i],k+1),Form("ZDC %s%s channel %d;TS [25 ns];Q [fC]",ttit[j],stit[i],k+1),10,-0.5,9.5,1,350000);
     chargetime_ZDC_PU[i][j][k] = new TProfile(Form("PU%s%schannel_prepost%d",ttit[j],stit[i],k+1),Form("ZDC %s%s channel %d;TS [25 ns];Q [fC]",ttit[j],stit[i],k+1),10,-0.5,9.5,1,350000);


      }

    


 for(int j = 0; j < 16; j++){
      adc_RPD[i][j] = new TH2F(Form("adc RPD%s channel %d",stit[i],j+1),Form("RPD%s channel %d;TS [25 ns];ADC",stit[i],j+1),10,-0.5,9.5,256,0,\
			       256);
      fC_RPD[i][j] = new TH2F(Form("fC RPD%s channel %d",stit[i],j+1),Form("RPD%s channel %d;TS [25 ns];Q [fC]",stit[i],j+1),10,-0.5,9.5,256,1\
			      ,350000);
    }

  }


TFile *fOut = TFile::Open("signalzdc.root","RECREATE");
TTree *fit = new TTree("fit","selected ntcltestle");
   if (!fOut) { return; }


  const int NTS=10;
  TLeaf* bxidLeaf = (TLeaf*) mumu->GetLeaf("bxid");
  TLeaf* zsideLeaf = (TLeaf*) mumu->GetLeaf("zside");
  TLeaf* sectionLeaf = (TLeaf*) mumu->GetLeaf("section");
  TLeaf* channelLeaf = (TLeaf*) mumu->GetLeaf("channel");
  TLeaf* ntrk = (TLeaf*) mumu->GetLeaf("nTrack");
  TLeaf* nHFneg = (TLeaf*) mumu->GetLeaf("nHFneg");
  TLeaf* nHFpos = (TLeaf*) mumu->GetLeaf("nHFpos");
  TLeaf* adcLeaf[NTS];
  TLeaf* fCleaf[NTS];

  for(int iTS = 0; iTS < NTS; iTS++){
    adcLeaf[iTS] = (TLeaf*) mumu->GetLeaf(Form("adc%d",iTS));
    fCleaf[iTS] = (TLeaf*) mumu->GetLeaf(Form("nfC%d",iTS));
  }
TRandom *r1 = new TRandom();
 std::cout<<"Event"<<"    "<<"Timeslice"<<"     "<<"Bin-values in each TS"<<"\n";
  //TH1D* hadcseven = new TH1D("hadcseven","adc_7 signal shape ; signal shape of 7 ;Events", 100 ,0, 4);
  //TH1D* hnfcseven = new TH1D("hnfcseven","nfc_7 signal shape ; signal shape of 7 ;Events", 100 ,0.0, 400.0);   
   for(int i = 0; i < mumu->GetEntries(); i++){
    //for(int i = 3452; i < 3502; i++){
    mumu->GetEntry(i);
     
      for(int n = 0; n < 50; n++){

	int side = (int)((zsideLeaf->GetValue(n)+1)/2.0);
	int type = (int)(sectionLeaf->GetValue(n))-1;
	int channel = (int)(channelLeaf->GetValue(n))-1;

                         /*
                         For each event
check if  ZDC+ Had2 has any signal in TS1 or TS7. If it does reject the event, then
check if  ZDC+ Had3 has any signal in TS1 or TS7. If it does reject the event, then
check if  ZDC+ Had4 has any signal in TS1 or TS7. If it does reject the event
                         */                         
	for(int iTS = 0; iTS < 10; iTS++){
	//if(iTS !=4){
	  if(type != 3){ // EM or HAD section                                                                                                    
	                   chargetime_ZDC_Tot[side][type][channel]->Fill(iTS,fCleaf[iTS]->GetValue(n),1);
	             if (side == 0 && type ==0)continue;
	                if (side ==1 && type ==0)continue;
	                if(side ==0 && type == 1)continue;
	                 if (side ==1 && type==1   && channel == 1 &&( (fCleaf[1]->GetValue(n))/(fCleaf[4]->GetValue(n)) >0.01)  && ((fCleaf[7]->GetValue(n))/(fCleaf[4]->GetValue(n)) >0.05))continue;
	                  if (side ==1 && type==1   && channel == 2 && (fCleaf[1]->GetValue(n))/(fCleaf[4]->GetValue(n)) >0.01  && (fCleaf[7]->GetValue(n))/(fCleaf[4]->GetValue(n)) >0.05)continue;
	                  if (side ==1 && type==1   && channel == 3 && (fCleaf[1]->GetValue(n))/(fCleaf[4]->GetValue(n)) >0.01  && (fCleaf[7]->GetValue(n))/(fCleaf[4]->GetValue(n)) >0.05)continue;
	                 if((fCleaf[1]->GetValue(n))/(fCleaf[4]->GetValue(n)) >0.01) continue;
	                 if((fCleaf[7]->GetValue(n))/(fCleaf[4]->GetValue(n)) >0.05) continue;
	                 adc_ZDC[side][type][channel]->Fill(iTS,adcLeaf[iTS]->GetValue(n));
	                 fC_ZDC[side][type][channel]->Fill(iTS,fCleaf[iTS]->GetValue(n));	     
			 chargetime_ZDC[side][type][channel]->Fill(iTS,fCleaf[iTS]->GetValue(n),1);
			 chargetime_ZDC_PU[side][type][channel]->Fill(iTS,(chargetime_ZDC_Tot[side][type][channel]->GetBinContent(iTS+1)- chargetime_ZDC[side][type][channel]->GetBinContent(iTS+1)),1);
			// chargetime_ZDC[side][type][channel]->Clone("EFC_ZDC[iTS]");           
			 Double_t mybin = chargetime_ZDC_PU[side][type][channel]->GetBinContent(iTS+1);
			
			//std::cout<<i<<"\t"<<side<<"\t"<<type<<"\t"<<channel<<"\t"<<iTS<<"      "<<mybin<<"\n";  
			//std::cout<<"Event"<<i<<"TimeSlice"<<iTS<<"Bin Value"<<mybin<<"\n";  
			//if(side==1 && type==1){
			//std::cout<<"z["<<iTS<<"]="<<mybin<<";"<<"\n";
			//}
			// std::cout<<"Charge"<<fCleaf[4]->GetValue(n)<<"\n";
			//if(side == 1&& type ==1){
			//std::cout<<i<<"\t"<<side<<"\t"<<type<<"\t"<<channel<<"\t"<<"z["<<iTS<<"]="<<mybin<<";"<<"\n";}
			//std::cout<<i<<"\t"<<channel<<"\t"<<"z["<<iTS<<"]="<<mybin<<";"<<"\n";}
                       //event->SetBinContent(iTS+1,fCleaf[iTS]->GetValue(n));                                                         
	                }//}
	             //   if(iTS !=4){
	  //if(type != 3){ // EM or HAD section                                                                                                    
	               			// chargetime_ZDC_SB[side][type][channel]->Fill(iTS,fCleaf[iTS]->GetValue(n),1);	                 
                                         //event->SetBinContent(iTS+1,fCleaf[iTS]->GetValue(n));                                                         
	      //          }}
	    else{ // RPD section                                                                                                               
           	    //if(adcLeaf[iTS]->GetValue(n) > max)                                                                    
	            //  max = adcLeaf[iTS]->GetValue(n);
                    adc_RPD[side][channel]->Fill(iTS,adcLeaf[iTS]->GetValue(n));
	            fC_RPD[side][channel]->Fill(iTS,fCleaf[iTS]->GetValue(n));
	        }
        	}

                }

           //if(i % 100000 == 0) 
       //  std::cout << i << " events are processed." << std::endl;
  }
  TCanvas *c1 = new TCanvas();      
  //c1->SetLogy();
  c1->SetLogz();

  for(int iside = 0; iside < 2; iside++)
    for(int itype = 0; itype < 2; itype++)
      for(int ich = 0; ich < 4; ich++){
        adc_ZDC[iside][itype][ich]->Draw("colz");
        //c1->SaveAs(Form("/home/cms/ZDCCalibration/pic/adc_ZDC_%s_%s_channel_%d.png",ttit[itype],stit2[iside],ich+1));
        //adc_ZDC[iside][itype][ich]->Write();
	// f2.Write();
      }


  for(int iside = 0; iside < 2; iside++)
    for(int itype = 0; itype < 2; itype++)
      for(int ich = 0; ich < 4; ich++){
        fC_ZDC[iside][itype][ich]->Draw("colz");
       //c1->SaveAs(Form("/home/cms/ZDCCalibration/image/fC_ZDC_%s_%s_channel_%d.png",ttit[itype],stit2[iside],ich+1));
        //fC_ZDC[iside][itype][ich]->Write();
	//  f2.Write();
      }
  for(int iside = 0; iside < 2; iside++)
    for(int itype = 0; itype < 2; itype++)
      for(int ich = 0; ich < 1; ich++){
        chargetime_ZDC[iside][itype][ich]->Draw("colz");
         chargetime_ZDC_PU[iside][itype][ich]->Draw("colz");
        //chargetime_ZDC_SB[iside][itype][ich]->Draw("ALP");
       // chargetime_ZDC_SB[iside][itype][ich]->Write();
        //chargetime_ZDC[iside][itype][ich]->Write();
        //c1->SaveAs(Form("/home/cms/ZDCCalibration/image/chargetime_ZDC_%s_%s_channel_%d.png",ttit[itype],stit2[iside],ich+1));
        c1->SaveAs(Form("/Users/ab/Documents/ZDCCalibration/newimg/chargetime_ZDC_%s_%s_channel_%d.png",ttit[itype],stit2[iside],ich+1));
            fit->Fill();                                                                                                                                                                          
      }





  // RPD- channels separately                                                                                                                  
  /*for(int i = 0; i < 16; i++){
    adc_RPD[0][i]->Draw("colz");
    c1->SaveAs(Form("/Users/md/Documents/ZDCCalibration/pic/adc_RPD_Minus_channel_%d.png",i));
    //adc_RPD[0][i]->Write();
    //f2.Write();
  }
auto cpu = new TCanvas();
cpu->cd();
for(int g=0;g<10;g++){
PU[g]->Draw("colz");
}
  // RPD+ channels separately                                                                                                                  
  for(int i = 0; i < 16; i++){
    adc_RPD[1][i]->Draw("colz");
    c1->SaveAs(Form("/Users/md/Documents/ZDCCalibration/pic/adc_RPD_Plus_channel_%d.png",i));
    //adc_RPD[1][i]->Write();
    //f2.Write();
  }

  // RPD- channels separately                                                                                                                  
  for(int i = 0; i < 16; i++){
    fC_RPD[0][i]->Draw("colz");
    c1->SaveAs(Form("/Users/md/Documents/ZDCCalibration/pic/fC_RPD_Minus_channel_%d.png",i));
    //fC_RPD[0][i]->Write();
    //f2.Write();
  }

  // RPD+ channels separately                                                                                                                  
  for(int i = 0; i < 16; i++){
    fC_RPD[1][i]->Draw("hist e");
    c1->SaveAs(Form("/Users/md/Documents/ZDCCalibration/pic/fC_RPD_Plus_channel_%d.png",i));
    //fC_RPD[1][i]->Write();
    // f2.Write();
  }*/

fit->Write();
   delete fit;

}
