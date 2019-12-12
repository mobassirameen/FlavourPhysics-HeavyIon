

#include "TMinuit.h"

Float_t x1[10], x2[10], x3[10], z[20], errorz[20] ;

//______________________________________________________________________________
Double_t func(float x1, float x2, float x3,  Double_t *par)
{
 Double_t value=par[0]*x1 + par[1]*x2 + par[2]*x3+ par[3];
 return value;
}

//______________________________________________________________________________
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
   const Int_t nbins = 10;
   Int_t i;

//calculate chisquare
   Double_t chisq = 0;
   Double_t delta;
   for (i=0;i<nbins; i++) {
     delta  = (z[i]-func(x1[i], x2[i], x3[i],par))/errorz[i];
     chisq += delta*delta;
   }
   f = chisq;
}

//______________________________________________________________________________
void templatefit_tmnt()
{

// The z values
z[0]=278.123;
z[1]=23067.7;
z[2]=1139.6;
z[3]=581.871;
z[4]=26731.9;
z[5]=1771.06;
z[6]=966.376;
z[7]=2375.01;
z[8]=591.682;
z[9]=655.143;
// The errors on z values 5% of the tprofile values of the data points

   errorz[0]=z[0]*0.05;
   errorz[1]=z[1]*0.05;
   errorz[2]=z[2]*0.05;
   errorz[3]=z[3]*0.05;
   errorz[4]=z[4]*0.05;
   errorz[5]=z[5]*0.05;
   errorz[6]=z[6]*0.05;
   errorz[7]=z[7]*0.05;
   errorz[8]=z[8]*0.05;
   errorz[9]=z[9]*0.05;
// the x1 values
   x1[0]=0;
   x1[1]=1.0;
   x1[2]=0.1;
   x1[3]=0.05;
   x1[4]=0;
   x1[5]=0;
   x1[6]=0;
   x1[7]=0;
   x1[8]=0;
   x1[9]=0;
// the x2 values
   x2[0]=0;
   x2[1]=0;
   x2[2]=0;
   x2[3]=0;
   x2[4]=1.0;
   x2[5]=0.1;
   x2[6]=0.05;
   x2[7]=0;
   x2[8]=0;
   x2[9]=0;
// the x3 values
   x3[0]=0;
   x3[1]=0;
   x3[2]=0;
   x3[3]=0;
   x3[4]=0;
   x3[5]=0;
   x3[6]=0;
   x3[7]=1.0;
   x3[8]=0.1;
   x3[9]=0.05;
   
   /*for(int l = 0; l<=10; l++)
   {
   std::cout<<*func(x1[l], x2[l], x3[l], par)<<"\n";
   } */
   TMinuit *gMinuit = new TMinuit(4);  //initialize TMinuit with a maximum of 5 params
   gMinuit->SetFCN(fcn);

  
   Double_t arglist[10];
   Int_t ierflg = 0;

   arglist[0] = 1;
   gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

// Set starting values and step sizes for parameters
   static Double_t vstart[4] = {3, 1 , 0.1 , 0.01};
   static Double_t step[4] = {0.1 , 0.1 , 0.01 , 0.001};
   gMinuit->mnparm(0, "a1", vstart[0], step[0], 0,0,ierflg);
   gMinuit->mnparm(1, "a2", vstart[1], step[1], 0,0,ierflg);
   gMinuit->mnparm(2, "a3", vstart[2], step[2], 0,0,ierflg);
   gMinuit->mnparm(3, "a4", vstart[3], step[3], 0,0,ierflg);

// Now ready for minimization step
   arglist[0] = 500;
   arglist[1] = 1.;
   gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

// Print results
   Double_t amin,edm,errdef;
   Int_t nvpar,nparx,icstat;
   gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
   //gMinuit->mnprin(3,amin);
   TCanvas *c1 = new TCanvas("c1","Parameter p0",200,10,700,700);
   c1->cd();
   gMinuit->Command("SCAn 0");
   TGraph *gr1 = (TGraph*)gMinuit->GetPlot();
   gr1->SetMarkerStyle(21);
   gr1->Draw("alp");
   TCanvas *c2 = new TCanvas("c2","parameter p1",200,10,700,700);
   c2->cd();
   gMinuit->Command("SCAn 1");
   TGraph *gr2 = (TGraph*)gMinuit->GetPlot();
   gr2->SetMarkerStyle(21);
   gr2->Draw("alp");
   TCanvas *c3 = new TCanvas("c3","parameter p2",200,10,700,700);
   c3->cd();
   gMinuit->Command("SCAn 2");
   TGraph *gr3 = (TGraph*)gMinuit->GetPlot();
   gr3->SetMarkerStyle(21);
   gr3->Draw("alp");
   TCanvas *c4 = new TCanvas("c4","parameter p3",200,10,700,700);
   c4->cd();
   gMinuit->Command("SCAn 3");
   TGraph *gr4 = (TGraph*)gMinuit->GetPlot();
   gr4->SetMarkerStyle(21);
   gr4->Draw("alp");
}

