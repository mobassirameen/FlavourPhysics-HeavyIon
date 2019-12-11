// vim: set sw=4 sts=4 filetype=cpp fdm=marker et: 

// -----------------------------------------------
//       Author: Xin Shi <Xin.Shi@cern.ch> 
//       Created:   [2013-08-15 Thu 14:54] 
// -----------------------------------------------
#include <stdio.h>
#include <sstream>
#include <sys/stat.h>
#include <math.h>
#include <string.h>
//#include <regex> //c++11 feature should be fine using gcc491.

#include <TSystem.h>
#include <TStyle.h>
#include <TF1.h>
#include <TF2.h>
#include <TH1.h>
#include <TH2.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>
#include <TMinuit.h>
#include <TRandom3.h>
#include <TFile.h>
#include <TPad.h> 
#include <TCanvas.h> 
#include <TChain.h> 
#include <TChainElement.h>
#include <TPaveText.h>
#include <TLatex.h>
#include <TString.h>
#include <TGraphAsymmErrors.h>
#include <TLorentzVector.h>

#include <RooConstVar.h>
#include <RooRealVar.h>
#include <RooAbsPdf.h>
#include <RooAddPdf.h>
#include <RooGaussian.h>
#include <RooBifurGauss.h>
#include <RooChebychev.h> 
#include <RooGenericPdf.h> 
#include <RooExponential.h>
#include <RooBernstein.h>
#include <RooPolynomial.h>
#include <RooExtendPdf.h>
#include <RooProdPdf.h>
#include <RooDataHist.h>
#include <RooDataSet.h>
#include <RooAbsData.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooAddition.h>
#include <RooProduct.h>

#include "tools.h" 

using namespace std; 
using namespace RooFit;

// Tags configration
bool is7TeVCheck = false; // Using 2011 efficiency map.
int isCDFcut = 0; // 0 for off, 1 for on;
TChain *ch=new TChain("tree");
//Constants, Fit results for efficiency, etc.. //{{{
char genQ2range[10][32] = {"genQ2 < 2.00 && genQ2 > 1.00",
                           "genQ2 < 4.30 && genQ2 > 2.00",
                           "genQ2 < 8.68 && genQ2 > 4.30",
                           "genQ2 <10.09 && genQ2 > 8.68",
                           "genQ2 <12.86 && genQ2 >10.09",
                           "genQ2 <14.18 && genQ2 >12.86",
                           "genQ2 <16.00 && genQ2 >14.18",
                           "genQ2 <19.00 && genQ2 >16.00",
                           "genQ2 <19.00 && genQ2 > 1.00",
                           "genQ2 < 6.00 && genQ2 > 1.00"};
char q2range[10][32] = {"Q2 < 2.00 && Q2 > 1.00",
                        "Q2 < 4.30 && Q2 > 2.00",
                        "Q2 < 8.68 && Q2 > 4.30",
                        "Q2 <10.09 && Q2 > 8.68",
                        "Q2 <12.86 && Q2 >10.09",
                        "Q2 <14.18 && Q2 >12.86",
                        "Q2 <16.00 && Q2 >14.18",
                        "Q2 <19.00 && Q2 >16.00",
                        "Q2 <19.00 && Q2 > 1.00",
                        "Q2 < 6.00 && Q2 > 1.00"};
double q2rangedn[10] = {1.00 , 2.00 , 4.30 , 8.68  , 10.09 , 12.86 , 14.18 , 16.00 ,  1.00 , 1.00};
double q2rangeup[10] = {2.00 , 4.30 , 8.68 , 10.09 , 12.86 , 14.18 , 16.00 , 19.00 , 19.00 , 6.00};
char mumuMassWindow[7][512] = { "Mumumass > 0",
                                "(Mumumass > 3.096916+3.5*Mumumasserr || Mumumass < 3.096916-5.5*Mumumasserr) && (Mumumass > 3.686109+3.5*Mumumasserr || Mumumass < 3.686109-3.5*Mumumasserr)",
                                "(Mumumass < 3.096916+3*Mumumasserr && Mumumass > 3.096916-5*Mumumasserr) || (Mumumass < 3.686109+3*Mumumasserr && Mumumass > 3.686109-3*Mumumasserr)",
                                "Mumumass < 3.096916+3*Mumumasserr && Mumumass > 3.096916-5*Mumumasserr",
                                "Mumumass < 3.686109+3*Mumumasserr && Mumumass > 3.686109-3*Mumumasserr",
                                "(Mumumass*Mumumass<8.68 && Bmass-Mumumass>2.182+0.16 && Bmass-Mumumass<2.182-0.16) || (Mumumass*Mumumass>10.09 && Mumumass*Mumumass<12.86 && Bmass-Mumumass>1.593+0.06 && Bmass-Mumumass<1.593-0.06) || (Mumumass*Mumumass>14.18 && Bmass-Mumumass>1.593+0.06 && Bmass-Mumumass<1.593-0.06)",
                                "(Mumumass*Mumumass < 8.68 && ( Bmass-Mumumass < 2.182+0.16 || Bmass-Mumumass > 2.182-0.16)) || (Mumumass*Mumumass>8.68 && Mumumass*Mumumass<10.09) || (Mumumass*Mumumass > 10.09 && Mumumass < 12.86 && ( Bmass-Mumumass < 1.593+0.06 || Bmass-Mumumass > 1.593-0.06)) || (Mumumass*Mumumass>12.86 && Mumumass*Mumumass<14.08) || (Mumumass*Mumumass > 14.08 && (Bmass-Mumumass > 1.593-0.06 || Bmass-Mumumass < 1.593+0.06))"
};//None, sig, bkg, Jpsi, Psi2S, CDF, anti-CDF
double genAfb[8]={-0.160,-0.066,0.182,0.317,0.374,0.412,0.421,0.376};
double genAfberr[8]={0.000434,0.000285,0.000223,0.000383,0.000277,0.000421,0.000395,0.000422};
double genFl [8]={0.705,0.791,0.649,0.524,0.454,0.399,0.369,0.341};
double genFlerr[8]={0.000568,0.000409,0.000271,0.000420,0.000287,0.000415,0.000369,0.000361};
double arrRecPar2011   [8][20] = {
            {0.0017687,-8.73499e-06,-0.000935262,-0.000787674,
            -0.00529644,0.00155626,0.00356542,0.000171319,
            0,0,0,0,
            0.00551706,-0.00339678,-0.00463866,0.00240596,
            -0.00197731,0.00184597,0.00200529,-0.00178549},
            {0.00157318,7.37446e-05,-0.000749493,-0.000781159,
            -0.00288927,0.000261854,0.00187545,0.000883533,
            0,0,0,0,
            0.00132593,-0.000331488,-0.00112174,-0.000108788,
            0,0,0,0},
            {0.00118453,-1.23943e-05,-0.000454793,-0.000541742,
            -2.12267e-05,0.000714674,-9.69142e-05,-0.000768923,
            -7.59472e-05,0.000229203,5.36592e-05,-6.0488e-05,
            -0.00106615,-0.000896725,0.000710503,0.00156238,
            0,0,0,0},
            {0.00106852,-7.35469e-05,-0.000403409,-0.00035673,
            0.000637779,0.000294361,-0.000303145,-0.000885132,
            5.60979e-05,8.91354e-05,-0.000103992,-5.73688e-05,
            -0.00112888,9.3696e-05,0.000728873,0.000650714,
            0,0,0,0},
            {0.00102525,-0.000118438,-0.000346987,-0.000315243,
            0.000741172,0.000462559,-0.000313513,-0.00064603,
            -0.0007614,0.00153595,0.000840568,-0.0014875,
            -1.66019e-05,-0.00190496,-0.000546402,0.00220778,
            0,0,0,0},
            {0.00108772,0.000130305,-0.000457955,-0.000514284,
            0.000225576,-0.00108671,0.000283481,0.000592928,
            4.66766e-05,-0.00017112,4.27933e-05,0.000247925,
            7.59011e-05,0.000852078,-0.000514188,-8.52169e-05,
            0,0,0,0},
            {0.000964546,4.77019e-05,-0.000249802,-0.00035122,
            0.00120212,-0.0013481,-0.00081191,0.00109161,
            -0.000838523,0.00118858,0.000579591,-0.00101538,
            0,0,0,0,
            0,0,0,0},
            {0.00103644,7.27679e-05,-0.000366081,-0.000162294,
            0.000571208,-0.00089218,-1.99667e-05,0.00066918,
            -0.000423471,0.000360459,0.000228834,-0.000207738,
            0,0,0,0,
            0,0,0,0}};

double arrRecParErr2011[8][20] = {
            {9.77352e-06,1.45679e-05,1.44267e-05,1.71613e-05,
            1.7851e-05,2.42707e-05,2.66439e-05,2.89112e-05,
            0,0,0,0,
            2.06968e-05,2.87918e-05,3.2287e-05,3.52259e-05,
            1.71744e-05,2.50134e-05,2.75925e-05,3.13877e-05},
            {1.93048e-05,4.54403e-05,3.19836e-05,6.34448e-05,
            6.97046e-05,0.000175846,0.000127026,0.000252818,
            0,0,0,0,
            5.60232e-05,0.000145372,0.000105324,0.000210834,
            0,0,0,0},
            {1.20764e-05,2.99275e-05,2.12433e-05,4.24699e-05,
            5.88911e-05,0.000155115,0.0001176,0.000229393,
            1.73594e-05,3.28355e-05,3.60904e-05,5.2698e-05,
            5.76856e-05,0.000147459,0.000126591,0.000229798,
            0,0,0,0},
            {1.07973e-05,2.65747e-05,1.89528e-05,3.78338e-05,
            8.26138e-05,0.000208233,0.000148573,0.000297846,
            1.46533e-05,4.07354e-05,3.25714e-05,6.24532e-05,
            0.000102431,0.000261492,0.00018865,0.000376784,
            0,0,0,0},
            {1.43812e-05,3.69822e-05,2.68571e-05,5.34828e-05,
            0.000109918,0.000289434,0.000212286,0.000422439,
            4.35946e-05,9.76997e-05,7.60719e-05,0.000139613,
            0.000142318,0.000374508,0.000278563,0.000548987,
            0,0,0,0},
            {4.57551e-05,0.000116528,7.84802e-05,0.000164459,
            0.000350888,0.000882504,0.000628837,0.00126522,
            6.95431e-05,0.000192566,0.000150432,0.000293712,
            0.000438904,0.00111066,0.000809841,0.00160928,
            0,0,0,0},
            {1.41088e-05,3.73896e-05,2.85699e-05,5.60619e-05,
            5.72502e-05,0.000142174,0.000109448,0.000209836,
            5.66788e-05,0.000139079,0.000107494,0.00020466,
            0,0,0,0,
            0,0,0,0},
            {1.32107e-05,3.6017e-05,2.71503e-05,5.4467e-05,
            3.99561e-05,0.000108688,8.82189e-05,0.000169816,
            3.68927e-05,0.000100146,8.35791e-05,0.000158447,
            0,0,0,0,
            0,0,0,0}};

std::string f_accXrecoEff_ord0[8] = { // default values
    "11627.982364*((1.195422e-04*exp(-0.5*((CosThetaL-(-1.727343e-01))/2.021796e-01)**2)+1.156964e-04*exp(-0.5*((CosThetaL-(2.507083e-01))/2.478225e-01)**2)+4.629809e-05*exp(-0.5*((CosThetaL-(-5.148565e-01))/1.407258e-01)**2))*(7.165504e-05-2.621913e-05*CosThetaK+1.453609e-04*CosThetaK**2+2.274953e-05*CosThetaK**3-2.398253e-04*CosThetaK**4-4.428545e-05*CosThetaK**5+9.677067e-05*CosThetaK**6))",
    "9916.540629*((7.278431e-05*exp(-0.5*((CosThetaL-(-4.905860e-01))/1.878949e-01)**2)+7.448700e-05*exp(-0.5*((CosThetaL-(5.058518e-01))/2.003984e-01)**2)+1.425194e-04*exp(-0.5*((CosThetaL-(1.313125e-02))/2.957232e-01)**2))*(8.311598e-05-2.316101e-05*CosThetaK+1.476586e-04*CosThetaK**2-2.367362e-05*CosThetaK**3-1.683845e-04*CosThetaK**4+5.042865e-06*CosThetaK**5+1.243843e-05*CosThetaK**6))",
    "8.353802e+03*((1.344021e-04+1.980409e-05*CosThetaL+4.029664e-05*CosThetaL**2+1.560540e-05*CosThetaL**3-2.131400e-04*CosThetaL**4-3.310795e-05*CosThetaL**5+5.462426e-05*CosThetaL**6)*(1.136958e-04-3.718097e-05*CosThetaK+6.443598e-05*CosThetaK**2+5.683602e-05*CosThetaK**3-4.802073e-05*CosThetaK**4-7.413557e-05*CosThetaK**5-3.107261e-05*CosThetaK**6))",
    "4.392777e+04*((2.114489e-05+2.400662e-06*CosThetaL+2.759247e-05*CosThetaL**2+1.100568e-06*CosThetaL**3-4.538219e-05*CosThetaL**4-2.412249e-06*CosThetaL**5+5.307765e-06*CosThetaL**6)*(2.406814e-05-7.583489e-06*CosThetaK-9.968329e-06*CosThetaK**2+1.463576e-05*CosThetaK**3+3.247851e-05*CosThetaK**4-1.619795e-05*CosThetaK**5-2.949584e-05*CosThetaK**6))",
    "6.506619e+03*((1.349742e-04+1.528919e-05*CosThetaL+8.605597e-05*CosThetaL**2+1.312572e-05*CosThetaL**3-2.948919e-05*CosThetaL**4-9.566140e-06*CosThetaL**5-5.879247e-05*CosThetaL**6)*(1.581494e-04-3.384666e-05*CosThetaK-1.447583e-05*CosThetaK**2+3.758161e-05*CosThetaK**3+6.777260e-05*CosThetaK**4-5.585069e-05*CosThetaK**5-8.495213e-05*CosThetaK**6))",
    "4.625695e+04*((1.803216e-05+6.423635e-07*CosThetaL+9.704679e-06*CosThetaL**2+1.065779e-05*CosThetaL**3-1.658277e-06*CosThetaL**4-1.799046e-05*CosThetaL**5+6.089049e-06*CosThetaL**6)*(2.270524e-05-9.322913e-06*CosThetaK-1.587276e-05*CosThetaK**2+2.152708e-05*CosThetaK**3+5.615584e-05*CosThetaK**4-1.901528e-05*CosThetaK**5-4.887378e-05*CosThetaK**6))",
    "5.118383e+03*((1.668313e-04+1.911185e-05*CosThetaL+1.716389e-05*CosThetaL**2-3.192265e-05*CosThetaL**3+2.000329e-04*CosThetaL**4+1.783316e-05*CosThetaL**5-1.334724e-04*CosThetaL**6)*(2.056593e-04-4.151040e-05*CosThetaK-6.658669e-05*CosThetaK**2+3.742139e-05*CosThetaK**3+1.666491e-04*CosThetaK**4-5.072888e-05*CosThetaK**5-1.492963e-04*CosThetaK**6))",
    "3.837453e+03*((2.362599e-04-4.438020e-06*CosThetaL+3.318080e-05*CosThetaL**2+1.313482e-05*CosThetaL**3+7.878926e-05*CosThetaL**4-3.939653e-06*CosThetaL**5-2.211163e-05*CosThetaL**6)*(2.669904e-04-4.272653e-05*CosThetaK+1.487773e-05*CosThetaK**2+1.983652e-05*CosThetaK**3-9.317172e-05*CosThetaK**4-3.937610e-05*CosThetaK**5+4.831201e-05*CosThetaK**6))"
};
// Lumi = Nreco/(cross section*branch factor*filter efficiency), cross section is 49.59e9 [pb] for 8TeV and 48.44e9 [pb] for 7TeV.
// BF_BuToK*MuMu = 1.07E-6, 1.12E-6(2014)
// BF_BuToK*Jpsi = 1.43E-3, 1.44E-3(2014)
// BF_BuToK*Psi2S = 6.7E-4, 6.7 E-4(2014)
// BF_K*ToK0Pi  = 1  (K* decays to Kpi always)
// BF_K0ToKs  = 1/2
// BF_KsToPiPi = 2/3
double datasetLumi[5] = {19.4,26070.9,114.868,124.823,10.};//data, BuToKstarMuMu(11361.5+14709.4), BuToKstarJpsi(42.305+72.563), BuToKstarPsi2S(42.996,81.826), JpsiX
//}}}

double readParam(int iBin, const char parName[], int iColumn, double defVal=0., double forceReturn=999.)
{//{{{
    // Remark: first value is at iColumn=0.
    if (forceReturn != 999.) return forceReturn;

    std::vector<double> output;
    char lineBuff[1024];
    char *valBuff;
    memset(lineBuff,' ',1024*sizeof(char));
    FILE *fp = fopen(TString::Format("fitParameters%d.txt",iBin),"r");
    if (!fp){
        printf("WARNING: readParam, missing parameter files, by default return %f.",defVal);
        return defVal;
    }
    while(fgets(lineBuff,1024,fp) != NULL ){
        valBuff = strtok(lineBuff," ");
        if ( strcmp(valBuff,parName) == 0 ){
            printf("INFO: readParam, matched %s!\n",valBuff);
            valBuff = strtok(NULL," ");
            while(valBuff != NULL){
                //output.push_back(stof(valBuff));//stof if c++11 function, use other function
                if (strcmp(valBuff,"nan") == 0 || strcmp(valBuff,"inf") == 0 ){
                    output.push_back(defVal);
                }else{
                    output.push_back(std::atof(valBuff));
                }
                valBuff = strtok(NULL," ");
            }
            break;
        }
        memset(lineBuff,' ',1024*sizeof(char));
    }
    fclose(fp);
    
    if (iColumn < output.size() ){
        printf("INFO: readParam, get %s[%d]=%e\n",parName,iColumn,output.at(iColumn));
        return output.at(iColumn);
    }else{
        printf("WARNING: readParam, empty column! Return %s[%d]=defVal=%f.\n",parName,iColumn,defVal);
        return defVal;
    }
}//}}}
std::string readParam(int iBin, const char parName[], string defVal="", string forceReturn="defaultForceReturn")
{//{{{
    // Remark: first value is at iColumn=0.
    if (forceReturn != "defaultForceReturn") return forceReturn;

    string output;
    char lineBuff[1024];
    char *valBuff;
    memset(lineBuff,' ',1024*sizeof(char));
    FILE *fp = fopen(TString::Format("fitParameters%d.txt",iBin),"r");
    if (!fp){
        printf("WARNING: readParam, missing parameter files, by default return %s.",defVal.c_str());
        return defVal;
    }
    while(fgets(lineBuff,1024,fp) != NULL ){
        valBuff = strtok(lineBuff," ");
        if ( strcmp(valBuff,parName) == 0 ){
            printf("INFO: readParam, matched %s!\n",valBuff);
            valBuff = strtok(NULL,"\n");
            output=string(valBuff);
            break;
        }
        memset(lineBuff,' ',1024*sizeof(char));
    }
    fclose(fp);
    
    if (output != ""){
        printf("INFO: readParam, get %s=%s\n",parName,output.c_str());
        return output;
    }else{
        printf("WARNING: readParam, empty item! Return %s=defVal=%s.\n",parName, defVal.c_str());
        return defVal;
    }
}//}}}
void writeParam(int iBin, const char parName[], double *val, int nVal=2, bool overwrite=true)
{//{{{
    if ( !overwrite ) return;

    struct stat fiBuff;
    FILE *fi = 0;
    if (stat(TString::Format("fitParameters%d.txt",iBin),&fiBuff) == 0){
        rename(TString::Format("fitParameters%d.txt",iBin),TString::Format("fitParameters%d.txt.temp",iBin));
        fi = fopen(TString::Format("fitParameters%d.txt.temp",iBin),"r");
    }else{
        fi = fopen(TString::Format("fitParameters%d.txt.temp",iBin),"w");
    }
    
    bool parExist = false;
    char lineBuff[1024];
    char *valBuff = 0;
    memset(lineBuff,' ',1024*sizeof(char));
    FILE *fp = fopen(TString::Format("fitParameters%d.txt",iBin),"w");
    while(fgets(lineBuff,1024,fi) != NULL ){
        valBuff = strtok(lineBuff," ");
        if ( strcmp(valBuff,parName) == 0 ){
            fprintf(fp,"%s",parName);
            int iVal = 0;
            while(iVal < nVal){
                fprintf(fp," %e",val[iVal]);
                iVal++;
            }
            fprintf(fp,"\n");
            parExist = true;
        }else{
            fprintf(fp,"%s",lineBuff);
            valBuff = strtok(NULL," ");
            while( valBuff != NULL ){
                fprintf(fp," %s",valBuff);
                valBuff = strtok(NULL," ");
            }
        }
        memset(lineBuff,' ',512*sizeof(char));
    }
    if (parExist == false){
        fprintf(fp,"%s",parName);
        int iVal = 0;
        while(iVal < nVal){
            fprintf(fp," %e",val[iVal]);
            iVal++;
        }
        fprintf(fp,"\n");
    }
    fclose(fp);
    fclose(fi);
    remove(TString::Format("fitParameters%d.txt.temp",iBin));
}//}}}
void writeParam(int iBin, const char parName[], string instring, bool overwrite=true)
{//{{{
    if ( !overwrite ) return;

    struct stat fiBuff;
    FILE *fi = 0;
    if (stat(TString::Format("fitParameters%d.txt",iBin),&fiBuff) == 0){
        rename(TString::Format("fitParameters%d.txt",iBin),TString::Format("fitParameters%d.txt.temp",iBin));
        fi = fopen(TString::Format("fitParameters%d.txt.temp",iBin),"r");
    }else{
        fi = fopen(TString::Format("fitParameters%d.txt.temp",iBin),"w");
    }
    
    bool parExist = false;
    char lineBuff[1024];
    char *valBuff = 0;
    memset(lineBuff,' ',1024*sizeof(char));
    FILE *fp = fopen(TString::Format("fitParameters%d.txt",iBin),"w");
    while(fgets(lineBuff,1024,fi) != NULL ){
        valBuff = strtok(lineBuff," ");
        if ( strcmp(valBuff,parName) == 0 ){
            fprintf(fp,"%s %s\n", parName, instring.c_str());
            parExist = true;
        }else{
            fprintf(fp,"%s",lineBuff);
            valBuff = strtok(NULL," ");
            while( valBuff != NULL ){
                fprintf(fp," %s",valBuff);
                valBuff = strtok(NULL," ");
            }
        }
        memset(lineBuff,' ',512*sizeof(char));
    }
    if (parExist == false){
        fprintf(fp,"%s %s\n", parName, instring.c_str());
    }
    fclose(fp);
    fclose(fi);
    remove(TString::Format("fitParameters%d.txt.temp",iBin));
}//}}}

TF2  *f2_fcn = NULL;
double model_2D(double *x, double *par)
{//{{{
    double xx = x[0];
    double yy = x[1];
    for (int i = 0; i < f2_fcn->GetNpar(); i++) f2_fcn->SetParameter(i,par[i]);
    return f2_fcn->Eval(xx,yy);
}//}}}

TH2F *h2_fcn = NULL;
// nParameters, ???, return fcn value, parameter array, strategy
void fcn_binnedChi2_2D(int &npar, double *gin, double &f, double *par, int iflag)
{//{{{
    f=0;
    for (int i = 1; i <= h2_fcn->GetNbinsX(); i++) {
        for (int j = 1; j <= h2_fcn->GetNbinsY(); j++) {
            int gBin = h2_fcn->GetBin(i,j);
            double x[2] = {h2_fcn->GetXaxis()->GetBinCenter(i),h2_fcn->GetYaxis()->GetBinCenter(j)};
            double measure  = h2_fcn->GetBinContent(gBin);
            double error    = h2_fcn->GetBinError(gBin);
            
            //// Naively test using center value
            //double func     = model_2D(x, par);//Take center value
            //double delta    = (measure-func)/error;
            //if (measure != 0) 
            //    f+=delta*delta;
            
            //// Real run using integral
            for (int k = 0; k < f2_fcn->GetNpar(); k++){//nPar MUST be the same value as f2_fcn
                f2_fcn->SetParameter(k,par[k]);
            }
            double xi = h2_fcn->GetXaxis()->GetBinLowEdge(i);
            double xf = h2_fcn->GetXaxis()->GetBinUpEdge(i);
            double yi = h2_fcn->GetYaxis()->GetBinLowEdge(j);
            double yf = h2_fcn->GetYaxis()->GetBinUpEdge(j);
            //f2_fcn->SetRange(xi,xf,yi,yf);
            //double minX, minY;
            //f2_fcn->GetMinimumXY(minX,minY);
            //if (f2_fcn->Eval(minX,minY) < 0){
            //    f += 100;
            //}else{
                f += pow( (f2_fcn->Integral(xi,xf,yi,yf)/(xf-xi)/(yf-yi)-measure)/error,2);
            //}
        }
    }
    //printf("FCN in calls = %f\n",f);
}//}}}

//_________________________________________________________________________________

void bmass(int iBin, const char outfile[] = "bmass")
{//{{{
  bool test = false; 
  
  ch->SetBranchStatus("*",0);
  ch->SetBranchStatus("Bmass"         , 1);
  ch->SetBranchStatus("Mumumass"      , 1);
  ch->SetBranchStatus("Mumumasserr"   , 1);
  ch->SetBranchStatus("Q2"            , 1);
  RooRealVar Bmass("Bmass", "B^{+/-} mass(GeV/c^{2})", 5.27953-0.28, 5.27953+0.28) ;
  RooRealVar Q2("Q2","q^{2}",0.5,20.);
  RooRealVar Mumumass("Mumumass","M^{#mu#mu}",0.,10.);
  RooRealVar Mumumasserr("Mumumasserr","Error of M^{#mu#mu}",0.,10.);
  RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(Q2, Bmass, Mumumass, Mumumasserr),TString::Format("(%s) && (%s)",q2range[iBin],mumuMassWindow[1]),0);

  data->Print();

  // Create model and dataset
  // -------------------------------------------------------------------------
  // Gaussian signal 
  RooRealVar mean("mean","mean of gaussians", 5.27, 5.23, 5.32) ;
  RooRealVar sigma1("sigma1","width of Gaussian1", 0.0285, 0.01, 0.05) ;
  RooRealVar sigma2("sigma2","width of Gaussian2", 0.08, 0.05, 0.35) ;
  RooRealVar sigM_frac("sigM_frac","fraction of Gaussians",0,0.,1.);
  RooGaussian sigGauss1("siggauss1","Signal component", Bmass, mean, sigma1) ;  
  RooGaussian sigGauss2("siggauss2","Signal component", Bmass, mean, sigma2) ;  
  RooAddPdf sig("sig","sig",RooArgList(sigGauss1,sigGauss2),RooArgList(sigM_frac));

  // Build Chebychev polynomial p.d.f.  
  RooRealVar a0("a0", "constant", 0.5, -1, 1) ;
  RooRealVar a1("a1", "linear", 0.6, -1, 1) ;
  RooRealVar a2("a2", "quadratic", 0.1, -1, 1) ;
  RooChebychev bkg("bkg", "Background", Bmass, RooArgSet(a0, a1, a2)) ;

  // Construct signal+background PDF
  RooRealVar nsig("nsig", "number of signal events", 4648, 0, 1E8); 
  RooRealVar nbkg("nbkg", "number of background events", 21472, 0, 1E8);
  RooAddPdf  model("model", "g+c", RooArgList(bkg, sig), RooArgList(nbkg, nsig)) ;
  
  // Print structure of composite p.d.f.
  model.Print("t") ;

  // Fit model to data, save fitresult 
  // ------------------------------------------------------------------------
  RooFitResult* fitres; 
  if (! test) {
    fitres = model.fitTo(*data, Extended(kTRUE), Minos(kTRUE), Save(kTRUE)) ;
    fitres->Print("v"); 
  }
  
  // Plot model 
  // ---------------------------------------------------------
  TString title = "B^{+/-} mass";
  int nbins = 20; 
  RooPlot* frame = Bmass.frame(Title(title), Bins(nbins));
  data->plotOn(frame) ;
  model.plotOn(frame) ;

  // Overlay the background component of model with a dashed line
  model.plotOn(frame,Components("bkg"), LineStyle(kDashed), LineColor(2)) ;

  // Draw the frame on the canvas
  TCanvas *c = new TCanvas("c"); 
  set_root_style(); 
  c->UseCurrentStyle() ;

  gPad->SetLeftMargin(0.15) ;
  frame->GetYaxis()->SetTitleOffset(1.7) ; 
  frame->Draw();

  TPaveText* paveText = new TPaveText(0.75, 0.82, 1., 1., "NDC"); 
  paveText->SetBorderSize(0);
  paveText->SetFillColor(kWhite);
  paveText->AddText(Form("nsig = %.0f #pm %.0f " , nsig.getVal()  , nsig.getError())); 
  paveText->AddText(Form("nbkg = %.0f #pm %.0f " , nbkg.getVal()  , nbkg.getError())); 
  paveText->AddText(Form("mean = %.3f #pm %.3f " , mean.getVal()  , mean.getError())); 
  paveText->AddText(Form("sigma1 = %.3f #pm %.3f ", sigma1.getVal(), sigma1.getError())); 
  paveText->AddText(Form("sigma2 = %.3f #pm %.3f ", sigma2.getVal(), sigma2.getError())); 
  paveText->AddText(Form("frac = %.3f #pm %.3f ", sigM_frac.getVal(), sigM_frac.getError())); 
  paveText->Draw(); 

  c->Print(TString::Format("./plots/%s_bin%d.pdf",outfile,iBin));
  
  // Persist fit result in root file 
  // -------------------------------------------------------------
  //TFile resf(TString::Format("./plots/%s.root",outfile), "RECREATE") ;
  //gPad->Write("plot"); 
  //if (! test) fitres->Write("fitres") ;
  //resf.Close() ;

  // In a clean ROOT session retrieve the persisted fit result as follows:
  // RooFitResult* r = gDirectory->Get("fitres") ;
  
  delete paveText; 
  delete c;

}//}}}

void angular3D_bin(int iBin, const char outfile[] = "angular3D")
{//{{{
    // Remark: You must use RooFit!! It's better in unbinned fit.
    //         Extended ML fit is adopted by Mauro, just follow!!
    //         This function is NOT designed to determine parameters in resonance region.
    if (iBin==3 || iBin==5) return;
    
    // Read data
    RooRealVar CosThetaK("CosThetaK"     , "cos#theta_{K}"       , -1. , 1.   ) ;
    RooRealVar CosThetaL("CosThetaL"     , "cos#theta_{L}"       , -1. , 1.   ) ;
    RooRealVar Bmass("Bmass"             , "M_{K^{*}#mu#mu}"     , 5.  , 5.56 ) ;
    RooRealVar Mumumass("Mumumass"       , "M^{#mu#mu}"          , 0.  , 10.  ) ;
    RooRealVar Mumumasserr("Mumumasserr" , "Error of M^{#mu#mu}" , 0.  , 10.  ) ;
    RooRealVar Q2("Q2"                   , "q^{2}"               , 0.5 , 20.  ) ;
    RooAddition Bmass_offset("Bmass_offset","Bmass_offset",RooArgSet(Bmass,RooConst(-5)));
    RooProduct Bmass_norm("Bmass_norm","Bmass_norm",RooArgSet(Bmass_offset,RooConst(1./0.56)));
        
    // Angular parameters
    RooRealVar afb("afb" , "A_{FB}" , genAfb[iBin] , -1. , 1.);
    RooRealVar fl("fl"   , "F_{L}"  , genFl[iBin]  , 0.  , 1.);
    RooRealVar fs("fs"   , "F_{S}"  , 0 , 0.  , 1.);//Derive from B0ToKstarJpsi , Bin3
    RooRealVar as("as"   , "A_{S}"  , 0 , -1. , 1.);//Derive from B0ToKstarJpsi , Bin3
    if (iBin != 3 && iBin != 5){
        // 2011 cross check
        //fs.setVal(0.0129254);
        //fs.setAsymError(-0.00898344,0.0101371);
        //as.setVal(-0.0975919);
        //as.setAsymError(-0.00490805,0.0049092);

        // read parameter from datacard
        fs.setVal(readParam(3,"fs",0, 0.0129254));
        fs.setAsymError(readParam(3,"fs",1, -0.00898344),readParam(3,"fs",2, 0.0101371));
        as.setVal(readParam(3,"as",0, -0.0975919));
        as.setAsymError(readParam(3,"as",1, -0.00490805),readParam(3,"as",2, 0.0049092));
    }

    // Create parameters and PDFs
        // Signal double gaussian
    RooRealVar sigGauss_mean("sigGauss_mean","M_{K*#Mu#Mu}",5.28,5.25,5.30);
    RooRealVar sigGauss1_sigma("sigGauss1_sigma","#sigma_{1}",readParam(iBin,"sigGauss1_sigma",0,0.02),.01,.05);
    RooRealVar sigGauss2_sigma("sigGauss2_sigma","#sigma_{2}",readParam(iBin,"sigGauss2_sigma",0,0.08),.05,.40);
    sigGauss1_sigma.setError(readParam(iBin,"sigGauss1_sigma",1));
    sigGauss2_sigma.setError(readParam(iBin,"sigGauss2_sigma",1));
    RooGaussian f_sigMGauss1("f_sigMGauss1","f_sigMGauss1", Bmass, sigGauss_mean, sigGauss1_sigma);//double gaussian with shared mean
    RooGaussian f_sigMGauss2("f_sigMGauss2","f_sigMGauss2", Bmass, sigGauss_mean, sigGauss2_sigma);//double gaussian with shared mean
    RooRealVar sigM_frac("sigM_frac","sigM_frac",readParam(iBin,"sigM_frac",0,0.5),0.,1.);
    sigM_frac.setError(readParam(iBin,"sigM_frac",1));
    RooAddPdf f_sigM("f_sigM","f_sigM", RooArgList(f_sigMGauss1, f_sigMGauss2), sigM_frac);

    // Efficiency, acceptance and angular distribution of signal (checked)
    RooRealVar recK0L0("recK0L0","recK0L0",readParam(iBin,"accXrecoEff2", 0));
    RooRealVar recK1L0("recK1L0","recK1L0",readParam(iBin,"accXrecoEff2", 1));
    RooRealVar recK2L0("recK2L0","recK2L0",readParam(iBin,"accXrecoEff2", 2));
    RooRealVar recK3L0("recK3L0","recK3L0",readParam(iBin,"accXrecoEff2", 3));
    RooRealVar recK0L2("recK0L2","recK0L2",readParam(iBin,"accXrecoEff2", 4));
    RooRealVar recK1L2("recK1L2","recK1L2",readParam(iBin,"accXrecoEff2", 5));
    RooRealVar recK2L2("recK2L2","recK2L2",readParam(iBin,"accXrecoEff2", 6));
    RooRealVar recK3L2("recK3L2","recK3L2",readParam(iBin,"accXrecoEff2", 7));
    RooRealVar recK0L3("recK0L3","recK0L3",readParam(iBin,"accXrecoEff2", 8));
    RooRealVar recK1L3("recK1L3","recK1L3",readParam(iBin,"accXrecoEff2", 9));
    RooRealVar recK2L3("recK2L3","recK2L3",readParam(iBin,"accXrecoEff2",10));
    RooRealVar recK3L3("recK3L3","recK3L3",readParam(iBin,"accXrecoEff2",11));
    RooRealVar recK0L4("recK0L4","recK0L4",readParam(iBin,"accXrecoEff2",12));
    RooRealVar recK1L4("recK1L4","recK1L4",readParam(iBin,"accXrecoEff2",13));
    RooRealVar recK2L4("recK2L4","recK2L4",readParam(iBin,"accXrecoEff2",14));
    RooRealVar recK3L4("recK3L4","recK3L4",readParam(iBin,"accXrecoEff2",15));
    RooRealVar recK0L6("recK0L6","recK0L6",readParam(iBin,"accXrecoEff2",16));
    RooRealVar recK1L6("recK1L6","recK1L6",readParam(iBin,"accXrecoEff2",17));
    RooRealVar recK2L6("recK2L6","recK2L6",readParam(iBin,"accXrecoEff2",18));
    RooRealVar recK3L6("recK3L6","recK3L6",readParam(iBin,"accXrecoEff2",19));
    recK0L0.setError(readParam(iBin,"accXrecoEff2Err", 0));
    recK1L0.setError(readParam(iBin,"accXrecoEff2Err", 1));
    recK2L0.setError(readParam(iBin,"accXrecoEff2Err", 2));
    recK3L0.setError(readParam(iBin,"accXrecoEff2Err", 3));
    recK0L2.setError(readParam(iBin,"accXrecoEff2Err", 4));
    recK1L2.setError(readParam(iBin,"accXrecoEff2Err", 5));
    recK2L2.setError(readParam(iBin,"accXrecoEff2Err", 6));
    recK3L2.setError(readParam(iBin,"accXrecoEff2Err", 7));
    recK0L3.setError(readParam(iBin,"accXrecoEff2Err", 8));
    recK1L3.setError(readParam(iBin,"accXrecoEff2Err", 9));
    recK2L3.setError(readParam(iBin,"accXrecoEff2Err",10));
    recK3L3.setError(readParam(iBin,"accXrecoEff2Err",11));
    recK0L4.setError(readParam(iBin,"accXrecoEff2Err",12));
    recK1L4.setError(readParam(iBin,"accXrecoEff2Err",13));
    recK2L4.setError(readParam(iBin,"accXrecoEff2Err",14));
    recK3L4.setError(readParam(iBin,"accXrecoEff2Err",15));
    recK0L6.setError(readParam(iBin,"accXrecoEff2Err",16));
    recK1L6.setError(readParam(iBin,"accXrecoEff2Err",17));
    recK2L6.setError(readParam(iBin,"accXrecoEff2Err",18));
    recK3L6.setError(readParam(iBin,"accXrecoEff2Err",19));
    RooArgSet f_sigA_argset(CosThetaL,CosThetaK);
    f_sigA_argset.add(RooArgSet(fl,afb,fs,as));
    f_sigA_argset.add(RooArgSet(recK0L0,recK1L0,recK2L0,recK3L0));
    f_sigA_argset.add(RooArgSet(recK0L2,recK1L2,recK2L2,recK3L2));
    TString f_sigA_format;
    TString f_ang_format = "9/16*((2/3*fs+4/3*as*CosThetaK)*(1-CosThetaL**2)+(1-fs)*(2*fl*CosThetaK**2*(1-CosThetaL**2)+1/2*(1-fl)*(1-CosThetaK**2)*(1+CosThetaL**2)+4/3*afb*(1-CosThetaK**2)*CosThetaL))";
    TString f_rec_ord0 = readParam(iBin,"f_accXrecoEff_ord0",f_accXrecoEff_ord0[iBin]);
    TString f_rec_format, f_rec_L0, f_rec_L2, f_rec_L3, f_rec_L4, f_rec_L6;
    f_rec_L0 = "(recK0L0+recK1L0*CosThetaK+recK2L0*(3*CosThetaK**2-1)/2+recK3L0*(5*CosThetaK**3-3*CosThetaK)/2)";
    f_rec_L2 = "(recK0L2+recK1L2*CosThetaK+recK2L2*(3*CosThetaK**2-1)/2+recK3L2*(5*CosThetaK**3-3*CosThetaK)/2)*CosThetaL**2";
    f_rec_L3 = "(recK0L3+recK1L3*CosThetaK+recK2L3*(3*CosThetaK**2-1)/2+recK3L3*(5*CosThetaK**3-3*CosThetaK)/2)*CosThetaL**3";
    f_rec_L4 = "(recK0L4+recK1L4*CosThetaK+recK2L4*(3*CosThetaK**2-1)/2+recK3L4*(5*CosThetaK**3-3*CosThetaK)/2)*CosThetaL**4";
    f_rec_L6 = "(recK0L6+recK1L6*CosThetaK+recK2L6*(3*CosThetaK**2-1)/2+recK3L6*(5*CosThetaK**3-3*CosThetaK)/2)*CosThetaL**6";

    if (iBin == 0) {
        f_sigA_argset.add(RooArgSet(recK0L4,recK1L4,recK2L4,recK3L4));
        f_sigA_argset.add(RooArgSet(recK0L6,recK1L6,recK2L6,recK3L6));
        f_sigA_format = TString::Format("%s*(1+%s+%s+%s+%s)*%s",f_rec_ord0.Data(),f_rec_L0.Data(),f_rec_L2.Data(),f_rec_L4.Data(),f_rec_L6.Data(),f_ang_format.Data());
    }else if (iBin == 1) {
        f_sigA_argset.add(RooArgSet(recK0L4,recK1L4,recK2L4,recK3L4));
        f_sigA_format = TString::Format("%s*(1+%s+%s+%s)*%s",f_rec_ord0.Data(),f_rec_L0.Data(),f_rec_L2.Data(),f_rec_L4.Data(),f_ang_format.Data());
    }else if (iBin > 1 && iBin < 6) {
        f_sigA_argset.add(RooArgSet(recK0L3,recK1L3,recK2L3,recK3L3));
        f_sigA_argset.add(RooArgSet(recK0L4,recK1L4,recK2L4,recK3L4));
        f_sigA_format = TString::Format("%s*(1+%s+%s+%s+%s)*%s",f_rec_ord0.Data(),f_rec_L0.Data(),f_rec_L2.Data(),f_rec_L3.Data(),f_rec_L4.Data(),f_ang_format.Data());
    }else{
        f_sigA_argset.add(RooArgSet(recK0L3,recK1L3,recK2L3,recK3L3));
        f_sigA_format = TString::Format("%s*(1+%s+%s+%s)*%s",f_rec_ord0.Data(),f_rec_L0.Data(),f_rec_L2.Data(),f_rec_L3.Data(),f_ang_format.Data());
    }
    RooGenericPdf f_sigA("f_sigA", f_sigA_format, f_sigA_argset);

        // merge mass and angular distro of signal
    RooProdPdf f_sig("f_sig","f_sig",RooArgSet(f_sigM,f_sigA));
    printf("INFO: f_sig prepared.\n");
    
    // Create combinatorial background distribution (to be checked)
    RooRealVar bkgCombM_c1("bkgCombM_c1","c1",readParam(iBin,"bkgCombM_c1",0),-100,100);
    RooRealVar bkgCombM_c2("bkgCombM_c2","c2",readParam(iBin,"bkgCombM_c2",0),-100,100);
    RooRealVar bkgCombM_c3("bkgCombM_c3","c3",readParam(iBin,"bkgCombM_c3",0),-100,100);
    RooRealVar bkgCombM_c4("bkgCombM_c4","c4",readParam(iBin,"bkgCombM_c4",0),-100,100);
    RooRealVar bkgCombM_c5("bkgCombM_c5","c5",readParam(iBin,"bkgCombM_c5",0),-100,100);
    RooRealVar bkgCombM_c6("bkgCombM_c6","c6",readParam(iBin,"bkgCombM_c6",0),-100,100);
    RooRealVar bkgCombM_c7("bkgCombM_c7","c7",readParam(iBin,"bkgCombM_c7",0),-100,100);
    RooArgSet f_bkgCombM_argset;
    switch (iBin){
        case 0:
        case 1:
        case 2:
        case 3:
        case 4:
        case 5:
        case 6:
        case 7:
            f_bkgCombM_argset.add(RooArgSet(bkgCombM_c1,bkgCombM_c2,bkgCombM_c3,bkgCombM_c4));
            break;
        default:
            f_bkgCombM_argset.add(RooArgSet(bkgCombM_c1,bkgCombM_c2,bkgCombM_c3,bkgCombM_c4));
            break;
    }
    RooPolynomial f_bkgCombM("f_bkgCombM","f_bkgCombM",Bmass_norm,f_bkgCombM_argset);
    RooRealVar bkgCombL_c1("bkgCombL_c1","c1",readParam(iBin,"bkgCombL_c1",0),-100,100);
    RooRealVar bkgCombL_c2("bkgCombL_c2","c2",readParam(iBin,"bkgCombL_c2",0.01),-100,100);
    RooRealVar bkgCombL_c3("bkgCombL_c3","c3",readParam(iBin,"bkgCombL_c3",0),-100,100);
    RooRealVar bkgCombL_c4("bkgCombL_c4","c4",readParam(iBin,"bkgCombL_c4",0.01),-100,100);
    RooRealVar bkgCombL_c5("bkgCombL_c5","c5",readParam(iBin,"bkgCombL_c5",0),-100,100);
    RooRealVar bkgCombL_c6("bkgCombL_c6","c6",readParam(iBin,"bkgCombL_c6",0.01),-100,100);
    RooRealVar bkgCombL_c7("bkgCombL_c7","c7",readParam(iBin,"bkgCombL_c7",0),-100,100);
    RooRealVar bkgCombL_c8("bkgCombL_c8","c8",readParam(iBin,"bkgCombL_c8",0),-100,100);
    RooArgSet f_bkgCombL_argset;
    switch (iBin) {
        case 0:
        case 1:
        case 2:
        case 3:
        case 4:
        case 5:
        case 6:
        case 7:
            f_bkgCombL_argset.add(RooArgSet(bkgCombL_c1,bkgCombL_c2,bkgCombL_c3,bkgCombL_c4));
            break;
        default:
            bkgCombL_c1.setVal(0.);
            bkgCombL_c2.setVal(0.);
            bkgCombL_c3.setVal(0.);
            bkgCombL_c4.setVal(0.);
            bkgCombL_c1.setConstant(kTRUE);
            bkgCombL_c2.setConstant(kTRUE);
            bkgCombL_c3.setConstant(kTRUE);
            bkgCombL_c4.setConstant(kTRUE);
            break;
    }
    RooPolynomial f_bkgCombL("f_bkgCombL","f_bkgCombL",CosThetaL,f_bkgCombL_argset); // RooPolynomial
//    RooGaussian f_bkgCombL_gaus1("f_bkgCombL_gaus1","f_bkgCombL_gaus1",CosThetaL,bkgCombL_c1,bkgCombL_c2);
//    RooGaussian f_bkgCombL_gaus2("f_bkgCombL_gaus2","f_bkgCombL_gaus2",CosThetaL,bkgCombL_c3,bkgCombL_c4);
//    RooGaussian f_bkgCombL_gaus3("f_bkgCombL_gaus3","f_bkgCombL_gaus3",CosThetaL,bkgCombL_c5,bkgCombL_c6);
//    RooAddPdf f_bkgCombL("f_bkgCombL","f_bkgCombL",RooArgList(f_bkgCombL_gaus1,f_bkgCombL_gaus2,f_bkgCombL_gaus3),RooArgList(bkgCombL_c7,bkgCombL_c8));//3 Gaussians
    RooRealVar bkgCombK_c1("bkgCombK_c1","c1",readParam(iBin,"bkgCombK_c1",0),-100,100);
    RooRealVar bkgCombK_c2("bkgCombK_c2","c2",readParam(iBin,"bkgCombK_c2",0),-100,100);
    RooRealVar bkgCombK_c3("bkgCombK_c3","c3",readParam(iBin,"bkgCombK_c3",0),-100,100);
    RooRealVar bkgCombK_c4("bkgCombK_c4","c4",readParam(iBin,"bkgCombK_c4",0),-100,100);
    RooRealVar bkgCombK_c5("bkgCombK_c5","c5",readParam(iBin,"bkgCombK_c5",0),-100,100);
    RooRealVar bkgCombK_c6("bkgCombK_c6","c6",readParam(iBin,"bkgCombK_c6",0),-100,100);
    RooRealVar bkgCombK_c7("bkgCombK_c7","c7",readParam(iBin,"bkgCombK_c7",0),-100,100);
    RooArgSet f_bkgCombK_argset;
    switch (iBin) {
        case 0:
        case 1:
        case 2:
        case 3:
        case 4:
        case 5:
        case 6:
        case 7:
            //f_bkgCombK_argset.add(RooArgSet(bkgCombK_c1,bkgCombK_c2,bkgCombK_c3,bkgCombK_c4,bkgCombK_c5,bkgCombK_c6,bkgCombK_c7));
            f_bkgCombK_argset.add(RooArgSet(bkgCombK_c1,bkgCombK_c2,bkgCombK_c3,bkgCombK_c4));
            break;
        default:
            bkgCombK_c1.setVal(0.);
            bkgCombK_c2.setVal(0.);
            bkgCombK_c3.setVal(0.);
            bkgCombK_c4.setVal(0.);
            bkgCombK_c1.setConstant(kTRUE);
            bkgCombK_c2.setConstant(kTRUE);
            bkgCombK_c3.setConstant(kTRUE);
            bkgCombK_c4.setConstant(kTRUE);
            break;
    }
    RooPolynomial f_bkgCombK("f_bkgCombK","f_bkgCombK",CosThetaK,f_bkgCombK_argset);
    RooProdPdf f_bkgCombA("f_bkgCombA", "f_bckCombA",RooArgSet(f_bkgCombK,f_bkgCombL));
    RooProdPdf f_bkgComb("f_bkgComb", "f_bckComb",RooArgSet(f_bkgCombA,f_bkgCombM));
    
    //TString f_bkgComb_format = "";
    //RooArgSet f_bkgComb_argset;
    //f_bkgComb_argset.add(RooArgSet(bkgCombK_c1,bkgCombK_c2,bkgCombK_c3,bkgCombK_c4,CosThetaK));
    //f_bkgComb_argset.add(RooArgSet(bkgCombL_c1,bkgCombL_c2,bkgCombL_c3,bkgCombL_c4,CosThetaL));
    //f_bkgComb_argset.add(RooArgSet(bkgCombM_c1,bkgCombM_c2,bkgCombM_c3,bkgCombM_c4,Bmass_norm));
    //RooGenericPdf f_bkgComb("f_bkgComb", f_bkgComb_format,f_bkgComb_argset);
    printf("INFO: f_bkgComb prepared.\n");
    
    // Create peak background distribution(jpsi/psi2s)
    RooRealVar bkgjpsiGauss1_mean1("bkgjpsiGauss1_mean1","M_{K*#Mu#Mu}",readParam(iBin,"bkgjpsiGauss1_mean1",0));
    RooRealVar bkgjpsiGauss1_mean2("bkgjpsiGauss1_mean2","M_{K*#Mu#Mu}",readParam(iBin,"bkgjpsiGauss1_mean2",0));
    RooRealVar bkgjpsiGauss1_sigma1("bkgjpsiGauss1_sigma1","#sigma_{11}",readParam(iBin,"bkgjpsiGauss1_sigma1",0));
    RooRealVar bkgjpsiGauss1_sigma2("bkgjpsiGauss1_sigma2","#sigma_{12}",readParam(iBin,"bkgjpsiGauss1_sigma2",0));
    RooRealVar bkgjpsiM_frac1("bkgjpsiM_frac1","bkgjpsiM_frac1",readParam(iBin,"bkgjpsiM_frac1",0));
    RooRealVar bkgjpsiM_fracC1("bkgjpsiM_fracC1"           , "bkgjpsiM_fracC1"  , readParam(iBin,"bkgjpsiM_fracC1",0));
    RooRealVar bkgjpsiGauss2_mean1("bkgjpsiGauss2_mean1","M_{K*#Mu#Mu}",readParam(iBin,"bkgjpsiGauss2_mean1",0));
    RooRealVar bkgjpsiGauss2_mean2("bkgjpsiGauss2_mean2","M_{K*#Mu#Mu}",readParam(iBin,"bkgjpsiGauss2_mean2",0));
    RooRealVar bkgjpsiGauss2_sigma1("bkgjpsiGauss2_sigma1","#sigma_{21}",readParam(iBin,"bkgjpsiGauss2_sigma1",0));
    RooRealVar bkgjpsiGauss2_sigma2("bkgjpsiGauss2_sigma2","#sigma_{22}",readParam(iBin,"bkgjpsiGauss2_sigma2",0));
    RooRealVar bkgjpsiM_frac2("bkgjpsiM_frac2","bkgjpsiM_frac2",readParam(iBin,"bkgjpsiM_frac2",0));
    RooRealVar bkgjpsiM_frac12("bkgjpsiM_frac12","bkgjpsiM_frac12",readParam(iBin,"bkgjpsiM_frac12",0));
    bkgjpsiGauss1_mean1.setError(readParam(iBin,"bkgjpsiGauss1_mean1",1));
    bkgjpsiGauss1_mean2.setError(readParam(iBin,"bkgjpsiGauss1_mean2",1));
    bkgjpsiGauss1_sigma1.setError(readParam(iBin,"bkgjpsiGauss1_sigma1",1));
    bkgjpsiGauss1_sigma2.setError(readParam(iBin,"bkgjpsiGauss1_sigma2",1));
    bkgjpsiM_frac1.setError(readParam(iBin,"bkgjpsiM_frac1",1));
    bkgjpsiM_fracC1.setError(readParam(iBin,"bkgjpsiM_fracC1",1));
    bkgjpsiGauss2_mean1.setError(readParam(iBin,"bkgjpsiGauss2_mean1",1));
    bkgjpsiGauss2_mean2.setError(readParam(iBin,"bkgjpsiGauss2_mean2",1));
    bkgjpsiGauss2_sigma1.setError(readParam(iBin,"bkgjpsiGauss2_sigma1",1));
    bkgjpsiGauss2_sigma2.setError(readParam(iBin,"bkgjpsiGauss2_sigma2",1));
    bkgjpsiM_frac2.setError(readParam(iBin,"bkgjpsiM_frac2",1));
    bkgjpsiM_frac12.setError(readParam(iBin,"bkgjpsiM_frac12",1));
    RooGaussian f_bkgjpsiPeakMGauss11("f_bkgjpsiPeakMGauss11","f_bkgjpsiPeakMGauss11", Bmass, bkgjpsiGauss1_mean1, bkgjpsiGauss1_sigma1);
    RooGaussian f_bkgjpsiPeakMGauss12("f_bkgjpsiPeakMGauss12","f_bkgjpsiPeakMGauss12", Bmass, bkgjpsiGauss1_mean2, bkgjpsiGauss1_sigma2);
    RooGaussian f_bkgjpsiPeakMGauss21("f_bkgjpsiPeakMGauss21","f_bkgjpsiPeakMGauss21", Bmass, bkgjpsiGauss2_mean1, bkgjpsiGauss2_sigma1);
    RooGaussian f_bkgjpsiPeakMGauss22("f_bkgjpsiPeakMGauss22","f_bkgjpsiPeakMGauss22", Bmass, bkgjpsiGauss2_mean2, bkgjpsiGauss2_sigma2);
    RooAddPdf f_bkgjpsiPeakM1("f_bkgjpsiPeakM1","f_bkgjpsiPeakM1", RooArgList(f_bkgjpsiPeakMGauss11, f_bkgjpsiPeakMGauss12), bkgjpsiM_frac1);
    RooAddPdf f_bkgjpsiPeakM2("f_bkgjpsiPeakM2","f_bkgjpsiPeakM2", RooArgList(f_bkgjpsiPeakMGauss21, f_bkgjpsiPeakMGauss22), bkgjpsiM_frac2);
    RooAddPdf f_bkgjpsiPeakM12("f_bkgjpsiPeakM12","f_bkgjpsiPeakM12", RooArgList(f_bkgjpsiPeakM1, f_bkgjpsiPeakM2), bkgjpsiM_frac12);
    RooRealVar bkgpsi2sGauss1_mean1("bkgpsi2sGauss1_mean1","M_{K*#Mu#Mu}",readParam(iBin,"bkgpsi2sGauss1_mean1",0));
    RooRealVar bkgpsi2sGauss1_mean2("bkgpsi2sGauss1_mean2","M_{K*#Mu#Mu}",readParam(iBin,"bkgpsi2sGauss1_mean2",0));
    RooRealVar bkgpsi2sGauss1_sigma1("bkgpsi2sGauss1_sigma1","#sigma_{11}",readParam(iBin,"bkgpsi2sGauss1_sigma1",0));
    RooRealVar bkgpsi2sGauss1_sigma2("bkgpsi2sGauss1_sigma2","#sigma_{12}",readParam(iBin,"bkgpsi2sGauss1_sigma2",0));
    RooRealVar bkgpsi2sM_frac1("bkgpsi2sM_frac1","bkgpsi2sM_frac1",readParam(iBin,"bkgpsi2sM_frac1",0));
    RooRealVar bkgpsi2sM_fracC1("bkgpsi2sM_fracC1"           , "bkgpsi2sM_fracC1"  , readParam(iBin,"bkgpsi2sM_fracC1",0));
    RooRealVar bkgpsi2sGauss2_mean1("bkgpsi2sGauss2_mean1","M_{K*#Mu#Mu}",readParam(iBin,"bkgpsi2sGauss2_mean1",0));
    RooRealVar bkgpsi2sGauss2_mean2("bkgpsi2sGauss2_mean2","M_{K*#Mu#Mu}",readParam(iBin,"bkgpsi2sGauss2_mean2",0));
    RooRealVar bkgpsi2sGauss2_sigma1("bkgpsi2sGauss2_sigma1","#sigma_{21}",readParam(iBin,"bkgpsi2sGauss2_sigma1",0));
    RooRealVar bkgpsi2sGauss2_sigma2("bkgpsi2sGauss2_sigma2","#sigma_{22}",readParam(iBin,"bkgpsi2sGauss2_sigma2",0));
    RooRealVar bkgpsi2sM_frac2("bkgpsi2sM_frac2","bkgpsi2sM_frac2",readParam(iBin,"bkgpsi2sM_frac2",0));
    RooRealVar bkgpsi2sM_frac12("bkgpsi2sM_frac12","bkgpsi2sM_frac12",readParam(iBin,"bkgpsi2sM_frac12",0));
    bkgpsi2sGauss1_mean1.setError(readParam(iBin,"bkgpsi2sGauss1_mean1",1));
    bkgpsi2sGauss1_mean2.setError(readParam(iBin,"bkgpsi2sGauss1_mean2",1));
    bkgpsi2sGauss1_sigma1.setError(readParam(iBin,"bkgpsi2sGauss1_sigma1",1));
    bkgpsi2sGauss1_sigma2.setError(readParam(iBin,"bkgpsi2sGauss1_sigma2",1));
    bkgpsi2sM_frac1.setError(readParam(iBin,"bkgpsi2sM_frac1",1));
    bkgpsi2sM_fracC1.setError(readParam(iBin,"bkgpsi2sM_fracC1",1));
    bkgpsi2sGauss2_mean1.setError(readParam(iBin,"bkgpsi2sGauss2_mean1",1));
    bkgpsi2sGauss2_mean2.setError(readParam(iBin,"bkgpsi2sGauss2_mean2",1));
    bkgpsi2sGauss2_sigma1.setError(readParam(iBin,"bkgpsi2sGauss2_sigma1",1));
    bkgpsi2sGauss2_sigma2.setError(readParam(iBin,"bkgpsi2sGauss2_sigma2",1));
    bkgpsi2sM_frac2.setError(readParam(iBin,"bkgpsi2sM_frac2",1));
    bkgpsi2sM_frac12.setError(readParam(iBin,"bkgpsi2sM_frac12",1));
    RooGaussian f_bkgpsi2sPeakMGauss11("f_bkgpsi2sPeakMGauss11","f_bkgpsi2sPeakMGauss11", Bmass, bkgpsi2sGauss1_mean1, bkgpsi2sGauss1_sigma1);
    RooGaussian f_bkgpsi2sPeakMGauss12("f_bkgpsi2sPeakMGauss12","f_bkgpsi2sPeakMGauss12", Bmass, bkgpsi2sGauss1_mean2, bkgpsi2sGauss1_sigma2);
    RooGaussian f_bkgpsi2sPeakMGauss21("f_bkgpsi2sPeakMGauss21","f_bkgpsi2sPeakMGauss21", Bmass, bkgpsi2sGauss2_mean1, bkgpsi2sGauss2_sigma1);
    RooGaussian f_bkgpsi2sPeakMGauss22("f_bkgpsi2sPeakMGauss22","f_bkgpsi2sPeakMGauss22", Bmass, bkgpsi2sGauss2_mean2, bkgpsi2sGauss2_sigma2);
    RooAddPdf f_bkgpsi2sPeakM1("f_bkgpsi2sPeakM1","f_bkgpsi2sPeakM1", RooArgList(f_bkgpsi2sPeakMGauss11, f_bkgpsi2sPeakMGauss12), bkgpsi2sM_frac1);
    RooAddPdf f_bkgpsi2sPeakM2("f_bkgpsi2sPeakM2","f_bkgpsi2sPeakM2", RooArgList(f_bkgpsi2sPeakMGauss21, f_bkgpsi2sPeakMGauss22), bkgpsi2sM_frac2);
    RooAddPdf f_bkgpsi2sPeakM12("f_bkgpsi2sPeakM12","f_bkgpsi2sPeakM12", RooArgList(f_bkgpsi2sPeakM1, f_bkgpsi2sPeakM2), bkgpsi2sM_frac12);
    switch (iBin) {// Should be constants.
        case 2:
            // double Gaussian
            bkgjpsiM_frac12.setVal(1.);
            bkgjpsiM_frac12.setConstant(kTRUE);
            bkgpsi2sM_frac12.setVal(1.);
            bkgpsi2sM_frac12.setConstant(kTRUE);
//            f = new RooExtendPdf("f","f",f_bkgPeakM1,nbkgPeak); // triple Gaussian at low
            break;
        case 4:
            // two double Gaussians
            bkgjpsiM_fracC1.setVal(0.);
            bkgjpsiM_fracC1.setConstant(kTRUE);
            bkgjpsiM_frac12.setVal(0.);
            bkgpsi2sM_fracC1.setVal(0.);
            bkgpsi2sM_fracC1.setConstant(kTRUE);
            bkgpsi2sM_frac12.setVal(1.);
            bkgpsi2sM_frac12.setConstant(kTRUE);
//            f = new RooExtendPdf("f","f",f_bkgPeakM12,nbkgPeak); // two double Gaussian
            break;
        case 6:
            // single Gaussian
            bkgjpsiM_frac12.setVal(0.);
            bkgjpsiM_frac12.setConstant(kTRUE);
            bkgpsi2sM_frac12.setVal(0.);
            bkgpsi2sM_frac12.setConstant(kTRUE);
//            f = new RooExtendPdf("f","f",f_bkgPeakM2,nbkgPeak); // single Gaussian at high
            break;
        case 3:
        case 5:
            bkgjpsiM_frac12.setVal(1.);
            bkgjpsiM_frac12.setConstant(kTRUE);
            bkgjpsiM_fracC1.setVal(0.);
            bkgjpsiM_fracC1.setConstant(kTRUE);
            bkgpsi2sM_frac12.setVal(1.);
            bkgpsi2sM_frac12.setConstant(kTRUE);
            bkgpsi2sM_fracC1.setVal(0.); 
            bkgpsi2sM_fracC1.setConstant(kTRUE);
//            f = new RooExtendPdf("f","f",f_bkgPeakM1,nbkgPeak); // double Gaussian peak around B+
            break;
        default:
            break;
    }
        // Angular distribution of peaking background
    RooRealVar bkgjpsiPeak_c1("bkgjpsiPeak_c1","bkgjpsiPeak_c1",readParam(3,"bkgjpsiPeak_c1",0));
    RooRealVar bkgjpsiPeak_c2("bkgjpsiPeak_c2","bkgjpsiPeak_c2",readParam(3,"bkgjpsiPeak_c2",0));
    RooRealVar bkgjpsiPeak_c3("bkgjpsiPeak_c3","bkgjpsiPeak_c3",readParam(3,"bkgjpsiPeak_c3",0));
    RooRealVar bkgjpsiPeak_c4("bkgjpsiPeak_c4","bkgjpsiPeak_c4",readParam(3,"bkgjpsiPeak_c4",0));
    RooRealVar bkgjpsiPeak_c5("bkgjpsiPeak_c5","bkgjpsiPeak_c5",readParam(3,"bkgjpsiPeak_c5",0));
    RooRealVar bkgjpsiPeak_c6("bkgjpsiPeak_c6","bkgjpsiPeak_c6",readParam(3,"bkgjpsiPeak_c6",0));
    RooRealVar bkgjpsiPeak_c7("bkgjpsiPeak_c7","bkgjpsiPeak_c7",readParam(3,"bkgjpsiPeak_c7",0));
    RooRealVar bkgjpsiPeak_c8("bkgjpsiPeak_c8","bkgjpsiPeak_c8",readParam(3,"bkgjpsiPeak_c8",0));
    bkgjpsiPeak_c1.setError(readParam(3,"bkgjpsiPeak_c1",1));
    bkgjpsiPeak_c2.setError(readParam(3,"bkgjpsiPeak_c2",1));
    bkgjpsiPeak_c3.setError(readParam(3,"bkgjpsiPeak_c3",1));
    bkgjpsiPeak_c4.setError(readParam(3,"bkgjpsiPeak_c4",1));
    bkgjpsiPeak_c5.setError(readParam(3,"bkgjpsiPeak_c5",1));
    bkgjpsiPeak_c6.setError(readParam(3,"bkgjpsiPeak_c6",1));
    bkgjpsiPeak_c7.setError(readParam(3,"bkgjpsiPeak_c7",1));
    bkgjpsiPeak_c8.setError(readParam(3,"bkgjpsiPeak_c8",1));
    RooRealVar bkgpsi2sPeak_c1("bkgpsi2sPeak_c1","bkgpsi2sPeak_c1",readParam(5,"bkgpsi2sPeak_c1",0));
    RooRealVar bkgpsi2sPeak_c2("bkgpsi2sPeak_c2","bkgpsi2sPeak_c2",readParam(5,"bkgpsi2sPeak_c2",0));
    RooRealVar bkgpsi2sPeak_c3("bkgpsi2sPeak_c3","bkgpsi2sPeak_c3",readParam(5,"bkgpsi2sPeak_c3",0));
    RooRealVar bkgpsi2sPeak_c4("bkgpsi2sPeak_c4","bkgpsi2sPeak_c4",readParam(5,"bkgpsi2sPeak_c4",0));
    RooRealVar bkgpsi2sPeak_c5("bkgpsi2sPeak_c5","bkgpsi2sPeak_c5",readParam(5,"bkgpsi2sPeak_c5",0));
    RooRealVar bkgpsi2sPeak_c6("bkgpsi2sPeak_c6","bkgpsi2sPeak_c6",readParam(5,"bkgpsi2sPeak_c6",0));
    RooRealVar bkgpsi2sPeak_c7("bkgpsi2sPeak_c7","bkgpsi2sPeak_c7",readParam(5,"bkgpsi2sPeak_c7",0));
    RooRealVar bkgpsi2sPeak_c8("bkgpsi2sPeak_c8","bkgpsi2sPeak_c8",readParam(5,"bkgpsi2sPeak_c8",0));
    bkgpsi2sPeak_c1.setError(readParam(5,"bkgjpsiPeak_c1",1));
    bkgpsi2sPeak_c2.setError(readParam(5,"bkgjpsiPeak_c2",1));
    bkgpsi2sPeak_c3.setError(readParam(5,"bkgjpsiPeak_c3",1));
    bkgpsi2sPeak_c4.setError(readParam(5,"bkgjpsiPeak_c4",1));
    bkgpsi2sPeak_c5.setError(readParam(5,"bkgjpsiPeak_c5",1));
    bkgpsi2sPeak_c6.setError(readParam(5,"bkgjpsiPeak_c6",1));
    bkgpsi2sPeak_c7.setError(readParam(5,"bkgjpsiPeak_c7",1));
    bkgpsi2sPeak_c8.setError(readParam(5,"bkgjpsiPeak_c8",1));
    switch (iBin) {// Should be fixed constants already.
        case 3:
        case 5:
            bkgjpsiPeak_c1.setConstant(kTRUE);
            bkgjpsiPeak_c2.setConstant(kTRUE);
            bkgjpsiPeak_c3.setConstant(kTRUE);
            bkgjpsiPeak_c4.setConstant(kTRUE);
            bkgjpsiPeak_c5.setConstant(kTRUE);
            bkgjpsiPeak_c6.setConstant(kTRUE);
            bkgjpsiPeak_c7.setConstant(kTRUE);
            bkgjpsiPeak_c8.setConstant(kTRUE);
            bkgpsi2sPeak_c1.setConstant(kTRUE);
            bkgpsi2sPeak_c2.setConstant(kTRUE);
            bkgpsi2sPeak_c3.setConstant(kTRUE);
            bkgpsi2sPeak_c4.setConstant(kTRUE);
            bkgpsi2sPeak_c5.setConstant(kTRUE);
            bkgpsi2sPeak_c6.setConstant(kTRUE);
            bkgpsi2sPeak_c7.setConstant(kTRUE);
            bkgpsi2sPeak_c8.setConstant(kTRUE);
            break;
        default:
            bkgjpsiPeak_c1.setConstant(kTRUE);
            bkgjpsiPeak_c2.setConstant(kTRUE);
            bkgjpsiPeak_c3.setConstant(kTRUE);
            bkgjpsiPeak_c4.setConstant(kTRUE);
            bkgjpsiPeak_c5.setConstant(kTRUE);
            bkgjpsiPeak_c6.setConstant(kTRUE);
            bkgjpsiPeak_c7.setConstant(kTRUE);
            bkgjpsiPeak_c8.setConstant(kTRUE);
            bkgpsi2sPeak_c1.setConstant(kTRUE);
            bkgpsi2sPeak_c2.setConstant(kTRUE);
            bkgpsi2sPeak_c3.setConstant(kTRUE);
            bkgpsi2sPeak_c4.setConstant(kTRUE);
            bkgpsi2sPeak_c5.setConstant(kTRUE);
            bkgpsi2sPeak_c6.setConstant(kTRUE);
            bkgpsi2sPeak_c7.setConstant(kTRUE);
            bkgpsi2sPeak_c8.setConstant(kTRUE);
            break;
    }
    TString f_bkgjpsiPeakA_format = "(1+\
                                bkgjpsiPeak_c1*CosThetaK+\
                                bkgjpsiPeak_c2*CosThetaK**2+\
                                bkgjpsiPeak_c3*CosThetaL**2+\
                                bkgjpsiPeak_c4*CosThetaL**4+\
                                bkgjpsiPeak_c5*CosThetaL**2*CosThetaK+\
                                bkgjpsiPeak_c6*CosThetaL**4*CosThetaK+\
                                bkgjpsiPeak_c7*CosThetaL**2*CosThetaK**2+\
                                bkgjpsiPeak_c8*CosThetaL**4*CosThetaK**2)*\
                                exp(-(CosThetaK+0.5*CosThetaL)**4-(CosThetaK-0.5*CosThetaL)**4+1.5*(CosThetaK+0.5*CosThetaL)**2+1.5*(CosThetaK-0.5*CosThetaL)**2)";
    TString f_bkgpsi2sPeakA_format = "(1+\
                                bkgpsi2sPeak_c1*CosThetaK+\
                                bkgpsi2sPeak_c2*CosThetaK**2+\
                                bkgpsi2sPeak_c3*CosThetaL**2+\
                                bkgpsi2sPeak_c4*CosThetaL**4+\
                                bkgpsi2sPeak_c5*CosThetaL**2*CosThetaK+\
                                bkgpsi2sPeak_c6*CosThetaL**4*CosThetaK+\
                                bkgpsi2sPeak_c7*CosThetaL**2*CosThetaK**2+\
                                bkgpsi2sPeak_c8*CosThetaL**4*CosThetaK**2)*\
                                exp(-(CosThetaK+0.5*CosThetaL)**4-(CosThetaK-0.5*CosThetaL)**4+1.5*(CosThetaK+0.5*CosThetaL)**2+1.5*(CosThetaK-0.5*CosThetaL)**2)";

    RooArgSet f_bkgjpsiPeakA_argset;
    f_bkgjpsiPeakA_argset.add(RooArgSet(CosThetaL,CosThetaK));
    f_bkgjpsiPeakA_argset.add(RooArgSet(bkgjpsiPeak_c1,bkgjpsiPeak_c2,bkgjpsiPeak_c3,bkgjpsiPeak_c4,bkgjpsiPeak_c5,bkgjpsiPeak_c6,bkgjpsiPeak_c7,bkgjpsiPeak_c8));
    RooGenericPdf f_bkgjpsiPeakA("f_bkgjpsiPeakA", "f_bkgjpsiPeakA",f_bkgjpsiPeakA_format,f_bkgjpsiPeakA_argset);

    RooArgSet f_bkgpsi2sPeakA_argset;
    f_bkgpsi2sPeakA_argset.add(RooArgSet(CosThetaL,CosThetaK));
    f_bkgpsi2sPeakA_argset.add(RooArgSet(bkgpsi2sPeak_c1,bkgpsi2sPeak_c2,bkgpsi2sPeak_c3,bkgpsi2sPeak_c4,bkgpsi2sPeak_c5,bkgpsi2sPeak_c6,bkgpsi2sPeak_c7,bkgpsi2sPeak_c8));
    RooGenericPdf f_bkgpsi2sPeakA("f_bkgpsi2sPeakA", "f_bkgpsi2sPeakA",f_bkgpsi2sPeakA_format,f_bkgpsi2sPeakA_argset);

        // merge mass with angular term
    RooProdPdf f_bkgjpsiPeak("f_bkgjpsiPeak", "f_bkgjpsiPeak",RooArgSet(f_bkgjpsiPeakA,f_bkgjpsiPeakM12));
    RooProdPdf f_bkgpsi2sPeak("f_bkgpsi2sPeak", "f_bkgpsi2sPeak",RooArgSet(f_bkgpsi2sPeakA,f_bkgpsi2sPeakM12));
    printf("INFO: f_bkgPeak prepared.\n");

    // Observed spectrum = model*fullEfficiency
    RooRealVar nsig("nsig","nsig",50,0,5E3);
    RooRealVar nbkgComb("nbkgComb","nbkgComb",100,0,1E8);
    RooRealVar nbkgjpsiPeak("nbkgjpsiPeak","nbkgjpsiPeak",readParam(iBin,"nbkgjpsiPeak",0)*datasetLumi[0]/datasetLumi[2],0,1E8);
    RooRealVar nbkgpsi2sPeak("nbkgpsi2sPeak","nbkgpsi2sPeak",readParam(iBin,"nbkgpsi2sPeak",0)*datasetLumi[0]/datasetLumi[3],0,1E8);
    switch (iBin){
        case 2:
        case 3:
            nbkgpsi2sPeak.setVal(0.);
            nbkgpsi2sPeak.setConstant(kTRUE);
            break;
        case 4:
            break;
        case 5:
        case 6:
            nbkgjpsiPeak.setVal(0.);
            nbkgjpsiPeak.setConstant(kTRUE);
            break;
        case 0:
        case 1:
        case 7:
        default:
            nbkgjpsiPeak.setVal(0.);
            nbkgpsi2sPeak.setVal(0.);
            nbkgjpsiPeak.setConstant(kTRUE);
            nbkgpsi2sPeak.setConstant(kTRUE);
            break;
    }
    RooAddPdf f("kernel","kernel",RooArgList(f_bkgComb,f_bkgjpsiPeak,f_bkgpsi2sPeak,f_sig),RooArgList(nbkgComb,nbkgjpsiPeak,nbkgpsi2sPeak,nsig));// no penalty term
//    RooAddPdf f("kernel","kernel",RooArgList(f_bkgCombM,f_bkgjpsiPeakM12,f_bkgpsi2sPeakM12,f_sigM),RooArgList(nbkgComb,nbkgjpsiPeak,nbkgpsi2sPeak,nsig));// no angular part
//    RooAddPdf f("kernel","kernel",RooArgList(f_bkgComb,f_bkgjpsiPeakM12,f_bkgpsi2sPeakM12,f_sig),RooArgList(nbkgComb,nbkgjpsiPeak,nbkgpsi2sPeak,nsig));// no angular part


    // Extra penalty term to confine As, Fs, Fl, Afb.
    //RooRealVar t_penalty("t_penalty","t",0.01);
    //RooGenericPdf f_penaltyAfb("f_penaltyAfb","(1-TMath::Erf((afb-0.75*(1-fl))/(1.5*t_penalty*(1-fl))))*(1-TMath::Erf((-afb-0.75*(1-fl))/(1.5*t_penalty*(1-fl))))",RooArgSet(afb,fl,t_penalty));
    //RooGenericPdf f_penaltyAs("f_penaltyAfb","(1-TMath::Erf((afb-2*(1-fl)/3)/(1.5*t_penalty*(1-fl))))*(1-TMath::Erf((-afb-0.75*(1-fl))/(1.5*t_penalty*(1-fl))))",RooArgSet(afb,fl,t_penalty));
    //RooProdPdf f_penalty("f_penalty","f_penalty",f_penaltyAfb,f_penaltyAs);
    //RooProdPdf f("f","f",f_model,f_penalty);
    printf("INFO: f_penalty NOT applied.\n");

    // Gaussian constraints
    RooGaussian gaus_sigGauss1_sigma("gaus_sigGauss1_sigma","gaus_sigGauss1_sigma",sigGauss1_sigma,RooConst(readParam(iBin,"sigGauss1_sigma",0)),RooConst(readParam(iBin,"sigGauss1_sigma",1)));
    RooGaussian gaus_sigGauss2_sigma("gaus_sigGauss2_sigma","gaus_sigGauss2_sigma",sigGauss2_sigma,RooConst(readParam(iBin,"sigGauss2_sigma",0)),RooConst(readParam(iBin,"sigGauss2_sigma",1)));
    RooGaussian gaus_sigM_frac("gaus_sigM_frac","gaus_sigM_frac",sigM_frac,RooConst(readParam(iBin,"sigM_frac",0)),RooConst(readParam(iBin,"sigM_frac",1)));
    RooGaussian gaus_nbkgjpsiPeak("gaus_nbkgjpsiPeak","gaus_nbkgjpsiPeak",nbkgjpsiPeak,RooConst(readParam(iBin,"nbkgjpsiPeak",0)*datasetLumi[0]/datasetLumi[2]),RooConst(readParam(iBin,"nbkgjpsiPeak",1)*datasetLumi[0]/datasetLumi[2]));
    RooGaussian gaus_nbkgpsi2sPeak("gaus_nbkgpsi2sPeak","gaus_nbkgpsi2sPeak",nbkgpsi2sPeak,RooConst(readParam(iBin,"nbkgpsi2sPeak",0)*datasetLumi[0]/datasetLumi[3]),RooConst(readParam(iBin,"nbkgpsi2sPeak",1)*datasetLumi[0]/datasetLumi[3]));
    RooBifurGauss gaus_fs("gaus_fs","gaus_fs",fs,RooConst(readParam(iBin,"fs",0,0.00129254)),RooConst(readParam(iBin,"fs",1)),RooConst(readParam(iBin,"fs",1,0.0101371)));
    RooBifurGauss gaus_as("gaus_as","gaus_as",as,RooConst(readParam(iBin,"as",0,-0.0975917)),RooConst(readParam(iBin,"as",1)),RooConst(readParam(iBin,"as",1,0.0049092)));
//    RooBifurGauss gaus_fs("gaus_fs","gaus_fs",fs,RooConst(0.0129254),RooConst(0.00898344),RooConst(0.0101371));// 2011 result
//    RooBifurGauss gaus_as("gaus_as","gaus_as",as,RooConst(-0.0975919),RooConst(0.00490805),RooConst(0.0049092));
    
    RooArgSet gausConstraints(gaus_sigGauss1_sigma,gaus_sigGauss2_sigma,gaus_sigM_frac);
    if (iBin<2 || iBin > 6) gausConstraints.add(RooArgSet(gaus_nbkgjpsiPeak,gaus_nbkgpsi2sPeak));
    if (iBin == 3 || iBin == 5) gausConstraints.add(RooArgSet(gaus_fs,gaus_as));
    printf("INFO: gausConstraints are settled.\n");
    
    // Get data and apply unbinned fit
    int mumuMassWindowBin = 1+4*isCDFcut;
    if (iBin==3 || iBin==5) mumuMassWindowBin = 0; // no cut
    RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(Q2, CosThetaK, CosThetaL, Bmass, Mumumass, Mumumasserr),TString::Format("(%s) && (%s)",q2range[iBin],mumuMassWindow[mumuMassWindowBin]),0);
    f.fitTo(*data,Extended(kTRUE),ExternalConstraints(gausConstraints),Minimizer("Minuit"),Minos(kTRUE));

    // Draw the frame on the canvas
    TCanvas* c = new TCanvas("c");
    RooPlot* framemass = Bmass.frame();
    data->plotOn(framemass,Binning(20));// TGraphPainter options are allowed by DrawOption()
    f.plotOn(framemass,LineColor(1));
    f.plotOn(framemass,Components(f_sig),LineColor(4),LineWidth(2));
    f.plotOn(framemass,Components(f_bkgComb),LineColor(2),LineWidth(2),LineStyle(2));
    f.plotOn(framemass,Components(f_bkgjpsiPeak),LineColor(6),LineWidth(2),LineStyle(2));
    f.plotOn(framemass,Components(f_bkgpsi2sPeak),LineColor(8),LineWidth(2),LineStyle(2));

    framemass->SetTitle("");
    framemass->SetMinimum(0);
    framemass->Draw();
    
    TLatex *t1 = new TLatex();
    t1->SetNDC();
    double fixNDC = 0;
    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",q2range[iBin]));
    t1->DrawLatex(.65,.79+fixNDC,TString::Format("Y_{Signal}=%.3f",nsig.getVal()));
    t1->DrawLatex(.65,.72+fixNDC,TString::Format("Y_{J/#Psi}=%.3f",nbkgjpsiPeak.getVal()));
    t1->DrawLatex(.65,.65+fixNDC,TString::Format("Y_{#Psi'}=%.3f",nbkgpsi2sPeak.getVal()));
    t1->DrawLatex(.65,.58+fixNDC,TString::Format("Y_{Comb}=%.3f",nbkgComb.getVal()));
    c->Print(TString::Format("./plots/%s_bin%d.pdf",outfile,iBin));

    // Draw projection to CosThetaK
    RooPlot* framecosk = CosThetaK.frame(); 
    data->plotOn(framecosk,Binning(20)); 
    f.plotOn(framecosk,LineColor(1)); 
    f.plotOn(framecosk,Components(f_sig),LineColor(4),LineWidth(2));
    f.plotOn(framecosk,Components(f_bkgComb),LineColor(2),LineWidth(2),LineStyle(2));
    f.plotOn(framecosk,Components(f_bkgjpsiPeak),LineColor(6),LineWidth(2),LineStyle(2));
    f.plotOn(framecosk,Components(f_bkgpsi2sPeak),LineColor(8),LineWidth(2),LineStyle(2));

    framecosk->SetTitle("");
    framecosk->SetMinimum(0);
    framecosk->Draw();

    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",q2range[iBin]));
    t1->DrawLatex(.35,.79+fixNDC,TString::Format("F_{L}=%.3f#pm%.3f",fl.getVal(),fl.getError()));
    c->Update();
    c->Print(TString::Format("./plots/%s_cosk_bin%d.pdf",outfile,iBin));

    // Draw projection to CosThetaL
    RooPlot* framecosl = CosThetaL.frame(); 
    data->plotOn(framecosl,Binning(20)); 
    f.plotOn(framecosl,LineColor(1)); 
    f.plotOn(framecosl,Components(f_sig),LineColor(4),LineWidth(2));
    f.plotOn(framecosl,Components(f_bkgComb),LineColor(2),LineWidth(2),LineStyle(2));
    f.plotOn(framecosl,Components(f_bkgjpsiPeak),LineColor(6),LineWidth(2),LineStyle(2));
    f.plotOn(framecosl,Components(f_bkgpsi2sPeak),LineColor(8),LineWidth(2),LineStyle(2));

    framecosl->SetTitle("");
    framecosl->SetMinimum(0);
    framecosl->Draw();

    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",q2range[iBin]));
    t1->DrawLatex(.35,.79+fixNDC,TString::Format("A_{FB}=%.3f#pm%.3f",afb.getVal(),afb.getError()));
    c->Update();
    c->Print(TString::Format("./plots/%s_cosl_bin%d.pdf",outfile,iBin));

    // clear
    delete t1;
    delete c;
    delete data;

    // write output
    double val[3]={0,0,0};
    val[0] = fl.getVal();val[1] = fl.getError();
    writeParam(iBin, "fl", val);
    val[0] = afb.getVal();val[1] = afb.getError();
    writeParam(iBin, "afb",val);
    val[0] = fs.getVal();val[1] = fs.getErrorHi();val[2] = fs.getErrorLo();
    writeParam(iBin, "fs", val, 3);
    val[0] = as.getVal();val[1] = as.getErrorHi();val[2] = fs.getErrorLo();
    writeParam(iBin, "as", val, 3);
}//}}}

void angular3D_bin_simplified(int iBin, const char outfile[]="testAngular3D")
{//{{{
    // Remark: You must use RooFit!! It's better in unbinned fit.
    //         Extended ML fit is adopted by Mauro, just follow!!
    //         This function is NOT designed to determine parameters in resonance region.
//    if (iBin==3 || iBin==5) return;
    
    // Read data
    RooRealVar CosThetaK("CosThetaK"     , "cos#theta_{K}"       , -1. , 1.   ) ;
    RooRealVar CosThetaL("CosThetaL"     , "cos#theta_{L}"       , -1. , 1.   ) ;
    RooRealVar Bmass("Bmass"             , "M_{K^{*}#mu#mu}"     , 5.  , 5.56 ) ;
    RooAddition Bmass_offset("Bmass_offset","Bmass_offset",RooArgSet(Bmass,RooConst(-5)));
    RooRealVar Mumumass("Mumumass"       , "M^{#mu#mu}"          , 0.  , 10.  ) ;
    RooRealVar Mumumasserr("Mumumasserr" , "Error of M^{#mu#mu}" , 0.  , 10.  ) ;
    RooRealVar Q2("Q2"                   , "q^{2}"               , 0.5 , 20.  ) ;
        
    // Create parameters and PDFs
    // Signal gaussian
    RooRealVar sigGauss_mean("sigGauss_mean","M_{K*#Mu#Mu}",5.274,5.25,5.30);
    RooRealVar sigGauss_sigma("sigGauss_sigma","#sigma",0.003,0.002,0.05);
    sigGauss_sigma.setError(0.0004);
    TString sigM_format = "exp(-0.5*((Bmass-sigGauss_mean)/sigGauss_sigma)**2)/2.50662827463100024/sigGauss_sigma";
    RooArgSet sigM_argset(Bmass,sigGauss_mean,sigGauss_sigma);
    RooGaussian f_sigM("f_sigM","f_sigM", Bmass, sigGauss_mean, sigGauss_sigma);
    RooGenericPdf f_sigM2("f_sigM2","f_sigM2",sigM_format,sigM_argset);

    RooRealVar sigA_c1("sigA_c1","",0.-1.,1);
    RooRealVar sigA_c2("sigA_c2","",0.-1.,1);
    TString sigK_format = "1+sigA_c1*CosThetaK";
    TString sigL_format = "1+sigA_c2*CosThetaL";
    TString sigA_format = TString::Format("(%s)*(%s)",sigK_format.Data(),sigL_format.Data());
    RooArgSet sigK_argset(CosThetaK,sigA_c1);
    RooArgSet sigL_argset(CosThetaL,sigA_c2);
    RooArgSet sigA_argset;
    sigA_argset.add(sigL_argset);
    sigA_argset.add(sigK_argset);
    RooGenericPdf f_sigL("f_sigL","f_sigL",sigL_format,sigL_argset);
    RooGenericPdf f_sigK("f_sigK","f_sigK",sigK_format,sigK_argset);
    RooGenericPdf f_sigA("f_sigA","f_sigA",sigA_format,sigA_argset);

    TString sigMK_format = TString::Format("(%s)*(%s)",sigM_format.Data(),sigK_format.Data());
    TString sigML_format = TString::Format("(%s)*(%s)",sigM_format.Data(),sigL_format.Data());
    TString sig_format = TString::Format("(%s)*(%s)",sigM_format.Data(),sigA_format.Data());
    RooArgSet sigMK_argset;
    sigMK_argset.add(sigM_argset);
    sigMK_argset.add(sigK_argset);
    RooArgSet sigML_argset;
    sigML_argset.add(sigM_argset);
    sigML_argset.add(sigL_argset);
    RooArgSet sig_argset;
    sig_argset.add(sigM_argset);
    sig_argset.add(sigA_argset);
    RooGenericPdf f_sigMK("f_sigMK","f_sigMK",sigMK_format,sigMK_argset);
    RooGenericPdf f_sigML("f_sigML","f_sigML",sigML_format,sigML_argset);

    RooProdPdf f_sig("f_sig","f_sig",RooArgSet(f_sigM2,f_sigA));
    RooGenericPdf f_sig2("f_sig2","f_sig2",sig_format,sig_argset);
    
        // background
//    RooRealVar bkgCombM_c1("bkgCombM_c1","",-2.5,-10,10);
    RooRealVar bkgCombM_c2("bkgCombM_c2","",2.6,-10,10);
//    TString bkgCombM_format = "(1+bkgCombM_c1*(Bmass-5)+bkgCombM_c2*(Bmass-5)**2)/(0.56+bkgCombM_c1*0.3136/2+bkgCombM_c2*0.175616/3)";
//    RooArgSet bkgCombM_argset(Bmass,bkgCombM_c1,bkgCombM_c2);
    RooRealVar bkgCombM_c1("bkgCombM_c1","",0.,-1.7,10);
    TString bkgCombM_format = "1+bkgCombM_c1*(Bmass-5)";
    RooArgSet bkgCombM_argset(Bmass,bkgCombM_c1);
    RooPolynomial f_bkgCombM("f_bkgCombM","f_bkgCombM",Bmass_offset,RooArgSet(bkgCombM_c1,bkgCombM_c2));
    RooGenericPdf f_bkgCombM2("f_bkgCombM2","f_bkgCombM2",bkgCombM_format,bkgCombM_argset);

    RooRealVar bkgCombA_c1("bkgCombA_c1","",-0.5,-1,10);
    RooRealVar bkgCombA_c2("bkgCombA_c2","",4,0.,20);
    RooRealVar bkgCombA_c3("bkgCombA_c3","", -0.006,  -1,1);
    RooRealVar bkgCombA_c4("bkgCombA_c4","", .5,1E-2,5);
    TString bkgCombK_format = "(1+bkgCombA_c1*CosThetaK+bkgCombA_c2*CosThetaK**2)";
    TString bkgCombL_format = "exp(-0.5*((CosThetaL-bkgCombA_c3)/bkgCombA_c4)**2)";
    TString bkgCombA_format = TString::Format("(%s)*(%s)",bkgCombK_format.Data(),bkgCombL_format.Data());
    TString bkgCombML_format = TString::Format("(%s)*(%s)",bkgCombM_format.Data(),bkgCombL_format.Data());
    TString bkgCombMK_format = TString::Format("(%s)*(%s)",bkgCombM_format.Data(),bkgCombK_format.Data());
    RooArgSet bkgCombK_argset(CosThetaK,bkgCombA_c1,bkgCombA_c2);
    RooArgSet bkgCombL_argset(CosThetaL,bkgCombA_c3,bkgCombA_c4);
    RooArgSet bkgCombA_argset;
    bkgCombA_argset.add(bkgCombK_argset);
    bkgCombA_argset.add(bkgCombL_argset);
    RooArgSet bkgCombML_argset;
    bkgCombML_argset.add(bkgCombM_argset);
    bkgCombML_argset.add(bkgCombL_argset);
    RooArgSet bkgCombMK_argset;
    bkgCombMK_argset.add(bkgCombM_argset);
    bkgCombMK_argset.add(bkgCombK_argset);
    RooGenericPdf f_bkgCombA("f_bkgCombA","f_bkgCombA",bkgCombA_format,bkgCombA_argset);
    RooGenericPdf f_bkgCombML("f_bkgCombML","f_bkgCombML",bkgCombML_format,bkgCombML_argset);
    RooGenericPdf f_bkgCombMK("f_bkgCombMK","f_bkgCombMK",bkgCombMK_format,bkgCombMK_argset);
    
    // Observed spectrum = model*fullEfficiency
//    RooRealVar nsig("nsig","nsig",7000,0,2E4);
//    RooRealVar nbkg("nbkg","nbkg",10000,1,2E4);
    RooRealVar nsig("nsig","nsig",10,0,2E4);
    RooRealVar nbkg("nbkg","nbkg",1000,1,2E4);
    TString f_bkg_format = TString::Format("(%s)*(%s)", bkgCombM_format.Data(), bkgCombA_format.Data());
    printf("f_bkg_format is '%s'\n",f_bkg_format.Data());
    RooProdPdf f_bkg("f_bkg","f_bkg",RooArgSet(f_bkgCombM,f_bkgCombA));
    
    RooArgList bkg_argset;
    bkg_argset.add(bkgCombM_argset);
    bkg_argset.add(bkgCombA_argset);
    RooGenericPdf f_bkg2("f_bkg2","f_bkg2",f_bkg_format,bkg_argset);

    RooAddPdf f("f","f",RooArgList(f_sig,f_bkg2),RooArgList(nsig,nbkg));
//    RooAddPdf f("f","f",RooArgList(f_sig2,f_bkg2),RooArgList(nsig,nbkg));
    
    // Get data and apply unbinned fit
    int mumuMassWindowBin = 1+4*isCDFcut;
    if (iBin==3 || iBin==5) mumuMassWindowBin = 0; // no cut
    RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(Q2, Bmass, CosThetaK, CosThetaL, Mumumass, Mumumasserr),TString::Format("(%s) && (%s)",q2range[iBin],mumuMassWindow[mumuMassWindowBin]),0);
    
//    f_bkgCombM.fitTo(*data,Minimizer("Minuit2"));
    f.fitTo(*data,Extended(kTRUE),Minimizer("Minuit"),Minos(kTRUE));

    // Draw the frame on the canvas
    TCanvas* c = new TCanvas("c");
    RooPlot* framemass = Bmass.frame();
    data->plotOn(framemass,Binning(20));// TGraphPainter options are allowed by DrawOption()
    f.plotOn(framemass,LineColor(1));
    f.plotOn(framemass,Components(f_sig),LineColor(4),LineWidth(2));
    f.plotOn(framemass,Components(f_bkg),LineColor(3),LineWidth(2));
    f.plotOn(framemass,Components(f_sig2),LineColor(6),LineWidth(2),LineStyle(7));
    f.plotOn(framemass,Components(f_bkg2),LineColor(5),LineWidth(2),LineStyle(7));

    framemass->SetTitle("");
    framemass->SetMinimum(0);
    framemass->Draw();
    
    TLatex *t1 = new TLatex();
    t1->SetNDC();
    double fixNDC = 0;
    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",q2range[iBin]));
    t1->DrawLatex(.65,.79+fixNDC,TString::Format("Y_{Sig}=%.3f",nsig.getVal()));
    t1->DrawLatex(.65,.72+fixNDC,TString::Format("Y_{Bkg}=%.3f",nbkg.getVal()));
    c->Print(TString::Format("./plots/test/%s_bin%d.pdf",outfile,iBin));
    
    // CosThetaK
    RooPlot* frameK = CosThetaK.frame();
    data->plotOn(frameK,Binning(20));// TGraphPainter options are allowed by DrawOption()
    f.plotOn(frameK,LineColor(1));
    f.plotOn(frameK,Components(f_sig),LineColor(4),LineWidth(2));
    f.plotOn(frameK,Components(f_bkg),LineColor(3),LineWidth(2));
    f.plotOn(frameK,Components(f_sig2),LineColor(6),LineWidth(2),LineStyle(7));
    f.plotOn(frameK,Components(f_bkg2),LineColor(5),LineWidth(2),LineStyle(7));

    frameK->SetTitle("");
    frameK->SetMinimum(0);
    frameK->Draw();
    
    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",q2range[iBin]));
    c->Update();
    c->Print(TString::Format("./plots/test/%s_cosk_bin%d.pdf",outfile,iBin));
    
    // CosThetaL
    RooPlot* frameL = CosThetaL.frame();
    data->plotOn(frameL,Binning(20));// TGraphPainter options are allowed by DrawOption()
    f.plotOn(frameL,LineColor(1));
    f.plotOn(frameL,Components(f_sig),LineColor(4),LineWidth(2));
    f.plotOn(frameL,Components(f_bkg),LineColor(3),LineWidth(2));
    f.plotOn(frameL,Components(f_sig2),LineColor(6),LineWidth(2),LineStyle(7));
    f.plotOn(frameL,Components(f_bkg2),LineColor(5),LineWidth(2),LineStyle(7));

    frameL->SetTitle("");
    frameL->SetMinimum(0);
    frameL->Draw();
    
    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",q2range[iBin]));
    c->Update();
    c->Print(TString::Format("./plots/test/%s_cosl_bin%d.pdf",outfile,iBin));

    // clear
    delete t1;
    delete c;
    delete data;

}//}}}

//_________________________________________________________________________________

void printListOfTChainElements(TChain *chain){
    TObjArray *fileElements=chain->GetListOfFiles();
    int nFiles = fileElements->GetEntries();
    TIter next(fileElements);
    TChainElement *chEl=0;
    for( int entry=0; entry < nFiles; entry++ ) {
        chEl=(TChainElement*)next();
        printf("%s\n",chEl->GetTitle());
    }
}

//_________________________________________________________________________________
int main(int argc, char** argv) {
    // Tags
    is7TeVCheck = false;   
    isCDFcut = 0;// 1 for true, 0 for false.
    
    // Help message
    if (argc <= 2) {
        printf("Usage       : ./test Function infile\n");
        printf("Functions   :\n");
        printf("    test                Target test function\n");
        printf("    bmass               Fit to mass spectrum using a double Gaussian signal and Chebyshev bkg.\n");
        printf("Remark      :\n");
        printf("    1. Outputs will be stored in ./plots, please keep the directory.\n");
        printf("    2. Wildcard is allowed for infile. But you must quote infile like \"inputData_Run*.root\"!\n");
        return 0;
    }

    // main
    TString func    = argv[1];
    TString infile  = argv[2];
  
    if (func == "bmass") {
        if (argc != 4){
            printf("./fit bmass infile binID\n");
            for (int i = 0; i < 10; i++) {
                printf("    Bin %d : %s\n",i,q2range[i]);
            }
            return 0;
        }
        int iBin = atoi(argv[3]);
        ch->Add(infile.Data());
        if (ch == NULL) return 1;
        const char outfile[]="bmass";
        bmass(iBin, outfile); 
    }else if (func == "test"){
        ch->Add(infile.Data());
        printListOfTChainElements(ch);
        if (ch == NULL) return 1;
        const char outfile[]="test";
        for (int iBin = 0; iBin < 8; iBin++) {
//            if (iBin != 3 && iBin != 5) continue;
//            if (iBin == 3 || iBin == 5) continue;
        }
//        angular3D_bin(2);
//        angular3D_bin_simplified(3);
        angular3D_bin_simplified(0);
        angular3D_bin_simplified(1);
    }else{ 
        cerr << "No function available for: " << func.Data() << endl; 
    }
    printf("%lld entries processed.\n",ch->GetEntries());

    return 0 ;
}
