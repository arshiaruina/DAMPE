// To compile, g++ -Wall `root-config --cflags --libs` -I$DMPSWSYS/include -L$DMPSWSYS/lib -lDmpEvent -lDmpService eta.cc -o ../submit/eta.exe

// example file
// /beegfs/users/ruina/VAequalisation/out/20181019/merged/merged_160919_104734.root

#include "TH1.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"

//C++
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <sstream>
#include <cmath>
#include <vector>
// ROOT
#include "TROOT.h"
#include "TFile.h"
#include "TMath.h"
#include "TF1.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TRandom.h"
#include "TStyle.h"

//#include "../inc/va_equalisation.h"
//#include "mylangaus.C"
#include "TPaveStats.h"
using namespace std;

Double_t langaufun(Double_t *x, Double_t *par) {
   //Fit parameters:
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //
   //In the Landau distribution (represented by the CERNLIB approximation),
   //the maximum is located at x=-0.22278298 with the location parameter=0.
   //This shift is corrected within this function, so that the actual
   //maximum is identical to the MP parameter.
      // Numeric constants
      Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
      Double_t mpshift  = -0.22278298;       // Landau maximum location
      // Control constants
      Double_t np = 100.0;      // number of convolution steps
      Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas
      // Variables
      Double_t xx;
      Double_t mpc;
      Double_t fland;
      Double_t sum = 0.0;
      Double_t xlow,xupp;
      Double_t step;
      Double_t i;
      // MP shift correction
      mpc = par[1] - mpshift * par[0];
      // Range of convolution integral
      xlow = x[0] - sc * par[3];
      xupp = x[0] + sc * par[3];
      step = (xupp-xlow) / np;
      // Convolution integral of Landau and Gaussian by sum
      for(i=1.0; i<=np/2; i++) {
         xx = xlow + (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);
         xx = xupp - (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);
      }
      return (par[2] * step * sum * invsq2pi / par[3]);
}


TF1 *langaufit(TH1D *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF)
{
   // Once again, here are the Landau * Gaussian parameters:
   //   par[0]=Width (scale) parameter of Landau density
   //   par[1]=Most Probable (MP, location) parameter of Landau density
   //   par[2]=Total area (integral -inf to inf, normalization constant)
   //   par[3]=Width (sigma) of convoluted Gaussian function
   //
   // Variables for langaufit call:
   //   his             histogram to fit
   //   fitrange[2]     lo and hi boundaries of fit range
   //   startvalues[4]  reasonable start values for the fit
   //   parlimitslo[4]  lower parameter limits
   //   parlimitshi[4]  upper parameter limits
   //   fitparams[4]    returns the final fit parameters
   //   fiterrors[4]    returns the final fit errors
   //   ChiSqr          returns the chi square
   //   NDF             returns ndf
   Int_t i;
   Char_t FunName[100];
   sprintf(FunName,"Fitfcn_%s",his->GetName());
   TF1 *ffitold = (TF1*)gROOT->GetListOfFunctions()->FindObject(FunName);
   if (ffitold) delete ffitold;
   TF1 *ffit = new TF1(FunName,langaufun,fitrange[0],fitrange[1],4);
   ffit->SetParameters(startvalues);
   ffit->SetParNames("Width","MP","Area","GSigma");
   for (i=0; i<4; i++) {
      ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);
   }
   his->Fit(FunName,"RB0");   // fit within specified range, use ParLimits, do not plot
   ffit->GetParameters(fitparams);    // obtain fit parameters
   for (i=0; i<4; i++) {
      fiterrors[i] = ffit->GetParError(i);     // obtain fit parameter errors
   }
   ChiSqr[0] = ffit->GetChisquare();  // obtain chi^2
   NDF[0] = ffit->GetNDF();           // obtain ndf
   return (ffit);              // return fit function
}
Int_t langaupro(Double_t *params, Double_t &maxx, Double_t &FWHM) {
   // Seaches for the location (x value) at the maximum of the
   // Landau-Gaussian convolute and its full width at half-maximum.
   //
   // The search is probably not very efficient, but it's a first try.
   Double_t p,x,fy,fxr,fxl;
   Double_t step;
   Double_t l,lold;
   Int_t i = 0;
   Int_t MAXCALLS = 10000;
   // Search for maximum
   p = params[1] - 0.1 * params[0];
   step = 0.05 * params[0];
   lold = -2.0;
   l    = -1.0;
   while ( (l != lold) && (i < MAXCALLS) ) {
      i++;
      lold = l;
      x = p + step;
      l = langaufun(&x,params);
      if (l < lold)
         step = -step/10;
      p += step;
   }
   if (i == MAXCALLS)
      return (-1);
   maxx = x;
   fy = l/2;
   // Search for right x location of fy
   p = maxx + params[0];
   step = params[0];
   lold = -2.0;
   l    = -1e300;
   i    = 0;
   while ( (l != lold) && (i < MAXCALLS) ) {
      i++;
      lold = l;
      x = p + step;
      l = TMath::Abs(langaufun(&x,params) - fy);
      if (l > lold)
         step = -step/10;
      p += step;
   }
   if (i == MAXCALLS)
      return (-2);
   fxr = x;
   // Search for left x location of fy
   p = maxx - 0.5 * params[0];
   step = -params[0];
   lold = -2.0;
   l    = -1e300;
   i    = 0;
   while ( (l != lold) && (i < MAXCALLS) ) {
      i++;
      lold = l;
      x = p + step;
      l = TMath::Abs(langaufun(&x,params) - fy);
      if (l > lold)
         step = -step/10;
      p += step;
   }
   if (i == MAXCALLS)
      return (-3);
   fxl = x;
   FWHM = fxr - fxl;
   return (0);
}

//int main(int argc, char** argv) {
int main() {


//    //std::ifstream histFile;
//    //histFile.open("histolist_new.txt");
//    std::string histName = "hVAEnergyX_154_2_0";
//    //std::string inFileName = argv[1];
//    std::string inFileName = "/beegfs/users/ruina/VAequalisation/out/20181019/merged/merged_160919_104734.root";
//    TFile *inFile = new TFile(inFileName.c_str());
//    //TH1D *hist0 = (TH1D*)inFile->Get(histName0.c_str());
//    TH1D *hSNR = (TH1D*)inFile->Get(histName.c_str());
//
//    TFile *outFile = new TFile("langaufit.root", "RECREATE");
//    TCanvas *c1 = new TCanvas("c1",histName.c_str(),800,600);
//
//
//    // Fill Histogram
//    //Int_t data[100] = {0,0,0,0,0,0,2,6,11,18,18,55,90,141,255,323,454,563,681,
//    //                 737,821,796,832,720,637,558,519,460,357,291,279,241,212,
//    //                 153,164,139,106,95,91,76,80,80,59,58,51,30,49,23,35,28,23,
//    //                 22,27,27,24,20,16,17,14,20,12,12,13,10,17,7,6,12,6,12,4,
//    //                 9,9,10,3,4,5,2,4,1,5,5,1,7,1,6,3,3,3,4,5,4,4,2,2,7,2,4};
//    //TH1F *hSNR = new TH1F("snr","Signal-to-noise",400,0,400);
//    //for (Int_t i=0; i<100; i++) hSNR->Fill(i,data[i]);
//    // Fitting SNR histo
//    printf("Fitting...\n");
//    // Setting fit range and start values
//    Double_t fr[2];
//    Double_t sv[4], pllo[4], plhi[4], fp[4], fpe[4];
//    //fr[0]=0.3*hSNR->GetMean();
//    //fr[1]=3.0*hSNR->GetMean();
//
//    //double fr[2], sv[4], pllo[4], plhi[4];
//    //fr[0] = 40.;      fr[1] = 150.;
//    //sv[0] = 2.;       sv[1] = 50.;       sv[2] = 2e5;      sv[3] = 3.;    
//    //pllo[0] = 0.;     pllo[1] = 10.;     pllo[2] = 1e4;    pllo[3] = 0.;  
//    //plhi[0] = 10.;    plhi[1] = 1e5;     plhi[2] = 1e6;    plhi[3] = 1e5; 
//
//
//    fr[0] = 0.;      fr[1] = 180.;
//    pllo[0]=0.5; pllo[1]=10.0; pllo[2]=1.0;         pllo[3]=0.4;
//    plhi[0]=5.0; plhi[1]=70.0; plhi[2]=1000000.0;   plhi[3]=5.0;
//    sv[0]=1.8;   sv[1]=50.0;   sv[2]=50000.0;       sv[3]=3.0;
//
//    //pllo[0]=0.5; pllo[1]=5.0; pllo[2]=1.0; pllo[3]=0.4;
//    //plhi[0]=5.0; plhi[1]=50.0; plhi[2]=1000000.0; plhi[3]=5.0;
//    //sv[0]=1.8; sv[1]=20.0; sv[2]=50000.0; sv[3]=3.0;
//    Double_t chisqr;
//    Int_t    ndf;
//    TF1 *fitsnr = langaufit(hSNR,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
//    Double_t SNRPeak, SNRFWHM;
//    langaupro(fp,SNRPeak,SNRFWHM);
//    printf("Fitting done\nPlotting results...\n");
//    // Global style settings
//    gStyle->SetOptStat(1111);
//    gStyle->SetOptFit(111);
//    gStyle->SetLabelSize(0.03,"x");
//    gStyle->SetLabelSize(0.03,"y");
//    hSNR->GetXaxis()->SetRange(0,70);
//    hSNR->Draw();
//    fitsnr->Draw("hist same");
//    c1->Print();
//    outFile->WriteTObject(c1);
//    outFile->Write();
//    outFile->Close();
//}


/* TODO:
 * Access the generated (merged) file.root which has the histograms of the energy distributions for all VAs
 * langaufit accessess each histogram and makes the fit
 * Fit results (fitted energy plots of all the VAs) stored in new root file
 * Fit results (MPVs of all VAs) stored in array (1152 VAs x 2 eta-regions) and plotted for the 2 eta regions
 * Mean of MPV distributions stored (2 values for the 2 eta regions)
 * Correction factor for each VA for each eta region = Mean MPV for that eta region / MPV of the VA for that eta region
 * */

	//TStopwatch sw;
	//sw.Start();

    //double MPV[1152][2] = {0.};
    TH1D *hMPV0 = new TH1D("hMPV0","MPV of all VAs in RO-strip region",100,40.,60.);
    TH1D *hMPV1 = new TH1D("hMPV1","MPV of all VAs in floating-strip region",100,30.,50.);


//-----------------------------------------------------------//
    std::ifstream histFile;
    histFile.open("histolist_new.txt");
    std::string histName;

    int count = 0;
    //std::string inFileName = argv[1];
    //std::string inFileName = "/beegfs/users/ruina/VAequalisation/out/20181019/merged/merged_160919_104734.root";
    std::string inFileName = "../out/20181001_20181009/merged/merged_231019_112429.root";
    TFile *inFile = new TFile(inFileName.c_str());
    // - TList *list = inFile->GetListOfKeys();
    //TFile *outFile = new TFile("../out/root_forum.root", "recreate");
    //TFile *outFile = new TFile("../out/20181019/langaufit.root", "RECREATE"); 
    //TFile *outFile = new TFile("../out/201810_wk1/langaufit.root", "RECREATE"); 
    TFile *outFile = new TFile("langaufit.root", "RECREATE"); 
    // - TKey *key;
    // - TIter iter(list); //or TIter iter(list->MakeIterator());
    // - static TString classname("TH1D");
    // - while((key = (TKey*)iter()) && count < 2) {
    // -     if (key->GetClassName() == classname) {
    // -         TH1D *hist = (TH1D*)key->ReadObj();
    // -         if (hist) {
    // -             const char* histName_c = hist->GetName();
    // -             std::string histName(histName_c);

    //while(std::getline(histFile,histName) && count < 5){ // for debug run
    while(std::getline(histFile,histName)){ // for analysis run
        
        //if(count%2 == 1) continue;

        //if(histName == "hEtaX" || histName == "hEtaY" || hist0->GetEntries() == 0 || hist1->GetEntries() == 0) continue;
        //if(!(histName == "hVAEnergyY_141_5" || histName == "hVAEnergyY_106_2")) continue;
        if(!(histName == "hVAEnergyX_154_2" || histName == "hVAEnergyY_33_5")) continue;
        //if(!(histName == "hVAEnergyX_154_2")) continue;
        string histName0 = histName + "_0";
        string histName1 = histName + "_1";
        //TH1D *hist = (TH1D*)inFile->Get(histName.c_str()); 
        TH1D *hist0 = (TH1D*)inFile->Get(histName0.c_str()); 
        TH1D *hist1 = (TH1D*)inFile->Get(histName1.c_str()); 

        /*------ Fitting start ------*/
   
        // --- Hist0 --- //             
                 
        std::cout << "hist0 std dev " << hist0->GetStdDev() << std::endl;
        std::cout << "hist0 mean " << hist0->GetMean() << std::endl;
        std::cout << "hist0 integral " << hist0->Integral() << std::endl;

        // Setting fit range and start values

        double fr0[2], sv0[4], pllo0[4], plhi0[4];
        //fr0[0] = 40.;      fr0[1] = 150.;
        //sv0[0] = 2.;       sv0[1] = 50.;       sv0[2] = 2e5;      sv0[3] = 3.;    
        //pllo0[0] = 0.;     pllo0[1] = 10.;     pllo0[2] = 1e4;    pllo0[3] = 0.;  
        //plhi0[0] = 10.;    plhi0[1] = 1e5;     plhi0[2] = 1e6;    plhi0[3] = 1e5; 
 
        fr0[0] = 40.;      fr0[1] = 150.;
        pllo0[0]=0.5; pllo0[1]=10.0; pllo0[2]=1.0;         pllo0[3]=0.4;
        plhi0[0]=9.0; plhi0[1]=70.0; plhi0[2]=1000000.0;   plhi0[3]=9.0;
        sv0[0]=1.8;   sv0[1]=50.0;   sv0[2]=50000.0;       sv0[3]=3.0;
              
        // Return values
        double fp0[4], fpe0[4];
        double chisqr0;
        int ndf0;
        TF1 *fitVAEnergy0 = langaufit(hist0,fr0,sv0,pllo0,plhi0,fp0,fpe0,&chisqr0,&ndf0);
        fitVAEnergy0->SetRange(0,200);

        double SNRPeak0, SNRFWHM0;
        langaupro(fp0,SNRPeak0,SNRFWHM0);
                
        // --- Hist1 --- //
        
        std::cout << "hist1 std dev " << hist1->GetStdDev() << std::endl;
        std::cout << "hist1 mean " << hist1->GetMean() << std::endl;
        std::cout << "hist1 integral " << hist1->Integral() << std::endl;
        
        // Setting fit range and start values
        double fr1[2], sv1[4], pllo1[4], plhi1[4];
        fr1[0] = 20.;   fr1[1] = 150.;
        sv1[0] = 2.;    sv1[1] = 40.;        sv1[2] = 1e5;       sv1[3] = 1.;
        pllo1[0] = 0.;  pllo1[1] = 10.;      pllo1[2] = 1.0;     pllo1[3] = 0.;
        plhi1[0] = 10.; plhi1[1] = 1e5;      plhi1[2] = 1e6;     plhi1[3] = 1e5;
               
        // Return values
        double fp1[4], fpe1[4];
        double chisqr1;
        int ndf1;
        TF1 *fitVAEnergy1 = langaufit(hist1,fr1,sv1,pllo1,plhi1,fp1,fpe1,&chisqr1,&ndf1);
        fitVAEnergy1->SetRange(0,200);

        double SNRPeak1, SNRFWHM1;
        langaupro(fp1,SNRPeak1,SNRFWHM1);
     
        //hMPV0->Fill(fp0[1]); // MPVs for eta region 0
        //hMPV1->Fill(fp1[1]); // MPVs for eta region 1 

        TCanvas *c1 = new TCanvas("c1",histName.c_str(),800,600);

        // Global style settings
        gStyle->SetOptStat(1111);
        gStyle->SetOptFit(111);
        gStyle->SetLabelSize(0.03,"x");
        gStyle->SetLabelSize(0.03,"y");
               
        hist0->SetTitle(histName.c_str());
 
        hist0->SetMaximum((hist0->GetMaximum())*1.1);
        hist0->Draw("hist");
        fitVAEnergy0->SetLineColor(kRed);
        fitVAEnergy0->Draw("hist same"); 
        
        hist1->Draw("hist sames");
        fitVAEnergy1->SetLineColor(kBlue);
        fitVAEnergy1->Draw("hist sames");
        
        gPad->Update(); 
        TPaveStats* sb0=(TPaveStats*)hist0->FindObject("stats");
        sb0->SetX1NDC(.65);
        sb0->SetX2NDC(.85);
        sb0->SetY1NDC(.65);
        sb0->SetY2NDC(.85);
        sb0->SetTextColor(kRed);
        TPaveStats* sb1=(TPaveStats*)hist1->FindObject("stats");
        sb1->SetX1NDC(.65);
        sb1->SetX2NDC(.85);
        sb1->SetY1NDC(.4);
        sb1->SetY2NDC(.6);
        sb1->SetTextColor(kBlue);
        
        c1->Print();
        //c1->Write();

        /**----- Fitting end ------*/


        //outFile->WriteTObject(hist);
        outFile->WriteTObject(c1);
        //outFile->WriteTObject(hMPV0);
        //outFile->WriteTObject(hMPV1);
        //delete hist;
        delete c1;
        count++;
        std::cout << "---------------------------------------------" << std::endl; 

     }
 
 
     //outFile->WriteTObject(hMPV0);
     //outFile->WriteTObject(hMPV1);
     outFile->cd();
     outFile->Write();
     outFile->Close();

    return 0;
} // end of main

