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
#include "TPaveStats.h"

//#include "../inc/va_equalisation.h"
//#include "mylangaus.C"
#include "../inc/correctionfactor.h"

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


std::pair <double, double> MPVGausFit(std::string f) {

    //std::string inFileName = "../out/20181001_20181009/langaufit_231019_112429.root";
    //std::string inFileName = "../out/20181001_20181009/langaufit_231019_112429_FMVAs.root";
    //TFile *inFile = new TFile(inFileName.c_str());
    TFile *inFile = new TFile(f.c_str());
    //TFile *outFile = new TFile("../out/20181001_20181009/MPV_langaufit_231019_112429.root", "RECREATE"); 
    //TFile *outFile = new TFile("../out/20181001_20181009/gausfitMPV_20181001_20181009_FMVAs.root", "RECREATE"); 
    TFile *outFile = new TFile("test_gausfitMPV_20181001_20181009_FMVAs.root", "RECREATE"); 
    //TFile *outFile = new TFile("gausfitMPV_20181001_20181009_FMVAs.root", "RECREATE"); 

    string histName0 = "hMPV0";
    string histName1 = "hMPV1";
    TH1D *hist0 = (TH1D*)inFile->Get(histName0.c_str());
    TH1D *hist1 = (TH1D*)inFile->Get(histName1.c_str());

    TF1 *fGaus0 = new TF1("fGaus0","gaus");
    fGaus0->SetParameters(1,0,1);
    hist0->Fit(fGaus0,"0");
    fGaus0->SetRange(40,60);

    TF1 *fGaus1 = new TF1("fGaus1","gaus");
    fGaus1->SetParameters(1,0,1);
    hist1->Fit(fGaus1,"0");
    fGaus1->SetRange(30,50);

    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(1111);
    gStyle->SetLabelSize(0.03,"x");
    gStyle->SetLabelSize(0.03,"y");

    TCanvas *c0 = new TCanvas("c0","c0");
    hist0->Draw("hist");
    fGaus0->SetLineColor(kBlack);
    fGaus0->Draw("hist sames");
    gPad->Update();
    TPaveStats* sb0=(TPaveStats*)hist0->FindObject("stats");
    sb0->SetX1NDC(.65);
    sb0->SetX2NDC(.85);
    sb0->SetY1NDC(.65);
    sb0->SetY2NDC(.85);
    sb0->SetTextColor(kBlack);

    TCanvas *c1 = new TCanvas("c1","c1");
    hist1->Draw("hist");
    fGaus1->SetLineColor(kBlack);
    fGaus1->Draw("hist sames");
    gPad->Update();
    TPaveStats* sb1=(TPaveStats*)hist1->FindObject("stats");
    sb1->SetX1NDC(.65);
    sb1->SetX2NDC(.85);
    sb1->SetY1NDC(.65);
    sb1->SetY2NDC(.85);
    sb1->SetTextColor(kBlack);

    c0->Write();
    c1->Write();

    outFile->Write();
    outFile->Close();

    delete c0;
    delete c1;

    double eqParam0 = fGaus0->GetParameter(1);
    double eqParam1 = fGaus1->GetParameter(1);

    std::cout << "eqParam0: " << eqParam0 << " eqParam1: " << eqParam1 << std::endl;

    return std::make_pair(eqParam0, eqParam1);

} 


//int main(int argc, char** argv) {
int main(){

    TH1D *hMPV0 = new TH1D("hMPV0","MPV of all VAs in RO-strip region",100,25.,85.);
    TH1D *hMPV1 = new TH1D("hMPV1","MPV of all VAs in floating-strip region",100,20.,60.);
    TH1D *hCorrFac0 = new TH1D("hCorrFac0","Correction factors in RO-strip region",100,0.,2.);
    TH1D *hCorrFac1 = new TH1D("hCorrFac1","Correction factors in floating-strip region",100,0,2.);
    TH1D *hCorrFacDiff = new TH1D("hCorrFacDiff","#Delta(Correction Factors)",100,-1.,1.);

    std::ifstream histFile;
    histFile.open("/beegfs/users/ruina/VAequalisation/resources/histolist_new.txt");
    std::string histName;

    //std::string inFileName = argv[1];
    std::string inFileName = "/beegfs/users/ruina/VAequalisation/out/20181001_20181009/merged/merged_231019_112429.root";
    TFile *inFile = new TFile(inFileName.c_str());
    // - TList *list = inFile->GetListOfKeys();
    //TFile *outFile = new TFile("../out/root_forum.root", "recreate");
    //TFile *outFile = new TFile("../out/20181019/langaufit.root", "RECREATE"); 
    //TFile *outFile = new TFile("../out/201810_wk1/langaufit.root", "RECREATE"); 
    //TFile *outFile = new TFile("../out/20181001_20181009/langaufit_231019_112429.root", "RECREATE"); 
    //TFile *outFile = new TFile("../out/20181001_20181009/langaufit_231019_112429_test.root", "RECREATE"); 
    
    std::string outFileName = "test_langaufit_231019_112429_FMVAs.root";
    //std::string outFileName = "langaufit_231019_112429_FMVAs.root";
    //std::string outFileName = "../out/20181001_20181009/langaufit_231019_112429_FMVAs.root";
    TFile *outFile = new TFile(outFileName.c_str(), "RECREATE"); 

    int iladder, iva = 0;
    int nVA = 2;
    int nLADDER = 2;
    int countVA = 0;
    TCanvas *c1[nLADDER][nVA];
    double MPV0[nLADDER][nVA] = {0.};
    double MPV1[nLADDER][nVA] = {0.};
    
    std::cout << "MPV0 and MPV1 elements... " << std::endl;
    for(int wt=0; wt<nLADDER; wt++){
        for(int wtf=0; wtf<nVA; wtf++){
            std::cout << MPV0[wt][wtf] << "\t" << MPV1[wt][wtf] << std::endl;
        }
    }

    double corrFac0[nLADDER][nVA] = {0.};
    double corrFac1[nLADDER][nVA] = {0.};
    std::vector<string> EQMladders = {"12","13","14","15","36","37","38","39","66","67","68"};

    while(std::getline(histFile,histName) && countVA < nVA){ // for debug run
    //while(std::getline(histFile,histName)){ // for analysis run
    
        std::size_t delim1 = histName.find("_");      
        std::size_t delim2 = histName.find("_",delim1+1);      
        iladder = std::stoi(histName.substr(delim1+1,delim2-delim1-1));
        iva = std::stoi(histName.substr(delim2+1));
        std::cout << "iladder ... " << iladder << std::endl;
        std::cout << "iva ... " << iva << std::endl;

 
        string histName0 = histName + "_0";
        string histName1 = histName + "_1";
        TH1D *hist0 = (TH1D*)inFile->Get(histName0.c_str()); 
        TH1D *hist1 = (TH1D*)inFile->Get(histName1.c_str()); 
        std::cout << histName << std::endl;
        if(histName == "hEtaX" || histName == "hEtaY" || hist0->GetEntries() == 0 || hist1->GetEntries() == 0) continue;
        //if(!(histName == "hVAEnergyY_141_5" || histName == "hVAEnergyY_106_2")) continue;
        //if(!(histName == "hVAEnergyX_154_2" || histName == "hVAEnergyY_33_5")) continue;
        //if(!(histName == "hVAEnergyX_154_2")) continue;
        //if(!(histName == "hVAEnergyX_154_2")) continue;
        c1[iladder][iva] = new TCanvas(histName.c_str(),histName.c_str(),800,600);
        if(countVA==0){ // open pdf
            std::cout << "plots.pdf created?" << std::endl;
            c1[iladder][iva]->Print("plots.pdf[","pdf");
        }

        std::cout << "debug 1" << std::endl;

        /**************************************************************
        
            Langau fits on the energy distributions of all VAs
            Generates 1152 plots with 2 distributions on each
            for the two eta regions
        
        **************************************************************/
    
        // --- Hist0 --- //             
                 
        std::cout << "hist0 std dev " << hist0->GetStdDev() << std::endl;
        std::cout << "hist0 mean " << hist0->GetMean() << std::endl;
        std::cout << "hist0 integral " << hist0->Integral() << std::endl;

        // Setting fit range and start values
        double fr0[2], sv0[4], pllo0[4], plhi0[4];
        fr0[0] = 40.;      fr0[1] = 150.;
        sv0[0]=1.8;   sv0[1]=50.0;   sv0[2]=50000.0;       sv0[3]=3.0;
        pllo0[0]=0.5; pllo0[1]=10.0; pllo0[2]=1.0;         pllo0[3]=0.;
        plhi0[0]=9.0; plhi0[1]=70.0; plhi0[2]=1000000.0;   plhi0[3]=15.;
              
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
        sv1[0] = 2.;    sv1[1] = 35.;        sv1[2] = 1e5;       sv1[3] = 1.;
        pllo1[0] = 0.;  pllo1[1] = 10.;      pllo1[2] = 1.0;     pllo1[3] = 0.;
        plhi1[0] = 10.; plhi1[1] = 1e5;      plhi1[2] = 1e6;     plhi1[3] = 10.;
               
        // Return values
        double fp1[4], fpe1[4];
        double chisqr1;
        int ndf1;
        TF1 *fitVAEnergy1 = langaufit(hist1,fr1,sv1,pllo1,plhi1,fp1,fpe1,&chisqr1,&ndf1);
        fitVAEnergy1->SetRange(0,200);

        double SNRPeak1, SNRFWHM1;
        langaupro(fp1,SNRPeak1,SNRFWHM1);
     
        //TCanvas *c1 = new TCanvas("c1",histName.c_str(),800,600);

    
        std::cout << "debug 2" << std::endl;

        // Global style settings
        gStyle->SetOptStat(1111);
        gStyle->SetOptFit(111);
        gStyle->SetLabelSize(0.03,"x");
        gStyle->SetLabelSize(0.03,"y");
               
        hist0->SetTitle(histName.c_str());
 
        hist0->SetMaximum((hist0->GetMaximum())*1.1);
        hist0->Draw("hist");
        hist0->SetLineColor(kRed);
        fitVAEnergy0->SetLineColor(kRed);
        fitVAEnergy0->Draw("hist same"); 
        
        hist1->Draw("hist sames");
        hist1->SetLineColor(kBlue);
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

        /************************************************
        
                   Saving the histograms and fits
        ************************************************/
        
        std::cout << "debug 3" << std::endl;
        if(countVA > 0 && countVA < (nVA-1))
        c1[iladder][iva]->Print("plots.pdf","pdf");
        std::cout << "debug 4" << std::endl;
      
        //c1[iva]->Write();
        //c1[iva]->Print(); //generates .ps files in current directory, but empty?
        //c1->Print();
        //c1->Write();

        //outFile->WriteTObject(hist);
        outFile->WriteTObject(c1[iladder][iva]);
        //outFile->WriteTObject(c1);
        //outFile->WriteTObject(hMPV0);
        //outFile->WriteTObject(hMPV1);
        //delete hist;
        //delete c1;

        std::cout << "debug 5" << std::endl;

        /***********************************************************************
        
            Mean of the fits stored in MPV0[1152] and MPV1[1152]
            to be used later to compute the correction factors for each VA
        
        ***********************************************************************/
        
        //MPV0[iva] = fp0[1];
        //MPV1[iva] = fp1[1];
        
        std::cout << "iladder " << iladder << std::endl;
        std::cout << "iva " << iva << std::endl;
        std::cout << "MPV0[iladder][iva] " << MPV0[iladder][iva] << std::endl;
        std::cout << "MPV1[iladder][iva] " << MPV1[iladder][iva] << std::endl;
        std::cout << "fp0[1] " << fp0[1] << std::endl;
        std::cout << "fp1[1] " << fp1[1] << std::endl;

        MPV0[iladder][iva] = fp0[1];
        MPV1[iladder][iva] = fp1[1];

        /***********************************************************************
    
            Mean of the fits (for FM VAs only) stored in hMPV0 and hMPV1
            to make a gaussian fit to them and use the mean as "eq. param."
      
        ***********************************************************************/
       
        std::cout << "debug 6" << std::endl;
 
        std::vector<string>::iterator it;
        it = std::find (EQMladders.begin(), EQMladders.end(), std::to_string(iladder));    
    
        std::cout << "debug 7" << std::endl;

        //it = std::find (EQMladders.begin(), EQMladders.end(), histName.substr(11,2));    
        //if(it != EQMladders.end() && histName.substr(13,1) == "_") //EQM ladder found
        if(it == EQMladders.end() || histName.substr(13,1) != "_") { //FM ladder found
        //    continue;
        //else {
            hMPV0->Fill(fp0[1]); // MPVs for eta region 0
            hMPV1->Fill(fp1[1]); // MPVs for eta region 1 
        }    

        std::cout << "-------  Finished working on VA " << countVA << "  ---------" << std::endl; 

        if(countVA==(nVA-1)){ //close pdf
            c1[iladder][iva]->Print("plots.pdf]","pdf");
        }

        countVA++;

        std::cout << "final debug" << std::endl;

    } // end of loop over all histnames

     
 
    /******************************************************
    
        Saving the hMPV0 and hMPV1 histograms and fits 
        to the same output file, making a gaussian fit
        and printing the values of the two eq. params.
        (which are the means of the two gaussian fits)
        for the two eta regions)
 
    ******************************************************/
        
    outFile->WriteTObject(hMPV0);
    outFile->WriteTObject(hMPV1);
    outFile->cd();
    outFile->Write();
    outFile->Close();

    std::pair<double, double> eqParams = MPVGausFit(outFileName);
    std::cout << eqParams.first << std::endl;
    std::cout << eqParams.second << std::endl;

    /******************************************************
    
        Computing the correction factors for the two
        eta regions and the difference between the values,
        saving them in histograms.

    ******************************************************/

    //TFile *outFileCorrFac = new TFile("corrFac.root", "RECREATE");
    TFile *outFileCorrFac = new TFile("test_corrFac.root", "RECREATE");

    //for(int iva = 0; iva < 1152; iva++){
    for(int iladder = 0; iva < nLADDER; iladder++){
        for(int iva = 0; iva < nVA; iva++){
            corrFac0[iladder][iva] = eqParams.first/MPV0[iladder][iva];
            corrFac1[iladder][iva] = eqParams.second/MPV1[iladder][iva];
            hCorrFac0->Fill(corrFac0[iladder][iva]);
            hCorrFac1->Fill(corrFac1[iladder][iva]);
            hCorrFacDiff->Fill(corrFac0[iladder][iva]-corrFac1[iladder][iva]);
        }
    }

    outFileCorrFac->WriteTObject(hCorrFac0);
    outFileCorrFac->WriteTObject(hCorrFac1);
    outFileCorrFac->WriteTObject(hCorrFacDiff);
    outFileCorrFac->cd();
    outFileCorrFac->Write();
    outFileCorrFac->Close();

    std::cout << " corr facs computed and stored" << std::endl;    

    return 0;
} // end of main
