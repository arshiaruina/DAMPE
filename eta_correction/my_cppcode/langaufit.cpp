#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"
#include "TGaxis.h"
#include "TRandom.h"
#include "TMath.h"
#include "TF1.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TPaveStats.h"

using namespace std;

//-----------------------------------------------------------------------
//
// Convoluted Landau and Gaussian Fitting Function
//         (using ROOT's Landau and Gauss functions)
//
//  Based on a Fortran code by R.Fruehwirth (fruhwirth@hephy.oeaw.ac.at)
//  Adapted for C++/ROOT by H.Pernegger (Heinz.Pernegger@cern.ch) and
//   Markus Friedl (Markus.Friedl@cern.ch)
//
//  to execute this example, do:
//  root > .x langaus.C
// or
//  root > .x langaus.C++
//
//-----------------------------------------------------------------------

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


TF1* langaufit(TH1D *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF, Int_t flag)
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
   //sprintf(FunName,"Fitfcn_%s",his->GetName());
   //Modified AR-28.10.2020
   if(flag==0) // fit func for proton
      sprintf(FunName,"Fitfcn_p_%s",his->GetName());
   if(flag==1) // fit func for helium
      sprintf(FunName,"Fitfcn_He_%s",his->GetName());
   TF1 *ffitold = (TF1*)gROOT->GetListOfFunctions()->FindObject(FunName);
   if (ffitold) delete ffitold;
   TF1 *ffit = new TF1(FunName,langaufun,fitrange[0],fitrange[1],4);
   ffit->SetParameters(startvalues);
   ffit->SetParNames("Width","MP","Area","GSigma");
   for (i=0; i<4; i++) {
      ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);
   }
   his->Fit(FunName,"RB0");   // fit within specified range, use ParLimits, do not plot
   //his->Fit(FunName,"RB");   // Removing "0" does not make a difference, still does not plot! AR-28.10.2020
   ffit->GetParameters(fitparams);    // obtain fit parameters
   for (i=0; i<4; i++) {
      fiterrors[i] = ffit->GetParError(i);     // obtain fit parameter errors
   }
   ChiSqr[0] = ffit->GetChisquare();  // obtain chi^2
   NDF[0] = ffit->GetNDF();           // obtain ndf
   return (ffit);              // return fit function
}


Int_t langaupro(Double_t *params, Double_t &maxx, Double_t &FWHM) {
   // Searches for the location (x value) at the maximum of the
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



void LangauFit(string inFileName, string outFileName){

    cout << "inside LangauFit" << endl;
    cout << "debug 1" << endl;

    float BGOenergybins[]={0.02, 0.10, 0.25, 0.50, 1., 5.}; // TeV, hard-coded according to skim ranges
    int BGOenergyNbins = 6;
    TH1D* hClustETot[BGOenergyNbins];
    TH1D* hClustETotCorr[BGOenergyNbins];
    TCanvas* cFit[BGOenergyNbins];
    TCanvas* cFitCorr[BGOenergyNbins];
    for(int ibin=0; ibin<BGOenergyNbins; ibin++){
        hClustETot[ibin] = new TH1D(Form("hClustETot_%f",BGOenergybins[ibin]),Form("Uncorrected STK cluster energy for BGO energy bin %f",BGOenergybins[ibin]),500,0.,500.);
        hClustETotCorr[ibin] = new TH1D(Form("hClustETotCorr_%f",BGOenergybins[ibin]),Form("Corrected STK cluster energy for BGO energy bin %f",BGOenergybins[ibin]),500,0.,500.);
        cFit[ibin] = new TCanvas(Form("cFitsToUncorrPeaks_%f",BGOenergybins[ibin]),Form("Fits to uncorrected energy peaks for BGO energy bin %f",BGOenergybins[ibin]),800,600);
        cFitCorr[ibin] = new TCanvas(Form("cFitsToCorrPeaks_%f",BGOenergybins[ibin]),Form("Fits to corrected energy peaks for BGO energy bin %f",BGOenergybins[ibin]),800,600);
    }
    
    TH1F *hMPVproton = new TH1F("hMPVproton","MPV of uncorrected proton peak",6,0.,6.);
    TH1F *hMPVhelium = new TH1F("hMPVhelium","MPV of uncorrected helium peak",6,0.,6.);
    TH1F *hMPVprotonCorr = new TH1F("hMPVprotonCorr","MPV of corrected proton peak",6,0.,6.);
    TH1F *hMPVheliumCorr = new TH1F("hMPVheliumCorr","MPV of corrected helium peak",6,0.,6.);
   
    TFile *inFile = new TFile(inFileName.c_str());
    TFile *outFile = new TFile(outFileName.c_str(), "RECREATE"); 

    TTree *tree = (TTree*) inFile->Get("T");
    tree->Print();

    float BGOenergy = 0.;
    float ClustETot = 0.;
    float ClustETotCorr = 0.;

    tree->SetBranchAddress("BGOenergy",&BGOenergy);
    tree->SetBranchAddress("ClustETot",&ClustETot);
    tree->SetBranchAddress("ClustETotCorr",&ClustETotCorr);

    // --- Histograms of (corrected and non-corrected) STK energy in BGO energy bins ---//
    cout << "debug 2" << endl;
    cout << "entries: " << tree->GetEntries() << endl;

    for (int ientry = 0; ientry < tree->GetEntries(); ientry++) {
    
        tree->GetEntry(ientry);

        for (int ibin=0; ibin<BGOenergyNbins; ibin++){
            if(ibin < BGOenergyNbins-1){
                if((BGOenergy > BGOenergybins[ibin]*1e6) && (BGOenergy < BGOenergybins[ibin+1]*1e6)){
                    cout << "what!!! " << BGOenergy << endl;
                    hClustETot[ibin]->Fill(ClustETot);
                    hClustETotCorr[ibin]->Fill(ClustETotCorr);
                }
            }
            else { //last bin
                if(BGOenergy > BGOenergybins[ibin]*10e6){
                    hClustETot[ibin]->Fill(ClustETot);
                    hClustETotCorr[ibin]->Fill(ClustETotCorr);
                }
            }
        }
    }

    cout << "debug 3" << endl;

    for (int ibin=0; ibin<BGOenergyNbins; ibin++){  
        cout << ibin << " " << BGOenergyNbins << endl;
        hClustETot[ibin]->Draw();
        hClustETotCorr[ibin]->Draw();
        outFile->WriteTObject(hClustETot[ibin]);
        outFile->WriteTObject(hClustETotCorr[ibin]);
    }

    cout << "debug 4" << endl;

    // --- Langau fit for the uncorrected p and He peaks --- //
    
    for (int ibin=0; ibin<BGOenergyNbins; ibin++){

        cout << "checkpoint 1" << endl;

        TH1D *hist = hClustETot[ibin];
        if(hist->GetEntries() == 0) continue;

        // --- Proton Peak --- //

        cout << "--------------------------" << endl;
        std::cout << "Fitting uncorrected proton peak in BGO energy bin " << BGOenergybins[ibin] << std::endl;
        cout << "--------------------------" << endl;
        //std::cout << "hist std dev " << hist->GetStdDev() << std::endl;
        //std::cout << "hist mean " << hist->GetMean() << std::endl;
        //std::cout << "hist integral " << hist->Integral() << std::endl;
    

        // Setting fit range and start values
        double fr0[2], sv0[4], pllo0[4], plhi0[4];
        fr0[0] = 20.;      fr0[1] = 120.;
        sv0[0]=40;      sv0[1]=60.0;    sv0[2]=20000.0;       sv0[3]=10.0;
        pllo0[0]=20;    pllo0[1]=40.0;  pllo0[2]=1.0;         pllo0[3]=0.;
        plhi0[0]=60;   plhi0[1]=80.0;  plhi0[2]=1000000.0;   plhi0[3]=20.;
        
        // Return values
        double fp0[4], fpe0[4];
        double chisqr0;
        int ndf0;
        int flag0=0;
        TF1 *fitFcnProton = langaufit(hist,fr0,sv0,pllo0,plhi0,fp0,fpe0,&chisqr0,&ndf0,flag0);
        fitFcnProton->SetRange(0,500);
        
        double SNRPeak0, SNRFWHM0;
        langaupro(fp0,SNRPeak0,SNRFWHM0);
       
        // --- Helium Peak --- //
        
        cout << "--------------------------" << endl;
        std::cout << "Fitting uncorrected He peak in BGO energy bin " << BGOenergybins[ibin] << std::endl;
        cout << "--------------------------" << endl;
        //std::cout << "hist std dev " << hist->GetStdDev() << std::endl;
        //std::cout << "hist mean " << hist->GetMean() << std::endl;
        //std::cout << "hist integral " << hist->Integral() << std::endl;
        
        // Setting fit range and start values
        double fr1[2], sv1[4], pllo1[4], plhi1[4];
        fr1[0] = 220.;   fr1[1] = 380.;
        sv1[0] = 80;        sv1[1] = 250.;        sv1[2] = 10000;       sv1[3] = 10.;
        pllo1[0] = 60.;    pllo1[1] = 230.;      pllo1[2] = 1.0;     pllo1[3] = 0.;
        plhi1[0] = 100.;    plhi1[1] = 270.;      plhi1[2] = 1e6;     plhi1[3] = 20.;
        
        // Return values
        double fp1[4], fpe1[4];
        double chisqr1;
        int ndf1;
        int flag1=1;
        TF1 *fitFcnHelium = langaufit(hist,fr1,sv1,pllo1,plhi1,fp1,fpe1,&chisqr1,&ndf1,flag1);
        fitFcnHelium->SetRange(0,500);
        
        double SNRPeak1, SNRFWHM1;
        langaupro(fp1,SNRPeak1,SNRFWHM1);

        cout << "checkpoint 2" << endl;

        //c1[iladder][iva] = new TCanvas(histName.c_str(),histName.c_str(),800,600);
        //if(countVA==0){ // open pdf
        //    c1[iladder][iva]->Print("plots.pdf[","pdf");
        //}
    
        cFit[ibin]->cd();

        // Global style settings
        gStyle->SetOptStat(1111);
        gStyle->SetOptFit(111);
        gStyle->SetLabelSize(0.03,"x");
        gStyle->SetLabelSize(0.03,"y");
               
        //hist->SetTitle(histName.c_str());

        hist->Sumw2(kFALSE); 
        hist->SetMaximum((hist->GetMaximum())*1.1);
        hist->Draw("hist");
        //hist->SetLineColor(kRed);
        fitFcnProton->Draw("sames"); 
        fitFcnProton->SetLineColor(kRed);
        //cFit[ibin]->Update();        
        //hist->Draw("hist sames");
        //hist->SetLineColor(kBlue);
        fitFcnHelium->Draw("sames");
        fitFcnHelium->SetLineColor(kBlue);
        //fitFcnProton->SetLineColor(kRed);
        //fitFcnProton->Draw("sames"); 
        cFit[ibin]->Update();        

        cout << "checkpoint 3" << endl;

        //hist->GetListOfFunctions()->ls();    

        //gPad->Update(); 
        //TPaveStats* sb0=(TPaveStats*)hist->FindObject("stats");
        //sb0->SetX1NDC(.65);
        //sb0->SetX2NDC(.85);
        //sb0->SetY1NDC(.65);
        //sb0->SetY2NDC(.85);
        //sb0->SetTextColor(kRed);
        //TPaveStats* sb1=(TPaveStats*)hist->FindObject("stats");
        //sb1->SetX1NDC(.65);
        //sb1->SetX2NDC(.85);
        //sb1->SetY1NDC(.4);
        //sb1->SetY2NDC(.6);
        //sb1->SetTextColor(kBlue);

        //c1[iladder][iva]->Print("plots.pdf","pdf");
        outFile->WriteTObject(cFit[ibin]);

        hMPVproton->SetBinContent(ibin,fp0[1]);
        hMPVhelium->SetBinContent(ibin,fp1[1]);

        //if(countVA==(nVA-1)){ //close pdf
        //    c1[iladder][iva]->Print("plots.pdf]","pdf");
        //}


        cout << "checkpoint 4" << endl;
    } 

    cout << "debug 5" << endl;

    // --- Langau fit for the corrected p and He peaks --- //
    
    for (int ibin=0; ibin<BGOenergyNbins; ibin++){

        TH1D *hist = hClustETotCorr[ibin];
        if(hist->GetEntries() == 0) continue;
        
        // --- Proton Peak --- //

        cout << "--------------------------" << endl;
        std::cout << "Fitting corrected proton peak in BGO energy bin " << BGOenergybins[ibin] << std::endl;
        cout << "--------------------------" << endl;
        //std::cout << "hist std dev " << hist->GetStdDev() << std::endl;
        //std::cout << "hist mean " << hist->GetMean() << std::endl;
        //std::cout << "hist integral " << hist->Integral() << std::endl;
        
        // Setting fit range and start values
        double fr0[2], sv0[4], pllo0[4], plhi0[4];
        fr0[0] = 20.;      fr0[1] = 120.;
        sv0[0]=40;      sv0[1]=60.0;    sv0[2]=20000.0;       sv0[3]=10.0;
        pllo0[0]=20;    pllo0[1]=40.0;  pllo0[2]=1.0;         pllo0[3]=0.;
        plhi0[0]=60;   plhi0[1]=80.0;  plhi0[2]=1000000.0;   plhi0[3]=20.;

        // Return values
        double fp0[4], fpe0[4];
        double chisqr0;
        int ndf0;
        int flag0=0;
        TF1 *fitFcnProton = langaufit(hist,fr0,sv0,pllo0,plhi0,fp0,fpe0,&chisqr0,&ndf0,flag0);
        fitFcnProton->SetRange(0,500);
        
        double SNRPeak0, SNRFWHM0;
        langaupro(fp0,SNRPeak0,SNRFWHM0);
        
        // --- Helium Peak --- //
        
        cout << "--------------------------" << endl;
        std::cout << "Fitting corrected He peak in BGO energy bin " << BGOenergybins[ibin] << std::endl;
        cout << "--------------------------" << endl;
        //std::cout << "hist std dev " << hist->GetStdDev() << std::endl;
        //std::cout << "hist mean " << hist->GetMean() << std::endl;
        //std::cout << "hist integral " << hist->Integral() << std::endl;
        
        // Setting fit range and start values
        double fr1[2], sv1[4], pllo1[4], plhi1[4];
        fr1[0] = 220.;   fr1[1] = 380.;
        sv1[0] = 80;       sv1[1] = 250.;        sv1[2] = 10000;       sv1[3] = 10.;
        pllo1[0] = 60.;    pllo1[1] = 230.;      pllo1[2] = 1.0;     pllo1[3] = 0.;
        plhi1[0] = 100.;    plhi1[1] = 270.;      plhi1[2] = 1e6;     plhi1[3] = 20.;
        
        // Return values
        double fp1[4], fpe1[4];
        double chisqr1;
        int ndf1;
        int flag1=1;
        TF1 *fitFcnHelium = langaufit(hist,fr1,sv1,pllo1,plhi1,fp1,fpe1,&chisqr1,&ndf1,flag1);
        fitFcnHelium->SetRange(0,500);
        
        double SNRPeak1, SNRFWHM1;
        langaupro(fp1,SNRPeak1,SNRFWHM1);

        //c1[iladder][iva] = new TCanvas(histName.c_str(),histName.c_str(),800,600);
        //if(countVA==0){ // open pdf
        //    c1[iladder][iva]->Print("plots.pdf[","pdf");
        //}
    
        cFitCorr[ibin]->cd();

        // Global style settings
        gStyle->SetOptStat(1111);
        gStyle->SetOptFit(111);
        gStyle->SetLabelSize(0.03,"x");
        gStyle->SetLabelSize(0.03,"y");
               
        //hist->SetTitle(histName.c_str());

        hist->Sumw2(kFALSE);
        hist->SetMaximum((hist->GetMaximum())*1.1);
        hist->Draw("hist");
        //hist->SetLineColor(kRed);
        fitFcnProton->SetLineColor(kRed);
        fitFcnProton->Draw("sames"); 
        //cFitCorr[ibin]->Update();
        //hist->Draw("hist sames");
        //hist->SetLineColor(kBlue);
        fitFcnHelium->SetLineColor(kBlue);
        fitFcnHelium->Draw("sames");
        cFitCorr[ibin]->Update();

        //hist->GetListOfFunctions()->ls();

        //gPad->Update(); 
        //TPaveStats* sb0=(TPaveStats*)hist->FindObject("stats");
        //sb0->SetX1NDC(.65);
        //sb0->SetX2NDC(.85);
        //sb0->SetY1NDC(.65);
        //sb0->SetY2NDC(.85);
        //sb0->SetTextColor(kRed);
        //TPaveStats* sb1=(TPaveStats*)hist->FindObject("stats");
        //sb1->SetX1NDC(.65);
        //sb1->SetX2NDC(.85);
        //sb1->SetY1NDC(.4);
        //sb1->SetY2NDC(.6);
        //sb1->SetTextColor(kBlue);

        //c1[iladder][iva]->Print("plots.pdf","pdf");
        outFile->WriteTObject(cFitCorr[ibin]);
        
        hMPVprotonCorr->SetBinContent(ibin+1,fp0[1]);
        hMPVheliumCorr->SetBinContent(ibin+1,fp1[1]);

        //if(countVA==(nVA-1)){ //close pdf
        //    c1[iladder][iva]->Print("plots.pdf]","pdf");
        //}

        
        cout << "debug 6" << endl;

    } 

    outFile->WriteTObject(hMPVproton);
    outFile->WriteTObject(hMPVhelium);

    outFile->WriteTObject(hMPVprotonCorr);
    outFile->WriteTObject(hMPVheliumCorr);
    
    TCanvas *cMPVproton = new TCanvas("cMPVproton","MPVs of proton peak",800,600);

    hMPVproton->SetMarkerStyle(8);
    hMPVproton->SetMarkerSize(0.7);
    hMPVproton->SetMarkerColor(1);

    hMPVprotonCorr->SetMarkerStyle(8);
    hMPVprotonCorr->SetMarkerSize(0.7);
    hMPVprotonCorr->SetMarkerColor(2);

    //hMPVproton->SetMaximum(1.5*hMPVproton->GetMaximum());
    //hMPVproton->SetMinimum(0.5*hMPVproton->GetMinimum());

    hMPVproton->Draw("PL");
    hMPVprotonCorr->Draw("PL same");

    outFile->WriteTObject(cMPVproton);
    
    TCanvas *cMPVhelium = new TCanvas("cMPVhelium","MPVs of helium peak",800,600);

    hMPVhelium->SetMarkerStyle(8);
    hMPVhelium->SetMarkerSize(0.7);
    hMPVhelium->SetMarkerColor(1);

    hMPVheliumCorr->SetMarkerStyle(8);
    hMPVheliumCorr->SetMarkerSize(0.7);
    hMPVheliumCorr->SetMarkerColor(2);

    //hMPVhelium->SetMaximum(1.5*hMPVhelium->GetMaximum());
    //hMPVhelium->SetMinimum(0.5*hMPVhelium->GetMinimum());
    
    hMPVhelium->Draw("PL");
    hMPVheliumCorr->Draw("PL same");

    outFile->WriteTObject(cMPVhelium);

    outFile->cd();
    outFile->Write();
    outFile->Close();

}

int main(int argc, char *argv[]){

    using namespace std;

    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <file_in.root> <file_out.root>" << endl;
        cerr << "\t <file_in.root> -- root file that contains the T tree with useful variables like STK cluster energy, BGO energy etc. stored in its branches" << endl;
        cerr << "\t <file_out.root> -- output root file that will contain the langau fitted p and He peaks in STK cluster energy as well as MPVs and sigmas of the fits" << endl;
        return 1;
    }   
    cout << "debug start" << endl;
    LangauFit(argv[1], argv[2]);
    cout << "debug finish" << endl;

    return 0;
 
}
