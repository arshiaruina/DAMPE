/* Code to compute eta values
 * Author: Arshia Ruina
 * To compile, g++ -Wall `root-config --cflags --libs` -I$DMPSWSYS/include -L$DMPSWSYS/lib -lDmpEvent -lDmpService eta.cc -o ../submit/eta.exe
 * */

// example file
// /beegfs/users/ruina/out/20181019/merged/merged_160919_104734.root

#include "../inc/va_equalisation.h"
#include "mylangaus.C"

using namespace std;

int main(int argc, char** argv) {

/* TODO:
 * Access the generated (merged) file.root which has the histograms of the energy distributions for all VAs
 * langaufit accessess each histogram and makes the fit
 * Fit results (fitted energy plots of all the VAs) stored in new root file
 * Fit results (MPVs of all VAs) stored in array (1152 VAs x 2 eta-regions) and plotted for the 2 eta regions
 * Mean of MPV distributions stored (2 values for the 2 eta regions)
 * Correction factor for each VA for each eta region = Mean MPV for that eta region / MPV of the VA for that eta region
 * */

	TStopwatch sw;
	sw.Start();

    double MPV[1152][2] = {0.};
    TH1D *hMPV0 = new TH1D("hMPV0","MPV of all VAs in RO-strip region",100,40.,60.);
    TH1D *hMPV1 = new TH1D("hMPV1","MPV of all VAs in floating-strip region",100,30.,50.);

    /* ------ Accessing the histo names from the histolist file ------ */

	std::string inFileName = argv[1];
	TFile *inFile = new TFile(inFileName.c_str());
    
    std::ifstream histFile;
    histFile.open("histolist.txt");
   
    std::string histName;
    
    TFile *outFile = new TFile("../out/20181019/langaufit.root", "RECREATE");
    
    int i = 0;
    //while(!histFile.eof()){
            //std::cout << "debug . #" << i << std::endl;
    while(std::getline(histFile,histName) && i < 10){ // for debug run
    //while(std::getline(histFile,histName)){ // for analysis run

    //...for(int iName = 0; iName < 1152; iName+=2){

        //float progress = 100.0 * ((float) i) / ((float) 2304);
        //if (!(i % 10)) cout << setprecision(3) << " [ " << progress << " % ] \r";
        //std::cout << std::setprecision(3) << "[ " << progress << " % ] \r";
    
        TH1D *hist = (TH1D*)inFile->Get(histName.c_str());
        //...TH1D *hist = (TH1D*)inFile->Get(histoNames.at(iName).c_str());
    
        if(hist->GetEntries() == 0) continue;
        /* ---- Fitting ---- */

        //TFile *f   = new TFile("../out/20181019/merged/merged_160919_104734.root");
        //TH1D *hVAEnergy = (TH1D*) f->Get("hVAEnergyX_71_0_0");
    
        std::cout << "Fitting histogram " << histName << std::endl;
    
        /*
         Once again, here are the Landau * Gaussian parameters:
          par[0]=Width (scale) parameter of Landau density
          par[1]=Most Probable (MP, location) parameter of Landau density
          par[2]=Total area (integral -inf to inf, normalization constant)
          par[3]=Width (sigma) of convoluted Gaussian function
       
        Variables for langaufit call:
          his             histogram to fit
          fitrange[2]     lo and hi boundaries of fit range
          startvalues[4]  reasonable start values for the fit
          parlimitslo[4]  lower parameter limits
          parlimitshi[4]  upper parameter limits
          fitparams[4]    returns the final fit parameters
          fiterrors[4]    returns the final fit errors
          ChiSqr          returns the chi square
          NDF             returns ndf
     
        */
    
        // Setting fit range and start values
        double fr[2], sv[4], pllo[4], plhi[4];
        fr[0] = 20.;         fr[1] = 200.;
        sv[0] = 2.;         sv[1] = 60.;        sv[2] = 50000.0;        sv[3] = 5.;
        pllo[0] = 0.5;      pllo[1] = 40.;      pllo[2] = 1.0;          pllo[3] = 1.;
        plhi[0] = 5.;       plhi[1] = 80.;      plhi[2] = 100000.0;     plhi[3] = 10.;
    
        // Return values
        double fp[4], fpe[4];
        double chisqr;
        int ndf;
        TF1 *fitVAEnergy = langaufit(hist,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
    
        double SNRPeak, SNRFWHM;
        langaupro(fp,SNRPeak,SNRFWHM);
    
        //std::cout << "Fitting done \n Plotting results..." << std::endl;
 
        if(i%4 == 0 || i%4 == 1){
            //MPV[i][0] = fp[1]; // MPVs for eta region 0
            hMPV0->Fill(fp[1]); // MPVs for eta region 0
        }
        if(i%4 == 2 || i%4 == 3){
            //MPV[i][1] = fp[1]; // MPVs for eta region 1
            hMPV1->Fill(fp[1]); // MPVs for eta region 0
        }

        TCanvas *c1 = new TCanvas("c1",histName.c_str(),800,600);
        //TCanvas can[i] = new TCanvas(Form(histName.c_str()),histName.c_str(),800,600);
 
        // Global style settings
        gStyle->SetOptStat(1111);
        gStyle->SetOptFit(111);
        gStyle->SetLabelSize(0.03,"x");
        gStyle->SetLabelSize(0.03,"y");
    
        hist->SetMaximum((hist->GetMaximum())*1.1);
        hist->Draw("hist");
        fitVAEnergy->Draw("hist same");

        /* ---- Fitting ends ---- */

        //hist->Draw(); 
        hist->Write();
        fitVAEnergy->Write();
        c1->Print();
        c1->Write();
        //outFile->Write();
        delete hist;
        delete fitVAEnergy;
        delete c1;
        i++;
    }
    
//    for(int j = 0; j < 1152; j++){
//        hMPV0->Fill(MPV[j][0]);
//        hMPV1->Fill(MPV[j][1]);
//    }
    hMPV0->Write();
    hMPV1->Write();
    std::cout << i << std::endl; 
	outFile->cd();
	outFile->Write();
	outFile->Close();
    histFile.close();
	sw.Stop();
    sw.Print();

    return 0;
} // end of main

