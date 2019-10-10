/* Code to compute eta values
 * Author: Arshia Ruina
 * To compile, g++ -Wall `root-config --cflags --libs` -I$DMPSWSYS/include -L$DMPSWSYS/lib -lDmpEvent -lDmpService eta.cc -o ../submit/eta.exe
 * */

// example file
// /beegfs/users/ruina/VAequalisation/out/20181019/merged/merged_160919_104734.root

#include "../inc/va_equalisation.h"
#include "mylangaus.C"
#include "TPaveStats.h"
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

	//TStopwatch sw;
	//sw.Start();

    //double MPV[1152][2] = {0.};
    TH1D *hMPV0 = new TH1D("hMPV0","MPV of all VAs in RO-strip region",100,40.,60.);
    TH1D *hMPV1 = new TH1D("hMPV1","MPV of all VAs in floating-strip region",100,30.,50.);


//-----------------------------------------------------------//
    std::ifstream histFile;
    histFile.open("histolist_new.txt");
    std::string histName;

    //int count = 0;
    std::string inFileName = argv[1];
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

    //while(std::getline(histFile,histName) && count < 20){ // for debug run
    while(std::getline(histFile,histName)){ // for analysis run
        
        //if(count%2 == 1) continue;

        string histName0 = histName + "_0";
        string histName1 = histName + "_1";
        //TH1D *hist = (TH1D*)inFile->Get(histName.c_str()); 
        TH1D *hist0 = (TH1D*)inFile->Get(histName0.c_str()); 
        TH1D *hist1 = (TH1D*)inFile->Get(histName1.c_str()); 
        if(histName == "hEtaX" || histName == "hEtaY" || hist0->GetEntries() == 0 || hist1->GetEntries() == 0) continue;

        /*------ Fitting start ------*/

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
    
                 // --- Hist0 --- Setting fit range and start values
                 double fr[2], sv[4], pllo[4], plhi[4];
                 fr[0] = 20.;         fr[1] = 200.;
                 sv[0] = 2.;         sv[1] = 50.;        sv[2] = 50000.0;        sv[3] = 5.;
                 pllo[0] = 0.5;      pllo[1] = 40.;      pllo[2] = 1.0;          pllo[3] = 1.;
                 plhi[0] = 5.;       plhi[1] = 80.;      plhi[2] = 100000.0;     plhi[3] = 10.;
                
                 // Return values
                 double fp0[4], fpe0[4];
                 double chisqr0;
                 int ndf0;
                 TF1 *fitVAEnergy0 = langaufit(hist0,fr,sv,pllo,plhi,fp0,fpe0,&chisqr0,&ndf0);

                 double SNRPeak0, SNRFWHM0;
                 langaupro(fp0,SNRPeak0,SNRFWHM0);
                
                 // --- Hist1 --- Setting fit range and start values
                 double fr[2], sv[4], pllo[4], plhi[4];
                 fr[0] = 20.;         fr[1] = 200.;
                 sv[0] = 2.;         sv[1] = 50.;        sv[2] = 50000.0;        sv[3] = 5.;
                 pllo[0] = 0.5;      pllo[1] = 40.;      pllo[2] = 1.0;          pllo[3] = 1.;
                 plhi[0] = 5.;       plhi[1] = 80.;      plhi[2] = 100000.0;     plhi[3] = 10.;
                
                 // Return values
                 double fp1[4], fpe1[4];
                 double chisqr1;
                 int ndf1;
                 TF1 *fitVAEnergy1 = langaufit(hist1,fr,sv,pllo,plhi,fp1,fpe1,&chisqr1,&ndf1);

                 double SNRPeak1, SNRFWHM1;
                 langaupro(fp1,SNRPeak1,SNRFWHM1);
                
                 //std::cout << "Fitting done \n Plotting results..." << std::endl;
 
                 //if(histName.back() == '0') {
                 //    hMPV0->Fill(fp[1]); // MPVs for eta region 0
                 //}
                 //else{
                 //    hMPV1->Fill(fp[1]); // MPVs for eta region 0
                 //}
     
                 hMPV0->Fill(fp0[1]); // MPVs for eta region 0
                 hMPV1->Fill(fp1[1]); // MPVs for eta region 1 
                 

                 TCanvas *c1 = new TCanvas("c1",histName.c_str(),800,600);
                 //TCanvas can[i] = new TCanvas(Form(histName.c_str()),histName.c_str(),800,600);
            
                 // Global style settings
                 gStyle->SetOptStat(1111);
                 gStyle->SetOptFit(111);
                 gStyle->SetLabelSize(0.03,"x");
                 gStyle->SetLabelSize(0.03,"y");
               
                 hist0->SetTitle(histName.c_str());
 
                 hist0->SetMaximum((hist0->GetMaximum())*1.1);
                 hist0->Draw("hist");
                 fitVAEnergy0->SetLineColor(kRed);
                 //hist0->GetFunction(fitVAEnergy0)->SetLineColour(kRed);
                 fitVAEnergy0->Draw("hist same"); 
                 
                 hist1->Draw("hist sames");
                 fitVAEnergy1->SetLineColor(kBlue);
                 //fitVAEnergy0->SetLineColor(kRed);
                 fitVAEnergy1->Draw("hist sames");
                 
                 gPad->Update(); 
                 //TPaveStats* sb2=(TPaveStats*)hist1->GetListOfFunctions()->FindObject("stats");
                 TPaveStats* sb0=(TPaveStats*)hist0->FindObject("stats");
                 //sb0->SetName(histName0.c_str());
                 sb0->SetX1NDC(.65);
                 sb0->SetX2NDC(.85);
                 sb0->SetY1NDC(.65);
                 sb0->SetY2NDC(.85);
                 sb0->SetTextColor(kRed);
                 TPaveStats* sb1=(TPaveStats*)hist1->FindObject("stats");
                 //sb1->SetName(histName1.c_str());
                 sb1->SetX1NDC(.65);
                 sb1->SetX2NDC(.85);
                 sb1->SetY1NDC(.4);
                 sb1->SetY2NDC(.6);
                 sb1->SetTextColor(kBlue);
                 //hist->Fit("gaus","V","same",0.,200.);
                 //hist->GetFunction("gaus")->SetLineColor(kBlue);
                 //hist->Draw("same");
                 //c1->Update();

                 //hist->Fit("landau");
                 //hist->GetFunction("landau")->SetLineColor(kGreen);
                 //hist->Draw("same");
                 //c1->Update();
                

                 /* ---- Fitting ends ---- */
            
                 //hist->Draw(); 
                 //hist->Write();
                 //fitVAEnergy->Write();
                 c1->Print();
                 //c1->Write();


                 /**----- Fitting end ------*/


                 //outFile->WriteTObject(hist);
                 outFile->WriteTObject(c1);
                 //outFile->WriteTObject(hMPV0);
                 //outFile->WriteTObject(hMPV1);
                 //delete hist;
                 delete c1;
                 //count++;
                 std::cout << "---------------------------------------------" << std::endl; 
            //}
         //}
     }
 
     //------------------------------------------------//
 
     outFile->WriteTObject(hMPV0);
     outFile->WriteTObject(hMPV1);
     outFile->cd();
     outFile->Write();
     outFile->Close();

//-----------------------------------------------------------//

   // std::ifstream histFile;
   // histFile.open("histolist.txt");
   //
   // std::string histName;
   // 
   // TFile *outFile = new TFile("../out/20181019/langaufit.root", "RECREATE");
   // 
   // int i = 0;
   // //while(!histFile.eof()){
   //         //std::cout << "debug . #" << i << std::endl;
   // while(std::getline(histFile,histName) && i < 10){ // for debug run
   // //while(std::getline(histFile,histName)){ // for analysis run

   // //...for(int iName = 0; iName < 1152; iName+=2){

   //     //float progress = 100.0 * ((float) i) / ((float) 2304);
   //     //if (!(i % 10)) cout << setprecision(3) << " [ " << progress << " % ] \r";
   //     //std::cout << std::setprecision(3) << "[ " << progress << " % ] \r";
   // 
   //     TH1D *hist = (TH1D*)inFile->Get(histName.c_str());
   //     //...TH1D *hist = (TH1D*)inFile->Get(histoNames.at(iName).c_str());
   // 
   //     if(hist->GetEntries() == 0) continue;
   //     /* ---- Fitting ---- */

   //     //TFile *f   = new TFile("../out/20181019/merged/merged_160919_104734.root");
   //     //TH1D *hVAEnergy = (TH1D*) f->Get("hVAEnergyX_71_0_0");
   // 
   //     std::cout << "Fitting histogram " << histName << std::endl;
   // 
   //     /*
   //      Once again, here are the Landau * Gaussian parameters:
   //       par[0]=Width (scale) parameter of Landau density
   //       par[1]=Most Probable (MP, location) parameter of Landau density
   //       par[2]=Total area (integral -inf to inf, normalization constant)
   //       par[3]=Width (sigma) of convoluted Gaussian function
   //    
   //     Variables for langaufit call:
   //       his             histogram to fit
   //       fitrange[2]     lo and hi boundaries of fit range
   //       startvalues[4]  reasonable start values for the fit
   //       parlimitslo[4]  lower parameter limits
   //       parlimitshi[4]  upper parameter limits
   //       fitparams[4]    returns the final fit parameters
   //       fiterrors[4]    returns the final fit errors
   //       ChiSqr          returns the chi square
   //       NDF             returns ndf
   //  
   //     */
   // 
   //     // Setting fit range and start values
   //     double fr[2], sv[4], pllo[4], plhi[4];
   //     fr[0] = 20.;         fr[1] = 200.;
   //     sv[0] = 2.;         sv[1] = 60.;        sv[2] = 50000.0;        sv[3] = 5.;
   //     pllo[0] = 0.5;      pllo[1] = 40.;      pllo[2] = 1.0;          pllo[3] = 1.;
   //     plhi[0] = 5.;       plhi[1] = 80.;      plhi[2] = 100000.0;     plhi[3] = 10.;
   // 
   //     // Return values
   //     double fp[4], fpe[4];
   //     double chisqr;
   //     int ndf;
   //     TF1 *fitVAEnergy = langaufit(hist,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
   // 
   //     double SNRPeak, SNRFWHM;
   //     langaupro(fp,SNRPeak,SNRFWHM);
   // 
   //     //std::cout << "Fitting done \n Plotting results..." << std::endl;
 
   //     if(i%4 == 0 || i%4 == 1){
   //         //MPV[i][0] = fp[1]; // MPVs for eta region 0
   //         hMPV0->Fill(fp[1]); // MPVs for eta region 0
   //     }
   //     if(i%4 == 2 || i%4 == 3){
   //         //MPV[i][1] = fp[1]; // MPVs for eta region 1
   //         hMPV1->Fill(fp[1]); // MPVs for eta region 0
   //     }

   //     TCanvas *c1 = new TCanvas("c1",histName.c_str(),800,600);
   //     //TCanvas can[i] = new TCanvas(Form(histName.c_str()),histName.c_str(),800,600);
 
   //     // Global style settings
   //     gStyle->SetOptStat(1111);
   //     gStyle->SetOptFit(111);
   //     gStyle->SetLabelSize(0.03,"x");
   //     gStyle->SetLabelSize(0.03,"y");
   // 
   //     hist->SetMaximum((hist->GetMaximum())*1.1);
   //     hist->Draw("hist");
   //     fitVAEnergy->Draw("hist same");

   //     /* ---- Fitting ends ---- */

   //     //hist->Draw(); 
   //     hist->Write();
   //     fitVAEnergy->Write();
   //     c1->Print();
   //     c1->Write();
   //     //outFile->Write();
   //     delete hist;
   //     delete fitVAEnergy;
   //     delete c1;
   //     i++;
   // }
   // 
// //   for(int j = 0; j < 1152; j++){
// //       hMPV0->Fill(MPV[j][0]);
// //       hMPV1->Fill(MPV[j][1]);
// //   }
   // hMPV0->Write();
   // hMPV1->Write();
   // std::cout << i << std::endl; 
   // outFile->cd();
   // outFile->Write();
   // outFile->Close();
   // histFile.close();
   // sw.Stop();
   // sw.Print();

    return 0;
} // end of main
