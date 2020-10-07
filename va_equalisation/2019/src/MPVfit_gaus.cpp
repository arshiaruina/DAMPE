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

//int main(int argc, char** argv) {
int main() {

    //std::string inFileName = "../out/20181001_20181009/langaufit_231019_112429.root";
    std::string inFileName = "../out/20181001_20181009/langaufit_231019_112429_FMVAs.root";
    TFile *inFile = new TFile(inFileName.c_str());
    //TFile *outFile = new TFile("../out/20181001_20181009/MPV_langaufit_231019_112429.root", "RECREATE"); 
    TFile *outFile = new TFile("../out/20181001_20181009/gausfitMPV_20181001_20181009_FMVAs.root", "RECREATE"); 

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

    double eqParam0 = fGaus0->GetParameter(1);
    double eqParam1 = fGaus1->GetParameter(1);

    std::cout << "eqParam0: " << eqParam0 << " eqParam1: " << eqParam1 << std::endl;

    /**----- Fitting end ------*/

    outFile->Write();
    outFile->Close();

    delete c0;
    delete c1;

    return 0;
} // end of main

