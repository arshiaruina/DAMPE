#ifndef CORRFAC_H
#define CORRFAC_H 1

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
#include "TH1.h"
#include "TH2.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TRandom.h"
#include "TStyle.h"
#include "TPaveStats.h"
//constants
#define nVA     6
#define nLADDER 192

class corrFac{

public:

    //corrFac();
    //~corrFac();

    std::pair <double, double> MPVGausFit(std::string inFileNameMPVGausFit);
    std::string VAEnergyLangauFit();
    void Compute(bool);
    
    std::string inFileNameLangauFit;   
    std::string outFileNameLangauFit;  
    std::string outFileNameMPVGausFit;
    std::string outFileNameCorrFac;

private:

    TH1D *hMPV0 = new TH1D("hMPV0","MPV of all VAs in RO-strip region",100,25.,85.);
    TH1D *hMPV1 = new TH1D("hMPV1","MPV of all VAs in floating-strip region",100,20.,60.);
    TH1D *hCorrFac0 = new TH1D("hCorrFac0","Correction factors in RO-strip region",100,0.,2.);
    TH1D *hCorrFac1 = new TH1D("hCorrFac1","Correction factors in floating-strip region",100,0,2.);
    TH1D *hCorrFacDiff = new TH1D("hCorrFacDiff","#Delta(Correction Factors)",100,-1.,1.);
    TH2D *hCorrFac = new TH2D("hCorrFac", "Correction factors", 192, 0, 192, 6, 0, 6);
    
    int countVA = 0;
    int iva; // 0-5
    int iladder; // 0-191
    //const int nVA = 6;
    //const int nLADDER = 192;
    //TCanvas *c1[50][2];
    //double MPV0[50][2];
    //double MPV1[50][2];
    //double corrFac0[50][2];
    //double corrFac1[50][2];
    TCanvas *c1[nLADDER][nVA];
    double MPV0[nLADDER][nVA];
    double MPV1[nLADDER][nVA];
    double corrFac0[nLADDER][nVA];
    double corrFac1[nLADDER][nVA]; // TODO declare outside main
    std::pair<double, double> eqParams;
    std::vector<std::string> EQMladders = {"12","13","14","15","36","37","38","39","66","67","68"};
    
    std::ifstream histFile;
    std::string histName;
    std::string histFileName = "/beegfs/users/ruina/VAequalisation/resources/histolist_new.txt";
    


    // some things need to be changed.... 11.02.2020
    
    //std::string startPeriodCompute  = "20181001";
    //std::string stopPeriodCompute   = "20181009";
    //std::string startPeriodApply    = "20181001"; /* = "20181101"; */
    //std::string stopPeriodApply     = "20181009"; /* = "20181109"; */

    //std::string dirComputeNotCorrected  = "/beegfs/users/ruina/VAequalisation/out/periodCompute/" + startPeriodCompute + "_" + stopPeriodCompute + "/not_corrected";
    //std::string dirComputeCorrected     = "/beegfs/users/ruina/VAequalisation/out/periodCompute/" + startPeriodCompute + "_" + stopPeriodCompute + "/corrected";
    //std::string dirApplyNotCorrected    = "/beegfs/users/ruina/VAequalisation/out/periodApply/" + startPeriodApply + "_" + stopPeriodApply + "/not_corrected";
    //std::string dirApplyCorrected       = "/beegfs/users/ruina/VAequalisation/out/periodApply/" + startPeriodApply + "_" + stopPeriodApply + "/corrected";

    //std::string mergeTag = "030220_174755"; //dirComputeNotCorrected 
    //std::string mergeTag = ""; // dirComputeCorrected
    //std::string mergeTag = "030220_112552"; //dirApplyNotCorrected 
    //std::string mergeTag = "050220_160930"; //dirApplyCorrected




    /////////////       To make langaufits, MPVgausfits and compute corrFac...
    
    //std::string startPeriodA = "20181001";
    //std::string stopPeriodA  = "20181009";
    //std::string dirA = "/beegfs/users/ruina/VAequalisation/out/periodA/data_selection_cuts/" + startPeriodA + "_" + stopPeriodA + "/";
    //std::string dirA = "/beegfs/users/ruina/VAequalisation/out/periodA/data_corrected/" + startPeriodA + "_" + stopPeriodA + "/";
    //
    ////std::string mergeTag = "030220_174755";
    ////std::string inFileNameLangauFit     = dirA + "merged/" + mergeTag + ".root";
    ////std::string outFileNameLangauFit    = dirA + "langaufit/" + mergeTag + ".root";
    ////std::string outFileNameMPVGausFit   = dirA + "gausfitMPV/" + mergeTag + ".root";
    ////std::string outFileNameCorrFac      = dirA + "corrFac/" + mergeTag + ".root";
    //
    //////////////        To make langaufits and MPVgausfits only...
    //
    //std::string startPeriodB = "20181101";
    //std::string stopPeriodB  = "20181109";
    //std::string dirB1 = "/beegfs/users/ruina/VAequalisation/out/periodB/data_selection_cuts/" + startPeriodB + "_" + stopPeriodB + "/";
    //std::string dirB2 = "/beegfs/users/ruina/VAequalisation/out/periodB/data_corrected/" + startPeriodB + "_" + stopPeriodB + "/";

    //////////////        ... for /data_selection_cuts
    //    
    ////std::string mergeTag = "030220_112552";
    ////std::string inFileNameLangauFit  = dirB1 + "merged/" + mergeTag + ".root";
    ////std::string outFileNameLangauFit = dirB1 + "langaufit/" + mergeTag + ".root";
    ////std::string outFileNameMPVGausFit = dirB1 + "gausfitMPV/" + mergeTag + ".root";
    ////std::string outFileNameCorrFac;
    //
    //////////////        ... for /data_corrected
    //    
    //std::string mergeTag = "050220_160930";
    //std::string inFileNameLangauFit  = dirB2 + "merged/" + mergeTag + ".root";
    //std::string outFileNameLangauFit = dirB2 + "langaufit/" + mergeTag + ".root";
    //std::string outFileNameMPVGausFit = dirB2 + "gausfitMPV/" + mergeTag + ".root";
    //std::string outFileNameCorrFac;
    
};
//global constants
//#define N_VA        1152
//#define N_LADDER    192
//double corrFac0[N_LADDER][N_VA];
//double corrFac1[N_LADDER][N_VA];

#endif
