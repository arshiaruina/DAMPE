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
//classes
//#include "landauGausConv.h"
//constants
#define nVA     6
#define nLADDER 192

class correctionFactors{

public:

    correctionFactors();
    ~correctionFactors();

    std::pair <double, double> MPVGausFit(std::string inFileNameMPVGausFit);
    std::string VAEnergyLangauFit();
    void Compute();
    
private:

    TH1D *hMPV0 = new TH1D("hMPV0","MPV of all VAs in RO-strip region",100,25.,85.);
    TH1D *hMPV1 = new TH1D("hMPV1","MPV of all VAs in floating-strip region",100,20.,60.);
    TH1D *hCorrFac0 = new TH1D("hCorrFac0","Correction factors in RO-strip region",100,0.,2.);
    TH1D *hCorrFac1 = new TH1D("hCorrFac1","Correction factors in floating-strip region",100,0,2.);
    TH1D *hCorrFacDiff = new TH1D("hCorrFacDiff","#Delta(Correction Factors)",100,-1.,1.);
    TH2D *hCorrFac = new TH2D("hCorrFac", "Correction factors", 192, 0, 191, 6, 0, 5);
    
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
    
    std::string computeDataPeriodStart = "20181001";
    std::string computeDataPeriodStop  = "20181009";
    std::string mergeTag = "231019_112429";
    std::string dir = "/beegfs/users/ruina/VAequalisation/out/" + computeDataPeriodStart + "_" + computeDataPeriodStop + "/";
    std::string inFileNameLangauFit  = dir + "merged/merged_" + mergeTag + ".root";
    std::string outFileNameLangauFit = dir + "test/langaufit_" + mergeTag + ".root";
    std::string outFileNameCorrFac   = dir + "test/corrFac_" + mergeTag + ".root";
};
//global constants
//#define N_VA        1152
//#define N_LADDER    192
//double corrFac0[N_LADDER][N_VA];
//double corrFac1[N_LADDER][N_VA];
