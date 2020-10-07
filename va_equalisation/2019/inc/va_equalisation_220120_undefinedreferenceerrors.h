//
// Created by Arshia Ruina on 13.01.20.
//

#ifndef VA_EQUALISATION_H
#define VA_EQUALISATION_H

// C++
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
#include "TTree.h"
#include "TMath.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TRandom.h"
#include "TStyle.h"
#include "TClonesArray.h"
#include "TStopwatch.h"
#include "TSystem.h"
// DAMPE
#include "DmpEvtHeader.h"
#include "DmpEvtBgoRaw.h"
#include "DmpEvtBgoHits.h"
#include "DmpEvtBgoRec.h"
#include "DmpStkSiCluster.h"
#include "DmpRootEvent.h"
#include "DmpChain.h"
#include "DmpEvtGlobTrack.h"
#include "DmpStkTrack.h"
#include "DmpStkTrackHelper.h"
#include "DmpStkSiCluster.h"
#include "DmpEvtPsdRec.h"
#include "DmpStkClusterCalibration.h"

#define N_LADDER     192
#define N_VA         6
#define N_ETAREG     2
#define MINTRACKHITS 6

class va_equalisation {

public:

    va_equalisation();
    ~va_equalisation();

////declared these in clusterEnergyGen.h
//int GetEtaRegion(double);
//int GetVANumber(int,int);
//bool IsLadderX1(int);
//bool IsLadderX2(int);
//bool IsLadderY1(int);
//bool IsLadderY2(int);
//double CalcEta(DmpStkSiCluster*);

//private:

    std::string computeDataPeriodStart = "20181001";
    std::string computeDataPeriodStop  = "20181009";
    std::string mergeTag = "231019_112429";
    std::string dir = "/beegfs/users/ruina/VAequalisation/out/" + computeDataPeriodStart + "_" + computeDataPeriodStop + "/";
    std::string inFileNameCorrFac = dir + "test/corrFac_" + mergeTag + ".root";
    
};

#endif //VA_EQUALISATION_H
