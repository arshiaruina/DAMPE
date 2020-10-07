#ifndef VAEQUALISATION_H
#define VAEQUALISATION_H 1
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

#define MAX_INCLINATION 60. // degrees
#define STEPS           12
#define INDEX           99

class VAequalisation {

    public:

        VAequalisation();
        ~VAequalisation();
        double CalcEta(DmpStkSiCluster* cluster);
        void InclInitialise();
        int CalcInclIndex(DmpStkTrack* track, std::string dir);
   
    private:

        double incl[STEPS];
};

#endif
