#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <string>
#include <stdio.h>
#include <sstream>
#include <math.h>
#include <cmath>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TChain.h"
#include "DmpEvtHeader.h"
#include "DmpEvtPsdRec.h"
#include "DmpRootEvent.h"
#include "DmpChain.h"
#include "DmpEvtGlobTrack.h"
#include "DmpSvcPsdEposCor.h"
#include "TClonesArray.h"
#include "TSystem.h"
#include "TROOT.h"

// #include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"

#include "DmpStkTrack.h"
#include "DmpStkTrackHelper.h"
#include "DmpStkSiCluster.h"
#include "DmpEvtPsdRec.h"
#include "TStopwatch.h"
#include "DmpStkClusterCalibration.h"

#include "track_selection.hpp"
#include "etacorr.hpp"
#include "analysis.hpp"

#include "definitions.hpp"
#include "helium_analysis.cpp" 
#include "simu_analysis.cpp"

using namespace std;
using namespace TMVA;


TH2D * hist_EvnEta(string t_name, string t_title) {
    TH2D * h = new TH2D(t_name.c_str(), t_title.c_str(), 50, 0., 1., 100, 0., 500.);
    h -> SetXTitle("#eta");
    h -> SetYTitle("STK #frac{dE}{dx}, [ADC]");
    h -> SetOption("colz");
    return h;
}


int main( int argc , char *argv[]){ 
    TStopwatch sw;
    sw.Start();

    // TMVA
    TMVA::Tools::Instance();
    map<string, int> Use;
    Use["BDT"] = 1; // uses Adaptive Boost

    int nToRun = 10; // set to negative to run all

    HeliumAnalysis * analysis = new HeliumAnalysis(argv[2], "RECREATE");
    analysis->setTChain(argv[1], true);
    analysis->addTTree();
 
    analysis->run(nToRun);
    delete analysis;

    SimuAnalysis * simu = new SimuAnalysis(argv[2], "UPDATE");
    simu->setTChain(argv[1], true);
    simu->openTTree("T");

    simu->run(nToRun);
    delete simu;
    return 0;
}
