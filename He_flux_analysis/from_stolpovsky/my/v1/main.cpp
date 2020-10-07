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

#include "definitions.hpp"
#include "reco_analysis.cpp"
#include "trig_analysis.cpp"

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

    int nToRun = -1; // set to negative to run all

    RecoAnalysis * reco = new RecoAnalysis(argv[2], "RECREATE");
    reco->setTChain(argv[1], true);
    reco->addTTree();

    reco->run(nToRun);
    delete reco;

    string filename(argv[2]);
    filename.insert(filename.find(".root"), "_trig");
    TrigAnalysis * trig = new TrigAnalysis(filename.c_str(), "RECREATE");
    trig->setTChain(argv[1], true);
    trig->addTTree();

    trig->run(nToRun);
    delete trig;
    return 0;
}
