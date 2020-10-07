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

#include "DmpStkTrack.h"
#include "DmpStkTrackHelper.h"
#include "DmpStkSiCluster.h"
#include "DmpEvtPsdRec.h"
#include "TStopwatch.h"
#include "DmpStkClusterCalibration.h"

#include "definitions.hpp"
#include "eta_corr.cpp"

using namespace std;


int main( int argc , char *argv[]){

    if (argc < 3) {
        cerr << "Usage: " << argv[0] << " <file_in.txt> <file_out.root> [Nev]" << endl;
	cerr << "\t <file_in.txt> -- name of file with the list of root files to analyse" << endl;
	cerr << "\t <file_out.root> -- name of the output root file" << endl;
	cerr << "\t [Nev] -- number of events to analyse (optional, useful for debug)" << endl;
	return 1;
    }

    char * filelist = argv[1];
    char * outfile = argv[2];
    int nToRun = argc==4 ? atoi(argv[3]) : -1; // negative = run all

    VaAnalysis * va_analysis = new VaAnalysis(outfile, "RECREATE");
    va_analysis->setTChain(filelist, true);
    va_analysis->addTTree();
    va_analysis->createHists();
    va_analysis->run(nToRun);

    delete va_analysis;
    return 0;
}
