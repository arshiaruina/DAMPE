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
#include "va_analysis.cpp"

using namespace std;


int main( int argc , char *argv[]){

    if (argc < 4) {
        cerr << "Usage: " << argv[0] << " <file_in.txt> <file_out.root> <0-1> <corr.root> [Nev]" << endl;
	cerr << "\t <file_in.txt> -- name of file with the list of root files to analyse" << endl;
	cerr << "\t <file_out.root> -- name of the output root file" << endl;
	cerr << "\t <corr.root> -- provide the name of file with correction coeffs, "
             << "or 0, if no correction" << endl;
	cerr << "\t [Nev] -- number of events to analyse (optional, useful for debug)" << endl;
	return 1;
    }

    char * filelist = argv[1];
    char * outfile = argv[2];
    char * corr = argv[3];
    bool calc_corr = strcmp(corr, "0") == 0;

    int nToRun = argc==5? atoi(argv[4]) : -1; // negative = run all

    VaAnalysis * va_analysis = new VaAnalysis(outfile, "RECREATE");
    va_analysis->setTChain(filelist, true);
    va_analysis->addTTree();
    va_analysis->createHists();
    if (!calc_corr) va_analysis->setCorr(corr);

    va_analysis->run(nToRun);

    delete va_analysis;
    return 0;
}
