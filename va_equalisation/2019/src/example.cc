/* Code to read EventHeader class from a root file 
 * From DAMPE SOftware documentation
 * Date: 09 July 2019
 * */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <sstream>
#include <cmath>
#include <vector>

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

#include "DmpEvtHeader.h"
#include "DmpEvtBgoRaw.h"
#include "DmpEvtBgoHits.h"
#include "DmpEvtBgoRec.h"

#include <filesystem>

using namespace std;
namespace fs = std::filesystem;

int main() {

	std::string path = "/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T232028_20181019T233706_00000_cu6fCNZZzBuq1VAgiuIw/";
	for (const auto & entry : fs::directory_iterator(path))
        std::cout << entry.path() << std::endl;

	gSystem->Load("./libDmpEvent.so");

	string inputfile = "/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T232028_20181019T233706_00000_cu6fCNZZzBuq1VAgiuIw/DAMPE_2A_OBS_20181019_20181019T232028_20181019T233706_00000_cu6fCNZZzBuq1VAgiuIw.root"; 
	//cout << "The accessed ROOT file is: " << inputfile << endl;
	
        //TChain *tc = new TChain("CollectionTree");
	//for(int i = 0; i < inputFiles.size(); i++) {
	//        tc->Add(inputFiles[i].c_str());
	//}

	TFile *f = new TFile(inputfile.c_str());
	TTree *t = (TTree*)f->Get("CollectionTree");

	DmpEvtHeader *fHeader = 0; 
	DmpEvtBgoRaw *fBgoRaw = 0;
	DmpEvtBgoHits *fBgoHits = 0;
	DmpEvtBgoRec *fBgoRec = 0;
	
	t->SetBranchAddress("EventHeader", &fHeader); 
	t->SetBranchAddress("DmpEvtBgoRaw", &fBgoRaw); 
	t->SetBranchAddress("DmpEvtBgoHits", &fBgoHits);
	t->SetBranchAddress("DmpEvtBgoRec", &fBgoRec);

	t->GetEntry(0);

	cout << "The timestamp of the first event in this ROOT file: ";
	cout << fHeader->GetSecond()<<endl; 

	cout << "The first ADC value of BGO of the first event in the ROOT file: ";
	cout << fBgoRaw->fADC.at(0) << endl;
	cout << fBgoRaw->fADC.size() << endl;

	cout << "The first dynode ID of BGO of the first event in the ROOT file: ";
	cout << fBgoRaw->fGlobalDynodeID.at(0) << endl;
	cout << fBgoRaw->fGlobalDynodeID.size() << endl;

	cout << "Energy [MeV] of the first hit of the first event in this ROOT file: ";
	cout << fBgoHits->fES0.at(0) << endl;
	cout << fBgoHits->fES0.size() << endl;

	cout << "Energy [MeV] deposited in the first layer of the first event in this ROOT file: ";
	cout << fBgoRec->GetELayer(0) << endl;

	return 0;
}

