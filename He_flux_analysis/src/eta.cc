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
#include "DmpStkSiCluster.h"

#include "DmpRootEvent.h"
#include "DmpChain.h"
#include "DmpEvtGlobTrack.h"
#include "TClonesArray.h"

#include "DmpStkTrack.h"
#include "DmpStkTrackHelper.h"
#include "DmpStkSiCluster.h"
#include "DmpEvtPsdRec.h"
#include "TStopwatch.h"
#include "DmpStkClusterCalibration.h"

//#include "track_selection.hpp"
//#include "etacorr.hpp"

//#include <filesystem>

using namespace std;
//namespace fs = std::filesystem;

int main() {

	gSystem->Load("./libDmpEvent.so");

	string inputfile = "/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T232028_20181019T233706_00000_cu6fCNZZzBuq1VAgiuIw/DAMPE_2A_OBS_20181019_20181019T232028_20181019T233706_00000_cu6fCNZZzBuq1VAgiuIw.root"; 

	TFile *f = new TFile(inputfile.c_str());
	TTree *t = (TTree*)f->Get("CollectionTree");

	//DmpStkSiCluster *fCluster = 0;
	//t->SetBranchAddress("DmpStkSiCluster", &fCluster);

TClonesArray* stkclusters  = new TClonesArray("DmpStkSiCluster");
t->SetBranchAddress("StkClusterCollection",stkclusters); // name o fthe branch

TClonesArray* stktracks = new TClonesArray("DmpStkTrack"); //name of the class
t->SetBranchAddress("StkKalmanTracks",stktracks); // name o fthe branch
//DmpStkTrackHelper(TClonesArray* stktracks, bool usebgo = false, DmpEvtBgoRec* bgorec=0, DmpEvtBgoHits* bgohits=0):



	int nEvents = t->GetEntries();
	for (int i = 0; i < nEvents; i++) {

		t->GetEntry(i);
		for(int itrack=0; itrack<=stktracks->GetLast();++itrack){
			DmpStkTrack* stktrack = (DmpStkTrack*) stktracks->ConstructedAt(itrack);
			for (int ipoint=0; ipoint<stktrack->GetNPoints(); ipoint++) {
			}
		}


            // Loop over clusters
            for (int ipoint=0; ipoint<stktrack->GetNPoints(); ipoint++) {
                for (int ixy=0; ixy<2; ixy<++) {
                    // Check if cluster i<s associated to a hit
                    DmpStkSiCluster* cluster;
                    if(ixy == 0)
                        cluster = stktrack -> GetClusterX(ipoint, stkclusters);
                    else
                        cluster = stktrack -> GetClusterY(ipoint, stkclusters);
                    if(!cluster) continue;


                    double e = cluster -> getEnergy() * costheta;
                    double eta = CalcEta(cluster);
                }
            }
        }


	return 0;
}


double CalcEta(DmpStkSiCluster* cluster){
  
	// Eta has to be calculated for each cluster associated with the selected track

	// initialization
	double sig = 0.;
	double sigMax1 = 0.;
	double sigMax2 = 0.;
	int nstrips = cluster->getNstrip(); // number of strips
	int stripMax1 = 0;
	int stripMax2 = 0;
	int ch1 = 0;
	int ch2 = 0;
	double eta = 0.;
	
	// finding highest signal
	for(int istrip = 0; istrip < nstrips; istrip++){
		sig = cluster->GetSignal(istrip);
		if(sig > sigMax1){
			sigMax1 = sig;
			stripMax1 = istrip;
		}
		ch1++;
	}
	
	// finding adjacent strip with second-highest signal
	if(cluster->GetSignal(stripMax1 - 1) > cluster->GetSignal(stripMax1 + 1)
		stripMax2 = stripMax1 - 1;
	else
		stripMax2 = stripMax1 + 1;
	
	sigMax2 = cluster->GetSignal(stripMax2);
	
	// compute eta
	ch1 = GetChannelID(stripMax1);
	ch2 = GetChannelID(stripMax2);
	if(ch1 < ch2)
		eta = sigMax1/(sigMax1 + sigMax2);
	else
		eta = sigMax2/(sigMax1 + sigMax2);
	
	return eta;
}

void PlotEta(){

        /* Input files */

        std::string inputFile1 = "/beegfs/dampe/prod/MC/reco/v6r0p0/allProton-v6r0p0_10TeV_100TeV-HP-p6/allProton-v6r0p0_10TeV_100TeV-HP-p6.noOrb.607416.reco.root";
        std::string inputFile2 = "/beegfs/dampe/prod/MC/reco/v6r0p0/allProton-v6r0p0_10TeV_100TeV-HP-p6/allProton-v6r0p0_10TeV_100TeV-HP-p6.noOrb.601082.reco.root";
        //TFile *ntuple1 = new TFile(inputFile1.c_str(), "READ");
        //TFile *ntuple2 = new TFile(inputFile2.c_str(), "READ");

        TChain *tc = new TChain("CollectionTree");
        tc->Add(inputFile1.c_str());
        tc->Add(inputFile2.c_str());

        /* Get entries */
        int nEvents = tc->GetEntries();

        /* Total energy of event from BGO calorimeter */
        std::vector<double> energyVec;
        TH1D *energyHist = new TH1D("energyHist", "allProton-v6r0p0_10TeV_100TeV-HP-p6; Reco energy; No. of events",22,10e5,10e7);

        DmpEvtBgoRec *bgoRec = new DmpEvtBgoRec();
        tc->SetBranchAddress("DmpEvtBgoRec",&bgoRec);


        /* Loop over all events and store values */


	for(int i = 0; i < nEvents; i++){
		tc->GetEntry(i);
		energyHist->Fill(bgoRec->GetTotalEnergy());
	}       
	
	TCanvas c("c", "c", 800, 600);
	c.SetLogy();
	c.SetLogx();
	energyHist->Draw();
	c.SaveAs("../out/energy_histogram.png");

}
