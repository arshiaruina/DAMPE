/* Code to draw a histogram of the BGO energy
 * Author: Arshia Ruina
 * Date: June 2019
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

#include "DmpEvtBgoRec.h"

int main(){

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

	return 0;
}
