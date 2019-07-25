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
#include "TSystem.h"

//#include "track_selection.hpp"
//#include "etacorr.hpp"

//#include <filesystem>

using namespace std;
//namespace fs = std::filesystem;

int main() {

        TStopwatch sw;
        sw.Start();

        gSystem->Load("./libDmpEvent.so");

        std::ifstream listOfFiles;
        listOfFiles.open("../resources/20181019.txt"); // one day data

        if(!listOfFiles.good()){ // Returns true if none of the stream's error state flags (eofbit, failbit and badbit) is set.
                std::cout << "Problem with list of files!" << std::endl;
                return 1;
        }

        int nFiles = 0;

        TH1D* hEta = new TH1D("hEta","Eta distribution; #eta; No. of events",50,0.,1.);
        TH12D* hEtaEnergy = new TH2D("hEtaEnergy","Cluster charge distribution; #eta; #frac{dE}{dx} [ADC counts]",50,0.,1.,50,0,200);

        while(!listOfFiles.eof()) { // uncomment for analysis run
        //while(!listOfFiles.eof() && nFiles<1) { // uncomment for debug run

                std::string fileName;
                std::getline(listOfFiles, fileName);

                TFile *f = new TFile(fileName.c_str());
                std::cout << "Reading file " << fileName.c_str() << std::endl;

                TTree *t = (TTree*)f->Get("CollectionTree");

                TClonesArray* stktracks = new TClonesArray("DmpStkTrack"); //name of the class
                t->SetBranchAddress("StkKalmanTracks",&stktracks); // name of the branch
                //DmpStkTrackHelper(TClonesArray* stktracks, bool usebgo = false, DmpEvtBgoRec* bgorec=0, DmpEvtBgoHits* bgohits=0):

                TClonesArray* stkclusters  = new TClonesArray("DmpStkSiCluster"); // name of the class
                t->SetBranchAddress("StkClusterCollection",&stkclusters); // name of the branch

                std::cout << std::endl;
                std::cout << "No. of events " << t->GetEntries() << std::endl;

                for (int i = 0; i < t->GetEntries(); i++){ // uncomment for analysis run
                //for (int i = 0; i < 100; i++){ // uncomment for debug run

                        std::cout << std::endl;
                        std::cout << "Processing event " << i+1 << std::endl;
                        t->GetEntry(i);

                        std::cout << std::endl;
                        std::cout << "No. of tracks in this event " << stktracks->GetLast()+1 << std::endl;

                        for(int itrack = 0; itrack <= stktracks->GetLast(); itrack++){

                                std::cout << std::endl;
                                std::cout << "Processing track " << itrack+1 << " of event " << i+1 << std::endl;
                                DmpStkTrack* stktrack = (DmpStkTrack*) stktracks->ConstructedAt(itrack);

                                std::cout << std::endl;
                                std::cout << "No. of clusters in this track " << stktrack->GetNPoints() << std::endl;

                                std::cout << "stktrack->getNhitXY() " << stktrack->getNhitXY() << std::endl;
                                std::cout << "stktrack->getNhitX() " << stktrack->getNhitX() << std::endl;
                                std::cout << "stktrack->getNhitY() " << stktrack->getNhitY() << std::endl;
                                std::cout << "stktrack->getImpactPoint() " << "[" << stktrack->getImpactPoint().X() << "," << stktrack->getImpactPoint().Y() << "," << stktrack->getImpactPoint().Z() << "]" << std::endl;
                                std::cout << "stktrack->getImpactPointHasX() " << stktrack->getImpactPointHasX() << std::endl;
                                std::cout << "stktrack->getImpactPointHasY() " << stktrack->getImpactPointHasY() << std::endl;
                                std::cout << "stktrack->isImpactPointLayerX() " << stktrack->isImpactPointLayerX() << std::endl;
                                std::cout << "stktrack->isImpactPointLayerY() " << stktrack->isImpactPointLayerY() << std::endl;

                                double clusterEnergy = 0.;
                                double clusterEta = 0.;
                                double cosTheta = stktrack->getDirection().CosTheta();

                                for (int icluster = 0; icluster < stktrack->GetNPoints(); icluster++) {

                                        std::cout << std::endl;
                                        std::cout << "Processing cluster " << icluster+1 << " of track " << itrack+1 << " of event " << i+1 << std::endl;

                                        DmpStkSiCluster* stkcluster = (DmpStkSiCluster*) stkclusters->ConstructedAt(icluster);
                                        
					std::cout << "debug 1" << std::endl;
                                        clusterEnergy = stkcluster->getEnergy()*cosTheta;
                                        
					std::cout << "debug 2" << std::endl;
                                        if(stkcluster->getNstrip()==0) continue;
                                        clusterEta = CalcEta(stkcluster);
                                        
					std::cout << "debug 3" << std::endl;
                                        hEta -> Fill(clusterEta);
                                        //... hEtaEnergy -> Fill(clusterEta, clusterEnergy);
                                        //... hEtaenergy->Draw("colz");

                                } // end of loop over clusters
                        } // end of loop over tracks
                } // end of loop over events
                nFiles++;
        } // end of loop over files     

        TCanvas c1("c1", "c1", 800, 600);
        hEta -> Draw();
        c1.SaveAs("../out/eta.pdf");

        //... TCanvas c2("c2", "c2", 800, 600);
     	//... hEtaEnergy -> Draw();
        //... c2.SaveAs("../out/etaenergy.pdf");

        sw.Stop();
        sw.Print();
} // end of main


double CalcEta(DmpStkSiCluster* cluster){

        // Eta has to be calculated for each cluster associated with the selected track

        //std::cout << "CalcEta::debug 1" << std::endl; 
        // initialization
        double sig = 0.;
        double sigMax1 = 0.;
        double sigMax2 = 0.;
        int nStrips = cluster->getNstrip(); // number of strips
        int stripMax1 = 0;
        int stripMax2 = 0;
        int ch1 = 0;
        int ch2 = 0;
        double eta = 0.;

        //std::cout << "CalcEta::debug 2" << std::endl; 
        // finding highest signal
        for(int istrip = 0; istrip < nStrips; istrip++){
                sig = cluster->GetSignal(istrip);
                if(sig > sigMax1){
                        sigMax1 = sig;
                        stripMax1 = istrip;
                }
        }

        //std::cout << "CalcEta::debug 3" << std::endl; 
        std::cout << "Size of cluster " << nStrips << std::endl;
        std::cout << "Strip with highest signal " << stripMax1 << std::endl;
        //std::cout << cluster->GetSignal(stripMax1-1) << std::endl; 
        std::cout << "Highest signal " << cluster->GetSignal(stripMax1) << std::endl;
        //std::cout << cluster->GetSignal(stripMax1+1) << std::endl; 
        // finding adjacent strip with second-highest signal
        if(nStrips==1){
                eta = 0.;//--> then eta is set to 0. (biased)
                std::cout << "debug me 1" << std::endl;
        }
        else {
                if(stripMax1==0){ // highest signal strip is first strip of the cluster
                        stripMax2 = 1;
                        sigMax2 = cluster->GetSignal(stripMax2);
                }
                else if(stripMax1==nStrips-1){ // highest signal strip is the last strip of the cluster
                        stripMax2 = stripMax1 - 1;
                        sigMax2 = cluster->GetSignal(stripMax2);
                }
                // so here I have assumed that the next neighbouring strip must contain the second higest signal

                //std::cout << "debug me 2" << std::endl;
                //if(nStrips==2) {
                //      std::cout << "debug me 3" << std::endl;
                //      stripMax2 = !stripMax1;
                //      sigMax2 = cluster->GetSignal(stripMax2);
                //}
                else {
                        std::cout << "debug me 4" << std::endl;
                        if(cluster->GetSignal(stripMax1 - 1) > cluster->GetSignal(stripMax1 + 1)){
                                stripMax2 = stripMax1 - 1;
                                std::cout << "debug me 5" << std::endl;
                        }
                        else {
                                stripMax2 = stripMax1 + 1;
                                std::cout << "debug me 6" << std::endl;
                        }
                        std::cout << "debug me 7" << std::endl;
                        sigMax2 = cluster->GetSignal(stripMax2);
                }
                //std::cout << "CalcEta::debug 4" << std::endl; 
                std::cout << "Strip with second highest signal " << stripMax2 << std::endl;
                std::cout << "Second highest signal " << cluster->GetSignal(stripMax2) << std::endl;
                // compute eta
                ch1 = cluster->GetChannelID(stripMax1);
                ch2 = cluster->GetChannelID(stripMax2);
                std::cout << "Channel ID of strip with highest signal " << ch1 << std::endl;
                std::cout << "Channel ID of strip with second highest signal " << ch2 << std::endl;
                if(ch1 > ch2)
                        eta = sigMax1/(sigMax1 + sigMax2);
                else
                        eta = sigMax2/(sigMax1 + sigMax2);
        }
        std::cout << "Eta for this cluster " << eta << std::endl;
        //std::cout << "CalcEta::debug 5" << std::endl; 
        return eta;
}

double CalcIncl(DmpStkTrack* track, std::string dir){

	double maxIncl = 60.;
	int steps = 12;
	double incl[steps] = {0.};
	double theta = 0.;
	double inclination = 99.;
 
	for (int i = 1; i < steps; i++) {
		incl[i] = incl[i-1] + (maxAngle/steps);
	}

	if (dir == "x")
		double theta = atan(track->getTrackParams().getSlopeX());
	else if (dir == "y")
		double theta = atan(track->getTrackParams().getSlopeY());
	
	theta *= 180./TMath::Pi(); // in deg
	
	for(int i=0; i<11; i++){
		if(theta >= incl[i] && theta < incl[i+1]) {
			inclination = incl[i];
			break;
		}
	}
	return inclination;

}

/*
void PlotEta(){

        // Input files

        std::string inputFile1 = "/beegfs/dampe/prod/MC/reco/v6r0p0/allProton-v6r0p0_10TeV_100TeV-HP-p6/allProton-v6r0p0_10TeV_100TeV-HP-p6.noOrb.607416.reco.root";
        std::string inputFile2 = "/beegfs/dampe/prod/MC/reco/v6r0p0/allProton-v6r0p0_10TeV_100TeV-HP-p6/allProton-v6r0p0_10TeV_100TeV-HP-p6.noOrb.601082.reco.root";
        //TFile *ntuple1 = new TFile(inputFile1.c_str(), "READ");
        //TFile *ntuple2 = new TFile(inputFile2.c_str(), "READ");

        TChain *tc = new TChain("CollectionTree");
        tc->Add(inputFile1.c_str());
        tc->Add(inputFile2.c_str());

        // Get entrie
        int nEvents = tc->GetEntries();

        // Total energy of event from BGO calorimeter 
        std::vector<double> energyVec;
        TH1D *energyHist = new TH1D("energyHist", "allProton-v6r0p0_10TeV_100TeV-HP-p6; Reco energy; No. of events",22,10e5,10e7);

        DmpEvtBgoRec *bgoRec = new DmpEvtBgoRec();
        tc->SetBranchAddress("DmpEvtBgoRec",&bgoRec);


        // Loop over all events and store values


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
*/
