/* Code to compute eta values
 * Author: Arshia Ruina
 * To compile, g++ -Wall `root-config --cflags --libs` -I$DMPSWSYS/include -L$DMPSWSYS/lib -lDmpEvent -lDmpService eta.cc -o ../submit/eta.exe
 * */

// example file
// /beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T001308_20181019T002610_00000_jEwJChikPnVApmeJU6og.root

#include <iostream>
#include <iomanip>
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
#include "TRandom.h"
#include "TStyle.h"

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

//#include "../inc/eta.h"

using namespace std;


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

        // finding highest signal
        for(int istrip = 0; istrip < nStrips; istrip++){
                sig = cluster->GetSignal(istrip);
                if(sig > sigMax1){
                        sigMax1 = sig;
                        stripMax1 = istrip;
                }
        }

        if(nStrips==1){
                eta = 0.;//--> then eta is set to 0. (biased)
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
                // Assuming that an adjoining strip must contain the second higest signal

                else {
                        if(cluster->GetSignal(stripMax1 - 1) > cluster->GetSignal(stripMax1 + 1)){
                                stripMax2 = stripMax1 - 1;
                        }
                        else {
                                stripMax2 = stripMax1 + 1;
                        }
                        sigMax2 = cluster->GetSignal(stripMax2);
                }
                
		// compute eta
                ch1 = cluster->GetChannelID(stripMax1);
                ch2 = cluster->GetChannelID(stripMax2);
                if(ch1 > ch2)
                        eta = sigMax1/(sigMax1 + sigMax2);
                else
                        eta = sigMax2/(sigMax1 + sigMax2);
        }
        return eta;
}

//double CalcIncl(DmpStkTrack* track, std::string dir){
int CalcInclIndex(DmpStkTrack* track, std::string dir){

	double maxIncl = 60.;
	int steps = 12;
	double incl[steps] = {0.};
	double theta = 0.;
	//double inclination = 99.;
	int index = 99; 

	for (int i = 0; i < steps; i++) {
		incl[i] = i * (maxIncl/steps);
	}

	if (dir == "x")
		theta = atan(track->getTrackParams().getSlopeX());
	else if (dir == "y")
		theta = atan(track->getTrackParams().getSlopeY());
	
	//double thetaX = atan(track->getTrackParams().getSlopeX());
	//double thetaY = atan(track->getTrackParams().getSlopeY());
	
	theta *= 180./TMath::Pi(); // in deg
	//thetaX *= 180./TMath::Pi(); // in deg
	//thetaY *= 180./TMath::Pi(); // in deg

	// use abs value of theta
	
	for(int i=0; i<12; i++){
		if(std::abs(theta) >= incl[i] && std::abs(theta) < incl[i]+5.) {
		//if(thetaX >= incl[i] && thetaX < incl[i+1] && thetaY >= incl[i] && thetaY < incl[i+1]) {
			//inclination = incl[i];
			index = i;
			break;
		}
	}
	//return inclination;
	return index;

}

int main(int argc, char** argv) {

	TStopwatch sw;
	sw.Start();

	//gSystem->Load("./libDmpEvent.so");

	std::string inFileName = argv[1];
	TFile *f = new TFile(inFileName.c_str());
	TTree *t = (TTree*)f->Get("CollectionTree");

	//DmpStkTrackHelper(TClonesArray* stktracks, bool usebgo = false, DmpEvtBgoRec* bgorec=0, DmpEvtBgoHits* bgohits=0):
	
	TClonesArray* stktracks = new TClonesArray("DmpStkTrack"); //name of the class
	t->SetBranchAddress("StkKalmanTracks",&stktracks); // name of the branch

	DmpStkTrackHelper* stkhelper = new DmpStkTrackHelper(stktracks,false,0,0);

    TClonesArray* stkclusters  = new TClonesArray("DmpStkSiCluster"); // name of the class
	t->SetBranchAddress("StkClusterCollection",&stkclusters); // name of the branch

	std::size_t found = inFileName.find_last_of("/");
	//std::string outFileName = "../out/20181019/" + inFileName.substr(found+1);
	std::string outFileName = "../out/201810/" + inFileName.substr(found+1);
	TFile *outFile = new TFile(outFileName.c_str(), "RECREATE");

    TH1D* hEtaX = new TH1D("hEtaX","Eta distribution for X planes; #eta; No. of events",50,0.,1.);
    TH1D* hEtaY = new TH1D("hEtaY","Eta distribution for Y planes; #eta; No. of events",50,0.,1.);
    //TH1D* hTrkChi2 = new TH1D("hTrkChi2","#chi^{2} ndof of tracks; #chi^{2}; No. of events",100,0.,50.);

	std::vector<TH2D*> hEtaEnergyVec;
	for(int ihist = 0; ihist < 12; ihist++){
		//std::string name = "hEtaEnergy_" + std::to_string(ihist);
		//std::string title = "Cluster charge distribution " + std::to_string(ihist*5) + "<|#theta_{x,y}|<" + std::to_string((ihist+1)*5);
		stringstream name, title;
		name << "hEtaEnergy_" << ihist;
		title << "Cluster charge distribution " << ihist*5 << "#leq|#theta_{x,y}|<" << (ihist+1)*5;
		TH2D* hist = new TH2D(name.str().c_str(), ";#eta; #frac{dE}{dx} [ADC counts]",50,0.,1.,50,0.,400.);
		hist->SetTitle(title.str().c_str());
		hist->SetOption("COLZ");
		hist->SetContour(30);
		hist->SetStats(0);
		hEtaEnergyVec.push_back(hist);
	}

	int nEntries = t->GetEntries();

	//for (int i = 0; i < nEntries; i++){ // uncomment for analysis run
	for (int i = 0; i < 100; i++){ // uncomment for debug run

		float progress = 100.0 * ((float) i) / ((float) nEntries);
        //if (!(i % 10)) cout << setprecision(3) << " [ " << progress << " % ] \r";
        std::cout << std::setprecision(3) << "[ " << progress << " % ] \r";

		t->GetEntry(i);

		stkhelper->SortTracks(4,false); // sorting tracks, trackquality=1,bgomatchcut=false
		if(stkhelper->GetSize() == 0) continue;

		//for(int itrack = 0; itrack <= stktracks->GetLast(); itrack++){
		for(int itrack = 0; itrack < stkhelper->GetSize(); itrack++){

			//DmpStkTrack* stktrack = (DmpStkTrack*) stktracks->ConstructedAt(itrack);
			DmpStkTrack* stktrack = (DmpStkTrack*) stkhelper->GetTrack(itrack);
			double cosTheta = stktrack->getDirection().CosTheta();
			//hTrkChi2 -> Fill(stktrack->getChi2NDOF());

			for (int ipoint = 0; ipoint < stktrack->GetNPoints(); ipoint++) {
	
				DmpStkSiCluster* stkcluster;
				double clusterEta = 0.;
				double clusterEnergy = 0.;			
				//double inclPerp = 0.;			
				int inclPerpIndex = 0;			

				for(int ixy = 0; ixy < 2; ixy++){
					if(ixy == 0){ // x cluster at this point
						stkcluster = stktrack -> GetClusterX(0,stkclusters);
						if(!stkcluster) continue;
						clusterEta = CalcEta(stkcluster);
						if(clusterEta == 0 || clusterEta == 1) continue;
						hEtaX -> Fill(clusterEta);
						inclPerpIndex = CalcInclIndex(stktrack,"y");
					}
					else{ // y cluster at this point
						stkcluster = stktrack -> GetClusterY(0,stkclusters);
						if(!stkcluster) continue;
						clusterEta = CalcEta(stkcluster);
						if(clusterEta == 0 || clusterEta == 1) continue;
						hEtaY -> Fill(clusterEta);
						inclPerpIndex = CalcInclIndex(stktrack,"x");
					}

					//if(clusterEta == 0. || clusterEta == 1.) continue;
					if(inclPerpIndex == 99) continue; //angle is >= 60deg
					clusterEnergy = stkcluster->getEnergy()*cosTheta;	
					hEtaEnergyVec.at(inclPerpIndex)->Fill(clusterEta, clusterEnergy);
					hEtaEnergyVec.at(inclPerpIndex)->Draw("colz");
				
				} // end of loop over clusters
			} // end of loop over points
		} // end of loop over tracks
	} // end of loop over entries

	//TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
	//..hEtaX->Draw();
	//..hEtaY->Write();
	//..for(int ihist = 0; ihist < 12; ihist++){
	//..	hEtaEnergyVec.at(ihist)->Draw();	
	//..	hEtaEnergyVec.at(ihist)->Write();	
	//..}
	//.. these lines were resulting in duplicate copies in the rootfile!
	//outFile->cd();
	outFile->Write();
	outFile->Close();
	std::cout << outFileName << " created." << std::endl;
	sw.Stop();
        sw.Print();
} // end of main

