/* Code to compute eta values
 * Author: Arshia Ruina
 * To compile, g++ -Wall `root-config --cflags --libs` -I$DMPSWSYS/include -L$DMPSWSYS/lib -lDmpEvent -lDmpService eta.cc -o ../submit/eta.exe
 * */

// example file
// /beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T001308_20181019T002610_00000_jEwJChikPnVApmeJU6og/DAMPE_2A_OBS_20181019_20181019T001308_20181019T002610_00000_jEwJChikPnVApmeJU6og.root

#include "../inc/va_equalisation.h"

using namespace std;

int GetEtaRegion(double eta){
    
    int etaReg = -99;
    if(eta < 0.2 || eta > 0.8)
        etaReg = 0;
    if(eta > 0.4 && eta < 0.6)
        etaReg = 1;
    return etaReg;
}

int GetVANumber(int firstStrip, int lastStrip){

    int vaNumber = -99;
    for(int iva = 0; iva < 6; iva++) {
        if(firstStrip >= iva*64 && lastStrip < (iva+1)*64)
            vaNumber = iva;
    }
    //if(firstStrip >= 0 &&  lastStrip < 64)  vaNumber = 0; 
    //if(firstStrip >= 64 &&  lastStrip < 128) vaNumber = 1;  
    //if(firstStrip >= 128 && lastStrip < 192) vaNumber = 2;
    //if(firstStrip >= 192 && lastStrip < 256) vaNumber = 3;
    //if(firstStrip >= 256 && lastStrip < 320) vaNumber = 4;
    //if(firstStrip >= 320 && lastStrip < 384) vaNumber = 5;
    return vaNumber;
}

// iladder          0 - 47          48 - 95
// X ladder IDs     48 - 95         144 - 191
// Y ladder IDs     0 - 47          96 - 143

bool IsLadderX1(int ladder){
    if(ladder >= 48 && ladder < 96)
        return true;
    else
        return false;
}

bool IsLadderX2(int ladder){
    if(ladder >= 144 && ladder < 192)
        return true;
    else
        return false;
}

bool IsLadderY1(int ladder){
    if(ladder >= 0 && ladder < 48)
        return true;
    else
        return false;
}

bool IsLadderY2(int ladder){
    if(ladder >= 96 && ladder < 144)
        return true;
    else
        return false;
}

double CalcEta(DmpStkSiCluster* cluster){

    double sig = 0.;
    double sigMax1 = 0.;
    double sigMax2 = 0.;
    int nStrips = cluster->getNstrip();
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

    if(nStrips==1){ // if cluster in one strip only
        eta = 0.;// then eta is set to 0. (biased)
    }
    else {
            // finding second highest signal
        if(stripMax1==0){ // if highest signal strip is first strip of the cluster
            stripMax2 = 1;
            sigMax2 = cluster->GetSignal(stripMax2);
        }
        else if(stripMax1==nStrips-1){ // if highest signal strip is the last strip of the cluster
            stripMax2 = stripMax1 - 1;
            sigMax2 = cluster->GetSignal(stripMax2);
        }
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



int main(int argc, char** argv) {

	TStopwatch sw;
	sw.Start();

	//gSystem->Load("./libDmpEvent.so");

	std::string inFileName = argv[1];
	TFile *f = new TFile(inFileName.c_str());
	TTree *t = (TTree*)f->Get("CollectionTree");
	
    TClonesArray* stktracks = new TClonesArray("DmpStkTrack"); //name of the class
	t->SetBranchAddress("StkKalmanTracks",&stktracks); // name of the branch

	//DmpStkTrackHelper(TClonesArray* stktracks, bool usebgo = false, DmpEvtBgoRec* bgorec=0, DmpEvtBgoHits* bgohits=0):
	DmpStkTrackHelper* stkhelper = new DmpStkTrackHelper(stktracks,false,0,0);

    TClonesArray* stkclusters  = new TClonesArray("DmpStkSiCluster"); // name of the class
	t->SetBranchAddress("StkClusterCollection",&stkclusters); // name of the branch

	std::size_t found = inFileName.find_last_of("/");
	//std::string outFileName = "../out/20181019/" + inFileName.substr(found+1);
	//std::string outFileName = "../out/201810/" + inFileName.substr(found+1);
	std::string outFileName = "test/" + inFileName.substr(found+1);
	TFile *outFile = new TFile(outFileName.c_str(), "RECREATE");

    hEtaX = new TH1D("hEtaX","Eta distribution for X planes; #eta; No. of events",50,0.,1.);
    hEtaY = new TH1D("hEtaY","Eta distribution for Y planes; #eta; No. of events",50,0.,1.);
    
    for(int iladder = 0; iladder < N_LADDER/2; iladder++) { 
        if(iladder < 48) {
            xLadder = iladder+48; // X ladders 48-95 
            yLadder = iladder;    // Y ladders 0-47
        }
        else {
            xLadder = iladder+96; // X ladders 144-191
            yLadder = iladder+48; // Y ladders 96-143
        }
        for(int iva = 0; iva < N_VA; iva++){
            for(int ietareg = 0; ietareg < 2; ietareg++){
                hVAEnergyX[iladder][iva][ietareg] = new TH1D(Form("hVAEnergyX_%d_%d_%d",xLadder,iva,ietareg),Form("Energy for ladder %d X VA %d #eta region %d",xLadder,iva,ietareg),200,0.,200.);
                hVAEnergyX[iladder][iva][ietareg] -> GetXaxis() -> SetTitle("Cluster energy");
                hVAEnergyX[iladder][iva][ietareg] -> GetYaxis() -> SetTitle("No. of events");
                hVAEnergyY[iladder][iva][ietareg] = new TH1D(Form("hVAEnergyY_%d_%d_%d",yLadder,iva,ietareg),Form("Energy for ladder %d Y VA %d #eta region %d",yLadder,iva,ietareg),200,0.,200.);
                hVAEnergyY[iladder][iva][ietareg] -> GetXaxis() -> SetTitle("Cluster energy");
                hVAEnergyY[iladder][iva][ietareg] -> GetYaxis() -> SetTitle("No. of events");
            }
        }
    }
    
	//.. std::vector<TH2D*> hEtaEnergyVec;
	//.. for(int ihist = 0; ihist < 12; ihist++){
	//.. 	//std::string name = "hEtaEnergy_" + std::to_string(ihist);
	//.. 	//std::string title = "Cluster charge distribution " + std::to_string(ihist*5) + "<|#theta_{x,y}|<" + std::to_string((ihist+1)*5);
	//.. 	stringstream name, title;
	//.. 	name << "hEtaEnergy_" << ihist;
	//.. 	title << "Cluster charge distribution " << ihist*5 << "#leq|#theta_{x,y}|<" << (ihist+1)*5;
	//.. 	TH2D* hist = new TH2D(name.str().c_str(), ";#eta; #frac{dE}{dx} [ADC counts]",50,0.,1.,50,0.,400.);
	//.. 	hist->SetTitle(title.str().c_str());
	//.. 	hist->SetOption("COLZ");
	//.. 	hist->SetContour(30);
	//.. 	hist->SetStats(0);
	//.. 	hEtaEnergyVec.push_back(hist);
	//.. }

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
			/*double*/ cosTheta = stktrack->getDirection().CosTheta();

			for (int ipoint = 0; ipoint < stktrack->GetNPoints(); ipoint++) {
	
				DmpStkSiCluster* stkcluster;
				//double clusterEta = 0.;
				//double clusterEnergy = 0.;			
				////double inclPerp = 0.;			
				////int inclPerpIndex = 0;			
                //int ladderNumber = 0;
                //int clusterFirstStrip = 0;
                //int clusterLastStrip = 0;
                //int clusterVA = 99;

                for(int ixy = 0; ixy < 2; ixy++){

					if(ixy == 0){
						stkcluster = stktrack -> GetClusterX(0,stkclusters);
				    }
                    else{
                        stkcluster = stktrack -> GetClusterY(0,stkclusters);
                    }
            		if(!stkcluster) continue;
					
                    clusterEnergy = stkcluster -> getEnergy()*cosTheta;
                    clusterEta = CalcEta(stkcluster);

                    if(ixy == 0 && clusterEta != 0 && clusterEta != 1) {hEtaX -> Fill(clusterEta); std::cout << "cluster X present" << std::endl;}
                    if(ixy == 1 && clusterEta != 0 && clusterEta != 1) {hEtaY -> Fill(clusterEta); std::cout << "cluster Y present" << std::endl;}
                    
                    ladderNumber = stkcluster -> getLadderHardware();
                    clusterFirstStrip = stkcluster -> getFirstStrip();
                    clusterLastStrip = stkcluster -> getLastStrip();
                    clusterVA = GetVANumber(clusterFirstStrip, clusterLastStrip);
                    clusterEtaReg = GetEtaRegion(clusterEta);

                    std::cout << "cluster ladder: " << ladderNumber << std::endl;
                    std::cout << "cluster VA: " << clusterVA << std::endl;

                    if(clusterVA < 0 || clusterEtaReg < 0) continue;
                    if(IsLadderX1(ladderNumber)) hVAEnergyX[ladderNumber-48][clusterVA][clusterEtaReg] -> Fill(clusterEnergy); 
                    if(IsLadderX2(ladderNumber)) hVAEnergyX[ladderNumber-96][clusterVA][clusterEtaReg] -> Fill(clusterEnergy); 
                    if(IsLadderY1(ladderNumber)) hVAEnergyY[ladderNumber][clusterVA][clusterEtaReg] -> Fill(clusterEnergy); 
                    if(IsLadderY2(ladderNumber)) hVAEnergyY[ladderNumber-48][clusterVA][clusterEtaReg] -> Fill(clusterEnergy); 
                   	//inclPerpIndex = CalcInclIndex(stktrack,"y");
					//inclPerpIndex = CalcInclIndex(stktrack,"x");

					//if(inclPerpIndex == 99) continue; //angle is >= 60deg
					//clusterEnergy = stkcluster->getEnergy()*cosTheta;	
					//hEtaEnergyVec.at(inclPerpIndex)->Fill(clusterEta, clusterEnergy);
					//hEtaEnergyVec.at(inclPerpIndex)->Draw("colz");
				
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

