/* Code to compute eta values
 * Author: Arshia Ruina
 * To compile, g++ -Wall `root-config --cflags --libs` -I$DMPSWSYS/include -L$DMPSWSYS/lib -lDmpEvent -lDmpService eta.cc -o ../submit/eta.exe
 * */

// example file
// /beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T001308_20181019T002610_00000_jEwJChikPnVApmeJU6og/DAMPE_2A_OBS_20181019_20181019T001308_20181019T002610_00000_jEwJChikPnVApmeJU6og.root

#include "va_equalisation.h"
//#include "../inc/track_selection.hpp"

using namespace std;


bool** read_bad_channels_file(const char* filename){
    bool** badchannels = new bool*[N_LADDER];
    for(int i=0; i<N_LADDER; i++){
        badchannels[i] = new bool[N_CHANNEL];
        for(int j=0; j<N_CHANNEL; j++){
            badchannels[i][j] = false;
        }
    }

    ifstream infile(filename);
    if (!infile)
    {
        std::cout << "[INFO] Can't open Bad Channels file: " << filename << " ==> throwing exception!" << std::endl;
        throw;
    }
    std::string line;
    while (getline(infile, line)) {
        if (line.c_str() == std::string("")) break;
        std::istringstream lineStream(line.c_str());
        std::string number;
        int i = 0;
        int ildr   = -1;
        int bdchnl = -1;
        while (getline(lineStream, number, ',')) {
            switch (i++) {
            case 0:
                ildr = atoi(number.c_str());
                break;
            case 1:
                bdchnl = atoi(number.c_str());
                break;
            }
        }
        if(ildr>=0 && bdchnl>=0 && ildr<N_LADDER && bdchnl<N_CHANNEL){
            badchannels[ildr][bdchnl] = true;
        }
    }
    std::cout << "[INFO] Done reading bad channels" << std::endl;
    return badchannels;
} // end read_bad_channels_file


bool is_cluster_bad_channel(DmpStkSiCluster * cluster, bool** badchannels) {
    int ladder = cluster->getLadderHardware(); //ladder num
    int minc = cluster->getIndex1();//address of the first strip in the cluster (0-363)
    int maxc = cluster->getIndex1() + cluster->getNstrip() -1;//last strip of the cluster

    for(int i = minc; i <= maxc; i++){
        if(badchannels[ladder][i]) return true;
    }
    return false;
} // end is_cluster_bad_channel



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

bool IsClusterAtVAEdge(int firstStrip, int lastStrip, int vaNumber){

    if(firstStrip == vaNumber*64 || lastStrip == (vaNumber+1)*64-1)
        return true;
    else
        return false;
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

    //int nClustersWithHotStrips = 0;

	//gSystem->Load("./libDmpEvent.so");

	std::string inFileName = argv[1];
	TFile *f = new TFile(inFileName.c_str(), "READ");
	TTree *t = (TTree*)f->Get("CollectionTree");

    // Bad channels list
    bool** badchannels =  read_bad_channels_file("/beegfs/users/ruina/VAequalisation/resources/bad_chan.txt");
	
    TClonesArray* stktracks = new TClonesArray("DmpStkTrack"); //name of the class
	t->SetBranchAddress("StkKalmanTracks",&stktracks); // name of the branch

	//DmpStkTrackHelper(TClonesArray* stktracks, bool usebgo = false, DmpEvtBgoRec* bgorec=0, DmpEvtBgoHits* bgohits=0):
	DmpStkTrackHelper* stkhelper = new DmpStkTrackHelper(stktracks,false,0,0);

    TClonesArray* stkclusters  = new TClonesArray("DmpStkSiCluster"); // name of the class
	t->SetBranchAddress("StkClusterCollection",&stkclusters); // name of the branch

    TClonesArray* stkladderadc = new TClonesArray("DmpStkLadderAdc"); // name of the class
	t->SetBranchAddress("DmpStkLadderAdcCollection",&stkladderadc); // name of the branch

    /* =====================================================
    
     Before submitting a job, check:
    
     -> input filename
     -> input file location
     -> output filename
     -> output file location
     -> number of entries looped over
     -> print statements
     
    ===================================================== */


	std::size_t found = inFileName.find_last_of("/");
    //std::string outFileName = "../out/20181001_20181009/" + inFileName.substr(found+1);
    //std::string outFileName = "../out/20181001_20181009_1/" + inFileName.substr(found+1);
    //std::string outFileName = "../out/20181019/" + inFileName.substr(found+1);
	//std::string outFileName = "../out/201810/" + inFileName.substr(found+1);
	//std::string outFileName = "../test/" + inFileName.substr(found+1);

    bool flagAppCorrFac = true;   //for corrFac application
    //bool flagAppCorrFac = false;    //without corrFac application

    //std::string outFileName;
    //TFile *inFileCorrFac = new TFile(inFileNameCorrFac.c_str());
    //std::string hCorrFacName = "hCorrFac";
    if(flagAppCorrFac) {
        inFileCorrFac = new TFile(inFileNameCorrFac.c_str(), "READ");
        hCorrFac = (TH2D*)inFileCorrFac->Get(hCorrFacName.c_str());
        outFileName = "out.root";
    }
    else {
        outFileName = dirA + "/unmerged/" + inFileName.substr(found+1);
        //outFileName = dirB + "/unmerged/" + inFileName.substr(found+1);
        //outFileName = "tmp/" + inFileName.substr(found+1);
    }
    TFile *outFile = new TFile(outFileName.c_str(), "RECREATE");

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
                histoNamesX.push_back("hVAEnergyX_" + std::to_string(xLadder) + "_" + std::to_string(iva) + "_" + std::to_string(ietareg));
                hVAEnergyY[iladder][iva][ietareg] = new TH1D(Form("hVAEnergyY_%d_%d_%d",yLadder,iva,ietareg),Form("Energy for ladder %d Y VA %d #eta region %d",yLadder,iva,ietareg),200,0.,200.);
                hVAEnergyY[iladder][iva][ietareg] -> GetXaxis() -> SetTitle("Cluster energy");
                hVAEnergyY[iladder][iva][ietareg] -> GetYaxis() -> SetTitle("No. of events");
                histoNamesY.push_back("hVAEnergyY_" + std::to_string(yLadder) + "_" + std::to_string(iva) + "_" + std::to_string(ietareg));
            }
        }
    }

    hEtaX = new TH1D("hEtaX","Eta distribution for X planes; #eta; No. of events",50,0.,1.);
    hEtaY = new TH1D("hEtaY","Eta distribution for Y planes; #eta; No. of events",50,0.,1.);
    //TH1F* hMeasCovXX 		= new TH1F("hMeasCovXX","MeasCovXX",100,-0.1,0.1);
    //TH1F* hMeasCovYY 		= new TH1F("hMeasCovYY","MeasCovYY",100,-0.1,0.1);
    //TH1F* hMeasCovXXSqrt 	= new TH1F("hMeasCovXXSqrt","Sqrt of MeasCovXX",100,-0.1,0.1);
    //TH1F* hMeasCovYYSqrt 	= new TH1F("hMeasCovYYSqrt","Sqrt of MeasCovYY",100,-0.1,0.1);
    //TH1F* hDistanceX 		= new TH1F("hDistanceX","|MeasHitX - FiltHitX|",100,0.,0.5);
    //TH1F* hDistanceY 		= new TH1F("hDistanceY","|MeasHitY - FiltHitY|",100,0.,0.5);
    //TH1I* hNoisyClusters 	= new TH1I("hNoisyClusters","Clusters with a noisy channel",1,0,1);    
    TH1I* hNTracksNoCuts 	= new TH1I("hNTracksNoCuts","No. of tracks with no selection cuts",10,0,10);    
    TH1I* hNTracks 			= new TH1I("hNTracks","No. of tracks that passed all cuts",10,0,10);    
    TH1D* hSmeanX 			= new TH1D("hSmeanX","hSmeanX",50,0.,5.);    
    TH1D* hSmeanY 			= new TH1D("hSmeanY","hSmeanY",50,0.,5.);    
    TH1D* hCheckEnergy 	    = new TH1D("hCheckEnergy","AdcValue - Energy",1000,-500.,500.);    
    
    std::stringstream oss;
    std::ofstream corrFacImplCheck;
    corrFacImplCheck.open("/beegfs/users/ruina/VAequalisation/tmp/corrFacImplCheck.txt");
    oss << "i |\t Ladder ID  |\t VA ID  |\t Adc Value  |\t Corr Fac  |\t Adc Val Corr \n";
    
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

    std::cout << "[INFO] variables set, starting loop over entries..." << std::endl;

	//for (int i = 0; i < nEntries; i++){ // uncomment for analysis run
    for (int i = 0; i < 100; i++){ // uncomment for debug run

		float progress = 100.0 * ((float) i) / ((float) nEntries);
        //if (!(i % 10)) cout << setprecision(3) << " [ " << progress << " % ] \r";
        std::cout << std::setprecision(3) << "[ " << progress << " % ] \r";

		t->GetEntry(i);

		hNTracksNoCuts->Fill(stkhelper->GetSize());
		stkhelper->SortTracks(4,false); // sorting tracks, trackquality=1,bgomatchcut=false		
        if(stkhelper->GetSize() == 0) continue;
        std::vector<int> vTrackIndex;

        //std::cout << "[DEBUG] track loop to make selections" << std::endl;

		//for(int itrack = 0; itrack <= stktracks->GetLast(); itrack++){
		for(int itrack = 0; itrack < stkhelper->GetSize(); itrack++){ // track loop to make selections

			//DmpStkTrack* stktrack = (DmpStkTrack*) stktracks->ConstructedAt(itrack);
			DmpStkTrack* stktrack = (DmpStkTrack*) stkhelper->GetTrack(itrack);
			/*double*/ cosTheta = stktrack->getDirection().CosTheta();

            int nClusterX = 0;
            int nClusterY = 0;
            double sMeanX = 0.;
            double sMeanY = 0.;

            // 2a. track made with 6 hits

            if(stktrack->getNhitXY() < MINTRACKHITS) continue;
    
			for (int ipoint = 0; ipoint < stktrack->GetNPoints(); ipoint++) {
	
				DmpStkSiCluster* stkcluster;

                for(int ixy = 0; ixy < 2; ixy++){

					if(ixy == 0){
						stkcluster = stktrack -> GetClusterX(ipoint,stkclusters);
                        //filling hDistanceX here without rejecting events with E < 20
                        //hDistanceX -> Fill(std::abs(stktrack->getHitMeasX(ipoint) - stktrack->getHitX(ipoint)));
				    }
                    else{
                        stkcluster = stktrack -> GetClusterY(ipoint,stkclusters);
                        //filling hDistanceY here without rejecting events with E < 20
                        //hDistanceY -> Fill(std::abs(stktrack->getHitMeasY(ipoint) - stktrack->getHitY(ipoint)));
                    }
            		if(!stkcluster) continue;
		
                    clusterFirstStrip = stkcluster -> getFirstStrip();
                    clusterLastStrip = stkcluster -> getLastStrip();
  					
                    // 1. remove clusters with noisy channels
                    if(is_cluster_bad_channel(stkcluster, badchannels)) continue;          
          
                    clusterEnergy = stkcluster->getEnergy()*cosTheta;
                    clusterEta = CalcEta(stkcluster);

                    if(ixy == 0 /*&& clusterEta != 0 && clusterEta != 1*/) {
                        nClusterX++;
                        sMeanX += clusterEnergy;
                        /*hEtaX -> Fill(clusterEta); std::cout << "cluster X present" << std::endl;*/
                    }
                    if(ixy == 1 /*&& clusterEta != 0 && clusterEta != 1*/) {
                        nClusterY++;
                        sMeanY += clusterEnergy;
                        /*hEtaY -> Fill(clusterEta); std::cout << "cluster Y present" << std::endl;*/}
	                
    			} // end of loop over x and y clusters
			} // end of loop over points
	    
			// 2b. has 6 clusters
			if(nClusterX < 6 || nClusterY < 6) continue;

			// 3. |Z| = 1
        	sMeanX = std::sqrt(sMeanX/nClusterX/60.);
        	sMeanY = std::sqrt(sMeanY/nClusterY/60.);

            hSmeanX->Fill(sMeanX);
            hSmeanY->Fill(sMeanY);

        	//if(!(sMeanX > 0.92 && sMeanX < 0.96 && sMeanY > 0.92 && sMeanY < 0.96)) continue;
        	if(sMeanX < .9 || sMeanX > 1.1) continue;
            if(sMeanY < .9 || sMeanY > 1.1) continue;

        	vTrackIndex.push_back(itrack);
        	hNTracks->Fill(vTrackIndex.size());


    	} // end of loop over tracks

        //std::cout << "[DEBUG] end of loop over tracks, selection" << std::endl;

    	if(vTrackIndex.empty()) continue;

        //std::cout << "[DEBUG] track loop " << std::endl;
        //std::cout << "[DEBUG] size of track vector " << vTrackIndex.size() << std::endl;

    	// loop over selected tracks
    	for(unsigned int itrack = 0; itrack < vTrackIndex.size(); itrack++){
    	//for(int itrack = 0; itrack < stkhelper->GetSize(); itrack++){

			//DmpStkTrack* stktrack = (DmpStkTrack*) stktracks->ConstructedAt(itrack);
			DmpStkTrack* stktrack = (DmpStkTrack*) stkhelper->GetTrack(itrack);
			/*double*/ cosTheta = stktrack->getDirection().CosTheta();

            //std::cout << "[DEBUG] points loop " << std::endl;
            
            //int nClusterX = 0;
            //int nClusterY = 0;
            //double sMeanX = 0.;
            //double sMeanY = 0.;

            // 2. Sel Cut
            //if(stktrack->getNhitXY() < MINTRACKHITS) continue;
    
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

                //hMeasCovXX->Fill(stktrack->getMeasCovXX(ipoint));
                //hMeasCovYY->Fill(stktrack->getMeasCovYY(ipoint));
                //hMeasCovXXSqrt->Fill(std::sqrt(stktrack->getMeasCovXX(ipoint)));
                //hMeasCovYYSqrt->Fill(std::sqrt(stktrack->getMeasCovYY(ipoint)));


                //std::cout << "[DEBUG] xy loop " << std::endl;
                
                for(int ixy = 0; ixy < 2; ixy++){

					if(ixy == 0){
						stkcluster = stktrack -> GetClusterX(ipoint,stkclusters);
                        //filling hDistanceX here without rejecting events with E < 20
                        //hDistanceX -> Fill(std::abs(stktrack->getHitMeasX(ipoint) - stktrack->getHitX(ipoint)));
				    }
                    else{
                        stkcluster = stktrack -> GetClusterY(ipoint,stkclusters);
                        //filling hDistanceY here without rejecting events with E < 20
                        //hDistanceY -> Fill(std::abs(stktrack->getHitMeasY(ipoint) - stktrack->getHitY(ipoint)));
                    }
            		if(!stkcluster) continue;
                                    
                    ladderNumber = stkcluster->getLadderHardware();
                    clusterFirstStrip = stkcluster -> getFirstStrip();
                    clusterLastStrip = stkcluster -> getLastStrip();
                    clusterVA = GetVANumber(clusterFirstStrip, clusterLastStrip);
                    clusterEta = CalcEta(stkcluster);
                    clusterEtaReg = GetEtaRegion(clusterEta);

                    if(clusterVA < 0 || clusterEtaReg < 0) continue;

                    clusterEnergyAdc = 0.;
                    for(int istrip = 0; istrip < stkcluster->getNstrip(); istrip++){
                        clusterEnergyAdc += stkcluster->GetAdcValue(istrip,stkladderadc);
                    }
                    
                    if(flagAppCorrFac){
                        oss << i << "\t" << ladderNumber << "\t" << clusterVA << "\t" << clusterEnergyAdc << "\t";
                        //clusterEnergy = stkcluster->getEnergy() * cosTheta * hCorrFac->GetBinContent(ladderNumber+1,clusterVA+1);
                        clusterEnergyAdc *= hCorrFac->GetBinContent(ladderNumber+1,clusterVA+1);
                        oss << hCorrFac->GetBinContent(ladderNumber+1,clusterVA+1) << "\t" << clusterEnergyAdc << "\n";
                    }
 
                    if(IsLadderX1(ladderNumber)) hVAEnergyX[ladderNumber-48][clusterVA][clusterEtaReg] -> Fill(clusterEnergyAdc); 
                    if(IsLadderX2(ladderNumber)) hVAEnergyX[ladderNumber-96][clusterVA][clusterEtaReg] -> Fill(clusterEnergyAdc); 
                    if(IsLadderY1(ladderNumber)) hVAEnergyY[ladderNumber][clusterVA][clusterEtaReg] -> Fill(clusterEnergyAdc); 
                    if(IsLadderY2(ladderNumber)) hVAEnergyY[ladderNumber-48][clusterVA][clusterEtaReg] -> Fill(clusterEnergyAdc); 
                   	//inclPerpIndex = CalcInclIndex(stktrack,"y");
					//inclPerpIndex = CalcInclIndex(stktrack,"x");

					//if(inclPerpIndex == 99) continue; //angle is >= 60deg
					//clusterEnergy = stkcluster->getEnergy()*cosTheta;	
					//hEtaEnergyVec.at(inclPerpIndex)->Fill(clusterEta, clusterEnergy);
					//hEtaEnergyVec.at(inclPerpIndex)->Draw("colz");
			
	                
    			} // end of loop over x and y clusters
			} // end of loop over points
	    
        	//sMeanX = std::sqrt(sMeanX/nClusterX/60.);
        	//sMeanY = std::sqrt(sMeanY/nClusterY/60.);

    	} // end of loop over tracks

        //std::cout << "[DEBUG] end of final loop over tracks" << std::endl;

	} // end of loop over entries

    corrFacImplCheck << oss.str();
    corrFacImplCheck.close();
    std::cout << "[INFO] Created tmp/corrFacImplCheck.txt" << std::endl;

    //std::cout << "[INFO] Creating overlaid histogram of energy distributions of all VAs" << std::endl;

    //TODO: ADC value distribution for all VAs together
/*    TCanvas* cAllVA = new TCanvas ("cAllVA","Energy distribution of all VAs",800,600);
    TH1D* hEmp = new TH1D ("hEmp","hEmp",200,0.,200.);
    hEmp->Draw("hist");
    cAllVA->Update();
    for(int iladder = 0; iladder < N_LADDER/2; iladder++) {
        for(int iva = 0; iva < N_VA; iva++){
            for(int ietareg = 0; ietareg < 2; ietareg++){
              
                if(IsLadderX1(iladder)) hVAEnergyX[iladder-48][iva][ietareg] -> Draw("hist same"); 
                if(IsLadderX2(iladder)) hVAEnergyX[iladder-96][iva][ietareg] -> Draw("hist same");
                if(IsLadderY1(iladder)) hVAEnergyY[iladder][iva][ietareg] -> Draw("hist same");
                if(IsLadderY2(iladder)) hVAEnergyY[iladder-48][iva][ietareg] -> Draw("hist same");

                cAllVA->Update();
            }
        }
    }
    
    cAllVA->Write();
*/
    
    //std::cout << "[DEBUG] end of loop over entries" << std::endl;

    //std::cout << "making checkfile" << std::endl;

    //stringstream ss;
    //std::ofstream checkfile;
    //checkfile.open("../out/20181019/check.txt");
    //ss << "making a check for VAs that have low energy events" << std::endl;
    //        std::cout << total[lad][va] << std::endl;
    //        std::cout << "doing things" << std::endl;
    //        frac = ( checkCounter[lad][va] / total[lad][va] ) * 100.;
    //        std::cout << "doing things" << std::endl;
    //        //if(frac > 1)
    //        //    ss << lad << "\t" << va << "\t" << checkCounter[lad][va] << "\t" << total[lad][va] << "\t" << frac << "  <<----" << std::endl;
    //        //else
    //        //    ss << lad << "\t" << va << "\t" << checkCounter[lad][va] << "\t" << total[lad][va] << "\t" << frac << std::endl;
    //    }
    //}
    //std::cout << "made it?" << std::endl;
    //checkfile << ss.str();
    //checkfile.close();


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
	if (outFile->IsZombie())
        std::cout << "[INFO] Error opening output file!" << std::endl;
    else
        std::cout << "[INFO] " << outFileName << " created." << std::endl;
	sw.Stop();
    sw.Print();
} // end of main

