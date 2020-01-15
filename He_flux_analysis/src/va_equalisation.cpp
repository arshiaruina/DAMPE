//
// Created by Arshia Ruina on 13.01.20.
//

#include "va_equalisation.h"
#include "clusterEnergyGen.h"

int main(int argc, char** argv) {
    TStopwatch sw;
    sw.Start();


    ////////    Initialise all variables    ////////


    clusterEnergyGen clenge;

    TH1D* hVAEnergyX[N_LADDER/2][N_VA][N_ETAREG];
    TH1D* hVAEnergyY[N_LADDER/2][N_VA][N_ETAREG];
    TH1D* hEtaX = new TH1D("hEtaX","Eta distribution for X planes; #eta; No. of events",50,0.,1.);
    TH1D* hEtaY = new TH1D("hEtaY","Eta distribution for Y planes; #eta; No. of events",50,0.,1.);
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
    int checkCounter[N_LADDER][N_VA] = {0};
    int total[N_LADDER][N_VA] = {0};
    //storage for names of the histograms
    std::vector<std::string> histoNamesX, histoNamesY;
    //std::stringstream ssX, ssY;
    int xLadder             = -99;
    int yLadder             = -99;
    double cosTheta         = -99.;
    double clusterEta       = -99.;
    double clusterEnergy    = -99.;
    //double inclPerp;
    // int inclPerpIndex;
    int ladderNumber        = -99;
    int clusterFirstStrip   = -99;
    int clusterLastStrip    = -99;
    int clusterVA           = -99;
    int clusterEtaReg       = -99;
    bool flagAppCorrFac     = true; //false, if not applying VA calibration


    ////////    Open the required files     ////////


    // Bad channels list
    bool** badchannels =  clenge.read_bad_channels_file("../resources/bad_chan.txt");

    std::string dataFileName = argv[1];
    TFile *dataFile = new TFile(dataFileName.c_str());
    TTree *dataTree = (TTree*)dataFile->Get("CollectionTree");

    TClonesArray* stktracks = new TClonesArray("DmpStkTrack"); //name of the class
    dataTree->SetBranchAddress("StkKalmanTracks",&stktracks); // name of the branch

    //DmpStkTrackHelper(TClonesArray* stktracks, bool usebgo = false, DmpEvtBgoRec* bgorec=0, DmpEvtBgoHits* bgohits=0):
    DmpStkTrackHelper* stkhelper = new DmpStkTrackHelper(stktracks,false,0,0);

    TClonesArray* stkclusters  = new TClonesArray("DmpStkSiCluster"); // name of the class
    dataTree->SetBranchAddress("StkClusterCollection",&stkclusters); // name of the branch

    /* =====================================================

    Before submitting a job, check:

    -> input filename
    -> input file location
    -> output filename
    -> output file location
    -> number of entries looped over
    -> print statements

    ===================================================== */

    std::size_t found = dataFileName.find_last_of("/");
    //std::string outFileName = "../out/20181001_20181009/" + dataFileName.substr(found+1);
    std::string outFileName = "../test/" + dataFileName.substr(found+1);
    TFile *outFile = new TFile(outFileName.c_str(), "RECREATE");


    ////////    Prepare the clusterEnergy histograms    ////////


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

    int nEntries = dataTree->GetEntries();


    ////////    Start loop over all events    ////////


    for (int i = 0; i < nEntries; i++){ // uncomment for analysis run
    // for (int i = 0; i < 100; i++){ // uncomment for debug run

        float progress = 100.0 * ((float) i) / ((float) nEntries);
        //if (!(i % 10)) cout << setprecision(3) << " [ " << progress << " % ] \r";
        std::cout << std::setprecision(3) << "[ " << progress << " % ] \r";

        dataTree->GetEntry(i);

        hNTracksNoCuts->Fill(stkhelper->GetSize());
        stkhelper->SortTracks(4,false); // sorting tracks, trackquality=1,bgomatchcut=false
        if(stkhelper->GetSize() == 0) continue;
        std::vector<int> vTrackIndex;

        ////////    Start loop over all tracks for selection    ////////

        //for(int itrack = 0; itrack <= stktracks->GetLast(); itrack++){
        for(int itrack = 0; itrack < stkhelper->GetSize(); itrack++) { // track loop to make selections

            //DmpStkTrack* stktrack = (DmpStkTrack*) stktracks->ConstructedAt(itrack);
            DmpStkTrack *stktrack = (DmpStkTrack *) stkhelper->GetTrack(itrack);
            cosTheta = stktrack->getDirection().CosTheta();

            int nClusterX = 0;
            int nClusterY = 0;
            double sMeanX = 0.;
            double sMeanY = 0.;

            // 2a. track made with 6 hits

            if (stktrack->getNhitXY() < MINTRACKHITS) continue;

            ////////    Start of loop over points    ////////

            for (int ipoint = 0; ipoint < stktrack->GetNPoints(); ipoint++) {
                DmpStkSiCluster *stkcluster;

                ////////    Start of loop over x and y clusters    ////////

                for (int ixy = 0; ixy < 2; ixy++) {
                    if (ixy == 0) {
                        stkcluster = stktrack->GetClusterX(ipoint, stkclusters);
                        //filling hDistanceX here without rejecting events with E < 20
                        //hDistanceX -> Fill(std::abs(stktrack->getHitMeasX(ipoint) - stktrack->getHitX(ipoint)));
                    } else {
                        stkcluster = stktrack->GetClusterY(ipoint, stkclusters);
                        //filling hDistanceY here without rejecting events with E < 20
                        //hDistanceY -> Fill(std::abs(stktrack->getHitMeasY(ipoint) - stktrack->getHitY(ipoint)));
                    }
                    if (!stkcluster) continue;
                    
                    // 1. remove clusters with noisy channels
                    if (is_cluster_bad_channel(stkcluster, badchannels)) continue;

                    clusterEnergy = stkcluster->getEnergy() * cosTheta;
                    //clusterEta = CalcEta(stkcluster);

                    if (ixy == 0 /*&& clusterEta != 0 && clusterEta != 1*/) {
                        nClusterX++;
                        sMeanX += clusterEnergy;
                        /*hEtaX -> Fill(clusterEta); std::cout << "cluster X present" << std::endl;*/
                    }
                    if (ixy == 1 /*&& clusterEta != 0 && clusterEta != 1*/) {
                        nClusterY++;
                        sMeanY += clusterEnergy;
                        /*hEtaY -> Fill(clusterEta); std::cout << "cluster Y present" << std::endl;*/
                    }
                } ////////    End of loop over x and y clusters    ////////

                // 2b. has 6 clusters
                if (nClusterX < 6 || nClusterY < 6) continue;

                // 3. |Z| = 1
                sMeanX = std::sqrt(sMeanX / nClusterX / 60.);
                sMeanY = std::sqrt(sMeanY / nClusterY / 60.);
                hSmeanX->Fill(sMeanX);
                hSmeanY->Fill(sMeanY);
                //if(!(sMeanX > 0.92 && sMeanX < 0.96 && sMeanY > 0.92 && sMeanY < 0.96)) continue;
                if (sMeanX < .9 || sMeanX > 1.1) continue;
                if (sMeanY < .9 || sMeanY > 1.1) continue;

                vTrackIndex.push_back(itrack);
                hNTracks->Fill(vTrackIndex.size());

            }   ////////    End of loop over points    ////////

        } ////////    End loop over all tracks for selection    ////////

        //if no tracks passed the selection criteria for this event, go to next event
        if (vTrackIndex.empty()) continue;

        ////////    Start loop over all selected tracks    ////////

        //for(int itrack = 0; itrack < stkhelper->GetSize(); itrack++){
        for(unsigned int itrack = 0; itrack < vTrackIndex.size(); itrack++){
            //DmpStkTrack* stktrack = (DmpStkTrack*) stktracks->ConstructedAt(itrack);
            DmpStkTrack* stktrack = (DmpStkTrack*) stkhelper->GetTrack(itrack);
            cosTheta = stktrack->getDirection().CosTheta();

            ////////    Start of loop over points    ////////

            for (int ipoint = 0; ipoint < stktrack->GetNPoints(); ipoint++) {

                DmpStkSiCluster* stkcluster;

                ////////    Start loop over x and y clusters    ////////

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
                    clusterVA = GetVANumber(clusterFirstStrip, clusterLastStrip);

                    ladderNumber = stkcluster->getLadderHardware();

                    if(flagAppCorrFac){
                        clusterEnergy = stkcluster->getEnergy() * cosTheta * corrFac[ladderNumber][clusterVA];
                    }
                    else {
                        clusterEnergy = stkcluster->getEnergy()*cosTheta;
                    }
    
                    clusterEta = CalcEta(stkcluster);
                    clusterEtaReg = GetEtaRegion(clusterEta);

                    if(clusterVA < 0 || clusterEtaReg < 0) continue;
                    if(IsLadderX1(ladderNumber)) hVAEnergyX[ladderNumber-48][clusterVA][clusterEtaReg] -> Fill(clusterEnergy);
                    if(IsLadderX2(ladderNumber)) hVAEnergyX[ladderNumber-96][clusterVA][clusterEtaReg] -> Fill(clusterEnergy);
                    if(IsLadderY1(ladderNumber)) hVAEnergyY[ladderNumber][clusterVA][clusterEtaReg] -> Fill(clusterEnergy);
                    if(IsLadderY2(ladderNumber)) hVAEnergyY[ladderNumber-48][clusterVA][clusterEtaReg] -> Fill(clusterEnergy);

                } ////    End of loop over x and y clusters   ////////

            } ////////    End of loop over points    ////////

        } ////////    End loop over all selected tracks    ////////

    } ////////    End loop over all events    ////////

    outFile->Write();
    outFile->Close();
    if (outFile->IsZombie())
        std::cout << "Error opening output file!" << std::endl;
    else
        std::cout << outFileName << " created." << std::endl;
    sw.Stop();
    sw.Print();

    return 0;
    } // end of main
