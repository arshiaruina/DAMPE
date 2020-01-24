//
// Created by Arshia Ruina on 13.01.20.
//

#include "../inc/va_equalisation.h"
#include "../inc/clusterEnergyGen.h"

int main(int argc, char** argv) {
    TStopwatch sw;
    sw.Start();

    clusterEnergyGen clenge;
    va_equalisation vaeq;

    bool flagAppCorrFac = true; //false, if not applying VA calibration

    ////////    Open the required files     ////////

    std::string dataFileName = argv[1];
    TFile *dataFile = new TFile(dataFileName.c_str());
    TTree *dataTree = (TTree*)dataFile->Get("CollectionTree");

    TClonesArray* stktracks = new TClonesArray("DmpStkTrack"); //name of the class
    dataTree->SetBranchAddress("StkKalmanTracks",&stktracks); // name of the branch

    //DmpStkTrackHelper(TClonesArray* stktracks, bool usebgo = false, DmpEvtBgoRec* bgorec=0, DmpEvtBgoHits* bgohits=0):
    DmpStkTrackHelper* stkhelper = new DmpStkTrackHelper(stktracks,false,0,0);

    TClonesArray* stkclusters  = new TClonesArray("DmpStkSiCluster"); // name of the class
    dataTree->SetBranchAddress("StkClusterCollection",&stkclusters); // name of the branch

    std::size_t found = dataFileName.find_last_of("/");
    //std::string outFileName = "../out/20181001_20181009/" + dataFileName.substr(found+1);
    std::string outFileName;
    TFile *inFileCorrFac = new TFile(vaeq.inFileNameCorrFac.c_str());
    std::string hCorrFacName = "hCorrFac";
    TH2D *hCorrFac = (TH2D*)inFileCorrFac->Get(hCorrFacName.c_str());  
    if(flagAppCorrFac) {
        outFileName = vaeq.dir + "test/AppCorrFac/" + dataFileName.substr(found+1);

  }
    else {
        outFileName = vaeq.dir + "test/" + dataFileName.substr(found+1);
    }
    TFile *outFile = new TFile(outFileName.c_str(), "RECREATE");

    // Bad channels list  
    bool** badchannels =  clenge.read_bad_channels_file("../resources/bad_chan.txt");

    /* ==========================================

        Before submitting a job, check:

        -> input filename
        -> input file location
        -> output filename
        -> output file location
        -> number of entries looped over
        -> print statements

    =========================================== */

    ////////    Initialize the clenge.clusterEnergy histograms    ////////

    for(int iladder = 0; iladder < N_LADDER/2; iladder++) {
        if(iladder < 48) {
            clenge.xLadder = iladder+48; // X ladders 48-95
            clenge.yLadder = iladder;    // Y ladders 0-47
        }
        else {
            clenge.xLadder = iladder+96; // X ladders 144-191
            clenge.yLadder = iladder+48; // Y ladders 96-143
        }
        for(int iva = 0; iva < N_VA; iva++){
            for(int ietareg = 0; ietareg < 2; ietareg++){
                clenge.hVAEnergyX[iladder][iva][ietareg] = new TH1D(Form("clenge.hVAEnergyX_%d_%d_%d",clenge.xLadder,iva,ietareg),Form("Energy for ladder %d X VA %d #eta region %d",clenge.xLadder,iva,ietareg),200,0.,200.);
                clenge.hVAEnergyX[iladder][iva][ietareg] -> GetXaxis() -> SetTitle("Cluster energy");
                clenge.hVAEnergyX[iladder][iva][ietareg] -> GetYaxis() -> SetTitle("No. of events");
                clenge.histoNamesX.push_back("clenge.hVAEnergyX_" + std::to_string(clenge.xLadder) + "_" + std::to_string(iva) + "_" + std::to_string(ietareg));
                clenge.hVAEnergyY[iladder][iva][ietareg] = new TH1D(Form("clenge.hVAEnergyY_%d_%d_%d",clenge.yLadder,iva,ietareg),Form("Energy for ladder %d Y VA %d #eta region %d",clenge.yLadder,iva,ietareg),200,0.,200.);
                clenge.hVAEnergyY[iladder][iva][ietareg] -> GetXaxis() -> SetTitle("Cluster energy");
                clenge.hVAEnergyY[iladder][iva][ietareg] -> GetYaxis() -> SetTitle("No. of events");
                clenge.histoNamesY.push_back("clenge.hVAEnergyY_" + std::to_string(clenge.yLadder) + "_" + std::to_string(iva) + "_" + std::to_string(ietareg));
            }
        }
    }

    ////////    Start loop over all events    ////////

    int nEntries = dataTree->GetEntries();

    //for (int i = 0; i < nEntries; i++){ // uncomment for analysis run
     for (int i = 0; i < 100; i++){ // uncomment for debug run

        float progress = 100.0 * ((float) i) / ((float) nEntries);
        //if (!(i % 10)) cout << setprecision(3) << " [ " << progress << " % ] \r";
        std::cout << std::setprecision(3) << "[ " << progress << " % ] \r";

        dataTree->GetEntry(i);

        clenge.hNTracksNoCuts->Fill(stkhelper->GetSize());
        stkhelper->SortTracks(4,false); // sorting tracks, trackquality=1,bgomatchcut=false
        if(stkhelper->GetSize() == 0) continue;
        std::vector<int> vTrackIndex;

        ////////    Start loop over all tracks for selection    ////////

        for(int itrack = 0; itrack < stkhelper->GetSize(); itrack++) { // track loop to make selections

            DmpStkTrack *stktrack = (DmpStkTrack *) stkhelper->GetTrack(itrack);
            clenge.cosTheta = stktrack->getDirection().CosTheta();

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
                    } else {
                        stkcluster = stktrack->GetClusterY(ipoint, stkclusters);
                    }
                    if (!stkcluster) continue;
                    
                    // 1. remove clusters with noisy channels
                    if (clenge.is_cluster_bad_channel(stkcluster, badchannels)) continue;

                    clenge.clusterEnergy = stkcluster->getEnergy() * clenge.cosTheta;
                    //clenge.clusterEta = CalcEta(stkcluster);

                    if (ixy == 0 /*&& clenge.clusterEta != 0 && clenge.clusterEta != 1*/) {
                        nClusterX++;
                        sMeanX += clenge.clusterEnergy;
                        /*hEtaX -> Fill(clenge.clusterEta); std::cout << "cluster X present" << std::endl;*/
                    }
                    if (ixy == 1 /*&& clenge.clusterEta != 0 && clenge.clusterEta != 1*/) {
                        nClusterY++;
                        sMeanY += clenge.clusterEnergy;
                        /*hEtaY -> Fill(clenge.clusterEta); std::cout << "cluster Y present" << std::endl;*/
                    }
                } ////////    End of loop over x and y clusters    ////////

                // 2b. has 6 clusters
                if (nClusterX < 6 || nClusterY < 6) continue;

                // 3. |Z| = 1
                sMeanX = std::sqrt(sMeanX / nClusterX / 60.);
                sMeanY = std::sqrt(sMeanY / nClusterY / 60.);
                clenge.hSmeanX->Fill(sMeanX);
                clenge.hSmeanY->Fill(sMeanY);
                //if(!(sMeanX > 0.92 && sMeanX < 0.96 && sMeanY > 0.92 && sMeanY < 0.96)) continue;
                if (sMeanX < .9 || sMeanX > 1.1) continue;
                if (sMeanY < .9 || sMeanY > 1.1) continue;

                vTrackIndex.push_back(itrack);
                clenge.hNTracks->Fill(vTrackIndex.size());

            }   ////////    End of loop over points    ////////

        } ////////    End loop over all tracks for selection    ////////

        //if no tracks passed the selection criteria for this event, go to next event
        if (vTrackIndex.empty()) continue;

        ////////    Start loop over all selected tracks    ////////

        for(unsigned int itrack = 0; itrack < vTrackIndex.size(); itrack++){
            DmpStkTrack* stktrack = (DmpStkTrack*) stkhelper->GetTrack(itrack);
            clenge.cosTheta = stktrack->getDirection().CosTheta();

            ////////    Start of loop over points    ////////

            for (int ipoint = 0; ipoint < stktrack->GetNPoints(); ipoint++) {

                DmpStkSiCluster* stkcluster;

                ////////    Start loop over x and y clusters    ////////

                for(int ixy = 0; ixy < 2; ixy++){

                    if(ixy == 0){
                        stkcluster = stktrack -> GetClusterX(ipoint,stkclusters);
                    }
                    else{
                        stkcluster = stktrack -> GetClusterY(ipoint,stkclusters);
                    }
                    if(!stkcluster) continue;

                    clenge.clusterFirstStrip = stkcluster -> getFirstStrip();
                    clenge.clusterLastStrip = stkcluster -> getLastStrip();
                    clenge.clusterVA = clenge.GetVANumber(clenge.clusterFirstStrip, clenge.clusterLastStrip);

                    clenge.ladderNumber = stkcluster->getLadderHardware();

                    if(flagAppCorrFac){
                        //clenge.clusterEnergy = stkcluster->getEnergy() * clenge.cosTheta * corrFac[clenge.ladderNumber][clenge.clusterVA];
                        clenge.clusterEnergy = stkcluster->getEnergy() * clenge.cosTheta * hCorrFac->GetBinContent(clenge.ladderNumber+1,clenge.clusterVA+1);
                    }
                    else {
                        clenge.clusterEnergy = stkcluster->getEnergy()*clenge.cosTheta;
                    }
    
                    clenge.clusterEta = CalcEta(stkcluster);
                    clenge.clusterEtaReg = clenge.GetEtaRegion(clenge.clusterEta);

                    if(clenge.clusterVA < 0 || clenge.clusterEtaReg < 0) continue;
                    if(clenge.IsLadderX1(clenge.ladderNumber)) clenge.hVAEnergyX[clenge.ladderNumber-48][clenge.clusterVA][clenge.clusterEtaReg] -> Fill(clenge.clusterEnergy);
                    if(clenge.IsLadderX2(clenge.ladderNumber)) clenge.hVAEnergyX[clenge.ladderNumber-96][clenge.clusterVA][clenge.clusterEtaReg] -> Fill(clenge.clusterEnergy);
                    if(clenge.IsLadderY1(clenge.ladderNumber)) clenge.hVAEnergyY[clenge.ladderNumber][clenge.clusterVA][clenge.clusterEtaReg] -> Fill(clenge.clusterEnergy);
                    if(clenge.IsLadderY2(clenge.ladderNumber)) clenge.hVAEnergyY[clenge.ladderNumber-48][clenge.clusterVA][clenge.clusterEtaReg] -> Fill(clenge.clusterEnergy);

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
