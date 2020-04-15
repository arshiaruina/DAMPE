// Date: 05.03.2020
// Source: https://www.nevis.columbia.edu/~seligman/root-class/files/MakeNtuple.C

// example file
// /beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T001308_20181019T002610_00000_jEwJChikPnVApmeJU6og/DAMPE_2A_OBS_20181019_20181019T001308_20181    019T002610_00000_jEwJChikPnVApmeJU6og.root

#include "../inc/ntupleMaker.h"

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
        //eta = 0.;// then eta is set to 0. (biased)
        eta = 1.;
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
    return vaNumber;
}

struct TRACK
{
    int trackIndex;
    int nClustersX;
    int nClustersY;
    int nClustersTotal;
};

int main(int argc, char** argv){
   
    std::string inputFileName = argv[1];
    TFile *inputFile = new TFile(inputFileName.c_str());
    if(inputFile == 0){
        std::cout << "[ERROR] File " << inputFileName << " not found!" << std::endl;
        return 0;
    }

    DmpChain *dmpch = new DmpChain("CollectionTree");
    dmpch->Add(inputFileName.c_str());
    dmpch->GetListOfFiles()->Print();
    int nEntries = dmpch->GetEntries();

    std::size_t found = inputFileName.find_last_of("/");
    std::string outputFileName = "/beegfs/users/ruina/VAequalisation/testNtuples/" + inputFileName.substr(found+1);
    //std::string outputFileName = "/beegfs/users/ruina/VAequalisation/ntuples_20181001_20181009/" + inputFileName.substr(found+1);
    //std::string outputFileName = "/beegfs/users/ruina/VAequalisation/ntuples_201801_201802/" + inputFileName.substr(found+1);
    TFile* outputFile = new TFile(outputFileName.c_str(),"recreate");
    TTree* outputTree = new TTree("MySelectionTree","Selections for VA calibration");

    int nClustersX;
    int nClustersY;
    int nClustersTotal;
    float clusterEnergy [MINTRACKCLUS] = {0.};
    float clusterEnergyAdc [MINTRACKCLUS] = {0.};
    float clusterEta [MINTRACKCLUS] = {0.};
    int clusterEtaRegion [MINTRACKCLUS] = {0};
    int clusterLadder [MINTRACKCLUS] = {0};
    int clusterVA [MINTRACKCLUS] = {0}; 

    outputTree->Branch("nClustersX", &nClustersX, "nClustersX/I");
    outputTree->Branch("nClustersY", &nClustersY, "nClustersY/I");
    outputTree->Branch("nClustersTotal", &nClustersTotal, "nClustersTotal/I");
    outputTree->Branch("clusterEnergy", clusterEnergy, "clusterEnergy[nClustersTotal]/F");
    outputTree->Branch("clusterEnergyAdc", clusterEnergyAdc, "clusterEnergyAdc[nClustersTotal]/F");
    outputTree->Branch("clusterEta", clusterEta, "clusterEta[nClustersTotal]/F");
    outputTree->Branch("clusterEtaRegion", clusterEtaRegion, "clusterEtaRegion[nClustersTotal]/I");
    outputTree->Branch("clusterLadder", clusterLadder, "clusterLadder[nClustersTotal]/I");
    outputTree->Branch("clusterVA", clusterVA, "clusterVA[nClustersTotal]/I");

    int nSel = 7;
    int nEntriesAfterSelection[nSel];
    int nEntriesAfterSelection1track[nSel];
    TH1I* hNEntriesAfterSelection = new TH1I("hNEntriesAfterSelection","Selections",nSel+1,0,nSel);
    TH1I* hNEntriesAfterSelection1track = new TH1I("hNEntriesAfterSelection1track","Selections (1 track events)",nSel+1,0,nSel);
    //TH1I* hNTracksNoCuts    = new TH1I("hNTracksNoCuts","No. of tracks with no selection cuts",10,0,10);    
    //TH1I* hNTracks          = new TH1I("hNTracks","No. of tracks that passed all cuts",10,0,10);    
    TH1D* hSmeanX           = new TH1D("hSmeanX","hSmeanX",50,0.,5.);    
    TH1D* hSmeanY           = new TH1D("hSmeanY","hSmeanY",50,0.,5.);    
    TH1F* hEnergy1stripCluster = new TH1F("hEnergy1stripCluster","Energy for 1-strip clusters",200,0.,200.);    

    for(int i = 0; i < nSel; i++){
        nEntriesAfterSelection[i] = 0;
        nEntriesAfterSelection1track[i] = 0;
    }
    
    //////////////////////////////////////////////////
    //             Loop over events
    //////////////////////////////////////////////////
    for(int ientry = 0; ientry < nEntries; ++ientry)
    //for(int ientry = 0; ientry < 20; ++ientry)
    {   

        float progress = 100.0 * ((float) ientry) / ((float) nEntries);
        std::cout << std::setprecision(3) << "[ " << progress << " % ] \r";

        std::vector<TRACK> vSelectedTracks;

        DmpEvent *pev = dmpch->GetDmpEvent();    
        int nTracks[nSel];
        for(int i = 0; i < nSel; i++){
            nTracks[i] = 0;
        }
        //nTracks[0] = pev->NStkKalmanTrack();
        int nHitsOnTrack = 0;
        //////////////////////////////////////////////////
        //              Loop over tracks
        //////////////////////////////////////////////////
        for(int itrack = 0; itrack < pev->NStkKalmanTrack(); itrack++){
        
            DmpStkTrack* stktrack = pev->pStkKalmanTrack(itrack);
            if(pev->pStkKalmanTrack(itrack))
                nTracks[0]++;
        
            float sMeanX = 0.;
            float sMeanY = 0.;
            int nClusX = 0.;
            int nClusY = 0.;
            int nholes = 0;
            int nonoverlaps = stktrack->getNhitX() + stktrack->getNhitY() - 2 * stktrack->getNhitXY();
            double chisqndof = stktrack->getChi2NDOF();
            float cosTheta = stktrack->getDirection().CosTheta();
        
            //////////////////////////////////////////////////
            //           Loop over points
            //////////////////////////////////////////////////
            for(int ipoint = 0; ipoint < stktrack->GetNPoints(); ipoint++){

                if(stktrack->getHitMeasX(ipoint) <= -9999. || stktrack->getHitMeasY(ipoint) <= -9999.)
                    nholes++;
    
                for(int ixy = 0; ixy < 2; ixy++){
                    
                    DmpStkSiCluster * stkcluster;
                    if(ixy == 0){
                        stkcluster = stktrack -> pClusterX(ipoint);
                    }
                    else{
                        stkcluster = stktrack -> pClusterY(ipoint);
                    }
                    if(!stkcluster) continue;
            
                    //float eta = CalcEta(stkcluster);
                    //int etaReg = GetEtaRegion(eta);
                    //int vaNumber = GetVANumber(stkcluster->getFirstStrip(), stkcluster->getLastStrip());
                    float energy = stkcluster->getEnergy()*cosTheta;

                    //if(eta == 0. || eta == 1.) continue;
                    //if(stkcluster->getNstrip() == 1) continue;
                        //hEnergy1stripCluster->Fill(energy);

                    //if(etaReg < 0) continue;
                    //if(vaNumber < 0) continue;

                    if(ixy == 0) {
                        nClusX++;
                        sMeanX += energy;
                    }
                    if(ixy == 1) {
                        nClusY++;
                        sMeanY += energy;
                    }
                } // end of loop over x and y clusters
            } // end of loop over points

            if(nClusX < MINTRACKHITS || nClusY < MINTRACKHITS) //continue;
            nTracks[1]++;

            ////// Tight criteria from track quality requirements /////
            if(chisqndof > 15.) continue;
            nTracks[2]++;
            
            if(!(stktrack->getImpactPointHasX()) || !(stktrack->getImpactPointHasY()) || nonoverlaps > 0) continue;
            nTracks[3]++;
            
            if(nholes > 0) continue;
            nTracks[4]++;

            /////////////////////

            if(stktrack->getNhitXY() < MINTRACKHITS) continue;
            nTracks[5]++;

            sMeanX = std::sqrt(sMeanX/nClusX/60.);
            sMeanY = std::sqrt(sMeanY/nClusY/60.);
            hSmeanX->Fill(sMeanX);
            hSmeanY->Fill(sMeanY);
            if(sMeanX < .9 || sMeanX > 1.1) continue;
            if(sMeanY < .9 || sMeanY > 1.1) continue;
            nTracks[6]++;
            
            vSelectedTracks.push_back(TRACK());
            vSelectedTracks.back().trackIndex = itrack;
            vSelectedTracks.back().nClustersX = nClusX;
            vSelectedTracks.back().nClustersY = nClusY;
            vSelectedTracks.back().nClustersTotal = nClusX + nClusY;

        } // end of loop over tracks

        if(nTracks[0]>0) 
            nEntriesAfterSelection[0]++;
        if(nTracks[0]==1)
            nEntriesAfterSelection1track[0]++;
        
        if(nTracks[1]>0) 
            nEntriesAfterSelection[1]++;
        if(nTracks[1]==1)
            nEntriesAfterSelection1track[1]++;
        
        if(nTracks[2]>0) 
            nEntriesAfterSelection[2]++;
        if(nTracks[2]==1)
            nEntriesAfterSelection1track[2]++;
        
        if(nTracks[3]>0) 
            nEntriesAfterSelection[3]++;
        if(nTracks[3]==1)
            nEntriesAfterSelection1track[3]++;
        
        if(nTracks[4]>0) 
            nEntriesAfterSelection[4]++;
        if(nTracks[4]==1)
            nEntriesAfterSelection1track[4]++;
        
        if(nTracks[5]>0) 
            nEntriesAfterSelection[5]++;
        if(nTracks[5]==1)
            nEntriesAfterSelection1track[5]++;
        
        if(nTracks[6]>0) 
            nEntriesAfterSelection[6]++;
        if(nTracks[6]==1)
            nEntriesAfterSelection1track[6]++;
       
        /////////////////////////////////////////////////////////////

        if(vSelectedTracks.size()==0) {
            //std::cout << "None of the tracks for event " << ientry << " passed the selection cuts!" << std::endl;
        }

        if(vSelectedTracks.size()>0) {
            //std::cout << "At least 1 track of event " << ientry << " passed all cluster requirements! :) " << std::endl;
        }
        
        if(vSelectedTracks.size()==1) {
            //std::cout << "Event " << ientry << " has exactly one track!!! " << std::endl;
        }
    /*
        // if we have a one-track event, we save it!
        nClustersX = vSelectedTracks.at(0).nClustersX;
        nClustersY = vSelectedTracks.at(0).nClustersY;
        nClustersTotal = vSelectedTracks.at(0).nClustersTotal;

        DmpStkTrack* selectedTrack = (DmpStkTrack*) stkhelper->GetTrack(vSelectedTracks.at(0).trackIndex);
        float cosTheta = selectedTrack->getDirection().CosTheta();

        for (int ipoint = 0; ipoint < selectedTrack->GetNPoints(); ipoint++) {

            for(int ixy = 0; ixy < 2; ixy++){

                DmpStkSiCluster* stkcluster;
                float energyAdc = 0.;

                if(ixy == 0){
                    stkcluster = selectedTrack -> GetClusterX(ipoint,stkclusters);
                }
                else{
                    stkcluster = selectedTrack -> GetClusterY(ipoint,stkclusters);
                }

               for(int istrip = 0; istrip < stkcluster->getNstrip(); istrip++){
                    energyAdc += stkcluster->GetAdcValue(istrip,stkladderadc);
                }
                
                float energy = stkcluster->getEnergy()*cosTheta;
                float eta = CalcEta(stkcluster);
                int etaReg = GetEtaRegion(eta);
                int ladNumber = stkcluster->getLadderHardware();
                int vaNumber = GetVANumber(stkcluster->getFirstStrip(),stkcluster->getLastStrip());

                int index = ipoint * 2 + ixy;
                clusterEnergy[index] = energy;
                clusterEnergyAdc[index] = energyAdc;
                clusterEta[index] = eta;
                clusterEtaRegion[index] = etaReg;
                clusterLadder[index] = ladNumber;
                clusterVA[index] = vaNumber;
                
            } // end of loop over x and y clusters
        } // end of loop over points
        */        
        outputTree->Fill();

    }// end of loop over entries

    std::cout << "nEntries: " << nEntries << std::endl;
    std::cout << "Selection 0: " << nEntriesAfterSelection[0] << "\t" << nEntriesAfterSelection1track[0] << std::endl;
    std::cout << "Selection 1: " << nEntriesAfterSelection[1] << "\t" << nEntriesAfterSelection1track[1] << std::endl;
    std::cout << "Selection 2: " << nEntriesAfterSelection[2] << "\t" << nEntriesAfterSelection1track[2] << std::endl;
    std::cout << "Selection 3: " << nEntriesAfterSelection[3] << "\t" << nEntriesAfterSelection1track[3] << std::endl;
    std::cout << "Selection 4: " << nEntriesAfterSelection[4] << "\t" << nEntriesAfterSelection1track[4] << std::endl;
    std::cout << "Selection 5: " << nEntriesAfterSelection[5] << "\t" << nEntriesAfterSelection1track[5] << std::endl;
    std::cout << "Selection 6: " << nEntriesAfterSelection[6] << "\t" << nEntriesAfterSelection1track[6] << std::endl;

    hSmeanX->Write();
    hSmeanY->Write();
    hEnergy1stripCluster->Write();
    //hNTracks->Write();
    //hNTracksNoCuts->Write();

    hNEntriesAfterSelection -> SetBinContent(1,nEntries);
    for(int i = 0; i < nSel; i++){
        hNEntriesAfterSelection -> SetBinContent(i+2,nEntriesAfterSelection[i]);
    }
    hNEntriesAfterSelection1track -> SetBinContent(1,nEntries);
    for(int i = 0; i < nSel; i++){
        hNEntriesAfterSelection1track -> SetBinContent(i+2,nEntriesAfterSelection1track[i]);
    }
    gStyle->SetOptStat(0);
    hNEntriesAfterSelection->Write();
    hNEntriesAfterSelection1track->Write();
    outputTree->Write();
    inputFile->Close();
    outputFile->Close();


    return 0;
}
