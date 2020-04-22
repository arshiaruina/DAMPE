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

//int main(int argc, char** argv){
int main() {

    std::cout << "hi" << std::endl;
    //std::string fileListName = "/beegfs/users/ruina/VAequalisation/fileList.txt";
    std::string line;
    std::ifstream fileList("/beegfs/users/ruina/VAequalisation/fileList.txt");

    //fileList.open(fileListName.c_str());
    //std::cout << fileListName << std::endl;

    //if(!fileList.good()){
    //    std::cout << "Problem with energy file!" << std::endl;
    //    return 1;
    //}

    while(getline(fileList,line)) {

        std::istringstream fileName(line);   
        //std::string inputFileName = argv[1];
        std::string inputFileName = fileName.str();
        std::cout << inputFileName << std::endl;
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
        std::cout << outputFileName << std::endl;
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
        int nEntriesAfterSelectionMulTrack[nSel];
        TH1I* hNEntriesAfterSelection = new TH1I("hNEntriesAfterSelection","Selections",nSel+1,0,nSel);
        TH1I* hNEntriesAfterSelection1track = new TH1I("hNEntriesAfterSelection1track","Selections (1 track events)",nSel+1,0,nSel);
        //TH1I* hNTracksNoCuts    = new TH1I("hNTracksNoCuts","No. of tracks with no selection cuts",10,0,10);    
        //TH1I* hNTracks          = new TH1I("hNTracks","No. of tracks that passed all cuts",10,0,10);    
        TH1I* hNClusX           = new TH1I("hNClusX","hNClusX",10,0.,10.);    
        TH1I* hNClusY           = new TH1I("hNClusY","hNClusY",10,0.,10.);    
        TH1D* hChiSqNDOF        = new TH1D("hChiSqNDOF","hChiSqNDOF",100,0.,100.);
        TH1D* hSmeanX           = new TH1D("hSmeanX","hSmeanX",50,0.,5.);    
        TH1D* hSmeanY           = new TH1D("hSmeanY","hSmeanY",50,0.,5.);    
        TH1D* hSmeanXAfterSel           = new TH1D("hSmeanXAfterSel","hSmeanXAfterSel",50,0.,5.);    
        TH1D* hSmeanYAfterSel           = new TH1D("hSmeanYAfterSel","hSmeanYAfterSel",50,0.,5.);    
        TH1F* hEnergy1stripCluster = new TH1F("hEnergy1stripCluster","Energy for 1-strip clusters",200,0.,200.);    

        for(int i = 0; i < nSel; i++){
            nEntriesAfterSelection[i] = 0;
            nEntriesAfterSelection1track[i] = 0;
            nEntriesAfterSelectionMulTrack[i] = 0;
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
                //std::cout << "entry: " << ientry << " track: " << itrack << " GetNPoints: " << stktrack->GetNPoints() << " getNhitXY: " << stktrack->getNhitXY() << std::endl;
            
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
               
                //std::cout << "debug 1" << std::endl; 
                sMeanX = std::sqrt(sMeanX/nClusX/60.);
                sMeanY = std::sqrt(sMeanY/nClusY/60.);
                //std::cout << "debug 2" << std::endl; 
                
                hNClusX->Fill(nClusX);
                hNClusY->Fill(nClusY);
                hChiSqNDOF->Fill(chisqndof);
                hSmeanX->Fill(sMeanX);
                hSmeanY->Fill(sMeanY);
                //std::cout << "debug 3" << std::endl; 
                
                ////// Selections //////

                //if(nClusX != MINTRACKHITS || nClusY != MINTRACKHITS) continue;
                //nTracks[1]++;

                if(chisqndof > 15.) continue;
                nTracks[2]++;
                
                //if(!(stktrack->getImpactPointHasX()) || !(stktrack->getImpactPointHasY()) || nonoverlaps > 0) continue;
                //nTracks[3]++;
                
                //if(nholes > 0) continue;
                //nTracks[4]++;

                //if(stktrack->GetNPoints() < MINTRACKHITS) continue;
                if(stktrack->getNhitXY() != MINTRACKHITS) continue;
                nTracks[5]++;

                hSmeanXAfterSel->Fill(sMeanX);
                hSmeanYAfterSel->Fill(sMeanY);
                if(sMeanX < .9 || sMeanX > 1.1) continue;
                if(sMeanY < .9 || sMeanY > 1.1) continue;
                nTracks[6]++;
                
                vSelectedTracks.push_back(TRACK());
                vSelectedTracks.back().trackIndex = itrack;
                vSelectedTracks.back().nClustersX = nClusX;
                vSelectedTracks.back().nClustersY = nClusY;
                vSelectedTracks.back().nClustersTotal = nClusX + nClusY;

            } // end of loop over tracks
          
            //if(nTracks[1]==1 && nTracks[2]==1) rach++;
            //if(nTracks[1]!=1 && nTracks[2]==1) rach1++;
            //if(nTracks[1]==1 && nTracks[2]!=1) rach2++;
 
            //if(nTracks[1]==1 && nTracks[2]==1)
            //    std::cout << "For event " << ientry << " nTracks[1] " << nTracks[1] << " nTracks[2] " << nTracks[2] << std::endl;
            for(int iSel = 0; iSel < nSel; iSel++){
                if(nTracks[iSel] > 0) {
                    nEntriesAfterSelection[iSel]++;
                }
                if(nTracks[iSel] == 1){
                    nEntriesAfterSelection1track[iSel]++;
                }
                if(nTracks[iSel] > 1){
                    nEntriesAfterSelectionMulTrack[iSel]++;
                }
            }

            //if(nTracks[0]>0) 
            //    nEntriesAfterSelection[0]++;
            //if(nTracks[0]==1)
            //    nEntriesAfterSelection1track[0]++;
            //
            //if(nTracks[1]>0) 
            //    nEntriesAfterSelection[1]++;
            //if(nTracks[1]==1)
            //    nEntriesAfterSelection1track[1]++;
            //
            //if(nTracks[2]>0) 
            //    nEntriesAfterSelection[2]++;
            //if(nTracks[2]==1)
            //    nEntriesAfterSelection1track[2]++;
            //
            //if(nTracks[3]>0) 
            //    nEntriesAfterSelection[3]++;
            //if(nTracks[3]==1)
            //    nEntriesAfterSelection1track[3]++;
            //
            //if(nTracks[4]>0) 
            //    nEntriesAfterSelection[4]++;
            //if(nTracks[4]==1)
            //    nEntriesAfterSelection1track[4]++;
            //
            //if(nTracks[5]>0) 
            //    nEntriesAfterSelection[5]++;
            //if(nTracks[5]==1)
            //    nEntriesAfterSelection1track[5]++;
            //
            //if(nTracks[6]>0) 
            //    nEntriesAfterSelection[6]++;
            //if(nTracks[6]==1)
            //    nEntriesAfterSelection1track[6]++;
           
            /////////////////////////////////////////////////////////////
/*
            if(vSelectedTracks.size()==0) {
                //std::cout << "None of the tracks for event " << ientry << " passed the selection cuts!" << std::endl;
            }

            if(vSelectedTracks.size()>0) {
                //std::cout << "At least 1 track of event " << ientry << " passed all cluster requirements! :) " << std::endl;
            }
            
            if(vSelectedTracks.size()==1) {
                //std::cout << "Event " << ientry << " has exactly one track!!! " << std::endl;
            }
        
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
        std::cout << "Selection 0: " << nEntriesAfterSelection[0] << "\t" << nEntriesAfterSelection1track[0] << "\t" << nEntriesAfterSelectionMulTrack[0] << std::endl;
        std::cout << "Selection 1: " << nEntriesAfterSelection[1] << "\t" << nEntriesAfterSelection1track[1] << "\t" << nEntriesAfterSelectionMulTrack[1] << std::endl;
        std::cout << "Selection 2: " << nEntriesAfterSelection[2] << "\t" << nEntriesAfterSelection1track[2] << "\t" << nEntriesAfterSelectionMulTrack[2] << std::endl;
        std::cout << "Selection 3: " << nEntriesAfterSelection[3] << "\t" << nEntriesAfterSelection1track[3] << "\t" << nEntriesAfterSelectionMulTrack[3] << std::endl;
        std::cout << "Selection 4: " << nEntriesAfterSelection[4] << "\t" << nEntriesAfterSelection1track[4] << "\t" << nEntriesAfterSelectionMulTrack[4] << std::endl;
        std::cout << "Selection 5: " << nEntriesAfterSelection[5] << "\t" << nEntriesAfterSelection1track[5] << "\t" << nEntriesAfterSelectionMulTrack[5] << std::endl;
        std::cout << "Selection 6: " << nEntriesAfterSelection[6] << "\t" << nEntriesAfterSelection1track[6] << "\t" << nEntriesAfterSelectionMulTrack[6] << std::endl;

        hNClusX->Write();
        hNClusY->Write();
        hChiSqNDOF->Write();
        hSmeanX->Write();
        hSmeanY->Write();
        hSmeanXAfterSel->Write();
        hSmeanYAfterSel->Write();
        hEnergy1stripCluster->Write();
        hNEntriesAfterSelection->Write();
        hNEntriesAfterSelection1track->Write();
        //hNTracks->Write();
        //hNTracksNoCuts->Write();

        TCanvas *c1 = new TCanvas("c1","c1",800,600);
        hNEntriesAfterSelection -> SetBinContent(1,nEntries);
        for(int i = 0; i < nSel; i++){
            hNEntriesAfterSelection -> SetBinContent(i+2,nEntriesAfterSelection[i]);
        }
        hNEntriesAfterSelection1track -> SetBinContent(1,nEntries);
        for(int i = 0; i < nSel; i++){
            hNEntriesAfterSelection1track -> SetBinContent(i+2,nEntriesAfterSelection1track[i]);
        }
        gStyle->SetOptStat(0);
        hNEntriesAfterSelection->SetLineColor(kBlack);
        hNEntriesAfterSelection->Draw("hist");
        hNEntriesAfterSelection1track->SetLineColor(kBlue);
        hNEntriesAfterSelection1track->Draw("hist same");
        c1->Write(); 
        outputTree->Write();
        inputFile->Close();
        outputFile->Close();

        std::cout << "end of this file, going to next ........" << std::endl;
    }
    //return 0;
}
