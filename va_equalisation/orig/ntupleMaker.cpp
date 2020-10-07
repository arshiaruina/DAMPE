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

//struct CLUSTER
//{
//    float energy;
//    float energyAdc;
//    float eta;
//    //float etaReg;
//    int va;
//    int ladder;
//};

struct TRACK
{
    int trackIndex;
    int nClustersX;
    int nClustersY;
    //int nClustersTotal;
    //float sMeanX;
    //float sMeanY;
    //float sMean;
    //float direction;
    //std::vector<CLUSTER> cluster;
};

int main(int argc, char** argv){
//int main() {

    std::string line;
    std::string fileListName = argv[1];
    std::ifstream fileList(fileListName.c_str());
    //std::ifstream fileList("/beegfs/users/ruina/VAequalisation/fileList.txt");

    std::cout << "INFO: Reading file list..." << std::endl;

    while(getline(fileList,line)) {

        std::istringstream fileName(line);   
        //std::string inputFileName = argv[1];
        std::string inputFileName = fileName.str();
        //std::string inputFileName = "/dpnc" + fileName.str(); // for running on baobab
        std::cout << inputFileName << std::endl;
        TFile *inputFile = new TFile(inputFileName.c_str());
        if(inputFile == 0){
            std::cout << "[ERROR] File " << inputFileName << " not found!" << std::endl;
            return 0;
        }

        std::size_t found = inputFileName.find_last_of("/");
        //std::string outputFileName = "/beegfs/users/ruina/VAequalisation/testNtuples/" + inputFileName.substr(found+1);
        std::string outputFileName = "/beegfs/users/ruina/VAequalisation/ntuples_201801_201802/" + inputFileName.substr(found+1);
        TFile* outputFile = new TFile(outputFileName.c_str(),"recreate");

        //-----------------------------------------------------//

        std::cout << "DEBUG: Initialising meta data tree..." << std::endl;
        //int ixk = 0;
        int const nSTKchannels = 384;
        int const nSTKladders  = 192;
        int nLadders = 192, nChannels = 384;
        float noiseInfo[nSTKladders][nSTKchannels];
        TH2I* hNoiseInfo = new TH2I("hNoiseInfo","Noise info",384,0,383,192,0,191);

        TTree *calib_tree = (TTree*) inputFile->Get("RunMetadataTree");
        TClonesArray* stk_ladder_sigmas      = new TClonesArray("DmpStkLadderCalibration");
        TClonesArray* stk_ladder_pedestals   = new TClonesArray("DmpStkLadderCalibration");
        TClonesArray* stk_ladder_badchannels = new TClonesArray("DmpStkLadderCalibration");
        TBranch *b_StkLadderSigmas;   //!
        TBranch *b_StkLadderPedestals;   //!
        TBranch *b_StkLadderBadchannels;   //!
        calib_tree->SetBranchAddress("StkLadderSigmas",     &stk_ladder_sigmas,       &b_StkLadderSigmas);
        calib_tree->SetBranchAddress("StkLadderPedestals",  &stk_ladder_pedestals ,   &b_StkLadderPedestals);
        calib_tree->SetBranchAddress("StkLadderBadchannels",&stk_ladder_badchannels , &b_StkLadderBadchannels);
        cout << "[INFO] Entries in RunMetadataTree " << calib_tree->GetEntries() << endl;

        TTree *outputTreeMeta = new TTree("MyMetaTree","Info on noisy channels");
        outputTreeMeta->Branch("nLadders", &nLadders, "nLadders/I");
        outputTreeMeta->Branch("nChannels", &nChannels, "nChannels/I");
        outputTreeMeta->Branch("noiseInfo", noiseInfo, "noiseInfo[nLadders][nChannels]/F");

        if(calib_tree->GetEntries() != 1) {
            cout << " calib_tree->GetEntries() = " << calib_tree->GetEntries() << ", not 1, stop!" << endl;
            exit (EXIT_FAILURE);
        }

        //-----------------------------------------------------//

        //std::stringstream ss;
        //std::ofstream noiseFile;
        //noiseFile.open("/beegfs/users/ruina/VAequalisation/out/noisyChannels.txt"); 
        //ss << "Ladder \t Channel \t Noise" << std::endl; 

        for (int ientry = 0; ientry < calib_tree->GetEntries(); ientry++) {

            b_StkLadderSigmas -> GetEntry(ientry);
            //b_StkLadderPedestals   ->GetEntry(ientry);    
            b_StkLadderBadchannels ->GetEntry(ientry);

            for (int iladsig = 0; iladsig < stk_ladder_sigmas->GetLast()+1; iladsig++){

                DmpStkLadderCalibration* sigma = (DmpStkLadderCalibration*) (stk_ladder_sigmas->ConstructedAt(iladsig));
                //int trbID    = sigma-> GetTRBNumber();  
                int ladderID = sigma->GetLadderID();
                //cout << " sigma: ladder/ladderID = " << iladsig << "/" << ladderID << endl;
                //std::cout << "iladsig: " << iladsig << " ladder id " << ladderID << std::endl;

                for(int ich = 0; ich < nSTKchannels; ich++) {

                    //int vaNumber = ich/64;

                    double noise = sigma->GetValue(ich);
                    //if(noise > 5.0) {
                    //    std::cout << iladsig << "\t" << vaNumber << "\t" << ich << "\t" << noise << std::endl;
                    //    //noisyChannels[iladsig][ich] = noise;
                    //}
                    //noiseInfo[iladsig][ich] = noise;
                    //if(noise > 5.0) 
                    //    std::cout << iladsig << "\t" << ich << "\t" << noiseInfo[iladsig][ich] << std::endl;
                    //std::cout << iladsig << "\t" << vaNumber << "\t" << ich << "\t" << noiseInfo[iladsig][ich] << std::endl;
                    //int bin= iladsig*384+ich+1;
                    //std::cout << iladsig << " " << ich << " " << bin << std::endl;
                    if(noise > 5.0){ 
                        hNoiseInfo->Fill(ich,iladsig);
                        //ss << iladsig << "\t" << ich << "\t" << noise << std::endl;
                        //hNoiseInfo->SetBinContent(bin,100);
                    }
                    if(noise <= 0) {
                        cout << " wrong STK noise for ladder " << ladderID << ", channel " << ich << ", noise " << noise << endl;
                        exit (EXIT_FAILURE);
                    }
                }
            }
        }
        hNoiseInfo->SetOption("COLZ");
        hNoiseInfo->SetContour(20);
        hNoiseInfo->SetStats(0);
        //hNoiseInfo->SetMinimum(0.9);
        //hNoiseInfo->SetMaximum(1.1);
        hNoiseInfo->Write();
        outputTreeMeta->Fill();
        outputTreeMeta->Write();

        //-----------------------------------------------------//

        DmpChain *dmpch = new DmpChain("CollectionTree");
        dmpch->Add(inputFileName.c_str());
        dmpch->GetListOfFiles()->Print();
        int nEntries = dmpch->GetEntries();
        //TClonesArray* stkladderadc = Get("DmpStkLadderAdc"); 
        //iinputTree->SetBranchAddress("DmpStkLadderAdcCollection",&stkladderadc);
        //std::size_t found = inputFileName.find_last_of("/");
        //std::string outputFileName = "/beegfs/users/ruina/VAequalisation/testNtuples/" + inputFileName.substr(found+1);
        //std::cout << outputFileName << std::endl;
        //std::string outputFileName = "/beegfs/users/ruina/VAequalisation/ntuples_20181001_20181009/" + inputFileName.substr(found+1);
        //std::string outputFileName = "/beegfs/users/ruina/VAequalisation/ntuples_201801_201802/" + inputFileName.substr(found+1);
        //std::string outputFileName = "/beegfs/users/ruina/VAequalisation/testNtuples/" + inputFileName.substr(found+1);
        //std::string outputFileName = "/dpnc/beegfs/users/ruina/VAequalisation/ntuples_201801_201802/" + inputFileName.substr(found+1); // for baobab
        //TFile* outputFile = new TFile(outputFileName.c_str(),"recreate");
        TTree* outputTree = new TTree("MySelectionTree","Selections for VA calibration");

        int nClustersX = 0;
        int nClustersY = 0;
        //int nClustersTotal = 0;
        int nSelTracks = 0;
        //std::vector <float> trackDir;
        //std::vector <float> trackSmeanX;
        //std::vector <float> trackSmeanY;
        //float clusterEnergy [MINTRACKCLUS] = {0.};
        //float clusterEnergyAdc [MINTRACKCLUS] = {0.};
        //float clusterEta [MINTRACKCLUS] = {0.};
        //int clusterEtaRegion [MINTRACKCLUS] = {0};
        //int clusterLadder [MINTRACKCLUS] = {0};
        //int clusterVA [MINTRACKCLUS] = {0}; 
        std::vector <float> clusterEnergy;
        std::vector <float> clusterEnergyAdc;
        std::vector <float> clusterEta;
        std::vector <int> clusterEtaRegion;
        std::vector <int> clusterLadder;
        std::vector <int> clusterVA; 
        std::vector <int> clusterFirstStrip; 
        std::vector <int> clusterLastStrip; 

        outputTree->Branch("nClustersX", &nClustersX, "nClustersX/I");
        outputTree->Branch("nClustersY", &nClustersY, "nClustersY/I");
        //outputTree->Branch("nClustersTotal", &nClustersTotal, "nClustersTotal/I");
        outputTree->Branch("nSelTracks", &nSelTracks, "nSelTracks/I");
        //outputTree->Branch("trackDir", &trackDir);
        //outputTree->Branch("trackSmeanX", &trackSmeanX);
        //outputTree->Branch("trackSmeanY", &trackSmeanY);
        //outputTree->Branch("clusterEnergy", clusterEnergy, "clusterEnergy[nClustersTotal]/F");
        //outputTree->Branch("clusterEnergyAdc", clusterEnergyAdc, "clusterEnergyAdc[nClustersTotal]/F");
        //outputTree->Branch("clusterEta", clusterEta, "clusterEta[nClustersTotal]/F");
        //outputTree->Branch("clusterEtaRegion", clusterEtaRegion, "clusterEtaRegion[nClustersTotal]/I");
        //outputTree->Branch("clusterLadder", clusterLadder, "clusterLadder[nClustersTotal]/I");
        //outputTree->Branch("clusterVA", clusterVA, "clusterVA[nClustersTotal]/I");
        outputTree->Branch("clusterEnergy", &clusterEnergy);
        outputTree->Branch("clusterEnergyAdc", &clusterEnergyAdc);
        outputTree->Branch("clusterEta", &clusterEta);
        outputTree->Branch("clusterEtaRegion", &clusterEtaRegion);
        outputTree->Branch("clusterLadder", &clusterLadder);
        outputTree->Branch("clusterVA", &clusterVA);
        outputTree->Branch("clusterFirstStrip", &clusterFirstStrip);
        outputTree->Branch("clusterLastStrip", &clusterLastStrip);

        int nSel = 7;
        int nEntriesAfterSelection[nSel];
        int nEntriesAfterSelection1track[nSel];
        int nEntriesAfterSelectionMulTrack[nSel];
        TH1I* hNEntriesAfterSelection = new TH1I("hNEntriesAfterSelection","Selections",nSel+1,0,nSel);
        TH1I* hNEntriesAfterSelection1track = new TH1I("hNEntriesAfterSelection1track","Selections (1 track events)",nSel+1,0,nSel);
        //TH1I* hNTracksNoCuts    = new TH1I("hNTracksNoCuts","No. of tracks with no selection cuts",10,0,10);    
        //TH1I* hNTracks          = new TH1I("hNTracks","No. of tracks that passed all cuts",10,0,10);    
        TH1I* hNClusX       = new TH1I("hNClusX","hNClusX",10,0.,10.);    
        TH1I* hNClusY       = new TH1I("hNClusY","hNClusY",10,0.,10.);    
        TH1I* hNClusXMulStrip       = new TH1I("hNClusXMulStrip","hNClusXMulStrip",10,0.,10.);    
        TH1I* hNClusYMulStrip       = new TH1I("hNClusYMulStrip","hNClusYMulStrip",10,0.,10.);    
        TH1D* hChiSqNDOF    = new TH1D("hChiSqNDOF","hChiSqNDOF",100,0.,100.);
        TH1D* hSmean        = new TH1D("hSmean","hSmean",200,0.,5.);    
        TH1D* hSmeanX       = new TH1D("hSmeanX","hSmeanX",200,0.,5.);    
        TH1D* hSmeanY       = new TH1D("hSmeanY","hSmeanY",200,0.,5.);    
        //TH1D* hSmeanXAfterSel           = new TH1D("hSmeanXAfterSel","hSmeanXAfterSel",200,0.,5.);    
        //TH1D* hSmeanYAfterSel           = new TH1D("hSmeanYAfterSel","hSmeanYAfterSel",200,0.,5.);    
        //TH1F* hEnergy1stripCluster = new TH1F("hEnergy1stripCluster","Energy for 1-strip clusters",200,0.,200.);    

        for(int i = 0; i < nSel; i++){
            nEntriesAfterSelection[i] = 0;
            nEntriesAfterSelection1track[i] = 0;
            nEntriesAfterSelectionMulTrack[i] = 0;
        }
        
        std::cout << "INFO: Start loop over entries" << std::endl;
        //////////////////////////////////////////////////
        //             Loop over events
        //////////////////////////////////////////////////
        for(int ientry = 0; ientry < nEntries; ++ientry)
        //for(int ientry = 0; ientry < 50; ++ientry)
        {   

            float progress = 100.0 * ((float) ientry) / ((float) nEntries);
            std::cout << std::setprecision(3) << "[ " << progress << " % ] \r";

            //std::cout << "Entry: " << ientry << std::endl;
 
            std::vector<TRACK> vSelectedTracks;

            DmpEvent *pev = dmpch->GetDmpEvent();    
          
            TClonesArray* stkladderadc = pev->GetStkLadderAdcCollection(); 
  
            //trackDir.clear();
            //trackSmeanX.clear();
            //trackSmeanY.clear();

            int nTracks[nSel];
            for(int i = 0; i < nSel; i++){
                nTracks[i] = 0;
            }
            //nTracks[0] = pev->NStkKalmanTrack();
 
            //Excluding events with more than 10 tracks
            if(pev->NStkKalmanTrack()>10) continue;

            //////////////////////////////////////////////////
            //              Loop over tracks
            //////////////////////////////////////////////////
            for(int itrack = 0; itrack < pev->NStkKalmanTrack(); itrack++){
            
                DmpStkTrack* stktrack = pev->pStkKalmanTrack(itrack);
                if(pev->pStkKalmanTrack(itrack))
                    nTracks[0]++;
            
                float sMean = 0.;
                float sMeanX = 0.;
                float sMeanY = 0.;
                int nClusXMulStrip = 0.;
                int nClusYMulStrip = 0.;
                int nClusX = 0.;
                int nClusY = 0.;
                //int nholes = 0;
                //int nonoverlaps = stktrack->getNhitX() + stktrack->getNhitY() - 2 * stktrack->getNhitXY();
                double chisqndof = stktrack->getChi2NDOF();
                float cosTheta = stktrack->getDirection().CosTheta();
            
                //////////////////////////////////////////////////
                //           Loop over points
                //////////////////////////////////////////////////
                //std::cout << "entry: " << ientry << " track: " << itrack << " GetNPoints: " << stktrack->GetNPoints() << " getNhitXY: " << stktrack->getNhitXY() << std::endl;
            
                for(int ipoint = 0; ipoint < stktrack->GetNPoints(); ipoint++){

                    //if(stktrack->getHitMeasX(ipoint) <= -9999. || stktrack->getHitMeasY(ipoint) <= -9999.)
                    //    nholes++;
        
                    for(int ixy = 0; ixy < 2; ixy++){
                        
                        DmpStkSiCluster* stkcluster;
                        if(ixy == 0){
                            stkcluster = stktrack -> pClusterX(ipoint);
                        }
                        else{
                            stkcluster = stktrack -> pClusterY(ipoint);
                        }
                        if(!stkcluster) continue;
                
                        if(ixy == 0) {
                            nClusX++;
                        }
                        if(ixy == 1) {
                            nClusY++;
                        }
                        
                        //float eta = CalcEta(stkcluster);
                        //int etaReg = GetEtaRegion(eta);
                        //int vaNumber = GetVANumber(stkcluster->getFirstStrip(), stkcluster->getLastStrip());
                        float energy = stkcluster->getEnergy()*cosTheta;


                        //if(eta == 0. || eta == 1.) continue;
                        if(stkcluster->getNstrip() == 1) continue;//not taking 1 strip clusters in the computation of Smean
                            //hEnergy1stripCluster->Fill(energy);

                        //if(etaReg < 0) continue;
                        //if(vaNumber < 0) continue;
                    
                        //NOTE: Here, nClusX and nClusY now become the number of multi strip clusters, not total clusters

                        if(ixy == 0) {
                            nClusXMulStrip++;
                            sMeanX += energy;
                        }
                        if(ixy == 1) {
                            nClusYMulStrip++;
                            sMeanY += energy;
                        }
                        sMean += energy;
                    } // end of loop over x and y clusters
                } // end of loop over points
               
                //std::cout << "debug 1" << std::endl; 
                sMean = std::sqrt(sMean/(nClusXMulStrip+nClusYMulStrip)/60.);
                sMeanX = std::sqrt(sMeanX/nClusXMulStrip/60.);
                sMeanY = std::sqrt(sMeanY/nClusYMulStrip/60.);
                //std::cout << "debug 2" << std::endl; 
                
                hNClusX->Fill(nClusX);
                hNClusY->Fill(nClusY);
                hNClusXMulStrip->Fill(nClusXMulStrip);
                hNClusYMulStrip->Fill(nClusYMulStrip);
                hChiSqNDOF->Fill(chisqndof);
                hSmean->Fill(sMean);
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

                //hSmeanXAfterSel->Fill(sMeanX);
                //hSmeanYAfterSel->Fill(sMeanY);
                //Stefania's cuts
                if(sMeanX < .72 || sMeanX > 1.16) continue;
                if(sMeanY < .72 || sMeanY > 1.16) continue;
                //if(sMeanX < .9 || sMeanX > 1.1) continue;
                //if(sMeanY < .9 || sMeanY > 1.1) continue;
                nTracks[6]++;
                
                vSelectedTracks.push_back(TRACK());
                vSelectedTracks.back().trackIndex = itrack;
                vSelectedTracks.back().nClustersX = nClusX;
                vSelectedTracks.back().nClustersY = nClusY;
                //vSelectedTracks.back().nClustersTotal = nClusX + nClusY;
                //vSelectedTracks.back().sMean = sMean;
                //vSelectedTracks.back().sMeanX = sMeanX;
                //vSelectedTracks.back().sMeanY = sMeanY;
                //vSelectedTracks.back().direction = cosTheta;

            } // end of loop over tracks
          
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

            //std::cout << "Filling ntuple..." << std::endl;

            ///////////////////////////////
            //        Filling ntuples
            ///////////////////////////////

            nSelTracks = vSelectedTracks.size();
            if(!nSelTracks) continue;

            //std::cout << "[DEBUG] Entry: " << ientry << " nSelTracks: " << nSelTracks << std::endl;

            int nClusXTot = 0;
            int nClusYTot = 0;
            //trackSmean.clear();
            //trackSmeanX.clear();
            //trackSmeanY.clear();
            clusterEnergy.clear();
            clusterEnergyAdc.clear();
            clusterEta.clear();
            clusterEtaRegion.clear();
            clusterLadder.clear();
            clusterVA.clear();
            clusterFirstStrip.clear();
            clusterLastStrip.clear();


            // Loop over selected tracks
            for(unsigned int itrack = 0; itrack < vSelectedTracks.size(); itrack++) {
    
                //std::cout << ientry << " " << itrack << " " << vSelectedTracks.at(itrack).nClustersX << std::endl;
                nClusXTot += vSelectedTracks.at(itrack).nClustersX;
                //std::cout << nClustersX << std::endl;
                nClusYTot += vSelectedTracks.at(itrack).nClustersY;
                //nClustersTotal += vSelectedTracks.at(itrack).nClustersTotal;

                //DmpStkTrack* selectedTrack = pev->pStkKalmanTrack(itrack);
                //std::cout << "[DEBUG] selected track index: " << vSelectedTracks.at(itrack).trackIndex << std::endl;
                DmpStkTrack* selectedTrack = pev->pStkKalmanTrack(vSelectedTracks.at(itrack).trackIndex);
                
                float cosTheta = selectedTrack->getDirection().CosTheta();
                //trackDir.push_back(vSelectedTracks.at(itrack).direction);
                //trackSmean.push_back(vSelectedTracks.at(itrack).sMean);
                //trackSmeanX.push_back(vSelectedTracks.at(itrack).sMeanX);
                //trackSmeanY.push_back(vSelectedTracks.at(itrack).sMeanY);

                for (int ipoint = 0; ipoint < selectedTrack->GetNPoints(); ipoint++) {

                    for(int ixy = 0; ixy < 2; ixy++){

                        DmpStkSiCluster* stkcluster;
                        float energyAdc = 0.;

                        if(ixy == 0){
                            stkcluster = selectedTrack -> pClusterX(ipoint);
                        }
                        else{
                            stkcluster = selectedTrack -> pClusterY(ipoint);
                        }

                       for(int istrip = 0; istrip < stkcluster->getNstrip(); istrip++){
                            energyAdc += stkcluster->GetAdcValue(istrip,stkladderadc);
                        }
                        
                        float energy = stkcluster->getEnergy()*cosTheta;
                        float eta = CalcEta(stkcluster);
                        int etaReg = GetEtaRegion(eta);
                        int ladNumber = stkcluster->getLadderHardware();
                        int firstStrip = stkcluster->getFirstStrip();
                        int lastStrip = stkcluster->getLastStrip();
                        int vaNumber = GetVANumber(firstStrip,lastStrip);

                        //int index = ipoint * 2 + ixy;
                        //clusterEnergy[index] = energy;
                        //clusterEnergyAdc[index] = energyAdc;
                        //clusterEta[index] = eta;
                        //clusterEtaRegion[index] = etaReg;
                        //clusterLadder[index] = ladNumber;
                        //clusterVA[index] = vaNumber;

                        clusterEnergy.push_back(energy);
                        clusterEnergyAdc.push_back(energyAdc);
                        clusterEta.push_back(eta);
                        clusterEtaRegion.push_back(etaReg);
                        clusterLadder.push_back(ladNumber);
                        clusterVA.push_back(vaNumber);   
                        clusterFirstStrip.push_back(firstStrip);   
                        clusterLastStrip.push_back(lastStrip);   
                                
                    } // end of loop over x and y clusters
                } // end of loop over points
            } // end of loop over selected tracks           
            
            nClustersX = nClusXTot;
            nClustersY = nClusYTot;

            //std::cout << "Entry: " << ientry << std::endl;
            //std::cout << "nSelTracks: " << nSelTracks << std::endl;
            //std::cout << "nClustersX: " << nClustersX << std::endl;
            //std::cout << "nClustersY: " << nClustersY << std::endl;
            //std::cout << "Size of cluster Energy: " << clusterEnergy.size() << std::endl;
            //std::cout << "Size of clusterEnergyAdc: " << clusterEnergyAdc.size() << std::endl;
            //std::cout << "Size of clusterEta: " << clusterEta.size() << std::endl;
            //std::cout << "Size of clusterEtaRegion: " << clusterEtaRegion.size() << std::endl;
            //std::cout << "Size of clusterLadder: " << clusterLadder.size() << std::endl;
            //std::cout << "Size of clusterVA: " << clusterVA.size() << std::endl;        

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
        hNClusXMulStrip->Write();
        hNClusYMulStrip->Write();
        hChiSqNDOF->Write();
        hSmean->Write();
        hSmeanX->Write();
        hSmeanY->Write();
        //hSmeanXAfterSel->Write();
        //hSmeanYAfterSel->Write();
        //hEnergy1stripCluster->Write();
        hNEntriesAfterSelection->Write();
        hNEntriesAfterSelection1track->Write();
        //hNTracks->Write();
        //hNTracksNoCuts->Write();

        //TCanvas *c1 = new TCanvas("c1","c1",800,600);
        //hNEntriesAfterSelection -> SetBinContent(1,nEntries);
        //for(int i = 0; i < nSel; i++){
        //    hNEntriesAfterSelection -> SetBinContent(i+2,nEntriesAfterSelection[i]);
        //}
        //hNEntriesAfterSelection1track -> SetBinContent(1,nEntries);
        //for(int i = 0; i < nSel; i++){
        //    hNEntriesAfterSelection1track -> SetBinContent(i+2,nEntriesAfterSelection1track[i]);
        //}
        //gStyle->SetOptStat(0);
        //hNEntriesAfterSelection->SetLineColor(kBlack);
        //hNEntriesAfterSelection->Draw("hist");
        //hNEntriesAfterSelection1track->SetLineColor(kBlue);
        //hNEntriesAfterSelection1track->Draw("hist same");
        //c1->Write(); 
 
        outputTree->Write();
        inputFile->Close();
        outputFile->Close();
        //c1->delete();
        //hNEntriesAfterSelection->delete();
        //hNEntriesAfterSelection1track->delete();
        //hNTracksNoCuts->delete(); 
        //hNTracks->delete();       
        //hNClusX->delete();
        //hNClusY->delete();
        //hNClusXMulStrip->delete();
        //hNClusYMulStrip->delete();
        //hChiSqNDOF->delete();
        //hSmean->delete();
        //hSmeanX->delete();
        //hSmeanY->delete();
        //hSmeanXAfterSel->delete();
        //hSmeanYAfterSel->delete();
        //hEnergy1stripCluster->delete();
        std::cout << "End of file ..." << std::endl;
    }
    //return 0;
}
