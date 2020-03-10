// Date: 05.03.2020
// Source: https://www.nevis.columbia.edu/~seligman/root-class/files/MakeNtuple.C

#include "../inc/makeNtuple.h"

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
}

bool is_cluster_bad_channel(DmpStkSiCluster * cluster, bool** badchannels) {
    int ladder = cluster->getLadderHardware(); //ladder num
    int minc = cluster->getIndex1();//address of the first strip in the cluster (0-363)
    int maxc = cluster->getIndex1() + cluster->getNstrip() -1;//last strip of the cluster

    for(int i = minc; i <= maxc; i++){
        if(badchannels[ladder][i]) return true;
    }
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
    TTree *inputTree = (TTree*)inputFile->Get("CollectionTree");
    if(inputTree == 0){
        std::cout << "[ERROR] CollectionTree not found!" << std::endl;
        return 0;
    }

    bool** badchannels =  read_bad_channels_file("/beegfs/users/ruina/VAequalisation/resources/bad_chan.txt");

    TClonesArray* stktracks = new TClonesArray("DmpStkTrack");
    inputTree->SetBranchAddress("StkKalmanTracks",&stktracks);
    DmpStkTrackHelper* stkhelper = new DmpStkTrackHelper(stktracks,false,0,0);

    TClonesArray* stkclusters  = new TClonesArray("DmpStkSiCluster");
    inputTree->SetBranchAddress("StkClusterCollection",&stkclusters);

    TClonesArray* stkladderadc = new TClonesArray("DmpStkLadderAdc"); 
    inputTree->SetBranchAddress("DmpStkLadderAdcCollection",&stkladderadc); 

    std::size_t found = inputFileName.find_last_of("/");
    //std::string outputFileName = "/beegfs/users/ruina/VAequalisation/testNtuples/" + inputFileName.substr(found+1);
    std::string outputFileName = "/beegfs/users/ruina/VAequalisation/ntuples_20181001_20181009/" + inputFileName.substr(found+1);
    TFile* outputFile = new TFile(outputFileName.c_str(),"recreate");
    TTree* outputTree = new TTree("MySelectionTree","Selections for VA calibration");

    int nClustersX;
    int nClustersY;
    int nClustersTotal;
    // size of the folllowing vectors will be the number of clusters in the event
    //std::vector<float> clusterEnergy;
    //std::vector<float> clusterEnergyAdc;
    //std::vector<float> clusterEta;
    //std::vector<int> clusterEtaRegion;
    //std::vector<int> clusterLadder;
    //std::vector<int> clusterVA;
    std::vector<float> clusterEnergy(MINTRACKHITS);
    std::vector<float> clusterEnergyAdc(MINTRACKHITS);
    std::vector<float> clusterEta(MINTRACKHITS);
    std::vector<int> clusterEtaRegion(MINTRACKHITS);
    std::vector<int> clusterLadder(MINTRACKHITS);
    std::vector<int> clusterVA(MINTRACKHITS);
    
    outputTree->Branch("nClustersX", &nClustersX);
    outputTree->Branch("nClustersY", &nClustersY);
    outputTree->Branch("nClustersTotal", &nClustersTotal);
    outputTree->Branch("clusterEnergy", &clusterEnergy);
    outputTree->Branch("clusterEnergyAdc", &clusterEnergyAdc);
    outputTree->Branch("clusterEta", &clusterEta);
    outputTree->Branch("clusterEtaRegion", &clusterEtaRegion);
    outputTree->Branch("clusterLadder", &clusterLadder);
    outputTree->Branch("clusterVA", &clusterVA);
   
    int nEntries = inputTree->GetEntries();
    for(int ientry = 0; ientry < nEntries; ++ientry)
    //for(int ientry = 0; ientry < 100; ++ientry)
    {
        float progress = 100.0 * ((float) ientry) / ((float) nEntries);
        std::cout << std::setprecision(3) << "[ " << progress << " % ] \r";

        inputTree->GetEntry(ientry);

        stkhelper->SortTracks(4,false);
        if(stkhelper->GetSize() == 0) continue;

        std::vector<TRACK> vSelectedTracks;

        for(int itrack = 0; itrack < stkhelper->GetSize(); itrack++){ 

            DmpStkTrack* stktrack = (DmpStkTrack*) stkhelper->GetTrack(itrack);
            float cosTheta = stktrack->getDirection().CosTheta();

            float sMeanX = 0.;
            float sMeanY = 0.;
            int nClusX = 0.;
            int nClusY = 0.;

            if(stktrack->getNhitXY() < MINTRACKHITS) continue;

            for (int ipoint = 0; ipoint < stktrack->GetNPoints(); ipoint++) {

                DmpStkSiCluster* stkcluster;
                float energy = 0.;
                float eta = 0.;
                int etaReg = 0;
                int vaNumber = 0;   

                for(int ixy = 0; ixy < 2; ixy++){

                    if(ixy == 0){
                        stkcluster = stktrack -> GetClusterX(ipoint,stkclusters);
                    }
                    else{
                        stkcluster = stktrack -> GetClusterY(ipoint,stkclusters);
                    }

                    if(!stkcluster) continue;
                    if(is_cluster_bad_channel(stkcluster, badchannels)) continue;
                    
                    eta = CalcEta(stkcluster);
                    etaReg = GetEtaRegion(eta);
                    vaNumber = GetVANumber(stkcluster->getFirstStrip(), stkcluster->getLastStrip());        
                    if(eta == 0. || eta == 1. || etaReg < 0 || vaNumber < 0) continue;

                    energy = stkcluster->getEnergy()*cosTheta;

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

            if(nClusX < MINTRACKHITS || nClusY < MINTRACKHITS) continue;

            sMeanX = std::sqrt(sMeanX/nClusX/60.);
            sMeanY = std::sqrt(sMeanY/nClusY/60.);

            if(sMeanX < .9 || sMeanX > 1.1) continue;
            if(sMeanY < .9 || sMeanY > 1.1) continue;

            vSelectedTracks.push_back(TRACK());
            vSelectedTracks.back().trackIndex = itrack;
            vSelectedTracks.back().nClustersX = nClusX;
            vSelectedTracks.back().nClustersY = nClusY;
            vSelectedTracks.back().nClustersTotal = nClusX + nClusY;

        } // end of loop over tracks

        if(vSelectedTracks.size()!=1) continue;

        // if we have a one-track event, we save it!
        nClustersX = vSelectedTracks.at(0).nClustersX;
        nClustersY = vSelectedTracks.at(0).nClustersY;
        nClustersTotal = vSelectedTracks.at(0).nClustersTotal;

        DmpStkTrack* selectedTrack = (DmpStkTrack*) stkhelper->GetTrack(vSelectedTracks.at(0).trackIndex);
        float cosTheta = selectedTrack->getDirection().CosTheta();

        for (int ipoint = 0; ipoint < selectedTrack->GetNPoints(); ipoint++) {

            for(int ixy = 0; ixy < 2; ixy++){

                DmpStkSiCluster* stkcluster;
                float energy = 0.;
                float energyAdc = 0.;
                float eta = 0.;
                int etaReg = 0; 
                int ladNumber = 0;
                int vaNumber = 0;

                if(ixy == 0){
                    stkcluster = selectedTrack -> GetClusterX(ipoint,stkclusters);
                }
                else{
                    stkcluster = selectedTrack -> GetClusterY(ipoint,stkclusters);
                }

               for(int istrip = 0; istrip < stkcluster->getNstrip(); istrip++){
                    energyAdc += stkcluster->GetAdcValue(istrip,stkladderadc);
                }
                
                energy = stkcluster->getEnergy()*cosTheta;
                eta = CalcEta(stkcluster);
                etaReg = GetEtaRegion(eta);
                ladNumber = stkcluster->getLadderHardware();
                vaNumber = GetVANumber(stkcluster->getFirstStrip(),stkcluster->getLastStrip());

                clusterEnergy.push_back(energy);
                clusterEnergyAdc.push_back(energyAdc);
                clusterEta.push_back(eta);
                clusterEtaRegion.push_back(etaReg);
                clusterLadder.push_back(ladNumber);
                clusterVA.push_back(vaNumber);                
            } // end of loop over x and y clusters
        } // end of loop over points
             
        outputTree->Fill();

    }// end of loop over entries

    outputTree->Write();
    inputFile->Close();
    outputFile->Close();

    return 0;
}
