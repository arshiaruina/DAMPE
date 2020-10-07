//
// Created by Arshia Ruina on 10.01.20.
//

/* Code to compute eta values
 * Author: Arshia Ruina
 * To compile, g++ -Wall `root-lconfig --cflags --libs` -I$DMPSWSYS/include -L$DMPSWSYS/lib -lDmpEvent -lDmpService eta.cc -o ../submit/eta.exe
 */

// example file
// /beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T001308_20181019T002610_00000_jEwJChikPnVApmeJU6og/DAMPE_2A_OBS_20181019_20181019T001308_20181019T002610_00000_jEwJChikPnVApmeJU6og.root

//#include "../inc/va_equalisation.h"
#include "../inc/clusterEnergyGen.h"
//#include "../inc/track_selection.hpp"

using namespace std;

bool** clusterEnergyGen::read_bad_channels_file(const char* filename) {
    bool **badchannels = new bool *[N_LADDER];
    for (int i = 0; i < N_LADDER; i++) {
        badchannels[i] = new bool[N_CHANNEL];
        for (int j = 0; j < N_CHANNEL; j++) {
            badchannels[i][j] = false;
        }
    }
    ifstream infile(filename);
    if (!infile)
    {
        std::cout<<"Can't open Bad Channels file: "<<filename<<" ==> throwing exception!"<<std::endl;
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
    std::cout<<"Done reading bad channels"<<std::endl;
    return badchannels;
} // end read_bad_channels_file


bool clusterEnergyGen::is_cluster_bad_channel(DmpStkSiCluster * cluster, bool** badchannels) {
    int ladder = cluster->getLadderHardware(); //ladder num
    int minc = cluster->getIndex1();//address of the first strip in the cluster (0-363)
    int maxc = cluster->getIndex1() + cluster->getNstrip() -1;//last strip of the cluster

    for(int i = minc; i <= maxc; i++){
        if(badchannels[ladder][i]) return true;
    }
    return false;
} // end is_cluster_bad_channel

int clusterEnergyGen::GetEtaRegion(double eta){

    int etaReg = -99;
    if(eta < 0.2 || eta > 0.8)
        etaReg = 0;
    if(eta > 0.4 && eta < 0.6)
        etaReg = 1;
    return etaReg;
}

int clusterEnergyGen::GetVANumber(int firstStrip, int lastStrip){

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

/*bool IsClusterAtVAEdge(int firstStrip, int lastStrip, int vaNumber){

    if(firstStrip == vaNumber*64 || lastStrip == (vaNumber+1)*64-1)
        return true;
    else
        return false;
}*/

// iladder          0 - 47          48 - 95
// X ladder IDs     48 - 95         144 - 191
// Y ladder IDs     0 - 47          96 - 143

bool clusterEnergyGen::IsLadderX1(int ladder){
    if(ladder >= 48 && ladder < 96)
        return true;
    else
        return false;
}

bool clusterEnergyGen::IsLadderX2(int ladder){
    if(ladder >= 144 && ladder < 192)
        return true;
    else
        return false;
}

bool clusterEnergyGen::IsLadderY1(int ladder){
    if(ladder >= 0 && ladder < 48)
        return true;
    else
        return false;
}

bool clusterEnergyGen::IsLadderY2(int ladder){
    if(ladder >= 96 && ladder < 144)
        return true;
    else
        return false;
}

double clusterEnergyGen::CalcEta(DmpStkSiCluster* cluster){

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
