//
// Created by Arshia Ruina on 13.01.20.
//

#ifndef VA_EQUALISATION_CLUSTERENERGYGEN_H
#define VA_EQUALISATION_CLUSTERENERGYGEN_H

#include "va_equalisation.h"

#define N_CHANNEL 384

//TODO: required include files
//TODO: required defines

class clusterEnergyGen
{

public:

    clusterEnergyGen();
    ~clusterEnergyGen();

    bool** read_bad_channels_file(const char*);
    bool is_cluster_bad_channel(DmpStkSiCluster *, bool**);
    int GetEtaRegion(double);
    int GetVANumber(int,int);
    bool IsLadderX1(int);
    bool IsLadderX2(int);
    bool IsLadderY1(int);
    bool IsLadderY2(int);
    double CalcEta(DmpStkSiCluster*);

//protected:

    TH1D* hVAEnergyX[N_LADDER/2][N_VA][N_ETAREG];
    TH1D* hVAEnergyY[N_LADDER/2][N_VA][N_ETAREG];
    TH1D* hEtaX = new TH1D("hEtaX","Eta distribution for X planes; #eta; No. of events",50,0.,1.);
    TH1D* hEtaY = new TH1D("hEtaY","Eta distribution for Y planes; #eta; No. of events",50,0.,1.);
    //TH1F* hMeasCovXX      = new TH1F("hMeasCovXX","MeasCovXX",100,-0.1,0.1);
    //TH1F* hMeasCovYY      = new TH1F("hMeasCovYY","MeasCovYY",100,-0.1,0.1);
    //TH1F* hMeasCovXXSqrt  = new TH1F("hMeasCovXXSqrt","Sqrt of MeasCovXX",100,-0.1,0.1);
    //TH1F* hMeasCovYYSqrt  = new TH1F("hMeasCovYYSqrt","Sqrt of MeasCovYY",100,-0.1,0.1);
    //TH1F* hDistanceX      = new TH1F("hDistanceX","|MeasHitX - FiltHitX|",100,0.,0.5);
    //TH1F* hDistanceY      = new TH1F("hDistanceY","|MeasHitY - FiltHitY|",100,0.,0.5);
    //TH1I* hNoisyClusters  = new TH1I("hNoisyClusters","Clusters with a noisy channel",1,0,1);
    TH1I* hNTracksNoCuts    = new TH1I("hNTracksNoCuts","No. of tracks with no selection cuts",10,0,10);
    TH1I* hNTracks          = new TH1I("hNTracks","No. of tracks that passed all cuts",10,0,10);
    TH1D* hSmeanX           = new TH1D("hSmeanX","hSmeanX",50,0.,5.);
    TH1D* hSmeanY           = new TH1D("hSmeanY","hSmeanY",50,0.,5.);
    //int checkCounter[N_LADDER][N_VA] = {0};
    //int total[N_LADDER][N_VA] = {0};
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
    
};


#endif //VA_EQUALISATION_CLUSTERENERGYGEN_H
