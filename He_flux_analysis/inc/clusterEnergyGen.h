//
// Created by Arshia Ruina on 13.01.20.
//

#ifndef VA_EQUALISATION_CLUSTERENERGYGEN_H
#define VA_EQUALISATION_CLUSTERENERGYGEN_H

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
    bool isLadderX1(int);
    bool isLadderX2(int);
    bool isLadderY1(int);
    bool isLadderY2(int);
    double CalcEta(DmpStkSiCluster*);
};


#endif //VA_EQUALISATION_CLUSTERENERGYGEN_H
