/* Code to check which VAs have a considerable low enregy peak
 * Author: Arshia Ruina
 * To compile, g++ -Wall `root-config --cflags --libs` -I$DMPSWSYS/include -L$DMPSWSYS/lib -lDmpEvent -lDmpService eta.cc -o ../submit/eta.exe
 * */

// example file
// /beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T001308_20181019T002610_00000_jEwJChikPnVApmeJU6og/DAMPE_2A_OBS_20181019_20181019T001308_20181019T002610_00000_jEwJChikPnVApmeJU6og.root

#include "../inc/va_equalisation.h"

//using namespace std;

int main(int argc, char** argv) {

	//gSystem->Load("./libDmpEvent.so");

	std::string inFileName = argv[1];
	TFile *f = new TFile(inFileName.c_str());
    
    std::ofstream checkLowEnergyPeakFile;
    checkLowEnergyPeakFile.open("checkLowEnergyPeak.txt");
    std::stringstream ss, ss00, ss01, ss10, ss11, ss20, ss21, ss30, ss31;
    ss00 << "VAs for which no. of events with E<20GeV is 0-2%" << std::endl;
    ss10 << "VAs for which no. of events with E<20GeV is 2-5%" << std::endl;
    ss20 << "VAs for which no. of events with E<20GeV is 5-10%" << std::endl;
    ss30 << "VAs for which no. of events with E<20GeV is >10%" << std::endl;
    ss << std::setw(15) << std::left << "Hist. no.";
    ss << std::setw(20) << std::left << "Hist. name";
    ss << std::setw(15) << std::left << "Total entries";
    ss << std::setw(15) << std::left << "Entries w/ E<20GeV";
    ss << std::setw(15) << std::left << "Fraction";
    ss << std::endl; 
    std::ifstream histoListFile;
    histoListFile.open("histolist.txt");   
    TH1D* myhist;
    std::string histName;
    
    int ihist = 0;
        int count[4][2] = {0};

    while(std::getline(histoListFile, histName)){
        myhist = (TH1D*)f->Get(histName.c_str());
        int nTotalEntries = myhist->GetEntries();
        int nMyEntries  = 0;
        float fraction = 0.;    
        for(int ibin = 0; ibin < 20; ibin++){
               nMyEntries += myhist->GetBinContent(ibin);  
        }
        if(nTotalEntries!=0 && nMyEntries!=0){
            fraction = (float)nMyEntries / (float)nTotalEntries;
            fraction *= 100.;    
        } 
        ss << std::setw(15) << std::left << ihist;
        ss << std::setw(20) << std::left << histName;
        ss << std::setw(15) << std::left << nTotalEntries;
        ss << std::setw(15) << std::left << nMyEntries;
        ss << setprecision(2) << std::setw(15) << std::left << fraction;
        ss << std::endl;

    
        if(fraction > 0. && fraction <= 2. && histName.back()=='0'){
            count[0][0]++;
            ss00 << histName << std::endl;
        }
        if(fraction > 0. && fraction <= 2. && histName.back()=='1'){
            count[0][1]++;
            ss01 << histName << std::endl;
        }
        if(fraction > 2. && fraction <= 5. && histName.back()=='0'){
            count[1][0]++;
            ss10 << histName << std::endl;
        }
        if(fraction > 2. && fraction <= 5. && histName.back()=='1'){
            count[1][1]++;
            ss11 << histName << std::endl;
        }
        if(fraction > 5. && fraction <= 10. && histName.back()=='0'){
            count[2][0]++;
            ss20 << histName << std::endl;
        }
        if(fraction > 5. && fraction <= 10. && histName.back()=='1'){
            count[2][1]++;
            ss21 << histName << std::endl;
        }
        if(fraction > 10. && histName.back()=='0'){
            count[3][0]++;
            ss30 << histName << std::endl;
        }
        if(fraction > 10. && histName.back()=='1'){
            count[3][1]++;
            ss31 << histName << std::endl;
        }


        ihist++;
    }
    
    ss00 << "Total: " << count[0][0] << std::endl;
    ss01 << "Total: " << count[0][1] << std::endl;
    ss10 << "Total: " << count[1][0] << std::endl;
    ss11 << "Total: " << count[1][1] << std::endl;
    ss20 << "Total: " << count[2][0] << std::endl;
    ss21 << "Total: " << count[2][1] << std::endl;
    ss30 << "Total: " << count[3][0] << std::endl;
    ss31 << "Total: " << count[3][1] << std::endl;

    checkLowEnergyPeakFile << ss00.str() << std::endl;
    checkLowEnergyPeakFile << ss01.str() << std::endl;
    checkLowEnergyPeakFile << ss10.str() << std::endl;
    checkLowEnergyPeakFile << ss11.str() << std::endl;
    checkLowEnergyPeakFile << ss20.str() << std::endl;
    checkLowEnergyPeakFile << ss21.str() << std::endl;
    checkLowEnergyPeakFile << ss30.str() << std::endl;
    checkLowEnergyPeakFile << ss31.str() << std::endl;
    checkLowEnergyPeakFile << ss.str() << std::endl;
    checkLowEnergyPeakFile.close();

    

//                std::cout << "hi2" << std::endl;  
//    for(int iladder = 0; iladder < N_LADDER/2; iladder++) { 
//                std::cout << "hi3" << std::endl;  
//        for(int iva = 0; iva < N_VA; iva++){
//                std::cout << "hi4" << std::endl;  
//            for(int ietareg = 0; ietareg < 2; ietareg++){
//                std::cout << "hi5" << std::endl;  
//                //s << "hVAEnergyX_" << xLadder << "_" << iva << "_" << ietareg;
//                std::cout << "hi6" << std::endl;  
//                myhist = (TH1D*)f->Get(s.str().c_str());
//                std::cout << "hi7" << std::endl;  
//                std::cout << myhist->GetNbinsX() << std::endl;  
//                std::cout << "hi8" << std::endl;  
//            }
//        }
//    }
    
	//outFile->Write();
	//outFile->Close();
	//std::cout << outFileName << " created." << std::endl;
	//sw.Stop();
    //sw.Print();

    return 0;
} // end of main

