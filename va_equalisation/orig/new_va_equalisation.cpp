/* Code to compute eta values
 * Author: Arshia Ruina
 * To compile, g++ -Wall `root-config --cflags --libs` -I$DMPSWSYS/include -L$DMPSWSYS/lib -lDmpEvent -lDmpService eta.cc -o ../submit/eta.exe
 * */

// example file
// /beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T001308_20181019T002610_00000_jEwJChikPnVApmeJU6og/DAMPE_2A_OBS_20181019_20181019T001308_20181019T002610_00000_jEwJChikPnVApmeJU6og.root

#include "../inc/new_va_equalisation.h"
//#include "../inc/track_selection.hpp"

// iladder          0 - 47          48 - 95
// X ladder IDs     48 - 95         144 - 191
// Y ladder IDs     0 - 47          96 - 143

int main(int argc, char** argv) {
//int main() {

	TStopwatch sw;
	sw.Start();

	std::string week = argv[1];
	//std::string inFileName = "/beegfs/users/ruina/VAequalisation/testNtuples/merged_180620.root";
	std::string inFileName = "/beegfs/users/ruina/VAequalisation/ntuples_201801_201802/merged/" + week + ".root";
	//std::string inFileName = "/beegfs/users/ruina/VAequalisation/ntuples_201801_201802/DAMPE_2A_OBS_20180102_20180102T064154_20180102T065209_00000_px5VV00AObfsPk9dsoKK.root"; // for example
    TFile *f = new TFile(inFileName.c_str());
    //f->ls();
	
    /*//---------------------------------------------//

    std::stringstream ss;
    std::ofstream badChannelList;
    badChannelList.open("/beegfs/users/ruina/VAequalisation/out/badChannelList.txt");

    ss << "Ladder \t Channel" << std::endl;
    TH2I *hNoise = (TH2I*)f->Get("hNoiseInfo");
    for(int iladd=0; iladd<192; iladd++){
        for(int ich=0; ich<384; ich++){
            if(hNoise->GetBinContent(ich+iladd*386)){
                //ss << iladd << "\t" << ich << "\t" << hNoise->GetBin(ich,iladd) << "\t" << ich+iladd*386 << "\t" << hNoise->GetBinContent(ich+iladd*386) << std::endl;
                ss << iladd << "\t" << ich << "\t" << std::endl;
            }
        }
    }
    //std::cout << "x axis bins, channels: " << hNoise->GetXaxis()->GetNbins() << std::endl;
    //std::cout << "y axis bins, ladders: " << hNoise->GetYaxis()->GetNbins() << std::endl;

    badChannelList << ss.str();
    badChannelList.close();   

    sw.Stop();
    sw.Print(); 
    return 0; 
    //---------------------------------------------//
    */


    TTree *t = (TTree*)f->Get("MySelectionTree");
    //t->Print();

	//std::size_t found = inFileName.find_last_of("/");
    
    //bool flagAppCorrFac = true;   //for corrFac application
    bool flagAppCorrFac = false;    //without corrFac application
    
    TFile *inFileCorrFac = new TFile(inFileNameCorrFac.c_str());
    hCorrFac = (TH2D*)inFileCorrFac->Get(hCorrFacName.c_str());

    int nClustersX = 0;
    int nClustersY = 0;
    int nSelTracks = 0;
    std::vector <float> *clusterEnergy = 0;
    std::vector <float> *clusterEnergyAdc = 0;
    std::vector <float> *clusterEta = 0;
    std::vector <int> *clusterEtaRegion = 0;
    std::vector <int> *clusterLadder = 0;
    std::vector <int> *clusterVA = 0;
    std::vector <int> *clusterFirstStrip = 0;
    std::vector <int> *clusterLastStrip = 0;

    t->SetBranchAddress("nClustersX", &nClustersX);
    t->SetBranchAddress("nClustersY", &nClustersY);
    t->SetBranchAddress("nSelTracks", &nSelTracks);
    t->SetBranchAddress("clusterEnergy", &clusterEnergy);
    t->SetBranchAddress("clusterEnergyAdc", &clusterEnergyAdc);
    t->SetBranchAddress("clusterEta", &clusterEta);
    t->SetBranchAddress("clusterEtaRegion", &clusterEtaRegion);
    t->SetBranchAddress("clusterLadder", &clusterLadder);
    t->SetBranchAddress("clusterVA", &clusterVA);
    t->SetBranchAddress("clusterFirstStrip", &clusterFirstStrip);
    t->SetBranchAddress("clusterLastStrip", &clusterLastStrip);

    //outFileName = "tmp/pre_correction_" + inFileName.substr(found+1);
    //outFileName = "tmp/201801_201802/pre_correction.root";
    //outFileName = "/beegfs/users/ruina/VAequalisation/ntuples_201801_201802/post_correction_oldfactors.root";
    //outFileName = "/beegfs/users/ruina/VAequalisation/ntuples_201801_201802/pre_correction_oldfactors.root";
    //outFileName = "/beegfs/users/ruina/VAequalisation/test_pre_correction_oldfactors.root";
    outFileName = "/beegfs/users/ruina/VAequalisation/pre_correction_oldfactors_201801_201802/test/" + week + ".root";
    TFile *outFile = new TFile(outFileName.c_str(), "RECREATE");
   
    gStyle->SetOptStat(0);
 
    for(int iladder = 0; iladder < N_LADDER; iladder++) { 
        hEnergyAdcLadders[iladder] = new TH1F(Form("hEnergyAdcLadder_%d",iladder),Form("Energy for ladder %d",iladder),200,0.,200.);
        hEnergyAdcLadders[iladder] -> GetXaxis() -> SetTitle("ADC counts");
        hEnergyAdcLadders[iladder] -> GetYaxis() -> SetTitle("No. of events");
        /*for(int iva = 0; iva < N_VA; iva++){
            for(int ietareg = 0; ietareg < 2; ietareg++){
                hVAEnergy[iladder][iva][ietareg] = new TH1F(Form("hVAEnergy_%d_%d_%d",iladder,iva,ietareg),Form("Energy for ladder %d VA %d #eta region %d",iladder,iva,ietareg),200,0.,200.);
                hVAEnergy[iladder][iva][ietareg] -> GetXaxis() -> SetTitle("ADC counts");
                hVAEnergy[iladder][iva][ietareg] -> GetYaxis() -> SetTitle("No. of events");
            }
        }*/
    }

    hCumulEnergyAdc = new TH1F("ADC for all VAs","ADC for all VAs",200,0.,200.);
    hCumulEnergyAdc -> GetXaxis() -> SetTitle("ADC counts");
    hCumulEnergyAdc -> GetYaxis() -> SetTitle("No. of events");
    hCumulEnergy = new TH1F("Energy for all VAs","Energy for all VAs",200,0.,200.);
    hCumulEnergy -> GetXaxis() -> SetTitle("Energy");
    hCumulEnergy -> GetYaxis() -> SetTitle("No. of events");

    int nEntries = t->GetEntries();
    std::cout << "[INFO] Number of entries " << nEntries << std::endl; 

    std::cout << "[INFO] variables set, starting loop over entries..." << std::endl;


    //Access bad channel list file
    std::ifstream badChannelList;
    int badChannel, badChannelLad; 


	for (int i = 0; i < nEntries; i++){ // uncomment for analysis run
    //for (int i = 0; i < 5000; i++){ // uncomment for debug run

		float progress = 100.0 * ((float) i) / ((float) nEntries);
        std::cout << std::setprecision(3) << "[ " << progress << " % ] \r";
    
		t->GetEntry(i);

        //int nClustersTotal = nClustersX + nClustersY;
        int nClustersTotal = 12; //TODO: should be changed to above
        //std::cout << "Entry: " << i << std::endl;
        //std::cout << "nSelTracks: " << nSelTracks << std::endl;
        //std::cout << "nClustersX: " << nClustersX << std::endl;
        //std::cout << "nClustersY: " << nClustersY << std::endl;
        ////std::cout << "nClustersTotal: " << nClustersTotal << std::endl;
        //std::cout << "Size of clusterEnergy: " << clusterEnergy->size() << std::endl;
        //std::cout << "Size of clusterEnergyAdc: " << clusterEnergyAdc->size() << std::endl;
        //std::cout << "Size of clusterEta: " << clusterEta->size() << std::endl;
        //std::cout << "Size of clusterEtaRegion: " << clusterEtaRegion->size() << std::endl;
        //std::cout << "Size of clusterLadder: " << clusterLadder->size() << std::endl;
        //std::cout << "Size of clusterVA: " << clusterVA->size() << std::endl;


        for(int index = 0; index < nClustersTotal; index++){

            //std::cout << "clusterLadder " << index << " " << clusterLadder->at(index) << std::endl;
            //std::cout << "clusterVA " << index << " " << clusterVA->at(index) << std::endl;
            //std::cout << "clusterEtaRegion " << index << " " << clusterEtaRegion->at(index) << std::endl;

            if(clusterVA->at(index) < 0 || clusterEtaRegion->at(index) < 0) continue;
            //TODO: only for eta region 0 for now...
            if(clusterEtaRegion->at(index)==1) continue;

            bool clusterContainsBadChannel = false;
            badChannelList.open("/beegfs/users/ruina/VAequalisation/out/badChannelList_correct.txt");
            while ( !badChannelList.eof() ) { // keep reading until end-of-file
                badChannelList >> badChannelLad >> badChannel;
                if(clusterLadder->at(index) == badChannelLad && badChannel > clusterFirstStrip->at(index)-1 && badChannel < clusterLastStrip->at(index)+1){
                    //cluster contains noisy channel
                    clusterContainsBadChannel = true;
                    break;
                }
            }
            badChannelList.close();

            if(clusterContainsBadChannel) continue;

            if (flagAppCorrFac){
                //hVAEnergy[clusterLadder->at(index)][clusterVA->at(index)][clusterEtaRegion->at(index)] -> Fill(clusterEnergyAdc->at(index) * hCorrFac->GetBinContent(clusterLadder->at(index)+1,clusterVA->at(index)+1));
                hCumulEnergyAdc -> Fill(clusterEnergyAdc->at(index) * hCorrFac->GetBinContent(clusterLadder->at(index)+1,clusterVA->at(index)+1));
                hCumulEnergy -> Fill(clusterEnergy->at(index) * hCorrFac->GetBinContent(clusterLadder->at(index)+1,clusterVA->at(index)+1));
            }
            else{
                //hVAEnergy[clusterLadder->at(index)][clusterVA->at(index)][clusterEtaRegion->at(index)] -> Fill(clusterEnergyAdc->at(index));        
                hCumulEnergyAdc -> Fill(clusterEnergyAdc->at(index));
                hCumulEnergy -> Fill(clusterEnergy->at(index));
    	    }
        
            //To check which VAs are causing the spikes...
            hEnergyAdcLadders[clusterLadder->at(index)] -> Fill(clusterEnergyAdc->at(index));
            
        }
    } 

    for(int ican = 0; ican < 16; ican++){
        int i1 = ican*12;
        int i2 = (ican+1)*12;
        myCanvas[ican] = new TCanvas(Form("cEnergyAdcLAdders_%d",ican),Form("Energy for ladders %d to %d",i1,i2-1),800,600);
        //auto myFrame = myCanvas[ican]->DrawFrame(-1., 0., 2., 2.);
        //myFrame->SetTitle("My Title");
        //Colours of overlaid histograms 
        int icol = 0;
        int colours[] = {1,2,3,4,5,6,7,8,9,11,12,29};
        //Replace title on canvas
        //gStyle->SetOptTitle(0);
        //TPaveLabel *title = new TPaveLabel(.11,.95,.35,.99,"new title","brndc");
        //title->Draw();
        //myCanvas[ican]->Update();
        //Change hist ymax
        float yMax = hEnergyAdcLadders[i1] -> GetMaximum();
        //Add legend on canvas
        TLegend *legend = new TLegend(0.75,0.55,0.88,0.88);
        for(int ilad = i1; ilad < i2; ilad++){
            if(yMax < hEnergyAdcLadders[ilad] -> GetMaximum())
                yMax = hEnergyAdcLadders[ilad] -> GetMaximum();
            hEnergyAdcLadders[i1] -> SetMaximum(1.1*yMax);
            hEnergyAdcLadders[i1] -> SetTitle(Form("Energy for ladders %d to %d",i1,i2-1));
            hEnergyAdcLadders[ilad] -> SetLineColor(colours[icol]);
            hEnergyAdcLadders[ilad] -> SetLineWidth(2);
            if(ilad%12==0)
                hEnergyAdcLadders[ilad] -> Draw();
            else
                hEnergyAdcLadders[ilad] -> Draw("same");
            icol++;
            legend->AddEntry(hEnergyAdcLadders[ilad],Form("ladder %d",ilad),"l");
        }
        legend->SetBorderSize(1);
        legend->SetTextSize(0.02);
        legend->Draw();
        myCanvas[ican]->Write();
    }

	outFile->Write();
	outFile->Close();
    f->Close();
    
    std::cout << "[INFO] " << outFileName << " created." << std::endl;
	sw.Stop();
    sw.Print();
    
} // end of main

