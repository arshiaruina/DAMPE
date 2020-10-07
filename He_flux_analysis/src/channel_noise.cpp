#include "../inc/new_va_equalisation.h"

// /beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T001308_20181019T002610_00000_jEwJChikPnVApmeJU6og/DAMPE_2A_OBS_20181019_20181019T001308_20181019T002610_00000_jEwJChikPnVApmeJU6og.root

//int main(int argc, char** argv) {
int main(){

    bool hasRunMetadataTree = true;
    int const nSTKchannels = 384; 
    int const nSTKladders  = 192; 
    int const nSTKTRBs  = 8; 
    //double Ped_calib[nSTKchannels][nSTKladders];
    double Sig_calib[nSTKchannels][nSTKladders];

    //std::string inFileName = argv[1];
    std::string inFileName = "/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T001308_20181019T002610_00000_jEwJChikPnVApmeJU6og/DAMPE_2A_OBS_20181019_20181019T001308_20181019T002610_00000_jEwJChikPnVApmeJU6og.root";
    TFile *f = new TFile(inFileName.c_str());
    TTree *calib_tree;
    if(hasRunMetadataTree) calib_tree = (TTree*) f->Get("RunMetadataTree");

    //if(hasRunMetadataTree && !isMC) {
    if(hasRunMetadataTree) {
    

        TClonesArray* stk_ladder_sigmas      = new TClonesArray("DmpStkLadderCalibration");
        TClonesArray* stk_ladder_pedestals   = new TClonesArray("DmpStkLadderCalibration");
        TClonesArray* stk_ladder_badchannels = new TClonesArray("DmpStkLadderCalibration");
        
        TBranch *b_StkLadderSigmas;   //!
        TBranch *b_StkLadderPedestals;   //!
        TBranch *b_StkLadderBadchannels;   //!
	
        calib_tree->SetBranchAddress("StkLadderSigmas",     &stk_ladder_sigmas,       &b_StkLadderSigmas);
    	calib_tree->SetBranchAddress("StkLadderPedestals",  &stk_ladder_pedestals ,   &b_StkLadderPedestals);
    	calib_tree->SetBranchAddress("StkLadderBadchannels",&stk_ladder_badchannels , &b_StkLadderBadchannels);
    	cout << " entries in RunMetadataTree " << calib_tree->GetEntries() << endl;
        
        if(calib_tree->GetEntries() != 1) {
        	cout << " calib_tree->GetEntries() = " << calib_tree->GetEntries() << ", not 1, stop!" << endl;
        	exit (EXIT_FAILURE);
        }
        
        for (Int_t i=0;i<calib_tree->GetEntries();i++) {
    	
    		b_StkLadderSigmas      ->GetEntry(i);    
    		//b_StkLadderPedestals   ->GetEntry(i);    
    		b_StkLadderBadchannels ->GetEntry(i);
    	
    		for (int ilad = 0 ; ilad <stk_ladder_sigmas->GetLast()+1 ;ilad++){
    	
    			DmpStkLadderCalibration* sigma = (DmpStkLadderCalibration*) (stk_ladder_sigmas->ConstructedAt(ilad));
    			//int trbID    = sigma-> GetTRBNumber();  
    			int ladderID = sigma->GetLadderID();  
    			//cout << " sigma: ladder/ladderID = " << ilad << "/" << ladderID << endl;
    	
    			for(int ich = 0; ich <nSTKchannels; ich++) {
    	
    				double noise = sigma->GetValue(ich);
    				int ibin = ich+ilad*384+1;
    				//h_sigmaCalib->SetBinContent(ibin, noise);
    	
    				if(ibin<=0) cout << "ladder = " <<  ladderID << ", noise " << noise << endl;
    				Sig_calib[ich][ilad] = noise;
    				
    				if(noise <= 0) {
    					cout << " wrong STK noise for ladder " << ladderID << ", channel " << ich << ", noise " << noise << endl;
    					exit (EXIT_FAILURE);
    				}
    				//if(ilad==1 && ich == 10) cout << "ladder = " <<  ladderID << ", channel " << ich << ", noise " << noise << endl;
    			}   
    		}
    	}
    }

    std::cout << "Worked, I guess?" << std::endl;   
    //then in channel loop
    
    //unsigned short noise = 0;
    //if(isMC) noise = 30;
    //else noise = (unsigned short) 10*Sig_calib[ich][ilad];
    //if(noise>255) noise = 255; //avoid overflow

}
