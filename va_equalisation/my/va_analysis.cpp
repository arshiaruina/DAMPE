#ifndef _VA
#define _VA

#include <string>
#include <stdio.h>
#include <cmath>
#include <iostream>
#include <sstream>

#include "analysis.hpp"
#include "definitions.hpp"
#include "event_selection.hpp"
#include "track_selection.hpp"
#include "etacorr.hpp"

//constants
#define NVA     6
#define NLADDER 192

using namespace std;
using namespace myDampeLib;

class VaAnalysis : public DmpAnalysis {
public :
    VaAnalysis(const char * filename, const char * option) : DmpAnalysis(filename, option) {

        mEventSelector = new DmpEventSelector();

        if(access( "bad_chan_201801.txt", F_OK ) != -1) {
            mTrackSelector = new DmpTrackSelector("bad_chan_201801.txt");
//            mTrackSelector->readBadChannelsFile();
//            mTrackSelector->addSelect(DmpTrackSelector::stk_bad_channel);
        }
        else
            mTrackSelector = new DmpTrackSelector();

//        mTrackSelector->addSelect(DmpTrackSelector::psd_match);
//        mTrackSelector->setXYnonoverlaps(0);//2); // max # of X-Y non-overlaps
//        mTrackSelector->addSelect(DmpTrackSelector::stk_xy_overlap);
//        mTrackSelector->addSelect(DmpTrackSelector::stk_missing_impact_point);
//        mTrackSelector->addSelect(DmpTrackSelector::stk_no_1strip_clusters);
//        mTrackSelector->setMaxAngStk2Bgo(0.15); // max angular distance between STK track and BGO shower
//        mTrackSelector->setMaxDistStk2Bgo(20.); // max distance between STK track and BGO shower
//        mTrackSelector->addSelect(DmpTrackSelector::stk_bgo_dist_high);
        mTrackSelector->setSTKtrackChi2Max(15.);
        mTrackSelector->addSelect(DmpTrackSelector::stk_chi2_cut);
        mTrackSelector->setMinNXYhits(6);
        mTrackSelector->addSelect(DmpTrackSelector::stk_nXYhits);
        mTrackSelector->setSTKProtonMipE(55.);
        mTrackSelector->setSTKXProtonMipRange(0.1);
        mTrackSelector->setSTKYProtonMipRange(0.1);
        mTrackSelector->addSelect(DmpTrackSelector::proton_mip);
        mEventSelector->setTrackSelector(mTrackSelector);
//      mEventSelector->addSelect(DmpEventSelector::bgo_skim_cuts);
//      mEventSelector->addSelect(DmpEventSelector::het);
        mEventSelector->addSelect(DmpEventSelector::has_STK_track);
        //mEventSelector->addSelect(DmpEventSelector::has_only_one_STK_track);

        mCorrection0.resize(NLADDER);
        mCorrection1.resize(NLADDER);
        
        for (int ilad=0; ilad<NLADDER; ilad++) {
	        mCorrection0[ilad].resize(NVA);
	        mCorrection1[ilad].resize(NVA);
	        for (int iva=0; iva<NVA; iva++) {
	            mCorrection0[ilad][iva] = 1.;
	            mCorrection1[ilad][iva] = 1.;
	        }
	    }

        mTree = NULL;
    }   

    ~VaAnalysis() {
        mOutputFile -> cd();
	 
        mEtaVsEnergy->Write();    

        for (int i=0; i<mVAhists.size(); i++) 
            mVAhists[i]->Write();
	
        //mBadChanClustE->Write();
	    //mGoodChanClustE->Write();
	    mAllChanClustE->Write();
	
        //mBadChanClustEnergy->Write();
	    //mGoodChanClustEnergy->Write();
	    mAllChanClustEnergy->Write();
	
        //mBadChanClustADC->Write();
	    //mGoodChanClustADC->Write();
	    mAllChanClustADC->Write();
	
        //mBadChanClustEcosT->Write();
	    //mGoodChanClustEcosT->Write();
	    mAllChanClustEcosT->Write();
	
        //mBadChanClustEnergycosT->Write();
	    //mGoodChanClustEnergycosT->Write();
	    mAllChanClustEnergycosT->Write();
	
        //mBadChanClustADCcosT->Write();
	    //mGoodChanClustADCcosT->Write();
	    mAllChanClustADCcosT->Write();
	
        mMultStripClustE->Write();
	    m1StripClustE->Write();
	    cout << "Number of selected events = " << mSMeanX->GetEntries() << endl;
        //cout << "vertical events: " << iVertical << endl;
	    mSMeanX->Write();
	    mSMeanY->Write();
	    mNStrip->Write();
	    mEta->Write();
	    mCosTheta->Write();
        mEventSelector->hSelect()->Write();
	    closeOutputFile();
	    delete mEventSelector;
    }

    void addTTree() {
	    
        mTree = new TTree("T", "Cluster information tree");

	    mTree->Branch("ClustFS", &mClustFS, "mClustFS/I"); // cluster first strip
	    mTree->Branch("ClustLS", &mClustLS, "mClustLS/I"); // cluster last strip
        mTree->Branch("ClustE", &mClustE, "mClustE[10]/F"); // array of strip energies in cluster, max 10 values
         
	    mTree->SetDirectory(0);
    } 

    void createHists() {
	// Histograms for each VA
     
        mEtaVsEnergy = new TH2F("mEtaVsEnergy", "Eta vs. Cluster Energy;Eta;STK cluster energy",50,0.5,1.,150,0.,150.);
   
        mVAhists.resize(NLADDER * NVA * 2);
	    for (int ilad = 0; ilad < NLADDER; ilad++) {
	        for (int iva = 0; iva < NVA; iva++) {
	            for (int ieta = 0; ieta <= 1; ieta++) {
	                stringstream name;
	                name << "hVAEnergy_" << ilad << "_" << iva << "_" << ieta;
                    mVAhists[ilad * NVA * 2 + iva * 2 + ieta] = new TH1F(name.str().c_str(), name.str().c_str(), 100, 0., 100.);
	            }
	        }
        }

    	//mGoodChanClustE = new TH1F("mGoodChanClustE", "Energy of clusters with good channels only;STK cluster energy;Counts", 100, 0., 100.);
	    //mBadChanClustE = new TH1F("mBadChanClustE", "Energy of clusters with at least one bad channel;STK cluster energy;Counts", 100, 0., 100.);
	    mAllChanClustE = new TH1F("mAllChanClustE", "Energy of clusters with all channels;STK cluster energy;Counts", 100, 0., 100.);
    	
    	//mGoodChanClustEcosT = new TH1F("mGoodChanClustEcosT", "EnergycosT of clusters with good channels only;STK cluster energycosT;Counts", 100, 0., 100.);
	    //mBadChanClustEcosT = new TH1F("mBadChanClustEcosT", "EnergycosT of clusters with at least one bad channel;STK cluster energycosT;Counts", 100, 0., 100.);
	    mAllChanClustEcosT = new TH1F("mAllChanClustEcosT", "EnergycosT of clusters with all channels;STK cluster energycosT;Counts", 100, 0., 100.);
    	
    	//mGoodChanClustEnergy = new TH1F("mGoodChanClustEnergy", "My energy of clusters with good channels only;STK cluster energy;Counts", 100, 0., 100.);
	    //mBadChanClustEnergy = new TH1F("mBadChanClustEnergy", "My energy of clusters with at least one bad channel;STK cluster energy;Counts", 100, 0., 100.);
	    mAllChanClustEnergy = new TH1F("mAllChanClustEnergy", "My energy of clusters with all channels;STK cluster energy;Counts", 100, 0., 100.);
    	
    	//mGoodChanClustEnergycosT = new TH1F("mGoodChanClustEnergycosT", "My energycosT of clusters with good channels only;STK cluster energycosT;Counts", 100, 0., 100.);
	    //mBadChanClustEnergycosT = new TH1F("mBadChanClustEnergycosT", "My energycosT of clusters with at least one bad channel;STK cluster energycosT;Counts", 100, 0., 100.);
	    mAllChanClustEnergycosT = new TH1F("mAllChanClustEnergycosT", "My energycosT of clusters with all channels;STK cluster energycosT;Counts", 100, 0., 100.);
    	
        //mGoodChanClustADC = new TH1F("mGoodChanClustADC", "ADC of clusters with good channels only;STK cluster ADC;Counts", 100, 0., 100.);
	    //mBadChanClustADC = new TH1F("mBadChanClustADC", "ADC of clusters with at least one bad channel;STK cluster ADC;Counts", 100, 0., 100.);
	    mAllChanClustADC = new TH1F("mAllChanClustADC", "ADC of clusters with all channels;STK cluster ADC;Counts", 100, 0., 100.);

        //mGoodChanClustADCcosT = new TH1F("mGoodChanClustADCcosT", "ADCcosT of clusters with good channels only;STK cluster ADCcosT;Counts", 100, 0., 100.);
	    //mBadChanClustADCcosT = new TH1F("mBadChanClustADCcosT", "ADCcosT of clusters with at least one bad channel;STK cluster ADCcosT;Counts", 100, 0., 100.);
	    mAllChanClustADCcosT = new TH1F("mAllChanClustADCcosT", "ADCcosT of clusters with all channels;STK cluster ADCcosT;Counts", 100, 0., 100.);

	    mMultStripClustE = new TH1F("mMultStripClustE", "Energy of multi-strip clusters;STK cluster energy;Counts", 100, 0., 100.);
	    m1StripClustE = new TH1F("m1StripClustE", "Energy of single-strip clusters;STK cluster energy;Counts", 100, 0., 100.);

	    mSMeanX = new TH1F("mSMeanX", "S_mean, XZ plane;STK cluster energy;Counts", 1000, 0., 10.);
	    mSMeanY = new TH1F("mSMeanY", "S_mean, YZ plane;STK cluster energy;Counts", 1000, 0., 10.);

	    mNStrip = new TH1F("mNStrip", "Number of strips;Number of strips;Counts", 100, 0., 100.);

	    mEta = new TH1F("mEta", "eta distribution;eta;Counts", 50, 0.5, 1.);
	    mCosTheta = new TH1F("mCosTheta", "cos theta distribution;cos theta;Counts", 100, -1, 1.);
    }

    void setCorr(const char * corrfile) {
    
        TFile * corrf = new TFile(corrfile, "READ");
	    TH2D * hcorr0 = (TH2D * ) corrf->Get("hCorrFac0");
	    TH2D * hcorr1 = (TH2D * ) corrf->Get("hCorrFac1");
        for (int ilad=0; ilad<NLADDER; ilad++) {
            for (int iva=0; iva<NVA; iva++) {
                mCorrection0[ilad][iva] = hcorr0->GetBinContent(ilad + 1, iva + 1);
                mCorrection1[ilad][iva] = hcorr1->GetBinContent(ilad + 1, iva + 1);
            }
        }
	    corrf->Close();
    }

    void analyseOneEvent(DmpEvent * pev) {   
    // Event selection
	
        mSelected = mEventSelector->selected(pev);
	    if (!mSelected) {
	        return;
	    }

        // STK track
	    DmpStkTrack * stktrack = mEventSelector->stkTrack();
	    TClonesArray * stkclusters = pev->GetStkSiClusterCollection();
        TClonesArray* stkladderadc = pev->GetStkLadderAdcCollection();

        // Build histogram of S_mean for X and Y
	    float s_mean_xy[2];
	    mTrackSelector->getSMean(stktrack, pev, s_mean_xy);
//        cout << s_mean_xy[0] << " " << s_mean_xy[1] << endl;
	    mSMeanX -> Fill(s_mean_xy[0]);
	    mSMeanY -> Fill(s_mean_xy[1]);

        // Loop over the STK cluster
	    double costheta = stktrack->getDirection().CosTheta();
        mCosTheta->Fill(costheta);

        //double theta_rad = stktrack->getDirection().Theta();
        //double theta_deg = theta_rad * 180 * (1/3.142);
        //if(fabs(theta_deg) < 10) iVertical++; 

//	    vector<int> layerBadChannel(mTrackSelector->layerBadChannel(stktrack, pev));

	    for (int ipoint=0; ipoint<stktrack->GetNPoints(); ipoint++) {
    
            for (int ixy=0; ixy<2; ixy++) {

                // Check if cluster is associated to a hit
		        DmpStkSiCluster* cluster;
		        float energyAdc = 0.;

                if(ixy == 0)
		            cluster = stktrack -> GetClusterX(ipoint, stkclusters);// my track access is different. TODO: CHECK
		        else
		            cluster = stktrack -> GetClusterY(ipoint, stkclusters);
                if(!cluster) continue;

                // get cluster ADC counts
                for(int istrip = 0; istrip < cluster->getNstrip(); istrip++){
                    energyAdc += cluster->GetAdcValue(istrip,stkladderadc);
                }

                // get cluster energy
		        //float e = cluster -> GetSignalTotal () * costheta; // GetSignalTotal gives the sum of adc counts 
                //float energy = cluster -> getEnergy () * costheta;            
		        float e = cluster -> GetSignalTotal (); // GetSignalTotal gives the sum of adc counts 
                float energy = cluster -> getEnergy ();             
    
		        if(cluster->getNstrip() == 1)
		            m1StripClustE -> Fill(e);
		        else
		            mMultStripClustE -> Fill(e);

		        if(cluster->getNstrip() == 1) continue;

		        mNStrip->Fill(cluster->getNstrip());

		        int ilad = cluster->getLadderHardware();
		        int iva;
		        mClustFS = cluster->getFirstStrip();
		        mClustLS = cluster->getLastStrip();
		        for(int i = 0; i < NVA; i++) {
		            if(mClustFS >= i * 64 && mClustLS < (i + 1) * 64) {
		        	    iva = i;
		        	    break;
		            }
		        }

		        // Fill the mClustE array
		        // First initialize it with zeros
		        for(int i = 0; i < 10; i++) {
		            mClustE[i] = 0.;
		        }
		        for(int i = 0; i < cluster->getNstrip(); i++) {
		            if (i < 10) {
		        	    mClustE[i] = cluster->GetSignal(i);
		            }
		        }

                mAllChanClustE -> Fill(e);
                mAllChanClustEnergy -> Fill(energy);
                mAllChanClustADC -> Fill(energyAdc);
                mAllChanClustEcosT -> Fill(e * costheta);
                mAllChanClustEnergycosT -> Fill(energy * costheta);
                mAllChanClustADCcosT -> Fill(energyAdc * costheta);

                //if (find(layerBadChannel.begin(), layerBadChannel.end(), ipoint) != layerBadChannel.end()){
		        //    mBadChanClustE -> Fill(e);
		        //    mBadChanClustEnergy -> Fill(energy);
                //    mBadChanClustADC -> Fill(energyAdc);
                //}
		        //else {
		        //    mGoodChanClustE -> Fill(e);
		        //    mGoodChanClustEnergy -> Fill(energy);
                //    mGoodChanClustADC -> Fill(energyAdc);
                //}    
		        double eta = mEtaCorr.calcEta(cluster);
		        mEta -> Fill(eta);


                // eta vs. cluster energy
                mEtaVsEnergy->Fill(eta,energyAdc*costheta*mCorrection0[ilad][iva]);
                mEtaVsEnergy->SetOption("COLZ");
                //mEtaVsEnergy->SetContour(20);
                //mEtaVsEnergy->SetStats(0);
                //mEtaVsEnergy->SetMinimum(0.9);
                //mEtaVsEnergy->SetMaximum(1.1);


		        // Define eta region
                int ieta;
		        if (eta > 0.5 && eta < 0.6) // float strip region
		            ieta = 1;
		        else if(eta > 0.8 && eta < 1.) // readout strip region
		            ieta = 0;
		        else continue; // if eta is not suitable for the VA analysis, continue with another cluster

                //cout << "corr fac for ladder " << ilad << " VA " << iva << ": " << mCorrection[ilad][iva] << endl;
		        mVAhists[ilad * NVA * 2 + iva * 2 + ieta] -> Fill(energyAdc * costheta * mCorrection0[ilad][iva]);
		        //if(ieta=0)
                //    mVAhists[ilad * NVA * 2 + iva * 2 + ieta] -> Fill(energyAdc * costheta * mCorrection0[ilad][iva]);
                //if(ieta=1)    
		        //    mVAhists[ilad * NVA * 2 + iva * 2 + ieta] -> Fill(energyAdc * costheta * mCorrection1[ilad][iva]);

		        this->mTree->Fill();
		        // mSelected is a special variable. Set it to false to not call mTree->Fill() again at the end of the event loop
		        mSelected=false;

		        // set mInterestingEvent to true for the events which you want to see on the event display
		        mInterestingEvent=false;
    	
            } // end of ixy loop
	    } // end of ipoint loop
    }

private :
    bool mMC;

    int mClustFS;
    int mClustLS;
    float mClustE[10];
    //int iVertical = 0;

    DmpEventSelector * mEventSelector;
    DmpTrackSelector * mTrackSelector;

    vector<vector<double> > mCorrection0;
    vector<vector<double> > mCorrection1;

    TH2F * mEtaVsEnergy;

    vector<TH1F *> mVAhists;

    //TH1F * mBadChanClustE;
    //TH1F * mGoodChanClustE;
    TH1F * mAllChanClustE;
    
    //TH1F * mBadChanClustEnergy;
    //TH1F * mGoodChanClustEnergy;
    TH1F * mAllChanClustEnergy;
    
    //TH1F * mBadChanClustADC;
    //TH1F * mGoodChanClustADC;
    TH1F * mAllChanClustADC;

    //TH1F * mBadChanClustEcosT;
    //TH1F * mGoodChanClustEcosT;
    TH1F * mAllChanClustEcosT;
    
    //TH1F * mBadChanClustEnergycosT;
    //TH1F * mGoodChanClustEnergycosT;
    TH1F * mAllChanClustEnergycosT;
    
    //TH1F * mBadChanClustADCcosT;
    //TH1F * mGoodChanClustADCcosT;
    TH1F * mAllChanClustADCcosT;

    TH1F * m1StripClustE;
    TH1F * mMultStripClustE;

    TH1F * mSMeanX;
    TH1F * mSMeanY;

    TH1F * mNStrip;

    TH1F * mEta;
    TH1F * mCosTheta;

    EtaCorr mEtaCorr;
};

#endif
