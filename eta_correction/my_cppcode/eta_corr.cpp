#ifndef _ETA
#define _ETA

#include <string>
#include <stdio.h>
#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <numeric>

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
    VaAnalysis(const char * filename,
               const char * option) :
        DmpAnalysis(filename, option) 
    {
        mEventSelector = new DmpEventSelector();

        if(access( "bad_chan_201801.txt", F_OK ) != -1) {
            mTrackSelector = new DmpTrackSelector("bad_chan_201801.txt");
            mTrackSelector->readBadChannelsFile();
//            mTrackSelector->addSelect(DmpTrackSelector::stk_bad_channel);
        }
        else
            mTrackSelector = new DmpTrackSelector();

	/*
	 * TrackSelector helps selecting the STK tracks
	 * General syntax is:
	 *     1.
	 *     mTrackSelector->setSelectionParameter(value);
	 *     to setup the selection threshold
	 *     2.
	 *     mTrackSelector->addSelect(DmpTrackSelector::selection_name);
	 *     to actually add the selection to the list
	 */

        mTrackSelector->addSelect(DmpTrackSelector::psd_match);
        mTrackSelector->setXYnonoverlaps(2); // max # of X-Y non-overlaps
        mTrackSelector->addSelect(DmpTrackSelector::stk_xy_overlap);
        mTrackSelector->addSelect(DmpTrackSelector::stk_missing_impact_point);
        // mTrackSelector->addSelect(DmpTrackSelector::stk_no_1strip_clusters);
        mTrackSelector->setMaxAngStk2Bgo(0.15); // max angular distance between STK track and BGO shower
        mTrackSelector->setMaxDistStk2Bgo(20.); // max distance between STK track and BGO shower
        mTrackSelector->addSelect(DmpTrackSelector::stk_bgo_dist_high);
        mTrackSelector->setSTKtrackChi2Max(15.);
        mTrackSelector->addSelect(DmpTrackSelector::stk_chi2_cut);
        // mTrackSelector->setMinNXYhits(6);
        // mTrackSelector->addSelect(DmpTrackSelector::stk_nXYhits);
        // mTrackSelector->setSTKProtonMipE(60.);
        // mTrackSelector->setSTKXProtonMipRange(0.1);
        // mTrackSelector->setSTKYProtonMipRange(0.1);
        // mTrackSelector->addSelect(DmpTrackSelector::proton_mip);


	/*
	 * EventSelector is used to trigger the events
	 * syntaxis is similar to that of the TrackSelector
	 * Also you can add the trackSelector instance to the EventSelector
	 * to have at least one STK track that corresponds to your TrackSelector
	 */


        mEventSelector->setTrackSelector(mTrackSelector);
        mEventSelector->addSelect(DmpEventSelector::bgo_skim_cuts);
        mEventSelector->addSelect(DmpEventSelector::het);
        mEventSelector->addSelect(DmpEventSelector::has_STK_track);
        // mEventSelector->addSelect(DmpEventSelector::has_only_one_STK_track);

	mTree = NULL;
    }

    ~VaAnalysis() {
	mOutputFile->cd();
	mEventSelector->hSelect()->Write();
	closeOutputFile();
	delete mEventSelector;
    }

    void addTTree() {
	mTree = new TTree("T", "Cluster information tree");

	mTree->Branch("ClustETot", &mClustETot, "mClustETot/F"); // cluster total signal, no theta correction
	mTree->Branch("CosTheta", &mCosTheta, "mCosTheta/F"); // cos(theta) of the track

	mTree->Branch("ClustFS", &mClustFS, "mClustFS/I"); // cluster first strip
	mTree->Branch("ClustLS", &mClustLS, "mClustLS/I"); // cluster last strip
    mTree->Branch("ClustE", &mClustE, "mClustE[10]/F"); // array of strip energies in cluster, max 10 values

	mTree->Branch("ClustEta", &mClustEta, "mClustEta/F"); // cluster eta
	mTree->Branch("ClustVA", &mClustVA, "mClustVA/I"); // cluster VA
 
	mTree->Branch("ThetaX", &mThetaX, "mThetaX/F"); // theta_x at a given point on the track
	mTree->Branch("ThetaY", &mThetaY, "mThetaY/F"); // theta_y at a given point on the track
	mTree->Branch("Direction", &mDirection, "mDirection[3]/F"); // global inclination of the track, x, y and z components
	
    mTree->SetDirectory(0);
    } 

    void createHists() {
	;
    }

    float radianToDegree(float rad){
        return (rad * 180./M_PI);
    }

    void analyseOneEvent(DmpEvent * pev)
    {   
        // Event selection
	    mSelected = mEventSelector->selected(pev);
	    if (!mSelected) {
	        return;
	    }

        // STK track
	    DmpStkTrack * stktrack = mEventSelector->stkTrack();
	    TClonesArray * stkclusters = pev->GetStkSiClusterCollection();

        // Loop over the STK clusters
	    // cout << "Before cluster loop" << endl;
	    
	    TVector3 direction = stktrack->getDirection();
	    mCosTheta = direction.CosTheta();

        for(int icomp=0; icomp<3; icomp++){
	        mDirection[icomp] = direction(icomp);
        }
        
        for (int ipoint=0; ipoint<stktrack->GetNPoints(); ipoint++) {
            for (int ixy=0; ixy<2; ixy++) {
                // Check if cluster is associated to a hit
	    	    DmpStkSiCluster* cluster;
	    	    if(ixy == 0)
	    	        cluster = stktrack -> GetClusterX(ipoint, stkclusters);
	    	    else
	    	        cluster = stktrack -> GetClusterY(ipoint, stkclusters);
	    	    if(!cluster) continue;

                mThetaX = radianToDegree(atan(stktrack->getSlopeX(ipoint)));
                mThetaY = radianToDegree(atan(stktrack->getSlopeY(ipoint)));

                // get cluster energy
                for(int i = 0; i < 10; i++) {
                    mClustE[i] = 0.;
                }
                for(int i = 0; i < cluster->getNstrip(); i++) {
                    if (i < 10) {
                        mClustE[i] = cluster->GetSignal(i);
                    }
                }

	    	    //ladder num
	    	    int ilad = cluster->getLadderHardware();
	    	    //VA num
	    	    int iva = -1;
	    	    mClustFS = cluster->getFirstStrip();
	    	    mClustLS = cluster->getLastStrip();
	    	    for(int i = 0; i < NVA; i++) {
	    	        if(mClustFS >= i * 64 && mClustLS < (i + 1) * 64) {
	    	    	iva = i;
	    	    	break;
	    	        }
	    	    }
	    	    mClustVA = iva;
	    	    if (iva > NVA || iva < 0) continue;

	    	    // Get cluster energy
	    	    mClustETot = cluster->getEnergy();

	    	    // Get eta
	    	    mClustEta = mEtaCorr.calcEta(cluster);

	    	    // mSelected is a special variable. Set it to false to do not call mTree->Fill() again at the end of the event loop
	    	    mSelected=false;

	    	    // set mInterestingEvent to true for the events which you want to see on the event display
	    	    mInterestingEvent=false;

	    	    // Fill the tree for each cluster separately
	    	    this->mTree->Fill();
	        }
	    }
    }

    private :
    
    bool mMC;

    float mClustETot;
    float mCosTheta;
    int mClustFS;
    int mClustLS;
    float mClustE[10];
    float mClustEta;
    int mClustVA;
    float mThetaX;
    float mThetaY;
    float mDirection[3];

    DmpEventSelector * mEventSelector;
    DmpTrackSelector * mTrackSelector;

    EtaCorr mEtaCorr;
};

/*TO DO list :
 *
 * make 2D distributions of the Cluster charge [ADC] vs eta in bins of theta
 * make slices of the 2D Cluster charge [ADC] vs eta with respect to eta bin width of 0.025
 * scan all eta-theta phase space, and apply LanGau fits for the proton peak selection, extract proton peak position
 * repeat the same for helium
 * plot proton peak position profiles with as function of eta in bins of theta
 * fit your plot proton peak position profiles with two linear functions (or find another function), extract constants of this fits
 * Repeat the same for helium, make sure helium peak is 4 times higher than proton
 * use the extracted constants to normalize and correct the proton peak position (we are expecting the flat profile distributions of the proton peaks afterwards)
 * Study the correction in bins of BGO energy
 * Study the correction for the proton and helium MC samples.
 *
 * More immediate todo list :
 *
 * Add following branches in the cppcode TTree:
 * - theta_xy (separate values for x and y)
 *   - psdcharge (can be used to select tracks)
 *   - bgo_energy
 *
 * */


#endif
