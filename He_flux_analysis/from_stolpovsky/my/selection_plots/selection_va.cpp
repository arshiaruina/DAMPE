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

//constants
#define NVA     6
#define NLADDER 192
#define ETALIM  0.25

using namespace std;
using namespace myDampeLib;

class SelectionAnalysis : public DmpAnalysis {
public :
    SelectionAnalysis(const char * filename,
               const char * option) :
        DmpAnalysis(filename, option) 
    {
        mEventSelector = new DmpEventSelector();

	if(access( "bad_chan.txt", F_OK ) != -1) {
            mTrackSelector = new DmpTrackSelector("bad_chan.txt");
//	    mTrackSelector->readBadChannelsFile();
//	    mTrackSelector->addSelect(DmpTrackSelector::stk_bad_channel);
	}
        else
            mTrackSelector = new DmpTrackSelector();


//	mTrackSelector->addSelect(DmpTrackSelector::psd_match);
        mTrackSelector->setXYnonoverlaps(0);//2); // max # of X-Y non-overlaps
        mTrackSelector->addSelect(DmpTrackSelector::stk_xy_overlap);
        mTrackSelector->addSelect(DmpTrackSelector::stk_missing_impact_point);
//	mTrackSelector->addSelect(DmpTrackSelector::stk_no_1strip_clusters);
        mTrackSelector->setMaxAngStk2Bgo(0.15); // max angular distance between STK track and BGO shower
	mTrackSelector->setMaxDistStk2Bgo(20.); // max distance between STK track and BGO shower
//        mTrackSelector->addSelect(DmpTrackSelector::stk_bgo_dist_high);
        mTrackSelector->setSTKtrackChi2Max(15.); //30.);
//        mTrackSelector->addSelect(DmpTrackSelector::stk_chi2_cut);
	mTrackSelector->setMinNPoints(6);
//	mTrackSelector->addSelect(DmpTrackSelector::stk_npoints);
	mTrackSelector->setSTKProtonMipE(60.);
	mTrackSelector->setSTKXProtonMipRange(0.1);
        mTrackSelector->setSTKYProtonMipRange(0.1);
//	mTrackSelector->addSelect(DmpTrackSelector::proton_mip);

        mEventSelector->setTrackSelector(mTrackSelector);

//	mEventSelector->addSelect(DmpEventSelector::bgo_skim_cuts);
//	mEventSelector->addSelect(DmpEventSelector::het);
	mEventSelector->addSelect(DmpEventSelector::has_STK_track);
//	mEventSelector->addSelect(DmpEventSelector::has_only_one_STK_track);
    }

    ~SelectionAnalysis() {
        mOutputFile -> cd();
	mSMeanX->Write();
	mSMeanY->Write();
	mNonOverlaps->Write();
	mSTK2BGOang->Write();
	mSTK2BGOdist->Write();
	mChi2->Write();
	mNpoints->Write();
	mEventSelector->hSelect()->Write();
	closeOutputFile();
	delete mEventSelector;
    }

    void addTTree() {
	mTree = new TTree("Trig", "Tree");
        mTree->Branch("Trig", &mTrig, "mTrig/I");
	mTree->Branch("TrackTrig", &mTrackTrig, "mTrackTrig/I");
	mTree->Branch("PrimInStk", &mPrimInStk, "mPrimInStk/B");
        mTree->Branch("BGOenergy", &mBGOenergy, "mBGOenergy/D");
        mTree->Branch("Eprim", &mEprim, "mEprim/D");

	mTree->Branch("mNonoverlaps", &mNonoverlaps, "mNonoverlaps/I");
	mTree->Branch("mSmeanX", &mSmeanX, "mSmeanX/F");
	mTree->Branch("mSmeanY", &mSmeanY, "mSmeanY/F");
	mTree->Branch("mAngle", &mAngle, "mAngle/F");
	mTree->Branch("mChisqndof", &mChisqndof, "mChisqndof/F");
	mTree->Branch("mNpts", &mNpts, "mNpts/I");

        mTree->SetDirectory(0);
    }

    void createHists() {
	mSMeanX = new TH1F("mSMeanX", "mSMeanX", 200, 0., 2.);
	mSMeanY = new TH1F("mSMeanY", "mSMeanY", 200, 0., 2.);
	mNonOverlaps = new TH1F("mNonOverlaps", "mNonOverlaps", 10, 0., 10.);
	mSTK2BGOang = new TH1F("mSTK2BGOang", "mSTK2BGOang", 100, 0., 0.5);
	mSTK2BGOdist = new TH1F("mSTK2BGOdist", "mSTK2BGOdist", 500, 0., 100.);
	mChi2 = new TH1F("mChi2", "mChi2", 200, 0., 100.);
	mNpoints = new TH1F("mNpoints", "mNpoints", 10, 0., 10.);
    }

    void analyseOneEvent(DmpEvent * pev)
    {   
        // Get primary energy
        if (mMC) {
            DmpEvtSimuPrimaries * prim = pev->pEvtSimuPrimaries();
            // float mEprim = (float)(prim->pvpart_ekin);

            // Does the primary track pass through STK?
            float x = prim->pv_x;
            float y = prim->pv_y;
            float z = prim->pv_z;

            float px = prim->pvpart_px;
            float py = prim->pvpart_py;
            float pz = prim->pvpart_pz;

            float stk_top_z = -210.0319976349; // z coord of the very first stk layer (farthest from bgo)
            float stk_bot_z = -41.96795840450001; // z coord of the last stk layer
	    float bgo_bot_z = 448; // last layer of BGO

            float x_top = x + (stk_top_z - z) * px/pz;
            float y_top = y + (stk_top_z - z) * py/pz;
            float x_bot = x + (bgo_bot_z - z) * px/pz;
            float y_bot = y + (bgo_bot_z - z) * py/pz;

            float xy_lim_top = 400.;
	    float xy_lim_bot = 280.;

            mPrimInStk = abs(x_top) < xy_lim_top && abs(x_bot) < xy_lim_bot &&
                         abs(y_top) < xy_lim_top && abs(y_bot) < xy_lim_bot;
        }
        else {
	    mPrimInStk = true;
        }

	// mMC && mPrimInStk;

        // Event selection
        mTrig = 0;
	mTrackTrig = 0;
        vector<DmpEventSelector::Select> selectTypes = mEventSelector->selectTypes();
        for (int i = 0; i<selectTypes.size();  i++) {
            if (mEventSelector->pass(pev, selectTypes[i])) mTrig++;
            else break;
        }

	// mInterestingEvent = (mTrig == 1) && mPrimInStk;

	// Track selection
	DmpStkTrack * track;
	DmpStkTrack * best_track;
	if (mTrig >= 2 && pev->NStkKalmanTrack() >= 1) { // If skim and het triggers passed
	    int n = 0;
	    TClonesArray * tracks = pev->GetStkKalmanTrackCollection();

	    // Choose the best track
	    for(int i=0; i < pev->NStkKalmanTrack(); i++)
	    {
		track = (DmpStkTrack*)tracks->ConstructedAt(i);
		if (n == 0) {
		    best_track = track;
		}
		else
		{
		    if(mTrackSelector->first_is_better(track, best_track, pev)) {
			best_track = track;
		    }
		}
		n++;
	    }
	    vector<DmpTrackSelector::Select> selectTrTypes = mTrackSelector->selectTypes();
	    for (int i = 0; i<selectTrTypes.size();  i++) {
		if (mTrackSelector->pass(best_track, pev, selectTrTypes[i])) mTrackTrig++;
		else break;
	    }
	}
        // if (mTrackTrig == 6) mInterestingEvent = true;

        // Get energy
        DmpEvtBgoRec * bgoRec = pev->pEvtBgoRec();
        mBGOenergy = bgoRec->GetTotalEnergy();
        if (mMC) {
            DmpEvtSimuPrimaries * prim = pev->pEvtSimuPrimaries();
            mEprim = prim->pvpart_ekin;
        }
        else {
            mEprim = 0.;
        }

        // Event selection
	bool selected = mEventSelector->selected(pev);
	if (!selected) {
	    return;
	}

        // Fill some histograms
	if (mEventSelector->Ntracks() > 0) {
            // STK track
	    DmpStkTrack * stktrack = mEventSelector->stkTrack();

	    mNonoverlaps = stktrack->getNhitX() + stktrack->getNhitY() - 2 * stktrack->getNhitXY();
	    mNonOverlaps->Fill(mNonoverlaps);

	    // Build histogram of S_mean for X and Y
	    float s_mean_xy[2];
	    mTrackSelector->getSMean(stktrack, pev, s_mean_xy);
	    mSmeanX = s_mean_xy[0];
	    mSmeanY = s_mean_xy[1];
	    mSMeanX -> Fill(s_mean_xy[0]);
	    mSMeanY -> Fill(s_mean_xy[1]);

            // Ang dist to BGO
	    mSTK2BGOang -> Fill(mTrackSelector->stk2bgoAngl(stktrack, pev));

	    // Dist to BGO
	    mSTK2BGOdist -> Fill(mTrackSelector->stk2bgoDist(stktrack, pev));

            // chi2
	    mChisqndof = stktrack -> getChi2NDOF();
	    mChi2->Fill(mChisqndof);

            // N points
	    mNpts = stktrack->GetNPoints();
	    mNpoints->Fill(stktrack->GetNPoints());
	}
    }

private :
    bool mPrimInStk;

    DmpEventSelector * mEventSelector;
    DmpTrackSelector * mTrackSelector;

    TH1F * mNonOverlaps;
    TH1F * mSTK2BGOang;
    TH1F * mSTK2BGOdist;
    TH1F * mChi2;
    TH1F * mNpoints;
    TH1F * mSMeanX;
    TH1F * mSMeanY;

    int mNonoverlaps;
    float mSmeanX;
    float mSmeanY;
    float mAngle;
    float mChisqndof;
    int mNpts;

    int mTrig;
    int mTrackTrig;
    double mBGOenergy;
    double mEprim;
};

#endif
