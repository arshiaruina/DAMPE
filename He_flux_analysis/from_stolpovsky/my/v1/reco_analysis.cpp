#ifndef _RECO
#define _RECO

#include <string>
#include <stdio.h>
#include <cmath>
#include <algorithm>
#include <iostream>

#include "analysis.hpp"
#include "definitions.hpp"
#include "etacorr.hpp"
#include "psd_charge.hpp"
#include "event_selection.hpp"
#include "track_selection.hpp"
#include "base_analysis.cpp"

#include "TMVA/Reader.h"
#include "TMVA/Tools.h"


using namespace std;
using namespace myDampeLib;
using namespace TMVA;

class RecoAnalysis : public BaseAnalysis {
public :
    RecoAnalysis(const char * filename,
                 const char * option) :
        BaseAnalysis(filename, option) 
    {
        mEtaCorr.setTargetParameters(64.27, 6., 257.75, 6.*257.75/64.27);

        //TMVA
        Tools::Instance();
        mReader = new Reader( "!Color:!Silent" );
        mReader -> AddVariable( "PSDenergy[0]", &mPSDenergy_0 );
        mReader -> AddVariable( "PSDenergy[1]", &mPSDenergy_1 );
        mReader -> AddVariable( "Ntracks", &mNtracks );
        mReader -> AddVariable( "STKenergy[0]", &mSTKenergy_0 );
        mReader -> AddVariable( "STKenergy[1]", &mSTKenergy_1 );
        mReader -> AddVariable( "STKenergy[2]", &mSTKenergy_2 );
        mReader -> AddVariable( "STKenergy[3]", &mSTKenergy_3 );
        mReader -> AddVariable( "BGOenergy", &mBGOenergy );

        mReader->BookMVA( "BDTG", "/atlas/users/stolpovs/DAMPE/helium_flux/ML/weights/TMVAMulticlass_BDTG.weights.xml" );
    }

    ~RecoAnalysis() {
        mOutputFile->cd();
        mEventSelector->hSelect() -> Write();
        closeOutputFile();
    }

    void addTTree()
    {
        mTree = new TTree("T", "Tree");

        mTree->Branch("PSDenergy", &mPSDenergy, "mPSDenergy[2]/F");
        mTree->Branch("STKenergy", &mSTKenergy, "mSTKenergy[12]/F");
        mTree->Branch("Ntracks", &mNtracks, "mNtracks/F");
        mTree->Branch("BGOenergy", &mBGOenergy, "mBGOenergy/F");
        mTree->Branch("BGOsigX", &mBGOsigX, "mBGOsigX[154]/F");
        mTree->Branch("BGOsigY", &mBGOsigY, "mBGOsigY[154]/F");
        mTree->Branch("Eprim", &mEprim, "Eprim/F");
        mTree->Branch("ThetaX", &mThetaX, "ThetaX/F");
        mTree->Branch("ThetaY", &mThetaY, "ThetaY/F");
	mTree->Branch("PrimInStk", &mPrimInStk, "mPrimInStk/B");
	mTree->Branch("BGOshowerProjCoord", &mBGOshowerProjCoord, "BGOshowerProjCoord[2]/F");
	mTree->Branch("PrimDistTop", &mPrimDistTop, "PrimDistTop/F");
	mTree->Branch("PrimDistBot", &mPrimDistBot, "PrimDistBot/F");
	mTree->Branch("STKtrackImpactPointPlane", &mSTKtrackImpactPointPlane, "STKtrackImpactPointPlane/F");
	mTree->Branch("STKchi2", &mSTKchi2, "STKchi2/F");

        mTree->Branch("MVA", &mMVA, "MVA[3]/F");

        mTree->SetDirectory(0);
    }

    int stkZ2L(float z) {
	// Transform z measurement to the layer number
	vector<int> layers({210, 206, 176, 173, 144, 140, 111, 107, 79, 74, 46, 41});
	vector<int>::iterator it = find(layers.begin(), layers.end(), (int)(-z));
	if (it != layers.end()) return distance(layers.begin(), it);
	else return -1;
    }

    void analyseOneEvent(DmpEvent * pev)
    {
        // Event selection
        mSelected = mEventSelector->selected(pev);
        if (!mSelected) return;

        DmpEvtBgoRec * bgoRec = pev->pEvtBgoRec();
        mBGOenergy = (float)(bgoRec->GetTotalEnergy());


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

            float psd_top_xz = -291.5; // z coord of the very first layer of PSD
            float psd_top_yz = -317.7;
            float stk_top_z = -210.0319976349; // z coord of the very first stk layer (farthest from bgo)
            float stk_bot_z = -41.96795840450001; // z coord of the last stk layer
            float bgo_bot_z = 448; // last layer of BGO

            float x_top = x + (psd_top_xz - z) * px/pz;
            float y_top = y + (psd_top_yz - z) * py/pz;
            float x_bot = x + (bgo_bot_z - z) * px/pz;
            float y_bot = y + (bgo_bot_z - z) * py/pz;

            float xy_lim_top = 400.;
            float xy_lim_bot = 280.;

            mPrimInStk = abs(x_top) < xy_lim_top && abs(x_bot) < xy_lim_bot &&
                         abs(y_top) < xy_lim_top && abs(y_bot) < xy_lim_bot;

//            x_top = x + (stk_top_z - z) * px/pz;
//            y_top = y + (stk_top_z - z) * py/pz;
//            x_bot = x + (bgo_bot_z - z) * px/pz;
//            y_bot = y + (bgo_bot_z - z) * py/pz;

//           xy_lim_top = 400.;
//            xy_lim_bot = 280.;

//            primInStk_old = abs(x_top) < xy_lim_top && abs(x_bot) < xy_lim_bot &&
//                            abs(y_top) < xy_lim_top && abs(y_bot) < xy_lim_bot;
        }
        else {
            mPrimInStk = true;
        }

        // Event selection
        mSelected = mEventSelector->selected(pev);
	// Use mInterestingEvent to printout the event (just set it true)
	// mInterestingEvent = mPrimInStk && !mSelected;
        if (!mSelected) return;

        // STK track
        DmpStkTrack * stktrack = mEventSelector->stkTrack();
        mNtracks = (float)(mEventSelector->Ntracks());
        TClonesArray * stkclusters = pev->GetStkSiClusterCollection();
        double costheta = stktrack->getDirection().CosTheta();
        for (int ipoint=0; ipoint<stktrack->GetNPoints(); ipoint++) {
            for (int ixy=0; ixy<2; ixy++) {
                // Check if cluster is associated to a hit
                DmpStkSiCluster* cluster;
                if(ixy == 0)
	            cluster = stktrack -> GetClusterX(ipoint, stkclusters);
	        else
	            cluster = stktrack -> GetClusterY(ipoint, stkclusters);
                if(!cluster) continue;
		// double z = cluster -> GetZ();
                double e = cluster -> getEnergy() * costheta;
                double eta = mEtaCorr.calcEta(cluster);
		int index = ixy + (cluster -> getPlane()) * 2;
		mSTKenergy[index] = mEtaCorr.corrEnergy(e, eta, costheta);
            }
        }
	cout << "Dist from prim = " << mTrackSelector->dist_from_prim(stktrack, pev) << endl;

	// MC only
	if (mMC) {
	    DmpEvtSimuPrimaries * prim = pev->pEvtSimuPrimaries();

            float x = prim->pv_x;
            float y = prim->pv_y;
            float z = prim->pv_z;

            float px = prim->pvpart_px;
            float py = prim->pvpart_py;
            float pz = prim->pvpart_pz;

	    float stk_x = stktrack->getImpactPoint().X();
	    float stk_y = stktrack->getImpactPoint().Y();
	    float stk_z = stktrack->getImpactPoint().Z();

            float x_ = x + (stk_z - z) * px/pz;
            float y_ = y + (stk_z - z) * py/pz;

	    mPrimDistTop = (x_ - stk_x) * (x_ - stk_x) + 
                           (y_ - stk_y) * (y_ - stk_y);
  
	    float slopeX = stktrack->getSlopeX(0);
	    float slopeY = stktrack->getSlopeY(0);

	    float stk_bot_z = -41.96795840450001; // z coord of the last stk layer

	    stk_x += (stk_bot_z - stk_z) * slopeX;
	    stk_y += (stk_bot_z - stk_z) * slopeY;

	    x_ = x + (stk_bot_z - z) * px/pz;
	    y_ = y + (stk_bot_z - z) * py/pz;

	    mPrimDistBot = (x_ - stk_x) * (x_ - stk_x) +
                           (y_ - stk_y) * (y_ - stk_y);
	}
	else {
	    mPrimDistTop = -1.;
	    mPrimDistBot = -1.;
	}
	mSTKtrackImpactPointPlane = float(stktrack->getImpactPointPlane());
	mSTKchi2 = float(stktrack->getChi2NDOF());

        // PSD
        DmpEvtPsdHits *psdHits = pev->pEvtPsdHits();
        int * ibar1 = new int(-1);
        int * ibar2 = new int(-1);
        mPSDenergy[0] = psdEnergy(stktrack, psdHits, ibar1, mMC, 1);
        mPSDenergy[1] = psdEnergy(stktrack, psdHits, ibar2, mMC, 2);
        delete ibar1;
        delete ibar2;

	// BGO to PSD match
	float psdXZ = -291.5;
	float psdYZ = -317.7;
	mBGOshowerProjCoord[0] = bgoRec->GetSlopeXZ()*psdXZ + bgoRec->GetInterceptXZ();
	mBGOshowerProjCoord[1] = bgoRec->GetSlopeYZ()*psdYZ + bgoRec->GetInterceptYZ();


        // Loop over the BGO bars
        for (int ilayer=0; ilayer<NLAYERSXY; ilayer++){
            for (int ibar=0; ibar<NBARSL; ibar++) {
                mBGOsigX[ilayer*22 + ibar] = bgoRec->GetEdep(ilayer * 2 + 1, ibar);
                mBGOsigY[ilayer*22 + ibar] = bgoRec->GetEdep(ilayer * 2, ibar);
            }
        }

        mThetaX = bgoRec -> GetSlopeXZ ();
        mThetaY = bgoRec -> GetSlopeYZ ();

        // TMVA
        bool cut = true; //mPSDenergy[0] < 100. && 
	//                   mPSDenergy[1] < 100. && 
	//                   mSTKenergy[0] < 3000. && 
	//                   mSTKenergy[1] < 3000. &&
	//                   mSTKenergy[2] < 3000. &&
	//                   mSTKenergy[3] < 3000.;

        if (cut) {
            mPSDenergy_0 = mPSDenergy[0];
            mPSDenergy_1 = mPSDenergy[1];
            mSTKenergy_0 = mSTKenergy[0];
            mSTKenergy_1 = mSTKenergy[1];
            mSTKenergy_2 = mSTKenergy[2];
            mSTKenergy_3 = mSTKenergy[3];
            vector<float> mva(mReader->EvaluateMulticlass( "BDTG" ));
            for(int i=0; i<3; i++) mMVA[i] = mva[i];
        }

    }

private :
    bool ** mBadChannels;

    EtaCorr mEtaCorr;

    float mPSDenergy[NPSDL]; // Note that TMVA doesn't work with double, only float
    float mSTKenergy[NSTKL];
    float mNtracks;
    float mBGOenergy;
    float mEprim;
    float mBGOsigX[NBARSXY];
    float mBGOsigY[NBARSXY];
    float mThetaX, mThetaY;
    float mBGOshowerProjCoord[2]; // Coordinates of the BGO shower projection on the level of PSD
    bool mPrimInStk;

    float mSTKchi2;
    float mPrimDistTop;
    float mPrimDistBot;
    float mSTKtrackImpactPointPlane; // TMVA doesn't like ints

    // TMVA
    Reader * mReader;
    float mMVA[3]; 
    float mPSDenergy_0;
    float mPSDenergy_1;
    float mSTKenergy_0;
    float mSTKenergy_1;
    float mSTKenergy_2;
    float mSTKenergy_3;
};

#endif
