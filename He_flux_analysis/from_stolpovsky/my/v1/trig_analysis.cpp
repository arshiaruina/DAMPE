#ifndef _TRIG
#define _TRIG

#include <string>
#include <stdio.h>
#include <cmath>
#include <iostream>

#include "analysis.hpp"
#include "base_analysis.cpp"
#include "definitions.hpp"
#include "etacorr.hpp"
#include "psd_charge.hpp"
#include "event_selection.hpp"
#include "track_selection.hpp"

using namespace std;
using namespace myDampeLib;

class TrigAnalysis : public BaseAnalysis {
public :
    TrigAnalysis(const char * filename,
                 const char * option) :
        BaseAnalysis(filename, option) 
    {;}

    ~TrigAnalysis()
    {
	mOutputFile -> cd();
        mEventSelector -> hSelect() -> Write();
	closeOutputFile();
    }

    void addTTree()
    {
        mTree = new TTree("Trig", "Tree");
	mTree->Branch("Trig", &mTrig, "mTrig/I");
        mTree->Branch("BGOenergy", &mBGOenergy, "mBGOenergy/D");
	mTree->Branch("Eprim", &mEprim, "mEprim/D");
        mTree->SetDirectory(0);
    }

    void analyseOneEvent(DmpEvent * pev)
    {
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

        // Fill hSelect
	mEventSelector->selected(pev);

        // Event selection
	mTrig = 0;
	vector<DmpEventSelector::Select> selectTypes = mEventSelector->selectTypes();
        for (int i = 0; i<selectTypes.size();  i++) {
	    if (mEventSelector->pass(pev, selectTypes[i])) mTrig++;
	    else break;
	}

//        if(mTrig == 2) mInterestingEvent = true;
    }

private :
    int mTrig;
    double mBGOenergy;
    double mEprim;
};

#endif
