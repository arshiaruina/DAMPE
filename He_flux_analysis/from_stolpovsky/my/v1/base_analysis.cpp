#ifndef _BASE
#define _BASE

#include "analysis.hpp"
#include "definitions.hpp"
#include "etacorr.hpp"
#include "psd_charge.hpp"
#include "event_selection.hpp"
#include "track_selection.hpp"
#include <sys/stat.h>
#include <unistd.h>

using namespace std;
using namespace myDampeLib;

/*
 * Specify the selection criteria here
 */
class BaseAnalysis : public DmpAnalysis {
public:
    BaseAnalysis(const char * filename,
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

        mTrackSelector->setNimpactPoints(1);
        mTrackSelector->addSelect(DmpTrackSelector::stk_missing_impact_point);
//	mTrackSelector->addSelect(DmpTrackSelector::stk_impact_first_layer);
        mTrackSelector->setMaxAngStk2Bgo(25. / 180. * M_PI); // max angular distance between STK track and BGO shower
        mTrackSelector->setMaxDistStk2Bgo(30.); // max distance between STK track and BGO shower
        mTrackSelector->addSelect(DmpTrackSelector::stk_bgo_dist_high);
        mTrackSelector->setSTKtrackChi2Max(25.);
        mTrackSelector->addSelect(DmpTrackSelector::stk_chi2_cut);
        mTrackSelector->addSelect(DmpTrackSelector::psd_match);
	mTrackSelector->setDistFromPrimCut(1.); // 10 mm limit (cut is set on logarithmic scale in mm)
	mTrackSelector->addSelect(DmpTrackSelector::match_primary);
	mTrackSelector->setComparison(DmpTrackSelector::prim_dist);

        mEventSelector->addSelect(DmpEventSelector::bgo_skim_cuts);
        mEventSelector->addSelect(DmpEventSelector::het);
//        mEventSelector->setBgoPsdDist(40.);
//        mEventSelector->addSelect(DmpEventSelector::bgo_psd_match);
        mEventSelector->setTrackSelector(mTrackSelector);
        mEventSelector->addSelect(DmpEventSelector::has_STK_track);
    }

    ~BaseAnalysis() {
	delete mEventSelector; // track selector is deleted inside event selector
    }

protected:
    DmpEventSelector * mEventSelector;
    DmpTrackSelector * mTrackSelector;
};

#endif
