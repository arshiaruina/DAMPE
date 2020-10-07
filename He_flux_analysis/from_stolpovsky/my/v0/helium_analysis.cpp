#import "analysis.hpp"
#import "definitions.hpp"

using namespace myDampeLib;
using namespace std;

class HeliumAnalysis : public DmpAnalysis {
public:
    HeliumAnalysis(const char * filename, const char * option) :
        DmpAnalysis(filename, option) {;}

    ~HeliumAnalysis() {;}

    void addTTree()
    {
        mTree = new TTree("T", "Tree");

        addBranch(&mBGOenergy, "BGOenergy");
        addBranch(mBGOsigX, NBARSXY, "BGOsigX");
	addBranch(mBGOsigY, NBARSXY, "BGOsigY");
    }

    void analyseOneEvent() {
        DmpEvent * pev = mChain->GetDmpEvent();
	DmpEvtBgoRec * bgoRec = pev->pEvtBgoRec();
        mBGOenergy = bgoRec->GetTotalEnergy();

	for (int ilayer=0; ilayer<NLAYERSXY; ilayer++){
            for (int ibar=0; ibar<NBARSL; ibar++) {
                mBGOsigX[ilayer*22 + ibar] = bgoRec->GetEdep(ilayer * 2 + 1, ibar);
		mBGOsigY[ilayer*22 + ibar] = bgoRec->GetEdep(ilayer * 2, ibar);
            }
	}
    }

private:
    double mBGOenergy;
    float mBGOsigX[NBARSXY];
    float mBGOsigY[NBARSXY];
};

