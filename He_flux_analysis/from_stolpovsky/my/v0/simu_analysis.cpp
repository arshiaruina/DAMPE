#import "analysis.hpp"
#import "definitions.hpp"

using namespace myDampeLib;
using namespace std;

class SimuAnalysis : public DmpAnalysis {
public:
    SimuAnalysis(const char * filename, const char * option) :
        DmpAnalysis(filename, option) {;}

    ~SimuAnalysis() {;}

    void addTTree() {;}

    void openTTree(const char * treename)
    {
        cout << "Ok1" << endl;
        DmpAnalysis::openTTree(treename);

        mSlopeXBranch = mTree->Branch("SlopeX", &mSlopeX, "SlopeX");
	mSlopeYBranch = mTree->Branch("SlopeY", &mSlopeY, "SlopeY");
	cout << "Ok2" << endl;
    }

    void analyseOneEvent() {
        mTree->GetEntry(mCurrentEvent);
        DmpEvent * pev = mChain->GetDmpEvent();
	DmpEvtBgoRec * bgoRec = pev->pEvtBgoRec();
	mSlopeX = bgoRec -> GetSlopeXZ ();
	mSlopeY = bgoRec -> GetSlopeYZ ();
    }

    void run(int n=-1)
    { 
        int nnn = (n < 0)? mNEvents : n;
	for(mCurrentEvent = 0; mCurrentEvent < nnn; mCurrentEvent++) {
            if ((mCurrentEvent % (nnn / 10 + 1)) == int(nnn / 10)) {
                int percentage = 10. * int((mCurrentEvent + 1) / int(nnn / 10));
                cout << "Processing percentage: " << percentage << "% \n";
            }

            this -> analyseOneEvent();
            mSlopeXBranch -> Fill();
	    mSlopeYBranch -> Fill();
	}
    }

private:
    float mSlopeX;
    float mSlopeY;

    TBranch * mSlopeXBranch;
    TBranch * mSlopeYBranch;
};

