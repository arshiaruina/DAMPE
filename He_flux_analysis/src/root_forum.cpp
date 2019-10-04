/* Code to compute eta values
 * Author: Arshia Ruina
 * To compile, g++ -Wall `root-config --cflags --libs` -I$DMPSWSYS/include -L$DMPSWSYS/lib -lDmpEvent -lDmpService eta.cc -o ../submit/eta.exe
 * */

// example file
// /beegfs/users/ruina/out/20181019/merged/merged_160919_104734.root

#include "../inc/va_equalisation.h"
#include "mylangaus.C"
#include "TKey.h"

using namespace std;

int main(int argc, char** argv) {

    std::string inFileName = argv[1];
    TFile *inFile = new TFile(inFileName.c_str());
    TList *list = inFile->GetListOfKeys();
    TFile *outFile = new TFile("../out/root_forum.root", "recreate");

    //------------------------------------------------//

    //string classname = "TH1D";
    //TH1D *hist;
    //TIter iter(list); //or TIter iter(list->MakeIterator());
    //while(hist = (TH1D*)iter()){
    //    if(hist->ClassName() == classname.c_str()){
    //        outFile->WriteTObject(hist);
    //    }
    //}


    //------------------------------------------------//
    
    TKey *key;
    TIter iter(list); //or TIter iter(list->MakeIterator());
    static TString classname("TH1D");
    while((key = (TKey*)iter())) {
        if (key->GetClassName() == classname) {
            TH1D *hist = (TH1D*)key->ReadObj();
            if (hist) {
                const char* histName_c = hist->GetName();
                std::string histName_s(histName_c);
                if(histName_s.back() == '0')
                std::cout << hist->GetName() << " " << hist->GetEntries() << std::endl;
                //outFile->WriteTObject(hist);
                delete hist;
            }
        }
    }

    //------------------------------------------------//
    
    outFile->cd();
    outFile->Write();
    outFile->Close();

    return 0;
}
