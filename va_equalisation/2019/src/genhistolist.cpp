/* Code to generate list of VA histogram names "histolist.text" 
 * Author: Arshia Ruina
 * */

// example file
// /beegfs/users/ruina/out/20181019/merged/merged_160919_104734.root

#include "../inc/va_equalisation.h"

using namespace std;

int main(int argc, char** argv) {

    stringstream name;
    ofstream of;
    //of.open("histolist.txt");
    of.open("histolist_new.txt");
    int xLadder, yLadder;

    /* Storing histo names in vectors */
    int iName = 0;

    for(int iladder = 0; iladder < N_LADDER/2; iladder++) {
        if(iladder < 48) {
            xLadder = iladder+48; // X ladders 48-95 
            yLadder = iladder;    // Y ladders 0-47
        }
        else {
            xLadder = iladder+96; // X ladders 144-191
            yLadder = iladder+48; // Y ladders 96-143
        }
        for(int iva = 0; iva < N_VA; iva++){
            //for(int ietareg = 0; ietareg < 2; ietareg++){
                
                histoNamesX.push_back("hVAEnergyX_" + to_string(xLadder) + "_" + to_string(iva));
                histoNamesY.push_back("hVAEnergyY_" + to_string(yLadder) + "_" + to_string(iva));

                //histoNamesX.push_back("hVAEnergyX_" + to_string(xLadder) + "_" + to_string(iva) + "_" + to_string(ietareg));
                //histoNamesY.push_back("hVAEnergyY_" + to_string(yLadder) + "_" + to_string(iva) + "_" + to_string(ietareg));

                //name << "hVAEnergyX_" << xLadder << "_" << iva << "_" << ietareg << endl;
                //name << "hVAEnergyY_" << yLadder << "_" << iva << "_" << ietareg << endl;
                iName++;
            //}
        }
    }

    /* Printing histo names in text file to check */

    //if(histoNamesX.empty() || histoNamesY.empty()) {
    //    cout << "histo name vectors empty!" << endl;
    //    exit(1);
    //}

    //of << "histoNamesX" << "\t" << "histoNamesY" << endl;
    //of << "No. of histoNamesX: " << histoNamesX.size() << endl;
    //of << "No. of histoNamesY: " << histoNamesY.size() << endl;

    for(unsigned int i = 0; i < histoNamesX.size(); i++){
    //for(int i = 0; i < 1152; i++){
        //cout << histoNamesX[i] << endl;
        of << histoNamesX.at(i) << endl;
        //of << histoNamesY.at(i) << endl;
    }
    for(unsigned int i = 0; i < histoNamesY.size(); i++){
        of << histoNamesY.at(i) << endl;
    }
    
    //of << name.str();
    of.close();
    
    return 0;
}
