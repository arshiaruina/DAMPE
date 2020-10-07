//C++
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <sstream>
#include <cmath>
#include <vector>
// ROOT
#include "TROOT.h"
#include "TFile.h"
#include "TMath.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TRandom.h"
#include "TStyle.h"
#include "TPaveStats.h"
//global constants
//#define N_VA        1152
//#define N_LADDER    192
//double corrFac0[N_LADDER][N_VA];
//double corrFac1[N_LADDER][N_VA];
TH2D *hCorrFac = new TH2D("hCorrFac", "Correction factors", 192, 0, 191, 6, 0, 5);
