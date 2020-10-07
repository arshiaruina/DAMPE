/// \file
/// \ingroup tutorial_tmva
/// \notebook -nodraw
/// This macro provides a simple example for the training and testing of the TMVA
/// multiclass classification
/// - Project   : TMVA - a Root-integrated toolkit for multivariate data analysis
/// - Package   : TMVA
/// - Root Macro: TMVAMulticlass
///
/// \macro_output
/// \macro_code
/// \author Andreas Hoecker

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"


#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
//#include "TMVA/DataLoader.h"
//#include "TMVA/TMVAMultiClassGui.h"


using namespace TMVA;

void TMVAMulticlass( TString myMethodList = "", TString energy_tag = "100GeV-10TeV" )
{

   // This loads the library
   TMVA::Tools::Instance();

   // to get access to the GUI and all tmva macros
   //
   //     TString tmva_dir(TString(gRootDir) + "/tmva");
   //     if(gSystem->Getenv("TMVASYS"))
   //        tmva_dir = TString(gSystem->Getenv("TMVASYS"));
   //     gROOT->SetMacroPath(tmva_dir + "/test/:" + gROOT->GetMacroPath() );
   //     gROOT->ProcessLine(".L TMVAMultiClassGui.C");


   //---------------------------------------------------------------
   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;
   Use["MLP"]             = 0;
   Use["BDTG"]            = 1;
   Use["DNN_CPU"]         = 0;
   Use["FDA_GA"]          = 0;
   Use["PDEFoam"]         = 0;
   //---------------------------------------------------------------

   std::cout << std::endl;
   std::cout << "==> Start TMVAMulticlass" << std::endl;

   if (myMethodList != "") {
      for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

      std::vector<TString> mlist = TMVA::gTools().SplitString( myMethodList, ',' );
      for (UInt_t i=0; i<mlist.size(); i++) {
         std::string regMethod(mlist[i]);

         if (Use.find(regMethod) == Use.end()) {
            std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
            for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
            std::cout << std::endl;
            return;
         }
         Use[regMethod] = 1;
      }
   }

   // Create a new root output file.
   TString outfileName = TString("TMVAMulticlass_") + energy_tag + TString(".root");
   std::cout << "Output File: " << outfileName << std::endl;
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

   TMVA::Factory *factory = new TMVA::Factory( "TMVAMulticlass", outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=multiclass" );
//   TMVA::DataLoader *dataloader = new TMVA::DataLoader("dataset");

  // factory->AddVariable( "HET", 'B' );
   factory->AddVariable( "PSDenergy[0]", 'F' );
   factory->AddVariable( "PSDenergy[1]", 'F' );
   factory->AddVariable( "Ntracks", 'F');
   factory->AddVariable( "STKenergy[0]", 'F' );
   factory->AddVariable( "STKenergy[1]", 'F' );
   factory->AddVariable( "STKenergy[2]", 'F' );
   factory->AddVariable( "STKenergy[3]", 'F' );
   factory->AddVariable( "BGOenergy", 'F');

   TFile * input_p_LE  (0);
   TFile * input_p_HE  (0);
   TFile * input_He_LE (0);
   TFile * input_He_HE (0);
   TFile * input_C_LE  (0);
   TFile * input_C_HE  (0);
   TString fname_p_LE   = TString("../../Results/allProton-100GeV-10TeV.root");
   TString fname_p_HE   = TString("../../Results/allProton-10TeV-100TeV.root");
   TString fname_He_LE  = TString("../../Results/allHe-100GeV-10TeV.root");
   TString fname_He_HE  = TString("../../Results/allHe-10TeV-100TeV.root");
   TString fname_C_LE   = TString("../../Results/allC-100GeV-10TeV.root");
   TString fname_C_HE   = TString("../../Results/allC-10TeV-100TeV.root");
   if (!gSystem->AccessPathName( fname_p_LE ) &&
       !gSystem->AccessPathName( fname_He_LE) &&
       !gSystem->AccessPathName( fname_C_LE )) {
      input_p_LE  = TFile::Open( fname_p_LE  ); // check if file in local directory exists
      input_He_LE = TFile::Open( fname_He_LE );
      input_C_LE  = TFile::Open( fname_C_LE  );
   }
   if (!gSystem->AccessPathName( fname_p_HE ) &&
       !gSystem->AccessPathName( fname_He_HE) &&
       !gSystem->AccessPathName( fname_C_HE )) {
       input_p_HE  = TFile::Open( fname_p_HE  ); // check if file in local directory exists
       input_He_HE = TFile::Open( fname_He_HE );
       input_C_HE  = TFile::Open( fname_C_HE  );
   }


//   if (!input_p || !input_He || /*!input_Li ||*/ !input_C ) {
//     std::cout << "ERROR: could not open data file" << std::endl;
//     std::cout << "File names are:" << std::endl;
//     std::cout << fname_p << "\t" << fname_He << "\t" << fname_C << std::endl;
//     exit(1);
//   }

//   std::cout << "--- TMVAClassification       : Using input file: " << input_p->GetName()
//	     << ", " << input_He->GetName()
//	     << ", " << input_C->GetName()
//	     << std::endl;

   // Register the training and test trees

   TTree *signalTree_LE  = (TTree*)input_He_LE ->Get("T");
   TTree *signalTree_HE  = (TTree*)input_He_HE ->Get("T");
   TTree *background0_LE = (TTree*)input_p_LE  ->Get("T");
   TTree *background0_HE = (TTree*)input_p_HE  ->Get("T");
   TTree *background1_LE = (TTree*)input_C_LE  ->Get("T");
   TTree *background1_HE = (TTree*)input_C_HE  ->Get("T");

   gROOT->cd( outfileName+TString(":/") );
   factory->AddSignalTree(signalTree_LE, 1. / signalTree_LE->GetEntries(), "Training and Testing");
   factory->AddSignalTree(signalTree_HE, 1. / signalTree_HE->GetEntries(), "Training and Testing");
//   factory->AddTree(signalTree_LE  , "Signal");
//   factory->AddTree(signalTree_HE  , "Signal");
   factory->AddTree(background0_LE , "bg0");
   factory->AddTree(background0_HE , "bg0");
   factory->AddTree(background1_LE , "bg1");
   factory->AddTree(background0_HE , "bg0");

   TCut cuts = ""; //PSDenergy[0]<100. && PSDenergy[1] < 100. && STKenergy[0]<3000. && STKenergy[1]<3000. && STKenergy[2]<3000. && STKenergy[3]<3000.";
   factory->PrepareTrainingAndTestTree( cuts, //"SplitMode=Random:NormMode=NumEvents:!V");
                                        "nTrain_Signal=20000:nTrain_bg0=20000:nTrain_bg1=20000:nTest_Signal=20000:nTest_bg0=20000:nTest_bg1=20000:SplitMode=Random:NormMode=NumEvents:!V" );

   if (Use["BDTG"]) // gradient boosted decision trees
      factory->BookMethod( TMVA::Types::kBDT, "BDTG", "!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.50:nCuts=20:MaxDepth=2");
   if (Use["MLP"]) // neural network
      factory->BookMethod( TMVA::Types::kMLP, "MLP", "!H:!V:NeuronType=tanh:NCycles=1000:HiddenLayers=N+5,5:TestRate=5:EstimatorType=MSE");
   if (Use["FDA_GA"]) // functional discriminant with GA minimizer
      factory->BookMethod( TMVA::Types::kFDA, "FDA_GA", "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:PopSize=300:Cycles=3:Steps=20:Trim=True:SaveBestGen=1" );
   if (Use["PDEFoam"]) // PDE-Foam approach
      factory->BookMethod( TMVA::Types::kPDEFoam, "PDEFoam", "!H:!V:TailCut=0.001:VolFrac=0.0666:nActiveCells=500:nSampl=2000:nBin=5:Nmin=100:Kernel=None:Compress=T" );

   if (Use["DNN_CPU"]) {
      TString layoutString("Layout=TANH|100,TANH|50,TANH|10,LINEAR");
      TString training0("LearningRate=1e-1, Momentum=0.5, Repetitions=1, ConvergenceSteps=10,"
                        " BatchSize=256, TestRepetitions=10, Multithreading=True");
      TString training1("LearningRate=1e-2, Momentum=0.0, Repetitions=1, ConvergenceSteps=10,"
                        " BatchSize=256, TestRepetitions=7, Multithreading=True");
      TString trainingStrategyString("TrainingStrategy=");
      trainingStrategyString += training0 + "|" + training1;
      TString nnOptions("!H:V:ErrorStrategy=CROSSENTROPY:VarTransform=N:"
                        "WeightInitialization=XAVIERUNIFORM:Architecture=CPU");
      nnOptions.Append(":");
      nnOptions.Append(layoutString);
      nnOptions.Append(":");
      nnOptions.Append(trainingStrategyString);
      factory->BookMethod( TMVA::Types::kDNN, "DNN_CPU", nnOptions);
   }

   // Train MVAs using the set of training events
   factory->TrainAllMethods();

   // Evaluate all MVAs using the set of test events
   factory->TestAllMethods();

   // Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();

   // --------------------------------------------------------------

   // Save the output
   outputFile->Close();

   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVAMulticlass is done!" << std::endl;

   if (!gROOT->IsBatch()) TMVAMultiClassGui( outputFile->GetName() );

   delete factory;

}

int main( int argc, char** argv )
{
   // Select methods (don't look at this code - not of interest)
   TString methodList;
   for (int i=1; i<argc; i++) {
      TString regMethod(argv[i]);
      if(regMethod=="-b" || regMethod=="--batch") continue;
      if (!methodList.IsNull()) methodList += TString(",");
      methodList += regMethod;
   }
   TMVAMulticlass(methodList, TString("100GeV-10TeV"));
   return 0;
}

