#include "corrFac.h"
#include "langauConv.h"

//using namespace std; //TODO: remove!


std::pair <double, double> corrFac::MPVGausFit(std::string inFileNameMPVGausFit) {

    //std::string inFileNameMPVGausFit = "../out/20181001_20181009/langaufit_231019_112429.root";
    //std::string inFileNameMPVGausFit = "../out/20181001_20181009/langaufit_231019_112429_FMVAs.root";
    //TFile *inFileMPVGausFit = new TFile(inFileName.c_str());
    TFile *inFileMPVGausFit = new TFile(inFileNameMPVGausFit.c_str());
    //TFile *outFileMPVGausFit = new TFile("../out/20181001_20181009/MPV_langaufit_231019_112429.root", "RECREATE"); 
    //TFile *outFileMPVGausFit = new TFile("../out/20181001_20181009/gausfitMPV_20181001_20181009_FMVAs.root", "RECREATE"); 
    //TFile *outFileMPVGausFit = new TFile("test_gausfitMPV_20181001_20181009_FMVAs.root", "RECREATE"); 
   
    //std::string outFileNameMPVGausFit = "gausfitMPV_20181001_20181009_FMVAs.root"; // this line shifted to the header file on 04.02.2020
    TFile *outFileMPVGausFit = new TFile(outFileNameMPVGausFit.c_str(), "RECREATE"); 

    std::string histName0 = "hMPV0";
    std::string histName1 = "hMPV1";
    TH1D *hist0 = (TH1D*)inFileMPVGausFit->Get(histName0.c_str());
    TH1D *hist1 = (TH1D*)inFileMPVGausFit->Get(histName1.c_str());

    TF1 *fGaus0 = new TF1("fGaus0","gaus");
    fGaus0->SetParameters(1,0,1);
    hist0->Fit(fGaus0,"0");
    fGaus0->SetRange(40,60);

    TF1 *fGaus1 = new TF1("fGaus1","gaus");
    fGaus1->SetParameters(1,0,1);
    hist1->Fit(fGaus1,"0");
    fGaus1->SetRange(30,50);

    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(1111);
    gStyle->SetLabelSize(0.03,"x");
    gStyle->SetLabelSize(0.03,"y");

    TCanvas *c0 = new TCanvas("c0","c0");
    hist0->Draw("hist");
    fGaus0->SetLineColor(kBlack);
    fGaus0->Draw("hist sames");
    gPad->Update();
    TPaveStats* sb0=(TPaveStats*)hist0->FindObject("stats");
    sb0->SetX1NDC(.65);
    sb0->SetX2NDC(.85);
    sb0->SetY1NDC(.65);
    sb0->SetY2NDC(.85);
    sb0->SetTextColor(kBlack);

    TCanvas *c1 = new TCanvas("c1","c1");
    hist1->Draw("hist");
    fGaus1->SetLineColor(kBlack);
    fGaus1->Draw("hist sames");
    gPad->Update();
    TPaveStats* sb1=(TPaveStats*)hist1->FindObject("stats");
    sb1->SetX1NDC(.65);
    sb1->SetX2NDC(.85);
    sb1->SetY1NDC(.65);
    sb1->SetY2NDC(.85);
    sb1->SetTextColor(kBlack);

    c0->Write();
    c1->Write();

    //I made some changes here on 04.02.2020
    //outFileMPVGausFit->Close();
    //outFileMPVGausFit->WriteTObject(hist0);
    //outFileMPVGausFit->WriteTObject(hist1);
    //outFileMPVGausFit->cd();
    outFileMPVGausFit->Write();
    outFileMPVGausFit->Close();

    delete c0;
    delete c1;

    double eqParam0 = fGaus0->GetParameter(1);
    double eqParam1 = fGaus1->GetParameter(1);

    std::cout << "eqParam0: " << eqParam0 << " eqParam1: " << eqParam1 << std::endl;

    return std::make_pair(eqParam0, eqParam1);

} 

std::string corrFac::VAEnergyLangauFit(){

    langauConv langauObj;
    
    histFile.open(histFileName.c_str());
    
    TFile *inFileLangauFit = new TFile(inFileNameLangauFit.c_str());
    TFile *outFileLangauFit = new TFile(outFileNameLangauFit.c_str(), "RECREATE"); 
    

    ////////    Start loop over all histograms      ////////


    while(std::getline(histFile,histName) && countVA < 12){ // for debug run
    //while(std::getline(histFile,histName)){ // for analysis run

        std::cout << outFileNameLangauFit << std::endl;

        std::size_t delim1 = histName.find("_");
        std::size_t delim2 = histName.find("_",delim1+1);
        iladder = std::stoi(histName.substr(delim1+1,delim2-delim1-1));
        iva = std::stoi(histName.substr(delim2+1));
        std::cout << "iladder ... " << iladder << std::endl;
        std::cout << "iva ... " << iva << std::endl;

        //iva = countVA;
        //iladder = 48;

        std::string histName0 = histName + "_0";
        std::string histName1 = histName + "_1";
        TH1D *hist0 = (TH1D*)inFileLangauFit->Get(histName0.c_str()); 
        TH1D *hist1 = (TH1D*)inFileLangauFit->Get(histName1.c_str()); 
        std::cout << histName << std::endl;
        if(histName == "hEtaX" || histName == "hEtaY" || hist0->GetEntries() == 0 || hist1->GetEntries() == 0) continue;
        //if(!(histName == "hVAEnergyY_141_5" || histName == "hVAEnergyY_106_2")) continue;
        //if(!(histName == "hVAEnergyX_154_2" || histName == "hVAEnergyY_33_5")) continue;
        //if(!(histName == "hVAEnergyX_154_2")) continue;
        //if(!(histName == "hVAEnergyX_154_2")) continue;

        //vecHis.push_back(Hist());
        //vecHis.back().hist0  = hist0;
        //vecHis.back().hist1  = hist1;
        //vecHis.back().name   = histName;
        //vecHis.back().ladder = iladder;
        //vecHis.back().va     = iva;

        c1[iladder][iva] = new TCanvas(histName.c_str(),histName.c_str(),800,600);
        if(countVA==0){ // open pdf
            c1[iladder][iva]->Print("plots.pdf[","pdf");
        }

        /**************************************************************
        
            Langau fits on the energy distributions of all VAs
            Generates 1152 plots with 2 distributions on each
            for the two eta regions
        
        **************************************************************/
    
        // --- Hist0 --- //             
                 
        std::cout << "hist0 std dev " << hist0->GetStdDev() << std::endl;
        std::cout << "hist0 mean " << hist0->GetMean() << std::endl;
        std::cout << "hist0 integral " << hist0->Integral() << std::endl;

        // Setting fit range and start values
        double fr0[2], sv0[4], pllo0[4], plhi0[4];
        fr0[0] = 40.;      fr0[1] = 150.;
        sv0[0]=1.8;   sv0[1]=50.0;   sv0[2]=50000.0;       sv0[3]=3.0;
        pllo0[0]=0.5; pllo0[1]=10.0; pllo0[2]=1.0;         pllo0[3]=0.;
        plhi0[0]=9.0; plhi0[1]=70.0; plhi0[2]=1000000.0;   plhi0[3]=15.;
              
        // Return values
        double fp0[4], fpe0[4];
        double chisqr0;
        int ndf0;
        TF1 *fitVAEnergy0 = langauObj.langaufit(hist0,fr0,sv0,pllo0,plhi0,fp0,fpe0,&chisqr0,&ndf0);
        fitVAEnergy0->SetRange(0,200);

        double SNRPeak0, SNRFWHM0;
        langauObj.langaupro(fp0,SNRPeak0,SNRFWHM0);
                
        // --- Hist1 --- //
        
        std::cout << "hist1 std dev " << hist1->GetStdDev() << std::endl;
        std::cout << "hist1 mean " << hist1->GetMean() << std::endl;
        std::cout << "hist1 integral " << hist1->Integral() << std::endl;
        
        // Setting fit range and start values
        double fr1[2], sv1[4], pllo1[4], plhi1[4];
        fr1[0] = 20.;   fr1[1] = 150.;
        sv1[0] = 2.;    sv1[1] = 35.;        sv1[2] = 1e5;       sv1[3] = 1.;
        pllo1[0] = 0.;  pllo1[1] = 10.;      pllo1[2] = 1.0;     pllo1[3] = 0.;
        plhi1[0] = 10.; plhi1[1] = 1e5;      plhi1[2] = 1e6;     plhi1[3] = 10.;
               
        // Return values
        double fp1[4], fpe1[4];
        double chisqr1;
        int ndf1;
        TF1 *fitVAEnergy1 = langauObj.langaufit(hist1,fr1,sv1,pllo1,plhi1,fp1,fpe1,&chisqr1,&ndf1);
        fitVAEnergy1->SetRange(0,200);

        double SNRPeak1, SNRFWHM1;
        langauObj.langaupro(fp1,SNRPeak1,SNRFWHM1);
     
        // Global style settings
        gStyle->SetOptStat(1111);
        gStyle->SetOptFit(111);
        gStyle->SetLabelSize(0.03,"x");
        gStyle->SetLabelSize(0.03,"y");
               
        hist0->SetTitle(histName.c_str());
 
        hist0->SetMaximum((hist0->GetMaximum())*1.1);
        hist0->Draw("hist");
        hist0->SetLineColor(kRed);
        fitVAEnergy0->SetLineColor(kRed);
        fitVAEnergy0->Draw("hist same"); 
        
        hist1->Draw("hist sames");
        hist1->SetLineColor(kBlue);
        fitVAEnergy1->SetLineColor(kBlue);
        fitVAEnergy1->Draw("hist sames");
        
        gPad->Update(); 
        TPaveStats* sb0=(TPaveStats*)hist0->FindObject("stats");
        sb0->SetX1NDC(.65);
        sb0->SetX2NDC(.85);
        sb0->SetY1NDC(.65);
        sb0->SetY2NDC(.85);
        sb0->SetTextColor(kRed);
        TPaveStats* sb1=(TPaveStats*)hist1->FindObject("stats");
        sb1->SetX1NDC(.65);
        sb1->SetX2NDC(.85);
        sb1->SetY1NDC(.4);
        sb1->SetY2NDC(.6);
        sb1->SetTextColor(kBlue);

        /************************************************
        
                   Saving the histograms and fits

        ************************************************/
        
        c1[iladder][iva]->Print("plots.pdf","pdf");
        outFileLangauFit->WriteTObject(c1[iladder][iva]);

        /***********************************************************************
        
            Mean of the fits stored in MPV0[1152] and MPV1[1152]
            to be used later to compute the correction factors for each VA
        
        ***********************************************************************/
        
        MPV0[iladder][iva] = fp0[1];
        MPV1[iladder][iva] = fp1[1];
        std::cout << "MPV0 " << iladder << " " << iva << " " << MPV0[iladder][iva] << std::endl; 

        /***********************************************************************
    
            Mean of the fits (for FM VAs only) stored in hMPV0 and hMPV1
            to make a gaussian fit to them and use the mean as "eq. param."
      
        ***********************************************************************/
        
        std::vector<std::string>::iterator it;
        it = std::find (EQMladders.begin(), EQMladders.end(), histName.substr(11,2));    
        //if(it != EQMladders.end() && histName.substr(13,1) == "_") //EQM ladder found
        if(it == EQMladders.end() || histName.substr(13,1) != "_") { //FM ladder found
        //    continue;
        //else {
            hMPV0->Fill(fp0[1]); // MPVs for eta region 0
            hMPV1->Fill(fp1[1]); // MPVs for eta region 1 
        }    

        std::cout << "-------  Finished working on VA " << countVA << "  ---------" << std::endl; 

        if(countVA==(nVA-1)){ //close pdf
            c1[iladder][iva]->Print("plots.pdf]","pdf");
        }

        countVA++;
        std::cout << "[INFO] Output file: " << outFileNameLangauFit << std::endl;

    } 


    ////////    End of loop over all histograms     ////////


    /******************************************************
    
        Saving the hMPV0 and hMPV1 histograms and fits 
        to the same output file, making a gaussian fit
        and printing the values of the two eq. params.
        (which are the means of the two gaussian fits)
        for the two eta regions)
 
     ******************************************************/
        
    outFileLangauFit->WriteTObject(hMPV0);
    outFileLangauFit->WriteTObject(hMPV1);

    outFileLangauFit->cd();
    outFileLangauFit->Write();
    outFileLangauFit->Close();

    return(outFileNameLangauFit);
    
}

void corrFac::Compute(bool computeFlag){
    
    eqParams = MPVGausFit(VAEnergyLangauFit());   
    
    std::stringstream oss;
    std::ofstream corrFacCompCheck;
    corrFacCompCheck.open("/beegfs/users/ruina/VAequalisation/tmp/corrFacCompCheck.txt");
    oss << "Eq. param. (0): " << eqParams.first << std::endl;
    oss << "Ladder ID  |\t VA ID  |\t MPV(0)  |\t Corr fac \n";
    
    if (!computeFlag) return;   
    
    else{

        TFile *outFileCorrFac = new TFile(outFileNameCorrFac.c_str(), "RECREATE");
    
        //for(int iiladder = 0; iiladder < nLADDER; iiladder++){
        for(int iiladder = 48; iiladder < 50; iiladder++){
        //int iiladder = 48;
            for(int iiva = 0; iiva < nVA; iiva++){
    
                //corrFac for eta region 0, stored in a 2D array
                corrFac0[iiladder][iiva] = eqParams.first/MPV0[iiladder][iiva];
                oss << iiladder << "\t" << iiva << "\t" << MPV0[iiladder][iiva] << "\t" << corrFac0[iiladder][iiva] << "\n";
                
                //corrFac for eta region 1, stored in a 2D array
                corrFac1[iiladder][iiva] = eqParams.second/MPV1[iiladder][iiva];
    
                //corrFac for eta region 0, stored in a histogram, saved to outFileCorrFac
                hCorrFac->SetBinContent(iiladder+1,iiva+1,corrFac0[iiladder][iiva]);
    
                //histograms to see the distribution / spread of the values of the corrFacs
                hCorrFac0->Fill(corrFac0[iiladder][iiva]);
                hCorrFac1->Fill(corrFac1[iiladder][iiva]);
                hCorrFacDiff->Fill(corrFac0[iiladder][iiva]-corrFac1[iiladder][iiva]);
            }
        }
   
        corrFacCompCheck << oss.str();
        corrFacCompCheck.close();
        std::cout << "[INFO] Created tmp/corrFacCompCheck.txt" << std::endl; 
    
        //gStyle->SetPalette(kBird); //can be used only in ROOT 6 and higher versions
        hCorrFac->SetOption("COLZ");
        hCorrFac->SetContour(20);
        hCorrFac->SetStats(0);
        hCorrFac->SetMinimum(0.9);
        hCorrFac->SetMaximum(1.1);
    
        outFileCorrFac->WriteTObject(hCorrFac);
        outFileCorrFac->WriteTObject(hCorrFac0);
        outFileCorrFac->WriteTObject(hCorrFac1);
        outFileCorrFac->WriteTObject(hCorrFacDiff);
        outFileCorrFac->cd();
        outFileCorrFac->Write();
        outFileCorrFac->Close();
    
        return;
    }
}

int main(){

    std::string startPeriodCompute  = "20181001";
    std::string stopPeriodCompute   = "20181009";
    std::string startPeriodApply    = /*"20181201"; "20180901";*/ "20181001"; /* = "20181101"; */
    std::string stopPeriodApply     = /*"20181209"; "20180909";*/ "20181009"; /* = "20181109"; */
    std::string dir, mergeTag;
    std::string task = "CompNotCorr";
    //std::string task = "CompCorr";
    //std::string task = "AppNotCorr";
    //std::string task = "AppCorr";
   
    bool flagComputeCorrFac;
    
    if(task=="CompNotCorr"){
        flagComputeCorrFac = true;
        dir = "";
        //mergeTag = "030220_174755";
        mergeTag = "250220_191643";
    }
    //if(task=="CompCorr"){
    //    dir = "/beegfs/users/ruina/VAequalisation/out/periodCompute/" + startPeriodCompute + "_" + stopPeriodCompute + "/corrected";
    //    mergeTag = "";
    //}
    if(task=="AppNotCorr"){ 
        flagComputeCorrFac = false;
        //201809
        //dir = "/beegfs/users/ruina/VAequalisation/out/periodApply/" + startPeriodApply + "_" + stopPeriodApply + "/not_corrected";
        //mergeTag = "120220_212503"; /*"030220_112552";*/
        //201812
        dir = "/beegfs/users/ruina/VAequalisation/out/periodApply/" + startPeriodApply + "_" + stopPeriodApply + "/not_corrected";
        mergeTag = "120220_213649"; /*"030220_112552";*/
    }
    if(task=="AppCorr"){ 
        flagComputeCorrFac = false;
        //201810
        dir = "/beegfs/users/ruina/VAequalisation/out/periodApply/" + startPeriodApply + "_" + stopPeriodApply + "/corrected";
        mergeTag = "270220_000215";
        //mergeTag = "120220_220045";
        //201809
        //dir = "/beegfs/users/ruina/VAequalisation/out/periodApply/" + startPeriodApply + "_" + stopPeriodApply + "/corrected";
        //mergeTag = "130220_164432";
        //201812
        //dir = "/beegfs/users/ruina/VAequalisation/out/periodApply/" + startPeriodApply + "_" + stopPeriodApply + "/corrected";
        //mergeTag = "130220_165835";
    }
    
    corrFac myCorrFacObj;
    myCorrFacObj.inFileNameLangauFit    = mergeTag + ".root";
    myCorrFacObj.outFileNameLangauFit   = mergeTag + "_out.root";
    myCorrFacObj.outFileNameMPVGausFit  = mergeTag + "_mpv.root";
    myCorrFacObj.outFileNameCorrFac     = mergeTag + "_corr.root"; 

    //TODO make langau fits to histograms where correction factors have been applied
    //TODO give the filename as an argument

    myCorrFacObj.Compute(flagComputeCorrFac);
    
    //std::cout << "Gaussian fit to the MPVs..." << std::endl;
    //std::pair<double, double> eqParams = myCorrFacObj.MPVGausFit(myCorrFacObj.VAEnergyLangauFit());
    //std::cout << eqParams.first << std::endl;
    //std::cout << eqParams.second << std::endl;

    /******************************************************
    
        Computing the correction factors for the two
        eta regions and the difference between the values,
        saving them in histograms.

    ******************************************************/

    return 0;
}

