#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"
#include "TGaxis.h"
#include "TRandom.h"

#include "/atlas/users/stolpovs/DAMPE/myDampeLib/include/etacorr.hpp"
 
void checkoveretabins()
{

    EtaCorr mEtaCorr;
 
    TFile *f1 = new TFile("Results/bkp_081020/all_100_000.root");
    TFile *f2 = new TFile("Results/corr_all_100_000.root");
    TTree *t1 = (TTree*)f1->Get("T");
    TTree *t2 = (TTree*)f2->Get("T");
    float ClustETot1, ClustETot2, BGOenergy1, BGOenergy2, ClustEta1, ClustEta2;
    t1->SetBranchAddress("ClustETot",&ClustETot1);
    t1->SetBranchAddress("BGOenergy",&BGOenergy1);
    t1->SetBranchAddress("ClustEta",&ClustEta1);
    t2->SetBranchAddress("ClustETotCorr",&ClustETot2);
    t2->SetBranchAddress("BGOenergy",&BGOenergy2);
    t2->SetBranchAddress("ClustEta",&ClustEta2);

    TH1F *ha1 = new TH1F("ha1","STK cluster energy for 0.5 < eta < 0.6",500,0,500);
    TH1F *ha2 = new TH1F("ha2","STK cluster energy for BGO energy range 1-2TeV",500,0,500);
    TH1F *hb1 = new TH1F("hb1","STK cluster energy for 0.6 < eta < 0.7",500,0,500);
    TH1F *hb2 = new TH1F("hb2","STK cluster energy for BGO energy range 2-4TeV",500,0,500);
    TH1F *hc1 = new TH1F("hc1","STK cluster energy for 0.7 < eta < 0.8",500,0,500);
    TH1F *hc2 = new TH1F("hc2","STK cluster energy for BGO energy range 4-8TeV",500,0,500);
    TH1F *hd1 = new TH1F("hd1","STK cluster energy for 0.8 < eta < 0.9",500,0,500);
    TH1F *hd2 = new TH1F("hd2","STK cluster energy for BGO energy range >8TeV",500,0,500);
    TH1F *he1 = new TH1F("he1","aSTK cluster energy for 0.9 < eta < 1",500,0,500);
    TH1F *he2 = new TH1F("he2","aSTK cluster energy for BGO energy range >8TeV",500,0,500);

    int n1 = t1->GetEntries();
    for (int i=0; i<n1; i++) {
        t1->GetEntry(i);
        if((ClustETot1>0)&&(ClustETot1<500)){
            if((ClustEta1>=0.5)&&(ClustEta1<0.6)) ha1->Fill(ClustETot1);       
            if((ClustEta1>=0.6)&&(ClustEta1<0.7)) hb1->Fill(ClustETot1);       
            if((ClustEta1>=0.7)&&(ClustEta1<0.8)) hc1->Fill(ClustETot1);       
            if((ClustEta1>=0.8)&&(ClustEta1<0.9)) hd1->Fill(ClustETot1);       
            if((ClustEta1>=0.9)&&(ClustEta1<1.0)) he1->Fill(ClustETot1);       
        }
    }
    int n2 = t2->GetEntries();
    for (int i=0; i<n2; i++) {
        t2->GetEntry(i);
        if((ClustETot2>0)&&(ClustETot2<500)){
            if((ClustEta2>=0.5)&&(ClustEta2<0.6)) ha2->Fill(ClustETot2);       
            if((ClustEta2>=0.6)&&(ClustEta2<0.7)) hb2->Fill(ClustETot2);       
            if((ClustEta2>=0.7)&&(ClustEta2<0.8)) hc2->Fill(ClustETot2);       
            if((ClustEta2>=0.8)&&(ClustEta2<0.9)) hd2->Fill(ClustETot2);       
            if((ClustEta2>=0.9)&&(ClustEta2<1.0)) he2->Fill(ClustETot2);       
        }
    }

    mEtaCorr.fit_pr(ha_1);

    TCanvas *c1 = new TCanvas("c1","c1",500,1000);
    c1.Divide(1,5);
    gStyle->SetOptStat(kFALSE);

    ha1->SetLineColor(kBlack);
    ha2->SetLineColor(kBlue);
    //ha2->SetMarkerStyle(8);
    //ha2->SetMarkerSize(0.3);
    //ha2->SetMarkerColor(kBlack);

    hb1->SetLineColor(kBlack);
    hb2->SetLineColor(kBlue);
    //hb2->SetMarkerStyle(8);
    //hb2->SetMarkerSize(0.3);
    //hb2->SetMarkerColor(kBlue);

    hc1->SetLineColor(kBlack);
    hc2->SetLineColor(kBlue);
    //hc2->SetMarkerStyle(8);
    //hc2->SetMarkerSize(0.3);
    //hc2->SetMarkerColor(kRed);

    hd1->SetLineColor(kBlack);
    hd2->SetLineColor(kBlue);
    //hd2->SetMarkerStyle(8);
    //hd2->SetMarkerSize(0.3);
    //hd2->SetMarkerColor(kGreen);
    
    he1->SetLineColor(kBlack);
    he2->SetLineColor(kBlue);
    //he2->SetMarkerStyle(8);
    //he2->SetMarkerSize(0.3);
    //he2->SetMarkerColor(6);
    
    c1_1->cd();
    ha1->Draw();
    ha2->Draw("same");
    
    c1_2->cd();
    hb1->Draw();
    hb2->Draw("same");
    
    c1_3->cd();
    hc1->Draw();
    hc2->Draw("same");
    
    c1_4->cd();
    hd1->Draw();
    hd2->Draw("same");
    
    c1_5->cd();
    he1->Draw();
    he2->Draw("same");

    c1->Update();
    c1->SaveAs("overeta.pdf");
 
}
