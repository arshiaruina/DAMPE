#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"
#include "TGaxis.h"
#include "TRandom.h"
 
void checkoverbgoenergybins()
{
 
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

    TH1F *ha1 = new TH1F("ha1","STK cluster energy for BGO energy range 1-2TeV",500,0,500);
    TH1F *ha2 = new TH1F("ha2","STK cluster energy for BGO energy range 1-2TeV",500,0,500);
    TH1F *hb1 = new TH1F("hb1","STK cluster energy for BGO energy range 2-4TeV",500,0,500);
    TH1F *hb2 = new TH1F("hb2","STK cluster energy for BGO energy range 2-4TeV",500,0,500);
    TH1F *hc1 = new TH1F("hc1","STK cluster energy for BGO energy range 4-8TeV",500,0,500);
    TH1F *hc2 = new TH1F("hc2","STK cluster energy for BGO energy range 4-8TeV",500,0,500);
    TH1F *hd1 = new TH1F("hd1","STK cluster energy for BGO energy range >8TeV",500,0,500);
    TH1F *hd2 = new TH1F("hd2","STK cluster energy for BGO energy range >8TeV",500,0,500);

    int n1 = t1->GetEntries();
    for (int i=0; i<n1; i++) {
        t1->GetEntry(i);
        if((ClustETot1>0)&&(ClustETot1<500)&&(ClustEta1<1)){
            if((BGOenergy1>1e6)&&(BGOenergy1<2e6)) ha1->Fill(ClustETot1);       
            if((BGOenergy1>2e6)&&(BGOenergy1<4e6)) hb1->Fill(ClustETot1);       
            if((BGOenergy1>4e6)&&(BGOenergy1<8e6)) hc1->Fill(ClustETot1);       
            if((BGOenergy1>8e6)) hd1->Fill(ClustETot1);       
        }
    }
    int n2 = t2->GetEntries();
    for (int i=0; i<n2; i++) {
        t2->GetEntry(i);
        if((ClustETot2>0)&&(ClustETot2<500)&&(ClustEta2<1)){
            ClustETot2 = ClustETot2 * 60./64.27; 
            //ClustETot2 = ClustETot2/2.; 
            if((BGOenergy2>1e6)&&(BGOenergy2<2e6)) ha2->Fill(ClustETot2);       
            if((BGOenergy2>2e6)&&(BGOenergy2<4e6)) hb2->Fill(ClustETot2);       
            if((BGOenergy2>4e6)&&(BGOenergy2<8e6)) hc2->Fill(ClustETot2);       
            if((BGOenergy2>8e6)) hd2->Fill(ClustETot2);       
        }
    }

    TCanvas *c1 = new TCanvas("c1","c1",800,600);
   
    gStyle->SetOptStat(kFALSE);

    ha1->SetMaximum(ha1->GetMaximum()*1.1);

    ha1->SetLineColor(kBlack);
    ha1->SetLineWidth(0.5);
    ha2->SetMarkerStyle(8);
    ha2->SetMarkerSize(0.4);
    ha2->SetMarkerColor(kBlack);

    hb1->SetLineColor(kBlue);
    hb1->SetLineWidth(0.5);
    hb2->SetMarkerStyle(8);
    hb2->SetMarkerSize(0.4);
    hb2->SetMarkerColor(kBlue);

    hc1->SetLineColor(kRed);
    hc1->SetLineWidth(0.5);
    hc2->SetMarkerStyle(8);
    hc2->SetMarkerSize(0.4);
    hc2->SetMarkerColor(kRed);

    hd1->SetLineColor(kGreen);
    hd1->SetLineWidth(0.5);
    hd2->SetMarkerStyle(8);
    hd2->SetMarkerSize(0.4);
    hd2->SetMarkerColor(kGreen);
    
    ha1->Draw();
    ha2->Draw("psame");
    hb1->Draw("same");
    hb2->Draw("psame");
    hc1->Draw("same");
    hc2->Draw("psame");
    hd1->Draw("same");
    hd2->Draw("psame");

    c1->Update();
    c1->SaveAs("cbgo2.pdf");
 
}
