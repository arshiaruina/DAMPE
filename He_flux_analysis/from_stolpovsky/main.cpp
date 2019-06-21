#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <string>
#include <stdio.h>
#include <sstream>
#include <math.h>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TChain.h"
#include "DmpEvtHeader.h"
#include "DmpEvtPsdRec.h"
#include "DmpRootEvent.h"
#include "DmpChain.h"
#include "DmpEvtGlobTrack.h"
#include "DmpSvcPsdEposCor.h"
#include "TClonesArray.h"

#include "DmpStkTrack.h"
#include "DmpStkTrackHelper.h"
#include "DmpStkSiCluster.h"
#include "DmpEvtPsdRec.h"
#include "TStopwatch.h"
#include "DmpStkClusterCalibration.h"

#include "track_selection.hpp"
#include "etacorr.hpp"

// See some global variables defined in track_selection.h

using namespace std;


double psdEnergy(DmpStkTrack * track, DmpEvtPsdRec *psdRec, 
                 int * ibar1, int * ibar2, bool mc=false, int layer=0) {
    // layer -- get energy from a specific psd layer

    TVector3 globImpact = track->getImpactPoint();
    TVector3 globDirection = track->getDirection();

    double e = 0.;
    int nhits = 0;

    // Loop over psd layers
    for (int ilayer = 0; ilayer < 2; ilayer++) {
        if (layer != 0 && ilayer != layer-1) continue;
        // Loop over psd bars
        for (int ibar = 0; ibar < 41; ibar++) {
            // psd energy, no correction
            double etemp = psdRec->GetEdep(ilayer, ibar);

            // Get hit point and path length
            double * len = new double[2];
            if(!mc && !gPsdECor->GetPathLengthPosition(ilayer, ibar, globDirection, globImpact, len))
                continue;
            if(mc && !gPsdECor->GetPathLPMC(ilayer, ibar, globDirection, globImpact, len))
                continue;
            // len[0] -- hit point
            // len[1] -- path length

            switch (ilayer) {
                case 0 : *ibar1 = ibar; break;
                case 1 : *ibar2 = ibar; break;
            }

            // correction for the path length
            etemp *= 10. / len[1]; // 10. is the psd bar thickness in mm

            // correction for the attenuation
            double corr = gPsdECor -> GetPsdECorSp3(ilayer, ibar, len[0]);
            etemp *= corr;

            e += etemp;
            nhits += 1;
        } // end bar loop
    } // end layer loop
    double ret = e / nhits;
    if(!std::isnormal(ret)) ret = 0.;
    return ret;
}

double psdCharge(double e, double proton_peak=2.07) {
    // proton_peak -- position of the proton peak in the dE/dx distribution
    return TMath::Sqrt(e / proton_peak);
}

double psdCharge(DmpStkTrack * track, DmpEvtPsdRec *psdRec, 
                 int * ibar1, int * ibar2, bool mc=false) {
    double e = psdEnergy(track, psdRec, ibar1, ibar2, mc);
    return psdCharge(e);
}

TH2D * hist_EvnEta(string t_name, string t_title) {
    TH2D * h = new TH2D(t_name.c_str(), t_title.c_str(), 50, 0., 1., 100, 0., 500.);
    h -> SetXTitle("#eta");
    h -> SetYTitle("STK #frac{dE}{dx}, [ADC]");
    h -> SetOption("colz");
    return h;
}


int main( int argc , char *argv[]){ 
    TStopwatch sw;
    sw.Start();

    bool test = false;

    // TChain
    DmpChain *t = new DmpChain("CollectionTree");
    TFile * f;

    ifstream runlist(argv[1]);
    string line;
    int nfilesChained = 0;

    while(getline(runlist, line)) {
        f = TFile::Open(line.c_str(), "READ");
        if(f) {
            cout << "File found " << line << endl;
            t->Add(line.c_str());
            nfilesChained++;
        }
    }
    runlist.close();
    // TChain done

    // Bad channels list
    bool** badchannels =  read_bad_channels_file(argv[3]);

    // EtaCorr object
    EtaCorr etacorr;

    string datfile_name("etacorr.dat");
    ifstream ifile(datfile_name);
    bool datfile_exists = (bool)ifile;
    ifile.close();

    int nbins;
    TH1D * slices[1000];
    if(!datfile_exists){
        TFile * f_etacorrh = TFile::Open("EvsEta_hist.root");
        TH2D * he_eta_ = nullptr;
        f_etacorrh->GetObject("EvsEta_He", he_eta_);
        etacorr.fit_pr_He(he_eta_);
        nbins = he_eta_->GetNbinsX();
        for (int ibin=1; ibin<=nbins; ibin++) {
            string hname = Form("EvsEta_%d", ibin);
            slices[ibin - 1] = he_eta_->ProjectionY(hname.c_str(), ibin, ibin);
            slices[ibin - 1]->GetListOfFunctions()->Add(etacorr.getFitFunction(ibin - 1));
            slices[ibin - 1]->SetDirectory(0);
        }
        f_etacorrh->Close();
        etacorr.writeCorrFile(datfile_name);
    }
    else {
        etacorr.readCorrFile(datfile_name);
    }
    // etacorr.setTargetParameters(64.27, 257.75); // from Vitillo's code
    etacorr.setTargetParameters(64.27, 6., 257.75, 6.*257.75/64.27);

    // Output file
    TString output_file_name = argv[2];
    TFile *newrootfile= new TFile(output_file_name, "RECREATE");

    // Loop over events
    int NEvents = test ? 1000 : t->GetEntries();
    cout << "Nentries "<< NEvents << endl;

    // Some histograms
    auto emax = 500.;
    auto enbin = 100;
    // auto etanbin = 50;
    auto psdnbin = 50;
    auto psdmax = 5.;

    TH2D * he_eta = hist_EvnEta("EvsEta", "dE/dx vs eta");
    TH2D * he_eta_1STK = hist_EvnEta("EvsEta_1STK", "dE/dx vs eta, STK 1st layer");
    TH2D * he_eta_2STK = hist_EvnEta("EvsEta_2STK", "dE/dx vs eta, STK 2 first layers");
    // TH2D * he_eta_corr = hist_EvnEta("ECorrvsEta", "dE/dx corrected vs eta");
    TH2D * he_eta_corr_my = hist_EvnEta("ECorrvsEtaMy", "dE/dx corrected vs eta (my)");
    TH2D * he_eta_He = hist_EvnEta("EvsEta_He", "dE/dx vs eta, 1.5<PSD<2.5");

    TH1D * hpsd_energy = new TH1D("PSDenergy", "PSD energy deposit", 100, 0., 10.);
    TH2D * hpsd_stk = new TH2D("PSDvsSTK", "STK first layer vs PSD mean", psdnbin, 0., psdmax, enbin, 0., emax);
    TH2D * hpsd_stk_corr = new TH2D("PSDvsSTKcorr", "STK first layer (corrected) vs PSD mean", psdnbin, 0., psdmax, enbin, 0., emax);

    // Create a TTree to use it later in TMVA
    TTree * tree = new TTree("T", "Tree for the TMVA analysis");
    Double_t psdcharge;
    Double_t psdenergy;
    Double_t psdenergy_1;
    Double_t psdenergy_2;
    tree->Branch("psd_1", &psdenergy_1, "psdenergy_1/D");
    tree->Branch("psd_2", &psdenergy_2, "psdenergy_2/D");
    Double_t stksig[12];
    Int_t ntracks;
    tree->Branch("ntracks", &ntracks, "ntracks/I");
    Int_t reac_psd;
    tree->Branch("reac_psd", &reac_psd, "reac_psd/I");
    for (int ilayer = 0; ilayer<12; ilayer++) {
        stringstream name, type;
        name << "stksig" << ilayer;
        type << "stksig" << ilayer << "/D";
        tree->Branch(name.str().c_str(), stksig + ilayer, type.str().c_str());
    }
    Double_t stop_z;
    tree->Branch("stop_z", &stop_z, "stop_z/D");
    Double_t max_dist; // Max distance btw impact points of tracks in first layer of STK
    tree->Branch("max_dist", &max_dist, "max_dist/D");
    Double_t bgo_energy;
    tree->Branch("bgo_energy", &bgo_energy, "bgo_energy/D");
    bool het;
    tree->Branch("het", &het, "het/b");

    // Write some event monitors
    bool interesting_event;

    // for (int event=0; event<10; event++) {
    for (int event=0; event<NEvents; event++) {
        interesting_event = false;

        if ((event % (NEvents / 10 + 1)) == int(NEvents / 10)) {
            int percentage = 10. * int((event + 1) / int(NEvents / 10));
            cout << "Processing percentage: " << percentage << "% \n";
        }

        DmpEvent *pev = t->GetDmpEvent();
        // is it flight data or Monte-Carlo?
        bool mc = (pev->pEvtSimuHeader() != nullptr);

        // High Energy Trigger
        het = pev->pEvtHeader()->GeneratedTrigger(3);
        // if (!het) continue;

        DmpEvtPsdRec *psdRec = pev->pEvtPsdRec();

        // STK helper
        DmpStkTrackHelper * stk_helper = new DmpStkTrackHelper(pev->GetStkKalmanTrackCollection (),
            true, pev->pEvtBgoRec(), pev->pEvtBgoHits());
        stk_helper -> SortTracks(4, true);
        ntracks = stk_helper->GetSize();

        if (ntracks == 0) continue;

        // There was reaction in PSD?
        // reac_psd = reaction_in_psd(stk_helper, psdRec);

        if(mc) {
            for (int st=0; st<pev->NSimuTrajectory(); st++) {
                DmpSimuTrajectory * simu_trajectory = pev->pSimuTrajectory(st);
                if(simu_trajectory -> parentID == 0)
                    stop_z = simu_trajectory -> stop_z;
            }
        }

        bgo_energy = pev->pEvtBgoRec()->GetTotalEnergy();

        max_dist = 0.;

        int ibar1_tck0, ibar2_tck0;

        //Loop over tracks
        for(int tck = 0; tck < stk_helper->GetSize(); tck++) {
            DmpStkTrack* stktrack = stk_helper->GetTrack(tck);

            if ((tck < stk_helper->GetSize() - 1) && (stktrack->getImpactPointPlane() == 0)) {
                for (int tck2 = tck+1; tck2 < stk_helper->GetSize(); tck2++) {
                    DmpStkTrack * stktrack2 = stk_helper->GetTrack(tck2);
                    if (stktrack->getImpactPointPlane() == 0) {
                        TVector3 ip1 = stktrack ->getImpactPoint();
                        TVector3 ip2 = stktrack2->getImpactPoint();
                        Double_t d = pow(ip1.x() - ip2.x(), 2) + pow(ip1.y() - ip2.y(), 2);
                        d = sqrt(d);
                        if(d > max_dist) max_dist = d;
                    }
                }
            }

            if (stktrack->getNhitXY() < MINTRACKHITS) continue;
            // cout << "\tTrack number " << tck << endl;

            double costheta = stktrack->getDirection().CosTheta();

            int * ibar1 = new int(0);
            int * ibar2 = new int(0);

            psdenergy = psdEnergy(stktrack, psdRec, ibar1, ibar2, mc);
            psdcharge = psdCharge(psdenergy);
            hpsd_energy -> Fill(psdenergy);

            if (tck == 0) {
                psdenergy_1 = psdEnergy(stktrack, psdRec, ibar1, ibar2, mc, 1);
                psdenergy_2 = psdEnergy(stktrack, psdRec, ibar1, ibar2, mc, 2);
            }

            delete ibar1, ibar2;



            TClonesArray * stkclusters = pev->GetStkSiClusterCollection();

            // Loop over clusters
            for (int ipoint=0; ipoint<stktrack->GetNPoints(); ipoint++) {
                for (int ixy=0; ixy<2; ixy++) {
                    // Check if cluster is associated to a hit
                    DmpStkSiCluster* cluster;
                    if(ixy == 0) 
                        cluster = stktrack -> GetClusterX(ipoint, stkclusters);
                    else
                        cluster = stktrack -> GetClusterY(ipoint, stkclusters);
                    if(!cluster) continue;
                    
                    if(is_cluster_bad_channel(cluster, badchannels)) continue;

                    double e = cluster -> getEnergy() * costheta;
                    double eta = CalcEta(cluster);
                    // double e_corr = stk_helper->CorrEnergy(stktrack, cluster);
                    double e_corr_my = etacorr.corrEnergy(stktrack, cluster);

                    stksig[ipoint*2 + ixy] = e_corr_my;

                    he_eta -> Fill(eta, e);
                    // he_eta_corr -> Fill(eta, e_corr);
                    he_eta_corr_my -> Fill(eta, e_corr_my);

                    if (ipoint == 0) {
                        hpsd_stk -> Fill(psdcharge, e);
                        hpsd_stk_corr -> Fill(psdcharge, e_corr_my);

                        he_eta_1STK -> Fill(eta, e);
                        he_eta_2STK -> Fill(eta, e);
                    }
                    if (ipoint == 1)
                        he_eta_2STK -> Fill(eta, e);

                    if (psdcharge > 1.5 && psdcharge < 2.5) {
                        he_eta_He -> Fill(eta, e);
                        // if (e_corr_my > 40. && e_corr_my < 90.) {
                        //     // we are at the proton peak in STK,
                        //     // but at the helium peak in PSD
                        //     interesting_event = true;
                        // }
                    }
                }
            } // end loop clusters
        } // end loop tracks
        if (ntracks == 1) reac_psd = 2;

        tree->Fill();
        delete stk_helper;

        if (interesting_event) {
            cout << "This event is interesting! ";
            if (!mc) {
                DmpEvtHeader * header = pev->pEvtHeader();
                cout << header->GetSecond() << 
                 " " << header->GetMillisecond() << endl;
            }
            else {
                cout << event << endl;
            }
        }
    } // end loop events

    // check fitting
    // EtaCorr etacorr_new;
    // etacorr_new.fit_pr_He(he_eta_corr_my);
    // TH1D * slices_after_corr[nbins];
    // for (int ibin=1; ibin<=nbins; ibin++) {
    //     string hname = Form("EvsEta_%d_corr", ibin);
    //     slices_after_corr[ibin - 1] = he_eta_corr_my->ProjectionY(hname.c_str(), ibin, ibin);
    //     slices_after_corr[ibin - 1]->GetListOfFunctions()->Add(etacorr_new.getFitFunction(ibin - 1));
    // }

    // Write output
    newrootfile->cd();
    newrootfile->Write();
    if(!datfile_exists) {
        for (int i=0; i<nbins; i++) {
            slices[i]->Write();
            // slices_after_corr[i]->Write();
        }
    }
    newrootfile->Close();
    return 0;
}