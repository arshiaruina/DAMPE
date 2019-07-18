{
   gSystem->Load("libDmpEvent.so");
   //TFile* f =new TFile("beam_test_file.root","READ");
   TFile* f =new TFile("/beegfs/dampe/prod/MC/reco/v6r0p0/allProton-v6r0p4_1GeV_100GeV_FTFP_CSUP05/allProton-v6r0p4_1GeV_100GeV_FTFP_CSUP05.noOrb.000601.reco.root","READ");
   TTree* t =f->Get("CollectionTree");
   
   // Register STK collections
   TClonesArray* stkclusters = new TClonesArray("DmpStkSiCluster");
   t->SetBranchAddress("StkClusterCollection",&stkclusters);
   DmpStkEventMetadata* stkMetadata = new DmpStkEventMetadata();
   t->SetBranchAddress("DmpStkEventMetadata", &stkMetadata);
   TClonesArray* stkladderadc = new TClonesArray("DmpStkLadderAdc");
   t->SetBranchAddress("DmpStkLadderAdcCollection", &stkladderadc);
       

   // Check if STK tracks collection exists 
   bool fStkKalmanTracksFound = false;
   for(int i=0; i<t->GetListOfBranches()->GetEntries(); i++){
      if (std::string(t->GetListOfBranches()->At(i)->GetName()) == std::string("StkKalmanTracks")){
         fStkKalmanTracksFound = true;
         break;
      }
   }

   // Register STK tracks collection
   stktracks = new TClonesArray("DmpStkTrack", 200);
   if(fStkKalmanTracksFound)
      t->SetBranchAddress("StkKalmanTracks",&stktracks);
       
    

   // Register STK Single-TRB collections (FM and EQM ladders)
   //TClonesArray* stkclusters_stktrb = new TClonesArray("DmpStkSiCluster");
   //t->SetBranchAddress("StkClusterCollection_STKTRB",&stkclusters_stktrb);   

   // Register BGO constainer
   DmpEvtBgoHits* bgohits  = new  DmpEvtBgoHits();
   t->SetBranchAddress("DmpEvtBgoHits",&bgohits);

   // Register BGO REC constainer
   DmpEvtBgoRec* bgorec  = new  DmpEvtBgoRec();
   t->SetBranchAddress("DmpEvtBgoRec",&bgorec);
   
   // Register PSD constainer
   DmpEvtPsdHits* psdhits  = new  DmpEvtPsdHits();
   t->SetBranchAddress("DmpEvtPsdHits",&psdhits);

   // Register DAMPE event header (event metadata information)
   DmpEvtHeader* evtheader = new DmpEvtHeader();
   t->SetBranchAddress("EventHeader",&evtheader);

   // Register AMS clusters
   TClonesArray* amsclusters = new TClonesArray("Cluster");
   t->SetBranchAddress("AmsClusterCollection",&amsclusters);

   // Register ANCILLARY container
   AncillaryEventIons* ancevent = new AncillaryEventIons();
   t->SetBranchAddress("AncillaryEventIons",&ancevent);



   // Event LOOP
   for(int entry=0; entry<t->GetEntries(); entry++){
   //for(int entry=0; entry<10; entry++){
      t->GetEntry(entry);
      
      // STK metadata
      printf("STK data mode = %d\n", stkMetadata->fRunMode );   
      // STK modes: 2 (COMPRESSED), 3(RAW), 5(ULD),  6(DLD)
      // ... for more information see: http://dpnc.unige.ch/SVNDAMPE/DAMPE/DmpSoftware/trunk/Event/Stk/include/DmpStkEventMetadata.h

      // Loop over STK clusters in event
      for(int i=0; i<stkclusters->GetLast()+1; i++){
         DmpStkSiCluster* cluster = (DmpStkSiCluster*)stkclusters->ConstructedAt(i);
         printf("\nSTK cluster info:\n");
         printf("   total ADC counts = %f\n",cluster->getEnergy());
         printf("   number of strips = %d\n",cluster->getNstrip());         
         printf("   signal in the first strip = %f\n", cluster->GetSignal(0));
         printf("   signal in the last  strip = %f\n", cluster->GetSignal(cluster->getNstrip()-1));
         printf("   noise  in the first strip = %f\n", cluster->GetNoise(0));
         printf("   noise  in the last  strip = %f\n", cluster->GetNoise(cluster->getNstrip()-1));
         printf("   center of gravity = %f\n",cluster->getLadderCentroid());
         printf("   ladder ID (hardware) = %d\n",cluster->getLadderHardware());
         //etc.
         // For more details on STK clusters see http://dpnc.unige.ch/SVNDAMPE/DAMPE/DmpSoftware/trunk/Event/Stk/include/DmpStkSiCluster.h
      }

      // Loop over STK adc 
      for(int i=0; i<stkladderadc->GetLast()+1; i++){
         DmpStkLadderAdc* laddderadc = (DmpStkLadderAdc*)stkladderadc->ConstructedAt(i);
         int ladderid = laddderadc->GetLadderID();
         for(int j=0; j<384; j++){
            double adc = laddderadc->GetChannelAdc(j);
         }
      }

      // Loop over STK tracks (if present)
      if(fStkKalmanTracksFound)
      {
         for(int i=0; i<stktracks->GetLast()+1; i++){
            DmpStkTrack* track = (DmpStkTrack* ) stktracks-> ConstructedAt(i);
            double track_x = track->getImpactPoint().x();
            double track_y = track->getImpactPoint().y();
            double track_z = track->getImpactPoint().z();
            double track_incl_x = track->getTrackParams().getSlopeX();
            double track_incl_y = track->getTrackParams().getSlopeY();

            // GET CLUSTERS FOR THE TRACK
            int ntrackpoints = track->GetNPoints();
            for(int p=0; p<ntrackpoints; p++){
               DmpStkSiCluster* clx = track->GetClusterX(p,stkclusters);
               DmpStkSiCluster* cly = track->GetClusterY(p,stkclusters);
     
               if(clx) printf("Cluster x exists for the point \n");
               else printf("No x cluter assigned to the track point in plane \n");
               if(cly) printf("Cluster y exists for the point \n");
               else printf("No y cluter assigned to the track point in plane \n");
            }
         }
      }

                

      // Get event time stamp
      //
      //   NOTE: due to the absence of information on the time zone where the file was recorded,
      //                the returned value of the time stamp can be one-hour different w.r.t. 
      //                original time stamp. Please cross-check the results of this funcntion 
      //                with the relevant run log file  
      //
      int sec = evtheader->GetSecond();                                        
      string timestamp = DmpStkHkeepHeader::TCtoString(sec);   
      std::cout<<"Human readable time stamp: "<<timestamp<<std::endl;
      int msec = evtheader->GetMillisecond(); 
      int timestamp_msec = sec*1000 + msec;


      //
      // Check STK trigger (internal VS external)
      //
      if( evtheader->GeneratedPeriodTrigger()) std::cout<<"STK trigger is internal! (Noise run)"<<std::endl;
      else     std::cout<<"STK trigger is NOT internal! (Cosmics run)"<<std::endl;
      // 
      //  .... the following is equivalent to the previous one
      if( evtheader->EnabledPeriodTrigger())    std::cout<<"STK trigger is internal! (Noise run)"<<std::endl;
      else   std::cout<<"STK trigger is NOT internal! (Cosmics run)"<<std::endl;
      // 
      //  .... the following is equivalent to the previous one
      if( evtheader->GeneratedTrigger(0))  std::cout<<"STK trigger is NOT internal! (Cosmics run)"<<std::endl;      
      else     std::cout<<"STK trigger is internal! (Noise run)"<<std::endl;      

      





      // Loop over STK single-trb clusters in event
      //
      //    WARNING:  The single-TRB data contains information of all 24-ladders, while only up to 6-real ladders (out of 24) are connected
      //                        Usually, only clusters that have   cluster->getLadderHardware() = 26 , 27 , 30, 31,   34, 35  come from the real ladders
      //
      //                        Please don't use GetX, GetY, GetZ methods for these clusters, the results will be wrong, since x,y,z are calculated 
      //                        assuming the STK ladders are integrated into DAMPE, while it's not the case for the considered single TRB  
      //
      for(int i=0; i<stkclusters_stktrb->GetLast()+1; i++){
         DmpStkSiCluster* cluster = (DmpStkSiCluster*)stkclusters_stktrb->ConstructedAt(i);
         printf("\nSTK cluster info:\n");
         printf("   total ADC counts = %f\n",cluster->getEnergy());
         printf("   number of strips = %d\n",cluster->getNstrip());         
         printf("   signal in the first strip = %f\n", cluster->GetSignal(0));
         printf("   signal in the last  strip = %f\n", cluster->GetSignal(cluster->getNstrip()-1));
         printf("   noise  in the first strip = %f\n", cluster->GetNoise(0));
         printf("   noise  in the last  strip = %f\n", cluster->GetNoise(cluster->getNstrip()-1));
         printf("   center of gravity = %f\n",cluster->getLadderCentroid());
         printf("   ladder ID (hardware) = %d\n",cluster->getLadderHardware());
         //etc.
         // For more details on STK clusters see http://dpnc.unige.ch/SVNDAMPE/DAMPE/DmpSoftware/trunk/Event/Stk/include/DmpStkSiCluster.h
      }

                

      // Loop over AMS clusters
      for(int i=0; i<amsclusters->GetLast()+1; i++){
         Cluster* amscluster = (Cluster*)amsclusters->ConstructedAt(i);
         printf("\nAMS cluster info:\n");
         printf("   number of strips = %f\n",amscluster->length);
         printf("   address of the first strip =%d\n", amscluster->address);
         printf("   cluster center of gravitu = %f\n", amscluster->GetCoG());
         // etc.
         // For more details on AMS clusters see http://dpnc.unige.ch/SVNDAMPE/DAMPE/DmpSoftware/trunk/Event/Ams/include/Cluster.hh


      }

      // Get BGO energy and direction
      printf("BGO total energy = %f\n", bgorec->GetTotalEnergy());
      // BGO trajectory 
      double x = bgorec->GetTrajectoryLocation2D().x();
      double y = bgorec->GetTrajectoryLocation2D().y();
      double z = bgorec->GetTrajectoryLocation2D().z();
      double incl_x = bgorec->GetSlopeXZ();
      double incl_y = bgorec->GetSlopeYZ();

      // Loop over the BGO (PSD) hits
      for(int i=0; i<bgohits->fEnergy.size(); i++){
         printf("\nBGO hit information:\n");
         printf("   hit energy: %f\n",bgohits->fEnergy[i]);
         printf("   hit x=%f\n", bgohits->GetHitX(i));
         printf("   hit y=%f\n", bgohits->GetHitY(i));
         printf("   hit z=%f\n", bgohits->GetHitZ(i));
         // etc.
         // For more details on BGO event class see http://dpnc.unige.ch/SVNDAMPE/DAMPE/DmpSoftware/trunk/Event/Bgo/include/DmpEvtBgoHits.h
      }
      

      // Ancillary detector
      printf("\nAnclillary detector information:\n");
      printf("   ADC s1      = %f\n",ancevent->s1);
      printf("   ADC s1_6db  = %f\n",ancevent->s1_6db);
      printf("   ADC s1_12db = %f\n",ancevent->s1_12db);
      printf("   ADC s2      = %f\n",ancevent->s2);
      printf("   ADC s2_6db  = %f\n",ancevent->s2_6db);
      printf("   ADC s2_12db = %f\n",ancevent->s2_12db);
      printf("   ADC si0_sg = %f\n",ancevent->si0_sg);
      printf("   ADC si0_lg = %f\n",ancevent->si0_lg);
      printf("   ADC si1_sg = %f\n",ancevent->si1_sg);
      printf("   ADC si1_lg = %f\n",ancevent->si1_lg);
      // etc.
      // For more information see http://dpnc.unige.ch/SVNDAMPE/DAMPE/DmpSoftware/trunk/Event/Ams/include/AncillaryEventIons.hh

   }
   // END OF EVENT LOOP
}
