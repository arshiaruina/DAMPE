 15 double CalcEta(DmpStkSiCluster* cluster){
 16         
 17         double sig = 0.;
 18         double sigMax1 = 0.;
 19         double sigMax2 = 0.;
 20         int nStrips = cluster->getNstrip();
 21         int stripMax1 = 0;
 22         int stripMax2 = 0; 
 23         int ch1 = 0;
 24         int ch2 = 0;
 25         double eta = 0.;
 26         
 27         // finding highest signal 
 28         for(int istrip = 0; istrip < nStrips; istrip++){
 29                 sig = cluster->GetSignal(istrip);
 30                 if(sig > sigMax1){
 31                         sigMax1 = sig; 
 32                         stripMax1 = istrip;
 33                 } 
 34         }
 35         
 36         if(nStrips==1){ // if cluster in one strip only
 37                 eta = 0.;// then eta is set to 0. (biased)
 38         }
 39         else {
 40                 // finding second highest signal
 41                 if(stripMax1==0){ // if highest signal strip is first strip of the cluster
 42                         stripMax2 = 1;
 43                         sigMax2 = cluster->GetSignal(stripMax2);
 44                 }
 45                 else if(stripMax1==nStrips-1){ // if highest signal strip is the last strip of the cluster
 46                         stripMax2 = stripMax1 - 1;
 47                         sigMax2 = cluster->GetSignal(stripMax2);
 48                 }
 49                 else {
 50                         if(cluster->GetSignal(stripMax1 - 1) > cluster->GetSignal(stripMax1 + 1)){
 51                                 stripMax2 = stripMax1 - 1;
 52                         }
 53                         else {
 54                                 stripMax2 = stripMax1 + 1;
 55                         }
 56                         sigMax2 = cluster->GetSignal(stripMax2);
 57                 }
 58                 
 59                 // compute eta
 60                 ch1 = cluster->GetChannelID(stripMax1);
 61                 ch2 = cluster->GetChannelID(stripMax2);
 62                 if(ch1 > ch2)
 63                         eta = sigMax1/(sigMax1 + sigMax2);
 64                 else    
 65                         eta = sigMax2/(sigMax1 + sigMax2);
 66         }
 67         return eta;
 68 }

 69 void InclInitialise() {
 70 
 71     for (int i = 0; i < STEPS; i++) {
 72         incl[i] = i * (MAX_INCLINATION / STEPS);
 73     }
 74 }
 75 
 76 int CalcInclIndex(DmpStkTrack* track, std::string dir){
 77 
 78     double theta = 0.;
 79 
 80     if (dir == "x")
 81         theta = atan(track->getTrackParams().getSlopeX());
 82     else if (dir == "y")
 83         theta = atan(track->getTrackParams().getSlopeY());
 84 
 85     theta *= 180./TMath::Pi(); // in deg
 86     // use abs value of theta
 87     for(int i = 0; i < STEPS; i++){
 88         if(std::abs(theta) >= incl[i] && std::abs(theta) < incl[i]+5.) {
 89             //inclination = incl[i];
 90             index = i;
 91             break;
 92         }
 93     }
 94     //return inclination;
 95     return index;
 96 }

 99     std::vector<TH2D*> hEtaEnergyVec;
100     for(int ihist = 0; ihist < 12; ihist++){
101         //std::string name = "hEtaEnergy_" + std::to_string(ihist);
102         //std::string title = "Cluster charge distribution " + std::to_string(ihist*5) + "<|#theta_{x,y}|<" + std::to_string((ihist+1)*5);
103         stringstream name, title;
104         name << "hEtaEnergy_" << ihist;
105         title << "Cluster charge distribution " << ihist*5 << "#leq|#theta_{x,y}|<" << (ihist+1)*5;
106         TH2D* hist = new TH2D(name.str().c_str(), ";#eta; #frac{dE}{dx} [ADC counts]",50,0.,1.,50,0.,400.);
107         hist->SetTitle(title.str().c_str());
108         hist->SetOption("COLZ");
109         hist->SetContour(30);
110         hist->SetStats(0);
111         hEtaEnergyVec.push_back(hist);
112     }
113 

