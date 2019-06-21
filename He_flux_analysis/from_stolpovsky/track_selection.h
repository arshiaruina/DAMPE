#include <string>
#include <stdio.h> 
#include <cmath>
#include <iostream>
#include "DmpStkSiCluster.h"
#include "DmpStkLadderAdc.hh"

#define NLADDERS     192
#define NCHANNELS    384
#define NLAYERS      6
#define VA           6
#define ETA          2

// Define track cuts
#define MINTRACKHITS 6 // 5
#define MINCHISQ     6.0
#define CHI2KALMANCUT 5.

int const h_nstrips = 5;
//int const h_nparticles = 6;
int const h_nparticles = 6;
int const h_nplanes = 6;
int const nlayers_psd = 4; 
//int const NANGLES = 12;
int const NANGLES = 12;

int const Inclination =24;
//int const NPARTICLES =7; //until O w/o N //for now only proton and Helium, only P to test VA eq
int const NPARTICLES =2; //until O w/o N //for now only proton and Helium, only P to test VA eq
int const ETABIN = 50;//we fixed the number of bins in eta to 50 (~5mum bin)
int const FillEneDist =1;//0 fill VA plots, 1 not
int const SelType =1;//to increase statistics 0 select with || in thetax, 1 && (select holed cones in thetax and thetay

//int nEBIN[NPARTICLES]  = {200,600,1000,3000,4000,6000};//low bin after Z =2 too speed up the job running, since the tree is filled anyhow
int nEBIN[NPARTICLES]  = {200,600};//low bin after Z =2 too speed up the job running, since the tree is filled anyhow
//int nEBIN[NPARTICLES]  = {200,600};//low bin after Z =2 too speed up the job running, since the tree is filled anyhow
//int nEBINS[NPARTICLES] = {200,600,1000,3000,4000,6000};//we want histos only until C
int nEBINS[NPARTICLES] = {200,600};//we want histos only until C
//int nEBINS[NPARTICLES] = {200,600};
int const MAXACC = 60;
int const NANGLESMAX = 12;
//acceptance

//FillEneDist 0 fills ene dist in each VA
//1 does not fill ene dist
//2 fills ene dist with VA + eta corr



using namespace std;

bool  is_cluster_bad_channel(int clusterladder, int firststrip, int laststrip, bool** badchannels){ 
	for(int i=firststrip; i<=laststrip; i++){
		if(badchannels[clusterladder][i]) return true;
	}
	return false;
}


bool** read_bad_channels_file(const char* filename){
	bool** badchannels = new bool*[NLADDERS];
	for(int i=0; i<NLADDERS; i++){
		badchannels[i] = new bool[NCHANNELS];
		for(int j=0; j<NCHANNELS; j++){
			badchannels[i][j] = false;
		}
	}

	ifstream infile(filename);	
	if (!infile)
	{
		std::cout<<"Can't open Bad Channels file: "<<filename<<" ==> throwing exception!"<<std::endl;
		throw;
	}
	std::string line;
	while (getline(infile, line)) {
		if (line.c_str() == std::string("")) break;
		std::istringstream lineStream(line.c_str());
		std::string number;
		int i = 0;
		int ildr   = -1;
		int bdchnl = -1;
		while (getline(lineStream, number, ',')) {
			switch (i++) {
			case 0:
				ildr = atoi(number.c_str());
				break;
			case 1:
				bdchnl = atoi(number.c_str());
				break;			
			}
		}
		if(ildr>=0 && bdchnl>=0 && ildr<NLADDERS && bdchnl<NCHANNELS){
			badchannels[ildr][bdchnl] = true;
		}
	}
	std::cout<<"Done reading bad channels"<<std::endl;
	return badchannels;
}


//void get_linear_fit(double* x, double* y, int npoints,  double& a){
void get_linear_fit(double* x, double* y, int npoints,  double& a, double& b){
        // #@
        // #@  y = ax + b
        // #@
        
        double n   = 0;
        double Ex  = 0;
        double Ex2 = 0;
        double Ey  = 0;
        double Exy = 0;
        for(int i=0; i<npoints; i++){
		double argument = x[i];
		double value    = y[i];
		Ex+= argument;
		//Ex2+=argument**2;
		Ex2+=argument*argument;
		Ey+=value;
		Exy+=argument*value;
		n++;
	}
        
        /* a = (n*Exy - Ex* Ey )/(n*Ex2 - Ex**2); */
        /* b = (Ex * Exy - Ey * Ex2) / (Ex**2 - n * Ex2); */

	a = (n*Exy - Ex* Ey )/(n*Ex2 - Ex*Ex);
        b = (Ex * Exy - Ey * Ex2) / (Ex*Ex - n * Ex2);

}

void get_residuals(double* args, double* vals ,int npoints, double a, double b, double* residuals){
	for(int i=0; i<npoints; i++){
		double arg = args[i];
		double val = vals[i];
		double projval = a*arg + b;
		double diff  = val - projval;
		residuals[i] = diff;
	}
}


int whichParticle(double sqrtS, int x_or_y){ //// int  = 0 for X and int = 1 for Y

  /////////////////// Fit results preliminary min6points no Kalmanchi2 //////////////////
  
  // X particle 1 chi2/ndf 1.94395 mean  1.13824 sigma 0.0966518
  // Y particle 1 chi2/ndf 1.48315 mean  1.12822 sigma 0.0926166
  // X particle 2 chi2/ndf 1.08495 mean  2.34962 sigma 0.28228
  // Y particle 2 chi2/ndf 0.818483 mean  2.29956 sigma 0.293626
  /////////////////////////////////////////////////////////////////



  /////////////////// Fit results preliminary min5points with Kalmanchi2 //////////////////

  // 0 EandP X: mean 0.941922 sigma 0.0555021 chi^2/NDF 121.183
  // 0 EandP Y: mean 0.93563 sigma 0.0507962 chi^2/NDF 170.469
  // 1 He X: mean 1.99252 sigma 0.15482 chi^2/NDF 10.6068
  // 1 He Y: mean 1.97135 sigma 0.144205 chi^2/NDF 17.9441
  // 2 Be X: mean 4.45259 sigma 0.361992 chi^2/NDF 1.71093 R: 
  // 2 Be Y: mean 4.40607 sigma 0.336888 chi^2/NDF 1.62
  // 3 B  X: mean 5.70487 sigma 0.371221 chi^2/NDF 1.19621
  // 3 B  Y: mean 5.58466 sigma 0.309454 chi^2/NDF 1.65823
  // 4 C  X: mean 6.85596 sigma 0.372864 chi^2/NDF 1.54203 !!!Note: C is computed on Smean_x_hc
  // 4 C  Y: mean 6.67094 sigma 0.370338 chi^2/NDF 2.24394
  // 5 O  X: mean 9.09751 sigma 0.570488 chi^2/NDF 1.30016
  // 5 O  Y: mean 8.92516 sigma 0.552746 chi^2/NDF 1.06529

  /////////////////////////////////////////////////////////////////

  int particle_type = -9;
  double nsigma = 2.0;


  /* double mean_x[h_nparticles] = {1.02417, 2.18264}; */
  /* double mean_y[h_nparticles] = {1.01769, 2.14549}; */
  
  /* double sigma_x[h_nparticles] = {0.0827143,0.261576}; */
  /* double sigma_y[h_nparticles] = {0.0764869,0.261596}; */
  //until O
  // double mean_x[h_nparticles] = {0.941922, 1.99252, 2.90279, 4.45259, 5.70487, 6.85596, 9.09751 };
  //  double mean_y[h_nparticles] = {0.93563,  1.97135, 2.86109, 4.40607, 5.58466, 6.67094, 8.92516 };
  //let's put 2sigma to detect Z =1 particles
  // double sigma_x[h_nparticles] = {2*0.0555021, 0.15482, 1.51186e-01,  0.361992, 0.371221, 0.372864, 0.570488};
  //double sigma_y[h_nparticles] = {2*0.0507962, 0.144205, 1.50847e-01, 0.336888, 0.309454, 0.370338, 0.552746};

  //only until C
  double mean_x[h_nparticles] = {0.941922, 1.99252, 2.90279, 4.45259, 5.70487, 6.85596 };
  double mean_y[h_nparticles] = {0.93563,  1.97135, 2.86109, 4.40607, 5.58466, 6.67094};
  //let's put 2sigma to detect Z =1 particles
  double sigma_x[h_nparticles] = {2*0.0555021, 0.15482, 1.51186e-01,  0.361992, 0.371221, 0.372864};
  double sigma_y[h_nparticles] = {2*0.0507962, 0.144205, 1.50847e-01, 0.336888, 0.309454, 0.370338};
 
  //  double sigma_x[h_nparticles] = {0.0555021, 0.15482,  0.361992, 0.371221, 0.372864, 0.570488};
  // double sigma_y[h_nparticles] = {0.0507962, 0.144205, 0.336888, 0.309454, 0.370338, 0.552746};

  //x
  if(x_or_y == 0){
    for(int ip = 0; ip < h_nparticles; ip++){
      if( (sqrtS > mean_x[ip] - nsigma*sigma_x[ip]) &&   (sqrtS < mean_x[ip] + nsigma*sigma_x[ip]))
	particle_type =ip+ 1;
      
    }
  }
  
  //y
  else if(x_or_y == 1){
    for(int ip = 0; ip < h_nparticles; ip++){
      if( (sqrtS > mean_y[ip] - nsigma*sigma_y[ip]) &&   (sqrtS < mean_y[ip] +  nsigma*sigma_y[ip])){
	particle_type =ip+ 1;
      }
      
    }
    
    
  }
  
  return particle_type;
  
  
  
}


double eta_calc(vector<double> ene_clu , vector<double> ch_clu){
  //// defined as R/(R+L)

  double eta_value = 0.0;

  // int clu_size = ene_clu.size();
  
  // if(clu_size >= 2){
  
  double max_ene =  *max_element(ene_clu.begin(),ene_clu.end());
  int max_pos = find(ene_clu.begin(),ene_clu.end(),max_ene) - ene_clu.begin();
  //cout<<" max ene "<<max_ene<<" max_pos "<<max_pos<<endl;
  
  /* cout<<endl;         */
  /* cout<<"in f "<<clu_size<<endl; */
  int ch_right = max_pos + 1;
  int ch_left  = max_pos - 1;
  double energy_r = 0.0;
  double energy_l = 0.0;
  /* double strip_r = 0.0; */
  /* double strip_l = 0.0; */
  
  
  
  /* cout<< "right ene "<<energy_r<<" ch  "<<strip_r<<"  is = to "<<ch_right<<endl; */
  /* cout<< "left ene "<<energy_l<<" ch  "<<strip_l<<"  is = to "<<ch_left<<endl; */
  
  
  /* for(int i = 0; i < clu_size; i++){ */
  /* 	cout<<ene_clu.at(i)<<"  "<<ch_clu.at(i)<<endl; */
  /* } */
  
  if( find(ch_clu.begin(), ch_clu.end(), ch_right) != ch_clu.end()){
    int pos_r = find(ch_clu.begin(), ch_clu.end(), ch_right) - ch_clu.begin(); 
    energy_r = ene_clu.at(pos_r);
    //	strip_r = ch_clu.at(pos_r);
    //	cout<< "right ene "<<energy_r<<" ch  "<<strip_r<<"  is = to "<<ch_right<<endl; 
  }
  
  
  if( find(ch_clu.begin(), ch_clu.end(), ch_left) != ch_clu.end()){
    int pos_l = find(ch_clu.begin(), ch_clu.end(), ch_left) - ch_clu.begin(); 
    energy_l = ene_clu.at(pos_l);
    // strip_l = ch_clu.at(pos_l);
    //	cout<< "left ene "<<energy_l<<" ch  "<<strip_l<<"  is = to "<<ch_left<<endl;
  }
  
  /* cout<< "right ene "<<energy_r<<" ch  "<<strip_r<<"  is = to "<<ch_right<<endl; */
  /* cout<< "left ene "<<energy_l<<" ch  "<<strip_l<<"  is = to "<<ch_left<<endl; */
  
  
  if(energy_r > energy_l){
    eta_value = energy_r/(energy_r+max_ene);
    //	cout<<" numeratore "<<energy_r<<endl;
  }
  else if(energy_l > energy_r){
    eta_value = max_ene/(energy_l+max_ene);
    //	cout<<" numeratore "<<max_ene<<endl;
  }
  // cout<<" eta "<<eta_value<<endl;
  
  //}//close if 2 clusters
// else
//  eta_value = -9;



return eta_value;
}//close etacalc


//////////////// TRY THE CLUSTER ALGORITHM: ONLY THE CLUSTER WITH HIGHEST CHARGE BELONGING TO THE TRACK IS TAKEN 
//////////////// IT'LL BE OBSOLATE WITH THE NEW REPROCESSING? 


vector<double>  find_sub_cluster(vector<double> event){ //
  int  nchannels = event.size();
  //// cluster parameter definition ////
 
  
  
  // cout<<" in loop event size "<<event.size()<<endl;
  vector<double> sub_cluster;
  vector<double> sub_cluster_ch;
  sub_cluster.clear();
  sub_cluster_ch.clear();
   
 
  double max_in_event;   ///// find max in the event
  max_in_event =  *max_element(event.begin(),event.end());
  int canale_max = find(event.begin(),event.end(),max_in_event) - event.begin();
    
  int high_strip = canale_max + 1;
  int low_strip  = canale_max;
  
  if(max_in_event >= 45. && event.size()!= 1){ /// max_in_event > T_cut * sigma[canale_max]  /// && sigma[canale_max] < t_sigma //// the clusters have been already seleceted and belong to the tracks
    //    cout<<" in loop max in event  "<<max_in_event<<endl;
    sub_cluster.push_back(max_in_event);		
    sub_cluster_ch.push_back(canale_max);	
    
    // cout<<"high strip "<<high_strip<<" low strip "<<low_strip<<endl;
    /* ///// loop on smaller channels */
    for(int lloop = 0; lloop < low_strip; lloop++){
      //cout<<" loop "<<lloop<<endl;
      int canale_low = low_strip - lloop -1;
      int canale_low_check = low_strip - lloop ;
      //      cout<< canale_low<<" ??? "<<canale_low_check<<endl;
      //   if(!event.at(canale_low)) continue;
      //    cout<<" canale low "<<canale_low<<endl;
      if (event.at(canale_low) >= 15. && event.at(canale_low) <= event.at(canale_low_check)){ //t_cut * sigma[canale_low] && event[canale_low] < event[canale_low_check]&& sigma[canale_low] < t_sigma
	  sub_cluster.push_back(event[canale_low]);
	  sub_cluster_ch.push_back(canale_low);
	  //  cout<<" canale low "<< event[canale_low]<<endl;
	  
      }/// if lloop
      else
	break;
    }/// for llopp
    
    /* cout<<" in loop event size again! "<<event.size()<<endl; */
    /* ///// loop on higher channels */
    int high_range = nchannels - high_strip;
    //  cout<<" high range "<<high_range<<endl;
    for(int hloop = 0; hloop < high_range; hloop++){
      // cout<<hloop<<endl;
      int canale_high = high_strip + hloop;
      int canale_high_check = high_strip + hloop -1;
      //  cout<<"canale high: "<<canale_high<<endl;
      //   if(!event.at(canale_high)) continue;
      if (event.at(canale_high) >= 15. && event.at(canale_high) <= event.at(canale_high_check) ){ // event[canale_high] > t_cut *sigma[canale_high] && event[canale_high] < event[canale_high_check] && sigma[canale_high] < t_sigma
      	sub_cluster.push_back(event[canale_high]);
      	sub_cluster_ch.push_back(canale_high);
	//	cout<<" canale high ADC "<< event[canale_high]<<endl;
	
      }/// if hloop
      else
      	break;
    }/// for hloop
  }
  else if(event.size() == 1)
    sub_cluster.push_back(event.at(0));
 
  event.clear();
  //  cout<<" in loop sub cluster size "<<sub_cluster.size()<<endl;
    /* else */
    /*   break; */
  
  return sub_cluster;

}

void findinclination(double thetaxy, int &incl){
  
  double angle[NANGLESMAX]={0.};
  for(int iangle =0; iangle <NANGLESMAX; iangle++){

    if(iangle ==0){
      angle[iangle]= MAXACC/NANGLESMAX;
      
      //  cout << "angle  " << angle[iangle] << endl;
    }//close if

    else if(iangle >0){
      angle[iangle] = angle[iangle-1] + angle[0];
      //  cout << "angle  " << angle[iangle] << endl;
      // cout << iangle << endl;
    }//close if
    
  }//close iangle for

  for(int inclin =0; inclin < NANGLESMAX; inclin++ ){
    if(inclin==0){
      if(thetaxy >= angle[0] && thetaxy <= angle[1]){
	incl = inclin;
      }//close if
    }//
    else if(inclin >0){
      if(thetaxy >= angle[inclin-1] && thetaxy <= angle[inclin] ){
	incl = inclin;
      }//close if thetaxy
    }//close if inclin
  }//close for inclin

}//close findinclination



double calcEta(DmpStkSiCluster* x_cluster, int &chout){

  // initialize
  double Sx;
  double Smax1x = 0;
  double Smax2x = 0;
  int GetNstripX = x_cluster->getNstrip();
  int chSmax1x = x_cluster->getIndex1();
  int ch1x = chSmax1x;
  int chSmax2x = x_cluster->getIndex1();
  int ch2x = chSmax2x;
  double EtaX;
  
  //find first higher signal
  for(int istripx = 0; istripx < GetNstripX; istripx++){
    Sx = x_cluster->GetSignal(istripx);
    if(Sx>Smax1x && istripx == 0){
      Smax1x = Sx;
      ch1x = chSmax1x;
      chSmax1x++;
    }//close if
    else if(Sx > Smax1x && istripx >0) {  
      Smax1x = Sx;
      ch1x = chSmax1x;
      chSmax1x++;
    }//close if
    else chSmax1x++;
  }//close istripx for
  
  //find second higher signal
  for(int istripx = 0; istripx < GetNstripX; istripx++){
    Sx = x_cluster->GetSignal(istripx);
    if(Sx>Smax2x && abs(Sx-Smax1x)>0.001 &&  istripx==0){
      Smax2x = Sx;
      ch2x= chSmax2x;
      chSmax2x++;
    }//close if
    else if(Sx>Smax2x && abs(Sx-Smax1x)>0.001 &&  istripx>0 ){
      Smax2x = Sx;
      ch2x= chSmax2x;
      chSmax2x++;
    }//close if
    else chSmax2x++;
  }//close istripx for
  
  //COMPUTE ETAx
  if(ch1x < ch2x){
    EtaX = Smax2x/(Smax1x + Smax2x);
    //  cout << "Eta X first "<< EtaX << endl;
    
  }//close if
  else {
    EtaX = Smax1x/(Smax1x + Smax2x);
    // cout << "Eta X second "<< EtaX << endl;
  }//close else
  
  //END CALC ETA X//////////////////
  chout = ch1x;
  return EtaX;

}//close calcEta


double calcEtaNoCorr(DmpStkSiCluster* x_cluster,TClonesArray* ladderColl,int &chout){

  // initialize
  double Sx;
  double Smax1x = 0;
  double Smax2x = 0;
  int GetNstripX = x_cluster->getNstrip();
  int chSmax1x = x_cluster->getIndex1();
  int ch1x = chSmax1x;
  int chSmax2x = x_cluster->getIndex1();
  int ch2x = chSmax2x;
  double EtaX;
  
  //find first higher signal
  for(int istripx = 0; istripx < GetNstripX; istripx++){
    Sx = x_cluster->GetSignal(istripx);
    if(Sx>Smax1x && istripx == 0){
      Smax1x = Sx;
      ch1x = chSmax1x;
      chSmax1x++;
    }//close if
    else if(Sx > Smax1x && istripx >0) {  
      Smax1x = Sx;
      ch1x = chSmax1x;
      chSmax1x++;
    }//close if
    else chSmax1x++;
  }//close istripx for
  
  //find second higher signal
  for(int istripx = 0; istripx < GetNstripX; istripx++){
    Sx = x_cluster->GetSignal(istripx);
    if(Sx>Smax2x && abs(Sx-Smax1x)>0.001 &&  istripx==0){
      Smax2x = Sx;
      ch2x= chSmax2x;
      chSmax2x++;
    }//close if
    else if(Sx>Smax2x && abs(Sx-Smax1x)>0.001 &&  istripx>0 ){
      Smax2x = Sx;
      ch2x= chSmax2x;
      chSmax2x++;
    }//close if
    else chSmax2x++;
  }//close istripx for
  
  //COMPUTE ETAx
  if(ch1x < ch2x){
    EtaX = Smax2x/(Smax1x + Smax2x);
    //  cout << "Eta X first "<< EtaX << endl;
    
  }//close if
  else {
    EtaX = Smax1x/(Smax1x + Smax2x);
    // cout << "Eta X second "<< EtaX << endl;
  }//close else
  
  //END CALC ETA X//////////////////
  chout = ch1x;
  return EtaX;

}//close calcEta
