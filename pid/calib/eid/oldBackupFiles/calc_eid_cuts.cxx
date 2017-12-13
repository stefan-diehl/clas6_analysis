//////////////////////////////////////
/*
 
 David Riser, University of Connecticut
 calc_eid_cuts.cxx -> test the usage of libraries
 July 9, 2016

 Edited : Nov 5, 2016
 
 */
//////////////////////////////////////

// c++ includes
#include <iostream>
#include <vector>
using namespace std; 

// my includes
#include "CommonTools.h"
#include "h22Event.h"
#include "h22Option.h"
#include "h22Reader.h"
#include "Parameters.h"
#include "ParameterSet.h"

// root includes
#include "TCanvas.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TLine.h"
#include "TVector3.h"

class EIDCutGenerator : public h22Reader{
public:
  EIDCutGenerator(){
  }
  ~EIDCutGenerator(){
  }

protected:
  const static int N_SLICES = 100;

  TH1F *h1_vz[6];
  TH1F *h1_nphe[6];
  TH1F *h1_ecu[6];
  TH1F *h1_ecv[6];
  TH1F *h1_ecw[6];
  TH1D *h1_slice[6][N_SLICES];
  TH2F *h2_ec_sampling[6];

  TF1 *f_ecSamplingMean_pol3[6];
  TF1 *f_sigma_pol3[6];
  TGraphErrors *g_ecSamplingMean[6];
  TGraphErrors *g_sigma[6];

  Parameters eidParameters; 

  double P_MIN, P_MAX, P_STEP;
  int P_BINS;
  
  double ecSamplingMean[N_SLICES];
  double sigma[N_SLICES];
  double ecSamplingMeanError[N_SLICES];
  double sigma_e[N_SLICES];
  double x[N_SLICES];
  double dx[N_SLICES];
  
public:
  void Initialize();
  void DoECSamplingCut();
  void Loop(int numberOfEvents);
  void ProcessEvent();
  void Print(); 
  void Save();
};

void EIDCutGenerator::Initialize(){
  
  P_MIN  = 0.5;
  P_MAX  = 2.80;
  P_BINS    = N_SLICES;
  P_STEP = (P_MAX-P_MIN)/(N_SLICES-1);
  
  for (int s=0; s<6; s++){
      h1_vz[s]   = new TH1F(Form("h1_vz_%d",s),Form(" Z Vertex Sector %d",s+1),400,-85,0);
      h1_nphe[s] = new TH1F(Form("h1_nphe_%d",s),Form(" N. Photo-electrons Sector %d",s+1),400,0,100);
      h1_ecu[s] = new TH1F(Form("h1_ecu_%d",s),Form(" EC U Coord. Sector %d",s+1),400,0,500);
      h1_ecv[s] = new TH1F(Form("h1_ecv_%d",s),Form(" EC V Coord. Sector %d",s+1),400,0,500);
      h1_ecw[s] = new TH1F(Form("h1_ecw_%d",s),Form(" EC W Coord. Sector %d",s+1),400,0,500);
      h2_ec_sampling[s] = new TH2F(Form("h2_ec_sampling_%d",s),Form("EC Sampling Fraction Sector %d",s+1),P_BINS,P_MIN,P_MAX,100,0.05,0.5);
    }
  
  Init();
}

void EIDCutGenerator::Loop(int numberOfEvents){

    for (int iev=0; iev<numberOfEvents; iev++){
        GetEntry(iev);
	ProcessEvent();
        if (iev%10000 == 0) { std::cout << "\r Doing event " << iev << std::flush; }
    }

}

void EIDCutGenerator::ProcessEvent(){
  TVector3 uvw = event.GetUVWVector(0);
  int sector = event.ec_sect[0]-1;
  if (sector > -1 && event.q[0] < 0){
      h1_vz[sector]          ->Fill(event.vz[0]);
      h1_nphe[sector]        ->Fill(event.nphe[0]);
      h1_ecu[sector]         ->Fill(uvw.X());
      h1_ecv[sector]         ->Fill(uvw.Y());
      h1_ecw[sector]         ->Fill(uvw.Z());
      h2_ec_sampling[sector] ->Fill(event.p[0],event.etot[0]/event.p[0]);
    }
}

void EIDCutGenerator::Save(){
  eidParameters.saveParameters("pars.dat");
}

void EIDCutGenerator::DoECSamplingCut(){
    
    TF1 *f_gauss = new TF1("f_gauss","gaus",0.25,0.45);    
    
    // For passing to TGraph constructor
    for (int b=0; b<N_SLICES; b++) {x[b] = (double) b*P_STEP + P_MIN; dx[b] = 0.0;}
    
    std::cout << "\n Done event loop! " << std::endl;
    
    // Do the fitting and slicing
    for (int s=0; s<6; s++){
        for (int b=0; b<N_SLICES; b++){
            std::string name       = Form("h1_slice_%d_%d",s,b);
            h1_slice[s][b]    = new TH1D(name.c_str(),name.c_str(),100,0.05,0.5);
            h2_ec_sampling[s] ->ProjectionY(name.c_str(),b+1,b+2);
            h1_slice[s][b]    ->Fit("f_gauss","RQ");
            
            ecSamplingMean[b] = f_gauss->GetParameter(1);
            sigma[b] = f_gauss->GetParameter(2);
            ecSamplingMeanError[b] = f_gauss->GetParError(1);
            sigma_e[b] = f_gauss->GetParError(2);
        }
        
        f_ecSamplingMean_pol3[s]  = new TF1(Form("f_ecSamplingMean_pol3_%d",s),"pol3",0,N_SLICES);
        f_sigma_pol3[s] = new TF1(Form("f_sigma_pol3_%d",s),"pol3",0,N_SLICES);
        g_ecSamplingMean[s] = new TGraphErrors(N_SLICES,x,ecSamplingMean,dx,ecSamplingMeanError);
        g_sigma[s] = new TGraphErrors(N_SLICES,x,sigma,dx,sigma_e);
        g_ecSamplingMean[s] ->Fit(Form("f_ecSamplingMean_pol3_%d",s),"RQ");
        g_sigma[s] ->Fit(Form("f_sigma_pol3_%d",s),"RQ");

	cout << "ec sampling mean d = " << f_sigma_pol3[s]->GetParameter(0) << endl; 
	cout << "ec sampling mean c = " << f_sigma_pol3[s]->GetParameter(1) << endl; 
	cout << "ec sampling mean b = " << f_sigma_pol3[s]->GetParameter(2) << endl; 
	cout << "ec sampling mean a = " << f_sigma_pol3[s]->GetParameter(3) << endl; 

    }
}

void EIDCutGenerator::Print(){

    TCanvas * c1 = new TCanvas("c1", "", 1200, 400);
    c1->Divide(3,1);
    c1->Print("test.pdf[");
    

    TLine line;
    line.SetLineColor(kRed);
    

    for (int s=0; s<6; s++){
        c1->cd(1);
        h1_ecu[s]->Draw();
	//        line.DrawLine(pars.ECUMIN[s],0,pars.ECUMIN[s],h1_ecu[s]->GetMaximum());
	//        line.DrawLine(pars.ECUMAX[s],0,pars.ECUMAX[s],h1_ecu[s]->GetMaximum());
        c1->cd(2);
        h1_ecv[s]->Draw();
	//        line.DrawLine(pars.ECVMIN[s],0,pars.ECVMIN[s],h1_ecv[s]->GetMaximum());
	//        line.DrawLine(pars.ECVMAX[s],0,pars.ECVMAX[s],h1_ecv[s]->GetMaximum());
        c1->cd(3);
        h1_ecw[s]->Draw();
	//        line.DrawLine(pars.ECWMIN[s],0,pars.ECWMIN[s],h1_ecw[s]->GetMaximum());
	//        line.DrawLine(pars.ECWMAX[s],0,pars.ECWMAX[s],h1_ecw[s]->GetMaximum());
        c1->Print("test.pdf");
    }
    
    c1->cd(1);
    h2_ec_sampling[0]->Draw("colz"); g_ecSamplingMean[0]->Draw("same");
    c1->cd(2);
    h2_ec_sampling[1]->Draw("colz"); g_ecSamplingMean[1]->Draw("same");
    c1->cd(3);
    h2_ec_sampling[2]->Draw("colz"); g_ecSamplingMean[2]->Draw("same");
    c1->Print("test.pdf");
    
    c1->cd(1);
    h2_ec_sampling[3]->Draw("colz"); g_ecSamplingMean[3]->Draw("same");
    c1->cd(2);
    h2_ec_sampling[4]->Draw("colz"); g_ecSamplingMean[4]->Draw("same");
    c1->cd(3);
    h2_ec_sampling[5]->Draw("colz"); g_ecSamplingMean[5]->Draw("same");
    c1->Print("test.pdf");
    
    c1->Clear();
    c1->Divide(6,1);    
    
    for (int s=0; s<6; s++) { c1->cd(s+1); h1_nphe[s]->Draw(); } c1->Print("test.pdf");
    for (int s=0; s<6; s++) { c1->cd(s+1); h1_vz[s]->Draw(); } c1->Print("test.pdf");
    
    // Print Slices
    for (int b=0; b<N_SLICES; b++){
        for (int s=0; s<6; s++){
            c1->cd(s+1);
            h1_slice[s][b]->Draw();
        }
        c1->Print("test.pdf");
    }
    c1->Print("test.pdf]");
}

int main(int argc, char * argv[]){
    
    // Setup Options
    h22Options opts;
    opts.set(argc,argv);

    int numberOfEvents = opts.args["N"].arg; 
    numberOfEvents = smallest(numberOfEvents, opts.ifiles.size());

    EIDCutGenerator *cutGenerator = new EIDCutGenerator(); 
    
    vector<string>::iterator it; 
    for (it=opts.ifiles.begin(); it<opts.ifiles.end(); it++) { cutGenerator->AddFile(*it); }

    cutGenerator->Initialize(); 
    cutGenerator->Loop(numberOfEvents); 
    cutGenerator->DoECSamplingCut(); 
    //    cutGenerator->Print(); 
    cutGenerator->Save(); 
    
    return 0;
}
