#include <iostream>
using std::cout; 
using std::endl; 

// Mine
#include "h22Option.h"
#include "Parameters.h"

// ROOT
#include "TApplication.h"
#include "TF1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TStyle.h"

int main(int argc, char *argv[]){

  // Setup options
  h22Options opts; 
  opts.args["INPUT"].args = "data.root";  
  opts.args["INPUT"].type = 1;  
  opts.args["INPUT"].name = "Input filename";  
  opts.args["PARS"].args = "data.pars";
  opts.args["PARS"].type = 1;
  opts.args["PARS"].name = "Parameter file";
  opts.set(argc, argv);

  // Get parameters
  Parameters pars;
  pars.loadParameters(opts.args["PARS"].args); 

  // Histos 
  TFile *file = TFile::Open(opts.args["INPUT"].args.c_str());
  TH2I *h[6];
  TF1  *upper[6];
  TF1  *lower[6];

  // Load 
  for(int i=0; i<6; i++){
    h[i]  = (TH2I*) file->Get(Form("ecSamplingFraction_%d",i)); 

    double min = h[i]->GetXaxis()->GetXmin();
    double max = h[i]->GetXaxis()->GetXmax();
    upper[i] = new TF1(Form("upper%d",i),"pol3", min, max); 
    lower[i] = new TF1(Form("lower%d",i),"pol3", min, max); 

    double am     = pars.getParameter("EL_SF_MU_A")   .getValue(i);
    double bm     = pars.getParameter("EL_SF_MU_B")   .getValue(i);
    double cm     = pars.getParameter("EL_SF_MU_C")   .getValue(i);
    double dm     = pars.getParameter("EL_SF_MU_D")   .getValue(i);
    double as     = pars.getParameter("EL_SF_SIGMA_A").getValue(i);
    double bs     = pars.getParameter("EL_SF_SIGMA_B").getValue(i);
    double cs     = pars.getParameter("EL_SF_SIGMA_C").getValue(i);
    double ds     = pars.getParameter("EL_SF_SIGMA_D").getValue(i);
    double NSIGMA = pars.getParameter("EL_EC_NSIGMA") .getValue(0);

    upper[i]->SetParameter(0,dm+NSIGMA*ds);
    upper[i]->SetParameter(1,cm+NSIGMA*cs);
    upper[i]->SetParameter(2,bm+NSIGMA*bs);
    upper[i]->SetParameter(3,am+NSIGMA*as);
    lower[i]->SetParameter(0,dm-NSIGMA*ds);
    lower[i]->SetParameter(1,cm-NSIGMA*cs);
    lower[i]->SetParameter(2,bm-NSIGMA*bs);
    lower[i]->SetParameter(3,am-NSIGMA*as);

    upper[i]->SetLineStyle(8); 
    lower[i]->SetLineStyle(8); 
    upper[i]->SetLineWidth(2); 
    lower[i]->SetLineWidth(2); 

    cout << "[Pol3 Params] UPPER, Sector=" << i << " A=" << upper[i]->GetParameter(3) << " B=" << upper[i]->GetParameter(2) << " C=" << upper[i]->GetParameter(1) << " D=" << upper[i]->GetParameter(0) << endl; 

  }

  // Launch 
  //  TApplication *app = new TApplication("MyApp", &argc, argv); 
  TCanvas *can      = new TCanvas("can","",1200,800); 

  can->Divide(3,2);
  for (int i=0; i<6; i++){
    can->cd(i+1);
    h[i]->Draw("colz");
    upper[i]->Draw("same");
    lower[i]->Draw("same");

    gPad->SetGridx(); 
    gPad->SetGridy(); 

    gStyle->SetOptStat(0); 
  }
  
  can->Print("img/SamplingFraction.png");
  
  //  app->Run();

  return 0; 
}
