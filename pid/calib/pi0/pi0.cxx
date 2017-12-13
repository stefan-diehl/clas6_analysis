////////////////////////////////////////
/*
 David Riser, University of Connecticut
 
 */
////////////////////////////////////////

// c++ includes
#include <iostream>
using namespace std;

// my includes
#include "CommonTools.h"
#include "Corrections.h"
#include "h22Event.h"
#include "h22Option.h"
#include "h22Reader.h"
#include "GenericAnalysis.h"
#include "ParticleFilter.h"
#include "Parameters.h"
#include "ParameterSet.h"
#include "Pars.h"

// root includes
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"

class Pi0Calibration : public GenericAnalysis {

public:
  Pi0Calibration(h22Options *opts, Parameters *params) : GenericAnalysis(opts) { pars = params; }
  ~Pi0Calibration(){ }

public:
  Parameters     *pars; 
  ParticleFilter *filter;
  Corrections     corr;
  
  // Histograms
  TH1F *openingAngle[6];
  TH1F *invMass[6];  
  TH1F *energyDiff[6];  
  TH2F *openEnergy[6];
  TH1F *numberPhotons;

  void DoFits();
  void Initialize();
  void ProcessEvent();
  void Save(string outputFilename);
  void WriteHardCodedParameters();
};

void Pi0Calibration::Initialize(){

    filter = new ParticleFilter(pars);
    filter->set_info(GSIM, GetRunNumber());

    numberPhotons = new TH1F("numberPhotons","",6,0,5); 
	  
    for(int s=0; s<6; ++s){
      openingAngle[s]  = new TH1F(Form("openingAngle%d",s),"",100, 0, 180); 
      invMass[s]       = new TH1F(Form("invMass%d",s),"",     100,-0.2,0.8); 
      energyDiff[s]    = new TH1F(Form("energyDiff%d",s),"",  100,-1.5,3.5); 
      openEnergy[s]    = new TH2F(Form("openEnergy%d",s),"",100,-1.5,3.5,100,0,180); 
    }
}

void Pi0Calibration::ProcessEvent(){

  vector<int> photons = filter->getVectorOfParticleIndices(event, 22); 
  
  
  numberPhotons->Fill((double)photons.size()); 

  if (photons.size() >= 2){

    for (int iphot=0; iphot<photons.size(); iphot++){
      TLorentzVector firstPhoton = event.GetTLorentzVector(photons[iphot], 22); 
      for(int jphot=iphot+1; iphot<photons.size(); iphot++){
	
	
	TLorentzVector otherPhoton = event.GetTLorentzVector(photons[jphot], 22);
	TLorentzVector pion        = firstPhoton+otherPhoton; 
	
	int sector = floor((pion.Phi()*to_degrees+180.0)/60.0); 
	
	invMass[sector]     ->Fill(pion.Mag());
	openingAngle[sector]->Fill(firstPhoton.Angle(otherPhoton.Vect())*to_degrees); 
	energyDiff[sector]  ->Fill(firstPhoton.E()-otherPhoton.E());
	openEnergy[sector]  ->Fill(firstPhoton.E()-otherPhoton.E(), firstPhoton.Angle(otherPhoton.Vect())*to_degrees); 
      }
    }
    
  }

}

void Pi0Calibration::DoFits(){
  TF1 *fitGauss = new TF1("fitGauss","gaus",0.1,0.2); 
  
  double mu[6], sigma[6]; 
  ParameterSet muPar, sigmaPar, nsigmaPar; 
  muPar.setName("PI0_INVMASS_MU"); 
  sigmaPar.setName("PI0_INVMASS_SIGMA");
  nsigmaPar.setName("PI0_INVMASS_NSIGMA"); 

  for(int s=0; s<6; s++){
    invMass[s]->Fit("fitGauss","RQ");
    mu[s]    = fitGauss->GetParameter(1);
    sigma[s] = fitGauss->GetParameter(2);

    muPar   .addValueAndError(mu[s],    fitGauss->GetParError(1)); 
    sigmaPar.addValueAndError(sigma[s], fitGauss->GetParError(2)); 

    cout << "[Pi0Calibration::DoFits] Mu = " << mu[s] << " Sigma = " << sigma[s] << " 1.5 sigma cut = [" << mu[s]-1.5*sigma[s] << "," << mu[s]+1.5*sigma[s] << "]" << endl; 
  }

  nsigmaPar.addValueAndError(1.5, 0.0); 
  pars->addParameterSet(muPar); 
  pars->addParameterSet(sigmaPar); 
  pars->addParameterSet(nsigmaPar); 

}

void Pi0Calibration::Save(string outputFilename){
  TFile *outputFile = new TFile(outputFilename.c_str(), "recreate"); 
  
  for (int s=0; s<6; s++){
    invMass[s]     ->Write(); 
    openingAngle[s]->Write(); 
    energyDiff[s]  ->Write();
    openEnergy[s]  ->Write();

  }
    numberPhotons->Write();   
  outputFile->Write();
  outputFile->Close();

}

int main(int argc, char * argv[]){
  
  // Setup Options
  h22Options opts;
  opts.args["PARS"].args = "/u/home/dmriser/mydoc/analysis/root_scripts/Analysis_v2/lists/data.pars";
  opts.args["PARS"].type = 1;
  opts.args["PARS"].name = "Parameter file";
  opts.set(argc,argv);
  
  long int nev    = opts.args["N"].arg;
  string parfile  = opts.args["PARS"].args;
  
  Parameters pars;
  pars.loadParameters(opts.args["PARS"].args);

  Pi0Calibration Analysis(&opts, &pars); 
  for (std::vector<std::string>::iterator it=opts.ifiles.begin(); it<opts.ifiles.end(); it++) { Analysis.AddFile(*it); }
  Analysis.RunAnalysis(); 
  Analysis.DoFits(); 
  Analysis.Save(opts.args["OUT"].args); 

  pars.saveParameters(opts.args["PARS"].args);  
  
  return 0;
}
