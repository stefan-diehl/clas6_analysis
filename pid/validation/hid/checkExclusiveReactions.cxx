////////////////////////////////////////
/*
 David Riser, University of Connecticut
 
 July 13, 2016
 
 Modified:
 March 8, 2017
 April 14, 2017 

 meson.cxx -> fill histograms for mesons 
 
 */
////////////////////////////////////////

#include <iostream>

#include "CommonTools.h"
#include "Corrections.h"
#include "DataEventCut.h"
#include "h22Event.h"
#include "h22Option.h"
#include "h22Reader.h"
#include "GenericAnalysis.h"
#include "MesonHistograms.h"
#include "MomCorr.h"
#include "Parameters.h"
#include "ParameterSet.h"
#include "Pars.h"
#include "ParticleFilter.h"
#include "PhysicsEvent.h"
#include "PhysicsEventBuilder.h"
#include "StandardHistograms.h"

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"

class HIDCalibration : public GenericAnalysis {

public:
    HIDCalibration(h22Options *opts, Parameters *params) : GenericAnalysis(opts), pars(params) { }
    ~HIDCalibration(){ }

public:
  Parameters          *pars;
  ParticleFilter      *filter;
  Corrections          corr;
  PhysicsEventBuilder *builder;   
  MomCorr_e1f         *momCorr; 

  MesonHistograms *kpHistos; 
  MesonHistograms *kmHistos; 
  StandardHistograms *kpStand; 
  StandardHistograms *kmStand; 

  TH1D *h1_mm2_withProton[7];

  void Initialize();
  void ProcessEvent();
  void Save(string outputFilename);
  void InitHistos();
  bool CurrentParticleIsNotElectronCandidate(std::vector<int> &electronCandidates, int index);

};

void HIDCalibration::InitHistos() {
  kpHistos = new MesonHistograms("kp",  321);
  kmHistos = new MesonHistograms("km",  -321);
  kpStand  = new StandardHistograms("kp", 0); 
  kmStand  = new StandardHistograms("km", 0); 

  for (int i=0; i<7; ++i){
    h1_mm2_withProton[i] = new TH1D(Form("h1_mm2_withProton_sect%d",i), Form("h1_mm2_withProton_sect%d",i), 40, -1.5, 1.5); 
  }
}

void HIDCalibration::Initialize(){
  InitHistos();
  
  filter = new ParticleFilter(pars);
  filter->set_info(GSIM, GetRunNumber());

  builder = new PhysicsEventBuilder(); 

  // not really safe for use on the farm because I don't set vars 
  // so not sure what to do about this.  maybe passing in the path
  // as a parameter would be ok? 
  std::string path = Global::Environment::GetAnalysisPath(); 
  std::string momCorrPath = Form("%s/momCorr/",path.c_str());
  momCorr = new MomCorr_e1f("/u/home/dmriser/Analysis_v2/momCorr/");
}

void HIDCalibration::ProcessEvent(){

  vector<int> electronCandidates = filter->getVectorOfParticleIndices(event, 11);
  if ( !electronCandidates.empty() ){
    int electronIndex = electronCandidates[0];
    event.SetElectronIndex(electronIndex); 
    corr.correctEvent(&event, GetRunNumber(), GSIM); 
    
    TLorentzVector electron = event.GetTLorentzVector(electronIndex, 11); 
    electron = momCorr->PcorN(electron, -1, 11);
    
    std::vector<int> kpIndices = filter->getVectorOfParticleIndices(event, 321); 
    std::vector<int> kmIndices = filter->getVectorOfParticleIndices(event, -321); 
    std::vector<int> protIndices = filter->getVectorOfParticleIndices(event, -211);
 
    PhysicsEvent candidateEvent = builder->getPhysicsEvent(electron); 

      // ------------------------------------------------------
      //        fill tracks that pass everything 
      // ------------------------------------------------------

    if (candidateEvent.qq > 1.0 && candidateEvent.w > 2.00) {

      if (!kpIndices.empty()){
	TLorentzVector meson = event.GetTLorentzVector(kpIndices[0], 321);
	meson = momCorr->PcorN(meson, 1, 321);
	
	PhysicsEvent physicsEvent = builder->getPhysicsEvent(electron, meson); 
	if (1){//(physicsEvent.mm2 > 1.0){
	  kpHistos->Fill(event, physicsEvent, kpIndices[0]);
	  kpStand->Fill(event, electronIndex, kpIndices[0], physicsEvent); 
	}
 
	if (!protIndices.empty()){
	  TLorentzVector proton     = event.GetTLorentzVector(protIndices[0], -211); 
	  PhysicsEvent protonEvent  = builder->getPhysicsEvent(electron, proton, meson); 
	  int sect                  = event.dc_sect[electronIndex]; 
	  
	  h1_mm2_withProton[0]   ->Fill(protonEvent.mm2); 
	  h1_mm2_withProton[sect]->Fill(protonEvent.mm2); 
	}
	
     }

      if (!kmIndices.empty()){
	TLorentzVector meson = event.GetTLorentzVector(kmIndices[0], -321);
	meson = momCorr->PcorN(meson, -1, -321);
	
	PhysicsEvent physicsEvent = builder->getPhysicsEvent(electron, meson); 
	kmHistos->Fill(event, physicsEvent, kmIndices[0]);
	kmStand->Fill(event, electronIndex, kmIndices[0], physicsEvent); 
      }
    }
  }
}

bool HIDCalibration::CurrentParticleIsNotElectronCandidate(std::vector<int> &electronCandidates,int index){
    return (std::find(electronCandidates.begin(), electronCandidates.end(), index) == electronCandidates.end());
}

void HIDCalibration::Save(string outputFilename){
    TFile *outputFile = new TFile(outputFilename.c_str(), "recreate");

    kpHistos->Save(outputFile);
    kmHistos->Save(outputFile);
    kpStand->Save(outputFile); 
    kmStand->Save(outputFile); 

    for (int i=0; i<7; ++i){
      h1_mm2_withProton[i]->Write(); 
    }

    outputFile->Write();
    outputFile->Close();
}

int main(int argc, char * argv[]){

    // Setup Options
    h22Options opts;
    opts.args["PARS"].args = "/u/home/dmriser/Analysis_v2/lists/data.pars";
    opts.args["PARS"].type = 1;
    opts.args["PARS"].name = "Parameter file";
    opts.args["LIST"].args = "UNSET";
    opts.args["LIST"].type = 1;
    opts.args["LIST"].name = "List of files";
    opts.set(argc,argv);

    int GSIM        = opts.args["MC"].arg;
    long int nev    = opts.args["N"].arg;
    string parfile  = opts.args["PARS"].args;

    Parameters pars;
    pars.loadParameters(opts.args["PARS"].args);

    std::vector<std::string> files; 
    if (opts.args["LIST"].args != "UNSET"){
      files = loadFilesFromList(opts.args["LIST"].args, 10000);
    } else {
      files = loadFilesFromCommandLine(&opts, 10000);
    }

    HIDCalibration Analysis(&opts, &pars);
    for (std::vector<std::string>::iterator it=files.begin(); it<files.end(); it++) { Analysis.AddFile(*it); }
    Analysis.RunAnalysis();
    Analysis.Save(opts.args["OUT"].args);

    Analysis.filter->GetSelector(321)->PrintSummary(); 
    Analysis.filter->GetSelector(-321)->PrintSummary(); 

    return 0;
}
