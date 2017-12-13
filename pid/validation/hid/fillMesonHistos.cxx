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

  MesonHistograms *posHistos; 
  MesonHistograms *negHistos; 

  MesonHistograms *posCutDVzHistos; 
  MesonHistograms *negCutDVzHistos; 

  MesonHistograms *posCutDCHistos; 

  MesonHistograms *posCutDISHistos; 
  MesonHistograms *negCutDISHistos; 

  MesonHistograms *posCutDISMMHistos; 
  MesonHistograms *negCutDISMMHistos; 

  MesonHistograms *kpHistos; 
  MesonHistograms *kmHistos; 

  StandardHistograms *posStand; 
  StandardHistograms *negStand; 
  StandardHistograms *kpStand; 
  StandardHistograms *kmStand; 

  void Initialize();
  void ProcessEvent();
  void Save(string outputFilename);
  void InitHistos();
  bool CurrentParticleIsNotElectronCandidate(std::vector<int> &electronCandidates, int index);

};

void HIDCalibration::InitHistos() {
  posHistos = new MesonHistograms("pos",  321);
  negHistos = new MesonHistograms("neg", -321);
  kpHistos = new MesonHistograms("kp",    321);
  kmHistos = new MesonHistograms("km",   -321);

  posCutDVzHistos = new MesonHistograms("posCutDVz",    321);
  negCutDVzHistos = new MesonHistograms("negCutDVz",   -321);

  posCutDISHistos = new MesonHistograms("posCutDIS",    321);
  negCutDISHistos = new MesonHistograms("negCutDIS",   -321);

  posCutDISMMHistos = new MesonHistograms("posCutDISMM",    321);
  negCutDISMMHistos = new MesonHistograms("negCutDISMM",   -321);

  posCutDCHistos = new MesonHistograms("posCutDC",    321);

  posStand = new StandardHistograms("pos", 0); 
  negStand = new StandardHistograms("neg", 0); 
  kpStand = new StandardHistograms("kp", 0); 
  kmStand = new StandardHistograms("km", 0); 
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

    PhysicsEvent candidateEvent = builder->getPhysicsEvent(electron); 

      // ------------------------------------------------------
      //        fill tracks that pass everything 
      // ------------------------------------------------------

    if (candidateEvent.qq > 1.0 && candidateEvent.w > 2.00) {

      if (!kpIndices.empty()){
	TLorentzVector meson = event.GetTLorentzVector(kpIndices[0], 321);
	meson = momCorr->PcorN(meson, 1, 321);
	
	PhysicsEvent physicsEvent = builder->getPhysicsEvent(electron, meson); 
	if (event.p[kpIndices[0]] <= 2.3){
	  kpHistos->Fill(event, physicsEvent, kpIndices[0]);
	  kpStand->Fill(event, electronIndex, kpIndices[0], physicsEvent); 
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

      // ------------------------------------------------------
      //        fill each cut 
      // ------------------------------------------------------
            
      for(int ipart=0; ipart<event.gpart; ipart++){
	double dvz = fabs(event.corr_vz[electronIndex]-event.corr_vz[ipart]); 
	
	// k- candidate 
	if (event.q[ipart] < 0 && CurrentParticleIsNotElectronCandidate(electronCandidates, ipart)){
	  TLorentzVector part1 = event.GetTLorentzVector(ipart, event.id[ipart]);
	  part1 = momCorr->PcorN(part1, event.q[ipart], event.id[ipart]);
	  PhysicsEvent physicsEvent = builder->getPhysicsEvent(electron, part1); 
 
	  negHistos->Fill(event, physicsEvent, ipart); 

	  if (dvz < 4.0){
	    negCutDVzHistos->Fill(event, physicsEvent, ipart); 
	  }

	  if (physicsEvent.qq > 1.0 && physicsEvent.w){ 
	    negCutDISHistos->Fill(event, physicsEvent, ipart);
	  }

	  if (physicsEvent.qq > 1.0 && physicsEvent.w > 2.0 && physicsEvent.mm2 > 1.0){ 
	    negCutDISMMHistos->Fill(event, physicsEvent, ipart);
	  }
	} 
	
	// k+ candidate 
	else if (event.q[ipart] > 0){
	  TLorentzVector part1 = event.GetTLorentzVector(ipart, event.id[ipart]);
	  part1 = momCorr->PcorN(part1, event.q[ipart], event.id[ipart]);
	  PhysicsEvent physicsEvent = builder->getPhysicsEvent(electron, part1); 
 
	  posHistos->Fill(event, physicsEvent, ipart); 

	  if (dvz < 4.0){
	    posCutDVzHistos->Fill(event, physicsEvent, ipart); 
	  }

	  if (filter->GetSelector(321)->GetCut("DC Region 1 Fid Cut")->IsPassed(event, ipart)){
	    posCutDCHistos->Fill(event, physicsEvent, ipart); 
	  }

	  if (physicsEvent.qq > 1.0 && physicsEvent.w){ 
	    posCutDISHistos->Fill(event, physicsEvent, ipart);
	  }

	  if (physicsEvent.qq > 1.0 && physicsEvent.w > 2.0 && physicsEvent.mm2 > 1.0){ 
	    posCutDISMMHistos->Fill(event, physicsEvent, ipart);
	  }


	}
      }
  }
}

bool HIDCalibration::CurrentParticleIsNotElectronCandidate(std::vector<int> &electronCandidates,int index){
    return (std::find(electronCandidates.begin(), electronCandidates.end(), index) == electronCandidates.end());
}

void HIDCalibration::Save(string outputFilename){
    TFile *outputFile = new TFile(outputFilename.c_str(), "recreate");

    posHistos->Save(outputFile);
    negHistos->Save(outputFile);

    posCutDVzHistos->Save(outputFile); 
    negCutDVzHistos->Save(outputFile); 

    posCutDCHistos->Save(outputFile); 

    posCutDISHistos->Save(outputFile); 
    negCutDISHistos->Save(outputFile); 

    posCutDISMMHistos->Save(outputFile); 
    negCutDISMMHistos->Save(outputFile); 

    kpHistos->Save(outputFile);
    kmHistos->Save(outputFile);

    posStand->Save(outputFile); 
    negStand->Save(outputFile); 
    kpStand->Save(outputFile); 
    kmStand->Save(outputFile); 

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

