#include <iostream>

// h22 libs 
#include "CommonTools.h"
#include "Corrections.h"
#include "CheckPoints.h"
#include "DBins.h"
#include "h22Event.h"
#include "h22Reader.h"
#include "MesonHistograms.h"
#include "MomCorr.h"
#include "ParticleFilter.h"
#include "Parameters.h"
#include "PhysicsEvent.h"
#include "PhysicsEventBuilder.h"
#include "SimpleNTupleWriter.h"
#include "StatusBar.h"
#include "StandardHistograms.h"

// root 
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "TStopwatch.h"
#include "TVector3.h"

class Analysis : public h22Reader {
public:  
  Analysis() {

    std::string path        = Global::Environment::GetAnalysisPath(); 
    std::string momCorrPath = Form("%s/momCorr/",path.c_str());
    momCorr                 = new MomCorr_e1f(momCorrPath);

    // setup reader options 
    GSIM = false; 
    Init();

    // needs parameters 
    params = new Parameters(); 
    params->loadParameters(Form("%s/lists/data_tofmass.pars", path.c_str())); 
    filter = new ParticleFilter(params);

    // ntuple 
    tupleWriter.addInt("helicity");
    tupleWriter.addInt("topology");
    tupleWriter.addTLorentzVector("electron"); 
    tupleWriter.addTLorentzVector("proton"); 
    tupleWriter.addTLorentzVector("kaon_pos"); 
    tupleWriter.addTLorentzVector("kaon_neg"); 
  }

   ~Analysis(){
    // total destruction 
  }

  void Loop(){
    
    // setup particle filter 
    filter->set_info(GSIM, GetRunNumber());

    StatusBar stat; 

    TStopwatch timer; 
    timer.Reset();
    timer.Start();

    // for every event
    for(int ientry=0; ientry<GetEntries(); ientry++){
      GetEntry(ientry);
      
      std::vector<int> electronIndices = filter->getVectorOfParticleIndices(event, 11);
      if(!electronIndices.empty()){
	
	int electronIndex = electronIndices[0];
	event.SetElectronIndex(electronIndex);
	Corrections::correctEvent(&event, GetRunNumber(), GSIM);
	
	// search mofo 
	std::vector<int> kpIndices     = filter->getVectorOfParticleIndices(event, 321);
	std::vector<int> kmIndices     = filter->getVectorOfParticleIndices(event, -321);
	std::vector<int> protonIndices = filter->getVectorOfParticleIndices(event, 2212);

	// something will be written 
	if ( (kpIndices.size() > 0)&&(kmIndices.size() > 0) || (kpIndices.size() > 0) && (protonIndices.size() > 0) || (protonIndices.size() > 0) && (kmIndices.size() > 0) ) {
	  int top = -1; 

	  TLorentzVector electron = event.GetTLorentzVector(electronIndex, 11); 
	  electron = momCorr->PcorN(electron, -1, 11);
	  
	  tupleWriter.setInt("helicity", event.corr_hel); 

	  if ( (kpIndices.size() > 0)&&(kmIndices.size() > 0) ){
	    top = 2;
	    tupleWriter.setInt("topology", top); 
	    
	    TLorentzVector kp = event.GetTLorentzVector(kpIndices[0], 321); 
	    TLorentzVector km = event.GetTLorentzVector(kmIndices[0], -321); 
	    
	    PhysicsEvent physicsEvent = builder.getPhysicsEvent(electron, kp, km); 
	    TLorentzVector proton = physicsEvent.finalState; 
	    tupleWriter.setTLorentzVector("electron", electron); 
	    tupleWriter.setTLorentzVector("proton", proton); 
	    tupleWriter.setTLorentzVector("kaon_pos", kp); 
	    tupleWriter.setTLorentzVector("kaon_neg", km); 
	    
	    tupleWriter.writeEvent(); 
	}
	  
	  if ( (kmIndices.size() > 0)&&(protonIndices.size() > 0) ){
	    top = 1;
	    tupleWriter.setInt("topology", top); 
	    
	    TLorentzVector proton = event.GetTLorentzVector(protonIndices[0], 2212); 
	    TLorentzVector km = event.GetTLorentzVector(kmIndices[0], -321); 
	    
	    PhysicsEvent physicsEvent = builder.getPhysicsEvent(electron, proton, km); 
	    TLorentzVector kp = physicsEvent.finalState; 
	    tupleWriter.setTLorentzVector("electron", electron); 
	    tupleWriter.setTLorentzVector("proton", proton); 
	    tupleWriter.setTLorentzVector("kaon_pos", kp); 
	    tupleWriter.setTLorentzVector("kaon_neg", km); 
	    
	    tupleWriter.writeEvent(); 
	}
	  
	  if ( (kpIndices.size() > 0)&&(protonIndices.size() > 0) ){
	    top = 0;
	    tupleWriter.setInt("topology", top); 
	    
	    TLorentzVector kp = event.GetTLorentzVector(kpIndices[0], 321); 
	    TLorentzVector proton = event.GetTLorentzVector(protonIndices[0], 2212); 
	    
	    PhysicsEvent physicsEvent = builder.getPhysicsEvent(electron, kp, proton); 
	    TLorentzVector km = physicsEvent.finalState; 
	    tupleWriter.setTLorentzVector("electron", electron); 
	    tupleWriter.setTLorentzVector("proton", proton); 
	    tupleWriter.setTLorentzVector("kaon_pos", kp); 
	    tupleWriter.setTLorentzVector("kaon_neg", km); 
	    
	    tupleWriter.writeEvent(); 
	}
	  
	  if ( (kpIndices.size() > 0)&&(kmIndices.size() > 0)&&(protonIndices.size() > 0) ){
	    top = 3;
	    tupleWriter.setInt("topology", top); 
	    
	    TLorentzVector kp     = event.GetTLorentzVector(kpIndices[0], 321); 
	    TLorentzVector km     = event.GetTLorentzVector(kmIndices[0], -321); 
	    TLorentzVector proton = event.GetTLorentzVector(protonIndices[0], 2212); 

	    tupleWriter.setTLorentzVector("electron", electron); 
	    tupleWriter.setTLorentzVector("proton", proton); 
	    tupleWriter.setTLorentzVector("kaon_pos", kp); 
	    tupleWriter.setTLorentzVector("kaon_neg", km); 
	    tupleWriter.writeEvent(); 
	}
	  
	  

	}	    
      }

      if (ientry%10000 == 0){
	stat.PrintStatus(ientry, GetEntries()); 
      }
    }

    std::cout << std::endl;
    double loopTime  = timer.RealTime(); 
    double eventRate = GetEntries()/loopTime; 
    std::cout << "[GenericAnalysis::Loop] Total time: " 
	      << loopTime << " seconds, Event rate: " 
	      << eventRate << " events/sec. " << std::endl;

  }

  bool IsMisidentifiedPion(TLorentzVector & electron, int mesonIndex){

    // Dan Carman method for rejecting some pi+ events in the 
    // kaon sample. 
    //
    // Classify an event and if it's a kaon you pretend 
    // it's a pion.  Then check the missing mass spectrum 
    // and if the event is near the neutron peak it must 
    // really be a pion, so reassign it to be.  
    //
    // Currently all it does it throw out the pion, 
    // if doesn't actually add it to the pion sample 
    // because the code never makes it that far in the 
    // pion case it would just return 0 pions and not enter.
    
    TLorentzVector pion = event.GetTLorentzVector(mesonIndex, 211); 
    PhysicsEvent pionEvent = builder.getPhysicsEvent(electron,pion); 
    
    // reject and dont process event 
    if (sqrt(pionEvent.mm2) < 1.07){
      return true; 
    }
    
    return false; 
  }

  void Save(std::string outFile){
    TFile *out = new TFile(outFile.c_str(), "recreate");

    tupleWriter.save(out); 

    out->Close();
  }


protected:
  MomCorr_e1f             *momCorr; 
  PhysicsEventBuilder      builder; 
  ParticleFilter          *filter; 
  Parameters              *params; 

  SimpleNTupleWriter      tupleWriter; 

};

int main(int argc, char *argv[]){

  if (argc > 1){
    // setup analysis
    Analysis analysis;

    // Setup Options
    h22Options opts;
    opts.args["PARS"].args = "/u/home/dmriser/Analysis_v2/lists/data.pars";
    opts.args["PARS"].type = 1;
    opts.args["PARS"].name = "Parameter file";
    opts.args["LIST"].args = "UNSET";
    opts.args["LIST"].type = 1;
    opts.args["LIST"].name = "List of files";
    opts.set(argc,argv);

    std::vector<std::string> files;
    if (opts.args["LIST"].args != "UNSET"){
      files = loadFilesFromList(opts.args["LIST"].args, 10000);
    } else {
      files = loadFilesFromCommandLine(&opts, 10000);
    }

    for (std::string f : files){
      analysis.AddFile(f);
    }

    // run analysis loop
    analysis.Loop();
    analysis.Save(opts.args["OUT"].args);

  } else {
    std::cerr << "No files found." << std::endl;
  }

  return 0;
}
