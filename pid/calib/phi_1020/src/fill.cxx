#include <iostream>
#include <map>

// h22 libs 
#include "CommonTools.h"
#include "Corrections.h"
#include "CheckPoints.h"
#include "DataEventCut.h"
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

    // needs parameters 
    params = new Parameters(); 
    //    params->loadParameters(Form("%s/lists/data_tofmass.pars", path.c_str())); 
    params->loadParameters(Form("%s/lists/parameters/data/dataLoose.pars", path.c_str())); 
    filter = new ParticleFilter(params);
    //    filter->positiveKaonSelector->DisableByRegex("Fid"); 

    // ntuple 
    tupleWriter.addInt("helicity");
    tupleWriter.addFloat("im_kk"); 
    tupleWriter.addFloat("im_pk"); 
    tupleWriter.addFloat("missing_mass_ekx_mk"); 
    tupleWriter.addFloat("missing_mass_epkx_mk"); 
    tupleWriter.addFloat("missing_mass_ekx_mpi"); 
    tupleWriter.addFloat("missing_mass_epkx_mpi"); 
    tupleWriter.addFloat("missing_mass_epx");
    tupleWriter.addFloat("x");
    tupleWriter.addFloat("w");
    tupleWriter.addFloat("q2");
    tupleWriter.addFloat("theta_h_p");
    tupleWriter.addFloat("theta_kk_lab");
    tupleWriter.addFloat("p_kp"); 
    tupleWriter.addFloat("p_km"); 
    tupleWriter.addFloat("p_prot"); 
    tupleWriter.addFloat("alpha_kp");
    tupleWriter.addFloat("alpha_prot");
    tupleWriter.addFloat("theta_h_lambda");
    tupleWriter.addFloat("phi_h_lambda"); 
  }

  ~Analysis(){
    // total destruction 
  }

  void Loop(){

    // setup reader options 
    GSIM = false; 
    Init();
    
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
	
	std::vector<int> kpIndices = filter->getVectorOfParticleIndices(event, 321);
	std::vector<int> protonIndices = filter->getVectorOfParticleIndices(event, 2212);

	if (!protonIndices.empty()){

	  // building event 
	  int protonIndex         = protonIndices[0]; 
	  TLorentzVector electron = event.GetTLorentzVector(electronIndex, 11);  
	  TLorentzVector proton   = event.GetTLorentzVector(protonIndex, 2212);
	  
	  // momentum correction done here 
	  // because if we dont find proton it's useless 
	  // to waste cpu doing it above 
	  electron = momCorr->PcorN(electron, -1, 11);

	  // build physics event with only proton 
	  PhysicsEvent ev = builder.getPhysicsEvent(electron, proton); 	  

	  //  now use kp if it exists 
	  if (!kpIndices.empty()){
	    int kaonIndex = kpIndices[0];
	    
	    // stop doing things if its actually a pion 
	    /* 
	    if( IsMisidentifiedPion(electron, kaonIndex) ){
	      continue; 
	    }
	    */
	    
	    // get four-vectors 
	    TLorentzVector kp               = event.GetTLorentzVector(kaonIndex, 321);
	    TLorentzVector kpWithMassOfPion = event.GetTLorentzVector(kaonIndex, 211);
	    
	    PhysicsEvent epkEvent = builder.getPhysicsEvent(electron, proton, kp);
	    PhysicsEvent epkEventWithMassOfPion = builder.getPhysicsEvent(electron, proton, kpWithMassOfPion);
	    

	    // event with proton left out 
	    PhysicsEvent ekEvent =  builder.getPhysicsEvent(electron, kp);
	    PhysicsEvent ekEventWithMassOfPion =  builder.getPhysicsEvent(electron, kpWithMassOfPion);

	    
	    // the final state could be a km 
	    TLorentzVector km = epkEvent.finalState; 
	    TLorentzVector phi = km+kp; 
	    TLorentzVector lambda = km+proton; 

	    // get info on our hadron pid 
	    DataEventCut_BetaPLikelihood *betaPLike_2212 = (DataEventCut_BetaPLikelihood*) filter->protonSelector->GetCut("Beta P Likelihood Cut 2212");
	    DataEventCut_BetaPLikelihood *betaPLike_321  = (DataEventCut_BetaPLikelihood*) filter->positiveKaonSelector->GetCut("Beta P Likelihood Cut 321");


	    //  write event 
	    tupleWriter.setFloat("missing_mass_ekx_mk",   sqrt(ekEvent.mm2));  
	    tupleWriter.setFloat("missing_mass_epkx_mk",  sqrt(epkEvent.mm2));  
	    tupleWriter.setFloat("missing_mass_ekx_mpi",  sqrt(ekEventWithMassOfPion.mm2));  
	    tupleWriter.setFloat("missing_mass_epkx_mpi", sqrt(epkEventWithMassOfPion.mm2));  
	    tupleWriter.setFloat("missing_mass_epx", sqrt(ev.mm2));  
	    tupleWriter.setFloat("im_kk", sqrt(phi.M2()));  
	    tupleWriter.setFloat("im_pk", sqrt(lambda.M2()));  
	    

	    tupleWriter.setFloat("p_prot", proton.P()); 
	    tupleWriter.setFloat("p_kp", kp.P()); 
	    tupleWriter.setFloat("p_km", km.P()); 

	    tupleWriter.setFloat("x", ev.x);  
	    tupleWriter.setFloat("w", ev.w);  
	    tupleWriter.setFloat("q2", ev.qq);  
	    tupleWriter.setFloat("theta_h_p", ev.thetaHadron);  
	    tupleWriter.setFloat("theta_kk_lab", to_degrees*kp.Angle(km.Vect()));  
	    

	    tupleWriter.setFloat("theta_h_lambda", ekEvent.thetaHadron); 
	    tupleWriter.setFloat("phi_h_lambda", ekEvent.phiHadron); 

	    // these have to be called again because they are setup in a bad way that 
	    // the local variable (protected) inside is set fConfidnce based on the last 
	    // call. 
	    
	    if(betaPLike_321->IsPassed(event,kpIndices[0])){
	      tupleWriter.setFloat("alpha_kp", betaPLike_321 ->GetConfidence());
	    }
	    if(betaPLike_2212->IsPassed(event,protonIndices[0])){
	      tupleWriter.setFloat("alpha_prot", betaPLike_2212->GetConfidence());
	    }
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
