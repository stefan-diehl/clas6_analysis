#include <iostream>

// h22 libs 
#include "CommonTools.h"
#include "Corrections.h"
#include "CheckPoints.h"
#include "DBins.h"
#include "h22Event.h"
#include "h22Reader.h"
#include "HadronID.h"
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
    params->loadParameters(Form("%s/lists/parameters/data/final.pars", path.c_str())); 

    // particle filter for doing pid
    filter      = new ParticleFilter(params);

    // setup structure of ntuple 
    tupleWriter.addInt("helicity");
    tupleWriter.addFloat("missing_mass"); 
    tupleWriter.addFloat("x"); 
    tupleWriter.addFloat("q2"); 
    tupleWriter.addFloat("z"); 
    tupleWriter.addFloat("pt"); 
    tupleWriter.addFloat("w"); 
    tupleWriter.addFloat("eta"); 
    tupleWriter.addFloat("phi_h"); 
    tupleWriter.addFloat("theta_h"); 
    tupleWriter.addFloat("p_ele"); 
    tupleWriter.addFloat("p_pip"); 
    tupleWriter.addFloat("phi_ele"); 
    tupleWriter.addFloat("phi_pip");
    tupleWriter.addFloat("theta_ele"); 
    tupleWriter.addFloat("theta_pip");
    tupleWriter.addFloat("dvz");
    tupleWriter.addFloat("dist_ecsf");
    tupleWriter.addFloat("dist_ec_edep");
    tupleWriter.addFloat("dist_vz");
    tupleWriter.addFloat("dist_cc_theta");
  }

  ~Analysis(){
  }

  void Loop(){

    // setup reader options 
    GSIM = false; 
    Init();
    
    // setup particle filter 
    filter->set_info(GSIM, GetRunNumber());

    StatusBar  stat; 
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
	  
	std::vector<int> pipIndices = filter->getVectorOfParticleIndices(event, 211); 

	if(!pipIndices.empty()){

	  int pipIndex       = pipIndices[0];
	  TLorentzVector pip = event.GetTLorentzVector(pipIndex, 211); 

	  // momentum correction done here 
	  // because if we dont find proton it's useless 
	  // to waste cpu doing it above 
	  TLorentzVector electron = event.GetTLorentzVector(electronIndex, 11); 
	  electron = momCorr->PcorN(electron, -1, 11);
	  
	  // build event 
	  PhysicsEvent ev = builder.getPhysicsEvent(electron, pip);

	  if ( sqrt(ev.mm2) > 0.8 && sqrt(ev.mm2) < 1.08 
	       && ev.w > 2.0 && ev.qq > 1.0) {

	    tupleWriter.setInt("helicity",       event.corr_hel);
	    tupleWriter.setFloat("missing_mass", sqrt(ev.mm2));
	    tupleWriter.setFloat("x",            ev.x);
	    tupleWriter.setFloat("q2",           ev.qq);
	    tupleWriter.setFloat("z",            ev.z);
	    tupleWriter.setFloat("pt",           ev.pT);
	    tupleWriter.setFloat("w",            ev.w);
	    tupleWriter.setFloat("eta",          ev.eta);
	    tupleWriter.setFloat("theta_h",      ev.thetaHadron);
	    tupleWriter.setFloat("phi_h",        ev.phiHadron);
	    tupleWriter.setFloat("p_ele",        electron.P()); 
	    tupleWriter.setFloat("p_pip",        pip.P()); 
	    tupleWriter.setFloat("theta_ele",    to_degrees*electron.Theta()); 
	    tupleWriter.setFloat("theta_pip",    to_degrees*pip.Theta()); 
	    tupleWriter.setFloat("phi_ele",      to_degrees*electron.Phi()); 
	    tupleWriter.setFloat("phi_pip",      to_degrees*pip.Phi()); 
	    tupleWriter.setFloat("dvz",          event.vz[electronIndex]-event.vz[pipIndex]); 
	    tupleWriter.writeEvent();
	  }
	}
      }
      
      if (ientry%10000 == 0){
	stat.PrintStatus(ientry, GetEntries()); 
      }
    }

    // summarize results 
    std::cout << std::endl;
    double loopTime  = timer.RealTime(); 
    double eventRate = GetEntries()/loopTime; 
    std::cout << "[GenericAnalysis::Loop] Total time: " 
	      << loopTime << " seconds, Event rate: " 
	      << eventRate << " events/sec. " << std::endl;

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
  SimpleNTupleWriter       tupleWriter; 
};

int main(int argc, char *argv[]){

  if (argc > 1){
  // setup analysis 
  Analysis analysis; 

  // Setup Options
  h22Options opts;
  opts.args["PARS"].args = "/u/home/dmriser/analysis-main/lists/data.pars";
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
