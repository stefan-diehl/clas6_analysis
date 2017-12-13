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


    /////////////////////////////////////////////////////
    // setup entries of the output tree

    tupleWriter.addInt("helicity");
    tupleWriter.addFloat("missing_mass2"); 
    tupleWriter.addFloat("x"); 
    tupleWriter.addFloat("q2"); 
    tupleWriter.addFloat("z"); 
    tupleWriter.addFloat("pt"); 
    tupleWriter.addFloat("w"); 
    tupleWriter.addFloat("eta"); 
    tupleWriter.addFloat("phi_h"); 
    tupleWriter.addFloat("theta_h"); 
    tupleWriter.addFloat("E_ele");
    tupleWriter.addFloat("px_ele"); 
    tupleWriter.addFloat("py_ele"); 
    tupleWriter.addFloat("pz_ele");  
    tupleWriter.addFloat("E_prot"); 
    tupleWriter.addFloat("px_prot"); 
    tupleWriter.addFloat("py_prot"); 
    tupleWriter.addFloat("pz_prot"); 
    tupleWriter.addFloat("p_ele"); 
    tupleWriter.addFloat("p_prot");
    tupleWriter.addFloat("phi_ele"); 
    tupleWriter.addFloat("phi_prot");
    tupleWriter.addFloat("theta_ele"); 
    tupleWriter.addFloat("theta_prot");
    tupleWriter.addFloat("dvz");

   //////////////////////////////////////////////////////

  }

////////////////////////////////////////////////////////
///  destructor

  ~Analysis(){
  }


//////////////////////////////////////////////////////////////////////////////////////////////
///  main event loop
//////////////////////////////////////////////////////////////////////////////////////////////

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

      std::vector<int> electronIndices = filter->getVectorOfParticleIndices(event, 11);   // get electrons from the event

      if(!electronIndices.empty()){    // if event has an electron
	
	int electronIndex = electronIndices[0];
	event.SetElectronIndex(electronIndex);

	Corrections::correctEvent(&event, GetRunNumber(), GSIM);
	  
	std::vector<int> pipIndices = filter->getVectorOfParticleIndices(event, 2212);    // get proton from the event

	if(!pipIndices.empty()){    // if event has proton (and electron)

	  int protIndex       = pipIndices[0];
	  TLorentzVector prot = event.GetTLorentzVector(protIndex, 2212); 

	  // momentum correction done here 
	  // because if we dont find proton it's useless 
	  // to waste cpu doing it above 

	  TLorentzVector electron = event.GetTLorentzVector(electronIndex, 11);  
	  electron = momCorr->PcorN(electron, -1, 11);    // corrections for the electron
	  prot = momCorr->PcorN(prot, 1, 2212);           // corrections for the electron

	  // build event 

	  PhysicsEvent ev = builder.getPhysicsEvent(electron, prot);

	  //if(ev.mm2 > -0.4 && ev.mm2 < 1) {

	    tupleWriter.setInt("helicity",        event.corr_hel);
	    tupleWriter.setFloat("missing_mass2", ev.mm2);
	    tupleWriter.setFloat("x",             ev.x);
	    tupleWriter.setFloat("q2",            ev.qq);
	    tupleWriter.setFloat("z",             ev.z);
	    tupleWriter.setFloat("pt",            ev.pT);
	    tupleWriter.setFloat("w",             ev.w);
	    tupleWriter.setFloat("eta",           ev.eta);
	    tupleWriter.setFloat("theta_h",       ev.thetaHadron);
	    tupleWriter.setFloat("phi_h",         ev.phiHadron);
	    tupleWriter.setFloat("E_ele",         electron.E());
	    tupleWriter.setFloat("px_ele",        electron.Px());
	    tupleWriter.setFloat("py_ele",        electron.Py());
            tupleWriter.setFloat("pz_ele",        electron.Pz());
	    tupleWriter.setFloat("E_prot",        prot.E()); 
            tupleWriter.setFloat("px_prot",       prot.Px()); 
            tupleWriter.setFloat("py_prot",       prot.Py()); 
            tupleWriter.setFloat("pz_prot",       prot.Pz()); 
            tupleWriter.setFloat("p_ele",    	  electron.P()); 
	    tupleWriter.setFloat("p_prot",    	  prot.P()); 
	    tupleWriter.setFloat("theta_ele",     to_degrees*electron.Theta()); 
	    tupleWriter.setFloat("theta_prot",    to_degrees*prot.Theta()); 
	    tupleWriter.setFloat("phi_ele",       to_degrees*electron.Phi()); 
	    tupleWriter.setFloat("phi_prot",      to_degrees*prot.Phi()); 
	    tupleWriter.setFloat("dvz",           event.vz[electronIndex]-event.vz[protIndex]); 
	    tupleWriter.writeEvent();

	  //}
	}
      }
      
      if (ientry%10000 == 0){
	stat.PrintStatus(ientry, GetEntries()); 
      }
    }

    ////////////////////////////////////////////////////////////////
    // summarize results 
    std::cout << std::endl;
    double loopTime  = timer.RealTime(); 
    double eventRate = GetEntries()/loopTime; 
    std::cout << "[GenericAnalysis::Loop] Total time: " 
	      << loopTime << " seconds, Event rate: " 
	      << eventRate << " events/sec. " << std::endl;

    ///////////////////////////////////////////////////////////////

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
