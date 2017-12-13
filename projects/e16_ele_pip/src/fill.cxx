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
    tupleWriter.addFloat("missing_mass"); 
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
    tupleWriter.addFloat("E_pip"); 
    tupleWriter.addFloat("px_pip"); 
    tupleWriter.addFloat("py_pip"); 
    tupleWriter.addFloat("pz_pip"); 
    tupleWriter.addFloat("p_ele"); 
    tupleWriter.addFloat("p_pip");
    tupleWriter.addFloat("phi_ele"); 
    tupleWriter.addFloat("phi_pip");
    tupleWriter.addFloat("theta_ele"); 
    tupleWriter.addFloat("theta_pip");
    tupleWriter.addFloat("dvz");


    // Define standard histograms for different particle types

    negatives = new StandardHistograms("negatives", 0); 
    electrons = new StandardHistograms("electrons", 0);
    //pim       = new StandardHistograms("pim", 0);
    //kminus    = new StandardHistograms("kam", 0);

    positives = new StandardHistograms("positives", 0); 
    protons   = new StandardHistograms("protons", 0);
    pip       = new StandardHistograms("pip", 0);
    //kplus     = new StandardHistograms("kap", 0);

    neutrals = new StandardHistograms("neutrals", 0);
    gammas = new StandardHistograms("gammas", 0);


    // disable specific cuts

    filter->electronSelector->DisableByName("Z-Vertex Cut");   // kill z vertex cut (different target position)
    filter->electronSelector->DisableByRegex("CCT");  // kill all cherenkov counter based cuts (contain CCT in their name)


    // set beam and target properties of the run

    TLorentzVector beam(0,0,5.754,5.754);
    TLorentzVector target(0,0,0,proton_mass);
    builder = PhysicsEventBuilder(beam, target);

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
 

 /////////////////////////////////////////////////////////
    /// Get entry (event) from tree and set helicity

      GetEntry(ientry);
      //SetHel();

    /////////////////////////////////////////////////////////
    /// fill standrad histograms
  /*
    for(int ipart=0; ipart<event.gpart; ipart++){
	if(event.q[ipart] < 0){ 
		negatives->Fill(event, ipart);
	}
        if(event.q[ipart] > 0){ 
		positives->Fill(event, ipart);
	}
        if(event.q[ipart] == 0){ 
		neutrals->Fill(event, ipart);
	}
    }
*/
    //////////////////////////////////////////////////////////


      std::vector<int> electronIndices = filter->getVectorOfParticleIndices(event, 11);

      if(!electronIndices.empty()){

        //electrons->Fill(event, electronIndices.front());
	
	int electronIndex = electronIndices[0];
	event.SetElectronIndex(electronIndex);

	// Corrections::correctEvent(&event, GetRunNumber(), GSIM);    //  --> TO be modified for e16
	  
	std::vector<int> pipIndices = filter->getVectorOfParticleIndices(event, 211); 

	if(!pipIndices.empty()){

          //pip->Fill(event, pipIndices.front());

	  int pipIndex       = pipIndices[0];
	  TLorentzVector pip = event.GetTLorentzVector(pipIndex, 211);

	  // momentum correction done here 
	  // because if we dont find proton it's useless 
	  // to waste cpu doing it above 
	  
          TLorentzVector electron = event.GetTLorentzVector(electronIndex, 11); 

	  //electron = momCorr->PcorN(electron, -1, 11);
	  

	  // build event
 
	  PhysicsEvent ev = builder.getPhysicsEvent(electron, pip);

	  //if(ev.w > 2.0 && ev.qq > 1.0) {

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
	    tupleWriter.setFloat("E_pip",        pip.E()); 
            tupleWriter.setFloat("px_pip",       pip.Px()); 
            tupleWriter.setFloat("py_pip",       pip.Py()); 
            tupleWriter.setFloat("pz_pip",       pip.Pz()); 
            tupleWriter.setFloat("p_ele",    	  electron.P()); 
	    tupleWriter.setFloat("p_pip",    	  pip.P()); 
	    tupleWriter.setFloat("theta_ele",     to_degrees*electron.Theta()); 
	    tupleWriter.setFloat("theta_pip",    to_degrees*pip.Theta()); 
	    tupleWriter.setFloat("phi_ele",       to_degrees*electron.Phi()); 
	    tupleWriter.setFloat("phi_pip",      to_degrees*pip.Phi()); 
	    tupleWriter.setFloat("dvz",           event.vz[electronIndex]-event.vz[pipIndex]); 

	    tupleWriter.writeEvent();

	 // }

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

    std::cout << "\nElectron cut summary: \n";
    filter->electronSelector->PrintSummary();    // infos on which cuts are passed and which failed

    //std::cout << "\nPi plus cut summary: \n";
    //filter->pipSelector->PrintSummary();    // infos on which cuts are passed and which failed


  }

  void Save(std::string outFile){
    TFile *out = new TFile(outFile.c_str(), "recreate");
    tupleWriter.save(out);
    negatives->Save(out);
    electrons->Save(out); 
    //pim->Save(out);
    //kminus->Save(out);
    positives->Save(out); 
    protons->Save(out);
    pip->Save(out);
    //kplus->Save(out);
    neutrals->Save(out);
    gammas->Save(out);
    out->Close();



  }


protected:

  MomCorr_e1f             *momCorr; 
  PhysicsEventBuilder      builder; 
  ParticleFilter          *filter;
  Parameters              *params;
  SimpleNTupleWriter       tupleWriter; 

  //////////////////////////////////////////////////////
  // create standard histograms

  StandardHistograms *negatives;  
  StandardHistograms *electrons; 
  //StandardHistograms *pim;
  //StandardHistograms *kminus;
  StandardHistograms *positives; 
  StandardHistograms *protons;
  StandardHistograms *pip;
  //StandardHistograms *kplus;
  StandardHistograms *neutrals;
  StandardHistograms *gammas;

  //////////////////////////////////////////////////////


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

