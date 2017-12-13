/* 

   Writing out Electron ID 
   histograms 
   
   March 22, 2016

   Updated: Dec 16, 2016 
   Updating to Analysis_v2 classes

 */

// C++ Libraries
#include <iostream>
#include <cstdlib>
#include <map>

// My Libraries
#include "CommonTools.h"
#include "h22Event.h"
#include "h22Reader.h"
#include "h22Option.h"
#include "Parameters.h"
#include "ParticleFilter.h"
#include "StatusBar.h"

// CERN Root Libraries
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TMVA/Reader.h"

class PIDHistograms {

public:
  PIDHistograms();

 ~PIDHistograms(){ 
   std::cout << "[PIDHistograms::~PIDHistograms] Destructor has been called. " << std::endl; 
 }
  
  const static int NSECT = 7;
  const static int NTYPE = 12;

  // 1-D 
  TH1D *h1_nphe[NTYPE][NSECT];
  TH1F *h1_ec_edep_inner[NTYPE][NSECT];
  TH1F *h1_ec_edep_outer[NTYPE][NSECT];
  TH1F *h1_p[NTYPE][NSECT];
  TH1F *h1_z_vertex[NTYPE][NSECT];

  // 2-D 
  TH2F *h2_cc_theta[NTYPE][NSECT];
  TH2F *h2_etot_p[NTYPE][NSECT];
  TH2F *h2_ang_fid[NTYPE][NSECT];
  TH2F *h2_ec_edep[NTYPE][NSECT];
  TH2F *h2_dcr1_fid[NTYPE][NSECT];
  TH2F *h2_dcr3_fid[NTYPE][NSECT];
  TH2F *h2_ec_fid[NTYPE][NSECT];
  

  void Fill(h22Event & event, int ipart, int cutType);
  void FillAllOthers(h22Event & event, int ipart, int cutType);
  void Save(string outputFilename);
};

PIDHistograms::PIDHistograms(){

  std::string type[NTYPE]   = {"allNegatives", "cuts","Z_VERTEX",
			       "CC_FID","CC_PHI","CC_THETA",
			       "DC_R1_FID","DC_R3_FID","EC_FID",
			       "EC_IN_OUT","EC_SAMPLING", "allOthers"};
  
  std::string sect[NSECT]  = {"all", "s1", "s2", "s3", "s4", "s5", "s6"};
  
  // initialize 
  for (int itype = 0; itype < NTYPE; itype++)
    for(int isect = 0; isect < NSECT; isect++)
      {
	// 1d
	h1_nphe[itype][isect]          = new TH1D(Form("h1_nphe_%s_%s",type[itype].c_str(),sect[isect].c_str()),Form("h1_nphe_%s_%s",type[itype].c_str(),sect[isect].c_str()),200,0,100);
	h1_ec_edep_inner[itype][isect] = new TH1F(Form("h1_ec_edep_inner_%s_%s",type[itype].c_str(),sect[isect].c_str()),Form("h1_ec_edep_inner_%s_%s",type[itype].c_str(),sect[isect].c_str()),200,0,4);
	h1_ec_edep_outer[itype][isect] = new TH1F(Form("h1_ec_edep_outer_%s_%s",type[itype].c_str(),sect[isect].c_str()),Form("h1_ec_edep_outer_%s_%s",type[itype].c_str(),sect[isect].c_str()),200,0,4);
	h1_p[itype][isect]             = new TH1F(Form("h1_p_%s_%s",type[itype].c_str(),sect[isect].c_str()),Form("h1_p_%s_%s",type[itype].c_str(),sect[isect].c_str()),200,0,5);
	h1_z_vertex[itype][isect]      = new TH1F(Form("h1_z_vertex_%s_%s",type[itype].c_str(),sect[isect].c_str()),Form("h1_z_vertex_%s_%s",type[itype].c_str(),sect[isect].c_str()),200,-35,-15);
	
	// 2d
	h2_cc_theta[itype][isect] = new TH2F(Form("h_cc_theta_%s_%s",type[itype].c_str(),sect[isect].c_str()),Form("h_cc_theta_%s_%s",type[itype].c_str(),sect[isect].c_str()),18,0,17,200,0,60);
	h2_etot_p[itype][isect]   = new TH2F(Form("h_etot_p_%s_%s",type[itype].c_str(),sect[isect].c_str()),Form("h_etot_p_%s_%s",type[itype].c_str(),sect[isect].c_str()),200,0,5,200,0.05,0.5);
	h2_ang_fid[itype][isect]  = new TH2F(Form("h_ang_fid_%s_%s",type[itype].c_str(),sect[isect].c_str()),Form("h_ang_fid_%s_%s",type[itype].c_str(),sect[isect].c_str()),200,-30,30,200,5,65);
	h2_ec_edep[itype][isect]  = new TH2F(Form("h_ec_edep_%s_%s",type[itype].c_str(),sect[isect].c_str()),Form("h_ec_edep_%s_%s",type[itype].c_str(),sect[isect].c_str()),200,0.01,0.5,200,0.01,0.5);
	h2_dcr1_fid[itype][isect] = new TH2F(Form("h_dcr1_fid_%s_%s",type[itype].c_str(),sect[isect].c_str()),Form("h_dcr1_fid_%s_%s",type[itype].c_str(),sect[isect].c_str()),200,-100,100,200,-100,100);
	h2_dcr3_fid[itype][isect] = new TH2F(Form("h_dcr3_fid_%s_%s",type[itype].c_str(),sect[isect].c_str()),Form("h_dcr3_fid_%s_%s",type[itype].c_str(),sect[isect].c_str()),200,50,400,200,-200,200);
	h2_ec_fid[itype][isect]   = new TH2F(Form("h_ec_fid_%s_%s",type[itype].c_str(),sect[isect].c_str()),Form("h_ec_fid_%s_%s",type[itype].c_str(),sect[isect].c_str()),200,-500,500,200,-500,500);
      }
}

void PIDHistograms::Fill(h22Event & event, int ipart, int cutType){
  
  int sector = event.dc_sect[ipart];
  
  // filling histograms for all negatives hN_abc[0][x]
  // 1-D 
  h1_nphe[cutType][0]          ->Fill(event.nphe[ipart]/10);
  h1_ec_edep_inner[cutType][0] ->Fill(event.ec_ei[ipart]);
  h1_ec_edep_outer[cutType][0] ->Fill(event.ec_eo[ipart]);
  h1_p[cutType][0]             ->Fill(event.p[ipart]);
  h1_z_vertex[cutType][0]      ->Fill(event.vz[ipart]);
  
  // 2-D 
  h2_cc_theta[cutType][0] ->Fill((event.cc_segm[ipart]%1000)/10, event.GetThetaCC(ipart));
  h2_etot_p[cutType][0]   ->Fill(event.p[ipart], event.etot[ipart]/event.p[ipart]);
  h2_ang_fid[cutType][0]  ->Fill(event.GetRelativePhi(ipart), event.GetTheta(ipart));
  h2_ec_edep[cutType][0]  ->Fill(event.ec_ei[ipart], event.ec_eo[ipart]);
  h2_dcr1_fid[cutType][0] ->Fill(event.tl1_x[ipart], event.tl1_y[ipart]);
  h2_dcr3_fid[cutType][0] ->Fill(event.tl3_x[ipart], event.tl3_y[ipart]);
  h2_ec_fid[cutType][0]   ->Fill(event.ech_x[ipart], event.ech_y[ipart]);
  
  if(sector > 0) {
    h1_nphe[cutType][sector]          ->Fill(event.nphe[ipart]/10);
    h1_ec_edep_inner[cutType][sector] ->Fill(event.ec_ei[ipart]);
    h1_ec_edep_outer[cutType][sector] ->Fill(event.ec_eo[ipart]);
    h1_p[cutType][sector]             ->Fill(event.p[ipart]);
    h1_z_vertex[cutType][sector]      ->Fill(event.vz[ipart]);
    
    h2_cc_theta[cutType][sector] ->Fill((event.cc_segm[ipart]%1000)/10, event.GetThetaCC(ipart));
    h2_etot_p[cutType][sector]   ->Fill(event.p[ipart], event.etot[ipart]/event.p[ipart]);
    h2_ang_fid[cutType][sector]  ->Fill(event.GetRelativePhi(ipart), event.GetTheta(ipart));
    h2_ec_edep[cutType][sector]  ->Fill(event.ec_ei[ipart], event.ec_eo[ipart]);
    h2_dcr1_fid[cutType][sector] ->Fill(event.tl1_x[ipart], event.tl1_y[ipart]);
    h2_dcr3_fid[cutType][sector] ->Fill(event.tl3_x[ipart], event.tl3_y[ipart]);
    h2_ec_fid[cutType][sector]   ->Fill(event.ech_x[ipart], event.ech_y[ipart]);
  } 
}

void PIDHistograms::FillAllOthers(h22Event & event, int ipart, int cutType){
  
  int sector = event.dc_sect[ipart];

  /* 
   0 - "allNegatives", 
   1 - "cuts",
   2 - "Z_VERTEX",
   3 - "CC_FID",
   4 - "CC_PHI",
   5 - "CC_THETA",
   6 - "DC_R1_FID",
   7 - "DC_R3_FID",
   8 - "EC_FID",
   9 - "EC_IN_OUT",
  10 - "EC_SAMPLING", 
  11 - "allOthers"
  */
  
  if(cutType == 2){
    h1_z_vertex[NTYPE-1][0]->Fill(event.vz[ipart]);
    if(sector > 0) { h1_z_vertex[NTYPE-1][sector]->Fill(event.vz[ipart]); }
  }
  else if (cutType == 3){
    h2_ang_fid[NTYPE-1][0]->Fill(event.GetRelativePhi(ipart), event.GetTheta(ipart));    
    if (sector > 0) { h2_ang_fid[NTYPE-1][sector]->Fill(event.GetRelativePhi(ipart), event.GetTheta(ipart)); }
  }
  else if (cutType == 4){
    // none
  }
  else if (cutType == 5){
    h2_cc_theta[NTYPE-1][0] ->Fill((event.cc_segm[ipart]%1000)/10, event.GetThetaCC(ipart));
    if (sector > 0) { h2_cc_theta[NTYPE-1][sector] ->Fill((event.cc_segm[ipart]%1000)/10, event.GetThetaCC(ipart));} 
  }
  else if (cutType == 6){
  h2_dcr1_fid[NTYPE-1][0] ->Fill(event.tl1_x[ipart], event.tl1_y[ipart]);
  if (sector > 0) {h2_dcr1_fid[NTYPE-1][sector] ->Fill(event.tl1_x[ipart], event.tl1_y[ipart]); }
  }
  else if (cutType == 7){
  h2_dcr3_fid[NTYPE-1][0] ->Fill(event.tl3_x[ipart], event.tl3_y[ipart]);
  if (sector > 0) { h2_dcr3_fid[NTYPE-1][sector] ->Fill(event.tl3_x[ipart], event.tl3_y[ipart]);}
  }
  else if (cutType == 8){
  h2_ec_fid[NTYPE-1][0]   ->Fill(event.ech_x[ipart], event.ech_y[ipart]);
  if (sector > 0) { h2_ec_fid[NTYPE-1][sector]   ->Fill(event.ech_x[ipart], event.ech_y[ipart]); }
  }
  else if (cutType == 9){
  h2_ec_edep[NTYPE-1][0]  ->Fill(event.ec_ei[ipart], event.ec_eo[ipart]);
  if (sector > 0) { h2_ec_edep[NTYPE-1][sector]  ->Fill(event.ec_ei[ipart], event.ec_eo[ipart]); } 
  }
  else if (cutType == 10){
    h2_etot_p[NTYPE-1][0]   ->Fill(event.p[ipart], event.etot[ipart]/event.p[ipart]);
    if(sector > 0){ h2_etot_p[NTYPE-1][sector]   ->Fill(event.p[ipart], event.etot[ipart]/event.p[ipart]); }
  }

}

void PIDHistograms::Save(string outputFilename){
  
  TFile outfile(outputFilename.c_str(),"recreate");  
  
  for (int itype = 0; itype < NTYPE; itype++){
    for(int isect = 0; isect < NSECT; isect++){
      // 1d
      h1_nphe[itype][isect]->Write();          
      h1_ec_edep_inner[itype][isect]->Write();
      h1_ec_edep_outer[itype][isect]->Write();
      h1_p[itype][isect]->Write();             
      h1_z_vertex[itype][isect]->Write();      
      
      // 2d
      h2_cc_theta[itype][isect]->Write(); 
      h2_etot_p[itype][isect]->Write();   
      h2_ang_fid[itype][isect]->Write();  
      h2_ec_edep[itype][isect]->Write();  
      h2_dcr1_fid[itype][isect]->Write(); 
      h2_dcr3_fid[itype][isect]->Write(); 
      h2_ec_fid[itype][isect]->Write();   
    }
  }
  
  outfile.Write(); 
  outfile.Close(); 
}

int main (int argc, char * argv[]){

  h22Options opts; 
  opts.args["OUT"].args = "PID.root";

  opts.args["PARS"].args = "/u/home/dmriser/Analysis_v2/lists/parameters/data/final.pars";
  opts.args["PARS"].type = 1;
  opts.args["PARS"].name = "Parameter file";

  opts.args["LIST"].type = 1; 
  opts.args["LIST"].args = "none"; 
  opts.args["LIST"].name = "File list"; 

  opts.set(argc, argv);
  
  int nFiles        = opts.args["N"].arg;
  string files      = opts.args["LIST"].args;
  string outputFile = opts.args["OUT"].args;

  if (nFiles == 0 || opts.args["LIST"].args == "none"){
    std::cerr << " >> Please input files using -LIST=files.dat and -N=1e6. " << std::endl; 
    exit(0);  
  }

  std::vector<std::string> fileVector = loadFilesFromList(files, nFiles);

  // setup file reader and add files
  h22Reader reader;
  for(std::string file : fileVector){
    reader.AddFile(file); 
  }

  reader.Init(); // set branch addresses
  int nEvents = reader.GetEntries(); 

  // Load the correct parameters, if using 
  // E1-6, you need to pass in parameters for 
  // that dataset (z-vertex is much different for example). 
  Parameters *pars = new Parameters(); 
  pars->loadParameters(opts.args["PARS"].args);
  
  ParticleFilter filter(pars);
  filter.set_info(reader.GSIM, reader.GetRunNumber()); 
  //  filter.getCut("Track Quality Cut")->disable();

  PIDHistograms histos;
  StatusBar     statusBar; 
  
  // loop over events
  for (int iEvent = 0; iEvent < nEvents; iEvent++){
      reader.GetEntry(iEvent);
      h22Event event = reader.GetEvent();

      // loop over all negatives in the event 
      for(int ipart = 0; ipart < event.gpart; ipart++){     
	  if (event.q[ipart] < 0){

	    //  holds the result of all cuts, intensive 
	    map<string, bool> eID_Status = filter.eid_map(event, ipart);
	    
	      histos.Fill(event, ipart, 0);
	      if(eID_Status["Z_VERTEX"]){    histos.Fill(event, ipart, 2); }
	      if(eID_Status["CC_FID"]){      histos.Fill(event, ipart, 3); }
	      if(eID_Status["CC_PHI"]){      histos.Fill(event, ipart, 4); }
	      if(eID_Status["CC_THETA"]){    histos.Fill(event, ipart, 5); }
	      if(eID_Status["DC_R1_FID"]){   histos.Fill(event, ipart, 6); }
	      if(eID_Status["DC_R3_FID"]){   histos.Fill(event, ipart, 7); }
	      if(eID_Status["EC_FID"]){      histos.Fill(event, ipart, 8); }
	      if(eID_Status["EC_IN_OUT"]){   histos.Fill(event, ipart, 9); }
	      if(eID_Status["EC_SAMPLING"]){ histos.Fill(event, ipart, 10); }
	      /*
	      if(eID_Status["CC_FID"] && 
		 eID_Status["EC_FID"] && 
		 eID_Status["DC_R1_FID"] && 
		 histos.TMVAReader->EvaluateMVA("SVM") > 0.75){ histos.Fill(event, ipart, 11); }
	      */

	      // we have to start checking the annoying status things 
	      // This is the ugliest thing ever written in the history of the universe.

	      if (eID_Status["CC_FID"]     && 
		  eID_Status["CC_PHI"]     && 
		  eID_Status["CC_THETA"]   && 
		  eID_Status["DC_R1_FID"]  && 
		  eID_Status["DC_R3_FID"]  && 
		  eID_Status["EC_FID"]     && 
		  eID_Status["EC_IN_OUT"]  &&
		  eID_Status["EC_SAMPLING"]){

		histos.FillAllOthers(event, ipart, 2);
	      }

	      if (eID_Status["Z_VERTEX"]   && 
		  eID_Status["CC_PHI"]     && 
		  eID_Status["CC_THETA"]   && 
		  eID_Status["DC_R1_FID"]  && 
		  eID_Status["DC_R3_FID"]  && 
		  eID_Status["EC_FID"]     && 
		  eID_Status["EC_IN_OUT"]  &&
		  eID_Status["EC_SAMPLING"]){

		histos.FillAllOthers(event, ipart, 3);
	      }

	      if (eID_Status["Z_VERTEX"]   && 
		  eID_Status["CC_FID"]     && 
		  eID_Status["CC_PHI"]     && 
		  eID_Status["DC_R1_FID"]  && 
		  eID_Status["DC_R3_FID"]  && 
		  eID_Status["EC_FID"]     && 
		  eID_Status["EC_IN_OUT"]  &&
		  eID_Status["EC_SAMPLING"]){

		histos.FillAllOthers(event, ipart, 5);
	      }

	      if (eID_Status["Z_VERTEX"]   && 
		  eID_Status["CC_FID"]     && 
		  eID_Status["CC_PHI"]     && 
		  eID_Status["CC_THETA"]   && 
		  eID_Status["DC_R1_FID"]  && 
		  eID_Status["DC_R3_FID"]  && 
		  eID_Status["EC_FID"]     && 
		  eID_Status["EC_IN_OUT"]  &&
		  eID_Status["EC_SAMPLING"]){

		histos.FillAllOthers(event, ipart, 6);
	      }

	      if (eID_Status["Z_VERTEX"]   && 
		  eID_Status["CC_FID"]     && 
		  eID_Status["CC_PHI"]     && 
		  eID_Status["CC_THETA"]   && 
		  eID_Status["DC_R3_FID"]  && 
		  eID_Status["EC_FID"]     && 
		  eID_Status["EC_IN_OUT"]  &&
		  eID_Status["EC_SAMPLING"]){

		histos.FillAllOthers(event, ipart, 7);
	      }

	      if (eID_Status["Z_VERTEX"]   && 
		  eID_Status["CC_FID"]     && 
		  eID_Status["CC_PHI"]     && 
		  eID_Status["CC_THETA"]   && 
		  eID_Status["DC_R1_FID"]  && 
		  eID_Status["EC_IN_OUT"]  &&
		  eID_Status["EC_SAMPLING"]){

		histos.FillAllOthers(event, ipart, 8);
	      }

	      if (eID_Status["Z_VERTEX"]   && 
		  eID_Status["CC_FID"]     && 
		  eID_Status["CC_PHI"]     && 
		  eID_Status["CC_THETA"]   && 
		  eID_Status["DC_R1_FID"]  && 
		  eID_Status["DC_R3_FID"]  && 
		  eID_Status["EC_FID"]     && 
		  eID_Status["EC_SAMPLING"]){

		histos.FillAllOthers(event, ipart, 9);
	      }

	      if (eID_Status["Z_VERTEX"]   && 
		  eID_Status["CC_FID"]     && 
		  eID_Status["CC_PHI"]     && 
		  eID_Status["CC_THETA"]   && 
		  eID_Status["DC_R1_FID"]  && 
		  eID_Status["DC_R3_FID"]  && 
		  eID_Status["EC_FID"]     && 
		  eID_Status["EC_IN_OUT"]){

		histos.FillAllOthers(event, ipart, 10);
	      }


	  }
	} // end ipart loop 
      

      // look for electron in event
      int e_index = filter.getByPID(event, 11);
      if (e_index > -123){ histos.Fill(event, e_index, 1); }

      if (iEvent%10000 == 0) {
	statusBar.PrintStatus(iEvent, nEvents);
      }
    }  // end loop over events

  histos.Save(outputFile); 


  return 0;
}
 

