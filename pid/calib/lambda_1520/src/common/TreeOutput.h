#ifndef tree_output_h 
#define tree_output_h 

// standard libs
#include <iostream>

// root includes 
#include "TBranch.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TTree.h"


class TreeOutput {
  
 public:
  TreeOutput(){
    tree = new TTree("events",""); 
    tree->SetDirectory(0); 

    b_p_prot       = tree->Branch("p_prot",       &p_prot); 
    b_p_kp         = tree->Branch("p_kp",         &p_kp); 
    b_p_km         = tree->Branch("p_km",         &p_km); 
    b_p_lambda     = tree->Branch("p_lambda",     &p_lambda);

    b_theta_kp     = tree->Branch("theta_kp", &theta_kp);
    b_theta_km     = tree->Branch("theta_km", &theta_km);

    b_phi_kp = tree->Branch("phi_kp", &phi_kp);
    
    b_angle_kp_km        = tree->Branch("angle_kp_km",   &angle_kp_km);
    b_angle_kp_prot      = tree->Branch("angle_kp_prot", &angle_kp_prot);
    b_angle_km_prot      = tree->Branch("angle_km_prot",  &angle_km_prot);
    b_angle_kp_lambda    = tree->Branch("angle_kp_lambda", &angle_kp_lambda);
    b_angle_km_lambda    = tree->Branch("angle_km_lambda", &angle_km_lambda);
    b_angle_prot_lambda  = tree->Branch("angle_prot_lambda", &angle_prot_lambda);
    b_angle_kp_lambda_cm = tree->Branch("angle_kp_lambda_cm", &angle_kp_lambda_cm);

    b_im_km_prot = tree->Branch("im_km_prot", &im_km_prot);
    b_im_km_kp   = tree->Branch("im_km_kp", &im_km_kp);

    b_x = tree->Branch("x", &x);
    b_w = tree->Branch("w", &w);
    b_q2 = tree->Branch("q2", &q2);
    b_z = tree->Branch("z", &z);
    b_pt = tree->Branch("pt", &pt);
    b_hel = tree->Branch("hel", &hel); 
    b_missing_mass = tree->Branch("missing_mass", &missing_mass);
  }

  TTree *tree; 

  // data 
  float p_prot; 
  float p_kp; 
  float p_km; 
  float p_lambda;
  int hel; 

  float theta_kp; 
  float theta_km; 
  float phi_kp; 
  float sig_kp, sig_prot; 

  float im_km_prot; 
  float im_km_kp; 

  float angle_kp_km;
  float angle_kp_prot; 
  float angle_km_prot; 
  float angle_kp_lambda_cm; 
  float angle_kp_lambda; 
  float angle_km_lambda; 
  float angle_prot_lambda; 

  float x, w, q2, z, pt; 
  float missing_mass; 

  // branches 
  TBranch *b_hel; 
  TBranch *b_p_prot; 
  TBranch *b_p_kp; 
  TBranch *b_p_km; 
  TBranch *b_p_lambda; 
  TBranch *b_theta_kp; 
  TBranch *b_theta_km; 
  TBranch *b_phi_kp; 
  TBranch *b_im_km_prot; 
  TBranch *b_im_km_kp; 
  TBranch *b_angle_kp_km; 
  TBranch *b_angle_kp_prot; 
  TBranch *b_angle_km_prot; 
  TBranch *b_angle_kp_lambda_cm; 
  TBranch *b_angle_kp_lambda;
  TBranch *b_angle_km_lambda;
  TBranch *b_angle_prot_lambda;
  TBranch *b_x; 
  TBranch *b_w; 
  TBranch *b_q2; 
  TBranch *b_z; 
  TBranch *b_pt; 
  TBranch *b_missing_mass; 

  void Save(TFile *out){
    if (out->IsOpen()){
      tree->Write(); 
    }
  }

};

#endif 
