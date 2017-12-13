{
  
  // -------------------------------------------------
  // User Parameters 
  // -------------------------------------------------
  const int NCONF         = 1; 
  string inputFile[NCONF] = {"/volatile/clas12/dmriser/rootFiles/pid/eid/data.root"}; 
  string variation[NCONF] = {"Data"};
  string imagePath        = "/volatile/clas12/dmriser/plots/pid/eid/";


  // -------------------------------------------------
  // Functionality 
  // -------------------------------------------------
  gROOT->LoadMacro("utils.C"); 

  TFile *files[NCONF];  
  TH1D *zVertex[2][NCONF];
  TH2D *ecEdep[2][NCONF];
  TH2D *samplingFraction[2][NCONF];

  loadHistograms(); 

}
