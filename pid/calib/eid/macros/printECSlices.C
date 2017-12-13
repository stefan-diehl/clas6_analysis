{

  gROOT->LoadMacro("utils.C");
  
  TFile *inputFile = TFile::Open("../out.root");

  const int numberSectors = 6; 
  const int numberSlices  = 16;

  TH1D *slices[numberSectors][numberSlices];
  
  loadSlices(inputFile, slices, numberSectors, numberSlices, "ecSlices");
  printSlices(slices, numberSectors, numberSlices, 0);
}
