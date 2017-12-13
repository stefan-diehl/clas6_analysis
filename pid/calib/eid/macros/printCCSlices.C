{

  gROOT->LoadMacro("utils.C");
  
  TFile *inputFile = TFile::Open("MCInelastic.root");

  const int numberSectors = 6; 
  const int numberSlices  = 18;

  TH1D *slices[numberSectors][numberSlices];
  
  loadSlices(inputFile, slices, numberSectors, numberSlices, "ccSlices");
  printSlices(slices, numberSectors, numberSlices, 1);
}
