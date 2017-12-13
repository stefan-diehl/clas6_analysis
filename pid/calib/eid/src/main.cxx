#include <iostream>
using namespace std;

#include "TApplication.h"
#include "h22Option.h"
#include "findElectronIDParameters.h"
#include "findElectronIDParameters.cxx"


int main(int argc, char *argv[]){

  findElectronIDParameters Analyzer;
  if (argc < 2) { cout << "No files?" << endl; return 0; }

  h22Options opts; 
  opts.args["PARS"].args = "eid.pars";
  opts.args["PARS"].type = 1;
  opts.args["PARS"].name = "Output parameter file name"; 
  opts.set(argc, argv);

  for(std::string filename : opts.ifiles){
    Analyzer.AddFile(filename);

    std::cout << "Adding file " << filename << std::endl; 
  }
        
  Analyzer.Init(); 
  Analyzer.Loop(opts.args["N"].arg); 
  Analyzer.WriteHardCodedParameters(); 
  Analyzer.CalculateSamplingFractionCut();
  Analyzer.CalculateCCThetaCut();
  Analyzer.SaveParameters(opts.args["PARS"].args); 
  Analyzer.SaveHistograms(opts.args["OUT"].args); 


  return 0; 
}
