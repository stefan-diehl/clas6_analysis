#include <iostream>
using namespace std;

#include "Photon.h"
#include "Photon.cxx"

int main(int argc, char * argv[]){

        Photon Analyzer;
        if (argc < 2) { return 0; }

        for (int ifile=1; ifile<argc; ifile++){ Analyzer.AddFile(argv[ifile]); } 
        Analyzer.Init(); 
        Analyzer.Loop(); 
	Analyzer.WriteParameters(); 
	Analyzer.histos.Save("photonIDHistograms.root");
        return 0; 
}
