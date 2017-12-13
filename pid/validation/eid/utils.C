void loadHistograms(){

  for(int c=0; c<NCONF; c++){
    files[c] = TFile::Open(inputFile[c].c_str()); 
    
    if(files[c]){
      cout << "[loadHistograms] Opened file " << inputFile[c] << endl; 
      zVertex[0][c]          = (TH1F*) files[c]->Get("h1_z_vertex_allNegatives_all"); 
      zVertex[1][c]          = (TH1F*) files[c]->Get("h1_z_vertex_cuts_all"); 
      samplingFraction[0][c] = (TH2F*) files[c]->Get("h_ec_edep_allNegatives_all"); 
      samplingFraction[1][c] = (TH2F*) files[c]->Get("h_ec_edep_cuts_all"); 
      ecEdep[0][c]           = (TH2F*) files[c]->Get("h_inner_outer_allNegatives_all"); 
      ecEdep[1][c]           = (TH2F*) files[c]->Get("h_inner_outer_cuts_all"); 
    } 

    else {
      cout << "[loadHistograms] Trouble opening file " << inputFile[c] << endl; 
    }

  }

}
