////////////////////////////////////////
/*
 David Riser, University of Connecticut
 
 July 17, 2016
 
 hid_histos.cxx -> description here.
 
 */
////////////////////////////////////////

// c++ includes
#include <iostream>
using namespace std;

// my includes from h22 library
#include "CommonTools.h"
#include "Corrections.h"
#include "h22Event.h"
#include "h22Option.h"
#include "h22Reader.h"
#include "pars.h"
#include "ParticleFilter.h"

// root includes
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLine.h"
#include "TStopwatch.h"

int main(int argc, char * argv[]){
    
    // Setup Options
    h22Options opts;
    opts.set(argc,argv);
    int GSIM     = opts.args["MC"].arg;
    long int nev = opts.args["N"].arg;
    string eparfile = opts.args["EPARS"].args;
    string hparfile = opts.args["HPARS"].args;
    string outfile  = opts.args["OUT"].args;
    
    hpars hpars; epars epars;
    epars.load(eparfile);
    hpars.load(hparfile);
    
    
    // Setup Reader
    h22Reader fReader(GSIM);
    for (auto it=opts.ifiles.begin(); it<opts.ifiles.end(); it++) { fReader.AddFile(*it); }
    fReader.Init();
    TFile out(outfile.c_str(),"recreate",outfile.c_str(),1);
    
    // If you ask for more the loop will analyze garbage
    nev = smallest(nev, fReader.GetEntries());
    
    // Setting important constants
    int runno = fReader.GetRunNumber();
    
    // Setting up PID and Corrections
    ParticleFilter filter(eparfile);
    filter.set_info(GSIM, runno);
    Corrections corr;
    
    // Histograms
    // [3][7] = [0:raw, 1:this, 2:all][0:All Sectors, 1-6:Sector #]
    
    //! 1-D
    TH1F * h1_dvz[3][7];
    
    //! 2-D
    TH2F * h2_dcr1[3][7];
    TH2F * h2_p_beta_prot[3][7];
    TH2F * h2_p_beta_pip[3][7];
    TH2F * h2_p_beta_pim[3][7];
    TH2F * h2_p_db_prot[3][7];
    TH2F * h2_p_db_pip[3][7];
    TH2F * h2_p_db_pim[3][7];
    
    
    // Initialize
    string cut[3] = {"raw","this","all"};
    
    for (int c=0; c<3; c++)
        for (int s=0; s<7; s++)
        {
            //! 1-D
            h1_dvz[c][s] = new TH1F(Form("h1_dvz_%s_%d",cut[c].c_str(),s),Form(" #Delta Z-Vertex %s Sector %d ",cut[c].c_str(),s),200,-15,15);
            
            //! 2-D
            h2_dcr1[c][s]        = new TH2F(Form("h2_dcr1_%s_%d",cut[c].c_str(),s),Form(" DC Region 1 Fid. %s Sector %d ",cut[c].c_str(),s),         100,15,60,100,-50,50);
            h2_p_beta_prot[c][s] = new TH2F(Form("h2_p_beta_prot_%s_%d",cut[c].c_str(),s), Form(" #beta vs P (prot) %s Sector %d ",cut[c].c_str(),s),100,0.5,5,100,0,1);
            h2_p_beta_pip[c][s]  = new TH2F(Form("h2_p_beta_pip_%s_%d",cut[c].c_str(),s), Form(" #beta vs P (#pi+) %s Sector %d ",cut[c].c_str(),s), 100,0.5,5,100,0,1);
            h2_p_beta_pim[c][s]  = new TH2F(Form("h2_p_beta_pim_%s_%d",cut[c].c_str(),s), Form(" #beta vs P (#pi-) %s Sector %d ",cut[c].c_str(),s), 100,0.5,5,100,0,1);
            
            h2_p_db_prot[c][s]   = new TH2F(Form("h2_p_db_prot_%s_%d",cut[c].c_str(),s), Form(" #Delta #beta vs P (prot) %s Sector %d ",cut[c].c_str(),s),100,0.5,5,100,-0.6,0.2);
            h2_p_db_pip[c][s]    = new TH2F(Form("h2_p_db_pip_%s_%d",cut[c].c_str(),s), Form(" #Delta #beta vs P (#pi+) %s Sector %d ",cut[c].c_str(),s), 100,0.5,5,100,-0.2,0.6);
            h2_p_db_pim[c][s]    = new TH2F(Form("h2_p_db_pim_%s_%d",cut[c].c_str(),s), Form(" #Delta #beta vs P (#pi-) %s Sector %d ",cut[c].c_str(),s), 100,0.5,5,100,-0.2,0.2);
            
        }
    
    // Keep track of time.
    TStopwatch timer;
    timer.Reset();
    timer.Start();
    
    for (long int iev=0; iev<nev; iev++)
    {
        fReader.GetEntry(iev);
        h22Event event = fReader.GetEvent();
        
        // keeping track of run number and making sure particle filter knows about it as well.
        if ( runno != fReader.GetRunNumber() ){ runno = fReader.GetRunNumber(); filter.set_info(GSIM, runno); }
        
        // Load up hadrons if we've electron.
        if (filter.has_electron(event))
        {
            double start_time = corr.electron_sct(event,0,runno,GSIM) - event.sc_r[0]/speed_of_light;
            
            for (int ipart=1; ipart<event.gpart; ipart++)
            {
                double beta = (event.sc_r[ipart]/(corr.hadron_sct(event,ipart,runno,GSIM)-start_time))/speed_of_light;
                int sector  = event.dc_sect[ipart];
                
                // Doing this cut
                map<string,bool> results = filter.hid_map(event, ipart);
                
                if (results["PROT_DVZ"])   { h1_dvz[1][0]->Fill(corr.vz(event,0,runno,GSIM)-corr.vz(event,ipart,runno,GSIM)); h1_dvz[1][sector]->Fill(corr.vz(event,0,runno,GSIM)-corr.vz(event,ipart,runno,GSIM)); }
                if (results["PROT_DCFID"]) { h2_dcr1[1][0]->Fill(event.GetRotatedDCR1PosX(ipart), event.GetRotatedDCR1PosY(ipart)); h2_dcr1[1][sector]->Fill(event.GetRotatedDCR1PosX(ipart), event.GetRotatedDCR1PosY(ipart)); }
                
                //! Fill Raw
                h1_dvz[0][0]  ->Fill(corr.vz(event,0,runno,GSIM) - corr.vz(event,ipart,runno,GSIM));
                h2_dcr1[0][0] ->Fill(event.GetRotatedDCR1PosX(ipart), event.GetRotatedDCR1PosY(ipart));
                
                h1_dvz[0][sector]  ->Fill(corr.vz(event,0,runno,GSIM) - corr.vz(event,ipart,runno,GSIM));
                h2_dcr1[0][sector] ->Fill(event.GetRotatedDCR1PosX(ipart), event.GetRotatedDCR1PosY(ipart));
                
                if (event.q[ipart] > 0)
                {
                    h2_p_db_prot[0][0]->Fill(event.p[ipart],(event.p[ipart]/sqrt(event.p[ipart]*event.p[ipart] + pid_to_mass(2212)*pid_to_mass(2212)))-beta);
                    h2_p_db_pip[0][0] ->Fill(event.p[ipart],(event.p[ipart]/sqrt(event.p[ipart]*event.p[ipart] + pid_to_mass(211)*pid_to_mass(211)))-beta);
                    h2_p_beta_prot[0][0]->Fill(event.p[ipart],beta);
                    h2_p_beta_pip[0][0] ->Fill(event.p[ipart],beta);
                    
                    h2_p_db_prot[0][sector]->Fill(event.p[ipart],(event.p[ipart]/sqrt(event.p[ipart]*event.p[ipart] + pid_to_mass(2212)*pid_to_mass(2212)))-beta);
                    h2_p_db_pip[0][sector] ->Fill(event.p[ipart],(event.p[ipart]/sqrt(event.p[ipart]*event.p[ipart] + pid_to_mass(211)*pid_to_mass(211)))-beta);
                    h2_p_beta_prot[0][sector]->Fill(event.p[ipart],beta);
                    h2_p_beta_pip[0][sector] ->Fill(event.p[ipart],beta);
                    
                    if (results["PROT_DBETA"]) { h2_p_db_prot[1][0]->Fill(event.p[ipart],(event.p[ipart]/sqrt(event.p[ipart]*event.p[ipart] + pid_to_mass(2212)*pid_to_mass(2212)))-beta);}
                    if (results["PIP_DBETA"])  { h2_p_db_pip[1][0]->Fill(event.p[ipart],(event.p[ipart]/sqrt(event.p[ipart]*event.p[ipart] + pid_to_mass(211)*pid_to_mass(211)))-beta); }
                    if (results["PROT_DBETA"]) { h2_p_beta_prot[1][0]->Fill(event.p[ipart],beta);}
                    if (results["PIP_DBETA"])  { h2_p_beta_pip[1][0] ->Fill(event.p[ipart],beta);}
                    
                    if (results["PROT_DBETA"]) { h2_p_db_prot[1][sector]->Fill(event.p[ipart],(event.p[ipart]/sqrt(event.p[ipart]*event.p[ipart] + pid_to_mass(2212)*pid_to_mass(2212)))-beta);}
                    if (results["PIP_DBETA"])  { h2_p_db_pip[1][sector]->Fill(event.p[ipart],(event.p[ipart]/sqrt(event.p[ipart]*event.p[ipart] + pid_to_mass(211)*pid_to_mass(211)))-beta); }
                    if (results["PROT_DBETA"]) { h2_p_beta_prot[1][sector]->Fill(event.p[ipart],beta);}
                    if (results["PIP_DBETA"])  { h2_p_beta_pip[1][sector] ->Fill(event.p[ipart],beta);}
                    
                }
                
                if (event.q[ipart] < 0)
                {
                    h2_p_db_pim[0][0]   ->Fill(event.p[ipart],(event.p[ipart]/sqrt(event.p[ipart]*event.p[ipart] + pid_to_mass(-211)*pid_to_mass(-211)))-beta);
                    h2_p_beta_pim[0][0] ->Fill(event.p[ipart],beta);
                    h2_p_db_pim[0][sector]   ->Fill(event.p[ipart],(event.p[ipart]/sqrt(event.p[ipart]*event.p[ipart] + pid_to_mass(-211)*pid_to_mass(-211)))-beta);
                    h2_p_beta_pim[0][sector] ->Fill(event.p[ipart],beta);
                    
                    if (results["PIM_DBETA"]){ h2_p_db_pim[1][0]     ->Fill(event.p[ipart],(event.p[ipart]/sqrt(event.p[ipart]*event.p[ipart] + pid_to_mass(-211)*pid_to_mass(-211)))-beta); }
                    if (results["PIM_DBETA"]){ h2_p_beta_pim[1][0]   ->Fill(event.p[ipart],beta); }
                    if (results["PIM_DBETA"]){ h2_p_db_pim[1][sector]->Fill(event.p[ipart],(event.p[ipart]/sqrt(event.p[ipart]*event.p[ipart] + pid_to_mass(-211)*pid_to_mass(-211)))-beta); }
                    if (results["PIM_DBETA"]){ h2_p_beta_pim[1][sector]->Fill(event.p[ipart],beta); }
                    
                }
                
                
            }
            
            
            // Doing All Cuts Passed
            int prot_index = filter.getByPID(event,2212);
            int pip_index = filter.getByPID(event,211);
            int pim_index = filter.getByPID(event,-211);
            
            if (prot_index > -123)
            {
                double beta = (event.sc_r[prot_index]/(corr.hadron_sct(event,prot_index,runno,GSIM)-start_time))/speed_of_light;
                int sector  = event.dc_sect[prot_index];
                
                h1_dvz[2][0]              ->Fill(corr.vz(event,0,runno,GSIM) - corr.vz(event,prot_index,runno,GSIM));
                h2_dcr1[2][0]             ->Fill(event.GetRotatedDCR1PosX(prot_index), event.GetRotatedDCR1PosY(prot_index));
                h1_dvz[2][sector]         ->Fill(corr.vz(event,0,runno,GSIM) - corr.vz(event,prot_index,runno,GSIM));
                h2_dcr1[2][sector]        ->Fill(event.GetRotatedDCR1PosX(prot_index), event.GetRotatedDCR1PosY(prot_index));
                h2_p_db_prot[2][0]        ->Fill(event.p[prot_index],(event.p[prot_index]/sqrt(event.p[prot_index]*event.p[prot_index] + pid_to_mass(2212)*pid_to_mass(2212)))-beta);
                h2_p_beta_prot[2][0]      ->Fill(event.p[prot_index],beta);
                h2_p_db_prot[2][sector]   ->Fill(event.p[prot_index],(event.p[prot_index]/sqrt(event.p[prot_index]*event.p[prot_index] + pid_to_mass(2212)*pid_to_mass(2212)))-beta);
                h2_p_beta_prot[2][sector] ->Fill(event.p[prot_index],beta);
                
            }
            
            if (pip_index > -123)
            {
                double beta = (event.sc_r[pip_index]/(corr.hadron_sct(event,pip_index,runno,GSIM)-start_time))/speed_of_light;
                int sector  = event.dc_sect[pip_index];
                
                h1_dvz[2][0]              ->Fill(corr.vz(event,0,runno,GSIM) - corr.vz(event,pip_index,runno,GSIM));
                h2_dcr1[2][0]             ->Fill(event.GetRotatedDCR1PosX(pip_index), event.GetRotatedDCR1PosY(pip_index));
                h1_dvz[2][sector]         ->Fill(corr.vz(event,0,runno,GSIM) - corr.vz(event,pip_index,runno,GSIM));
                h2_dcr1[2][sector]        ->Fill(event.GetRotatedDCR1PosX(pip_index), event.GetRotatedDCR1PosY(pip_index));
                h2_p_db_pip[2][0]         ->Fill(event.p[pip_index],(event.p[pip_index]/sqrt(event.p[pip_index]*event.p[pip_index] + pid_to_mass(211)*pid_to_mass(211)))-beta);
                h2_p_beta_pip[2][0]       ->Fill(event.p[pip_index],beta);
                h2_p_db_pip[2][sector]    ->Fill(event.p[pip_index],(event.p[pip_index]/sqrt(event.p[pip_index]*event.p[pip_index] + pid_to_mass(211)*pid_to_mass(211)))-beta);
                h2_p_beta_pip[2][sector]  ->Fill(event.p[pip_index],beta);
            }
            
            if (pim_index > -123)
            {
                double beta = (event.sc_r[pim_index]/(corr.hadron_sct(event,pim_index,runno,GSIM)-start_time))/speed_of_light;
                int sector  = event.dc_sect[pim_index];
                
                h1_dvz[2][0]              ->Fill(corr.vz(event,0,runno,GSIM) - corr.vz(event,pim_index,runno,GSIM));
                h2_dcr1[2][0]             ->Fill(event.GetRotatedDCR1PosX(pim_index), event.GetRotatedDCR1PosY(pim_index));
                h1_dvz[2][sector]         ->Fill(corr.vz(event,0,runno,GSIM) - corr.vz(event,pim_index,runno,GSIM));
                h2_dcr1[2][sector]        ->Fill(event.GetRotatedDCR1PosX(pim_index), event.GetRotatedDCR1PosY(pim_index));
                h2_p_db_pim[2][0]         ->Fill(event.p[pim_index],(event.p[pim_index]/sqrt(event.p[pim_index]*event.p[pim_index] + pid_to_mass(-211)*pid_to_mass(-211)))-beta);
                h2_p_beta_pim[2][0]       ->Fill(event.p[pim_index],beta);
                h2_p_db_pim[2][sector]    ->Fill(event.p[pim_index],(event.p[pim_index]/sqrt(event.p[pim_index]*event.p[pim_index] + pid_to_mass(-211)*pid_to_mass(-211)))-beta);
                h2_p_beta_pim[2][sector]  ->Fill(event.p[pim_index],beta);
                
            }
            
        }
        
        // Tell the user
        if ( (int)iev%(int)opts.args["PRINT"].arg == 0) cout << "\r done " << iev << " of " << nev << flush;
    }
    
    //  Performance readout
    double loop_time  = timer.RealTime();
    double event_rate = (double)nev/loop_time;
    cout << " Event loop finished in " << loop_time << " seconds w/ rate " << event_rate << " events/sec " << endl;
    
    out.Write();
    out.Close();
    
    
    
}
