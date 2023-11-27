#include "DC_fiducial_cuts.cxx"
#include <cstdlib>
#include <iostream>
#include <chrono>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include "clas12reader.h"
#include <optional>

using namespace clas12; 
// Generalise what to add to tree-> Could have a file that initialises the tree before the macro is run.
// This class should contain all functions needed to perform the Energy Corrections for the Forward Tagger.

// Preliminary Electron correction from Geraint's calculations. Needs to be calculated again for pass 2 data.
TLorentzVector  Correct_Electron(TLorentzVector x){

  Double_t E_new, Px_el, Py_el, Pz_el;
  TLorentzVector el_new;

  E_new = -1.1935+2.297*x.E()-0.4889*pow(x.E(),2)+0.0805*pow(x.E(),3)-0.0048*pow(x.E(),4);

  Px_el = E_new*(x.Px()/x.Rho());
  Py_el = E_new*(x.Py()/x.Rho());
  Pz_el = E_new*(x.Pz()/x.Rho());

  el_new.SetXYZM(Px_el, Py_el, Pz_el, 0.000511);

  return el_new;
}

TLorentzVector Correct_z_position(TLorentzVector e, Double_t x, Double_t y, Double_t z, Double_t r){

  Double_t px, py, pz;
  TLorentzVector el_corZ;
  Double_t mom=e.E()*e.E()-0.000511*0.000511;
  mom=TMath::Sqrt(mom);

  px = mom*(x/r);
  py = mom*(y/r);
  pz = mom*(z/r);

  el_corZ.SetXYZM(px, py, pz, 0.000511);

  return el_corZ;
}

Double_t Correct_Theta(TLorentzVector e){

  Double_t Th_cor, Corect_Th;

  Th_cor = exp(1.797 - 4.485*e.E()) + exp(-0.8671 - 1.078*e.E());
  Th_cor = TMath::DegToRad()*Th_cor;

  Corect_Th = e.Theta() + Th_cor;

  return Corect_Th;
}

Double_t Correct_Phi(TLorentzVector e){

  Double_t Ph_cor, Corect_Ph;

  Ph_cor = exp(4.918 - 3.828*e.E()) + exp(3.841 - 1.256*e.E()) + exp(2.874 - 0.2195*e.E());
  Ph_cor = TMath::DegToRad()*Ph_cor;

  Corect_Ph = e.Phi() + Ph_cor;

  return Corect_Ph;
}
// Function to do fiducial cuts.
void FidCuts(int part_DC_sector, Double_t part_DC_c1x, Double_t part_DC_c1y, Double_t part_DC_c1z, Double_t part_DC_c2x, Double_t part_DC_c2y, Double_t part_DC_c2z, Double_t part_DC_c3x, Double_t part_DC_c3y, Double_t part_DC_c3z, Int_t part_pid, bool outbending){
  
  while(false==outbending){
         if(!DC_fiducial_cut_theta_phi(part_DC_c1x, part_DC_c1y, part_DC_c1z, part_DC_c2x, part_DC_c2y, part_DC_c2z, part_DC_c3x, part_DC_c3y, part_DC_c3z, part_pid, part_DC_sector,1) ||
            !DC_fiducial_cut_theta_phi(part_DC_c1x, part_DC_c1y, part_DC_c1z, part_DC_c2x, part_DC_c2y, part_DC_c2z, part_DC_c3x, part_DC_c3y, part_DC_c3z, part_pid, part_DC_sector,2) ||
            !DC_fiducial_cut_theta_phi(part_DC_c1x, part_DC_c1y, part_DC_c1z, part_DC_c2x, part_DC_c2y, part_DC_c2z, part_DC_c3x, part_DC_c3y, part_DC_c3z, part_pid, part_DC_sector,3)) continue;
         }
         
      
  while(true==outbending){
         if(!DC_fiducial_cut_XY(part_DC_c1x, part_DC_c1y, part_DC_c1z, part_DC_c2x, part_DC_c2y, part_DC_c2z, part_DC_c3x, part_DC_c3y, part_DC_c3z, part_pid, part_DC_sector,1) &&
            !DC_fiducial_cut_XY(part_DC_c1x, part_DC_c1y, part_DC_c1z, part_DC_c2x, part_DC_c2y, part_DC_c2z, part_DC_c3x, part_DC_c3y, part_DC_c3z, part_pid, part_DC_sector,2) &&
            !DC_fiducial_cut_XY(part_DC_c1x, part_DC_c1y, part_DC_c1z, part_DC_c2x, part_DC_c2y, part_DC_c2z, part_DC_c3x, part_DC_c3y, part_DC_c3z, part_pid, part_DC_sector,3)) continue;
         }
          
}
// Generalised function to get hipo files from runs and do the cuts.
// FilePath: Directory of the runs.
void MMall_Slurm(TString infile, TString outfile){
  // Records the number of events that are put in the tree
  Int_t Tree_Events=0; 
  // Lorentz Vectors of particles.
  TLorentzVector el, el_corZ, el_new;
  TLorentzVector pip;
  TLorentzVector pim;
  TLorentzVector pr;
  TLorentzVector miss_el;
  // Exclusive missing mass.
  TLorentzVector miss_all;
  // Exclusive missing mass as a vector to write to tree.
  vector<double> miss_all_m2;
  TLorentzVector beam(0,0,10.6,10.6);
  // Set 4 vector for target.
  TLorentzVector target(0,0,0,0.938272); 
  // To hold number of events.
  Int_t eventno; 
  // Create ROOT file to save the tree to be filled.
  TFile tree_file(outfile,"recreate");
  // Create new Tree.
  TTree *Tree = new TTree("Tree","tree"); 
  // Branches to relevant variables.
  Tree->Branch("eventno",&eventno);
  Tree->Branch("miss_all_m2",&miss_all_m2);

  Double_t x, y, z, r;
  Int_t part_pid;
  Double_t part_DC_c1x,part_DC_c1y,part_DC_c1z;
  Double_t part_DC_c2x,part_DC_c2y,part_DC_c2z;
  Double_t part_DC_c3x,part_DC_c3y,part_DC_c3z;

  // Initialise HIPO chain.
  TChain fake("hipo");
  TString inputFile(infile);
  // Looping over all the hipo files in the range.
  //for(int i=range[0]; i<range[1]; i++){
    //TString inputFile(FilePath+ Form("%i", i) + "/*.hipo");
    // Add everything to a HIPO chain.
  fake.Add(infile.Data());
  //}
  // Get the list of all the HIPO files from the chain.
  auto files=fake.GetListOfFiles();
  // Particle information Database.
  auto db=TDatabasePDG::Instance();
  
  for(Int_t i=0;i<files->GetEntries();i++){
    // Counts events.
    int counter=0;
    //create the event reader.
    clas12reader c12(files->At(i)->GetTitle());
    // Print the file in loop.
    cout<<"File in loop: "<<files->At(i)->GetTitle()<<endl;
    //loop over all events in the file
    while(c12.next()==true){
      counter++;
      // Clear the vectors from the previous event.
      miss_all_m2.clear();
      // Getting the event number.
      eventno = c12.runconfig()->getEvent(); 
      // Continue to the next if no particles found.
      if(c12.getDetParticles().empty())
      	continue;
      // Get particle information from the HIPO files.
      auto electron=c12.getByID(11);
      auto proton=c12.getByID(2212);
      auto piplus=c12.getByID(211);
      auto piminus=c12.getByID(-211);
      // getStatus finds where the particle has been detected. {0,2000}-> only Forward Tagger, {2000,4000}-> only Forward Detector.
      if (electron.size()!=1  || proton.size()!=1 || piplus.size()!=1 || piminus.size()!=1
      || abs(electron[0]->par()->getStatus())>1999
      || abs(piplus[0]->par()->getStatus())<2000 || abs(piplus[0]->par()->getStatus())>3999
      || abs(piminus[0]->par()->getStatus())<2000 || abs(piminus[0]->par()->getStatus())>3999
      || abs(proton[0]->par()->getStatus())<2000 || abs(proton[0]->par()->getStatus())>3999) continue;
      // Set 4-Vector values from particle information obtained.
      el.SetXYZM(electron[0]->par()->getPx(),electron[0]->par()->getPy(),electron[0]->par()->getPz(),0.000511);
      pip.SetXYZM(piplus[0]->par()->getPx(),piplus[0]->par()->getPy(),piplus[0]->par()->getPz(),0.139600);
      pim.SetXYZM(piminus[0]->par()->getPx(),piminus[0]->par()->getPy(),piminus[0]->par()->getPz(),0.139600);
      pr.SetXYZM(proton[0]->par()->getPx(),proton[0]->par()->getPy(),proton[0]->par()->getPz(),0.938300);
      // Events in all 3 layers of drift chamber, for the fiducial cuts.
      if(piplus[0]->traj(DC,6)->getDetector() == 6 && piplus[0]->traj(DC,6)->getLayer() == 6){
        part_DC_c1x = piplus[0]->traj(DC,6)->getX();
        part_DC_c1y = piplus[0]->traj(DC,6)->getY();
        part_DC_c1z = piplus[0]->traj(DC,6)->getZ();
      }
      if(piplus[0]->traj(DC,18)->getDetector() == 6 && piplus[0]->traj(DC,18)->getLayer() == 18){
        part_DC_c2x = piplus[0]->traj(DC,18)->getX();
        part_DC_c2y = piplus[0]->traj(DC,18)->getY();
        part_DC_c2z = piplus[0]->traj(DC,18)->getZ();
      }
      if(piplus[0]->traj(DC,36)->getDetector() == 6 && piplus[0]->traj(DC,36)->getLayer() == 36){
        part_DC_c3x = piplus[0]->traj(DC,36)->getX();
        part_DC_c3y = piplus[0]->traj(DC,36)->getY();
        part_DC_c3z = piplus[0]->traj(DC,36)->getZ();
      }
       // Positive particles asssigned to pi^+, negative to pi^-.
      part_pid = 0;
      part_pid = 211;
      // Fiducial cuts for pi^+.
      int part_DC_sector_pip = determineSector(part_DC_c2x, part_DC_c2y,part_DC_c2z);

      // Fiducial cuts for inbending or outbending data, based on user input.
      FidCuts(part_DC_sector_pip, part_DC_c1x, part_DC_c1y, part_DC_c1z, part_DC_c2x, part_DC_c2y, part_DC_c2z, part_DC_c3x, part_DC_c3y, part_DC_c3z,  part_pid, true);
        //pi^-.
      if(piminus[0]->traj(DC,6)->getDetector() == 6 && piminus[0]->traj(DC,6)->getLayer() == 6){
        part_DC_c1x = piminus[0]->traj(DC,6)->getX();
        part_DC_c1y = piminus[0]->traj(DC,6)->getY();
        part_DC_c1z = piminus[0]->traj(DC,6)->getZ();
      }
      if(piminus[0]->traj(DC,18)->getDetector() == 6 && piminus[0]->traj(DC,18)->getLayer() == 18){
        part_DC_c2x = piminus[0]->traj(DC,18)->getX();
        part_DC_c2y = piminus[0]->traj(DC,18)->getY();
        part_DC_c2z = piminus[0]->traj(DC,18)->getZ();
      }
      if(piminus[0]->traj(DC,36)->getDetector() == 6 && piminus[0]->traj(DC,36)->getLayer() == 36){
        part_DC_c3x = piminus[0]->traj(DC,36)->getX();
        part_DC_c3y = piminus[0]->traj(DC,36)->getY();
        part_DC_c3z = piminus[0]->traj(DC,36)->getZ();
      }
      // Positive particles asssigned to pi^+, negative to pi^-
      part_pid = 0;
      part_pid = -211;
      // Fiducial cuts for pi^-.
      int part_DC_sector_pim = determineSector(part_DC_c2x, part_DC_c2y,part_DC_c2z);

      FidCuts(part_DC_sector_pim, part_DC_c1x, part_DC_c1y, part_DC_c1z, part_DC_c2x, part_DC_c2y, part_DC_c2z, part_DC_c3x, part_DC_c3y, part_DC_c3z,  part_pid, true);
      // Proton fiducial cut.
      if(proton[0]->traj(DC,6)->getDetector() == 6 && proton[0]->traj(DC,6)->getLayer() == 6){
        part_DC_c1x = proton[0]->traj(DC,6)->getX();
        part_DC_c1y = proton[0]->traj(DC,6)->getY();
        part_DC_c1z = proton[0]->traj(DC,6)->getZ();
      }
      if(proton[0]->traj(DC,18)->getDetector() == 6 && proton[0]->traj(DC,18)->getLayer() == 18){
        part_DC_c2x = proton[0]->traj(DC,18)->getX();
        part_DC_c2y = proton[0]->traj(DC,18)->getY();
        part_DC_c2z = proton[0]->traj(DC,18)->getZ();
      }
      if(proton[0]->traj(DC,36)->getDetector() == 6 && proton[0]->traj(DC,36)->getLayer() == 36){
        part_DC_c3x = proton[0]->traj(DC,36)->getX();
        part_DC_c3y = proton[0]->traj(DC,36)->getY();
        part_DC_c3z = proton[0]->traj(DC,36)->getZ();
      }
      // Positive particles asssigned to pi^+, negative to pi^-
      part_pid = 0;
      part_pid = 2212;
      int part_DC_sector_pr = determineSector(part_DC_c2x, part_DC_c2y,part_DC_c2z);

      FidCuts(part_DC_sector_pr, part_DC_c1x, part_DC_c1y, part_DC_c1z, part_DC_c2x, part_DC_c2y, part_DC_c2z, part_DC_c3x, part_DC_c3y, part_DC_c3z, part_pid, true);

      while(TMath::RadToDeg()*el.Theta() > 2.5 && TMath::RadToDeg()*el.Theta() < 4.5){
         miss_all = beam + target - el - pip - pim - pr;
         miss_all_m2.push_back(miss_all.M2());
         Tree->Fill();
         Tree_Events++;
        
         }

        }

      }

    Tree->Write(); 
    // Write information to the TTree
    cout<<"events in tree "<<Tree_Events<<endl;

}



