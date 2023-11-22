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

using namespace clas12;

TLorentzVector  Correct_Electron(TLorentzVector x){ //corrections can vary over beam time -different current. May need to look up if the correction chnages for different cases.

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

void MMpr_Tree(){
// Make list of good files to run over. 
Int_t Tree_Events=0; // Records the number of events that are put in the tree
Int_t eventno; // Event Number.

TLorentzVector el, el_corZ, el_new;
TLorentzVector pip;
TLorentzVector pim;
TLorentzVector pr;
TLorentzVector miss_el;
TLorentzVector miss_all;
TLorentzVector miss_pr, miss_pr_corrected;
vector<double> miss_pr_m;
TLorentzVector beam(0,0,10.6,10.6);
TLorentzVector target(0,0,0,0.938272); // Set 4 vector for target, stationary so no momentum
// Assign a branch to each measurable and name it
TChain fake("hipo");
TString FileName = "/lustre19/expphy/cache/clas12/rg-a/production/recon/fall2018/torus+1/pass2/dst/recon/0054" ;//, "/lustre19/expphy/cache/clas12/rg-a/production/recon/fall2018/torus+1/pass2/dst/recon/005425", "/lustre19/expphy/cache/clas12/rg-a/production/recon/fall2018/torus+1/pass2/dst/recon/005428", "/lustre19/expphy/cache/clas12/rg-a/production/recon/fall2018/torus+1/pass2/dst/recon/005430", "/lustre19/expphy/cache/clas12/rg-a/production/recon/fall2018/torus+1/pass2/dst/recon/005432", "/lustre19/expphy/cache/clas12/rg-a/production/recon/fall2018/torus+1/pass2/dst/recon/005434", "/lustre19/expphy/cache/clas12/rg-a/production/recon/fall2018/torus+1/pass2/dst/recon/005436", "/lustre19/expphy/cache/clas12/rg-a/production/recon/fall2018/torus+1/pass2/dst/recon/005440", "/lustre19/expphy/cache/clas12/rg-a/production/recon/fall2018/torus+1/pass2/dst/recon/005441", "/lustre19/expphy/cache/clas12/rg-a/production/recon/fall2018/torus+1/pass2/dst/recon/005442"};
//Telling which files to run over
for(int j=23; j<28; j++){
TString inputFile(FileName + Form("%i", j) + "/*.hipo");
std::cout<<inputFile<<endl;
//Creating a chain for the data from different files
fake.Add(inputFile.Data());
//get the hipo data
  
}
  auto files=fake.GetListOfFiles();
  auto db=TDatabasePDG::Instance();
  Int_t el_status_min[2] = {0,2000};
  Int_t el_status_max[2] = {1999,3999};

  Double_t x, y, z, r;

  bool inbending = true;
  bool outbending = false;

  Int_t part_pid;

  Double_t part_DC_c1x,part_DC_c1y,part_DC_c1z;
  Double_t part_DC_c2x,part_DC_c2y,part_DC_c2z;
  Double_t part_DC_c3x,part_DC_c3y,part_DC_c3z;

  TFile tree_file("hipo_tree_mmpr.root","recreate");
  TTree *Tree = new TTree("Tree","tree"); 
  Tree->Branch("eventno",&eventno);
  Tree->Branch("miss_pr_m",&miss_pr_m);
  //Going over all files listed above and all events inside those
  for(Int_t i=0;i<files->GetEntries();i++){
    auto start = std::chrono::high_resolution_clock::now();
    gBenchmark->Start("timer");
    int counter=0;
    //create the event reader
    clas12reader c12(files->At(i)->GetTitle());
    cout<<"File in loop: "<<files->At(i)->GetTitle()<<endl;


    c12.useFTBased(); //and use the Pids from RECFT

    // c12.addExactPid(11,1);    //exactly 1 electron
    // c12.addExactPid(211,1);    //exactly 1 pi+
    // c12.addExactPid(-211,1);    //exactly 1 pi-
    // c12.addExactPid(2212,1);    //exactly 1 proton

    //loop over all events in the file
    while(c12.next()==true){
      counter ++;
      // Clear the vectors from the previous event
      miss_pr_m.clear();
      // Define how to access the information for each event
      eventno = c12.runconfig()->getEvent(); // Getting the event number
      if(c12.getDetParticles().empty())
      	continue;
      

      auto electron=c12.getByID(11);
      auto proton=c12.getByID(2212);
      auto piplus=c12.getByID(211);
      auto piminus=c12.getByID(-211);

      if (electron.size()!=1  || proton.size()!=1 || piplus.size()!=1 || piminus.size()!=1
          || abs(electron[0]->par()->getStatus())>1999 //getStatus related to the status of the particle, where the particle has been detected. Electron only forward tagger.
          || abs(piplus[0]->par()->getStatus())<2000 || abs(piplus[0]->par()->getStatus())>3999 // only forward detector 
          || abs(piminus[0]->par()->getStatus())<2000 || abs(piminus[0]->par()->getStatus())>3999
          || abs(proton[0]->par()->getStatus())<2000 || abs(proton[0]->par()->getStatus())>3999) continue;

      el.SetXYZM(electron[0]->par()->getPx(),electron[0]->par()->getPy(),electron[0]->par()->getPz(),0.000511);
      pip.SetXYZM(piplus[0]->par()->getPx(),piplus[0]->par()->getPy(),piplus[0]->par()->getPz(),0.139600);
      pim.SetXYZM(piminus[0]->par()->getPx(),piminus[0]->par()->getPy(),piminus[0]->par()->getPz(),0.139600);
      pr.SetXYZM(proton[0]->par()->getPx(),proton[0]->par()->getPy(),proton[0]->par()->getPz(),0.938300);
      //cout<<"Ger1c"<<endl;

      //pip
      if(piplus[0]->traj(DC,6)->getDetector() == 6 && piplus[0]->traj(DC,6)->getLayer() == 6){ //events in all 3 layers of drift chamber, for the fiducial cuts.
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

      part_pid = 0;
      // Positive particles asssigned to pi^+, negative to pi^-

      part_pid = 211;

      int part_DC_sector_pip = determineSector(part_DC_c2x, part_DC_c2y,part_DC_c2z);

      // Use this cut if looking at inbending data
      // if(!DC_fiducial_cut_theta_phi(part_DC_c1x,part_DC_c1y,part_DC_c1z,
      //   part_DC_c2x, part_DC_c2y, part_DC_c2z,
      //   part_DC_c3x,part_DC_c3y, part_DC_c3z,
      //   part_pid,part_DC_sector,1) ||
      //   !DC_fiducial_cut_theta_phi(part_DC_c1x,part_DC_c1y,part_DC_c1z,
      //     part_DC_c2x, part_DC_c2y, part_DC_c2z,
      //     part_DC_c3x,part_DC_c3y, part_DC_c3z,
      //     part_pid,part_DC_sector,2) ||
      //     !DC_fiducial_cut_theta_phi(part_DC_c1x,part_DC_c1y,part_DC_c1z,
      //       part_DC_c2x, part_DC_c2y, part_DC_c2z,
      //       part_DC_c3x, part_DC_c3y, part_DC_c3z,
      //       part_pid,part_DC_sector,3))continue;


      // Use this cut if looking at outbending data
      if(!DC_fiducial_cut_XY(part_DC_c1x,part_DC_c1y,part_DC_c1z,
        part_DC_c2x, part_DC_c2y, part_DC_c2z,
        part_DC_c3x, part_DC_c3y, part_DC_c3z, part_pid,part_DC_sector_pip,1) &&
        !DC_fiducial_cut_XY(part_DC_c1x,part_DC_c1y,part_DC_c1z,
          part_DC_c2x, part_DC_c2y, part_DC_c2z,
          part_DC_c3x, part_DC_c3y, part_DC_c3z, part_pid,part_DC_sector_pip,2) &&
          !DC_fiducial_cut_XY(part_DC_c1x,part_DC_c1y,part_DC_c1z,
            part_DC_c2x, part_DC_c2y, part_DC_c2z,
            part_DC_c3x, part_DC_c3y, part_DC_c3z, part_pid,part_DC_sector_pip,3)) continue;
      //pip
//////////////////////////////////
      //pim
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

      part_pid = 0;
      // Positive particles asssigned to pi^+, negative to pi^-

      part_pid = -211;
      //determineSector->The sector of charged particles is extracted from the hit position in the local phi coordinate.
      int part_DC_sector_pim = determineSector(part_DC_c2x, part_DC_c2y,part_DC_c2z);

      // Use this cut if looking at inbending data
      // if(!DC_fiducial_cut_theta_phi(part_DC_c1x,part_DC_c1y,part_DC_c1z,
      //   part_DC_c2x, part_DC_c2y, part_DC_c2z,
      //   part_DC_c3x,part_DC_c3y, part_DC_c3z,
      //   part_pid,part_DC_sector,1) ||
      //   !DC_fiducial_cut_theta_phi(part_DC_c1x,part_DC_c1y,part_DC_c1z,
      //     part_DC_c2x, part_DC_c2y, part_DC_c2z,
      //     part_DC_c3x,part_DC_c3y, part_DC_c3z,
      //     part_pid,part_DC_sector,2) ||
      //     !DC_fiducial_cut_theta_phi(part_DC_c1x,part_DC_c1y,part_DC_c1z,
      //       part_DC_c2x, part_DC_c2y, part_DC_c2z,
      //       part_DC_c3x, part_DC_c3y, part_DC_c3z,
      //       part_pid,part_DC_sector,3))continue;


      // Use this cut if looking at outbending data
      if(!DC_fiducial_cut_XY(part_DC_c1x,part_DC_c1y,part_DC_c1z,
        part_DC_c2x, part_DC_c2y, part_DC_c2z,
        part_DC_c3x, part_DC_c3y, part_DC_c3z, part_pid,part_DC_sector_pim,1) &&
        !DC_fiducial_cut_XY(part_DC_c1x,part_DC_c1y,part_DC_c1z,
          part_DC_c2x, part_DC_c2y, part_DC_c2z,
          part_DC_c3x, part_DC_c3y, part_DC_c3z, part_pid,part_DC_sector_pim,2) &&
          !DC_fiducial_cut_XY(part_DC_c1x,part_DC_c1y,part_DC_c1z,
            part_DC_c2x, part_DC_c2y, part_DC_c2z,
            part_DC_c3x, part_DC_c3y, part_DC_c3z, part_pid,part_DC_sector_pim,3)) continue;
      //pim
      //////////////////////////////////
      //pr
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

      part_pid = 0;
      // Positive particles asssigned to pi^+, negative to pi^-

      part_pid = 2212;

      int part_DC_sector_pr = determineSector(part_DC_c2x, part_DC_c2y,part_DC_c2z);

      // Use this cut if looking at inbending data
      // if(!DC_fiducial_cut_theta_phi(part_DC_c1x,part_DC_c1y,part_DC_c1z,
      //   part_DC_c2x, part_DC_c2y, part_DC_c2z,
      //   part_DC_c3x,part_DC_c3y, part_DC_c3z,
      //   part_pid,part_DC_sector,1) ||
      //   !DC_fiducial_cut_theta_phi(part_DC_c1x,part_DC_c1y,part_DC_c1z,
      //     part_DC_c2x, part_DC_c2y, part_DC_c2z,
      //     part_DC_c3x,part_DC_c3y, part_DC_c3z,
      //     part_pid,part_DC_sector,2) ||
      //     !DC_fiducial_cut_theta_phi(part_DC_c1x,part_DC_c1y,part_DC_c1z,
      //       part_DC_c2x, part_DC_c2y, part_DC_c2z,
      //       part_DC_c3x, part_DC_c3y, part_DC_c3z,
      //       part_pid,part_DC_sector,3))continue;


      // Use this cut if looking at outbending data
      if(!DC_fiducial_cut_XY(part_DC_c1x,part_DC_c1y,part_DC_c1z,
        part_DC_c2x, part_DC_c2y, part_DC_c2z,
        part_DC_c3x, part_DC_c3y, part_DC_c3z, part_pid,part_DC_sector_pr,1) &&
        !DC_fiducial_cut_XY(part_DC_c1x,part_DC_c1y,part_DC_c1z,
          part_DC_c2x, part_DC_c2y, part_DC_c2z,
          part_DC_c3x, part_DC_c3y, part_DC_c3z, part_pid,part_DC_sector_pr,2) &&
          !DC_fiducial_cut_XY(part_DC_c1x,part_DC_c1y,part_DC_c1z,
            part_DC_c2x, part_DC_c2y, part_DC_c2z,
            part_DC_c3x, part_DC_c3y, part_DC_c3z, part_pid,part_DC_sector_pr,3)) continue;
      //pr

      if(TMath::RadToDeg()*el.Theta() > 2.5 && TMath::RadToDeg()*el.Theta() < 4.5){

        el_new = Correct_Electron(el);

        // If no final particles were missing. We would expect the invariant mass of miss_all to be zero, no?
        miss_all = beam + target - el - pip - pim - pr;
        // If the proton was missing from the final state particles. With the same logic, we would expect the invariant mass of the
        // system to be the proton mass.
        miss_pr = beam + target - el - pip - pim;
        // I don't really understand this cut. Why 0.009? 3 sigma of the width of the delta function. 
        if(miss_all.M2()>=-0.0102 && miss_all.M2()<=0.0102){
          miss_pr_m.push_back(miss_pr.M());
          // He probably fits a Gaussian to this from the terminal.
          //hmiss_pr->Fill(miss_pr.M());
          Tree->Fill();
          Tree_Events++;
        }

        // miss_all = beam + target - el_new - pip - pim - pr;
        //
        // miss_pr = beam + target - el_new - pip - pim;
        //
        // if(miss_all.M2()>=-0.006 && miss_all.M2()<=0.006){
        //   hmiss_pr_corrected->Fill(miss_pr.M());
        // }
    }

  }
auto finish = std::chrono::high_resolution_clock::now();
std::chrono::duration<double> elapsed = finish - start;
}
Tree->Write(); // Write information to the TTree
cout<<"events in tree "<<Tree_Events<<endl;
//fileOutput1.Write();  

}
