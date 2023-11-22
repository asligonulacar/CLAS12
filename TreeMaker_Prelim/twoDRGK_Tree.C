#include "Fiducial_cuts.cxx"
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
#include <TMath.h>
#include <TBenchmark.h>
#include "clas12reader.h"

using namespace clas12;

TLorentzVector  Correct_Electron(TLorentzVector x){

  Double_t E_new, Px_el, Py_el, Pz_el;
  TLorentzVector el_new;

  E_new = -1.412+2.415*x.E()-0.5048*pow(x.E(),2)+0.0805*pow(x.E(),3)-0.0048*pow(x.E(),4);

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

int getComponent(Double_t x, Double_t y) {

  int nCrystalX = 22;

  double crystal_size = 1.5;

  int ix = (int) ((x+11*crystal_size)/crystal_size);
  int iy = (int) ((y+11*crystal_size)/crystal_size);

  return iy*nCrystalX+ix;
}

void twoDRGK_Tree(){
  Int_t Tree_Events=0; // Records the number of events that are put in the tree
  Int_t eventno; // Event Number.
  


  TLorentzVector beam(0,0,10.6,10.6);
  TLorentzVector target(0,0,0,0.938300);

  TLorentzVector el;
  Double_t E_new, Px_el, Py_el, Pz_el;
  TLorentzVector el_new, el_corZ;
  TLorentzVector pip;
  TLorentzVector pim;
  TLorentzVector pr;
  TLorentzVector miss_el;
  TLorentzVector miss_all;
  TLorentzVector miss_pr;
  Double_t Delta_E_el;
  Double_t Delta_TH_el;
  Double_t Delta_PH_el;
  Double_t deltaPhi, deltaTheta;

  vector<double> miss_el_E;
  vector<double> el_E;
  vector<double> el_Phi;
  vector<double> el_Theta;
  vector<double> Delta_E_El;
  vector<double> DeltaPhi;
  vector<double> DeltaTheta;
  vector<double> pr_Rho;
  vector<double> pr_Theta;
  vector<double> pr_Phi;
  vector<double> pip_Rho;
  vector<double> pip_Theta;
  vector<double> pip_Phi;
  vector<double> pim_Rho;
  vector<double> pim_Theta;
  vector<double> pim_Phi;
  vector<double> X;
  vector<double> Y;
  vector<double> XYcomp;

  Double_t x, y, z, r;

  bool inbending = true;
  bool outbending = false;

  Int_t part_pid;

  Double_t part_DC_c1x,part_DC_c1y,part_DC_c1z;
  Double_t part_DC_c2x,part_DC_c2y,part_DC_c2z;
  Double_t part_DC_c3x,part_DC_c3y,part_DC_c3z;

  Int_t myeventnumber;

  TChain fake("hipo");
  TFile fileOutput1("S_w_RGK_all_75_2d_newCal_12_10_21.root","recreate");

  ofstream myfile;
  myfile.open ("Momentum_e_pPipPim_RGK_FTel.txt");
  TString FileName = "/lustre19/expphy/cache/clas12/rg-a/production/recon/fall2018/torus+1/pass2/dst/recon/0054" ;
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
  TFile tree_file("hipo_tree_2d.root","recreate");
  TTree *Tree = new TTree("Tree","tree"); 
  Tree->Branch("eventno",&eventno);
  Tree->Branch("miss_el_E",&miss_el_E);
  Tree->Branch("el_E",&el_E);
  Tree->Branch("el_Phi",&el_Phi);
  Tree->Branch("el_Theta",&el_Theta);
  Tree->Branch("Delta_E_El",& Delta_E_El);
  Tree->Branch("DeltaPhi",&DeltaPhi);
  Tree->Branch("DeltaTheta",&DeltaTheta);
  Tree->Branch("pr_Rho",&pr_Rho);
  Tree->Branch("pr_Theta",&pr_Theta);
  Tree->Branch("pr_Phi",&pr_Phi);
  Tree->Branch("pip_Rho",&pip_Rho);
  Tree->Branch("pip_Theta",&pip_Theta);
  Tree->Branch("pip_Phi",&pip_Phi);
  Tree->Branch("eventno",&eventno);
  Tree->Branch("pim_Rho",&pip_Rho);
  Tree->Branch("pim_Theta",&pim_Theta);
  Tree->Branch("pim_Phi",&pim_Phi);
  Tree->Branch("X",&X);
  Tree->Branch("Y",&Y);
  Tree->Branch("XYcomp", &XYcomp);
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
    myeventnumber=0;
    //open the file
    while(c12.next()==true){
      counter ++;
      miss_el_E.clear();
      el_E.clear();
      el_Phi.clear();
      el_Theta.clear();
      Delta_E_El.clear();
      DeltaPhi.clear();
      DeltaTheta.clear();
      pr_Rho.clear();
      pr_Theta.clear();
      pr_Phi.clear();
      pip_Rho.clear();
      pip_Theta.clear();
      pip_Phi.clear();
      pim_Rho.clear();
      pim_Theta.clear();
      pim_Phi.clear();
      X.clear();
      Y.clear();
      XYcomp.clear();
      eventno = c12.runconfig()->getEvent(); // Getting the event number
      if(c12.getDetParticles().empty()) continue;

      myeventnumber++;

      auto electron=c12.getByID(11);
      auto proton=c12.getByID(2212);
      auto piplus=c12.getByID(211);
      auto piminus=c12.getByID(-211);


      if (electron.size()!=1  || proton.size()!=1 || piplus.size()!=1 || piminus.size()!=1
          || abs(electron[0]->par()->getStatus())>1999
          || abs(piplus[0]->par()->getStatus())<2000 || abs(piplus[0]->par()->getStatus())>3999
          || abs(piminus[0]->par()->getStatus())<2000 || abs(piminus[0]->par()->getStatus())>3999
          || abs(proton[0]->par()->getStatus())<2000 || abs(proton[0]->par()->getStatus())>3999) continue;

      el.SetXYZM(electron[0]->par()->getPx(),electron[0]->par()->getPy(),electron[0]->par()->getPz(),0.000511);
      pip.SetXYZM(piplus[0]->par()->getPx(),piplus[0]->par()->getPy(),piplus[0]->par()->getPz(),0.139600);
      pim.SetXYZM(piminus[0]->par()->getPx(),piminus[0]->par()->getPy(),piminus[0]->par()->getPz(),0.139600);
      pr.SetXYZM(proton[0]->par()->getPx(),proton[0]->par()->getPy(),proton[0]->par()->getPz(),0.938300);
      // Getting FTCAL positions of electron?
      x = electron[0]->ft(FTCAL)->getX();
      y = electron[0]->ft(FTCAL)->getY();
      z = electron[0]->ft(FTCAL)->getZ();
//CReate 2 two-dimensional histograms for y vs x
//One histograms is filled with (x,y, Deltae) A
//The second histograms is filled with (x,y) B
//Once the loop over all vevents is done, you take the ratio of the two (A/B `and you get the average Deltae  function of y vs x



      //pip
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

      el_new = Correct_Electron(el);

      miss_el = beam + target - pip - pim - pr;

      miss_all = beam + target - el - pip - pim - pr;
      miss_pr = beam + target - el - pip - pim;
      // What are these cuts? 
      if(miss_all.M2()<-0.0102 || miss_all.M2()>0.0102 || miss_pr.M()<0.73272 || miss_pr.M()>1.2536278) continue;

      if(TMath::RadToDeg()*el.Theta() > 2.5 && TMath::RadToDeg()*el.Theta() < 4.5){

        Delta_E_el = miss_el.E() - el.E();
        deltaPhi = (miss_el.Phi() - el.Phi())*TMath::RadToDeg();
        deltaTheta = (miss_el.Theta() - el.Theta())*TMath::RadToDeg();

        Delta_E_El.push_back(Delta_E_el);
        DeltaTheta.push_back(deltaTheta);
        DeltaPhi.push_back(deltaPhi);
        el_E.push_back(el.E());
        el_Phi.push_back(el.Phi());
        miss_el_E.push_back(miss_el.E());
        el_Theta.push_back(el.Theta());
        pr_Rho.push_back(pr.Rho());
        pr_Theta.push_back(pr.Theta());
        pr_Phi.push_back(pr.Phi());
        pip_Rho.push_back(pip.Rho());
        pip_Theta.push_back(pip.Theta());
        pip_Phi.push_back(pip.Phi());
        pim_Rho.push_back(pim.Rho());
        pim_Theta.push_back(pim.Theta());
        pim_Theta.push_back(pim.Theta());
        XYcomp.push_back(getComponent(x,y));

        //hel_d_e_VS_elE->Fill(el.E(),Delta_E_el);
        //hel_d_e_VS_missE->Fill(miss_el.E(),Delta_E_el);
        //hel_delE_VS_elPH->Fill(el.Phi()*TMath::RadToDeg(),Delta_E_el);
        //hel_delE_VS_elTH->Fill(el.Theta()*TMath::RadToDeg(),Delta_E_el);
        //hel_d_e_VS_crysN->Fill(getComponent(x,y),Delta_E_el);
        //hel_delPH_VS_elE->Fill(el.E(),deltaPhi);
        //hel_delTH_VS_elE->Fill(el.E(),deltaTheta);
        //hel_delPH_VS_elPH->Fill(el.Phi()*TMath::RadToDeg(), deltaPhi);
        //hel_delPH_VS_elTH->Fill(el.Theta()*TMath::RadToDeg(), deltaPhi);
        //hel_delTH_VS_elTH->Fill(el.Theta()*TMath::RadToDeg(), deltaTheta);
        //hel_delTH_VS_elPH->Fill(el.Phi()*TMath::RadToDeg(), deltaTheta);
        //hprTH_VS_prMom->Fill(pr.Rho(),pr.Theta()*TMath::RadToDeg());
        //hprPH_VS_prMom->Fill(pr.Rho(),pr.Phi()*TMath::RadToDeg());
        //hpipTH_VS_pipMom->Fill(pip.Rho(),pip.Theta()*TMath::RadToDeg());
        //hpipPH_VS_pipMom->Fill(pip.Rho(),pip.Phi()*TMath::RadToDeg());
        //hpimTH_VS_pimMom->Fill(pim.Rho(),pim.Theta()*TMath::RadToDeg());
        //hpimPH_VS_pimMom->Fill(pim.Rho(),pim.Phi()*TMath::RadToDeg());


        myfile<<myeventnumber<<endl;
        myfile<<11<<" "<<el.Px()<<" "<<el.Py()<<" "<<el.Pz()<<endl;
        myfile<<2212<<" "<<pr.Px()<<" "<<pr.Py()<<" "<<pr.Pz()<<endl;
        myfile<<211<<" "<<pip.Px()<<" "<<pip.Py()<<" "<<pip.Pz()<<endl;
        myfile<<-211<<" "<<pim.Px()<<" "<<pim.Py()<<" "<<pim.Pz()<<endl;
        //myfile<<"\n"<<endl;

      }
      X.push_back(x);
      Y.push_back(y);
      Tree->Fill();
      Tree_Events++;
      //hX_Y_dE->Fill(x,y,Delta_E_el);
      //hX_Y->Fill(x,y);

  }
auto finish = std::chrono::high_resolution_clock::now();
std::chrono::duration<double> elapsed = finish - start;
std::cout << "Elapsed time: " << elapsed.count()<< " events = "<<counter<< " s\n";
}
Tree->Write(); // Write information to the TTree
cout<<"events in tree "<<Tree_Events<<endl;
//fileOutput1.Write();  
//hX_Y_ratio->Divide(hX_Y_dE, hX_Y);
}
