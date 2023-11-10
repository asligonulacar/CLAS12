#include "DC_Fiducial_Cuts_CLAS12.cxx"
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
void HipoRead(){
// Create the file to store the tree.
TFile tree_file("tree_v3.root","recreate");
// Create the tree
TTree *tree = new TTree("tree","tree");

//Telling which files to run over
TString inputFile("/u/home/gclash/12_10_21_RGK/*.hipo");

//Creating a chain for the data from different files
TChain fake("hipo");
fake.Add(inputFile.Data());

//get the hipo data
auto files=fake.GetListOfFiles();
auto db=TDatabasePDG::Instance();

TLorentzVector beam(0,0,10.6,10.6);
TLorentzVector target(0,0,0,db->GetParticle(2212)->Mass());
TLorentzVector el(0,0,0,db->GetParticle(11)->Mass());
TLorentzVector pr(0,0,0,db->GetParticle(2212)->Mass());
TLorentzVector pip(0,0,0,db->GetParticle(211)->Mass());
TLorentzVector pim(0,0,0,db->GetParticle(-211)->Mass());

//TLorentzVector miss_el; // scattered electron missing
//TLorentzVector miss_all; // all of the final state particles missing 
//TLorentzVector miss_pr; // proton missing 

//Going over all files listed above and all events inside those
for(Int_t i=0;i<files->GetEntries();i++){
    //create the event reader
    clas12reader c12(files->At(i)->GetTitle());
    cout<<"File in loop: "<<files->At(i)->GetTitle()<<endl;
    c12.useFTBased(); //and use the Pids from RECFT
    c12.addExactPid(11,1);    //exactly 1 electron
    c12.addExactPid(211,1);    //exactly 1 pi+
    c12.addExactPid(-211,1);    //exactly 1 pi-
    c12.addExactPid(2212,1);    //exactly 1 proton

    //loop over all events in the file
    myeventnumber=0;
    //open the file
    while(c12.next()==true){    
      if(c12.getDetParticles().empty()) continue;   
      myeventnumber++;
      auto electron=c12.getByID(11);
      auto proton=c12.getByID(2212);
      auto piplus=c12.getByID(211);
      auto piminus=c12.getByID(-211);

      el.SetXYZM(electron[0]->par()->getPx(),electron[0]->par()->getPy(),electron[0]->par()->getPz(),0.000511);
      pip.SetXYZM(piplus[0]->par()->getPx(),piplus[0]->par()->getPy(),piplus[0]->par()->getPz(),0.139600);
      pim.SetXYZM(piminus[0]->par()->getPx(),piminus[0]->par()->getPy(),piminus[0]->par()->getPz(),0.139600);
      pr.SetXYZM(proton[0]->par()->getPx(),proton[0]->par()->getPy(),proton[0]->par()->getPz(),0.938300);


    }
tree->Fill()
}
tree->Write();
}
