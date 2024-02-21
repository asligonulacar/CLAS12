#include <cstdlib>
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TH1.h>
#include <TF1.h>
#include <TH2.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <vector>
#include <TLorentzVector.h>
#include <math.h>
#include <TMath.h>
#include <string>
#include <sstream>
#include <iostream>
#include <vector>
using namespace std;
 

// Now I want to plot the missing mass (obviously with more than 5 runs but I think I already have
// that data for the FD. I need to add trajectory information as a branch in my tree so I can do 
// a nice cut to select where the particle has been detected. What I want to do is sort the correct 
// momenta to the correct particles. The useFTbased(); command seems to be working, and I am actually 
// getting more electrons in the FT. 
// Get Region better than Get Status.
// Fix IGUANA stuff.

TLorentzVector  Correct_Electron(TLorentzVector x){
  
  Double_t E_new, Px_el, Py_el, Pz_el;
  TLorentzVector el_new;
   
  //E_new = x.E()-1.26573+0.865347*x.E()-0.201174*pow(x.E(),2)+0.0200699*pow(x.E(),3)-0.000724783*pow(x.E(),4);
  E_new = x.E() + 0.0208922 + 0.050158*x.E() - 0.0181107*pow(x.E(),2) + 0.00305671*pow(x.E(),3) - 0.000178235*pow(x.E(),4);
  Px_el = E_new*(x.Px()/x.Rho());
  Py_el = E_new*(x.Py()/x.Rho());
  Pz_el = E_new*(x.Pz()/x.Rho());

  el_new.SetXYZM(Px_el, Py_el, Pz_el, 0.000511);

  return el_new;
}
// Alter this so it gets the histograms from the tree.
void ana_2(){
  // Output
  TFile fileOutput1("elFTvsRECFT.root","recreate");
  TString Filename = "Trialpt2_0.root";

  TH1D *Pxhist = new TH1D("Pxhist", "Electron Theta; #theta ; Counts", 1000, -0.5, 0.5);
  TH1D *PxRECFThist = new TH1D("PxRECFThist", "Electron Theta; #theta ; Counts", 1000, -0.5, 0.5);

  TH2D *PxPxRECFT = new TH2D("PxPxRECFT", "momenta comparison", 1000, -9, 9, 1000, -9,9);
  TH2D *PyPyRECFT = new TH2D("PyPyRECFT", "momenta comparison", 1000, -9, 9, 1000, -9,9);
  TH2D *PzPzRECFT = new TH2D("PzPzRECFT", "momenta comparison", 1000, -9, 9, 1000, -9,9);
  TH2D *PidPidRECFT = new TH2D("PidPidRECFT", "momenta comparison", 1000, -3000, 3000, 1000, -3000, 3000);

  TLorentzVector beam(0,0,10.6,10.6);
  TLorentzVector target(0,0,0,0.938272);
  TLorentzVector el;
  TLorentzVector Kaon;
  TLorentzVector FourMomentum; 

  vector<double> *Pid_FT = 0;
  vector<double> *Px_FT = 0;
  vector<double>  PxFT;
  vector<double>  PyFT;
  vector<double>  PzFT;
  vector<double>  PidFT;
  vector<double> *Py_FT = 0;
  vector<double> *Pz_FT = 0;
  vector<double> *Beta_FT = 0;
  vector<double> *Chi2Pid_FT = 0;

  vector<double> *Pid_RECFT = 0;
  vector<double> *Px_RECFT = 0;
  vector<double>  PidRECFT;
  vector<double>  PxRECFT;
  vector<double>  PyRECFT;
  vector<double>  PzRECFT;
  vector<double> *Py_RECFT = 0;
  vector<double> *Pz_RECFT = 0;
  vector<double> *Beta_RECFT = 0;
  vector<double> *Chi2Pid_RECFT = 0;

  vector<TLorentzVector> *electrons = 0;
  vector<TLorentzVector> *Kaons = 0;
  //std::cout<<num_of_events<<endl;
  int event_count=0;
  int eventFT_count=0;
  int el_count=0;
  int elFD_count=0;
  int K_count=0;

  // Read Original Tree
  TFile *f = TFile::Open(Filename);
  // Get name of Tree.
  TTree *Tree_ft = (TTree*) f->Get(f->GetListOfKeys()->At(2)->GetName());
  TTree *Tree_Recft = (TTree*) f->Get(f->GetListOfKeys()->At(4)->GetName());

  vector <double>  a = {321, -321, 2212, -2212, 11, -211, 211};

  Tree_ft->SetBranchAddress("Pid_FT",&Pid_FT);
  Tree_ft->SetBranchAddress("Px_FT",&Px_FT);
  Tree_ft->SetBranchAddress("Py_FT",&Py_FT);
  Tree_ft->SetBranchAddress("Pz_FT",&Pz_FT);
  Tree_ft->SetBranchAddress("Beta_FT",&Beta_FT);
  Tree_ft->SetBranchAddress("Chi2Pid_FT",&Chi2Pid_FT);

  Tree_Recft->SetBranchAddress("Pid_RECFT",&Pid_RECFT);
  Tree_Recft->SetBranchAddress("Px_RECFT",&Px_RECFT);
  Tree_Recft->SetBranchAddress("Py_RECFT",&Py_RECFT);
  Tree_Recft->SetBranchAddress("Pz_RECFT",&Pz_RECFT);
  Tree_Recft->SetBranchAddress("Beta_RECFT",&Beta_RECFT);
  Tree_Recft->SetBranchAddress("Chi2Pid_RECFT",&Chi2Pid_RECFT);

  double num_of_eventsFTREC = Tree_Recft->GetEntries();
  double num_of_eventsFT = Tree_ft->GetEntries();

  //std::cout<< num_of_eventsFT << " " << num_of_eventsFTREC<<endl;
  // Loop over each event

  //Beta cut value
  double beta_cut_val = 20;

  for(int j=0; j<num_of_eventsFT; j++){
    //std::cout<< (j/num_of_events)*100 << '%'<<endl;
    // Increase event count by 1
    event_count++;

    // Set the loop to an event
    Tree_ft->GetEvent(j);

    for(int i=0; i<Px_FT->size(); i++){
            if(-beta_cut_val<Beta_FT->at(i)<beta_cut_val) {
            //std::cout<<Px->at(i)<<endl;
              //std::cout<<"check"<<endl;
              if ((std::find(a.begin(), a.end(), Pid_FT->at(i)) == a.end())) continue;
              //std::cout<<Pid_FT->at(i)<<endl;
              Pxhist->Fill(Px_FT->at(i));
              PxFT.push_back(Px_FT->at(i));
              PyFT.push_back(Py_FT->at(i));
              PzFT.push_back(Pz_FT->at(i));
              PidFT.push_back(Pid_FT->at(i));
          }  
        }
    eventFT_count++;
     // Set the loop to an event
    Tree_Recft->GetEvent(j);
    for(int k=0; k<Px_RECFT->size(); k++){
            //std::cout<<"check"<<endl;
            //std::cout<<Px_RECFT->at(k)<<endl;
            if(-beta_cut_val<Beta_RECFT->at(k)<beta_cut_val){
              if ((std::find(a.begin(), a.end(), Pid_RECFT->at(k)) == a.end())) continue;
              PxRECFThist->Fill(Px_RECFT->at(k));
              PxRECFT.push_back(Px_RECFT->at(k));
              PyRECFT.push_back(Py_RECFT->at(k));
              PzRECFT.push_back(Pz_RECFT->at(k));
              PidRECFT.push_back(Pid_RECFT->at(k));
              //if(!(std::find(a.begin(), a.end(), PidRECFT.at(k)) == a.end())) continue;
          }   
        } 

  // Loop over each event

    for(int l=0; l < PxFT.size(); l++){
      //std::cout<<"check"<<endl;
      if(PxFT.size()==PxRECFT.size()){

            std::sort(PidFT.begin(), PidFT.end());
            std::sort(PidRECFT.begin(), PidRECFT.end());
            
            if (PidFT==PidRECFT){
            //std::cout<<l<<endl;
            std::cout<<PxFT.at(l)/PxRECFT.at(l)<<endl;
            PxPxRECFT->Fill(PxFT.at(l), PxRECFT.at(l));
            PyPyRECFT->Fill(PyFT.at(l), PyRECFT.at(l));
            PzPzRECFT->Fill(PzFT.at(l), PzRECFT.at(l));
            PidPidRECFT->Fill(PidFT.at(l), PidRECFT.at(l));
            }
    }
}
    PidFT.clear();
    PxFT.clear();
    PyFT.clear();
    PzFT.clear();
    PidRECFT.clear();
    PxRECFT.clear();
    PyRECFT.clear();
    PzRECFT.clear();
    

}
fileOutput1.Write();  
}
// Plot all delta beta not just for one particle. Plot all for one mass see how particles move, plot for given mass see if clas12 assigns proper mass.

