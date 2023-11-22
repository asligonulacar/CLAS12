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
// Alter this so it gets the histograms from the tree.
void Draw2d(){
  // Read Tree
  TFile fout("draw2d.root","recreate");
  TFile *f = TFile::Open("hipo_tree_2d_2.root");
  TTree *Tree = (TTree*) f->Get("Tree;1");
  // Set Pointers for Branch Addresses 
  vector<double> *miss_el_E = 0;
  vector<double> *el_E = 0;
  vector<double> *el_Phi = 0;
  vector<double> *el_Theta = 0;
  vector<double> *Delta_E_El = 0;
  vector<double> *DeltaPhi = 0;
  vector<double> *DeltaTheta = 0;
  vector<double> *pr_Rho = 0;
  vector<double> *pr_Theta = 0;
  vector<double> *pr_Phi = 0;
  vector<double> *pip_Rho = 0;
  vector<double> *pip_Theta = 0;
  vector<double> *pip_Phi = 0;
  vector<double> *pim_Rho = 0;
  vector<double> *pim_Theta = 0;
  vector<double> *pim_Phi = 0;
  vector<double> *X = 0;
  vector<double> *Y = 0;
  vector<double> *XYcomp = 0;
  // Set Branch Addresses
  Tree->SetBranchAddress("miss_el_E",&miss_el_E);
  Tree->SetBranchAddress("el_E",&el_E);
  Tree->SetBranchAddress("el_Phi",&el_Phi);
  Tree->SetBranchAddress("el_Theta",&el_Theta);
  Tree->SetBranchAddress("Delta_E_El",&Delta_E_El);
  Tree->SetBranchAddress("DeltaPhi",&DeltaPhi);
  Tree->SetBranchAddress("DeltaTheta",&DeltaTheta);
  Tree->SetBranchAddress("pr_Rho",&pr_Rho);
  Tree->SetBranchAddress("pr_Theta",&pr_Theta);
  Tree->SetBranchAddress("pr_Phi",&pr_Phi);
  Tree->SetBranchAddress("pip_Rho",&pip_Rho);
  Tree->SetBranchAddress("pip_Theta",&pip_Theta);
  Tree->SetBranchAddress("pip_Phi",&pip_Phi);
  Tree->SetBranchAddress("pim_Rho",&pip_Rho);
  Tree->SetBranchAddress("pim_Theta",&pim_Theta);
  Tree->SetBranchAddress("pim_Phi",&pim_Phi);
  Tree->SetBranchAddress("X",&X);
  Tree->SetBranchAddress("Y",&Y);
  Tree->SetBranchAddress("XYcomp", &XYcomp);
  // Initialise histograms
  auto* hel_d_e_VS_elE=new TH2D("hel_d_e_VS_elE","Difference Between 'Missing' e' Energy and Detected e' Energy VS Detected e' Energy in FT;E [GeV];#DeltaE(e) [GeV]",100,0,10,200,-0.4,1);
  hel_d_e_VS_elE->GetXaxis()->SetLabelSize(0.05);
  hel_d_e_VS_elE->GetYaxis()->SetLabelSize(0.05);

  auto* hel_d_e_VS_missE=new TH2D("hel_d_e_VS_missE","Difference Between 'Missing' e' Energy and Detected e' Energy VS Missing e' Energy in FT;E [GeV];#DeltaE(e) [GeV]",1000,0,10,2000,-10,10);
  hel_d_e_VS_missE->GetXaxis()->SetLabelSize(0.05);
  hel_d_e_VS_missE->GetYaxis()->SetLabelSize(0.05);

  auto* hel_d_e_VS_crysN=new TH2D("hel_d_e_VS_crysN","Difference Between 'Missing' e' Energy and Detected e' Energy VS Crystsl Number;N_{crystal};#DeltaE(e) [GeV]",500,0,500,2000,-10,10);
  hel_d_e_VS_crysN->GetXaxis()->SetLabelSize(0.05);
  hel_d_e_VS_crysN->GetYaxis()->SetLabelSize(0.05);

  auto* hel_delE_VS_elPH=new TH2D("hel_delE_VS_elPH","Difference Between 'Missing' e' Energy and Detected e' Energy VS Detected e' #phi;#phi(e) [deg];#DeltaE(e) [GeV]",7200,-360,360,2000,-10,10);
  hel_delE_VS_elPH->GetXaxis()->SetLabelSize(0.05);
  hel_delE_VS_elPH->GetYaxis()->SetLabelSize(0.05);

  auto* hel_delE_VS_elTH=new TH2D("hel_delE_VS_elTH","Difference Between 'Missing' e' Energy and Detected e' Energy VS Detected e' #theta;#theta(e) [deg];#DeltaE(e) [GeV]",100,0,10,2000,-10,10);
  hel_delE_VS_elTH->GetXaxis()->SetLabelSize(0.05);
  hel_delE_VS_elTH->GetYaxis()->SetLabelSize(0.05);

  auto* hel_delPH_VS_elE=new TH2D("hel_delPH_VS_elE","Difference Between 'Missing' e' #Phi and Detected e' #Phi VS Detected e' Energy in FT;E [GeV];#Delta#phi(e) [deg]",1000,0,10,7200,-180,180);
  hel_delPH_VS_elE->GetXaxis()->SetLabelSize(0.05);
  hel_delPH_VS_elE->GetYaxis()->SetLabelSize(0.05);

  auto* hel_delTH_VS_elE=new TH2D("hel_delTH_VS_elE","Difference Between 'Missing' e' #Theta and Detected e' #Theta VS Detected e' Energy in FT;E [GeV];#Delta#theta(e) [deg]",1000,0,10,7200,-360,360);
  hel_delTH_VS_elE->GetXaxis()->SetLabelSize(0.05);
  hel_delTH_VS_elE->GetYaxis()->SetLabelSize(0.05);

  auto* hel_delPH_VS_elPH=new TH2D("hel_delPH_VS_elPH","Difference Between 'Missing' e' #phi and Detected e' #phi VS Detected e' #phi in FT;#phi [deg];#Delta#phi(e) [deg]",7200,-180,180,7200,-180,180);
  hel_delPH_VS_elPH->GetXaxis()->SetLabelSize(0.05);
  hel_delPH_VS_elPH->GetYaxis()->SetLabelSize(0.05);

  auto* hel_delPH_VS_elTH=new TH2D("hel_delPH_VS_elTH","Difference Between 'Missing' e' #phi and Detected e' #phi VS Detected e' #theta in FT;#theta [deg];#Delta#phi(e) [deg]",7200,-360,360,7200,-180,180);
  hel_delPH_VS_elTH->GetXaxis()->SetLabelSize(0.05);
  hel_delPH_VS_elTH->GetYaxis()->SetLabelSize(0.05);

  auto* hel_delTH_VS_elTH=new TH2D("hel_delTH_VS_elTH","Difference Between 'Missing' e' #theta and Detected e' #theta VS Detected e' #theta in FT;#theta [deg];#Delta#theta(e) [deg]",7200,-360,360,7200,-360,360);
  hel_delTH_VS_elTH->GetXaxis()->SetLabelSize(0.05);
  hel_delTH_VS_elTH->GetYaxis()->SetLabelSize(0.05);

  auto* hel_delTH_VS_elPH=new TH2D("hel_delTH_VS_elPH","Difference Between 'Missing' e' #theta and Detected e' #theta VS Detected e' #phi in FT;#phi [deg];#Delta#theta(e) [deg]",7200,-180,180,7200,-360,360);
  hel_delTH_VS_elPH->GetXaxis()->SetLabelSize(0.05);
  hel_delTH_VS_elPH->GetYaxis()->SetLabelSize(0.05);

  auto* hprTH_VS_prMom=new TH2D("hprTH_VS_prMom","Proton #theta Vs Momentum;P_{p} [GeV];#theta_{p} [deg]",1000,0,10,1800,0,180);
  hprTH_VS_prMom->GetXaxis()->SetLabelSize(0.05);
  hprTH_VS_prMom->GetYaxis()->SetLabelSize(0.05);

  auto* hprPH_VS_prMom=new TH2D("hprPH_VS_prMom","Proton #phi Vs Momentum;P_{p} [GeV];#phi_{p} [deg]",1000,0,10,3600,-180,180);
  hprPH_VS_prMom->GetXaxis()->SetLabelSize(0.05);
  hprPH_VS_prMom->GetYaxis()->SetLabelSize(0.05);

  auto* hpipTH_VS_pipMom=new TH2D("hpipTH_VS_pipMom","#pi^{+} #theta Vs Momentum;P_{#pi^{+}} [GeV];#theta_{#pi^{+}} [deg]",1000,0,10,1800,0,180);
  hpipTH_VS_pipMom->GetXaxis()->SetLabelSize(0.05);
  hpipTH_VS_pipMom->GetYaxis()->SetLabelSize(0.05);

  auto* hpipPH_VS_pipMom=new TH2D("hpipPH_VS_pipMom","#pi^{+} #phi Vs Momentum;P_{#pi^{+}} [GeV];#phi_{#pi^{+}} [deg]",1000,0,10,3600,-180,180);
  hpipPH_VS_pipMom->GetXaxis()->SetLabelSize(0.05);
  hpipPH_VS_pipMom->GetYaxis()->SetLabelSize(0.05);

  auto* hpimTH_VS_pimMom=new TH2D("hpimTH_VS_pimMom","#pi^{-} #theta Vs Momentum;P_{#pi^{-}} [GeV];#theta_{#pi^{-}} [deg]",1000,0,10,1800,0,180);
  hpimTH_VS_pimMom->GetXaxis()->SetLabelSize(0.05);
  hpimTH_VS_pimMom->GetYaxis()->SetLabelSize(0.05);

  auto* hpimPH_VS_pimMom=new TH2D("hpimPH_VS_pimMom","#pi^{-} #phi Vs Momentum;P_{#pi^{-}} [GeV];#phi_{#pi^{-}} [deg]",1000,0,10,3600,-180,180);
  hpimPH_VS_pimMom->GetXaxis()->SetLabelSize(0.05);
  hpimPH_VS_pimMom->GetYaxis()->SetLabelSize(0.05);

  auto* hX_Y_dE=new TH2D("hX_Y_dE","y VS x Position of the FT Electron (With #DeltaE Weight);x [cm];y [cm]",4000,-20,20,4000,-20,20);
  hX_Y_dE->GetXaxis()->SetLabelSize(0.05);
  hX_Y_dE->GetYaxis()->SetLabelSize(0.05);

  auto* hX_Y=new TH2D("hX_Y","y VS x Position of the FT Electron (With #DeltaE Weight);x [cm];y [cm]",4000,-20,20,4000,-20,20);
  hX_Y->GetXaxis()->SetLabelSize(0.05);
  hX_Y->GetYaxis()->SetLabelSize(0.05);

  auto* hX_Y_ratio=new TH2D("hX_Y_ratio","y VS x Position of the FT Electron (With #DeltaE Weight);x [cm];y [cm]",4000,-20,20,4000,-20,20);
  hX_Y_ratio->GetXaxis()->SetLabelSize(0.05);
  hX_Y_ratio->GetYaxis()->SetLabelSize(0.05);
  // Loop over the entries in the tree
  Int_t nentries = (Int_t)Tree->GetEntries();
  std::cout<<nentries<<endl;
  for (Int_t i=0; i<nentries; i++) {
      Tree->GetEntry(i);
      // Problem loop: I have 1350 entries in the tree but I only have 1070 entries for some vectors. idk how.
        hel_d_e_VS_elE->Fill(el_E->at(0),Delta_E_El->at(0));
        hel_d_e_VS_missE->Fill(miss_el_E->at(0),Delta_E_El->at(0));
        hel_delE_VS_elPH->Fill(el_Phi->at(0)*TMath::RadToDeg(),Delta_E_El->at(0));
        hel_delE_VS_elTH->Fill(el_Theta->at(0)*TMath::RadToDeg(),Delta_E_El->at(0));
        hel_delPH_VS_elE->Fill(el_E->at(0),DeltaPhi->at(0));
        hel_delTH_VS_elE->Fill(el_E->at(0),DeltaTheta->at(0));
        hel_delPH_VS_elPH->Fill(el_Phi->at(0)*TMath::RadToDeg(), DeltaPhi->at(0));
        hel_delPH_VS_elTH->Fill(el_Theta->at(0)*TMath::RadToDeg(), DeltaPhi->at(0));
        hel_delTH_VS_elTH->Fill(el_Theta->at(0)*TMath::RadToDeg(), DeltaTheta->at(0));
        hel_delTH_VS_elPH->Fill(el_Phi->at(0)*TMath::RadToDeg(), DeltaTheta->at(0));
      }
      //
      /*
      hel_d_e_VS_crysN->Fill(XYcomp->at(0),Delta_E_El->at(0));
      
      hprTH_VS_prMom->Fill(pr_Rho->at(0),pr_Theta->at(0)*TMath::RadToDeg());
      hprPH_VS_prMom->Fill(pr_Rho->at(0),pr_Phi->at(0)*TMath::RadToDeg());
      hpipTH_VS_pipMom->Fill(pip_Rho->at(0),pip_Theta->at(0)*TMath::RadToDeg());
      hpipPH_VS_pipMom->Fill(pip_Rho->at(0),pip_Phi->at(0)*TMath::RadToDeg());
      //hpimTH_VS_pimMom->Fill(pim_Rho->at(i),pim_Theta->at(i)*TMath::RadToDeg());
      //hpimPH_VS_pimMom->Fill(pim_Rho->at(i),pim_Phi->at(i)*TMath::RadToDeg());
      */
  
  fout.cd();
  hel_d_e_VS_elE->Write();
  hel_d_e_VS_missE->Write();
  hel_delE_VS_elPH->Write();
  hel_delE_VS_elTH->Write();
  hel_delPH_VS_elE->Write();
  hel_delPH_VS_elPH->Write();
  hel_delTH_VS_elE->Write();
  hel_delPH_VS_elTH->Write();
  hel_delTH_VS_elTH->Write();
  hel_delTH_VS_elPH->Write();
  fout.Close();

}
