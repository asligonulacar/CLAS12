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

void Fit_seprate_bins_RGK_125(){

  TFile *f1=new TFile("/Users/asligonulacar/Desktop/PhD/CLAS/Very_Strange/HipoRead/GeraintsCode/draw2d.root");

  TFile fileOutput1("fitBins.root","recreate");

  auto* hel_d_e_VS_elE_mean = new TH1F("hel_d_e_VS_elE_mean","Mean of Difference Between 'Missing' e' Energy and Detected e' Energy VS Detected e' Energy in FT;#DeltaE(e) [GeV];Counts",125,0,125);
  hel_d_e_VS_elE_mean->GetXaxis()->SetLabelSize(0.05);
  hel_d_e_VS_elE_mean->GetYaxis()->SetLabelSize(0.05);

  Double_t el_d_e_VS_elE_mean[125], el_d_e_VS_elE_mean_error[125], E[125];

  //Delta E Vs E
  auto *hel_d_e_VS_elE  = (TH2F*) f1->Get("hel_d_e_VS_elE");


  hel_d_e_VS_elE->Rebin2D(8,1);

  //hel_d_e_VS_elE->GetXaxis()->SetRangeUser(2, 7);
  //hel_d_e_VS_elE->GetYaxis()->SetRangeUser(-0.1, 0.6);

  TH1D** hel_d_e_VS_elE_bins = new TH1D*[125];

  for(Int_t i=0; i<125; i++){
    hel_d_e_VS_elE_bins[i] = hel_d_e_VS_elE->ProjectionY(Form("hel_d_e_VS_elE_bins_%d",i),i,i+1);

    TF1 *func1 = new TF1("func1","gaus(0)",hel_d_e_VS_elE_bins[i]->GetMean()-1*hel_d_e_VS_elE_bins[i]->GetRMS(), hel_d_e_VS_elE_bins[i]->GetMean()+1*hel_d_e_VS_elE_bins[i]->GetRMS());

    if(hel_d_e_VS_elE_bins[i]->GetEntries() > 1){

      hel_d_e_VS_elE_bins[i]->Fit("func1","R","B");

      hel_d_e_VS_elE_mean->SetBinContent(i,func1->GetParameter(1));

      el_d_e_VS_elE_mean[i] = func1->GetParameter(1);
      el_d_e_VS_elE_mean_error[i] = func1->GetParError(1);

    }

    /*else{hel_d_e_VS_elE_mean->SetBinContent(i,0);

      el_d_e_VS_elE_mean[i] = 0;
      el_d_e_VS_elE_mean_error[i] = 0;

    }*/

    E[i] = 0.08*i;

  }

  auto dE_VS_E_graph = new TGraphErrors(125, E, el_d_e_VS_elE_mean, 0, el_d_e_VS_elE_mean_error);

  dE_VS_E_graph->Draw("AP");
  gPad->Modified();
  dE_VS_E_graph->GetXaxis()->SetLimits(0,5);
  dE_VS_E_graph->GetYaxis()->SetLimits(-1,1);
  fileOutput1.Write();

}
