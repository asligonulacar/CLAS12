enum ParIndex_t {p0, p1, p2, p3, Mean, Sigma, Scale, n_PAR};
// Use this map to (re-)name parameters for the plot
const std::map<ParIndex_t,std::string> parNames{{p0, "P0"},{p1, "P1"}, {p2, "P2"}, {p3, "P3"}, {Mean, "Mean"}, {Sigma, "Sigma"}, {Scale, "Scale"}};
// Polynomial background function.
Double_t background(Double_t *x, Double_t *par) {
    return par[p0] + par[p1]*x[0] + par[p2]*x[0]*x[0] + par[p3]*x[0]*x[0]*x[0];
}
Double_t sig(Double_t *x, Double_t *par) {
    return TMath::Gaus(x[0], par[Mean], par[Sigma], true)*par[Scale];
}
// Sum of background and peak function
Double_t fitFunction(Double_t *x, Double_t *par) {
      return background(x,par) + sig(x,par);
}
void FitHistogram(){
    TMinuit *gMinuit = new TMinuit(7);
    TFile *f = TFile::Open("5423_2000.root");
    //f->ls();
    TH1F *h1;
    f->GetObject("hmiss_all;1", h1);
    h1->SetAxisRange(-0.5, 0.5,"X");
    h1->SetAxisRange(0, 900,"Y");
    h1->GetXaxis()->SetTickLength(0.1);
    auto FitFNCTN = new TF1("FitFNCTN", fitFunction, -0.08, 0.04, n_PAR);
    for (auto& idx_name : parNames) {
     FitFNCTN->SetParName(idx_name.first, idx_name.second.c_str());
    }
    FitFNCTN->SetNpx(500);
    FitFNCTN->SetParameters(-234.743, -4000.24 ,14378.8, 161920 , -0.0009, 0.00344, 2);
    h1->Fit("FitFNCTN","R+", "ep");   
    TF1 *backFcn = new TF1("backFcn",background, -0.08, 0.04, n_PAR);
    //backFcn->SetLineColor(kGreen);
    TF1 *signalFcn = new TF1("signalFcn",sig,-0.003, 0.003, n_PAR);
    signalFcn->SetLineColor(kMagenta);
    signalFcn->SetNpx(700);
    double par[n_PAR];
    FitFNCTN->GetParameters(par);
    //backFcn->SetParameters(par);
    //backFcn->Draw("same");
    signalFcn->SetParameters(par);
    //h1->Fit(signalFcn, "VR+", "ep");
    //signalFcn->Draw("same");
    // draw the legend
}
