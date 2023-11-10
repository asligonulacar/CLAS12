{
  //Reads Trees.
  gROOT->Reset();
  gROOT->SetStyle("Plain");

  // Open file 
  TFile f("tree_v3.root");

  // Get tree
  TTree* tree = (TTree*)f.Get("tree;1");

  Float_t Px[5];
  Float_t Py[5];
  Float_t Pz[5];

  tree->SetBranchAddress("Px",&Px);
  tree->SetBranchAddress("Py",&Py);
  tree->SetBranchAddress("Pz",&Pz);
  
  tree->GetEvent(3);

  for(int i=0;i<5;i++){
      cout << Px[i] << endl;
  } 

}