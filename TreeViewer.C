{
// Views Trees.
gROOT->Reset();
gROOT->SetStyle("Plain");

// Open file
TFile f("tree_v2.root");

// Get tree
TTree* tree = (TTree*)f.Get("tree;1");

tree->Show(0);

tree->Print();

tree->Scan();


}