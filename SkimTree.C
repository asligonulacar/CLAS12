{
// Fills Trees.
// Create the file to store the tree
TFile tree_file("tree_v3.root","recreate");

// Create the tree
TTree *tree = new TTree("tree","tree");

// Define the number of events
Int_t numofevn = 10;

// Define the number of final state particles 
const Int_t numofpart = 5;

// Create branches
Float_t Px[5];
Float_t Py[5];
Float_t Pz[5];
Float_t Beta[5];
Int_t Pid[5];

//auto branch_numofpart = tree.Branch("numofpart",&numofpart,"numofpart/I");
tree->Branch("Px",Px,"Px[5]/F");
tree->Branch("Py",Py,"Py[5]/F");
tree->Branch("Pz",Pz,"Pz[5]/F");
tree->Branch("Beta",Beta,"Beta[5]/F");
tree->Branch("Pid",Pid,"Pid[5]/I");


// Create loop for filling branches
TRandom3 rndgen;
Int_t i = 0;
for(Int_t i=0;i<numofevn;i++){

    Int_t j = 0;
    for(Int_t j=0;j<5;j++){
        Px[j]=rndgen.Gaus(-2,1);
        Py[j]=rndgen.Gaus(-2,1);
        Pz[j]=rndgen.Gaus(-2,1);
        Beta[j]=j*i;
        Pid[j]=100+i;           
    }

    tree->Fill(); 

}

tree->Write();

}
