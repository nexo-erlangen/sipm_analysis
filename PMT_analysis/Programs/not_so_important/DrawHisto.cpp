void DrawHisto(){
	unsigned int evt_n;
	double peak_integral;
	double peak_maximum;
	double time_maximum;
	//vector<double>  peakfinder_maximum;
	//vector<double>  peakfinder_max_time;
	double baseline_mean;
	double baseline_rms;
	
	//~ TCanvas *c1= new TCanvas ("c1","c1",800,800);
	//~ TH1F *h1= new TH1F("h1","h1",1000,0,600);
	//~ TFile *myfile= new TFile("datatree.root","READ");
	//~ TTree *treeData= (TTree* myfile)->Get("treeDatal");
	//~ treeData->Branch("evt_n", &evt_n, "evt_n/i");
	//~ treeData->Branch("peak_integral",&peak_integral,"peak_integral/D");
	//~ treeData->Branch("peak_maximum",&peak_maximum,"peak_maximum/D");
	//~ treeData->Branch("time_maximum",&time_maximum,"time_maximum/D");
	//~ treeData->Branch("baseline_mean",&baseline_mean,"baseline_mean/D");
	//~ treeData->Branch("baseline_rms",&baseline_rms,"baseline_rms/D");
	//~ unsigned int entries= treeData->GetEntries();
	//~ for (int i=0; i<entries; i++)
	//~ {
		//~ treeData->GetEntry(i);
		//~ h1->Fill(peak_integral);
	//~ }
	//~ h1->Draw();
	//~ TFile * file = TFile::Open("Run00700/datatree.root");
	//~ TTree *treeData= (TTree* file)->Get("treeData");
	//~ treeData->Draw("peak_integral >> h1(400, -20, 100)");
	
	
	
	TFile *file = new TFile("Run00700/datatree.root","READ","",5);
	if(!file->IsOpen() || file->IsZombie()) {cout << "Error: cannot open " << "datatree.root !\n";	exit(-1);}
	treeInput = (TTree*) file->Get("treeData");
	//~ treeInput->SetBranchAddress("evt_n",&evt_n);
	treeInput->SetBranchAddress("evt_n", &evt_n);
	treeInput->SetBranchAddress("peak_integral",&peak_integral);
	treeInput->SetBranchAddress("peak_maximum",&peak_maximum);
	treeInput->SetBranchAddress("time_maximum",&time_maximum);
	treeInput->SetBranchAddress("baseline_mean",&baseline_mean);
	treeInput->SetBranchAddress("baseline_rms",&baseline_rms);
	
	treeData->Draw("peak_integral >> h1(400, -40, 100)");
}