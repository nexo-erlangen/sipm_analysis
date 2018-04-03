void apply_cuts(){
	string runnumber[15] = {"Run00695", "Run00696", "Run00697", "Run00698", "Run00699", "Run00700", "Run00701", "Run00702", "Run00703", "Run00704", "Run00705", "Run00706", "Run00707", "Run00708", "Run00709"};
	double baseline_cut[15] = {0.196875, 0.196875, 0.196875, 0.196875, 0.241875, 0.241875, 0.241875, 0.241875, 0.23625, 0.23625, 0.23625, 0.20625, 0.20625, 0.20625, 0.20625};
	for ( int j = 0; j < 15; j++ ){
		string path = "/home/vault/capm/sn0527/PMT_LED_Temp/";
		string runname = runnumber[j];
		
		TFile *file1 = TFile::Open(( path + runname+"/datatree.root" ).c_str());
		TTree *treeInput = (TTree*) file1->Get("treeData");
		
		unsigned int evt_n;
		double peak_integral;
		double peak_maximum;
		double time_maximum;
		double baseline_mean;
		double baseline_rms;
		
		treeInput->SetBranchAddress("evt_n", &evt_n);
		treeInput->SetBranchStatus("evt_n",0);
		treeInput->SetBranchAddress("peak_integral",&peak_integral);
		treeInput->SetBranchStatus("peak_integral",1);
		treeInput->SetBranchAddress("peak_maximum",&peak_maximum);
		treeInput->SetBranchStatus("peak_maximum",0);
		treeInput->SetBranchAddress("time_maximum",&time_maximum);
		treeInput->SetBranchStatus("time_maximum",0);
		treeInput->SetBranchAddress("baseline_mean",&baseline_mean);
		treeInput->SetBranchStatus("baseline_mean",0);
		treeInput->SetBranchAddress("baseline_rms",&baseline_rms);
		treeInput->SetBranchStatus("baseline_rms",1);
		
		double x = baseline_cut[j];
		TFile *file2 = new TFile( ( path + runname+"/fit_tree.root" ).c_str(), "recreate");
		TTree *tree_fit = treeInput->CopyTree(TString::Format("baseline_rms < %f", x)); 
		tree_fit->Write("", TObject::kOverwrite); // sonst werden wieder 2 Trees erzeugt
		delete tree_fit;
		file2->Close();
		delete file2;   
		
		delete treeInput;
		file1->Close();
		delete file1;
	}	
}