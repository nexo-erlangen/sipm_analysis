void get_baseline_cuts(){
	TFile *file1 = TFile::Open("datatree.root");
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
	
	
	string runnumber[15] = {"Run00695", "Run00696", "Run00697", "Run00698", "Run00699", "Run00700", "Run00701", "Run00702", "Run00703", "Run00704", "Run00705", "Run00706", "Run00707", "Run00708", "Run00709"};
	double baseline_cut[15] = {0};
	int number = 0;
	double cut[4] = {0};
	for ( int j = 0; j < 15; j++ ){
		//~ string runname = runnumber[i];
		//~ double end = number + 900000;
		TH1F *hist1 = new TH1F ("hist", "hist", 400, 0, 1 );
		for( int i = number; i < ( 900000 + number ); i++){
		//~ for( int i = 5400000; i < 6300000; i++){
			treeInput->GetEntry(i);
			hist1->Fill(baseline_rms);
			//~ cout << baseline_rms<<endl;
		}
		number += 900000;
		cout << number << endl;
		int binmax = hist1->GetMaximumBin();
		baseline_cut[j] = hist1->GetXaxis()->GetBinCenter(binmax);
		//~ cout << baseline_cut[j] << endl;
		//~ double x = hist1->GetXaxis()->GetBinCenter(binmax);
		//~ cout << x << endl;
		//~ cout << binmax <<endl;
		//~ hist1->Draw();
/*		double x = baseline_cut[i];
		TFile *file2 = new TFile( ( runname+"/fit_tree.root" ).c_str(), "recreate");
		TTree *tree_integral = treeInput->CopyTree(TString::Format("Entry$ > %i - 900000 && Entry$ < %i && baseline_rms < %f", number, number, x)); 
		tree_integral->Write();
		delete tree_integral;
		file2->Close();
		delete file2;                 // not necessary here since it has to be done again, is only cut determination
*/		
		delete hist1;
	}
	
	double sum = 0;
	for ( int i = 0; i < 4; i++){
		sum += baseline_cut[i];
	}	
	cut[0] = sum/4.0;
	
	sum = 0;
	for ( int i = 4; i < 8; i++){
		sum += baseline_cut[i];
	}	
	cut[1] = sum/4.0;	

	sum = 0;
	for ( int i = 8; i < 11; i++){
		sum += baseline_cut[i];
	}	
	cut[2] = sum/3.0;	

	sum = 0;
	for ( int i = 11; i < 15; i++){
		sum += baseline_cut[i];
	}	
	cut[3] = sum/4.0;	
	
	cout << "final cuts:\t" << cut[0] << endl;
	cout << "\t\t" << cut[1] << endl;
	cout << "\t\t" << cut[2] << endl;
	cout << "\t\t" << cut[3] << endl;
}
