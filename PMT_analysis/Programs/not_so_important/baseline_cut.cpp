void baseline_cut(){
	
	string path = "/home/vault/capm/sn0527/PMT_LED_Temp/";
	string runnumber[15] = {"Run00695", "Run00696", "Run00697", "Run00698", "Run00699", "Run00700", "Run00701", "Run00702", "Run00703", "Run00704", "Run00705", "Run00706", "Run00707", "Run00708", "Run00709"};
	double maximum[15] = {0};
	double baseline_cut_gaus[15] = {0};
	double cut[4] = {0};

	//~ for ( int j = 0; j < 15; j++ ){
		// reading of tree for fitting
		int j = 13;
		TFile *file1 = TFile::Open( ( path + runnumber[j] + "/datatree_900000.root" ).c_str() );
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

		// fill baseline_rms into histogram
		TCanvas *c1 = new TCanvas( "c1", "c1", 800, 800 );
		TH1F *hist1 = new TH1F ("hist", "hist", 400, 0, 1 );
		for( int i = 0; i < 900000 ; i++){
			treeInput->GetEntry(i);
			hist1->Fill(baseline_rms);
		}

		int binmax = hist1->GetMaximumBin();
		maximum[j] = hist1->GetXaxis()->GetBinCenter(binmax);
		cout << "maximum:\t" << maximum[j] << endl;
		cout << "Binmax:\t:" << binmax <<endl;	
		
		// fit gaus to baseline:
		double left_end = maximum[j] - 0.02;
		double right_end = maximum[j] + 0.015;
		
		//~ double left_end = maximum[j] - 0.035;
		//~ double right_end = maximum[j] + 0.025;
		TF1 *g1 = new TF1("g1", "gaus", left_end, right_end);
		hist1->Fit(g1, "R") ;
		
		c1->SetLogy(0);
		hist1->Draw();
		g1->Draw("same");
		
		double amplitude = g1->GetParameter(0);
		double mean = g1->GetParameter(1);
		double sigma = g1->GetParameter(2);
		
		cout << "amplitude:\t" << amplitude << endl;
		cout << "mean:\t" << mean << endl;
		cout << "sigma:\t" << sigma << endl;
		
		baseline_cut_gaus[j] = mean + sigma;
		cout << "baseline_cut_gaus:\t" << baseline_cut_gaus[j] << endl;
	
		//~ delete g1;
		//~ delete hist1;
		//~ delete treeInput;
		//~ file1->Close();
		//~ delete file1;
		// delete treeInput; // segmentation violation wenn an dieser Stelle einkommentiert :(
	//~ }

	// calculate cut:
	cout << "\n\n";
	cout << "----------------------------------------final cuts----------------------------------------" << endl;
	cout << "\n\n";
		
	// loops over different runs
	double sum = 0;
	for ( int i = 0; i < 4; i++){
		sum += baseline_cut_gaus[i];
	}	
	cut[0] = sum/4.0;
	//~ cout << "cut:\t" << cut[0][0] << endl;
	
	sum = 0;
	for ( int i = 4; i < 8; i++){
		sum += baseline_cut_gaus[i];
	}
	cut[1] = sum/4.0;	

	sum = 0;
	for ( int i = 8; i < 11; i++){
		sum += baseline_cut_gaus[i];
	}	
	cut[2] = sum/3.0;	

	sum = 0;
	for ( int i = 11; i < 15; i++){
		sum += baseline_cut_gaus[i];
	}	
	cut[3] = sum/4.0;	
	
	cout << "Cut Nr.:\t" << j << endl;
	cout << "cut[0]:\t" << cut[0] << endl;
	cout << "cut[1]:\t" << cut[1] << endl;
	cout << "cut[2]:\t" << cut[2] << endl;
	cout << "cut[3]:\t" << cut[3] << endl;
		
	
	double cuts_all[15] = { 0 };
	for ( int i = 0; i < 4; i++ ){
		cuts_all[i] = cut[0];
	}
	for ( int i = 4; i < 8; i++){
		cuts_all[i] = cut[1];
	}
	for ( int i = 8; i < 11; i++){
		cuts_all[i] = cut[2];

	}	
	for ( int i = 11; i < 15; i++){
		cuts_all[i] = cut[3];
	}
	
	//~ // apply cuts:
	//~ cout << "\n\n";
	//~ cout << "----------------------------------------writing trees with cuts----------------------------------------" << endl;
	//~ cout << "\n\n";
	//~ for ( int j = 0; j < 15; j++ ){
		//~ // string path = "/home/vault/capm/sn0527/PMT_LED_Temp/";
		//~ string runname = runnumber[j];
		
		//~ TFile *file2 = TFile::Open(( path + runname+"/datatree.root" ).c_str());
		//~ TTree *treeInput2 = (TTree*) file2->Get("treeData");
		
		//~ unsigned int evt_n;
		//~ double peak_integral;
		//~ double peak_maximum;
		//~ double time_maximum;
		//~ double baseline_mean;
		//~ double baseline_rms;
		
		//~ treeInput2->SetBranchAddress("evt_n", &evt_n);
		//~ treeInput2->SetBranchStatus("evt_n",0);
		//~ treeInput2->SetBranchAddress("peak_integral",&peak_integral);
		//~ treeInput2->SetBranchStatus("peak_integral",1);
		//~ treeInput2->SetBranchAddress("peak_maximum",&peak_maximum);
		//~ treeInput2->SetBranchStatus("peak_maximum",0);
		//~ treeInput2->SetBranchAddress("time_maximum",&time_maximum);
		//~ treeInput2->SetBranchStatus("time_maximum",0);
		//~ treeInput2->SetBranchAddress("baseline_mean",&baseline_mean);
		//~ treeInput2->SetBranchStatus("baseline_mean",0);
		//~ treeInput2->SetBranchAddress("baseline_rms",&baseline_rms);
		//~ treeInput2->SetBranchStatus("baseline_rms",1);
		
		
		//~ double x = cuts_all[j];
		//~ cout << "run: " << runname << "   " << "cut: " << x << endl;
		//~ string m = to_string(1);
		//~ TFile *file3 = new TFile( ( path + runname+"/fit_tree_cut_" + m + "_sigma.root" ).c_str(), "recreate");
		//~ TTree *tree_fit = treeInput2->CopyTree(TString::Format("baseline_rms < %f", x)); 
		//~ tree_fit->Write("", TObject::kOverwrite); // sonst werden wieder 2 Trees erzeugt
		//~ delete tree_fit;
		//~ file3->Close();
		//~ delete file3;   
			
		//~ delete treeInput2;
		//~ file2->Close();
		//~ delete file2;
	//~ }
}
