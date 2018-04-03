void get_baseline_cuts_gaus(){
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
	double maximum[15] = {0};
	double baseline_cut_gaus[15][4] = {0};
	int number = 0;
	double cut[4][4] = {0};
	// Open a new ROOT file for the histograms
	TFile *filehistos = TFile::Open("myhistos.root", "RECREATE");
	for ( int j = 0; j < 3; j++ ){
		//~ int j = 1;
		number += 900000;
		//~ number += 900000;
		//~ number += 900000;
		//~ number += 900000;
		
		TH1F *hist1 = new TH1F ("hist", "hist", 400, 0, 1 );
		for( int i = number; i < ( 900000 + number ); i++){
			treeInput->GetEntry(i);
			hist1->Fill(baseline_rms);
		}
		hist1->Write();
		
		number += 900000;
		cout << "number:\t" << number << endl;
		int binmax = hist1->GetMaximumBin();
		maximum[j] = hist1->GetXaxis()->GetBinCenter(binmax);
		cout << "maximum:\t" << maximum[j] << endl;
		cout << "Binmax:\t:" << binmax <<endl;	
		
		// fit gaus to baseline:
		double left_end = maximum[j] - 0.02;
		double right_end = maximum[j] + 0.015;
		TF1 *g1 = new TF1("g1", "gaus", left_end, right_end);
		hist1->Fit(g1, "R") ;
		
		hist1->Draw();
		g1->Draw("same");
		
		double amplitude = g1->GetParameter(0);
		double mean = g1->GetParameter(1);
		double sigma = g1->GetParameter(2);
		
		cout << "amplitude:\t" << amplitude << endl;
		cout << "mean:\t" << mean << endl;
		cout << "sigma:\t" << sigma << endl;
		
		
		for ( int n = 0; n < 4; n++){
		//~ int n = 1;
			baseline_cut_gaus[j][n-1] = mean + n * sigma;
			cout << "baseline_cut_gaus:\t" << baseline_cut_gaus[j][n-1] << endl;
		}
		delete g1;
		delete hist1;
	}
	treeInput->Scan("baseline_rms", "", "", 100, 270000);
	filehistos->Close();
	//~ delete file1;
	/*
	// calculate cuts:
	cout << "----------------------------------------final cuts----------------------------------------" << endl;
	for ( int j = 0; j < 4; j++ ){
	
		double sum = 0;
		for ( int i = 0; i < 4; i++){
			sum += baseline_cut_gaus[i][j];
		}	
		cut[0][j] = sum/4.0;
		//~ cout << "cut:\t" << cut[0][0] << endl;
		
		sum = 0;
		for ( int i = 4; i < 8; i++){
			sum += baseline_cut_gaus[i][j];
		}
		cut[1][j] = sum/4.0;	

		sum = 0;
		for ( int i = 8; i < 11; i++){
			sum += baseline_cut_gaus[i][j];
		}	
		cut[2][j] = sum/3.0;	

		sum = 0;
		for ( int i = 11; i < 15; i++){
			sum += baseline_cut_gaus[i][j];
		}	
		cut[3][j] = sum/4.0;	
		
		cout << "Cut Nr.:\t" << j << endl;
		cout << "cut[0][j]:\t" << cut[0][j] << endl;
		cout << "cut[1][j]:\t" << cut[1][j] << endl;
		cout << "cut[2][j]:\t" << cut[2][j] << endl;
		cout << "cut[3][j]:\t" << cut[3][j] << endl;
		
	}
	
	double cuts_all[15][4] = { 0 };
	for ( int j = 0; j < 4; j++ ){
		for ( int i = 0; i < 4; i++ ){
			cuts_all[i][j] = cut[0][j];
		}
		for ( int i = 4; i < 8; i++){
			cuts_all[i][j] = cut[1][j];
		}
		for ( int i = 8; i < 11; i++){
			cuts_all[i][j] = cut[2][j];

		}	
		for ( int i = 11; i < 15; i++){
			cuts_all[i][j] = cut[3][j];
		}
	}
	
	// apply cuts:
	for ( int j = 0; j < 15; j++ ){
		string path = "/home/vault/capm/sn0527/PMT_LED_Temp/";
		string runname = runnumber[j];
		
		TFile *file2 = TFile::Open(( path + runname+"/datatree.root" ).c_str());
		TTree *treeInput = (TTree*) file2->Get("treeData");
		
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
		
		
		for ( int n = 0; n < 4; n++ ){
			double x = cuts_all[j][n];
			cout << "cut:\t" << x << endl;
			string m = to_string(n);
			TFile *file3 = new TFile( ( path + runname+"/fit_tree_" + m + "_sigma.root" ).c_str(), "recreate");
			TTree *tree_fit = treeInput->CopyTree(TString::Format("baseline_rms < %f", x)); 
			tree_fit->Write("", TObject::kOverwrite); // sonst werden wieder 2 Trees erzeugt
			delete tree_fit;
			file3->Close();
			delete file3;   
		}
		delete treeInput;
		file2->Close();
		delete file2;
	}*/
}
