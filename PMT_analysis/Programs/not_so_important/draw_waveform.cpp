void draw_waveform(){

	TFile *file1 = new TFile("/home/vault/capm/sn0527/PMT_LED_Temp/Run00702/input.root","READ");
	TTree *treeInput2;
	unsigned int pulses_length;
	
	double *time = new double [1002];
	double *voltage = new double [1002];
	unsigned int evt_n;
	treeInput2 = (TTree*) file1->Get("treeInput");
	treeInput2->SetBranchAddress("time",time);
	treeInput2->SetBranchAddress("voltage",voltage);
	treeInput2->SetBranchAddress("evt_n",&evt_n);
	treeInput2->SetBranchAddress("length",&pulses_length);

	TFile *file2 = new TFile("/home/vault/capm/sn0527/PMT_LED_Temp/Run00702/datatree.root","READ");
	//~ c1 = ROOT.TCanvas("c1","",800,800)
	//~ TCanvas *c1 = new TCanvas("c1","",800,800);
	TH1F* hist = new TH1F("hist", "hist", 1000,0,2);
	//TTree *Tree = file->Get("treeData");
	TTree *treeInput;
	double baseline_rms;
	treeInput = (TTree*) file2->Get("treeData");
	treeInput->SetBranchAddress("baseline_rms",&baseline_rms);
	
	unsigned int n = treeInput->GetEntries();
	
	// fill baseline_rms into histogram
	TH1F *hist1 = new TH1F ("hist", "hist", 400, 0, 1 );
	for( int i = 0; i < n; i++){
		treeInput->GetEntry(i);
		hist1->Fill(baseline_rms);
	}
	
	int binmax = hist1->GetMaximumBin();
	double maximum = hist1->GetXaxis()->GetBinCenter(binmax);
	cout << "maximum:\t" << maximum << endl;
	cout << "Binmax:\t:" << binmax <<endl;	
	
	// fit gaus to baseline:
	double left_end = maximum[j] - 0.02;
	double right_end = maximum[j] + 0.015;
	TCanvas *c2 = new TCanvas("c2","",800,800);
	TF1 *g1 = new TF1("g1", "gaus", left_end, right_end);
	hist1->Fit(g1, "R") ;
	c2->SetLogy(1);
	//~ hist1->Draw();
	//~ g1->Draw("same");
	
	double amplitude = g1->GetParameter(0);
	double mean = g1->GetParameter(1);
	double sigma = g1->GetParameter(2);
	
	//~ cout << "mean: " << mean << endl;
	//~ cout << "sigma: " << sigma << endl;
	//~ cout << "mean + sigma: " << mean + sigma << endl;	
	//~ cout << "left end: " << mean + sigma - 0.05*sigma << endl;	
	//~ cout << "right end: " << mean + sigma + 0.05*sigma << endl;	
	
	//~ cout << "mean - sigma: " << mean - sigma << endl;	
	//~ cout << "left end: " << mean - sigma - 0.05*sigma << endl;	
	//~ cout << "right end: " << mean - sigma + 0.05*sigma << endl;
	
	// Draw waveform
	for( int i =126; i < n; i++){
		treeInput->GetEntry(i);
		//~ double baseline =  baseline_rms;
		cout << "i: " << i << endl;
		if ( baseline_rms > ( mean - 1.0 * sigma - 0.05*sigma ) && baseline_rms < ( mean - 1.0 * sigma + 0.05*sigma ) ){
		//~ if ( i == 8 ){
			treeInput2->GetEntry(i);
			TCanvas *c3 = new TCanvas("c3","",800,800);
			TGraph *graph1 = new TGraph( pulses_length, time, voltage );
			//~ graph1->GetYaxis()->SetLimits(-1.,1.);
			graph1->Draw();
			cout << "test" << endl;
			break;
		}
	}

}