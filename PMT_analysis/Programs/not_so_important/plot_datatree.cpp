void plot_datatree(){
	TFile *file = new TFile("/home/vault/capm/sn0527/PMT_LED_Temp/Run00700/datatree.root","READ","",5);
	TH1F* hist = new TH1F("hist", "hist", 500000,-20,70);
	//TTree *Tree = file->Get("treeData");
	TTree *treeInput;
	double peak_integral;
	treeInput = (TTree*) file->Get("treeData");
	treeInput->SetBranchAddress("peak_integral",&peak_integral);
	int entries = treeInput->GetEntries();
	double y = 0;
	for ( int i=0; i < entries; i++ ) {
		//y = entry->peak_integral;
		y = treeInput->GetEvent(i);
		hist->Fill(y);
	}
	hist->Rebin(4);
	hist->SetMinimum(1);
	hist->Draw();
}

/*
#Get histogram from ROOT file
c1 = ROOT.TCanvas("c1","",2500,2200)
file = ROOT.TFile("%s/%s.root" % (path, basename), "read") # hier passiert ansch der Fehler
#print "test"
Tree = file.Get("treeData")

print Tree.GetEntries()


hist = ROOT.TH1F("hist", "hist", 1000,-20,70)
#hist.Sumw2(1)

for entry in Tree:
	y = entry.peak_integral
	#print y 
	hist.Fill(y)
#print "test"#, ( "%s.root" % file )
hist.Rebin(4)
hist.SetMinimum(1)
hist.Draw();

*/
