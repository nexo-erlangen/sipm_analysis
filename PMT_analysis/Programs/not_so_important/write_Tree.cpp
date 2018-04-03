void write_Tree(){
	stringstream TERMINAL;
	TERMINAL << "ls -d Run0* > list_runs.txt" << endl;
	system(TERMINAL.str().c_str());
	TERMINAL.str("");
	
	
	ifstream runfile;
	runfile.open("list_runs.txt");
/*	
	ofstream intfile;
	intfile.open( "list_runs_dark_root.txt" );
	intfile << "# " << "run" << "\t" << "integral pmt" << endl;
*/
	
	string folder = "";
	
	while( !runfile.eof() ) {
		runfile >> folder;
		//~ runfile >> folder;
		//~ runfile >> folder;
		//~ runfile >> folder;
		//~ runfile >> folder;
		//~ runfile >> folder;
		//~ runfile >> folder;
		//~ runfile >> folder;
		//~ runfile >> folder;
		//~ runfile >> folder;
		//~ runfile >> folder;
	
		ifstream brightfile;
		brightfile.open( ( folder +"/Auswertung/integral_15-35ns/peak_integral.dat" ).c_str() ) ;
		
		if ( brightfile.is_open() ) {
			cout << "File is open" << endl;
		}
		else{
			cout << "Failed to open file" << endl;
		}
		
		Float_t integral = 0.0, mean_y, delta_y;
		Int_t nlines = 0;
		
		TFile *f = new TFile( ( folder +"/Auswertung/integral_15-35ns/muell.root" ).c_str() , "RECREATE" );
		TTree *tree = new TTree( "histogramm", "Histogrammdaten" );
		tree->Branch( "integral", &integral );
		//~ TCanvas *c1 = new TCanvas();
		//~ TH1F *h1 = new TH1F( "h1", "y distribution", 1000, -200000, 500000 );
		
		while ( 1 ){
			//~ cout << "test" << endl;
			brightfile >> integral;
			//~ cout << integral << endl;
			tree->Fill();
			if ( !brightfile.good() ){
				break;
			}
			//~ if ( nlines < 5 ) printf ( "x = %8f, \t y = %8f\n", x , y );
			//~ h1->Fill(y);
			//~ cout << integral << endl;
			nlines++;
		}
		//~ printf("found %d points\n", nlines);
		cout << "Zeilen:\t" << nlines << endl;
		brightfile.close();
		tree->Write();
		f->Close();
		//~ c1->SetLogy();
		//~ h1->Draw();
		//~ mean_y = h1->GetMean();
		//~ delta_y = h1->GetStdDev();
		//~ std::cout << folder << ":\t" << mean_y << "\t" << delta_y << std::endl; // is  3907.41 wert passt noch nicht ganz, ansch mal spaltendreher in anticorr plots, aufpassen!!!
		//~ infy tfile << folder << "\t" << mean_y << endl;
		//~ delete tree;
	}
	
//	intfile.close();
	
	runfile.close();
}