/*
 * FOLDER.cpp
 *
 *  Created on: Apr 15, 2016
 *      Author: sn0515
 */

#include "FOLDER.h"

FOLDER::FOLDER() {}

FOLDER::~FOLDER() {
	delete treeInput;
	delete x;
	delete voltage;
//	delete treeAnalysis;
//	delete treeFit;
}

unsigned int FOLDER::counterWaves=0;

FOLDER::FOLDER(MAININPUT* input, const string folder) {
	FOLDER::counterWaves=0;
	WAVE::counter=0;
	info=input;
	filepath=folder;
	cout << endl << "=========Folder====================================================" << endl;
	cout << filepath << endl;
	setStructure();
}

void FOLDER::setStructure() const {
	stringstream TERMINAL;
	TERMINAL.str("");
	TERMINAL << "mkdir -m 770 -p " << filepath << "Auswertung/" << endl;
	TERMINAL << "mkdir -m 770 -p " << filepath << "Auswertung/temp" << endl;
	TERMINAL << "mkdir -m 770 -p " << filepath << "Auswertung/frames" << endl;
	if(info->bool_root || info->bool_fit || info->bool_awesome) {
		TERMINAL << "mkdir -m 770 -p " << filepath << "Auswertung/frames/1pe" << endl;
		TERMINAL << "mkdir -m 770 -p " << filepath << "Auswertung/frames/2pe" << endl;
		TERMINAL << "mkdir -m 770 -p " << filepath << "Auswertung/frames/3pe" << endl;
		TERMINAL << "mkdir -m 770 -p " << filepath << "Auswertung/frames/4pe" << endl;
		TERMINAL << "mkdir -m 770 -p " << filepath << "Auswertung/frames/5pe" << endl;
	}
	system(TERMINAL.str().c_str());
	TERMINAL.str("");
}

void FOLDER::readConditions() {
//	size_t pos1 = filepath.find("C_");
//	size_t pos2 = filepath.find("V_");
//	size_t pos3 = filepath.find("mm_");
//	if(pos1!=std::string::npos)	{
//		size_t pos11 =  filepath.rfind("_",pos1);
//		stringstream(filepath.substr(pos11+1,pos1-pos11-1)) >> data.temp;
//	}
//	if(pos2!=string::npos) {
//		size_t pos22 =  filepath.rfind("_",pos2);
//		stringstream(filepath.substr(pos22+1,pos2-pos22-1)) >> data.volt;
//	}
//	if(pos3!=string::npos) {
//		size_t pos33 =  filepath.rfind("_",pos3);
//		stringstream(filepath.substr(pos33+1,pos3-pos33-1)) >> data.dist;
//	}

//	TFile out((filepath+"conditions.root").c_str(),"RECREATE");
//	if(!out.IsOpen()) {
//		cout << "Error: cannot open " << filepath << "conditions.root !" << endl;
//	}
//	tree0 = new TTree("tree0","tree0");
//	tree0->Branch("run",&data.run,"run/i");
//	tree0->Branch("sipmBias",&data.sipmBias,"sipmBias/D");
//	tree0->Branch("pmtBias",&data.pmtBias,"pmtBias/D");
//	tree0->Branch("dist",&data.dist,"dist/D");
//	tree0->Branch("temp",&data.temp,"press/D");
////	tree0->Branch("sign",&sign,"sign/I");
//	//reading conditions ...
//	tree0->Fill();
//	tree0->Write();
//	out.Close();
}

void FOLDER::readFolder() {
	stringstream TERMINAL;
	ifstream INPUT;
	TERMINAL.str("");
	TERMINAL << "(cd " << filepath << "; find *.txt *.dat *.trc 2> /dev/null > Auswertung/temp/input.log)" << endl;
	system(TERMINAL.str().c_str());
	TERMINAL.str("");

	INPUT.open((filepath+"Auswertung/temp/input.log").c_str());
	while(!INPUT.eof())	{
		string a;
		INPUT >> a;
		files.push_back(a);
	};
	files.pop_back();
	INPUT.close();
	cout << "Number of Files:\t\t" << files.size() << endl;
}

void FOLDER::readFiles() {
	double start2=getTime(0,0);
	cout << "Reading Data:";

	TFile *outFile = new TFile((filepath+"input.root").c_str(),"RECREATE","",5);
	if(!outFile->IsOpen() || outFile->IsZombie()) {cout << "Error: cannot open " << filepath << "input.root !" << endl;}
	TChain *chain = new TChain("treeInput");
	vector<string> files_temp;

	for(unsigned int m=0; m<files.size(); m++) {//loop over all files with data
		CFILE* file= new CFILE(this,m);
		file->readFile();
		chain->Add((file->getRootFile()).c_str());
		files_temp.push_back(file->getRootFile());
		delete file;
		loadbar(m, files.size(), 100, 100, start2);
	}
	cout << endl;
	chain->Merge((filepath+"input.root").c_str());

	stringstream TERMINAL;
	TERMINAL.str("");
	for(unsigned int m=0; m<files_temp.size(); m++) {//delete all temporary root files
		TERMINAL << "(rm " << files_temp.at(m) << ")" << endl;
	}
	system(TERMINAL.str().c_str());
	TERMINAL.str("");

	delete chain;
	outFile->Close();
	delete outFile;
}

void FOLDER::readRootInput() {
	double start2=getTime(0,0);
	cout << "Opening Data:\r";

	TFile *file = new TFile((filepath+"input_old.root").c_str(),"READ","",5);
	cout << "path: "<<  filepath+"input_old.root" << endl;
	if(!file->IsOpen() || file->IsZombie()) {cout << "Error: cannot open " << filepath << "input_old.root !\n";	exit(-1);}
	treeInput = (TTree*) file->Get("treeInput");
	treeInput->SetBranchAddress("evt_n",&evt_n);
	treeInput->SetBranchAddress("length",&length);
	treeInput->SetBranchAddress("trigger",&trigger);
	treeInput->SetBranchAddress("time",x);
	treeInput->SetBranchAddress("voltage",voltage);
//	treeInput->Scan("*");
	FOLDER::counterWaves = (Int_t)treeInput->GetEntries();

//	tree0 = (TTree*) file->Get("tree0");
//	tree0->SetBranchAddress("run",&data.run);
//	tree0->SetBranchAddress("sipmBias",&data.sipmBias);
//	tree0->SetBranchAddress("pmtBias",&data.pmtBias);
//	tree0->SetBranchAddress("dist",&data.dist);
//	tree0->SetBranchAddress("temp",&data.temp);

	loadbar(0, 1, 1, 100, start2);
}

void FOLDER::analysisWaves() {
	double start3=getTime(0,0);
	cout << "\nAnalysing Data:\r";

//	TFile *out = new TFile((filepath+"analysisWave.root").c_str(),"RECREATE");
//	if(!out->IsOpen() || out->IsZombie()) {
//		cout << "Error: cannot open " << filepath << "analysisWave.root !" << endl;
//	}
//	treeAnalysis = new TTree("treeAnalysis","treeAnalysis");
//	treeAnalysis->Branch("evt_n",&evt_n,"evt_n/i");
//	treeAnalysis->Branch("trigger_time",&trigger_time,"trigger_time/D");
//	treeAnalysis->Branch("peaks_n",&peaks_n,"peaks_n/i");
//	treeAnalysis->Branch("peak_time",&peak_time,"peak_time/D");
//	treeAnalysis->Branch("peak_value",&peak_value,"peak_value/D");
//	treeAnalysis->Branch("base",&base,"base/D");
//	treeAnalysis->Branch("fluct",&fluct,"fluct/D");
	TFile *out2 = new TFile((filepath+"analysisFit.root").c_str(),"RECREATE");
	if(!out2->IsOpen() || out2->IsZombie()) {cout << "Error: cannot open " << filepath << "analysisFit.root !" << endl; exit(-1);}

	treeFit = new TTree("treeFit","treeFit");
	treeFit->Branch("evt_n",&evt_n,"evt_n/i");
	treeFit->Branch("peaks_n",&peaks_n,"peaks_n/i");
	treeFit->Branch("peak_i",&peak_i,"peak_i/i");
//	treeFit->Branch("TriggerTime",&dPulse.TriggerTime,"TriggerTime/D");
	treeFit->Branch("TriggerTime",&trigger,"TriggerTime/D");
	treeFit->Branch("Baseline",&dPulse.Baseline,"Baseline/D");
	treeFit->Branch("Amp",&dPulse.Amp,"Amp/D");
	treeFit->Branch("Time",&dPulse.Time,"Time/D");
	treeFit->Branch("FitLowLimit",&dPulse.FitLowLimit,"FitLowLimit/D");
	treeFit->Branch("FitHighLimit",&dPulse.FitHighLimit,"FitHighLimit/D");
	treeFit->Branch("FitBaseline",&dPulse.FitBaseline,"FitBaseline/D");
	treeFit->Branch("FitTime",&dPulse.FitTime,"FitTime/D");
	treeFit->Branch("FitAmp",&dPulse.FitAmp,"FitAmp/D");
	treeFit->Branch("FitRiseTime",&dPulse.FitRiseTime,"FitRiseTime/D");
	treeFit->Branch("FitFallTime",&dPulse.FitFallTime,"FitFallTime/D");
//	treeFit->Branch("FitTime2Frac",&dPulse.FitTime2Frac,"FitTime2Frac/D");
//	treeFit->Branch("FitFallTime2",&dPulse.FitFallTime2,"FitFallTime2/D");
	treeFit->Branch("FitChi2",&dPulse.FitChi2,"FitChi2/D");
	treeFit->Branch("FitNDF",&dPulse.FitNDF,"FitNDF/D");
	treeFit->Branch("FuncTime",&dPulse.FuncTime,"FuncTime/D");
	treeFit->Branch("FuncAmp",&dPulse.FuncAmp,"FuncAmp/D");

	unsigned int limit=FOLDER::counterWaves;
	if(info->bool_frames && info->frames<FOLDER::counterWaves) {
		limit=info->frames;
	}
	for(unsigned int i=0; i<limit; ++i) { //loop over all waves
		WAVE *wave = new WAVE(this, i);
		wave->analysis();
		if(WAVE::counter<=info->single) {
			wave->writeWave("","","");
		}
		if(info->bool_fit) {
			int nr_peaks=wave->getNrPeaks("fit");
			mPulse=wave->getPulses();
			for(int i=0; i<nr_peaks; i++) {
				evt_n=WAVE::counter;
				peak_i=i+1;
				peaks_n=nr_peaks;
				dPulse=mPulse[i];
				treeFit->Fill();
			}
		}
//		if(info->bool_find) {
//			int nr_peaks = wave->getNrPeaks("find") ;
//			for(int i=nr_peaks-1; i>=0; i--) {
//				evt_n=WAVE::counter;
//				int location = data.dummy_time.size()-1-i;
//				peaks_n=data.dummy_time.at(location);
//				peak_value=data.maximum_root.at(location);
//				peak_time=data.maximum_root_time.at(location);
//				trigger_time=data.trigger_time.at(location);
//				treeAnalysis->Fill();
//			}
//		}
		delete wave;
		loadbar(i, limit, 100, 100, start3);
	}
//	treeAnalysis->Scan("*");
//	treeAnalysis->Write();
//	delete treeAnalysis;
//	out->Close();
//	delete out;
//	treeFit->Scan("*");

	treeFit->Write();
	delete treeFit;
	out2->Close();
	delete out2;
}

void FOLDER::readRootFit() {
	double start2=getTime(0,0);
	TFile *file = TFile::Open((filepath+"analysisFit.root").c_str(),"READ");
	if(!file->IsOpen()) {
		cout << "Error: cannot open " << filepath << "analysisFit.root !\n";
		exit(-1);
	}
	cout << "Opening Analysis:\r";
	treeFit = (TTree*) file->Get("treeFit");
	treeFit->SetBranchAddress("evt_n",&evt_n);
	treeFit->SetBranchAddress("peaks_n",&peaks_n);
	treeFit->SetBranchAddress("peak_i",&peak_i);
	treeFit->SetBranchAddress("TriggerTime",&dPulse.TriggerTime);
	treeFit->SetBranchAddress("Baseline",&dPulse.Baseline);
	treeFit->SetBranchAddress("Amp",&dPulse.Amp);
	treeFit->SetBranchAddress("Time",&dPulse.Time);
	treeFit->SetBranchAddress("FitLowLimit",&dPulse.FitLowLimit);
	treeFit->SetBranchAddress("FitHighLimit",&dPulse.FitHighLimit);
	treeFit->SetBranchAddress("FitBaseline",&dPulse.FitBaseline);
	treeFit->SetBranchAddress("FitTime",&dPulse.FitTime);
	treeFit->SetBranchAddress("FitAmp",&dPulse.FitAmp);
	treeFit->SetBranchAddress("FitRiseTime",&dPulse.FitRiseTime);
	treeFit->SetBranchAddress("FitFallTime",&dPulse.FitFallTime);
//	treeFit->SetBranchAddress("FitTime2Frac",&dPulse.FitTime2Frac);
//	treeFit->SetBranchAddress("FitFallTime2",&dPulse.FitFallTime2);
	treeFit->SetBranchAddress("FitChi2",&dPulse.FitChi2);
	treeFit->SetBranchAddress("FitNDF",&dPulse.FitNDF);
	treeFit->SetBranchAddress("FuncTime",&dPulse.FuncTime);
	treeFit->SetBranchAddress("FuncAmp",&dPulse.FuncAmp);

	FOLDER::counterWaves = (Int_t)treeFit->GetEntries();
//	treeFit->Scan("*");
	loadbar(0, 1, 1, 100, start2);
}

void FOLDER::analysisNasty() {
	double start2=getTime(0,0);
	cout << "\nAnalyse Folder:\r";

//	writePulseData("all", -50.0, 9000.0, 1,2,2,100);
	if(info->bool_2dWave) {
		make2DHistogramWave();
	}
	if(info->bool_prompt || info->bool_crosstalk) {
		stringstream TERMINAL;
		TERMINAL.str("");
		TERMINAL << "mkdir -p " << filepath << "Auswertung/prompt" << endl;
		system(TERMINAL.str().c_str());
		TERMINAL.str("");

		string prompt_histo="prompt/prompt";
		const unsigned int counter = info->nr_peaks;
		Data max[counter+1];
		Data min[counter+1];
		vector<Data> HISTO;

		writePulseData(prompt_histo, -1000.0, 9000.0, 1,1,20,100);
		if(info->bool_prompt) {
			makingHistogram(filepath, prompt_histo, prompt_histo+"_histo",0.15,2);
//			readGnuHistogram(filepath,"Auswertung/"+prompt_histo+"_histo.dat", &HISTO);
//			findHistogramPeaks(filepath, prompt_histo+"_extremes", counter, &HISTO, max ,min);
//			writeGnuScriptSpec(filepath, prompt_histo, prompt_histo+"_histo", max, min, counter);
//			writeGnuScriptGain(filepath,prompt_histo+"_histo-param", prompt_histo+"-linear", counter);
		}
		if(info->bool_crosstalk) {
			calcCrosstalkLinear(filepath, prompt_histo, info->pe1, info->pe2);
			calcCrosstalkQuadratic(filepath, prompt_histo, info->gain);
		}
	}
	if(info->bool_dcr) {
		stringstream TERMINAL;
		TERMINAL.str("");
		TERMINAL << "mkdir -p " << filepath << "Auswertung/dcr" << endl;
		TERMINAL << "mkdir -p " << filepath << "Auswertung/rectime" << endl;
		system(TERMINAL.str().c_str());
		TERMINAL.str("");

		string output = "dcr/dcr";
		writeTimeDiff(output,"rectime/rec",-10, 20, info->gain, info->recTime, info->dcr);
		writeGnuScriptDCR(filepath, output, output+"_plot",info->dcrBegin) ;
	}
	if(info->bool_rectime) {
		stringstream TERMINAL;
		TERMINAL.str("");
		TERMINAL << "mkdir -p " << filepath << "Auswertung/rectime" << endl;
		system(TERMINAL.str().c_str());
		TERMINAL.str("");

		string output = "rectime/rec";
//		writePulseDataRec(output,info->gain,0.6);
		writeTimeDiff("dcr/dcr",output,-10, 20, info->gain, info->recTime, info->dcr);
		writeGnuScriptRecTime(filepath, output , output+"_plot", 0.3, 1);
		writeTimeDiff2("heightvstime",-10, 20, info->gain);
	}
	if(info->bool_2dHisto) {
		stringstream TERMINAL;
		TERMINAL.str("");
		TERMINAL << "mkdir -p " << filepath << "Auswertung/2dHisto" << endl;
		system(TERMINAL.str().c_str());
		TERMINAL.str("");

		string output="2dHisto/2dHisto";
		double xlow = -100.0;
		double xhigh = 900.0;
		double ylow = 0.0;
		double yhigh = 45.0;
		int xbins = 1000;
		int ybins = 450;
		writePulseData(output, xlow, xhigh, 1,20,20,100);
		make2DHistogram(output, output+"-matrix", xlow, xhigh, xbins, ylow, yhigh, ybins);
		writeGnuScript2DHisto(filepath, output+"-matrix", output+"-plot", xlow, xhigh, xbins, ylow, yhigh, ybins);
	}
	loadbar(0, 1, 100, 100, start2);
}

void FOLDER::writePulseData(string file, double bound_lower, double bound_upper, unsigned int nr_peak_lower=1, unsigned int nr_peak_upper=20, unsigned int nr_peaks=20, double chi2_max=10.0) {
	ofstream OUTPUT;
	OUTPUT.open((filepath+"Auswertung/"+file+".dat").c_str());
	int nr_entries=(Int_t)treeFit->GetEntries();
	for(int i=0; i<nr_entries; i++) {
		treeFit->GetEntry(i);
		if(peaks_n<=nr_peaks && dPulse.FitChi2/dPulse.FitNDF<chi2_max && dPulse.FuncAmp > 0.1) {
			if(dPulse.FuncTime>bound_lower && dPulse.FuncTime<bound_upper) {
				if(peak_i>=nr_peak_lower && peak_i<=nr_peak_upper) {
					OUTPUT << dPulse.FuncTime << "\t" <<  dPulse.FuncAmp << "\t" << dPulse.FitTime << "\t" <<  dPulse.FitAmp << "\t" << peak_i <<  endl;
				}
			}
		}
	}
	OUTPUT.close();
}

void FOLDER::writeTimeDiff(string file, string file2, double bound_lower, double bound_upper, double gain, double recTime, double dcr) {
	TCanvas *c2 = new TCanvas("c2","c2",600,400);
    c2->SetLogx();
    c2->SetLogy();
	const Int_t nbins = 50;
	Double_t xmin = 1e0;
	Double_t xmax = 1e10;
	Double_t logxmin = log10(xmin);
	Double_t logxmax = log10(xmax);
	Double_t binwidth = (logxmax-logxmin)/nbins;
	Double_t xbins[nbins+1];
	xbins[0] = xmin;
	for (Int_t i=1;i<=nbins;i++) {
	  xbins[i] = xmin + TMath::Power(10,logxmin+i*binwidth);
	}
	TH1D *h1 = new TH1D("h1","h1",nbins,xbins);
	ofstream OUTPUT1;
	OUTPUT1.open((filepath+"Auswertung/"+file2+".dat").c_str());
	int nr_entries=(Int_t)treeFit->GetEntries();
//	cout << "nr_entries: " << nr_entries << endl;
	for(int i=0; i<nr_entries-1; i++) {
		treeFit->GetEntry(i);
		if(peak_i==1) {
			if(dPulse.FuncAmp>(0.5*gain) && dPulse.FuncAmp<(1.5*gain) && dPulse.FuncTime>bound_lower && dPulse.FuncTime<bound_upper) {
				double time_temp_i = dPulse.FuncTime + dPulse.TriggerTime;
				treeFit->GetEntry(i+1);
				if(!(dPulse.TriggerTime==0.0 && peak_i==1)) {
					if(dPulse.FuncTime>-30 && dPulse.FuncAmp/gain > 0.32) {
						double time_temp_ii = dPulse.FuncTime + dPulse.TriggerTime;
						OUTPUT1 << (time_temp_ii-time_temp_i) << "\t" << (dPulse.FuncAmp/gain) << endl;
						h1->Fill(time_temp_ii-time_temp_i);
//						cout << "(time_temp_ii-time_temp_i): " << (time_temp_ii-time_temp_i) << endl;
//						cout << "amplitude: " << dPulse.FuncAmp/gain << endl;
					}
				}
			}
		}
	}
	OUTPUT1.close();

	h1->Sumw2();
	double scale = h1->GetEntries();
//	cout << "entries: " << h1->GetEntries() << endl;
//	cout << "bins: " << h1->GetNbinsX() << endl;
	for(int iBin=1; iBin<=h1->GetNbinsX(); iBin++){
		h1->SetBinContent(iBin,h1->GetBinContent(iBin)/h1->GetBinWidth(iBin)/scale);
		h1->SetBinError(iBin,h1->GetBinError(iBin)/h1->GetBinWidth(iBin)/scale);
	}
	ofstream OUTPUT;
	OUTPUT.open((filepath+"Auswertung/"+file+".dat").c_str());
	for(int i =0; i<nbins; i++) {
		OUTPUT << h1->GetBinCenter(i) << "\t"<< h1->GetBinContent(i) << "\t" << h1->GetBinError(i) << endl;
	}

	int gNAP=5;
	double * NAP = new double [gNAP];
	double * TauAP = new double [gNAP];
	double * recFrac = new double [gNAP];
	double * NAP_err = new double [gNAP];
	double * TauAP_err = new double [gNAP];
	double * recFrac_err = new double [gNAP];
	NAP[0]=0.5;TauAP[0]=100.0;recFrac[0]=0.95;
	NAP[1]=0.05;TauAP[1]=1000.0;recFrac[1]=1.0;
	NAP[2]=0.01;TauAP[2]=10000.0;recFrac[2]=1.0;
	NAP[3]=0.01;TauAP[3]=15000.0;recFrac[3]=1.0;
	NAP[4]=0.01;TauAP[4]=16000.0;recFrac[4]=1.0;
//	NAP[0]=0.5;TauAP[0]=100.0;recFrac[0]=1.0;
//	NAP[1]=0.3;TauAP[1]=150.0;recFrac[1]=1.0;
//	NAP[2]=0.2;TauAP[2]=300.0;recFrac[2]=1.0;
//	NAP[3]=0.1;TauAP[3]=1000.0;recFrac[3]=1.0;
//	NAP[4]=0.1;TauAP[4]=2000.0;recFrac[4]=1.0;
//	NAP[5]=0.1;TauAP[5]=3000.0;recFrac[5]=1.0;
//	NAP[6]=0.1;TauAP[6]=4000.0;recFrac[6]=1.0;
//	NAP[7]=0.05;TauAP[7]=10000.0;recFrac[7]=1.0;

	double fDeadTime = 30.0;
	double fReadoutGapStart = 9e3;
	double fReadoutGapEnd = 17e3;
	double fMinTimeWhenRec = 5.0;
	double recoveryTime=recTime;

	TF1* FDNAPWGapRec = new TF1("FDNAPWGapRec",funcDNAPWGapRec,2,1e10,7+3*gNAP);
	FDNAPWGapRec->FixParameter(0,dcr*1.0e-9);
	FDNAPWGapRec->FixParameter(1,gNAP);
	FDNAPWGapRec->FixParameter(2,fDeadTime);
	FDNAPWGapRec->FixParameter(3,fReadoutGapStart);
	FDNAPWGapRec->FixParameter(4,fReadoutGapEnd);
	FDNAPWGapRec->FixParameter(5,fMinTimeWhenRec);
	FDNAPWGapRec->FixParameter(6,recoveryTime);
	for(int iAP=0; iAP<gNAP; iAP++){
		FDNAPWGapRec->SetParameter(7+3*iAP,NAP[iAP]);
		FDNAPWGapRec->SetParLimits(7+3*iAP,0.0, 1.0);
		FDNAPWGapRec->SetParameter(8+3*iAP,TauAP[iAP]);
//		FDNAPWGapRec->FixParameter(8+3*iAP,TauAP[iAP]);
		FDNAPWGapRec->SetParameter(9+3*iAP,recFrac[iAP]);
		if(iAP>0) {FDNAPWGapRec->FixParameter(9+3*iAP,recFrac[iAP]);}
		else {FDNAPWGapRec->SetParLimits(9+3*iAP,0.0, 1.0);}
	}
	TVirtualFitter::SetMaxIterations(10000);
	h1->Fit("FDNAPWGapRec","RM");
	c2->SaveAs((filepath+"Auswertung/dcr/AP.pdf").c_str());
	for(int iAP=0; iAP<gNAP; iAP++){
		NAP[iAP]=FDNAPWGapRec->GetParameter(7+3*iAP);
		NAP_err[iAP]=FDNAPWGapRec->GetParError(7+3*iAP);
		TauAP[iAP]=FDNAPWGapRec->GetParameter(8+3*iAP);
		TauAP_err[iAP]=FDNAPWGapRec->GetParError(8+3*iAP);
		recFrac[iAP]=FDNAPWGapRec->GetParameter(9+3*iAP);
		recFrac_err[iAP]=FDNAPWGapRec->GetParError(9+3*iAP);
	}

	// calculate parameters from function
	double tAPtot=0.;
	double tAP1us=0.;
	for(int iAP=0; iAP<gNAP; iAP++){
		double tRecFrac=1.;
		if(iAP==0) tRecFrac = recFrac[iAP];
		double tLa = TauAP[iAP]*recoveryTime/(TauAP[iAP]+recoveryTime);
		double tAP = NAP[iAP]*(exp(-fDeadTime/TauAP[iAP])*(1-tRecFrac)+tRecFrac*(exp(-fMinTimeWhenRec/TauAP[iAP])-tLa/TauAP[iAP]*exp(-fMinTimeWhenRec/tLa)));
		cout << iAP
		 << " Time constant (mus) " << (TauAP[iAP])*1e-3
		 << " us | Rate " << tAP << endl;
		tAPtot+=tAP;
		tAP1us+=tAP;
		tAP1us-=(NAP[iAP]*(exp(-1000./TauAP[iAP])-tLa/TauAP[iAP]*exp(-1000./tLa)));
	}
	cout << "AP tot " << tAPtot << " | AP 1us " << tAP1us << endl;


	OUTPUT.close();
	OUTPUT.open((filepath+"Auswertung/dcr/AP_func_data.dat").c_str());
	for(int i = 1 ; i<=100*nbins; i++) {
		int x = xmin + TMath::Power(10,logxmin+i*binwidth/100.0);
		OUTPUT << x << "\t" << FDNAPWGapRec->Eval(x) << endl;
	}
	OUTPUT.close();
	OUTPUT.open((filepath+"Auswertung/dcr/AP_param.dat").c_str());
	OUTPUT << filepath << endl;
	for(int i = 0 ; i<gNAP; i++) {
		OUTPUT << i << "\t" << NAP[i] << "\t" << NAP_err[i] << "\t" << TauAP[i] << "\t" << TauAP_err[i] << "\t" << recFrac[i] << "\t" << recFrac_err[i] << "\t" << tAPtot << "\t" << tAP1us << endl;
	}
	OUTPUT.close();

	delete h1;
	delete c2;
}

void FOLDER::writeTimeDiff2(string file, double bound_lower, double bound_upper, double gain) {
	ofstream OUTPUT;
	OUTPUT.open((filepath+"Auswertung/"+file+".dat").c_str());
	int nr_entries=(Int_t)treeFit->GetEntries();
	for(int i=0; i<nr_entries-1; i++) {
		treeFit->GetEntry(i);
		if(peak_i==1) {
			if(dPulse.FuncAmp>(0.5*gain) && dPulse.FuncAmp<(5.5*gain) && dPulse.FuncTime>bound_lower && dPulse.FuncTime<bound_upper) {
				double time_temp_i = dPulse.FuncTime + dPulse.TriggerTime;
				OUTPUT << dPulse.FuncTime << "\t" << (dPulse.FuncAmp/gain) << endl;
				treeFit->GetEntry(i+1);
				if(!(dPulse.TriggerTime==0.0 && peak_i==1)) {
					if(dPulse.FuncTime>-30) {
						double time_temp_ii = dPulse.FuncTime + dPulse.TriggerTime;
						OUTPUT << (time_temp_ii-time_temp_i) << "\t" << (dPulse.FuncAmp/gain) << endl;
					}
				}
			}
		}
	}
	OUTPUT.close();
}

void FOLDER::makingHistogram(const string filepath, string file_in, string file_out, double binwidth, int column) {
	stringstream TERMINAL;
	TERMINAL.str("");
	TERMINAL << "echo -e -n \""
			<< "reset \n"
			<< "bin(x,width)=width*floor(x/width) \n "
			<< "binwidth= " << binwidth << "\n"
//			<< "set output \\\"| tail +5 | head -n -2 > " << filepath << "Auswertung/" << file_out << ".dat\\\" \n"
			<< "set output \\\"| head -n -2 > " << filepath << "Auswertung/" << file_out << ".dat\\\" \n"
			<< "set table \n"
			<< "plot \\\"" << filepath << "Auswertung/" << file_in << ".dat\\\" u (bin(\\$" << column << ",binwidth)):(1.0) smooth frequency \n"
			<< "unset table \n"
			<< "reset \n"
			<< "\" | gnuplot";
	system(TERMINAL.str().c_str());
	TERMINAL.str("");
}

void FOLDER::make2DHistogramWave() {
	TFile *file = new TFile((filepath+"input.root").c_str(),"READ");
	if(!file->IsOpen() || file->IsZombie()) {
		cout << "Error: cannot open " << filepath << "input.root !\n";
		exit(-1);
	}
	cout << "Opening Data:\r";
	unsigned int length;
	unsigned int nr_waves;
	unsigned int evt_input;

	treeInput = (TTree*) file->Get("treeInput");
	treeInput->SetBranchAddress("evt_n",&evt_input);
	treeInput->SetBranchAddress("length",&length);
	treeInput->SetBranchAddress("time",x);
	treeInput->SetBranchAddress("voltage",voltage);
	nr_waves = (Int_t)treeInput->GetEntries();

	double start2=getTime(0,0);
	//int ybins=128;
	int ybins=295;
	double ymin=-105.0;
	double ymax=25.0;
//	TH2D *h2 = new TH2D("h2","h2",25002,-1000.2,9000.2,ybins,ymin,ymax);
	TH2D *h2 = new TH2D("h2","h2",25002,-1000.2,9000.2,ybins+1,ymin-((ymax-ymin)/ybins/2),ymax+((ymax-ymin)/ybins/2));

	for (unsigned int i = 0; i < nr_waves ; i++) {
		treeInput->GetEntry(i);
		if(treeFit->GetEntries(TString::Format("evt_n ==%i",evt_input))==1) {
			for(unsigned int j=0; j<length; ++j) {
				h2->Fill(x[j],voltage[j]);
			}
		}
//		else{cout << i << "\t" << evt_input << "\t" << treeFit->GetEntries(TString::Format("evt_n ==%i",evt_input)) << endl;}
		loadbar(i, nr_waves, 100, 100, start2);
	}
	ofstream OUTPUT;
	OUTPUT.open((filepath+"Auswertung/persistence.dat").c_str());
	for(int j=0; j<ybins+1; ++j) {
		for (unsigned int i = 0; i < length ; i++) {
			if(h2->GetBinContent(i,j)<5) {OUTPUT << "0" << "\t";}
			else {OUTPUT << h2->GetBinContent(i,j) << "\t";}
		}
		OUTPUT << "\n" ;
	}
	OUTPUT.close();
	delete h2;
}

void FOLDER::findHistogramPeaks(const string filepath, string file_out, unsigned const int counter, vector<Data>* HISTO, Data* max, Data* min) {
	const unsigned int boundary = 5;
	int max_i[counter+1];
	bool bool_max[counter+1];
	int min_i[counter+1];
	bool bool_min[counter+1];

	for(unsigned int k=0; k<=counter; k++)	{
		max[k].value=0.0;
		max[k].x=0.0;
		max_i[k]=0;
		bool_max[k]=false;
		min[k].value=0.0;
		min[k].x=0.0;
		min_i[k]=0;
		bool_min[k]=false;
	}
	ofstream OUTPUT;
	OUTPUT.open((filepath+"Auswertung/"+file_out+".dat").c_str());

	for(unsigned int k=1; k<=counter; k++) {
//		cout << "------" << k << "------" << endl;
		for(unsigned int i=min_i[k-1]+boundary; (i<HISTO->size()-boundary && !bool_max[k]); i++){
			unsigned int j=1;
			while(HISTO->at(i).value>1 && HISTO->at(i).value >= HISTO->at(i+j).value && HISTO->at(i).value >= HISTO->at(i-j).value && j<=boundary)	{
				if(j==boundary)	{
					max[k].value = HISTO->at(i).value;
					max[k].x = HISTO->at(i).x;
					max_i[k] = i;
					bool_max[k]=true;
//					cout << "max:\t" << max_i[k] << "\t" << max[k].x << "\t" << max[k].value << endl;
				}
				j++;
			}
		}
		for(unsigned int i=max_i[k]+1; (i<HISTO->size()-boundary && !bool_min[k]); i++)		{
			unsigned int j=1;
			while(HISTO->at(i).value <= HISTO->at(i+j).value && HISTO->at(i).value <= HISTO->at(i-j).value && j<=boundary) {
				if(j==boundary)	{
					min[k].value = HISTO->at(i).value;
					min[k].x = HISTO->at(i).x;
					min_i[k] = i;
					bool_min[k]=true;
//					cout << "min:\t" << min_i[k] << "\t" << min[k].x << "\t" << min[k].value << endl;
				}
				j++;
			}
		}
		if (max[k].x==0.0) {
			max[k].x=max[k-1].x;
			max[k].value=max[k-1].value;
		}
		if (min[k].x==0.0) {
			min[k].x=min[k-1].x;
			min[k].value=min[k-1].value;
		}
		OUTPUT << max[k].x << "\t" <<  max[k].value << "\t" << min[k].x << "\t" << min[k].value << "\n";
	}
	OUTPUT.close();
}

void FOLDER::writeGnuScriptSpec(const string folder, string file_in, string file_out, Data* max, Data* min, const unsigned int counter) {
	ofstream OUTPUT;
	OUTPUT.open((folder+"Auswertung/"+file_out+"-script.gnu").c_str());
	OUTPUT 	<< "reset \n"
			<< "load \"/home/hpc/capm/sn0515/Master/template.gnu\" \n"
			<< "set xlabel \"maximum [V]\" \n"
			<< "set ylabel \"counts #\" \n"
			<< "g(x,a,b,c)=a*exp(-0.5*(x-b)**2/c**2) \n";
	for(unsigned int i=1; i<=counter ;i++)	{
		OUTPUT << "b"<< i <<"=" << max[i].x <<" \n"
			   << "c"<< i <<"= 1.0\n"
			   << "a"<< i <<"=" << max[i].value <<" \n";
//		for(int j=1; j<=5 ;j++)
//		{
//			OUTPUT	<< "fit [(b"<< i <<"-(2)):(b"<< i <<"+(2))] g(x,a"<< i <<",b"<< i <<",c"<< i <<") \"" << folder << "Auswertung/"+file_out+".dat\" via a"<< i <<",b"<< i <<",c"<< i <<" \n";
//		}
		for(int j=1; j<=5 ;j++)
		{
			OUTPUT	<< "fit [(b"<< i <<"-(2*abs(c"<< i <<"))):(b"<< i <<"+(2*abs(c"<< i <<")))] g(x,a"<< i <<",b"<< i <<",c"<< i <<") \"" << folder << "Auswertung/"+file_out+".dat\" via a"<< i <<",b"<< i <<",c"<< i <<" \n";
		}
	}
	OUTPUT << " \n"
		   << "set print \"" << folder << "Auswertung/"+file_out+"-param.dat\" \n"; //append
	for(unsigned int i=1; i<=counter ;i++) {
		OUTPUT << "print a"<< i <<", a"<< i <<"_err, b"<< i <<",  b"<< i <<"_err, c"<< i <<", c"<< i <<"_err \n";
	}
	OUTPUT <<  "set yr [1:] \n"
			<<  "set log y \n"
			<<  "set xr [0:50] \n"
			<<  "set terminal dumb \n"
			<< "plot \"" << folder << "Auswertung/"+file_out+".dat\" with histeps t \"data\" lt 1\n";
			for(unsigned int i=1; i<=counter ;i++)
			{
				OUTPUT << "replot [b"<< i <<"-2*abs(c"<< i <<"):b"<< i <<"+2*abs(c"<< i <<")] g(x,a"<< i <<",b"<< i <<",c"<< i <<") t \""<< i <<"\" lt 1 lc "<< i+1 <<"\n";
			}
			OUTPUT << "set terminal cairolatex pdf \n"
					<< "set output \"" << folder << "Auswertung/"+file_out+".tex\" \n"
			<< "replot \n"
			<< "set terminal png \n"
			<< "set output \"" << folder << "Auswertung/"+file_out+".png\" \n"
			<< "replot \n"
			<< "set terminal x11 \n"
			<< "replot \n"
			<< "reset \n";
	OUTPUT.close();

	stringstream TERMINAL;
	TERMINAL.str("");
	TERMINAL << "rm -f " << folder << "Auswertung/"+file_out+".pdf";
	system(TERMINAL.str().c_str());
	TERMINAL.str("");
	TERMINAL << "gnuplot " << folder << "Auswertung/"+file_out+"-script.gnu &> " << folder << "Auswertung/"+file_out+"-temp.dat" ;
	system(TERMINAL.str().c_str());
	TERMINAL.str("");
}

void FOLDER::writeGnuScriptGain(const string folder, string file_in, string file_out, int counter) {
	ofstream OUTPUT;
	OUTPUT.open((folder+"Auswertung/"+file_out+"-script.gnu").c_str());
	OUTPUT 	<< "reset \n"
			<< "load \"/home/hpc/capm/sn0515/Master/template.gnu\" \n"
			<< "set xlabel \"photon equivalents [#]\" \n"
			<< "set ylabel \"mean pulse height [mV]\" \n"
			<< "g(x,m,t)=m*x+t \n"
			<< "fit g(x,m,t) \"" << folder << "Auswertung/"+file_in+".dat\" u ($0+1):($3):($4) via m,t \n"
			<< "set print \"" << folder << "Auswertung/"+file_out+"-param.dat\" \n" //append
			<< "print m, m_err, t, t_err \n"
			<< "set yr [0:] \n"
			<< "set xr [0:" << counter+1 << "] \n"
			<< "set terminal cairolatex pdf \n"
			<< "set output \"" << folder << "Auswertung/"+file_out+".tex\" \n"
			<< "plot \"" << folder << "Auswertung/"+file_in+".dat\" u ($0+1):($3):($4) wi yerrorbars t \"data\" lt 1 ,\\\n"
			<< "g(x,m,t) t \"fit\" lt 1 lc 2\n"
			<< "set terminal png \n"
			<< "set output \"" << folder << "Auswertung/"+file_out+".png\" \n"
			<< "replot \n"
			<< "set terminal x11 \n"
			<< "replot \n"
			<< "reset \n";
	OUTPUT.close();

	stringstream TERMINAL;
	TERMINAL.str("");
	TERMINAL << "rm -f " << folder << "Auswertung/"+file_out+".pdf";
	system(TERMINAL.str().c_str());
	TERMINAL.str("");
	TERMINAL << "gnuplot " << folder << "Auswertung/"+file_out+"-script.gnu &> " << folder << "Auswertung/"+file_out+"-temp.dat" ;
	system(TERMINAL.str().c_str());
	TERMINAL.str("");
}

void FOLDER::writeGnuScriptDCR(const string folder, string file_in, string file_out, double bound_lower) {
	ofstream OUTPUT;
	OUTPUT.open((folder+"Auswertung/"+file_out+"-script.gnu").c_str());
	OUTPUT 	<< "reset \n"
			<< "load \"/home/hpc/capm/sn0515/Master/template_big.gnu\" \n"
			<< "set xlabel \"time to next pulse [ns]\" \n"
			<< "set ylabel \"events [#]\" \n"
			<< "f(x)=w*(R/1e9)*exp(-1*R*x/1e9)\n"
			<< "R=30; w=0.5 \n"
			<< "fit [" << bound_lower<< ":1e10] f(x) \"" << folder << "Auswertung/"+file_in+".dat\" u 1:2:($3>0.0 ? $3 : 1/0) yerrors via R,w \n"
			<< "set print \"" << folder << "Auswertung/"+file_out+"-param.dat\" \n" //append
			<< "print R, R_err, w, w_err \n"
			<< "set yr [1e-12:1e-7] \n"
			<< "set xr [1e5:1e10] \n"
			<< "set log xy \n"
			<< "set terminal png \n"
			<< "set output \"" << folder << "Auswertung/"+file_out+".png\" \n"
			<< "plot \"" << folder << "Auswertung/"+file_in+".dat\" u 1:2:3 wi errorbars t \"data\" pt 7 ps 0.5 ,\\\n"
			<< "f(x) t \"fit\" lt 1 lc 2\n"
			<< "set terminal x11 \n"
			<< "replot \n"
			<< "reset \n";
	OUTPUT.close();

	stringstream TERMINAL;
	TERMINAL.str("");
	TERMINAL << "rm -f " << folder << "Auswertung/"+file_out+".pdf";
	system(TERMINAL.str().c_str());
	TERMINAL.str("");
	TERMINAL << "gnuplot " << folder << "Auswertung/"+file_out+"-script.gnu &> " << folder << "Auswertung/"+file_out+"-temp.dat" ;
	system(TERMINAL.str().c_str());
	TERMINAL.str("");
}

void FOLDER::writePulseDataRec(string file, double gain, double chi2_max=10.0) {
	ofstream OUTPUT;
	OUTPUT.open((filepath+"Auswertung/"+file+".dat").c_str());
	int nr_entries=(Int_t)treeFit->GetEntries();
//	for(int i=0; i<nr_entries; i++) {
//		treeFit->GetEntry(i);
//		if (peak_i==1 || peak_i==2) {
//			if(dPulse.FitChi2/dPulse.FitNDF<chi2_max) {
//				if(peak_i==1){OUTPUT << dPulse.FuncTime << "\t" << "0.0" << "\t" <<  dPulse.FuncAmp/gain  << "\t" << "1" <<  endl;}
//				if(peak_i==2){
//					double time_ii = dPulse.FuncTime;
//					double amp = dPulse.FuncAmp;
//					treeFit->GetEntry(i-1);
//					double time_i=dPulse.FuncTime;
//					OUTPUT << time_ii << "\t" << time_ii-time_i << "\t" <<  amp/gain << "\t" << "2" <<  endl;
//				}
//			}
//		}
//	}

	for(int i=0; i<nr_entries-1; i++) {
		treeFit->GetEntry(i);
		if(peak_i==1) {
			if(!(dPulse.TriggerTime==0.0)) {
				double time_temp_ii = dPulse.FuncTime + dPulse.TriggerTime;
				treeFit->GetEntry(i+1);
				double time_temp_i = dPulse.FuncTime + dPulse.TriggerTime;
				OUTPUT << time_temp_i-time_temp_ii << "\t" << dPulse.FuncAmp/gain << endl;
			}
		}
	}


	OUTPUT.close();
}

void FOLDER::writeGnuScriptRecTime(const string folder, string file_in, string file_out, double bound, int counter) {
	ofstream OUTPUT;
	OUTPUT.open((folder+"Auswertung/"+file_out+"-script.gnu").c_str());
	OUTPUT 	<< "reset \n"
			<< "load \"/home/hpc/capm/sn0515/Master/template_big.gnu\" \n"
			<< "set xlabel \"time[ns]\" \n"
			<< "set ylabel \"pulse height [mV]\" \n"
			<< "f(x,t,c)=(1-exp(-1*(x/t)))+c \n";
	for(int i=1; i<=counter ;i++)	{
		OUTPUT << "t"<< i <<"= 60\n"
			   << "c"<< i <<"= "<< i-1 <<"\n";
		for(int j=1; j<=5 ;j++)
		{
			OUTPUT	<< "fit f(x,t"<< i <<",c"<< i <<") \"" << folder << "Auswertung/"+file_in+".dat\" u 1:(($2<f($1,t1,c"<< i <<"+"<< bound <<") && $2>f($1,t1,-"<< bound <<"+c"<< i <<")) ? $2 : 1/0) via t"<< i <<"\n";
		}
	}
	OUTPUT << " \n"
		   << "set print \"" << folder << "Auswertung/"+file_out+"-param.dat\" \n"; //append
	for(int i=1; i<=counter ;i++) {
		OUTPUT << "print t"<< i <<", t"<< i <<"_err\n";
	}
	OUTPUT << "set yr [0:4.5] \n"
			<< "set xr [1:1e10] \n"
			<< "set log x \n"
			<< "set terminal dumb \n"
			<< "plot \"" << folder << "Auswertung/"+file_in+".dat\" u 1:(($2<f($1,t1,"<< bound <<") && $2>f($1,t1,-"<< bound <<")) ? $2 : 1/0) w d t \"data\" lt 1 \n"
			<< "replot \"" << folder << "Auswertung/"+file_in+".dat\" u 1:(($2<f($1,t1,"<< bound <<") && $2>f($1,t1,-"<< bound <<")) ? 1/0 : $2) w d notitle lt 1 lc -1 \n";
	for(int i=1; i<=counter ;i++)
	{
		OUTPUT << "replot f(x,t"<< i <<",c"<< i <<") t \"fit\" lt 1 lc 2,\\\n"
			<< "f(x,t"<< i <<",c"<< i <<"+"<<bound<<") t \"fit-limit\" lt 1 lc 3,\\\n"
			<< "f(x,t"<< i <<",c"<< i <<"-"<<bound<<") t \"fit-limit\" lt 1 lc 3 \n";
	}
	OUTPUT << "set terminal pdf \n"
			<< "set output \"" << folder << "Auswertung/"+file_out+".pdf\" \n"
			<< "replot \n"
			<< "set terminal png \n"
			<< "set output \"" << folder << "Auswertung/"+file_out+".png\" \n"
			<< "replot \n"
			<< "set terminal x11 \n"
			<< "replot \n"
			<< "reset \n";
	OUTPUT.close();

	stringstream TERMINAL;
	TERMINAL.str("");
	TERMINAL << "rm -f " << folder << "Auswertung/"+file_out+".pdf";
	system(TERMINAL.str().c_str());
	TERMINAL.str("");
	TERMINAL << "gnuplot " << folder << "Auswertung/"+file_out+"-script.gnu &> " << folder << "Auswertung/"+file_out+"-temp.dat" ;
	system(TERMINAL.str().c_str());
	TERMINAL.str("");
}

void FOLDER::make2DHistogram(string file_in, string file_out, double xlow, double xhigh, int xbins, double ylow, double yhigh, int ybins) {
	double x,y;
	TH2D *h2 = new TH2D("h2","h2",xbins,xlow,xhigh,ybins,ylow,yhigh);

	ifstream INPUT;
	INPUT.open((filepath+"Auswertung/"+file_in+".dat").c_str());
	if(!INPUT.is_open()) {
		cout << "Error: cannot open " << filepath << "Auswertung/" << file_in << ".dat !" << endl;
	}
	double dummy;
	double check;
	while(!INPUT.eof()) {
		INPUT >> x;
		INPUT >> y;
		INPUT >> check;
		INPUT >> dummy;
		INPUT >> dummy;
		if(check != 1.0 && check!=2.0 && check!=3.0) {
			h2->Fill(x,y);
		}
	}
	INPUT.close();

	ofstream OUTPUT;
	OUTPUT.open((filepath+"Auswertung/"+file_out+".dat").c_str());
	for(int j=0; j < ybins; ++j) {
		for (int i = 0; i < xbins ; i++) {
			OUTPUT << h2->GetBinContent(i,j) << "\t";
		}
		OUTPUT << "\n" ;
	}
	OUTPUT.close();
	delete h2;
}

void FOLDER::writeGnuScript2DHisto(const string folder, string file_in, string file_out, double xlow, double xhigh, int xbins, double ylow, double yhigh, int ybins) {
	ofstream OUTPUT;
	OUTPUT.open((folder+"Auswertung/"+file_out+"-script.gnu").c_str());
	OUTPUT 	<< "reset \n"
			<< "load \"/home/hpc/capm/sn0515/Master/template.gnu\" \n"
			<< "unset key \n"
			<< "unset grid \n"
			<< "set xlabel \"time [ns]\" \n"
			<< "set ylabel \"pulse height [mV]\" \n"
			<< "set yr [" << ylow << ":" << yhigh << "] \n"
			<< "set xr [" << xlow << ":" << xhigh << "] \n"
			<< "set log cb\n"
			<< "set pal mod RGB funct (gray>0.?sqrt(gray):1), (gray>0.?gray**3:1), (gray>0.?sin(gray*2.*pi):1) \n"
			<< "set terminal pdf \n"
			<< "set output \"" << folder << "Auswertung/"+file_out+".pdf\" \n"
			<< "plot \"" << folder << "Auswertung/"+file_in+".dat\" matrix u (($1*" << (xhigh-xlow)/xbins << ")+" << xlow << "):(($2*" << (yhigh-ylow)/ybins << ")+" << ylow << "):3 wi image \n"
			<< "set terminal png \n"
			<< "set output \"" << folder << "Auswertung/"+file_out+".png\" \n"
			<< "replot \n"
			<< "set terminal x11 \n"
			<< "replot \n"
			<< "reset \n";
	OUTPUT.close();

	stringstream TERMINAL;
	TERMINAL.str("");
	TERMINAL << "rm -f " << folder << "Auswertung/"+file_out+".pdf";
	system(TERMINAL.str().c_str());
	TERMINAL.str("");
	TERMINAL << "gnuplot " << folder << "Auswertung/"+file_out+"-script.gnu &> " << folder << "Auswertung/"+file_out+"-temp.dat" ;
	system(TERMINAL.str().c_str());
	TERMINAL.str("");
}

//void DATA::fitGnuMultigauss(const string folder, Data* max, Data* min, unsigned const int counter) {
//	stringstream TERMINAL;
//	TERMINAL << "rm -f " << folder << "Auswertung/temp/Hist_spectrum-par.dat";
//	system(TERMINAL.str().c_str());
//	TERMINAL.str("");
//	//***********************************Create Spectrum Multigauss************************************
//	TERMINAL << "echo -e -n \""
//			<< "set fit errorvariables \n"
//			<< "set sample 5000\n"
//			<< "FIT_LIMIT = 1e-6\n"
//			<< "FIT_MAXITER=1000\n"
//			<< "g(x,a,b,c)=a*exp(-0.5*(x-b)**2/c**2) \n";
//	for(unsigned int i=1; i<=counter ;i++)
//	{
//		TERMINAL << "b"<< i <<"=" << max[i].x <<" \n"
//				<< "c"<< i <<"=" << (min[i].x-min[i-1].x)/3 <<" \n"
//				<< "a"<< i <<"=" << max[i].value <<" \n";
//		for(int j=1; j<=5 ;j++)
//		{
//			TERMINAL << "fit ["<< min[i-1].x <<":"<< min[i].x <<"] g(x,a"<< i <<",b"<< i <<",c"<< i <<") \\\"" << folder << "Auswertung/temp/Hist_spectrum.dat\\\" via a"<< i <<",b"<< i <<",c"<< i <<" \n";
//		}
//		for(int j=1; j<=5 ;j++)
//		{
//			TERMINAL << "fit [(b-(2)):(b+(2))] g(x,a"<< i <<",b"<< i <<",c"<< i <<") \\\"" << folder << "Auswertung/temp/Hist_spectrum.dat\\\" via b"<< i <<",a"<< i <<" \n"
//					<< "fit [(b-(2)):(b+(2))] g(x,a"<< i <<",b"<< i <<",c"<< i <<") \\\"" << folder << "Auswertung/temp/Hist_spectrum.dat\\\" via b"<< i <<" \n"
//					<< "fit [(b-(2)):(b+(2))] g(x,a"<< i <<",b"<< i <<",c"<< i <<") \\\"" << folder << "Auswertung/temp/Hist_spectrum.dat\\\" via a"<< i <<" \n"
//					<< "fit [(b-(2)):(b+(2))] g(x,a"<< i <<",b"<< i <<",c"<< i <<") \\\"" << folder << "Auswertung/temp/Hist_spectrum.dat\\\" via c"<< i <<" \n"
//					<< "fit [(b-(2)):(b+(2))] g(x,a"<< i <<",b"<< i <<",c"<< i <<") \\\"" << folder << "Auswertung/temp/Hist_spectrum.dat\\\" via a"<< i <<",b"<< i <<",c"<< i <<" \n";
//		}
//		for(int j=1; j<=5 ;j++)
//		{
//			TERMINAL << "fit [(b-(c)):(b+(c))] g(x,a"<< i <<",b"<< i <<",c"<< i <<") \\\"" << folder << "Auswertung/temp/Hist_spectrum.dat\\\" via b"<< i <<",a"<< i <<" \n"
//					<< "fit [(b-(c)):(b+(c))] g(x,a"<< i <<",b"<< i <<",c"<< i <<") \\\"" << folder << "Auswertung/temp/Hist_spectrum.dat\\\" via b"<< i <<" \n"
//					<< "fit [(b-(c)):(b+(c))] g(x,a"<< i <<",b"<< i <<",c"<< i <<") \\\"" << folder << "Auswertung/temp/Hist_spectrum.dat\\\" via a"<< i <<" \n"
//					<< "fit [(b-(c)):(b+(c))] g(x,a"<< i <<",b"<< i <<",c"<< i <<") \\\"" << folder << "Auswertung/temp/Hist_spectrum.dat\\\" via c"<< i <<" \n"
//					<< "fit [(b-(c)):(b+(c))] g(x,a"<< i <<",b"<< i <<",c"<< i <<") \\\"" << folder << "Auswertung/temp/Hist_spectrum.dat\\\" via a"<< i <<",b"<< i <<",c"<< i <<" \n";
//		}
//		for(int j=1; j<=5 ;j++)
//		{
//			TERMINAL << "fit [(b-(c)):(b+(c))] g(x,a"<< i <<",b"<< i <<",c"<< i <<") \\\"" << folder << "Auswertung/temp/Hist_spectrum.dat\\\" via a"<< i <<",b"<< i <<",c"<< i <<" \n";
//		}
//	}
////	TERMINAL << "multiG(x)=g(x,a1,b1,c1)";
////	for(unsigned int i=2; i<=counter ;i++)
////	{
////		TERMINAL << "+g(x,a"<< i <<",b"<< i <<",c"<< i <<")";
////	}
////	TERMINAL << " \n"
////			<< "fit [0:"<< min[counter].x <<"] multiG(x) \\\"" << folder << "Auswertung/temp/Hist_spectrum.dat\\\" via a1,b1,c1";
////	for(unsigned int i=2; i<=counter ;i++)
////	{
////		TERMINAL << ",a"<< i <<",b"<< i << ",c"<< i << "";
////	}
//	TERMINAL << " \n"
//			<< "set print \\\"" << folder << "Auswertung/temp/Hist_spectrum-par.dat\\\" append \n";
//	for(unsigned int i=1; i<=counter ;i++)
//	{
//		TERMINAL << "print a"<< i <<", a"<< i <<"_err, b"<< i <<",  b"<< i <<"_err, c"<< i <<", c"<< i <<"_err \n";
//	}
//	TERMINAL <<  "set yr [1:] \n"
//			<<  "set log y \n"
//			<<  "set terminal pdf \n"
//			<< "set output \\\"" << folder << "Auswertung/Hist_spectrum.pdf\\\" \n"
//			<< "set samples = 10000\n"
//			<< "plot \\\"" << folder << "Auswertung/temp/Hist_spectrum.dat\\\" with histeps \n";
//			for(unsigned int i=1; i<=counter ;i++)
//			{
//				TERMINAL << "replot [b"<< i <<"-2*c"<< i <<":b"<< i <<"+2*c"<< i <<"] g(x,a"<< i <<",b"<< i <<",c"<< i <<") lw 2\n";
//			}
//			TERMINAL << "\" | gnuplot&>" << folder << "Auswertung/temp/Hist_spectrum-fit.dat";
//	system(TERMINAL.str().c_str());
//
////	ifstream INPUT;
////	INPUT.open((folder+"Auswertung/temp/Hist_spectrum-par.dat").c_str());
////	double pos[counter];
////	double pos_err[counter];
////	double dummy;
////	for(unsigned int i=0; i<counter; i++)
////	{
////		INPUT >> dummy;
////		INPUT >> dummy;
////		INPUT >> pos[i];
////		INPUT >> pos_err[i];
////		INPUT >> dummy;
////		INPUT >> dummy;
////	}
//}

void FOLDER::readGnuHistogram(string filepath, string file, vector<Data>* histo) {
	string dummy;
	ifstream INPUT_HISTO;
	INPUT_HISTO.open((filepath+file).c_str());
	if(INPUT_HISTO.is_open()) {
		while(!INPUT_HISTO.eof())	{
			Data a;
			INPUT_HISTO >> a.x;
			INPUT_HISTO >> a.value;
			INPUT_HISTO >> dummy;
			if(dummy=="u") {break;};
			histo->push_back(a);
		}
	}
	else {cout<<"histo could not be opened" << endl;}
	INPUT_HISTO.close();
	histo->pop_back();
//	cout << histo->size() << endl;
//	for (int i=0; i<histo->size(); i++) {
//	cout << histo->at(i).x << "\t" << histo->at(i).value << endl;
//	}
}

void FOLDER::calcFitParameter() {
	TCanvas *c2 = new TCanvas("c2","c2",600,400);
	TH1D *RiseTime = new TH1D("RiseTime","RiseTime histogram",400,0.0,20.0);
	TH1D *FallTime = new TH1D("FallTime","FallTime histogram",400,0.0,400.0);
//	vector<double> rise;
//	vector<double> fall;

	int nr_entries=(Int_t)treeFit->GetEntries();
	for(int i=0; i<nr_entries; i++) {
		treeFit->GetEntry(i);
		if((std::isfinite(dPulse.FitAmp) &&
				std::isfinite(dPulse.FitTime) &&
				std::isfinite(dPulse.FitBaseline) &&
				std::isfinite(dPulse.FitRiseTime) &&
				std::isfinite(dPulse.FitFallTime) &&
				std::isfinite(dPulse.FitChi2) &&
				std::isfinite(dPulse.FitNDF) &&
				std::isfinite(dPulse.Baseline) &&
				std::isfinite(dPulse.Time) &&
				std::isfinite(dPulse.Amp) &&
				std::isfinite(dPulse.FitLowLimit) &&
				std::isfinite(dPulse.FitHighLimit))) {
//				std::isfinite(dPulse.FuncTime) &&
//				std::isfinite(dPulse.FuncAmp))) {
			RiseTime->Fill(dPulse.FitRiseTime);
			FallTime->Fill(dPulse.FitFallTime);
//			rise.push_back(dPulse.FitRiseTime);
//			fall.push_back(dPulse.FitFallTime);
		}
	}
	TF1 *funcRise = new TF1("funcRise","gaus");
	TF1 *funcFall = new TF1("funcFall","gaus");
	funcRise->SetNpx(10000);
	funcFall->SetNpx(10000);

	funcRise->SetParameter(1,2.0);
	funcRise->SetParameter(2,1.5);
	funcFall->SetParameter(1,80.0);
	funcFall->SetParameter(2,20.0);

	ofstream OUTPUT1;
	OUTPUT1.open((filepath+"Auswertung/fit2/histo_Rise.dat").c_str());
	for(int i=0; i<400 ; ++i) {
		OUTPUT1 << RiseTime->GetBinCenter(i) << "\t"<< RiseTime->GetBinContent(i) << endl;
	}
	OUTPUT1.close();
	ofstream OUTPUT2;
	OUTPUT2.open((filepath+"Auswertung/fit2/histo_Fall.dat").c_str());
	for(int i=0; i<400 ; ++i) {
		OUTPUT2 << FallTime->GetBinCenter(i) << "\t"<< FallTime->GetBinContent(i) << endl;
	}
	OUTPUT2.close();

	double sigma=1.0;
	int repeat=5;

	for(int i=0; i<repeat; i++) {
		funcRise->SetRange(funcRise->GetParameter(1)-(sigma*funcRise->GetParameter(2)),funcRise->GetParameter(1)+(sigma*funcRise->GetParameter(2)));
		funcFall->SetRange(funcFall->GetParameter(1)-(sigma*funcFall->GetParameter(2)),funcFall->GetParameter(1)+(sigma*funcFall->GetParameter(2)));
		if(i!=repeat-1) {
			RiseTime->Fit(funcRise,"EMQR0");
			FallTime->Fit(funcFall,"EMQR0");
		}
		else {
			RiseTime->Fit(funcRise,"EMQR");
			FallTime->Fit(funcFall,"EMQR");
		}
	}
	RiseTime->Draw();
	c2->SaveAs((filepath+"Auswertung/fit2/RiseTime-Histogram.pdf").c_str());
	FallTime->Draw();
	c2->SaveAs((filepath+"Auswertung/fit2/FallTime-Histogram.pdf").c_str());

	cout << endl << endl;
	cout << "Entries in TTree Fit:\t" << nr_entries << endl;
	cout << "RiseTime:\t" << funcRise->GetParameter(1) << "\t +- \t" << funcRise->GetParError(1) << "\t\tsigma:\t" << funcRise->GetParameter(2) << endl;
	cout << "FallTime:\t" << funcFall->GetParameter(1) << "\t +- \t" << funcFall->GetParError(1) << "\t\tsigma:\t" << funcFall->GetParameter(2) << endl;

//	double rise_mean=funcRise->GetParameter(1);
//	double fall_mean=funcFall->GetParameter(1);
//	double rise_sig=0.0;
//	double fall_sig=0.0;
//	int counter_rise=0;
//	int counter_fall=0;
//	for(int i=0; i<rise.size(); i++) {
//		if(rise.at(i)>0.3 && rise.at(i)< 19.8) {
//			rise_sig += pow((rise.at(i)-rise_mean),2.0);
//			counter_rise++;
//		}
//		if(fall.at(i)>9.0 && fall.at(i)<398.0) {
//			fall_sig += pow((fall.at(i)-fall_mean),2.0);
//			counter_fall++;
//		}
//	}
//	rise_sig=sqrt(rise_sig/(counter_rise-1));
//	fall_sig=sqrt(fall_sig/(counter_fall-1));


	ofstream OUTPUT;
	OUTPUT.open((filepath+"Auswertung/fit2/parameter.dat").c_str());
	OUTPUT << "#NrEntries\tRiseTime\tRiseTime_sig\tFallTime\tFallTime_sig\n";
	OUTPUT << nr_entries << "\t"
			<< funcRise->GetParameter(1) << "\t"
			<< funcRise->GetParError(1) << "\t"
			<< funcRise->GetParameter(2) << "\t"
			<< funcFall->GetParameter(1) << "\t"
			<< funcFall->GetParError(1) << "\t"
			<< funcFall->GetParameter(2) << endl;
	OUTPUT.close();

	delete funcRise;
	delete funcFall;
	delete RiseTime;
	delete FallTime;
	delete c2;
}

void FOLDER::calcCrosstalkLinear(string filepath, string file_in, double pe1 , double pe2) {
	double peak;
	int counter1=0;
	int counterM=0;
	double dummy;
	double time;

	ifstream INPUT;
	INPUT.open((filepath+"Auswertung/"+file_in+".dat").c_str());
	if(!INPUT.is_open()) {
		cout << "Error: cannot open " << filepath << "Auswertung/" << file_in << ".dat !" << endl;
	}

	double limit = (pe1+pe2)/2.0;
//	cout << "limit: " << limit << endl;

	while(!INPUT.eof()) {
		INPUT >> time;
		INPUT >> peak;
		INPUT >> dummy;
		INPUT >> dummy;
		INPUT >> dummy;
//		if(time >-20 && time < 20) {
			if(peak>0.0 && peak < 100.0) {
				if(peak<limit) {counter1++;}
				else {counterM++;}
			}
//		}
	}
	INPUT.close();

	cout << "counter1: " << counter1 << endl;
	cout << "counterM: " << counterM << endl;

	cout << endl << "Crosstalk (linear) is:\t\t" << (100.0*counterM)/(counter1+counterM) << "\t%"<< endl;
}

void FOLDER::calcCrosstalkQuadratic(string filepath, string file_in, double gain) {
	TH1D *peakhis = new TH1D("peakhis","pulsheighthisto",500,0,100);  //histogramm aus analysis.fit

	double dummy;
	double peak;
	double time;

	ifstream INPUT;
	INPUT.open((filepath+"Auswertung/"+file_in+".dat").c_str());
	if(!INPUT.is_open()) {
		cout << "Error: cannot open " << filepath << "Auswertung/" << file_in << ".dat !" << endl;
	}

	while(!INPUT.eof()) {
		INPUT >> time;
		INPUT >> peak;
		INPUT >> dummy;
		INPUT >> dummy;
		INPUT >> dummy;
		if(time >-20 && time < 20) {
			if(peak>0.0 && peak < 100.0) {
				peakhis->Fill(peak);
			}
		}
	}
	INPUT.close();

	cout << "Crosstalk (quadratic) is:\t"<< 100*((peakhis->GetMean()/gain)-1) << "\t%"<< endl;

	delete peakhis;
}

void FOLDER::analysisRun() {
	/*
	 * variable bin size:
	 * https://root.cern.ch/root/roottalk/roottalk00/2318.html
	 */

//	TCanvas *c2 = new TCanvas("c2","c2");
//	c2->SetLogy(1);
//	data.dummyHist->Draw();
//	data.maxima->SetLineColor(kRed);
//	data.maxima->Draw("SAME");
//	c2->SaveAs((filepath+"Auswertung/Histogram_test.pdf").c_str());
//	c2->Clear();
//	delete c2;
//
//	TSpectrum *spec = new TSpectrum(20);
//	Int_t nfound = spec->Search(data.maxima,2.5);
//
//	Double_t par[3*nfound];
//	const Double_t *par_err;
//	Double_t *xpeaks = spec->GetPositionX();
//	std::sort(xpeaks,xpeaks+nfound);
//	char multigauss[255] = "gaus(0)";
//
//	for (int p=0;p<nfound;p++) {
//		Double_t lower_bound;
//		Double_t upper_bound;
//		Int_t bin = data.maxima->GetXaxis()->FindBin(xpeaks[p]);
//		Double_t yp = data.maxima->GetBinContent(bin);
//		if(p!=0) {lower_bound=(xpeaks[p]+xpeaks[p-1])/2;}
//		else {lower_bound=(3*xpeaks[p]-xpeaks[p+1])/2;}
//		if(p!=(nfound-1)) {upper_bound=(xpeaks[p+1]+xpeaks[p])/2;}
//		else {upper_bound=(3*xpeaks[p]-xpeaks[p-1])/2;}
//		par[3*p] = yp;
//		par[3*p+1] = xpeaks[p];
//		par[3*p+2] = (upper_bound-lower_bound)/5;
//
//		TF1 *g1 = new TF1("test","gaus",lower_bound,upper_bound);
//		data.maxima->Fit(g1,"Q0R");
//		g1->GetParameters(&par[3*p]);
//		delete g1;
//
//		if(p!=0) {sprintf(multigauss, "%s+gaus(%i)", multigauss, 3*p);}
//	}
//	Double_t sig=2.0;
//	Double_t lower_bound=par[1]-(sig*par[2]);
//	Double_t upper_bound=par[3*nfound-2]+(sig*par[3*nfound-1]);
//
//	TF1 *total = new TF1("total",multigauss,lower_bound,upper_bound);
//	total->SetParameters(par);
//	total->SetNpx(1000);
//	data.maxima->Fit(total,"Q0R");
//	TCanvas *c1 = new TCanvas("c1","c1");
//	c1->SetLogy(1);
//	data.maxima->Draw();
//	total->Draw("SAME");
//	c1->SaveAs((filepath+"Auswertung/Histogram_log.pdf").c_str());
//	c1->Clear();
//	c1->SetLogy(0);
//	c1->SaveAs((filepath+"Auswertung/Histogram.pdf").c_str());
//	c1->Clear();
//
//	total->GetParameters(&par[0]);
//	par_err=total->GetParErrors();
//	double pe[nfound], pe_err[nfound];
//	double pos[nfound], pos_err[nfound];
//	double fwhm[nfound], fwhm_err[nfound];
//	double res[nfound], res_err[nfound];
//	for(int i=0; i<nfound; i++) {
//		pe[i]=i+1;
//		pe_err[i]=0.0;
//		pos[i]=par[3*i+1];
//		pos_err[i]=par_err[3*i+1];
//		fwhm[i]=FWHM*par[3*i+2];
//		fwhm_err[i]=FWHM*par_err[3*i+2];
//		res[i]=FWHM*par[3*i+2]/par[3*i+1];
//		res_err[i]=res[i]*sqrt(pow(par_err[3*i+2]/par[3*i+2],2.0)+pow(par_err[3*i+1]/par[3*i+1],2.0));
//	}
//	TGraphErrors *gr1=new TGraphErrors(nfound,pe,pos,pe_err,pos_err);
//	TF1 *fline = new TF1("fline","pol1");
//	gr1->Fit("fline","Q");
//	data.gain=fline->GetParameter(1);
//	data.gain_err=fline->GetParError(1);
//	gr1->Draw();
//	c1->SaveAs((filepath+"Auswertung/Gain.pdf").c_str());
//	c1->Clear();
//
//	if(info->bool_verbose) {
//		//Histogram of all calculated baselines
//		int bxmin=data.baselines->FindFirstBinAbove(10.0);
//		double xmin=data.baselines->GetXaxis()->GetBinCenter(bxmin);
//		int bxmax=data.baselines->FindLastBinAbove(10.0);
//		double xmax=data.baselines->GetXaxis()->GetBinCenter(bxmax);
//		data.baselines->GetXaxis()->SetLimits(xmin,xmax);
//		data.baselines->Draw();
//		c1->SaveAs((filepath+"Auswertung/baselines.pdf").c_str());
//		c1->Clear();
//
//		//Histogram of all calculated thresholds
//		int bxmin2=data.thresholds->FindFirstBinAbove(10.0);
//		double xmin2=data.thresholds->GetXaxis()->GetBinCenter(bxmin2);
//		int bxmax2=data.thresholds->FindLastBinAbove(10.0);
//		double xmax2=data.thresholds->GetXaxis()->GetBinCenter(bxmax2);
//		data.thresholds->GetXaxis()->SetLimits(xmin2,xmax2);
//		data.thresholds->Draw();
//		c1->SaveAs((filepath+"Auswertung/thresholds.pdf").c_str());
//		c1->Clear();
//
//		//Plot of Peak FWHM
//		TGraphErrors *gr2=new TGraphErrors(nfound,pos,fwhm,pos_err,fwhm_err);
//		gr2->Draw();
//		c1->SaveAs((filepath+"Auswertung/FWHM.pdf").c_str());
//		c1->Clear();
//		delete gr2;
//
//		//Plot of Peak Resolution=FWHM/Pos
//		TGraphErrors *gr3=new TGraphErrors(nfound,pos,res,pos_err,res_err);
//		gr3->Draw();
//		c1->SaveAs((filepath+"Auswertung/Resolution.pdf").c_str());
//		c1->Clear();
//		delete gr3;
//	}
//	if(info->bool_dark) {
//		Double_t integral_all=data.maxima->Integral(1,data.maxima->GetSize()-1);
//		Double_t integral_left=data.maxima->Integral(1,data.maxima->FindFixBin((par[1]+par[4])/2,1));
//		Double_t prob_one=double(integral_left/integral_all);
//		data.crosstalk=1E2*(1.0-prob_one);
//		data.crosstalk_err=1E2*sqrt((prob_one+pow(prob_one,2.0))/integral_all);
//	}
//
//	delete spec;
//	delete total;
//	delete fline;
//	delete c1;
//	delete gr1;
}

unsigned int FOLDER::getEventLength(unsigned int i) const {
	treeInput->GetEntry(i);
	return length;
}

double FOLDER::getEventTrigger(unsigned int i) const {
	treeInput->GetEntry(i);
	return trigger;
}

double* FOLDER::getEventTime(unsigned int i) {
	treeInput->GetEntry(i);
	return x;
}

double* FOLDER::getEventValue(unsigned int i) {
	treeInput->GetEntry(i);
	return voltage;
}

string FOLDER::getFilepath() const {
	return filepath;
}

string FOLDER::getFile(unsigned int m) const {
	return files.at(m);
}

MAININPUT* FOLDER::getInfo_ptr() const {
	return info;
}

DATA* FOLDER::getData_ptr() {
	return &data;
}

TTree* FOLDER::getTreeInput() {
	return treeInput;
}

vector<string> FOLDER::getFolder() const {
	return files;
}

unsigned int FOLDER::getNumberFiles() const {
	return files.size();
}

/*
void FOLDER::writeData() {
	data.writeData(filepath);
}

void FOLDER::writeDataTree() {
	data.writeDataTree(filepath);
}
*/

void FOLDER::writeOutput() {
	if ( info->type == 4 ){
		data.writeDataTree(filepath);
	}
	else {
		data.writeData(filepath);
	}
}

void FOLDER::printInfo() {
	cout << endl;
	if(WAVE::counter!=0) {
		cout << "number of waves:\t\t" << WAVE::counter << endl;
		//cout << "avg. number of peaks per wave:\t" << (double(data.maxima->Integral(1,data.maxima->GetSize()-1))/WAVE::counter) << endl;
//		if(info->bool_run) {
//			cout << "Gain [mV]:\t\t\t" << data.gain << "  +-  " << data.gain_err << endl;
//			if(info->bool_dark) {
//				cout << "Crosstalk [%]:\t\t\t" << data.crosstalk << "  +-  " << data.crosstalk_err << endl;
//			}
//		}
	}
//	else {
//		cout << "!! empty Folder !! Unable to read Waves !!" << endl;
//	}
	cout << "===================================================================" << endl;
}

void FOLDER::readRootFit_ChangeTree() {
//	double start2=getTime(0,0);
//	TFile *file = TFile::Open((filepath+"analysisFit.root").c_str(),"READ");
//	if(!file->IsOpen()) {
//		cout << "Error: cannot open " << filepath << "analysisFit.root !\n";
//		exit(-1);
//	}
//	cout << "Opening Analysis:\r";
//	treeFit = (TTree*) file->Get("treeFit");
//	treeFit->SetBranchAddress("evt_n",&evt_n);
//	treeFit->SetBranchAddress("peaks_n",&peaks_n);
//	treeFit->SetBranchAddress("peak_i",&peak_i);
//	treeFit->SetBranchAddress("TriggerTime",&dPulse.TriggerTime);
//	treeFit->SetBranchAddress("Baseline",&dPulse.Baseline);
//	treeFit->SetBranchAddress("Amp",&dPulse.Amp);
//	treeFit->SetBranchAddress("Time",&dPulse.Time);
//	treeFit->SetBranchAddress("FitLowLimit",&dPulse.FitLowLimit);
//	treeFit->SetBranchAddress("FitHighLimit",&dPulse.FitHighLimit);
//	treeFit->SetBranchAddress("FitBaseline",&dPulse.FitBaseline);
//	treeFit->SetBranchAddress("FitTime",&dPulse.FitTime);
//	treeFit->SetBranchAddress("FitAmp",&dPulse.FitAmp);
//	treeFit->SetBranchAddress("FitRiseTime",&dPulse.FitRiseTime);
//	treeFit->SetBranchAddress("FitFallTime",&dPulse.FitFallTime);
//	treeFit->SetBranchAddress("FitChi2",&dPulse.FitChi2);
//	treeFit->SetBranchAddress("FitNDF",&dPulse.FitNDF);
//	treeFit->SetBranchAddress("FuncTime",&dPulse.FuncTime);
//	treeFit->SetBranchAddress("FuncAmp",&dPulse.FuncAmp);
//
////	cout << endl << endl << (Int_t)treeFit->GetEntries() << endl << endl;
////	treeFit->Scan("*");
////	exit(-1);
//
//
//	TFile *out2 = new TFile((filepath+"analysisFit1.root").c_str(),"RECREATE");
//	if(!out2->IsOpen() || out2->IsZombie()) {
//		cout << "Error: cannot open " << filepath << "analysisFit1.root !" << endl;
//	}
//	tPulse tPulse;
//	unsigned int evt;
//	unsigned int pea_n;
//	unsigned int pea_i;
//	TTree *treeFit1 = new TTree("treeFit","treeFit");
//	treeFit1->Branch("evt_n",&evt,"evt_n/i");
//	treeFit1->Branch("peaks_n",&pea_n,"peaks_n/i");
//	treeFit1->Branch("peak_i",&pea_i,"peak_i/i");
//	treeFit1->Branch("TriggerTime",&tPulse.TriggerTime,"TriggerTime/D");
//	treeFit1->Branch("Baseline",&tPulse.Baseline,"Baseline/D");
//	treeFit1->Branch("Amp",&tPulse.Amp,"Amp/D");
//	treeFit1->Branch("Time",&tPulse.Time,"Time/D");
//	treeFit1->Branch("FitLowLimit",&tPulse.FitLowLimit,"FitLowLimit/D");
//	treeFit1->Branch("FitHighLimit",&tPulse.FitHighLimit,"FitHighLimit/D");
//	treeFit1->Branch("FitBaseline",&tPulse.FitBaseline,"FitBaseline/D");
//	treeFit1->Branch("FitTime",&tPulse.FitTime,"FitTime/D");
//	treeFit1->Branch("FitAmp",&tPulse.FitAmp,"FitAmp/D");
//	treeFit1->Branch("FitRiseTime",&tPulse.FitRiseTime,"FitRiseTime/D");
//	treeFit1->Branch("FitFallTime",&tPulse.FitFallTime,"FitFallTime/D");
//	treeFit1->Branch("FitChi2",&tPulse.FitChi2,"FitChi2/D");
//	treeFit1->Branch("FitNDF",&tPulse.FitNDF,"FitNDF/D");
//	treeFit1->Branch("FuncTime",&tPulse.FuncTime,"FuncTime/D");
//	treeFit1->Branch("FuncAmp",&tPulse.FuncAmp,"FuncAmp/D");
//
////	cout << endl << endl << (Int_t)treeFit1->GetEntries() << endl << endl;
////	treeFit1->Scan("*");
//
//
//
//	TF1* FExpGaus = (TF1*) gROOT->FindObjectAny("FExpGaus");
//	if(FExpGaus) FExpGaus->Delete();
//	FExpGaus = new TF1("FExpGaus",FuncExpGausMulti,-1000,9000,6); // range reset later
//	FExpGaus->SetNpx(10000);
//
//
//int counter=0;
//	int nr_entries=(Int_t)treeFit->GetEntries();
//	for(int i=0; i<nr_entries; i++) {
//		treeFit->GetEntry(i);
//
//		if(!(std::isfinite(dPulse.FitAmp) &&
//				std::isfinite(dPulse.FitTime) &&
//				std::isfinite(dPulse.FitBaseline) &&
//				std::isfinite(dPulse.FitRiseTime) &&
//				std::isfinite(dPulse.FitFallTime) &&
//				std::isfinite(dPulse.FitChi2) &&
//				std::isfinite(dPulse.FitNDF) &&
//				std::isfinite(dPulse.Baseline) &&
//				std::isfinite(dPulse.Time) &&
//				std::isfinite(dPulse.Amp) &&
//				std::isfinite(dPulse.TriggerTime) &&
//				std::isfinite(dPulse.FitLowLimit) &&
//				std::isfinite(dPulse.FitHighLimit))) {
//			dPulse.FitTime = -100.0;
//			dPulse.FitAmp = -100.0;
//			dPulse.FitBaseline = -100.0;
//			dPulse.FitRiseTime = -100.0;
//			dPulse.FitFallTime = -100.0;
//			dPulse.FitChi2 = -100.0;
//			dPulse.FitNDF = -100.0;
//			dPulse.Time = -100.0;
//			dPulse.Amp = -100.0;
//			dPulse.Baseline = -100.0;
//			dPulse.FitLowLimit = -100.0;
//			dPulse.FitHighLimit = -100.0;
//			dPulse.FuncTime = -100.0;
//			dPulse.FuncAmp = -100.0;
//			dPulse.TriggerTime = -100.0;
//			cout << "is not normal\t\t" << i << "\t\t" << dPulse.FitAmp << endl;
//			counter++;
//		}
//		else {
//			FExpGaus->SetParameter(0,0.0);
//			FExpGaus->SetParameter(1,dPulse.FitRiseTime);
//			FExpGaus->SetParameter(2,dPulse.FitFallTime);
//			FExpGaus->SetParameter(3,1);
//			FExpGaus->SetParameter(4,dPulse.FitAmp);
//			FExpGaus->SetParameter(5,dPulse.FitTime);
//			FExpGaus->SetRange(dPulse.FitTime-10,dPulse.FitTime+10);
//			dPulse.FuncTime = FExpGaus->GetMaximumX(dPulse.FitTime-10,dPulse.FitTime+10);
//			dPulse.FuncAmp = FExpGaus->Eval(dPulse.FuncTime);
//		}
//
//		tPulse=dPulse;
//		evt=evt_n;
//		pea_n=peaks_n;
//		pea_i=peak_i;
//		treeFit1->Fill();
//
//		loadbar(i, nr_entries, 1, 100, start2);
//	}
//	delete FExpGaus;
//
////	treeFit1->Scan("*");
//
////	loadbar(0, 1, 1, 100, start2);
//	cout << "\n\n\n\n" << counter*1.0/nr_entries << "\n";
//	treeFit1->Write();
//	delete treeFit1;
//	out2->Close();
//	delete out2;
}

void FOLDER::readRootAnalysis() {
//	double start2=getTime(0,0);
//	TFile *file = TFile::Open((filepath+"analysisWave.root").c_str(),"READ");
//	if(!file->IsOpen()) {
//		cout << "Error: cannot open " << filepath << "analysisWave.root !\n";
//		exit(-1);
//	}
//	cout << "Opening Analysis:\r";
//	treeAnalysis = (TTree*) file->Get("treeAnalysis");
//	treeAnalysis->SetBranchAddress("evt_n",&evt_n);
//	treeAnalysis->SetBranchAddress("trigger_time",&trigger_time);
//	treeAnalysis->SetBranchAddress("peaks_n",&peaks_n);
//	treeAnalysis->SetBranchAddress("peak_time",&peak_time);
//	treeAnalysis->SetBranchAddress("peak_value",&peak_value);
////	treeAnalysis->SetBranchAddress("base",&base);
////	treeAnalysis->SetBranchAddress("fluct",&fluct);
//	FOLDER::counterWaves = (Int_t)treeAnalysis->GetEntries();
////	treeAnalysis->Scan("*");
//
//	loadbar(0, 1, 1, 100, start2);
}

void FOLDER::analysisRootFits() {
//	double start2=getTime(0,0);
//	TFile *file = TFile::Open((filepath+"analysisFit.root").c_str(),"READ");
//	if(!file->IsOpen()) {
//		cout << "Error: cannot open " << filepath << "analysisFit.root !\n";
//		exit(-1);
//	}
//	treeFit = (TTree*) file->Get("treeFit");
//	treeFit->SetBranchAddress("evt_n",&evt_n);
//	treeFit->SetBranchAddress("peaks_n",&peaks_n);
//	treeFit->SetBranchAddress("pulses",&mPulse.Baseline);
//	loadbar(0, 1, 100, 100, start2);
//
//	TH1D *h1 = new TH1D("h1", "h1 title", 100, 0, 10);
//	TH1D *h2 = new TH1D("h2", "h2 title", 50, 0, 200);
//	int nr_entries=(Int_t)treeFit->GetEntries();
//	for(int i=0; i<nr_entries; i++) {
//		treeFit->GetEntry(i);
//		cout << mPulse.FitRiseTime << "\t" << mPulse.FitFallTime << endl;
//		h1->Fill(mPulse.FitRiseTime);
//		h2->Fill(mPulse.FitFallTime);
//	}
//
//	cout << endl << endl << "h1 mean:   " << h1->GetMean() << "\t h2 mean   " << h2->GetMean() << endl;
//	exit(-1);
//
//	delete h1;
//	delete h2;
}



