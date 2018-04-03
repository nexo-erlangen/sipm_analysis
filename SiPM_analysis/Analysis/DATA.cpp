/*
 * DATA.cpp
 *
 *  Created on: Apr 21, 2016
 *      Author: sn0515
 */

#include "DATA.h"

//DATA::DATA(string path) {
DATA::DATA() {
	// TODO Auto-generated constructor stub
//	filepath = path;
	param.resize(5);
	for(int i=0; i<5; i++) {
		param.at(i).resize(11);
	}

//	out1 = new TFile((filepath+"Auswertung/peak.root").c_str(),"RECREATE");
//	if(!out1->IsOpen() || out1->IsZombie()) {
//		cout << "Error: cannot open " << filepath << "Auswertung/peak.root !\n";
//		exit(-1);
//	}
//	else {
//		cout << "\nOpening " << filepath << "Auswertung/peak.root ! \n\n";
//	}
//	treeFit = new TTree("treeFit","treeFit");
//	treeFit->Branch("evt_n",&evt_n,"evt_n/i");
//	treeFit->Branch("Amp",&Amp,"Amp/D");
//	treeFit->Branch("Time",&Time,"Time/D");
//	treeFit->Branch("Baseline",&Baseline,"Baseline/D");
//	treeFit->Branch("RiseTime",&RiseTime,"RiseTime/D");
//	treeFit->Branch("FallTime",&FallTime,"FallTime/D");
//	treeFit->Branch("Time2Frac",&Time2Frac,"Time2Frac/D");
//	treeFit->Branch("FallTime2",&FallTime2,"FallTime2/D");
//	treeFit->Branch("Chi2",&Chi2,"Chi2/D");
//	treeFit->Branch("NDF",&NDF,"NDF/D");
//
//	out2 = new TFile((filepath+"Auswertung/analysisWave.root").c_str(),"RECREATE");
//	if(!out2->IsOpen() || out2->IsZombie()) {
//		cout << "Error: cannot open " << filepath << "Auswertung/analysisWave.root !\n";
//		exit(-1);
//	}
//	else {
//		cout << "\nOpening " << filepath << "Auswertung/analysisWave.root works! \n\n";
//	}
//	treeAnalysis = (TTree*) gROOT->FindObjectAny("treeAnalysis");
//	if(treeAnalysis) treeAnalysis->Delete();
//	treeAnalysis = new TTree("treeAnalysis","treeAnalysis");
//	treeAnalysis->Branch("evt_n",&evt_n,"evt_n/i");
//	treeAnalysis->Branch("max_time",&max_time,"max_time/D");
//	treeAnalysis->Branch("max",&max,"max/D");
//	treeAnalysis->Branch("integral",&integral2,"integral/D");

//	treeAnalysis->Branch("FitNPulse",&FitNPulse,"FitNPulse/i");
//	treeAnalysis->Branch("FitAmp",&FitAmp,"FitAmp/D");
//	treeAnalysis->Branch("FitTime",&FitTime,"FitTime/D");
//	treeAnalysis->Branch("FitBaseline",&FitBaseline,"FitBaseline/D");
//	treeAnalysis->Branch("FitRiseTime",&FitRiseTime,"FitRiseTime/D");
//	treeAnalysis->Branch("FitFallTime",&FitFallTime,"FitFallTime/D");
//	treeAnalysis->Branch("FitTime2Frac",&FitTime2Frac,"FitTime2Frac/D");
//	treeAnalysis->Branch("FitFallTime2",&FitFallTime2,"FitFallTime2/D");
//	treeAnalysis->Branch("FitChi2",&FitChi2,"FitChi2/D");
//	treeAnalysis->Branch("FitNDF",&FitNDF,"FitNDF/D");
}

DATA::~DATA() {
	// TODO Auto-generated destructor stub
//	delete maxima;
	delete dummyHist;
	delete dummyHist2;
}

//tPulse DATA::fixParameters(int FitType) {
//	tPulse mPulse;
//
//	return mPulse;
//}


void DATA::writeData(string filepath) {
	ofstream OUTPUT;
	if(dummy.size()>0) {
		OUTPUT.open((filepath+"Auswertung/peak_dummy.dat").c_str());
		for (unsigned int i = 0; i < dummy.size() ; ++i) {
			OUTPUT << dummy_time.at(i) << "\t" <<  dummy.at(i) << endl;
		}
		OUTPUT.close();
	}
	if(maximum_root.size()>0) {
		OUTPUT.open((filepath+"Auswertung/peak_root.dat").c_str());
		for (unsigned int i = 0; i < maximum_root_time.size() ; ++i) {
			double time_difference = 0.0;
			double time_difference_max = 0.0;
			if(i!=0 && trigger_time[i]!=0.0) {
				time_difference=trigger_time.at(i)-trigger_time.at(i-1);
				time_difference_max=trigger_time.at(i)+maximum_root_time.at(i)-trigger_time.at(i-1)-maximum_root_time.at(i-1);
			}
			OUTPUT << maximum_root_time.at(i) << "\t" <<  maximum_root.at(i) << "\t" << time_difference << "\t" << time_difference_max <<  endl;
		}
		OUTPUT.close();
	}
	if(maximum.size()>0) {
		OUTPUT.open((filepath+"Auswertung/peak.dat").c_str());
		for (unsigned int i = 0; i < maximum_time.size() ; ++i) {
			double time_difference = 0.0;
			double time_difference_max = 0.0;
			if(i!=0 && trigger_time[i]!=0.0) {
				time_difference=trigger_time.at(i)-trigger_time.at(i-1);
				time_difference_max=trigger_time.at(i)+maximum_time.at(i)-trigger_time.at(i-1)-maximum_time.at(i-1);
			}
			OUTPUT << maximum_time.at(i) << "\t" <<  maximum.at(i) << "\t" << time_difference << "\t" << time_difference_max << endl;
//			OUTPUT << maximum_time.at(i) << "\t" <<  maximum.at(i) << endl;
		}
		OUTPUT.close();
	}
	if(waveIntegral.size()>0) {
		OUTPUT.open((filepath+"Auswertung/peak_integral.dat").c_str());
		for (unsigned int i = 0; i < waveIntegral.size() ; ++i) {
				OUTPUT << waveIntegral.at(i) << endl;
		}
		OUTPUT.close();
	}

	if(maximum.size()>0 && waveIntegral.size()>0) {
		OUTPUT.open((filepath+"Auswertung/peak.dat").c_str());
		if(waveIntegral.size()>0) {
			for (unsigned int i = 0; i < maximum_time.size() ; ++i) {
				OUTPUT << maximum_time.at(i) << "\t" << maximum.at(i) << "\t" << waveIntegral.at(i)  << endl;
			}
		}
		OUTPUT.close();
	}
	if(maximum_fit.size()>0) {
		OUTPUT.open((filepath+"Auswertung/peak_fit.dat").c_str());
		for (unsigned int i = 0; i < maximum_fit_time.size() ; ++i) {
			OUTPUT << maximum_fit_time.at(i) << "\t" <<  maximum_fit.at(i) << endl;
		}
		OUTPUT.close();
	}
	for(int i=0; i<5 ;i++) {
		if(param.at(i).at(1).size()>0) {
			stringstream CONVERT;
			CONVERT.str("");
			CONVERT << (i+1);
			string Peak=CONVERT.str();
			OUTPUT.open((filepath+"Auswertung/peak_param_fit_"+Peak+".dat").c_str());
			for (unsigned int j = 0; j < param.at(i).at(1).size() ; ++j) {
				for (int k=0 ; k<11; k++) {
					OUTPUT << param.at(i).at(k).at(j) << "\t";
				}
				OUTPUT << endl;
			}
			OUTPUT.close();
		}
	}

//	TFile *out1 = TFile::Open((filepath+"Auswertung/peak.root").c_str(),"RECREATE");
//	if(!out1->IsOpen()) {
//		cout << "Error: cannot open " << filepath << "Auswertung/peak.root !" << endl;
//	}
//	treeFit->Scan("*");
//	treeFit->Write();
//	out1->Close();

//	TFile *out2 = TFile::Open((filepath+"Auswertung/analysisWave.root").c_str(),"RECREATE");
//	if(!out2->IsOpen()) {
//		cout << "Error: cannot open " << filepath << "Auswertung/analysisWave.root !" << endl;
//	}
//	treeAnalysis->Scan("*");
//	treeAnalysis->Write();
//	out2->Close();
	//	delete treeFit;
	//	delete treeAnalysis;
//		delete out1;
//		delete out2;

}


void DATA::writeDataTree(string filepath) {
	TFile *file = new TFile((filepath+"datatree.root").c_str(),"RECREATE");
	if(!file->IsOpen() || file->IsZombie()) {cout << "Error: cannot open " << filepath << "datatree.root !" << endl; exit(-1);}
	unsigned int evt_n;
	double peak_integral;
	double peak_maximum;
	double time_maximum;
	//vector<double>  peakfinder_maximum;
	//vector<double>  peakfinder_max_time;
	double baseline_mean;
	double baseline_rms;

	TTree * treeData = new TTree("treeData","treeData");
	//treeData -> SetAutoSave();
	//gDirectory -> Purge();
	treeData->Branch("evt_n", &evt_n, "evt_n/i");
	treeData->Branch("peak_integral",&peak_integral,"peak_integral/D");
	treeData->Branch("peak_maximum",&peak_maximum,"peak_maximum/D");
	treeData->Branch("time_maximum",&time_maximum,"time_maximum/D");
	//treeData->Branch("peakfinder_maximum", &peakfinder_maximum, "peakfinder_maximum[15]/D");
	//treeData->Branch("peakfinder_max_time", &peakfinder_max_time, "peakfinder_max_time[15]/D");
	treeData->Branch("baseline_mean",&baseline_mean,"baseline_mean/D");
	treeData->Branch("baseline_rms",&baseline_rms,"baseline_rms/D");
	// Filling Tree
	for (unsigned int i = 0; i < waveIntegral.size() ; ++i) {
		evt_n = i;
		peak_integral = waveIntegral.at(i);
		peak_maximum = maximum.at(i);
		time_maximum = maximum_time.at(i);
		//peakfinder_maximum[i] = maximum_root[i];
		//peakfinder_max_time[i] = maximum_root_time[i];
		baseline_mean = mean_base.at(i);
		baseline_rms = rms_base.at(i);
		treeData->Fill();
	}
	//treeData->Write();
	treeData->Write("", TObject::kOverwrite); // ueberschreibt ersten tree und speichert nur zweiten
	delete treeData;
	file->Close();
	delete file;
}
