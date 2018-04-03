/*
 * WAVE.cpp
 *
 *  Created on: Apr 15, 2016
 *      Author: sn0515
 */

#include "WAVE.h"

WAVE::WAVE() {
	// TODO Auto-generated constructor stub
	cout << "!!WAVE: Standardkonstruktor called!!" << endl;
	WAVE::counter++;
}

WAVE::~WAVE() {
	WAVE::counter++;
	delete mPulse;
	delete PulseIndex ;
	delete PulseTime;
	// TODO Auto-generated destructor stub
}

unsigned int WAVE::counter=0;

WAVE::WAVE(FOLDER* object, unsigned int i) {
	pulses_length=object->getEventLength(i);
	trigger_time=object->getEventTrigger(i);
	info=object->getInfo_ptr();
	data=object->getData_ptr();
	filepath=object->getFilepath();
	Pulse_time = object->getEventTime(i);
	Pulse = object->getEventValue(i);
	data->evt_n=WAVE::counter;
	sampling = round(((Pulse_time[pulses_length-1]-Pulse_time[0])/(pulses_length-1))*10)/10;
}

string WAVE::get_counter_str() {
	stringstream CONVERT;
	CONVERT.str("");
	CONVERT << WAVE::counter;
	string Peak=CONVERT.str();
	return Peak;
}

void WAVE::writeWave(string path="", string note1="", string note2="") {
	ofstream OUTPUT;
	OUTPUT.open((filepath+"Auswertung/frames/"+path+"/"+note1+get_counter_str()+info->hist_name+note2+".dat").c_str());
	for(unsigned int j=0; j<pulses_length; j++) {
		OUTPUT << Pulse_time[j] << "\t" << Pulse[j] << endl;
	}
	OUTPUT.close();
}

void WAVE::analysis() {
	if(info->smooth!="") {
		analysisSmooth(Pulse,info->smooth_width);
	}
	setAnalysisLimits();
	if(info->bool_baseline) {
		analysisBaseline();
	}
	if(info->bool_find) {
		analysisPeakFinder();
		if(info->bool_fit) {
			analysisFitBaseline(-750,-150);
			analysisPeakFitter(info->fit_type,info->rise_time,info->fall_time);
			if(info->bool_refit) {
				analysisPeakReFitter(info->rise_time,info->fall_time);
			}
			calcPeak();
			checkSanity(0.0);
			writePulses();
		}
	}
	if(info->bool_peak) {
		analysisPeakHeight();
	}
	if(info->bool_integral) {
		analysisIntegral();
	}
	if(info->bool_datatree) {
		analysisPeakHeight();
		analysisIntegral();
		analysisBaseline();
	}
}

void WAVE::setAnalysisLimits() {
	limit[0] = Pulse_time[0];
	if(info->bound_lower!=0.0 && info->bound_lower>limit[0] && info->bound_lower<Pulse_time[pulses_length-1]) {
		limit[0]=info->bound_lower;
	}
	limit[1] = Pulse_time[pulses_length-1];
	if(info->bound_upper!=0.0 && info->bound_upper<limit[1] && info->bound_upper>Pulse_time[0]) {
		limit[1]=info->bound_upper;
	}
}

void WAVE::analysisPeakHeight() {
	unsigned int max=0;
	unsigned int counter=0;
	double peak_value=-50000.0;
	double base=0.0;

	for(unsigned int i=0; i<pulses_length && Pulse_time[i]<=limit[1]; i++) {//find peak value
		if(Pulse_time[i]>limit[0]) {
			if(Pulse[i] > peak_value) {
				peak_value = Pulse[i];
				max=i;
			}
		}
	}
	for(unsigned int i=0; Pulse_time[i]<Pulse_time[max]-10; i++) {//calculate rauschen before peak
			base += Pulse[i];
			counter++;
	}

	if(counter!=0)	{base = base/counter;}
	else {base=Pulse[0];}

	data->maximum.push_back(fabs(peak_value-base));
	data->maximum_time.push_back(Pulse_time[max]);
	data->trigger_time.push_back(trigger_time);
}

void WAVE::analysisIntegral() {
	unsigned int countera=0;
	unsigned int counterb=0;
	double base=0.0;
	double integral_value=0.0;
//	cout << "limit[0]: " << limit[0] << endl;
//	cout << "limit[1]: " << limit[1] << endl;
	for(unsigned int i=0; Pulse_time[i]<=limit[1]; ++i) {
//		if(Pulse_time[i] < limit[0]) {
		if(Pulse_time[i] < -20.0) {
			base += Pulse[i];
			countera++;
		}
		if(Pulse_time[i] > limit[0]) {
//		else {
			integral_value += Pulse[i];
			counterb++;
		}
	}
	if(countera!=0 || base!=0)	{
		base = base/countera;
		integral_value=integral_value-(base*counterb);
	}
//	cout << "countera: " << countera << endl;
//	cout << "counterb: " << counterb << endl;
//	cout << "sampling: " << sampling << endl;
//	cout << "integral: " << integral_value*sampling << endl;
//	getchar();

	data->waveIntegral.push_back(integral_value*sampling);
}

int WAVE::getNrPeaks(string type) {
	if(type=="fit") {
		return nr_peaks_fit;
	}
	else if(type=="find") {
		return nr_peaks;
	}
	else if(type=="peak") {
		return 1;
	}
	else {return 0;}
}

tPulse* WAVE::getPulses() {
	return mPulse;
}

void WAVE::analysisPeakFinder() {//check https://root.cern.ch/doc/master/classTSpectrumFit.html
	double * Pulse_new = new double [pulses_length];
	unsigned int averWindow = 1; //not used
	unsigned int decon = 3;
	unsigned int sigma = 40;
	double threshold = 20.0;

	TSpectrum *spec = new TSpectrum(NPulseMax);
	Int_t nfound = spec->SearchHighRes(Pulse,Pulse_new,pulses_length,sigma,threshold,true,decon,false,averWindow);
	Double_t *xpeaks_i = spec->GetPositionX();
	std::sort(xpeaks_i,xpeaks_i+nfound);

//	cout << "nfound: " << nfound << endl;

	unsigned int counter_smooth=0;
	double base_smooth=0.0;
	for(unsigned int i=0; Pulse_time[i]<150; i++) {//calculate bases before trigger
		if(Pulse_time[i]>-750){
			base_smooth += Pulse_new[i];
			counter_smooth++;
		}
	}
	if(counter_smooth!=0)	{base_smooth = base_smooth/counter_smooth;}
	else {base_smooth=0;}

	int boundary = 20;
	int width=5;
	vector<int> xpeak_smooth;
	vector<int> xpeak_real;
	for(int i=1;i<=nfound;i++) {//find only proper peaks
		Int_t xpos=round(xpeaks_i[i-1]);
		if(Pulse_new[xpos]-base_smooth>2) {
			if(Pulse_time[xpos]>-50 && Pulse_time[xpos]<8500) {
				xpeak_smooth.push_back(xpos);
			}
		}
	}

	for(unsigned int i=0;i<xpeak_smooth.size() ; i++) {//find real peak in unsmoothed data
		int xpos_smooth=xpeak_smooth.at(i);
		vector<double> vec_diff_left;
		vector<double> vec_local_max;
		for(int j=xpos_smooth-100; j<=xpos_smooth; j++) {
			int m=1;
			bool local_max=false;
			while(Pulse[j] > Pulse[j-m] && Pulse[j] >= Pulse[j+m]) {
				if(m==width) {
					local_max=true;
					break;
				}
				m++;
			}
			if(local_max) {
				double diff_left=0.0;
				for(int k=1; k<boundary; k++) {
					diff_left += (Pulse[j-k+1] - Pulse[j-k-1]);
				}
				vec_diff_left.push_back(diff_left);
				vec_local_max.push_back(j);
			}
		}
		if(vec_diff_left.size()>0) {
			std::vector<double>::iterator diff_biggest = std::max_element(vec_diff_left.begin(), vec_diff_left.end());
			int xpos_temp = vec_local_max.at(std::distance(vec_diff_left.begin(), diff_biggest));
			xpeak_real.push_back(xpos_temp);
		}
		else{
			xpeak_real.push_back(xpeak_smooth.at(i));
			cout <<  WAVE::counter << "\t" << "gradient in unsmoothed data zero\t\t" << xpeak_smooth.at(i) << "  of  " << xpeak_smooth.size() << "\t\t at " << Pulse_time[xpeak_smooth.at(i)] << endl;
			ofstream OUTPUT;
			OUTPUT.open((filepath+"Auswertung/frames/"+get_counter_str()+"_test.dat").c_str());
			for(unsigned int j=0; j<pulses_length; j++) {
				OUTPUT << Pulse_time[j] << "\t" << Pulse[j] << "\t" << Pulse_new[j] << endl;
			}
			OUTPUT.close();
		}
	}

	nr_peaks=xpeak_real.size();
//	cout << "nr_peaks: " << nr_peaks << endl;
/*	for(int i=0;i<nr_peaks;i++) {//write peaks to array
		int xpos = xpeak_real.at(i);
		mPulse[i].Amp=Pulse[xpos];
		mPulse[i].Time=Pulse_time[xpos];
	}
*/
	unsigned int begin_base=30;
	unsigned int length_base=20;

	for(unsigned  int i=0;i<xpeak_real.size();i++) {//find real peak heights and related bases
		unsigned int counter_each=0;
		double base_each=0.0;
		int xpos = xpeak_real.at(i);
		for(unsigned int i=xpos-begin_base; i>=0 && i>=xpos-begin_base-length_base; i--) {//calculate rauschen before peak
			base_each += Pulse[i];
			counter_each++;
		}
		if(counter_each!=0)	{base_each = base_each/counter_each;}
		else {base_each=Pulse[0];}

		data->maximum_root.push_back(Pulse[xpos]-base_each);
		data->maximum_root_time.push_back(Pulse_time[xpos]);
		data->dummy.push_back(xpeak_real.size());
		data->dummy_time.push_back(i+1);
		data->trigger_time.push_back(trigger_time);
		mPulse[i].Amp=Pulse[xpos]-base_each;
		mPulse[i].Time=Pulse_time[xpos];
	}
	delete Pulse_new;
	delete spec;
}

void WAVE::analysisFitBaseline(int limit_lower, int limit_upper) {
	if (limit_lower > limit_upper) { cout << "baseline limits awkward" << endl; exit(-1);}
	TGraph *gr=new TGraph(pulses_length,Pulse_time,Pulse);
	TF1* fbase = (TF1*) gROOT->FindObjectAny("fbase");
	if(fbase) fbase->Delete();
	fbase = new TF1("fbase","pol0");
	gr->Fit("fbase","Q0","",limit_lower, limit_upper);
	baseline=fbase->GetParameter(0);
	delete gr;
	delete fbase;
}

void WAVE::analysisPeakFitter(int FitType=2, double RiseTime=2.0, double FallTime=80.0) {
	TGraph *gr=new TGraph(pulses_length,Pulse_time,Pulse);
	TF1* FExpGaus = (TF1*) gROOT->FindObjectAny("FExpGaus");
	if(FExpGaus) FExpGaus->Delete();
	FExpGaus = new TF1("FExpGaus",FuncExpGausMulti,Pulse_time[0],Pulse_time[pulses_length-1],4+NFitPulseMax*2); // range reset later
	FExpGaus->SetNpx(1000);
	TCanvas *c2 = new TCanvas("c2","c2",600,400);
	int i=0;
	nr_peaks_fit=nr_peaks;
	while(i<nr_peaks_fit) {
		mPulse[i].FitLowLimit = mPulse[i].Time-20; //-n*risetime -nbin for baseline
		if(mPulse[i].FitLowLimit<Pulse_time[0]) mPulse[i].FitLowLimit=Pulse_time[0];

		mPulse[i].FitHighLimit = mPulse[i].Time+200; //+tau1*3 +tau2
		if(mPulse[i].FitHighLimit>Pulse_time[pulses_length-1]) mPulse[i].FitHighLimit=Pulse_time[pulses_length-1];

		mPulse[i].Baseline=baseline;

		int tNPulseInGroup=1; // total number of pulses fitted together
		while(i+tNPulseInGroup<nr_peaks_fit &&
				  mPulse[i].FitHighLimit>mPulse[i+tNPulseInGroup].Time){
			  mPulse[i].FitHighLimit = mPulse[i+tNPulseInGroup].Time+FallTime*3+RiseTime*5;
			  tNPulseInGroup++;
		}
		if(mPulse[i].FitHighLimit>Pulse_time[pulses_length-1]) {mPulse[i].FitHighLimit = Pulse_time[pulses_length-1];}
		for(int iFitPulse=i; iFitPulse<i+tNPulseInGroup; iFitPulse++){
		  mPulse[iFitPulse].FirstPulseInGroup = i;
		  mPulse[iFitPulse].FitLowLimit=mPulse[i].FitLowLimit;
		  mPulse[iFitPulse].FitHighLimit=mPulse[i].FitHighLimit;
		  mPulse[iFitPulse].Baseline=mPulse[i].Baseline;
		}
		FExpGaus->FixParameter(3,tNPulseInGroup); //number of pulses

		switch(FitType){
			case 1: // fixed time constants but free baseline
				FExpGaus->ReleaseParameter(0);
//				FExpGaus->FixParameter(0,mPulse[i].Baseline);
				FExpGaus->SetParameter(0,mPulse[i].Baseline);
				FExpGaus->FixParameter(1,RiseTime);
				FExpGaus->FixParameter(2,FallTime);
//				FExpGaus->FixParameter(3,mPulse_Parameters.FitTime2Frac);
//				FExpGaus->FixParameter(4,mPulse_Parameters.FitFallTime2);
				break;
			case 2: // free time constant
				FExpGaus->FixParameter(0,mPulse[i].Baseline);
				FExpGaus->SetParameter(0,mPulse[i].Baseline);
				FExpGaus->ReleaseParameter(1);
				FExpGaus->SetParameter(1,RiseTime);
				FExpGaus->SetParLimits(1,0.1*RiseTime,10*RiseTime);
				FExpGaus->ReleaseParameter(2);
				FExpGaus->SetParameter(2,FallTime);
				FExpGaus->SetParLimits(2,0.1*FallTime,5*FallTime);
//				if(mPulse_Parameters.FitTime2Frac>0.){
//					mFExpGaus->ReleaseParameter(3);
//					mFExpGaus->SetParameter(3,mPulse_Parameters.FitTime2Frac);
//					mFExpGaus->SetParLimits(3,0.0,1.0);
//					mFExpGaus->ReleaseParameter(4);
//					mFExpGaus->SetParameter(4,mPulse_Parameters.FitFallTime2);
//					mFExpGaus->SetParLimits(4,0.1*mPulse_Parameters.FitFallTime2,10*mPulse_Parameters.FitFallTime2);
//				}
//				else{
//					mFExpGaus->FixParameter(3,0.0);
//					mFExpGaus->FixParameter(4,0.0);
//				}
				break;
			case 3: // free time constant but with tight limits
				FExpGaus->FixParameter(0,mPulse[i].Baseline);
				FExpGaus->ReleaseParameter(1);
				FExpGaus->SetParameter(1,RiseTime);
				FExpGaus->SetParLimits(1,RiseTime*0.1, RiseTime*10);
				FExpGaus->ReleaseParameter(2);
				FExpGaus->SetParameter(2,FallTime);
				FExpGaus->SetParLimits(2,FallTime*0.1, FallTime*10);
//				if(mPulse_Parameters.FitTime2Frac>0.){
//					FExpGaus->ReleaseParameter(3);
//					FExpGaus->SetParameter(3,mPulse_Parameters.FitTime2Frac);
//					FExpGaus->SetParLimits(3,0.0,1.0);
//					FExpGaus->ReleaseParameter(4);
//					FExpGaus->SetParameter(4,mPulse_Parameters.FitFallTime2);
//					FExpGaus->SetParLimits(4,mPulse_Parameters.FitFallTime2*0.1, mPulse_Parameters.FitFallTime2*10);
//				}
//				else{
//					FExpGaus->FixParameter(3,0.0);
//					FExpGaus->FixParameter(4,0.0);
//				}
				break;
		}
		for(int iFitPulse=0; iFitPulse<tNPulseInGroup; iFitPulse++){
			FExpGaus->ReleaseParameter(4+2*iFitPulse);
			FExpGaus->SetParameter(4+2*iFitPulse,(mPulse[iFitPulse+i].Amp));
//			FExpGaus->SetParameter(4+2*iFitPulse,(mPulse[iFitPulse+i].Amp-mPulse[iFitPulse+i].Baseline));
//			if(FitType==1) {
			FExpGaus->SetParLimits(4+2*iFitPulse,0.0,100.0);
//			cout << "amplitude: " << (mPulse[iFitPulse+i].Amp) << endl;
/*			if ( (mPulse[iFitPulse+i].Amp) < 0 ){
				ofstream OUTPUT;
				OUTPUT.open((filepath+"Auswertung/fit2/"+get_counter_str()+"_test.dat").c_str());
				for(unsigned int j=0; j<pulses_length; j++) {
					OUTPUT << Pulse_time[j] << "\t" << Pulse[j] << endl;
				}
				OUTPUT.close();
			}*/
//			}
			FExpGaus->ReleaseParameter(5+2*iFitPulse);
			FExpGaus->SetParameter(5+2*iFitPulse,mPulse[iFitPulse+i].Time);
			FExpGaus->SetParLimits(5+2*iFitPulse,mPulse[i].FitLowLimit,mPulse[i].FitHighLimit); // bad things happen when pulses are allowed in front
//			cout << "time: " << mPulse[iFitPulse+i].Time << endl;
//			cout << "time_low_limit: " << mPulse[i].FitLowLimit << endl;
//			cout << "time_high_limit: " << mPulse[i].FitHighLimit << endl;
			if(iFitPulse>0 && iFitPulse<tNPulseInGroup-1) {
				FExpGaus->SetParLimits(5+2*iFitPulse,mPulse[iFitPulse+i-1].Time,mPulse[iFitPulse+i+1].Time); // bad things happen when pulses are allowed in front.
			}
			if(iFitPulse==0 && tNPulseInGroup>1) {
				FExpGaus->SetParLimits(5,mPulse[i].FitLowLimit,mPulse[i+1].Time);
			}
			if(iFitPulse==tNPulseInGroup-1 && tNPulseInGroup>1) {
				FExpGaus->SetParLimits(5+2*(tNPulseInGroup-1),mPulse[tNPulseInGroup+i-2].Time,mPulse[tNPulseInGroup+i-1].FitHighLimit);
			}
		}
		for(int iFitPulse=tNPulseInGroup; iFitPulse<NFitPulseMax; iFitPulse++){
			FExpGaus->FixParameter(4+2*iFitPulse,0);
			FExpGaus->FixParameter(5+2*iFitPulse,0);
		}
//		for(int iFitPulse=0; iFitPulse<tNPulseInGroup; iFitPulse++){
//			double mint,maxt,minv,maxv;
//			FExpGaus->GetParLimits(5+2*iFitPulse,mint,maxt);
//			FExpGaus->GetParLimits(4+2*iFitPulse,minv,maxv);
//			cout << WAVE::counter << "\t" << iFitPulse << "\t" << FExpGaus->GetParameter(5+2*iFitPulse) << "\t" << mint << "\t" << maxt << "\t" << FExpGaus->GetParameter(4+2*iFitPulse) << "\t" << minv << "\t" << maxv << endl;
//			cout << WAVE::counter << "\t" << iFitPulse << "\t" << FExpGaus->GetParameter(5+2*iFitPulse) << "\t" << FExpGaus->GetParameter(4+2*iFitPulse) << endl;
//		}
		FExpGaus->SetRange(mPulse[i].FitLowLimit,mPulse[i].FitHighLimit);

//		for(int j=0; j<tNPulseInGroup; j++){
//			cout << WAVE::counter << "\t" << j << "\t" << tNPulseInGroup << "\t"
//					<< mPulse[j+i].Time<< "\t"
//					<< mPulse[j+i].Amp<< "\t"
//					<< FExpGaus->GetParameter(0) << "\t"
//					<< FExpGaus->GetParameter(1) << "\t"
//					<< FExpGaus->GetParameter(2) << "\t"
//					<< mPulse[j+i].FitLowLimit << "\t"
//					<< mPulse[j+i].FitHighLimit << endl;
//		}

		Int_t fitStatus = gr->Fit("FExpGaus","QRM"); //option: M - "more" - improve fit results // B - "bounds" - use parameter limits
//		if (fitStatus!=0 && fitStatus!=4000) {cout << WAVE::counter << "\t" << i << "\tFit 1 has error:\t" << fitStatus << endl;}
		for(int iFitPulse=0; iFitPulse<tNPulseInGroup; iFitPulse++){
			mPulse[iFitPulse+i].FitAmp = FExpGaus->GetParameter(4+2*iFitPulse);
			mPulse[iFitPulse+i].FitTime = FExpGaus->GetParameter(5+2*iFitPulse);
			mPulse[iFitPulse+i].FitBaseline = FExpGaus->GetParameter(0);
			mPulse[iFitPulse+i].FitRiseTime = FExpGaus->GetParameter(1);
			mPulse[iFitPulse+i].FitFallTime = FExpGaus->GetParameter(2);
			mPulse[iFitPulse+i].FitChi2 = FExpGaus->GetChisquare();
			mPulse[iFitPulse+i].FitNDF = FExpGaus->GetNDF();
			mPulse[iFitPulse+i].RefitChi2 = -1.0;
			mPulse[iFitPulse+i].RefitNDF = -1.0;

//			cout << WAVE::counter << "\t" << i << "\t" << iFitPulse+i << "\t" << tNPulseInGroup << "\t"
//					<< mPulse[iFitPulse+i].FitTime<< "\t"
//					<< mPulse[iFitPulse+i].FitAmp<< "\t"
//					<< FExpGaus->GetParameter(0) << "\t"
//					<< FExpGaus->GetParameter(1) << "\t"
//					<< FExpGaus->GetParameter(2) << "\t"
//					<< mPulse[iFitPulse+i].FitLowLimit << "\t"
//					<< mPulse[iFitPulse+i].FitHighLimit << endl;
//
			if(WAVE::counter<info->single) {
				stringstream CONVERT;
				CONVERT.str("");
				CONVERT << i+1;
				string i_str=CONVERT.str();
				gr->Draw();
				c2->SaveAs((filepath+"Auswertung/frames/4pe/"+get_counter_str()+"-fit-"+i_str+".pdf").c_str());
			}
		}
		i+=tNPulseInGroup;
	}
	// >>> Sort
	std::sort(mPulse, mPulse+nr_peaks_fit, [](tPulse const &a, tPulse const &b){ return a.FitTime < b.FitTime; });

	delete FExpGaus;
	delete gr;
	c2->Close();
	delete c2;
}

void WAVE::analysisPeakReFitter(double RiseTime=2.0, double FallTime=80.0){
	TCanvas *c2 = new TCanvas("c2","c2",600,400);
	TGraph *gr=new TGraph(pulses_length,Pulse_time,Pulse);
	TF1* FExpGaus = (TF1*) gROOT->FindObjectAny("FExpGaus");
	if(FExpGaus) FExpGaus->Delete();
	FExpGaus = new TF1("FExpGaus",FuncExpGausMulti,Pulse_time[0],Pulse_time[pulses_length-1],4+NFitPulseMax*2); // range reset later
	FExpGaus->SetNpx(1000);

	int nPulseBeforeRefit=nr_peaks_fit;
	int newPulseAdded=0;
	double MinChi2ForRefit=0.6;
	for(int iPulse=0; iPulse<nPulseBeforeRefit; iPulse++){
		if(mPulse[iPulse].FirstPulseInGroup==iPulse &&
		   mPulse[iPulse].FitNDF>0 &&
		   //(mPulse[iPulse].mFitChi2/mPulse[iPulse].mFitNDF/mPulse[iPulse].mAmp*mMinChi2ForRefit)>3.5){ // refit this group
		(mPulse[iPulse].FitChi2/mPulse[iPulse].FitNDF)>MinChi2ForRefit) {
	// >>> Set function parameters that will not change
			FExpGaus->SetRange(mPulse[iPulse].FitLowLimit,mPulse[iPulse].FitHighLimit);
			FExpGaus->ReleaseParameter(0);
			FExpGaus->FixParameter(0,baseline);
			FExpGaus->FixParameter(1,RiseTime);
			FExpGaus->FixParameter(2,FallTime);

	// >>> Copy chi2 to refit
			mPulse[iPulse].RefitChi2 = mPulse[iPulse].FitChi2;
			mPulse[iPulse].RefitNDF = mPulse[iPulse].FitNDF;

	// >>> add pulses as long as the chi2 is too large and that the added pulse is not too close to a another pulse
			double tMinTimeDiff=5.; // minimum time between two pulses
			int tNFitPulse=1;
			int validRefit=1;// stop refitting if not much progress is being made
			while(tNFitPulse<(NFitPulseMax-1) && //can add one more pulse due to n pulse fit limit
			nr_peaks_fit<(NPulseMax-1) &&  //can add one more pulse due to total n pulse limit
			mPulse[iPulse].RefitNDF>0 &&
			//(mPulse[iPulse].mRefitChi2/mPulse[iPulse].mRefitNDF/mPulse[iPulse].mAmp*mMinChi2ForRefit)>3.5 &&
			(mPulse[iPulse].RefitChi2/mPulse[iPulse].RefitNDF)>MinChi2ForRefit &&
			validRefit){
	// >>> Set fit function starting parameters without additional pulse
				// Needed to figure out where to add the next pulse
				tNFitPulse=0;
				for(int iFitPulse=iPulse; iFitPulse<nr_peaks_fit; iFitPulse++){
					if(mPulse[iFitPulse].FirstPulseInGroup==iPulse){
						FExpGaus->ReleaseParameter(4+2*tNFitPulse);
						FExpGaus->SetParameter(4+2*tNFitPulse,mPulse[iFitPulse].FitAmp);
						FExpGaus->SetParLimits(4+2*tNFitPulse,0,100);
						//FExpGaus->SetParLimits(4+2*iFitPulse,0,100);
						FExpGaus->ReleaseParameter(5+2*tNFitPulse);
						FExpGaus->SetParameter(5+2*tNFitPulse,mPulse[iFitPulse].FitTime);
						FExpGaus->SetParLimits(5+2*tNFitPulse,mPulse[iPulse].FitLowLimit,mPulse[iPulse].FitHighLimit);
						tNFitPulse++;
					}
				}
				FExpGaus->FixParameter(3,tNFitPulse);

	// >>> Look for most likely position of next pulse
				int tLastBin = ((mPulse[iPulse].FitHighLimit-Pulse_time[0])/sampling)-1;//mWF->GetXaxis()->FindBin(mPulse[iPulse].mFitHighLimit)-1;
				int tFirstBin = ((mPulse[iPulse].FitLowLimit-Pulse_time[0])/sampling)+1;
				int iMaxDiff;
				double maxDiff=0.;
				double newPulseTime=0.;
				double newPulseAmp=0;
				double tDiff;
				for(int iBin=tFirstBin;iBin<=tLastBin; iBin++){
					tDiff=(Pulse[iBin]-FExpGaus->Eval(Pulse_time[iBin]));
					if(maxDiff<tDiff){
						maxDiff=tDiff;
						newPulseTime =  Pulse_time[iBin];
						newPulseAmp = Pulse[iBin] - FExpGaus->Eval(Pulse_time[iBin]);
					}
				}

	// >>> Abort if too close to an existing pulse
				tMinTimeDiff=5.;
				for(int iFitPulse=0; iFitPulse<nr_peaks_fit; iFitPulse++){
					if(mPulse[iFitPulse].FirstPulseInGroup==iPulse){
						tDiff=fabs(mPulse[iFitPulse].FitTime-newPulseTime);
						if(tMinTimeDiff>tDiff) {
							tMinTimeDiff=tDiff;
						}
					}
				}

				validRefit=0;
				if(tMinTimeDiff>3.0){
	// >>> Add one more pulse
					FExpGaus->ReleaseParameter(4+2*tNFitPulse);
					FExpGaus->SetParameter(4+2*tNFitPulse,newPulseAmp);
					FExpGaus->SetParLimits(4+2*tNFitPulse,0.0,100.0);
					FExpGaus->ReleaseParameter(5+2*tNFitPulse);
					FExpGaus->SetParameter(5+2*tNFitPulse,newPulseTime);
					FExpGaus->SetParLimits(5+2*tNFitPulse,mPulse[iPulse].FitLowLimit,mPulse[iPulse].FitHighLimit);
					tNFitPulse++;
					FExpGaus->FixParameter(3,tNFitPulse);
					mPulse[nr_peaks_fit].FirstPulseInGroup=iPulse;
					nr_peaks_fit++;

	// >>> Fit
					Int_t fitStatus = gr->Fit("FExpGaus","QRM");
//					if (fitStatus!=0) {	cout << WAVE::counter << "\t" << iPulse << "\tMulti Fit has error:\t" << fitStatus << endl;}

					double improvement_least=0.02;
					if((mPulse[iPulse].RefitChi2-FExpGaus->GetChisquare())/mPulse[iPulse].RefitChi2>improvement_least){//chi2 has decreased by more than 2%
						newPulseAdded=1;
						validRefit=1;
	// >>> Copy fit information
						int iFitParameter=0;
						for(int iFitPulse=iPulse; iFitPulse<nr_peaks_fit; iFitPulse++){
							if(mPulse[iFitPulse].FirstPulseInGroup==iPulse){
								mPulse[iFitPulse].FitAmp = FExpGaus->GetParameter(4+2*iFitParameter);
								mPulse[iFitPulse].FitTime = FExpGaus->GetParameter(5+2*iFitParameter);
								mPulse[iFitPulse].FitBaseline = FExpGaus->GetParameter(0);
								mPulse[iFitPulse].FitRiseTime = FExpGaus->GetParameter(1);
								mPulse[iFitPulse].FitFallTime = FExpGaus->GetParameter(2);
								mPulse[iFitPulse].FitChi2 = mPulse[iPulse].FitChi2;
								mPulse[iFitPulse].FitNDF = mPulse[iPulse].FitNDF;
								mPulse[iFitPulse].RefitChi2 = FExpGaus->GetChisquare();
								mPulse[iFitPulse].RefitNDF = FExpGaus->GetNDF();
								iFitParameter++;
							}
						}
//						if(WAVE::counter<info->single) {
//							stringstream CONVERT;
//							CONVERT.str("");
//							CONVERT << iPulse+1;
//							string i_str=CONVERT.str();
//							gr->Draw();
//							c2->SaveAs((filepath+"Auswertung/frames/4pe/"+get_counter_str()+"-fit-"+i_str+"-new.pdf").c_str());
//						}
					}
					else{
						nr_peaks_fit--;
					}
				}
			}
		}
		iPulse++;
	}
	// >>> Sort
	if(newPulseAdded) {
		std::sort(mPulse, mPulse+nr_peaks_fit, [](tPulse const &a, tPulse const &b){ return a.FitTime < b.FitTime; });
	}

	delete FExpGaus;
	delete gr;
	c2->Close();
	delete c2;
}

void WAVE::filterPulses() {
	if(nr_peaks_fit>1){
		vector<int> wrong_pulses_i;
		std::sort(mPulse, mPulse+nr_peaks_fit, [](tPulse const &a, tPulse const &b){ return a.FitTime < b.FitTime; });

		for(int i=0; i<nr_peaks_fit-1; i++) {
			if(mPulse[i+1].FitTime-mPulse[i].FitTime<5) {
				double Time_temp=(mPulse[i].FitTime+mPulse[i+1].FitTime)/2.0;
				mPulse[i].Time=Time_temp;
				mPulse[i+1].Time=Time_temp;
				double Amp_temp=(mPulse[i].FitAmp+mPulse[i+1].FitAmp)+mPulse[i].FitBaseline;
				mPulse[i].Amp=Amp_temp;
				mPulse[i+1].Amp=Amp_temp;
				wrong_pulses_i.push_back(i);
				cout << WAVE::counter << "\t" << i << "\t" << nr_peaks_fit << "\t"
						<< mPulse[i].FitChi2/mPulse[i].FitNDF << "\t"
						<< mPulse[i].FitTime<< "\t"
						<< mPulse[i].FitAmp<< "\t"
						<< mPulse[i].Time<< "\t"
						<< mPulse[i].Amp<< "\t"
						<< mPulse[i].FitBaseline<< "\t"
						<< "wrong pulse" << "\t" << i << endl;
			}
		}
		for(int i=0; i<nr_peaks_fit; i++) {
			if((mPulse[i].Amp-mPulse[i].FitBaseline)<0.5) {
				wrong_pulses_i.push_back(i);

				cout << WAVE::counter << "\t" << i << "\t" << nr_peaks_fit << "\t"
					<< mPulse[i].FitChi2/mPulse[i].FitNDF << "\t"
					<< mPulse[i].FitTime<< "\t"
					<< mPulse[i].FitAmp<< "\t"
					<< "wrong pulse --- amplitude" << "\t" << i << endl;
			}
		}

		sort( wrong_pulses_i.begin(), wrong_pulses_i.end() );
		wrong_pulses_i.erase( unique( wrong_pulses_i.begin(), wrong_pulses_i.end() ), wrong_pulses_i.end() );

		for(unsigned int i=0; i<wrong_pulses_i.size(); i++) {
			for (int j = wrong_pulses_i.at(i); j < nr_peaks_fit; ++j) {
				mPulse[j] = mPulse[j+1];
			}
		}
		nr_peaks_fit -=wrong_pulses_i.size();
	}
}

void WAVE::calcPeak() {
	TF1* FExpGaus = (TF1*) gROOT->FindObjectAny("FExpGaus");
	if(FExpGaus) FExpGaus->Delete();
	FExpGaus = new TF1("FExpGaus",FuncExpGausMulti,Pulse_time[0],Pulse_time[pulses_length-1],6); // range reset later
	FExpGaus->SetNpx(10000);
	for(int i=0; i<nr_peaks_fit; i++) {
		FExpGaus->SetParameter(0,0.0);
		FExpGaus->SetParameter(1,mPulse[i].FitRiseTime);
		FExpGaus->SetParameter(2,mPulse[i].FitFallTime);
		FExpGaus->SetParameter(3,1);
		FExpGaus->SetParameter(4,mPulse[i].FitAmp);
		FExpGaus->SetParameter(5,mPulse[i].FitTime);
		FExpGaus->SetRange(mPulse[i].FitTime-10,mPulse[i].FitTime+10);
		mPulse[i].FuncTime = FExpGaus->GetMaximumX(mPulse[i].FitTime-10,mPulse[i].FitTime+10);
		mPulse[i].FuncAmp = FExpGaus->Eval(mPulse[i].FuncTime);
	}
	delete FExpGaus;
}

void WAVE::checkSanity(double dummy=0.0) {
	for(int i=0; i<nr_peaks_fit; i++) {
		if(!(std::isfinite(mPulse[i].FitAmp) &&
				std::isfinite(mPulse[i].FitTime) &&
				std::isfinite(mPulse[i].FitBaseline) &&
				std::isfinite(mPulse[i].FitRiseTime) &&
				std::isfinite(mPulse[i].FitFallTime) &&
				std::isfinite(mPulse[i].FitChi2) &&
				std::isfinite(mPulse[i].FitNDF) &&
				std::isfinite(mPulse[i].Baseline) &&
				std::isfinite(mPulse[i].Time) &&
				std::isfinite(mPulse[i].Amp) &&
				std::isfinite(mPulse[i].FitLowLimit) &&
				std::isfinite(mPulse[i].FitHighLimit) &&
				std::isfinite(mPulse[i].FuncTime) &&
				std::isfinite(mPulse[i].FuncAmp))) {
		mPulse[i].FitTime = dummy;
		mPulse[i].FitAmp = dummy;
		mPulse[i].FitBaseline = dummy;
		mPulse[i].FitRiseTime = dummy;
		mPulse[i].FitFallTime = dummy;
		mPulse[i].FitChi2 = dummy;
		mPulse[i].FitNDF = dummy;
		mPulse[i].Time = dummy;
		mPulse[i].Amp = dummy;
		mPulse[i].Baseline = dummy;
		mPulse[i].FitLowLimit = dummy;
		mPulse[i].FitHighLimit = dummy;
		mPulse[i].FuncTime = dummy;
		mPulse[i].FuncAmp = dummy;
		}
	}
}

void WAVE::writePulses() {
	for(int i=0; i<nr_peaks_fit; i++) {
		data->param.at(0).at(0).push_back(mPulse[i].FitChi2/mPulse[i].FitNDF);
		data->param.at(0).at(1).push_back(mPulse[i].FitBaseline);
		data->param.at(0).at(2).push_back(mPulse[i].FitRiseTime);
		data->param.at(0).at(3).push_back(mPulse[i].FitFallTime);
		data->param.at(0).at(4).push_back(mPulse[i].Time);
		data->param.at(0).at(5).push_back(mPulse[i].Amp);
		data->param.at(0).at(6).push_back(mPulse[i].FitTime);
		data->param.at(0).at(7).push_back(mPulse[i].FitAmp);
		data->param.at(0).at(8).push_back(mPulse[i].FuncTime);
		data->param.at(0).at(9).push_back(mPulse[i].FuncAmp);
		data->param.at(0).at(10).push_back(i);
	}
}

void WAVE::analysisBaseline() {
	if ( info->type == 5 ){
		TH1D *distr = new TH1D("distr","distr. histogram",10000,-100.0,100);
		for(unsigned int i=0; i<limit[0]; i++) {
			distr->Fill(Pulse[i]);
		}
		int bxmin=distr->FindFirstBinAbove(10.0);
		double xmin=distr->GetXaxis()->GetBinCenter(bxmin);
		int bxmax=distr->FindLastBinAbove(10.0);
		double xmax=distr->GetXaxis()->GetBinCenter(bxmax);
		distr->GetXaxis()->SetLimits(xmin,xmax);

		TF1 *func1 = new TF1("func1","gaus",xmin,xmax);
		func1->SetNpx(10000);
		func1->SetParameter(0,distr->GetMaximumBin());
		func1->SetParameter(1,distr->GetXaxis()->GetBinCenter(distr->GetMaximumBin()));
		func1->SetParameter(2,0.1);
		distr->Fit(func1,"Q");
		baseline=func1->GetParameter(1);
		if(func1->GetParameter(2)*3>1.5) {
			threshold=baseline+(func1->GetParameter(2)*3.0); //threshold = 3*sigma + baseline
		}
		else {
			threshold=baseline+1.0; //threshold = 1.5mV + baseline
		}

		data->baselines->Fill(baseline);
		data->thresholds->Fill(threshold);
		delete func1;
		delete distr;
	}
	if ( info->type == 4 ){
		// Baseline histogram
		TH1* h1 = new TH1I("h1", "h1 title", 100, -10.0, 10.0);
		for(unsigned int i=0; Pulse_time[i] < limit[0]; i++) {
			h1->Fill(Pulse[i]);
		}
		data->mean_base.push_back(h1->GetMean());
		data->rms_base.push_back(h1->GetRMS());
		delete h1;
	}
}

void WAVE::analysisSmooth(double* pulse, unsigned int width) {
	if(info->smooth=="boxcar") {
		analysisSmoothBoxcarPlus(pulse,width);
		info->hist_name="_bxpl";
	}
	else if(info->smooth=="boxcar_pm") {
		analysisSmoothBoxcarPlusMinus(pulse,width);
		info->hist_name="_bxplmn";
	}
	else if(info->smooth=="gauss") {
		analysisSmoothGauss(pulse,width);
		info->hist_name="_gauss";
	}
}

void WAVE::analysisSmoothBoxcarPlus(double* pulse, unsigned int width) {//fix smoothing of last points
	double *PULSE_BOXCAR = new double[pulses_length];
	for(unsigned int i=0; i<(pulses_length-width); i++) {
		PULSE_BOXCAR[i]=0;
		for(unsigned int j=0; j<width; j++) {
			PULSE_BOXCAR[i] += pulse[i+j];
		}
		PULSE_BOXCAR[i] /= width;
	}
	for(unsigned int i=0; i<(pulses_length-width); i++) {
		pulse[i]=PULSE_BOXCAR[i];
	}
	delete PULSE_BOXCAR;
}

void WAVE::analysisSmoothBoxcarPlusMinus(double* pulse, unsigned int width) {
	double *PULSE_BOXCAR = new double[pulses_length];
	for(unsigned int i=0; i<pulses_length; i++) {
		PULSE_BOXCAR[i]=0;
		if(i<width) {
				PULSE_BOXCAR[i] += pulse[i];
		}
//		else if(i==1) {
//			for(int j=0; j<BOXCAR_smoothing-1; j++)	{
//				if(j==0) {
//					PULSE_BOXCAR[i] = PULSE_BOXCAR[i] + PULSE_VALUE[i+j];
//				}
//				else {
//					PULSE_BOXCAR[i] = PULSE_BOXCAR[i] + PULSE_VALUE[i+j] + PULSE_VALUE[i-j];
//				}
//			}
//			PULSE_BOXCAR[i] = PULSE_BOXCAR[i] + PULSE_VALUE[i+2];
//		}
		else {
			for(unsigned int j=0; j<width; j++) {
				if(j==0) {
					PULSE_BOXCAR[i] += pulse[i+j];
				}
				else {
					PULSE_BOXCAR[i] += pulse[i-j] + pulse[i+j];
				}
			}
		}
		PULSE_BOXCAR[i] = PULSE_BOXCAR[i]/(width*2-1);
	}
	for(unsigned int i=0; i<pulses_length; i++) {
		pulse[i]=PULSE_BOXCAR[i];
    }
	delete PULSE_BOXCAR;
}

void WAVE::analysisSmoothGauss(double* pulse, unsigned int width) { //does not work properly!!
	double *PULSE_BOXCAR = new double[pulses_length];
	//double gauss_sigma[3]={0.68, 0.274, 0.043};
	double gauss_sigma[3]={0.5, 0.5, 0.5};
	double gauss_number= gauss_sigma[0]+gauss_sigma[1]*2+gauss_sigma[2]*2;
	for(unsigned int i=0; i<pulses_length; i++) {
		PULSE_BOXCAR[i]=0;
		if(i<width)	{
			PULSE_BOXCAR[i] = pulse[i];
		}
		else {
			for(unsigned int j=0; j<width; j++) {
				if(j==0) {
					PULSE_BOXCAR[i] += pulse[i+j]*gauss_sigma[j];
				}
				else {
					PULSE_BOXCAR[i] += pulse[i+j]*gauss_sigma[j] + pulse[i-j]*gauss_sigma[j];
				}
			}
		}
		PULSE_BOXCAR[i] /= gauss_number;
	}
    for(unsigned int i=0; i<pulses_length; i++) {
    	pulse[i]=PULSE_BOXCAR[i];
    }
    delete PULSE_BOXCAR;
}
