/*
 * DATA.h
 *
 *  Created on: Apr 21, 2016
 *      Author: sn0515
 */

#ifndef DATA_H_
#define DATA_H_

#include "structurs.h"
using namespace std;

class DATA {
public:
	DATA();
	virtual ~DATA();
//	tPulse fixParameters(int);
	void writeData(std::string);
	void writeDataTree(std::string);

	vector<double> waveIntegral;
	vector<double> maximum;
	vector<double> maximum_time;
	vector<double> trigger_time;
	vector<double> mean_base;
	vector<double> rms_base;
	vector<double> maximum_time_abs;
	vector<double> maximum_fit;
	vector<double> maximum_fit_time;
	vector<double> maximum_fit_time_abs;
	vector<double> maximum_root;
	vector<double> maximum_root_time;
	vector<double> dummy;
	vector<double> dummy_time;
	vector<double> dummy3;
	vector<vector<vector<Double_t>>> param;

	/*
	 * use optimal binning (maybe: scotts_bin_width)
	 * http://www.astroml.org/modules/generated/astroML.density_estimation.scotts_bin_width.html#astroML.density_estimation.scotts_bin_width
	 */
//	TH1D *maxima = new TH1D("maxima","maximum histogram",120,0.0,100);
	TH1D *baselines= new TH1D("baselines","baseline histogram",1000,-100.0,100);
	TH1D *thresholds = new TH1D("thresholds","threshold histogram",1000,-100.0,100);
	TH1D *dummyHist= new TH1D("dummy","dummy histogram",1000,-100.0,100);
	TH1D *dummyHist2 = new TH1D("dummy2","dummy histogram2",1000,-100.0,100);
//	TFile *out1;
//	TTree *treeFit;
//	TFile *out2;
//	TTree *treeAnalysis;

	unsigned int run;
	double sipmBias=0.0;
	double pmtBias=0.0;
	double temp=0.0;
	double press=0.0;
	double dist=0.0;
	double gain=0.0;
	double gain_err=0.0;
	double crosstalk=0.0;
	double crosstalk_err=0.0;

	double max;
	double max_time;
	double integral2;
	unsigned int evt_n;

	double Amp ;
	double Time;
	double Baseline ;
	double RiseTime ;
	double FallTime ;
	double Time2Frac ;
	double FallTime2 ;
	double Chi2;
	double NDF;
//	double RefitChi2;
//	double RefitNDF;

private:
};

#endif /* DATA_H_ */
