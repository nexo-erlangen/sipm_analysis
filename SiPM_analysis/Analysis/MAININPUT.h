/*
 * MAININPUT.h
 *
 *  Created on: Apr 18, 2016
 *      Author: sn0515
 */

#ifndef MAININPUT_H_
#define MAININPUT_H_

#include "structurs.h"
using namespace std;

class MAININPUT {
public:
	MAININPUT();
	MAININPUT(int, char**);
	virtual ~MAININPUT();
	bool checkInput();
	bool checkRead();
	bool checkWave();
	bool checkDatatree();
	bool checkRun();
	bool checkAll();
	void printWarning();
	void printHelp();
	void setInputPath(std::string);
	std::string getInputPath();

	std::string input_path="";
	std::string smooth="";
	std::string hist_name="";
	vector <unsigned int> channel;
	unsigned int channel_pmt=0;
	int type=0;
	unsigned int single=0;
	unsigned int frames=0;
	int sign=0;
	double bound_lower=0.0;
	double bound_upper=0.0;
	unsigned int smooth_width=0;
	int nr_peaks=0;
	int fit_type=0;
	double rise_time=0.0;
	double fall_time=0.0;
	double pe1=0.0;
	double pe2=0.0;
	double gain=0.0;
	double dcrBegin=0.0;
	double recTime=0.0;
	double dcr=0.0;
	bool bool_frames=false;
	bool bool_integral=false;
	bool bool_peak=false;
	bool bool_datatree=false;
	bool bool_find=false;
	bool bool_run=false;
	bool bool_calcFitPar=false;
	bool bool_prompt=false;
	bool bool_crosstalk=false;
	bool bool_dcr=false;
	bool bool_rectime=false;
	bool bool_2dHisto=false;
	bool bool_2dWave=false;
	bool bool_verbose=false;
	bool bool_dark=true;
	bool bool_all=false;
	bool bool_root=false;
	bool bool_fit=false;
	bool bool_refit=false;
	bool bool_awesome=false;
	bool bool_read=false;
	bool bool_wave=false;
	bool bool_baseline=false;
//	bool bool_sipm=false;
//	bool bool_pmt=false;
};

#endif /* MAININPUT_H_ */
