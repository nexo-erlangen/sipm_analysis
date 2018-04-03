/*
 * WAVE.h
 *
 *  Created on: Apr 15, 2016
 *      Author: sn0515
 */

#ifndef WAVE_H_
#define WAVE_H_

#include "MAININPUT.h"
#include "FOLDER.h"
#include "CFILE.h"
#include "DATA.h"
#include "structurs.h"
using namespace std;

class WAVE {
public:
	WAVE();
	WAVE(FOLDER*, unsigned int);
	virtual ~WAVE();
	void analysis();
	void writeWave(std::string, std::string, std::string);
	int getNrPeaks(string);
	tPulse* getPulses();
	static unsigned int counter;

private:
	std::string filepath;
	MAININPUT* info;
	unsigned int pulses_length;
	double* Pulse;
	double* Pulse_time;
	double trigger_time;
	double sampling;
	DATA* data;
	double baseline;
	double threshold;
	double limit[2];
	int NFitPulseMax=4;
	int NPulseMax=15;
	tPulse *mPulse = new tPulse [NPulseMax];
	int * PulseIndex = new int[NPulseMax];
	double * PulseTime = new double[NPulseMax];
	int nr_peaks=0;
	int nr_peaks_fit=0;

	void analysisSmooth(double*,unsigned int);
	void analysisSmoothBoxcarPlus(double*,unsigned int);
	void analysisSmoothBoxcarPlusMinus(double*,unsigned int);
	void analysisSmoothGauss(double*,unsigned int);
	void analysisBaseline();
	void analysisFitBaseline(int,int);
	void analysisIntegral();
	void analysisPeakHeight();
	void analysisPeakFinder();
	void analysisPeakFitter(int, double, double);
	void analysisPeakReFitter(double, double);
	void filterPulses();
	void setAnalysisLimits();
	void writePulses();
	void checkSanity(double) ;
	void calcPeak();
	std::string get_counter_str();
};

#endif /* WAVE_H_ */
