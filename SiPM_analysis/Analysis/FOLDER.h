/*
 * FOLDER.h
 *
 *  Created on: Apr 15, 2016
 *      Author: sn0515
 */

#ifndef FOLDER_H_
#define FOLDER_H_

class CFILE;

#include "structurs.h"
#include "functions.h"
#include "MAININPUT.h"
#include "CFILE.h"
#include "DATA.h"
using namespace std;

class FOLDER {
public:
	FOLDER();
	FOLDER(MAININPUT*, const std::string);
	virtual ~FOLDER();
	friend class CFILE;					//???was tut das????
	void readConditions();
	void readFolder();
	void readFiles();
	void readRootInput();
	void readRootAnalysis();
	void readRootFit();
	void readRootFit_ChangeTree();
	void calcFitParameter();

	void analysisWaves();
	void analysisRootFits();
	void analysisRun();
	void analysisNasty();

	vector<std::string> getFolder() const;
	std::string getFilepath() const;
	std::string getFile(unsigned int) const;
	unsigned int getNumberFiles() const;
	MAININPUT* getInfo_ptr() const;
	DATA* getData_ptr();
	TTree* getTreeInput();
	double* getEventTime(unsigned int i);
	double* getEventValue(unsigned int i);
	unsigned int getEventLength(unsigned int i) const;
	double getEventTrigger(unsigned int i) const;
	void writeData();
	void writeDataTree();
	void writeOutput();
	void printInfo();

private:
	void setStructure() const;
	void makingHistogram(const std::string, string, string, double, int);
	void make2DHistogramWave();
	void make2DHistogram(string , string, double , double , int , double , double , int );
	void findHistogramPeaks(const string , string, unsigned const int,vector<Data>* , Data* , Data*);
	void writeGnuScriptSpec(const string , string, string, Data* , Data* , const unsigned int );
	void writeGnuScriptGain(const string , string , string , int) ;
	void writeGnuScriptDCR(const string , string , string , double);
	void writeGnuScriptRecTime(const string , string , string , double, int);
	void writeGnuScript2DHisto(const string , string , string , double , double , int , double , double , int ) ;
	void fitGnuMultigauss(const string , Data* , Data* , unsigned const int );
	void readGnuHistogram(string ,string, vector<Data>* );
	void writePulseData(string , double , double , unsigned int, unsigned int, unsigned int, double);
	void writePulseDataRec(string , double, double);
	void writeTimeDiff(string , string, double , double , double, double, double);
	void writeTimeDiff2(string , double , double , double);
	void calcCrosstalkLinear(string, string,double, double);
	void calcCrosstalkQuadratic(string, string, double);
	std::string filepath;
	vector<std::string> files;
	MAININPUT* info;
	TTree *treeInput;
	TTree *treeAnalysis;
	TTree *treeFit;
	tPulse *mPulse;
	tPulse dPulse;
//	TTree *treeInputInfo;
	DATA data;
	static const unsigned int maxWaves=100000;
	unsigned int evt_n;
	unsigned int length;
	double trigger;
	unsigned int peaks_n;
	unsigned int peak_i;
	double base;
	double trigger_time;
	double fluct;
	double peak_time;
	double peak_value;
	double *x = new double [maxWaves];
	double *voltage = new double [maxWaves];
	static unsigned int counterWaves;
};

#endif /* FOLDER_H_ */
