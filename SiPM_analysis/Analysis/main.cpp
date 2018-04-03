/*
 * main.cpp
 *
 *      Author: sn0515
 *      compile:
 *	g++ -std=c++0x -O3 -o SiPM main.cpp MAININPUT.cpp FOLDER.cpp CFILE.cpp WAVE.cpp DATA.cpp functions.cpp structures.cpp
 */

#include "structurs.h"
#include "functions.h"
#include "MAININPUT.h"
#include "FOLDER.h"
#include "CFILE.h"
#include "WAVE.h"
#include "DATA.h"

int main(int argc, char** argv) {
	double start = getTime(0,0);
//	gErrorIgnoreLevel = kWarning;
	gErrorIgnoreLevel=kError;
	MAININPUT info(argc, argv);

	if(info.checkInput()) {
		vector<string> FOLDER_LIST;
		FOLDER_LIST=INPUT_DATA(info.getInputPath());

		for(unsigned int k=0; k<FOLDER_LIST.size(); k++) {							//loop over all folders including files with data
			FOLDER* folder= new FOLDER(&info, FOLDER_LIST.at(k));
			if(info.checkRead()) {
				folder->readConditions();
				folder->readFolder();
				folder->readFiles();
			}
/*			if(info.checkWave()) {
				folder->readRootInput();
//				folder->analysisRootFits();
				folder->analysisWaves();
				folder->writeData();
			}*/
			if(info.checkWave()) {
				folder->readRootInput();
				folder->analysisWaves();
				folder->writeOutput();
			}
			if(info.checkRun()) {
//				folder->readRootAnalysis();
				folder->readRootFit();
				folder->analysisNasty();
//				if(info.bool_calcFitPar) {folder->calcFitParameter();}
//				else {folder->analysisNasty();}
//				folder->analysisRun();
			}
			folder->printInfo();
			delete folder;
		}
		if(FOLDER_LIST.size()>1 && info.checkAll()) {
			cout << "total analysis" << endl;
		}
		cout << "\r" << endl;
		cout << "\t\tElapsed time: " << getTime(1,start) << " seconds \n";
	}
	else {
		info.printWarning();
	}
	cout << endl << ">>>>>>>>>>>>>>>>>>>>>>>> Programm finished <<<<<<<<<<<<<<<<<<<<<<<<" << endl << endl;
	return 0;
}
