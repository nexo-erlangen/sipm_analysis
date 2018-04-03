/*
 * MAININPUT.cpp
 *
 *  Created on: Apr 18, 2016
 *      Author: sn0515
 */

#include "MAININPUT.h"

MAININPUT::MAININPUT() {
	// TODO Auto-generated constructor stub
	cout << "!!MAININPUT: Standardkonstruktor called!!" << endl;
}

MAININPUT::~MAININPUT() {
	// TODO Auto-generated destructor stub
}

MAININPUT::MAININPUT(int argc, char** argv) {
	for (int i=1; i<argc; i++)  {
		string tempstring=argv[i];
		if(tempstring == "-input") {
			input_path=argv[i+1];
		}
		if(tempstring == "-read") {
			bool_read=true;
		}
		if(tempstring == "-PSI")	{
			type=1;
		}
		if(tempstring == "-ascii") {
			type=2;
		}
		if(tempstring == "-binary") {
			type=3;
		}
		if(tempstring == "-sign") {
			sign=atof(argv[i+1]);
		}
		if(tempstring == "-wave") {
			bool_wave=true;
		}
		if(tempstring == "-single")	{
			single=atof(argv[i+1]);
		}
		if(tempstring == "-frames") {
			bool_frames=true;
			frames=atof(argv[i+1]);
		}
		if(tempstring == "-smooth") {//check if argv[i+2] does not start with "-", so if its not the next option
			smooth=argv[i+1];
			smooth_width=atof(argv[i+2]);
			if(smooth!="boxcar" && smooth!="boxcar_pm" && smooth!="gauss") {
				cout << endl << "\t\t!! invalide smooth option !! \n"
						<< "\t\t-- continue without smoothing --" << endl;
				smooth="";
			}
			else {
				if(smooth_width==0 || smooth_width>200) {
					smooth_width=25;
				}
			}
		}
		if(tempstring == "-integral") {
			bool_integral=true;
		}
		if(tempstring == "-datatree") {
			bool_datatree=true;
			type=4;
			//bool_integral=true;
			//bool_peak=true;
		}
		if(tempstring == "-find") {
			bool_find=true;
		}
		if(tempstring == "-fit") {
			bool_find=true;
			bool_fit=true;
			fit_type=atof(argv[i+1]);
			rise_time=atof(argv[i+2]);
			fall_time=atof(argv[i+3]);
			if(fit_type<1 || fit_type>3) {cout << "wrong fit type" << endl; exit(-1);}
			if(rise_time<=0.0 || fall_time <= 0.0) {cout << "wrong rise/fall times" << endl; exit(-1);}
		}
		if(tempstring == "-refit") {
			bool_refit=true;
		}
		if(tempstring == "-peak") {
			bool_peak=true;
		}
		if(tempstring == "-bounds") {
			bound_lower=atof(argv[i+1]);
			bound_upper=atof(argv[i+2]);
			if(bound_lower>bound_upper) {
				cout << "unsuitable parameters chosen for bounds" << endl;
				exit(-1);
			}
		}
		if(tempstring == "-baseline") {
			bool_baseline=true;
			type=5;
		}
		if(tempstring == "-run") {
			bool_run=true;
		}
		if(tempstring == "-calcFitPar") {
			bool_calcFitPar=true;
		}
		if(tempstring == "-prompt") {
			bool_prompt=true;
			nr_peaks=atof(argv[i+1]);
		}
		if(tempstring == "-crosstalk") {
			bool_crosstalk=true;
			pe1=atof(argv[i+1]);
			pe2=atof(argv[i+2]);
			gain=atof(argv[i+3]);
		}
		if(tempstring == "-dcr") {
			bool_dcr=true;
			gain=atof(argv[i+1]);
			dcr=atof(argv[i+2]);
			recTime=atof(argv[i+3]);
			dcrBegin=atof(argv[i+4]);
		}
		if(tempstring == "-2dhisto") {
			bool_2dHisto=true;
		}
		if(tempstring == "-2dWave") {
			bool_2dWave=true;
		}
		if(tempstring == "-rectime") {
			bool_rectime=true;
			gain=atof(argv[i+1]);
		}
		if(tempstring == "-verbose") {
			bool_run=true;
			bool_verbose=true;
		}
		if(tempstring == "-nodark") {
			bool_dark=false;
		}
		if(tempstring == "-all") {
			bool_all=true;
		}
		if(tempstring == "-root") {
			bool_root=true;
		}
		if(tempstring == "-awesome") {
			bool_awesome=true;
		}
		if(tempstring == "-help") {
			printHelp();
		}
	}
}

bool MAININPUT::checkInput() {
	if(input_path!="" && (checkRead() || checkWave() || checkRun() || checkAll())) {return true;}
	else {return false;}
}

bool MAININPUT::checkRead() {
	if(bool_read) {
		if((type==1 || type==2 || type==3) && fabs(sign)==1) {return true;}
		else {
			cout << "Warning: Input for data reading is invalid !" << endl;
			return false;
		}
	}
	else {return false;}
}

bool MAININPUT::checkWave() {
	if(bool_wave) {
		if (bool_peak || bool_integral || bool_fit || bool_find || bool_baseline || bool_datatree) {return true;}
		else {
			cout << "Warning: Input for wave analysis is invalid !" << endl;
			return false;
		}
	}
	else {return false;}
}

/*bool MAININPUT::checkDatatree() {
	if(bool_datatree) {return true;}
	else {
		cout << "Warning: Input for datatree analysis is invalid !" << endl;
		return false;
	}
}*/

bool MAININPUT::checkRun() {
	if(bool_run) {
		if(true) {return true;}
		else {
			cout << "Warning: .... !" << endl;
			return false;
		}
	}
	else {return false;}
}

bool MAININPUT::checkAll() {
	if(bool_all) {
		if(true) {return true;}
		else {
			cout << "Warning: Input for data reading is invalid !" << endl;
			return false;
		}
	}
	else {return false;}
}

void MAININPUT::printWarning() {
	cout << endl << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
	cout << "Not enough or suitable parameters!!" << endl;
	cout << "Program needs at least -input and (-read, -wave, -run or -all)" << endl;
	cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;

	printHelp();
}

void MAININPUT::printHelp() {
	cout 	<< endl << endl;
	cout 	<< "**************************************************************************" << endl;
	cout 	<< "\t\t\t" << "HELP FUNCTION"<< endl;
	cout 	<< "Find all possible options belonging to different analysis steps below." << endl;
	cout 	<< "**************************************************************************" << endl;
	cout 	<< "Option		||	Argument	||	Paramter(s) 	" << endl
			<< "--------------------------------------------------------------------------" << endl
			<< "Input path 	|| -input		|| 	< absolute path >		" << endl
			<< "Preprocess data	|| -read		|| 		" << endl
			<< "File type	||	-PSI		|| 		" << endl
			<< "File type	||	-ascii		|| 			" << endl
			<< "File type	||	-binary		|| 			" << endl
			<< "Sign of pulses	||	-sign		|| 	< -1/+1 >		" << endl
			<< "Wave Analysis	|| -wave		|| 		" << endl
			<< "Enable smoothing||	-smooth		|| 	< type > < binning(int) >	" << endl
			<< "Restrict # Waves||	-frames		||	< int>0 >	" << endl
			<< "Write # Waves	||	-single		|| 	< int>0 >		" << endl
			<< "Restrict Time	||	-bounds 	|| 	< double > < double >		" << endl
			<< "Analyze Baseline||	-baseline	||		" << endl
			<< "Peak height	||	-peak		|| 		" << endl
			<< "Integral	||	-integral	|| 		" << endl
			<< "Peakfinder	||	-find		||		" << endl
			<< "Finder + Fitter	||	-fit 		||	< type(int) > < rise(double) > < fall(double) >	"<< endl
			<< "Pulse Refit	||	-refit		|| 	 " << endl
			<< "Run Analysis	|| -run 		|| 		" << endl
			<< "PH Spectrum	||	-prompt 	|| 	< int >	" << endl
			<< "CT Analysis	||	-crosstalk	|| 	< 1pe(double) >	< 2pe(double) > < gain(double) >" << endl
			<< "DCR Analysis 	||	-dcr		|| 	< gain(double) > < gain(double) > < recTime(double) > < dcrBegin(double) >" << endl
			<< "2D Histo Plots	||	-2dHisto	||	 	" << endl
			<< "Persistence Plot||	-2dWave		||		" << endl
			<< "RecTime Analysis||	-rectime	||	< gain(double) >	" << endl
			<< "Final Analysis	|| -all			||		" << endl
			<< "**************************************************************************" << endl
			<< "Please note that all control parameters start with '-' sign!" << endl
			<< "Example: './Analysis -input /path-to-run/ -wave -single 10 -fit 2 2 80 -refit' " << endl
			<< "Good Luck :)" << endl
			<< endl;
	exit(1);
}

void MAININPUT::setInputPath(string path) {
	input_path=path;
	return;
}

string MAININPUT::getInputPath() {
	return input_path;
}



