/*
 * functions.cpp
 *
 *  Created on: Apr 21, 2016
 *      Author: sn0515
 */

#include "functions.h"

Double_t fpeaks(Double_t *x, Double_t *par) {
	Double_t result=0.0;
	result=par[4]*TMath::Exp(-1.0*(x[0]-par[0])/par[5]);
	result+=par[2]*TMath::Exp(-1.0*(x[0]-par[0])/par[3]);
	result*=1/(TMath::Exp(-1.0*(x[0]-par[0])/par[1])+1);
	return result;
}
Double_t base(Double_t *x, Double_t *par) {
	Double_t result;
	result=par[0];
	return result;
}
Double_t func_sum1(Double_t *x, Double_t *par) {
	return base(x,par) + fpeaks(x,&par[1]);
}
Double_t FuncExpGausMulti(Double_t *x, Double_t *p){
  // p[0]: baseline
  // p[1]: gaussian sig
  // p[2]: exponential decay constant
  // p[3]: number of pulses
  // p[4+2*i]: pulse[i] amplitude
  // p[5+2*i]: pulse[i] time
  // convolution of an exponential and a gaussian, formula found in a
  // chromatography paper on exponentially modified gaussians.
  // http://www.springerlink.com/content/yx7554182g164612/
	double val=p[0];
	for(int iPulse=0; iPulse<p[3]; iPulse++){
		double time = x[0]-p[5+2*iPulse];
		if(((p[1]*p[1]/p[2]/2.-time)/p[2])<700){ //otherwise exponential explodes
		  val+=p[4+2*iPulse]/2.*//p[2]/2.*
			  exp((p[1]*p[1]/p[2]/2.-time)/p[2])*
			  TMath::Erfc((p[1]*p[1]/p[2]-time)/sqrt(2)/p[1]);
		}
		//exp(1/2*p[1]*p[1]/p[2]/p[2]-time/p[2])*
		//(TMath::Erf(1/sqrt(2)*(p[5+2*iPulse]/p[1]+p[1]/p[2]))+
		// TMath::Erf(1/sqrt(2)*(time/p[1]-p[1]/p[2])));
	}
	return val;
}
Double_t FuncExpGausMulti2(double* x, double*p){
  // p[0]: baseline
  // p[1]: gaussian sig
  // p[2]: first  exponential decay constant
  // p[3]: second exponential fraction
  // p[4]: second exponential decay constant
  // p[5]: number of pulses
  // p[6+2*i]: pulse[i] amplitude
  // p[7+2*i]: pulse[i] time
  double val=0.;
  for(int iPulse=0; iPulse<p[5]; iPulse++){
    double time = x[0]-p[7+2*iPulse];
    if(((p[1]*p[1]/p[2]/2.-time)/p[2])<700){ //otherwise exponential explodes
      val+=(p[6+2*iPulse]/2.*exp((p[1]*p[1]/p[2]/2.-time)/p[2])*
	    TMath::Erfc((p[1]*p[1]/p[2]-time)/sqrt(2)/p[1]))*(1-p[3]);
      if(p[3]>0.){
	val+=(p[3]* p[6+2*iPulse]/2.*exp((p[1]*p[1]/p[4]/2.-time)/p[4])*
	      TMath::Erfc((p[1]*p[1]/p[4]-time)/sqrt(2)/p[1]));
      }
    }
  }
  val+=p[0];
  return val;
}

Double_t funcDNAPWGapRec(double* x, double* p){
	// p[0] = dark noise rate
	// p[1] = number of after-pulse time constant. Cannot float in fit
	// p[2] = fDeadTime
	// p[3] = fReadoutGapStart
	// p[4] = fReadoutGapEnd
	// p[5] = fMinTimeWhenRec
	// p[6] = fRecTime
	// p[2+i] = mean number of after-pulse for i time constant
	// p[3+i] = time constant
	// p[4+i] = fraction of pulses undergoing recovery
	double fDeadTime=p[2]; //15.0
	double fReadoutGapStart = p[3]; //2e4;
	double fReadoutGapEnd = p[4]; //4e4;
	double fMinTimeWhenRec = p[5]; //5;
	double fRecTime = p[6]; //49;
	if(x[0]<fDeadTime || (x[0]>=fReadoutGapStart && x[0]<=fReadoutGapEnd)) {return 0.;} // limit of pulse finder for all pulses
	double meanNPulseInTime = p[0];
	double meanNPulsePrev = p[0]*(x[0]-fDeadTime);
	//cout << meanNPulseInTime << " " << meanNPulsePrev << endl;
	if(x[0]>fReadoutGapEnd){
	meanNPulsePrev -= p[0]*(fReadoutGapEnd-fReadoutGapStart);
	}
	int nAP = (int) p[1];
	for(int iAP=0; iAP<nAP; iAP++){
	double tNPulseInTime =  p[7+3*iAP]*exp(-x[0]/p[8+3*iAP])/p[8+3*iAP];
	double tNPulsePrev = p[7+3*iAP]*(exp(-fDeadTime/p[8+3*iAP])-exp(-x[0]/p[8+3*iAP]));
	if(x[0]<fMinTimeWhenRec){
	  meanNPulseInTime += tNPulseInTime*(1-p[9+3*iAP]);
	  meanNPulsePrev += tNPulsePrev*(1-p[9+3*iAP]);
	}
	else{
	  double tauRecEff=p[7+3*iAP]*fRecTime/(p[7+3*iAP]+fRecTime);
	  meanNPulseInTime +=  tNPulseInTime*(1.-p[9+3*iAP]*exp(-x[0]/fRecTime));
	  meanNPulsePrev += tNPulsePrev-p[7+3*iAP]*p[9+3*iAP]*tauRecEff/p[7+3*iAP]*(exp(-fDeadTime/tauRecEff)-exp(-x[0]/tauRecEff));
	  if(x[0]>fReadoutGapEnd){ // assume 100% receovery by that time
	meanNPulsePrev -= p[7+3*iAP]*(exp(-fReadoutGapStart/p[8+3*iAP])-exp(-fReadoutGapEnd/p[8+3*iAP]));
	  }
	}
	}
	//cout << meanNPulseInTime << " " << meanNPulsePrev << endl;
	return exp(-meanNPulsePrev)*meanNPulseInTime;
}

Double_t funcDNR(double* x, double* p){
  // p[0] = dark noise rate
  // p[1] = prob for no AP
	p[0]/=1e9;
	return p[1]*p[0]*TMath::Exp(-1.0*x[0]*p[0]);
}

vector<string> INPUT_DATA(string filepath) {
	vector<string> FOLDER;
	size_t found1 = filepath.find(".dat",filepath.length()-4);
	size_t found2 = filepath.find(".txt",filepath.length()-4);
	if (found1!=string::npos || found2!=string::npos) {//input is file with list of folders
		ifstream INPUT;
		INPUT.open(filepath.c_str());
		while(!INPUT.eof()) {
			string folder="";
			INPUT >> folder;
			if(folder!="") {
				FOLDER.push_back(folder);
			}
		}
		INPUT.close();
	}
	else {//input is one folder
		size_t found = filepath.rfind("/");
		if(found!=filepath.length()-1) {
			filepath = filepath+"/";
		}

		stringstream TERMINAL;
		TERMINAL.str("");
		TERMINAL << "(cd " << filepath << "; ls -d $PWD/*_*/ 2> /dev/null > list.dat)" << endl;
		system(TERMINAL.str().c_str());
		TERMINAL.str("");

		ifstream INPUT;
		INPUT.open((filepath+"list.dat").c_str());
		while(!INPUT.eof())	{
			string folder="";
			INPUT >> folder;
			if(folder!="") {
				FOLDER.push_back(folder);
			}
		}
		if(FOLDER.size()==0) {
			TERMINAL << "(cd " << filepath << "; rm list.dat)" << endl;
			system(TERMINAL.str().c_str());
			TERMINAL.str("");
			FOLDER.push_back(filepath);
		}
		INPUT.close();
	}

	cout << endl << "=========Input=====================================================" << endl;
	cout << filepath << endl;
	cout << "Number of directories:\t" << FOLDER.size() << endl;
	cout << "===================================================================" << endl;

	return FOLDER;
}

//*******************************************************************************************************************************************
//*******************************************************************************************************************************************

void INPUT_HISTO(string file, vector<Data>* histo) {
	string dummy;
	ifstream INPUT_HISTO;

	INPUT_HISTO.open(file.c_str());
	while(!INPUT_HISTO.eof())
	{
		Data a;
		INPUT_HISTO >> a.x;
		INPUT_HISTO >> a.value;
		INPUT_HISTO >> dummy;
		histo->push_back(a);
	}
	INPUT_HISTO.close();
}

//*******************************************************************************************************************************************
//*******************************************************************************************************************************************

void OUTPUT_GNU_SCRIPT(string folder, Data* max, Data* min, const unsigned int counter) {
	ofstream OUTPUT;
	OUTPUT.open((folder+"Auswertung/temp/gnuplot_hist").c_str());
	OUTPUT 	<< "set xlabel \"maximum [V]\" \n"							//!!!write macro to use for multiple scipts!!!
			<< "set ylabel \"counts #\" \n"
			<< "unset title \n"
			<< "set xtics autofreq  norangelimit font \",15\" \n"
			<< "set ytics autofreq  norangelimit font \",15\" \n"
			<< "set cbtics autofreq  norangelimit font \",15\" \n"
			<< "set xlabel  offset character 0, 0, 0 font \"Arial,18\" textcolor lt -1 norotate \n"
			<< "set ylabel  offset character 0, 0, 0 font \"Arial,18\" textcolor lt -1 rotate by 90 \n"
			<< "set termoptions font \",11\" \n"
			<< "set fit errorvariables \n"
			<< "g(x,a,b,c)=a*exp(-0.5*(x-b)**2/c**2) \n"
			<< "multiG(x)=g(x,a1,b1,c1)";
	for(unsigned int i=2; i<=counter ;i++)
	{
		OUTPUT << "+g(x,a"<< i <<",b"<< i <<",c"<< i <<")";
	}
	OUTPUT 	<< "\n";
	for(unsigned int i=1; i<=counter ;i++)
	{
		OUTPUT << "b"<< i <<"=" << max[i].x <<" \n"
				<< "c"<< i <<"=" << (min[i].x-min[i-1].x)/3 <<" \n"
				<< "a"<< i <<"=" << max[i].value <<" \n";
	}
	OUTPUT 	<< "fit [0:"<< min[counter].x <<"] multiG(x) \"" << folder << "Auswertung/temp/Hist_spectrum.dat\" via a1,b1,c1";
	for(unsigned int i=2; i<=counter ;i++)
	{
		OUTPUT << ",a"<< i <<",b"<< i <<",c"<< i;
	}
	OUTPUT 	<< "\n";
	OUTPUT 	<< "plot \"" << folder << "Auswertung/temp/Hist_spectrum.dat\" w l lc 3, multiG(x) lc 1\n";
	OUTPUT.close();
}

//*******************************************************************************************************************************************
//*******************************************************************************************************************************************

double getTime(int t, int begin) {
    double time1=0;
    if (t==0) {
        time1=clock();
    }
    if (t==1) {
        double end=clock();
        time1=(end-begin)/CLOCKS_PER_SEC;
        time1=(int)(time1*100.0 + 0.5)/100.0;
    }
    return time1;
}

//********************************************************************
//-----Loadbar generieren---------------------------------------------
//********************************************************************

void loadbar(int x, int n, int r, int w, double start3) {
    x=x+1;
    static bool firstCall=true;
    if (firstCall)
    {
        firstCall=false;
        cout << "\r";
    }
    if ((x != n) && (x % (n / 100+1) != 0)) return;  // Only update r times
    float ratio = x/(float)n;
    int c=int(ratio * w);            // Calculuate the ratio of complete-to-incomplete.
    cout << "\t\t\t" << setw(4) << (int)(ratio*w) << "%|";
    /*
     for (int x=0;x<c;x++)                  // Show the load bar.
     {
     cout « "\u2588";
     }

     for (int x=c;x<w;x++)
     {
     cout « "\u2591";
     }
     */
    int mod=c%4;
    int counterL=c-mod;
    counterL=counterL/4;
    for(int i=0;i<counterL;i++)
    {
        cout << "\u25cf";
    }
    string circ="";
    if(mod==1)
    {
        cout << "\u25dc";
        circ = "\u25f7";
    }
    if(mod==2)
    {
        cout << "\u25dd";
        circ = "\u25f6";
    }
    if(mod==3)
    {
        cout << "\u25de";
        circ = "\u25f5";
    }
    if(mod==0)
    {
        cout << " ";
        circ = "\u25f4";
    }
    for(int i=counterL;i<w/4.0;i++)
    {
        cout << " ";
    }
    cout << "\b";
    cout << "|  " << getTime(1,start3) << " sec \r" << flush; // Move to the first column
}




