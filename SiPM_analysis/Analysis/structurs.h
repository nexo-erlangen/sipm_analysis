/*
 * structurs.h
 *
 *  Created on: 15.11.2015
 *      Author: patrick
 */

#ifndef STRUCTURS_H_
#define STRUCTURS_H_

#include <iostream>
#include <fstream>
#include <cmath>
#include <fcntl.h>
#include <limits>
#include <stdio.h>
//#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <sstream>
#include <time.h>
#include <complex>
#include <iomanip>
#include <stdlib.h>
#include <sys/time.h>
#include <vector>
#include <algorithm>

#include "TROOT.h"
#include "TTree.h"
#include "TChain.h"
//#include "TTreeReader.h"
//#include "TTreeReaderValue.h"
//#include "TTreeReaderArray.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
//#include "TStyle.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TRandom.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"
#include "TMath.h"
#include "TPolyMarker.h"

using namespace std;
using std::string;

#define _USE_MATH_DEFINES
#define CLOCKS_PER_MS (CLOCKS_PER_SEC / 1000)
#define FWHM 2.354820045030949382023138652919399275494771378771641077045

struct tPulse{
  double Baseline;
//  double AbsAmp;
  double Amp;
  double Time;
//  double Q;
//  double Width;
//  double SPTemplateChi2;
  double FitLowLimit;
  double FitHighLimit;
  double FitBaseline;
  double FitTime;
  double FitAmp;
  double FitRiseTime;
  double FitFallTime;
  double FitTime2Frac;
  double FitFallTime2;
  double FitChi2;
  double FitNDF;
  double RefitChi2;
  double RefitNDF;
  double FuncTime;
  double FuncAmp;
  double TriggerTime;
  int FirstPulseInGroup;
};

struct Data {
	double value;
	double x;
};

typedef struct {
   char           time_header[4];
   char           bn[2];
   unsigned short board_serial_number;
} THEADER;

typedef struct {
   char           event_header[4];
   unsigned int   event_serial_number;
   unsigned short year;
   unsigned short month;
   unsigned short day;
   unsigned short hour;
   unsigned short minute;
   unsigned short second;
   unsigned short millisecond;
   unsigned short range;
   char           bs[2];
   unsigned short board_serial_number;
   char           tc[2];
   unsigned short trigger_cell;
} EHEADER;

#endif /* STRUCTURS_H_ */
