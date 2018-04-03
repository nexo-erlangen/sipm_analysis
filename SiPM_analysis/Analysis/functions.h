/*
 * functions.h
 *
 *  Created on: Apr 21, 2016
 *      Author: sn0515
 */

#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

#include "structurs.h"
#include <iomanip>

using namespace std;

Double_t fpeaks(Double_t *, Double_t*);
Double_t base(Double_t *, Double_t*);
Double_t func_sum1(Double_t *, Double_t*);
Double_t FuncExpGausMulti(Double_t *, Double_t *);
Double_t FuncExpGausMulti2(double*, double*);
Double_t funcDNAPWGapRec(double* , double* );
Double_t funcDNR(double* , double*) ;
vector<string> INPUT_DATA(string);
void INPUT_HISTO(string, vector<Data>*);

void loadbar(int x, int n, int r, int w, double start3);
double getTime(int t, int begin);
void OUTPUT_GNU_SCRIPT(string, Data*, Data*, const unsigned int);

#endif /* FUNCTIONS_H_ */
