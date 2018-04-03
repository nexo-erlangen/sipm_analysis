/*
 * CFILE.h
 *
 *  Created on: Apr 18, 2016
 *      Author: sn0515
 */

#ifndef CFILE_H_
#define CFILE_H_

class FOLDER;

#include "structurs.h"
#include "MAININPUT.h"
#include "FOLDER.h"
#include "WAVE.h"
#include "DATA.h"
using namespace std;

class CFILE {
public:
	CFILE();
	CFILE(FOLDER*, unsigned int);
	virtual ~CFILE();
	void readFile();
	std::string getFile();
	std::string getRootFile();

private:
	FOLDER* folder;
	DATA* data;
	MAININPUT* info;
	TFile *out;
	TTree *treeInput;
	std::string filepath;
	std::string filename;
	std::string filenameTree;
	void readBinary();
	void readOszi();
	void readOszi_binary();
};

#endif /* CFILE_H_ */
