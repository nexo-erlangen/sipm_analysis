
/*
 * CFILE.cpp
 *
 *  Created on: Apr 18, 2016
 *      Author: sn0515
 */

#include "CFILE.h"

CFILE::CFILE() {
	// TODO Auto-generated constructor stub
	cout << "!!CFILE: Standardkonstruktor called!!" << endl;
}

CFILE::~CFILE() {
	// TODO Auto-generated destructor stub
	delete treeInput;
	out->Close();
	delete out;
}

CFILE::CFILE(FOLDER* object, unsigned int m) {
	folder=object;
	data=object->getData_ptr();
	treeInput=object->getTreeInput();
	info=object->getInfo_ptr();
	filepath=object->getFilepath();
	filename=object->getFile(m);
//	filenameTree=filepath+"Auswertung/temp/"+filename+".root";

	out = new TFile((filepath+filename+".root").c_str(),"RECREATE","",5);
	if(!out->IsOpen() || out->IsZombie()) {cout << "Error: cannot open " << filepath << "Auswertung/temp/input_" << m << ".root!" << endl;}
	treeInput = new TTree("treeInput","treeInput");
};

string CFILE::getFile() {
	return filename;
}

string CFILE::getRootFile() {
	if(out->IsOpen() && !out->IsZombie()) {return out->GetName();}
	else{return "None";}
}

void CFILE::readFile() {
	switch(info->type) {
		case 1 : readBinary(); break;
		case 2 : readOszi(); break;
		case 3 : readOszi_binary(); break;
	}
}

void CFILE::readOszi() {
	unsigned int pulses_file_num=0;
	unsigned int pulses_length=0;
	string unused;
	ifstream INPUT;
	INPUT.open((filepath+filename).c_str());
	if(!INPUT.is_open()) {
		cout << "Error: cannot open " << filepath << filename << " !" << endl;
	}
	getline(INPUT,unused);
	INPUT >> unused;
	INPUT >> pulses_file_num;
	INPUT >> unused;
	INPUT >> pulses_length;
	getline(INPUT,unused);
	getline(INPUT,unused);

	unsigned int evt_n;
	double trigger_time;
	double *x=new double [pulses_length];
	double *voltage = new double[pulses_length];
	double *trigger = new double [pulses_file_num];

	for(unsigned int n=0; n<pulses_file_num; n++)	{
		INPUT >> unused;
		INPUT >> unused;
		INPUT >> unused;
		INPUT >> trigger[n];
	}
	getline(INPUT,unused);
	getline(INPUT,unused);

	treeInput->Branch("evt_n",&evt_n,"evt_n/i");
	treeInput->Branch("length",&pulses_length,"length/i");
	treeInput->Branch("trigger",&trigger,"trigger/D");
	treeInput->Branch("time",x,"x[length]/D");
	treeInput->Branch("voltage",voltage,"voltage[length]/D");

	for(unsigned int n=0; n<pulses_file_num; n++) {
		for(unsigned int m=0; m<pulses_length; m++)	{						//loop over one pulses
			INPUT >> x[m];
			INPUT >> voltage[m];
			x[m]=1E9*x[m];
			voltage[m]=1E3*voltage[m]*info->sign;
		}
		evt_n=FOLDER::counterWaves;
		trigger_time=1E9*trigger[n];
		treeInput->Fill();
		++FOLDER::counterWaves;
	}
	INPUT.close();
	treeInput->Write();
//	tree->Scan("*");

	delete trigger;
	delete x;
	delete voltage;
}

void CFILE::readBinary() {
	//	THEADER th;
	//	EHEADER eh;
	//	char hdr[4];
	//	unsigned short voltage[1024];
	//	float bin_width[4][1024];
	//	int i, j, ch, n, chn_index;
	//	double t1, t2, dt;
	//	vector<pulses> PULSE;
	//	PULSE.resize(4);
	//	for(int i=0; i<4; i++) {
	//		PULSE.at(i).x.resize(1024);
	//		PULSE.at(i).value.resize(1024);
	//	}
	//	FILE *f = fopen((folder->filepath+filename).c_str(), "r");								// open the binary waveform file
	//	if (f == NULL) {
	//		printf("Cannot find file \'%s\'\n", (folder->filepath+filename).c_str());
	//		exit(-1);
	//	}
	//
	//	fread(&th, sizeof(th), 1, f);												// read time header
	//	printf("Found data for board #%d\n", th.board_serial_number);
	//
	//	memset(bin_width, sizeof(bin_width), 0);									// read time bin widths
	//	for (ch=0 ; ch<5 ; ch++) {
	//		fread(hdr, sizeof(hdr), 1, f);
	//		if (hdr[0] != 'C')
	//		{
	//			fseek(f, -4, SEEK_CUR);		// event header found
	//			break;
	//		}
	//		i = hdr[3] - '0' - 1;
	//		printf("Found timing calibration for channel #%d\n", i+1);
	//		fread(&bin_width[i][0], sizeof(float), 1024, f);
	//	}
	//
	//	for (n= 0 ; ; n++) {												// loop over all events in the data file
	//		i = (int)fread(&eh, sizeof(eh), 1, f);									// read event header
	//		if (i < 1) {
	//			break;
	//		}
	//
	//		for (ch=0 ; ch<5 ; ch++) {											// reach channel data
	//			i = (int)fread(hdr, sizeof(hdr), 1, f);
	//			if (i < 1) {
	//				break;
	//			}
	//			if (hdr[0] != 'C') {
	//				fseek(f, -4, SEEK_CUR);			// event header found
	//				break;
	//			}
	//			chn_index = hdr[3] - '0' - 1;
	//			fread(voltage, sizeof(short), 1024, f);
	//			for (i=0 ; i<1024 ; i++) {
	//				PULSE.at(chn_index).value.at(i) = input_info->sign*(voltage[i] / 65536. + eh.range/1000.0 - 0.5);		// convert data to volts
	//				for (j=0,PULSE.at(chn_index).x.at(i)=0 ; j<i ; j++)	{									// calculate time for this cell
	//					PULSE.at(chn_index).x.at(i) += bin_width[chn_index][(j+eh.trigger_cell) % 1024];
	//				}
	//				/*
	//				PULSE[chn_index][1][i] = MainInput->sign*(voltage[i] / 65536. + eh.range/1000.0 - 0.5);		// convert data to volts
	//				for (j=0,PULSE[chn_index][0][i]=0 ; j<i ; j++)											// calculate time for this cell
	//				{
	//					PULSE[chn_index][0][i] += bin_width[chn_index][(j+eh.trigger_cell) % 1024];
	//				}
	//				*/
	//			}
	//		}
	//		t1 =PULSE.at(1).x.at((1024-eh.trigger_cell) % 1024);						// align cell #0 of all channels
	//		for (ch=0 ; ch<3 ; ch++) {								 //!!change!!
	//			t2 = PULSE.at(ch).x.at((1024-eh.trigger_cell) % 1024);
	//			dt = t1 - t2;
	//			for (i=0 ; i<1024 ; i++) {
	//				PULSE.at(ch).x.at(i) += dt;
	//			}
	//		}
	//		for(unsigned int u=1; u<2 ; u++) {
	//			input_info->time_per_index=(PULSE.at(u).x.back()-PULSE.at(u).x.front())/(PULSE.at(u).x.size()-1);
	////			ANALYSIS_FRAME(folder, &PULSE.at(u), VALUES, MainInput);
	//		}
	//	}
}

void CFILE::readOszi_binary() {
	///////////////////////////////////////////////////////////////////////////
	// Deklarationen
	///////////////////////////////////////////////////////////////////////////

	// waveform display variables
	double pedmin=20;
	double pedmax=40;
	int maxtime0=200;

	FILE *pFile;
	int i=0, n=0, j=0, k=0, m=0, num=0,  number_points=0, seq_count=1, seq, seq_num;
	unsigned long number_in=0;
	int block0_length=0, data_address;
	float v_gain, v_offset, v_gain0, v_offset0;

	int ped0, pedn, ped02, pedn2;				// gate in channel

	double ped, ped2, signal; // charge, pC
	char c;
	string wdesc;
	float timestep, timestep0;
	double h_offset;
	short comm_type;

	float* timing;
	float* amplitude;

	int filename_length = (info->input_path+filename).length();

	int eventcount=-1;

	///////////////////////////////////////////////////////////////////////////
	// opening of files
	//////////////////////////////////////////////////////////////////////////

	pFile = fopen((info->input_path+filename).c_str(),"rb");

	if (pFile==NULL){
		cout << (info->input_path+filename).c_str() << " can not be opened!" << endl<<endl;
		exit(1);
	}
	/////////////////////////////
	//// reading of block0
	////////////////////////////
	fseek (pFile,0,SEEK_SET );
	fread (&c,1,1,pFile);
	if(c=='#') {
		fread (&c,1,1,pFile);
		if((c>='0')&&(c<='9'))
			block0_length=c-48+2;
		else cout<<"Warning: unknown format of block0"<<endl;
	}
	fseek (pFile,block0_length,SEEK_SET );
	fread (&wdesc[0],8,1,pFile);

	fseek (pFile,block0_length+32,SEEK_SET );
	fread (&comm_type, sizeof(comm_type),1,pFile);
	fread (&comm_type, sizeof(comm_type),1,pFile);

	/////////////////////////////////////////////
	//// calculation of data address
	/////////////////////////////////////////////
	data_address=block0_length;
	fread (&num, 4,1,pFile);
	data_address+=num;
	//cout << "length of WAVEDESC: " << num << endl;
	fread (&num, 4,1,pFile);
	data_address+=num;
	//cout << "length of USERTEXT: " << num << endl;
	fread (&num, 4,1,pFile);
	data_address+=num;
	//cout << "length of RES_DESC1: " << num << endl;
	fread (&num, 4,1,pFile);
	//cout << "data address of trigtime array: " << data_address << endl;
	// speichern von Anfangspunkt trigtime array:
	int startpoint_trigtime_array =  data_address;
	data_address+=num;
	//cout << "length of trigtime_array: " << num << endl;
	int length_trig_time_array = num;
	fread (&num, 4,1,pFile);
	data_address+=num;
	//cout << "length of ristime_array: " << num << endl;

	//fseek (pFile, block0_length+92 , SEEK_SET );
	//fread (&num, 4,1,pFile);
	//cout << "instrument number: " << num << endl;


	fseek (pFile,block0_length+32,SEEK_SET );
	fread (&comm_type, sizeof(comm_type),1,pFile);
	if(comm_type == 1)
	{
		cout << "Wrong format: 16-bit" << endl;
		fclose(pFile);
		exit(1);
	}

	//////////////////////////////////////////////////////////////////////////
	// Aus dem Waveform-Header entnimmt man die Werte fuer Vertical_Gain
	// und Vertical_Offset, die zur Berechnung der Amplituden aus den
	// raw data benoetigt werden.
	// Wieder muss der Offset von 11 Byte beruecksichtigt werden.
	//////////////////////////////////////////////////////////////////////////
	fseek (pFile, block0_length+156 , SEEK_SET );
	fread (&v_gain, sizeof(v_gain),1,pFile);
	//cout << "Vertical_gain  : " << v_gain;// << endl;

	fseek (pFile, block0_length+160 , SEEK_SET );
	fread (&v_offset, sizeof(v_offset),1,pFile);
	//cout << "    Vertical_offset: " << v_offset << endl;

	/////////////////////////////////////////////////////////////////////////////
	// number of sequences
	/////////////////////////////////////////////////////////////////////////////
	fseek (pFile, block0_length+174 , SEEK_SET );
	fread (&seq_count, 2,1,pFile);
	if(seq_count==0) {
		fseek (pFile, block0_length+112, SEEK_SET );
		fread (&seq_num, 2,1,pFile);
		fread (&seq_count, 2,1,pFile);
		seq_count*=seq_num;
	}

	/////////////////////////////////////////////////////////////////////////////
	// Total number of data points                                             //
	/////////////////////////////////////////////////////////////////////////////
	fseek (pFile, block0_length+116 , SEEK_SET );
	fread (&number_in, 4, 1, pFile);
	//cout << "Number of data points " << number_in << "!!!" << endl<< endl;

	/////////////////////////////////////////////////////////////////////////////
	// number of points in a single sequences                                  //
	/////////////////////////////////////////////////////////////////////////////
	number_points=number_in/seq_count;
	//cout<<"# of points in a single sequence="<<number_points<<endl;
	unsigned long num_total=number_points*(unsigned long)seq_count;
	//cout<<"# of points ="<<num_total<<endl;
	if (number_in!=num_total) {
		cout<<"Error: number of point in sequence is not digital or brocken sequence"<<endl;
		cout<<"Seq number="<<seq_count<<"  Total # of points="<<number_in<<"  # of points in sequence (truncated)="<<number_points<<endl;
		exit(1);
	}

	//////////////////////////////////////////////////////////////////////////////
	// Sampling interval for time between points                                //
	//////////////////////////////////////////////////////////////////////////////

	fseek (pFile, block0_length+176 , SEEK_SET);
	fread(&timestep, sizeof(timestep),1,pFile);
	//cout << "Time steps " << timestep << " s" << endl;

	fseek (pFile, block0_length+180, SEEK_SET);
	fread(&h_offset, sizeof(h_offset),1,pFile);
	//cout << "Horizontal_offset: " << h_offset << endl;

	//////////////////////////////////////////////////////////
	// Speichern der Zeit von ersten                        //
	// Trigger bis zum aktuellen                            //
	//////////////////////////////////////////////////////////

	fseek (pFile, startpoint_trigtime_array , SEEK_SET); // springt an Start von TRIGTIME ARRAY

	float rel_trigger[length_trig_time_array]={0};
	double trigger_time;
	treeInput->Branch("trigger",&trigger_time,"trigger/D");

	for( int i = 0; i < seq_count; i++){
		double memory_block_triggarray=0;
		fread(&memory_block_triggarray, 8,1,pFile);
		rel_trigger[i] = memory_block_triggarray * 1E9; // in ns umrechnen fuer weitere Verarbeitung
		//cout << "rel time to first trigger: " << rel_trigger[i] << endl;
		fread(&memory_block_triggarray, 8,1,pFile);
	}

	/////////////////////////////////////////////////////////////////////////////
	// Jetzt werden die Daten aus dem DATA_ARRAY_1 eingelesen. Da es sich nur um
	// die Werte fuer die Amplitude handelt muessen die zugehoerigen Zeiten ueber
	// die verwendete zeitliche Abtastrate berrechnet werden. Aus den raw data
	// berechnet man die Amplitude mittels vertical_gain und vertical_offset.
	/////////////////////////////////////////////////////////////////////////////
	fseek (pFile, data_address , SEEK_SET);

	double maxtime=(number_points+1)*timestep*1e9;
	//cout<<"maxtime0 "<<maxtime0<<endl;
	if(maxtime0>maxtime) maxtime0=maxtime;
	m=0;

	timing= new float[number_points];
	amplitude = new float[number_points];

	v_gain0=v_gain; v_offset0=v_offset;
	//cout<<v_gain<<endl;

	const int total_num_points=number_in;
	//cout <<"total number of points in file:\t" <<  total_num_points <<endl;

	// branches for trees:
	double *x=new double [number_points];
	double *voltage = new double[number_points];
	unsigned int evt_n;
	unsigned int pulses_length;

	//cout << pulses_length << "   " << number_points << endl ;

	treeInput->Branch("evt_n",&evt_n,"evt_n/i");
	treeInput->Branch("length",&pulses_length,"length/i");
	treeInput->Branch("time",x,"x[length]/D");
	treeInput->Branch("voltage",voltage,"voltage[length]/D");

	pulses_length = number_points;

	////////////////////////////////////////////
	//// Loop over sequences
	////////////////////////////////////////////
	for(seq=0; seq<seq_count; seq++){
		eventcount++;

		///////////////////////////////////////////
		//// Loop over waveform
		///////////////////////////////////////////
		double PED=0;
		for (n=0; n<number_points; n++)	{
			double amp,tim;
			timing[n]=0;
			amplitude[n]=0;
			char xx=0;
			fread(&xx, sizeof(xx),1,pFile);
			amplitude[n]= v_gain*xx - v_offset;

			amp=xx;
			tim=timing[n]= n*timestep*1e9;// + h_offset;

			x[n]=timing[n]+h_offset*1e9;
			voltage[n]=1E3*amplitude[n]*info->sign;
		}
		evt_n=FOLDER::counterWaves;
		trigger_time=rel_trigger[seq];
		treeInput->Fill();
		++FOLDER::counterWaves;
	}
	treeInput->Write();

	delete x;
	delete voltage;
	fclose(pFile);

	//~ treeInput->Scan("*");
	//cout<<"Number of evetnts "<<eventcount+1<<endl;
	//cout << "File is closed \n \n" << endl;
	//cout << "number of segments:\t" << seq_count << endl;
}
