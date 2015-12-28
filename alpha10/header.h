#pragma once

#include "stdafx.h"

using namespace std;

/*functions*/
int physio();
void Bsector(const vector<vector<float>>& env, float dangle);
void cairo(const vector<vector<float>>& env, float dangle);
//void BSector2(const vector<vector<float>>& env, float dangle, float fs);
//void Bsector3(const vector<vector<float>>& env, float dangle);


/*class of fileopen.cpp*/
class file {
	ifstream fin;
public:	
	file(string filename);
	~file();

	void open(string filename);
	void start();
	void warp(int pos);
	void go(int pos);
};


class a10 : file{
	ifstream fin;
public:
	unsigned short len_record, frame, line, sample, ch;
	char* probe_name = (char*)malloc(8 + 1);
	float frq_probe, max_angle, offset_from_center, rad_of_cuv,
		frq_t, frq_r, frq_s, acq_start, acq_end, range;
	unsigned short probe_type, pole, wave, burst, line_start, line_end, max_beam;
	unsigned short focus_num, PRT;
	float focus_first, FR;
	double RF_size;
	vector<vector<vector<vector<short>>>> RF;
	vector<vector<vector<short>>> RF0;

	//function
	a10(string filename);
	~a10();
	
	void loadheader();
	void printheader();
	void loadRF();
	void loadRF0(int frame);
	void freeRF();
	void freeRF0();
	void rmbias();
	short eledat(int frame, int line, int ch, int sample);
};

/*class of physio.cpp*/
class physio : file{
	ifstream fin;
public:
	string fn;
	vector<short> ECG;
	vector<short> PCG_min;
	vector<short> PCG_max;

	physio(string filename);
	~physio();
	int extract(int offset);
};