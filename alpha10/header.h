#pragma once

#include "stdafx.h"

using namespace std;

/*functions*/
int physio();
void Bsector(const vector<vector<float>>& env, float dangle);
void cairo(const vector<vector<float>>& env, float dangle);
void BSector2(const vector<vector<float>>& env, float dangle, float fs, int frm, float max);
float BSector22(const vector<vector<float>>& env, float dangle, float fs, int frm);
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


class a10{
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
	void open(string filename);
	void start();
	void warp(int pos);
	void go(int pos);
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

/*functions of DSP.cpp*/
void writev(vector<float> &v, string &str);
void writev(vector<float> &v, float &sc, string &str);
vector<float> NormCC(const vector<float> &x1, const vector<float> &x2, int lag, int offset);
vector<float> pwrspe(const vector<float> &v, const int &order);

/*class of GN.cpp*/
class GN{
	vector<vector<float>> jacob;
	vector<vector<float>> jacobt;
	vector<vector<float>> jj;
	vector<vector<float>> jjinv;
	vector<float> beta;
	vector<double> res;
	int m, n;
	Ipp32f *x, *y;

	float tmp;

public:
	GN(const vector<double> &x_ini, const vector<double> &y_ini);
	//~GN();
	void setj();
	void trans();
	void loadpoint();
	void mul(const vector<vector<float>> &a, const vector<vector<float>> &b);
	void mul(const vector<vector<float>> &a, const vector<float> &b);
	vector<double> solve(const double &r_ini, const double &c_ini);
};