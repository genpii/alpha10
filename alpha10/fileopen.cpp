// fileopen.cpp

#include "stdafx.h"
#include "fileopen.h"

using namespace std;

a10::a10(string filename)
{
	cout << "open file=" << filename << "\n";
	fin.open(filename, ios_base::in | ios_base::binary);
	if (!fin){
		cout << "couldn't load file.\n";
	}
}

a10::~a10()
{
	cout << "delete object.\n";
}

void a10::start()
{
	fin.clear();
	fin.seekg(0, ios_base::beg);
}

void a10::warp(int pos)
{
	fin.seekg(pos, ios_base::beg);
}

void a10::go(int pos)
{
	fin.seekg(pos, ios_base::cur);
}

void a10::loadheader()
{
	a10::start();
	fin.seekg(32, ios_base::beg);

	fin.read((char*)&len_record, sizeof(unsigned short));
	fin.read((char*)&frame, sizeof(unsigned short));
	frame = frame - 1;
	fin.read((char*)&line, sizeof(unsigned short));
	fin.read((char*)&sample, sizeof(unsigned short));
	fin.read((char*)&ch, sizeof(unsigned short));

	fin.seekg(24, ios_base::cur);

	fin.read(probe_name, 8);
	probe_name[8] = '\0';
	fin.read((char*)&probe_type, sizeof(unsigned short));
	fin.read((char*)&frq_probe, sizeof(float));
	fin.read((char*)&pole, sizeof(unsigned short));
	fin.read((char*)&wave, sizeof(unsigned short));
	fin.read((char*)&max_angle, sizeof(float));
	fin.read((char*)&offset_from_center, sizeof(float));
	fin.read((char*)&rad_of_cuv, sizeof(float));
	fin.read((char*)&frq_t, sizeof(float));
	fin.read((char*)&frq_r, sizeof(float));
	fin.read((char*)&frq_s, sizeof(float));
	fin.read((char*)&burst, sizeof(unsigned short));
	fin.read((char*)&acq_start, sizeof(float));
	fin.read((char*)&acq_end, sizeof(float));
	fin.read((char*)&line_start, sizeof(unsigned short));
	fin.read((char*)&line_end, sizeof(unsigned short));
	fin.read((char*)&max_beam, sizeof(unsigned short));
	fin.read((char*)&range, sizeof(float));

	fin.seekg(24, ios_base::cur);

	fin.read((char*)&focus_num, sizeof(unsigned short));
	fin.read((char*)&focus_first, sizeof(float));
	fin.read((char*)&PRT, sizeof(unsigned short));
	fin.read((char*)&FR, sizeof(float));

	fin.seekg(352, ios_base::beg);

	fin.read((char*)&RF_size, sizeof(double));
}

void a10::printheader()
{
	cout << "-----RF Data Information-----\n" << "record length:" << len_record << "\n";
	cout << "number of frames:" << frame << "\n";
	cout << "number of lines:" << line << "\n";
	cout << "samples per line:" << sample << "\n";
	cout << "number of channel:" << ch << "\n";
	cout << "probe name:UST-" << probe_name << "\n";
	switch (probe_type)
	{
	case 1:
		cout << "probe type:linear\n";
		break;
	case 2:
		cout << "probe type:convex\n";
		break;
	case 3:
		cout << "probe type:sector\n";
		break;
	case 4:
		cout << "probe type:annular\n";
		break;
	default:
		cout << "probe type:unknown\n";
		break;
	}
	cout << "probe frequency[MHz]:" << frq_probe << "\n";
	cout << "transmit pole:" << pole << "\n";
	cout << "wave pattern:" << wave << "\n";
	cout << "max angle of probe:" << max_angle << "\n";
	cout << "offset from center[mm]:" << offset_from_center << "\n";
	cout << "radius of curvature[mm]:" << rad_of_cuv << "\n";
	cout << "transmit frequency[MHz]:" << frq_t << "\n";
	cout << "receiving frequency[MHz]:" << frq_r << "\n";
	cout << "sampling frequency[MHz]:" << frq_s << "\n";
	cout << "burst cycle:" << burst << "\n";
	cout << "top of ROI[mm]:" << acq_start << "\n";
	cout << "bottom of ROI[mm]:" << acq_end << "\n";
	cout << "beam number of acquire start:" << line_start << "\n";
	cout << "beam number of acquire end:" << line_end << "\n";
	cout << "max number of beams per frame:" << max_beam << "\n";
	cout << "display range[mm]:" << range << "\n";
	cout << "number of focusing:" << focus_num << "\n";
	cout << "first transmit focus[mm]:" << focus_first << "\n";
	cout << "PRT[us]:" << PRT << "\n";
	cout << "frame rate[Hz]:" << FR << "\n";
	cout << "RF data size:" << RF_size << endl;
}

void a10::loadRF()
{
	fin.seekg(360, ios_base::beg);
	fin.seekg((line * ch * (sample + 3))* sizeof(short), ios_base::cur);

	RF = vector<vector<vector<vector<short>>>>(frame, vector<vector<vector<short>>>(line,
		vector<vector<short>>(ch, vector<short>(sample - 1, 0))));

	short tmp;
	cout << "loading RF...\n";
	for (int i = 0; i < frame; ++i)
		for (int j = 0; j < line; ++j){
			for (int k = 0; k < ch - 16; ++k){ // back of 80 elements
				fin.seekg(8, ios_base::cur); //attribute 6byte channel number 2byte
				for (int l = 0; l < sample - 1; ++l){
					fin.read((char*)&tmp, sizeof(short));
					RF[i][j][k + 16][l] = tmp - 2048;
				}
			}
			for (int k = 0; k < 16; ++k){ // front of 16 elements
				fin.seekg(8, ios_base::cur);
				for (int l = 0; l < sample - 1; ++l){
					fin.read((char*)&tmp, sizeof(short));
					RF[i][j][k][l] = tmp - 2048;
				}
			}
		}
}

void a10::loadRF0(int frame)
{
	fin.seekg(360, ios_base::beg);
	fin.seekg((line * ch * (sample + 3))* sizeof(short), ios_base::cur); //ignore first frame

	if (RF0.empty())
		RF0 = vector<vector<vector<short>>>(line, vector<vector<short>>(ch, vector<short>(sample - 1, 0)));

	for (int i = 0; i < frame - 1;++i)
		fin.seekg((line * ch * (sample + 3))* sizeof(short), ios_base::cur);

	short tmp;
	for (int i = 0; i < line; ++i){
		for (int j = 0; j < ch - 16; ++j){ // back of 80 elements
			fin.seekg(8, ios_base::cur); //attribute 6byte channel number 2byte
			for (int k = 0; k < sample - 1; ++k){
				fin.read((char*)&tmp, sizeof(short));
				RF0[i][j + 16][k] = tmp - 2048;
			}
		}
		for (int j = 0; j < 16; ++j){ // front of 16 elements
			fin.seekg(8, ios_base::cur);
			for (int k = 0; k < sample - 1; ++k){
				fin.read((char*)&tmp, sizeof(short));
				RF0[i][j][k] = tmp - 2048;
			}
		}
	}

}

void a10::freeRF()
{
	if (!RF.empty())
		vector<vector<vector<vector<short>>>>().swap(RF);
}

void a10::freeRF0()
{
	if (!RF0.empty())
		vector<vector<vector<short>>>().swap(RF0);
}

void a10::rmbias()
{
	int bias = 0;
	for (int i = 0; i < frame; ++i)
		for (int j = 0; j < line; ++j)
			for (int k = 0; k < ch; ++k){
				bias = accumulate(RF[i][j][k].begin(), RF[i][j][k].end(), 0);
				bias = bias / (sample - 1);
				for (int l = 0; l < sample - 1; ++l)
					RF[i][j][k][l] -= bias;
			}
}