//main.cpp

//author: Gen Onodera

#include "stdafx.h"

using namespace std;

int _tmain(int argc, _TCHAR* argv[])
{
	/* open data */
	cout << "Load started.\n";
	string filename = "./RF/52101_ssd.crf";
	ifstream fin(filename, ios_base::in | ios_base::binary);
	if (!fin){
		cout << "couldn't load file.\n";
		return 1;
	}
	cout << "Success!\n" << "Loaded " << filename << "\n";
	
	/* load header */
	fin.clear();
	fin.seekg(32, ios_base::beg);

	unsigned short len_record, frame, line, sample, ch;
	fin.read((char*)&len_record, sizeof(unsigned short));
	fin.read((char*)&frame, sizeof(unsigned short));
	const unsigned short frame2 = frame;
	fin.read((char*)&line, sizeof(unsigned short));
	fin.read((char*)&sample, sizeof(unsigned short));
	fin.read((char*)&ch, sizeof(unsigned short));
	cout << "-----RF Data Information-----\n" << "record length:" << len_record << "\n";
	cout << "number of frames:" << frame << "\n";
	cout << "number of lines:" << line << "\n";
	cout << "samples per line:" << sample << "\n";
	cout << "number of channel:" << ch << "\n";

	fin.seekg(24, ios_base::cur);
	char* probe_name = (char*)malloc(8+1);
	float frq_probe, max_angle, offset_from_center, rad_of_cuv,
		frq_t, frq_r, frq_s, acq_start, acq_end, range;
	unsigned short probe_type, pole, wave, burst, line_start, line_end, max_beam;
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
	cout << "probe name:" << probe_name << "\n";
	cout << "probe type:" << probe_type << "\n";
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

	fin.seekg(24, ios_base::cur);
	unsigned short focus_num, PRT;
	float focus_first, FR;
	fin.read((char*)&focus_num, sizeof(unsigned short));
	fin.read((char*)&focus_first, sizeof(float));
	fin.read((char*)&PRT, sizeof(unsigned short));
	fin.read((char*)&FR, sizeof(float));
	cout << "number of focusing:" << focus_num << "\n";
	cout << "first transmit focus[mm]:" << focus_first << "\n";
	cout << "PRT[us]:" << PRT << "\n";
	cout << "frame rate[Hz]:" << FR << "\n";

	fin.seekg(352, ios_base::beg);
	double RF_size;
	fin.read((char*)&RF_size, sizeof(double));
	cout << "RF data size:" << RF_size << endl;

	/* load RF data */
	// RF[frame][line][ch][sample]
	cout << "initializing array...\n";
	short tmp = 0;


	// random access method
	vector<vector<vector<vector<short>>>> RF(frame,
		vector<vector<vector<short>>>(line, vector<vector<short>>(ch, vector<short>(sample - 1, 0))));

	/*short ****RF = new short***[frame];
	for (int i = 0; i < frame; ++i){
		RF[i] = new short**[line];
		for (int j = 0; j < line; ++j){
			RF[i][j] = new short*[ch];
			for (int k = 0; k < ch; ++k){
				RF[i][j][k] = new short[sample - 1];
			}
		}
	}*/
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
	
	// push back method
	//vector<vector<vector<vector<short>>>> RF;
	/*for (int i = 0; i < frame; ++i){
		RF.push_back(vector<vector<vector<short>>>());
		for (int j = 0; j < line; ++j){
			RF[i].push_back(vector<vector<short>>());
			for (int k = 0; k < ch; ++k){
				fin.seekg(8, ios_base::cur);
				RF[i][j].push_back(vector<short>());
				for (int l = 0; l < sample - 1; ++l){
					fin.read((char*)&tmp, sizeof(short));
					if (tmp >= 2048) tmp -= 4096;
					RF[i][j][k].push_back(tmp);
				}
			}
		}
	}*/
	cout << "finished loading RF!\n";
	
	/* bias removal */
	int bias = 0;
	for (int i = 0; i < frame; ++i)
		for (int j = 0; j < line; ++j)
			for (int k = 0; k < ch; ++k){
				bias = accumulate(RF[i][j][k].begin(), RF[i][j][k].end(), 0);
				bias = bias / (sample - 1);
				for (int l = 0; l < sample - 1; ++l)
					RF[i][j][k][l] -= bias;	
			}
	cout << "finished removing bias!\n";

	///* channel RF draw */
	//string out = "element.dat";
	//ofstream fout(out, ios_base::out);
	//for (int i = 0; i < ch; ++i){
	//	for (int j = 0; j < sample - 1; ++j)
	//		fout << j << " " << RF[5][40][i][j] << "\n";
	//	fout << "\n";
	//}

	/* analize only 1 frame under here*/
	int nf = 5;

	/* interpolation */
	cout << "interpolating...\n";

	//spec and buffer setting
	int iporder = (int)(log((double)sample) / log(2.0));
	IppsFFTSpec_C_32fc *spec = 0;
	Ipp8u *specbuf, *initbuf, *workbuf;
	int size_spec, size_init, size_work;
	ippsFFTGetSize_C_32fc(iporder, IPP_FFT_NODIV_BY_ANY, ippAlgHintNone, &size_spec, &size_init, &size_work);
	//allocate buffer
	Ipp32fc *ipsrc = ippsMalloc_32fc((int)(4 * sample));
	for (int i = 0; i < sample - 1; ++i)
		ipsrc[i].re = RF[nf][40][40][i];
	Ipp32fc *ipdst = ippsMalloc_32fc((int)(4 * sample));
	specbuf = ippsMalloc_8u(size_spec);
	initbuf = ippsMalloc_8u(size_init);
	workbuf = ippsMalloc_8u(size_work);
	//initialize FFT
	ippsFFTInit_C_32fc(&spec, iporder, IPP_FFT_NODIV_BY_ANY, ippAlgHintNone, specbuf, initbuf);
	if (initbuf) ippFree(initbuf);
	//do FFT
	ippsFFTFwd_CToC_32fc(ipsrc, ipdst, spec, workbuf);
	//delete negative part
	for (int i = 0; i < sample / 2; ++i){
		ipdst[i + sample / 2].re = 0.0;
		ipdst[i + sample / 2].im = 0.0;
	}
	

	string out3 = "fft.dat";
	ofstream fout3(out3, ios_base::out);
	for (int j = 0; j < 4 * sample; ++j)
		fout3 << j << " " << ipdst[j].re << "\n";

	return 0;
}

