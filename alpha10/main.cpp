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

	//delete first frame
	fin.seekg((((sample - 1) + 8) * ch * line) * 2, ios_base::cur);
	frame = frame - 1;

	/* load RF data */
	// RF[frame][line][ch][sample]
	cout << "initializing array...\n";
	short tmp = 0;


	frame = 1;
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

	cout << "removing bias...\n";
	
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

	///* channel RF draw */
	//string out = "element.dat";
	//ofstream fout(out, ios_base::out);
	//for (int i = 0; i < ch; ++i){
	//	for (int j = 0; j < sample - 1; ++j)
	//		fout << j << " " << RF[5][40][i][j] << "\n";
	//	fout << "\n";
	//}

	/* do FFT and IFFT */
	cout << "creating analytic signal...\n";
	
	//initialize focused RF array
	vector<vector<vector<vector<float>>>> elere(frame,
		vector<vector<vector<float>>>(line, vector<vector<float>>(ch, vector<float>(4 * sample, 0))));
	vector<vector<vector<vector<float>>>> eleim(frame,
		vector<vector<vector<float>>>(line, vector<vector<float>>(ch, vector<float>(4 * sample, 0))));

	//spec and buffer setting for FFT
	Ipp8u *specbuff, *initbuff, *workbuff;
	Ipp8u *specbufi, *initbufi, *workbufi;
	int size_specf, size_initf, size_workf;
	int size_speci, size_initi, size_worki;
	IppsFFTSpec_C_32fc *specf = 0;
	IppsFFTSpec_C_32fc *speci = 0;
	Ipp32fc *ipsrc = ippsMalloc_32fc((int)(4 * sample));
	Ipp32fc *ipdst = ippsMalloc_32fc((int)(4 * sample));
	Ipp32fc *ipdst2 = ippsMalloc_32fc((int)(4 * sample));
	const int fftorder = (int)(log((double)sample) / log(2.0));
	const int ifftorder = (int)(log((double)(4 * sample)) / log(2.0));
	ippsFFTGetSize_C_32fc(fftorder, IPP_FFT_NODIV_BY_ANY, ippAlgHintNone, &size_specf, &size_initf, &size_workf);
	ippsFFTGetSize_C_32fc(ifftorder, IPP_FFT_NODIV_BY_ANY, ippAlgHintNone, &size_speci, &size_initi, &size_worki);
	specbuff = ippsMalloc_8u(size_specf);
	specbufi = ippsMalloc_8u(size_speci);
	initbuff = ippsMalloc_8u(size_initf);
	initbufi = ippsMalloc_8u(size_initi);
	workbuff = ippsMalloc_8u(size_workf);
	workbufi = ippsMalloc_8u(size_worki);
	ippsFFTInit_C_32fc(&specf, fftorder, IPP_FFT_NODIV_BY_ANY, ippAlgHintNone, specbuff, initbuff);
	ippsFFTInit_C_32fc(&speci, ifftorder, IPP_FFT_NODIV_BY_ANY, ippAlgHintNone, specbufi, initbufi);

	for (int i = 0; i < frame; ++i){
		for (int j = 0; j < line; ++j){
			for (int k = 0; k < ch; ++k){
				//set
				for (int l = 0; l < sample - 1; ++l)
					ipsrc[l].re = RF[i][j][k][l];

				//do FFT
				ippsFFTFwd_CToC_32fc(ipsrc, ipdst, specf, workbuff);
				
				//double positive part and delete negative part
				for (int l = 0; l < sample / 2; ++l){
					ipdst[l].re = ipdst[l].re * 2 / sample;
					ipdst[l].im = ipdst[l].im * 2 / sample;
					ipdst[l + sample / 2].re = 0.0;
					ipdst[l + sample / 2].im = 0.0;
				}
				ipdst[0].re /= 2;
				ipdst[0].im /= 2;

				//do IFFT
				ippsFFTInv_CToC_32fc(ipdst, ipdst2, speci, workbufi);
				
				//save
				for (int l = 0; l < 4 * sample; ++l){
					elere[i][j][k][l] = ipdst2[l].re;
					eleim[i][j][k][l] = ipdst2[l].im;
				}

			}
		}
	}
	//free IPP array
	ippsFree(ipsrc);
	ippsFree(ipdst);
	ippsFree(ipdst2);
	//free RF
	vector<vector<vector<vector<short>>>>().swap(RF);
	

	/* interpolation */
	cout << "interpolating...\n";

	vector<vector<vector<float>>> RFre(frame, vector<vector<float>>(line, vector<float>(sample, 0)));
	vector<vector<vector<float>>> RFim(frame, vector<vector<float>>(line, vector<float>(sample, 0)));

	//calculate delay
	const float c0 = 1500.0;
	float ratio;
	int dpoint;
	float shift, cendep;
	float smpt = 1.0 / frq_s; //us
	const complex<float> imunit(0.0f, 1.0f);
	vector<vector<float>> delay(line, vector<float>(ch, 0));
	vector<float> xi(ch, 0);
	for (int i = 0; i < ch; ++i)
		xi[i] = 0.2 * (47.5 - i);
	vector<float> theta(line, 0);
	for (int i = 0; i < line; ++i)
		theta[i] = max_angle * ((line - 1) / 2 - i) * (M_PI / 180.0);
	for (int i = 0; i < line; ++i)
		for (int j = 0; j < ch; ++j)
			delay[i][j] = 4 * frq_s * 1e+6 * (sqrt(pow(focus_first, 2) + pow(xi[j], 2)
			- 2 * focus_first * xi[j] * sin(theta[i])) - focus_first) / (c0*1e+3); //sampling point
	vector<float> eledep(sample, 0);
	for (int i = 0; i < sample; ++i)
		eledep[i] = i * (c0 / 2 / frq_s); //um

	//addition
	for (int i = 0; i < frame; ++i){
		for (int j = 0; j < line; ++j){
			for (int k = 0; k < sample; ++k){
				for (int l = 0; l < ch; ++l){
					
				}
			}
		}
	}
	
	
	
	//free eledata
	vector<vector<vector<vector<float>>>>().swap(elere);
	vector<vector<vector<vector<float>>>>().swap(eleim);

	return 0;
}

