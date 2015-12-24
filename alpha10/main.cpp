//main.cpp

//author: Gen Onodera

#include "stdafx.h"
#include "fileopen.h"

using namespace std;
int physio();
void Bsector(const vector<vector<float>>& env, float dangle);
void cairo(const vector<vector<float>>& env, float dangle);
//void BSector2(const vector<vector<float>>& env, float dangle, float fs);
//void Bsector3(const vector<vector<float>>& env, float dangle);
vector<short> ECG;
vector<short> PCG_min;
vector<short> PCG_max;
//const vector<vector<float>>& env, float dangle
int _tmain(int argc, _TCHAR* argv[])
{
	//cairo();
	
	
	physio();

	/* open data */
	cout << "Load started.\n";
	//string filename = "D:/RFdata/study/20151120/sector/RF20151120150735.crf";
	string filename = "D:/RFdata/study/20151026/sample1026.crf";
	a10 raw(filename);
	raw.loadheader();
	raw.printheader();

	raw.frame = 1;
	unsigned short frame = raw.frame;
	unsigned short line = raw.line;
	unsigned short sample = raw.sample;
	unsigned short ch = raw.ch;
	float max_angle = raw.max_angle;
	float frq_t = raw.frq_t;
	float frq_r = raw.frq_r;
	float frq_s = raw.frq_s;
	float FR = raw.FR;


	int physio_offset = 1000 / FR;
	ECG.erase(ECG.begin(), ECG.begin() + physio_offset);
	PCG_min.erase(PCG_min.begin(), PCG_min.begin() + physio_offset);
	PCG_max.erase(PCG_max.begin(), PCG_max.begin() + physio_offset);
	ofstream pecg("./ECG.dat", ios_base::out);
	ofstream ppcgmin("./PCG_min.dat", ios_base::out);
	ofstream ppcgmax("./PCG_max.dat", ios_base::out);
	ofstream framepoint("./fp.dat", ios_base::out);
	for (int i = 0; i < ECG.size(); ++i){
		pecg << i * (1000 / 998) << " " << ECG[i] << "\n";
		ppcgmin << i * (1000 / 998) << " " << PCG_min[i] << "\n";
		ppcgmax << i * (1000 / 998) << " " << PCG_max[i] << "\n";
	}
	for (int i = 0; i < frame; ++i)
		framepoint << i * (1000 / FR) << " 250\n";
	pecg.close();
	ppcgmin.close();
	ppcgmax.close();
	framepoint.close();
	/* load RF data */
	// RF[frame][line][ch][sample]
	cout << "initializing array...\n";
	short tmp = 0;
	
	///*frequency analytic*/
	//float han[60];
	//for (int i = 0; i < 60; ++i)
	//	han[i] = 0.5 - 0.5 * cos(2.0 * M_PI * i / 60);

	//Ipp8u *sbuf, *ibuf, *wbuf;
	//int s_s, s_i, s_w;
	//IppsFFTSpec_C_32fc *sp = 0;
	//Ipp32fc *src = ippsMalloc_32fc(1024);
	//ippsFFTGetSize_C_32fc(10, IPP_FFT_NODIV_BY_ANY, ippAlgHintNone, &s_s, &s_i, &s_w);
	//sbuf = ippsMalloc_8u(s_s);
	//ibuf = ippsMalloc_8u(s_i);
	//wbuf = ippsMalloc_8u(s_w);
	//ippsFFTInit_C_32fc(&sp, 10, IPP_FFT_NODIV_BY_ANY, ippAlgHintNone, sbuf, ibuf);
	//for (int i = 0; i < 60; ++i)
	//	src[i].re = RF[60][17][90][i + 1270] * han[i];
	//for (int i = 0; i < 1024 - 60; ++i)
	//	src[i + 60].re = 0.0;
	//ippsFFTFwd_CToC_32fc(src, src, sp, wbuf);
	//float freqi, tmpfc;
	//ofstream spect("./spectrum.dat", ios_base::out);
	//for (int i = 0; i < 512; ++i){
	//	freqi = static_cast<float>(i) * 15 / 512;
	//	tmpfc = 10.0 * log10(pow(src[i].re, 2) + pow(src[i].im, 2));
	//	spect << freqi << " " << tmpfc << "\n";
	//}
	//spect.close();

	///*calculate differences between frames before and after*/
	//short tmps, th;
	//th = 500; // threshold
	//for (int i = 0; i < 10; ++i){
	//	for (int j = 0; j < line; ++j){
	//		for (int k = 0; k < ch; ++k){
	//			for (int l = 0; l < sample - 1; ++l){
	//				tmps = abs(RF[i + 55][j][k][l] - RF[i + 54][j][k][l]);
	//				if (tmps > th && i==5)
	//					cout << "over " << th << ":(" << i + 55 << " " << j << " " << k << " " << l << ")\n";
	//			}
	//		}
	//	}
	//}

	raw.loadRF();
	//raw.loadRF0(60);
	/* do FFT and IFFT */
	/*string fname2 = "./raw.dat";
	ofstream fout2(fname2, ios_base::out);
	for (int i = 0; i < ch; ++i){
		for (int j = 0; j < sample - 1; ++j)
			fout2 << j << " " << raw.RF0[30][i][j] - i * 4096 << "\n";
		fout2 << "\n";
	}
	fout2.close();*/



	cout << "creating analytic signal...\n";
	
	//initialize focused RF array(vector)
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
	Ipp32fc *ipsrc = ippsMalloc_32fc((int)sample);
	Ipp32fc *ipdst = ippsMalloc_32fc((int)sample);
	Ipp32fc *ipsrc2 = ippsMalloc_32fc((int)(4 * sample));
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
				ippsZero_32fc(ipsrc, sample);
				ippsZero_32fc(ipdst, sample);
				ippsZero_32fc(ipsrc2, 4 * sample);
				ippsZero_32fc(ipdst2, 4 * sample);
				//set
				for (int l = 0; l < sample - 1; ++l)
					ipsrc[l].re = raw.RF[i][j][k][l];
					//ipsrc[l].re = raw.RF[60][j][k][l] - raw.RF[59][j][k][l];
 
				//do FFT
				ippsFFTFwd_CToC_32fc(ipsrc, ipdst, specf, workbuff);
				ippsZero_8u(workbuff, size_workf);
				//double positive part and delete negative part
				for (int l = 0; l < sample / 2; ++l){
					ipdst[l].re = ipdst[l].re * 2 / sample;
					ipdst[l].im = ipdst[l].im * 2 / sample;
					ipdst[l + sample / 2].re = 0.0;
					ipdst[l + sample / 2].im = 0.0;
				}
				for (int l = 0; l < 34; ++l){
					ipdst[l].re = 0.0;
					ipdst[l].im = 0.0;
				}
				
				
				//ipdst[0].re /= 2;
				//ipdst[0].im /= 2;

				for (int l = 0; l < sample; ++l){
					ipsrc2[l].re = ipdst[l].re;
					ipsrc2[l].im = ipdst[l].im;
				}


				//do IFFT
				ippsFFTInv_CToC_32fc(ipsrc2, ipdst2, speci, workbufi);
				ippsZero_8u(workbufi, size_worki);
				
				//save
				for (int l = 0; l < 4 * sample; ++l){
					elere[i][j][k][l] = 4 * sample * ipdst2[l].re;
					eleim[i][j][k][l] = 4 * sample * ipdst2[l].im;
				}
				ippsZero_32fc(ipdst2, 4 * sample);
			}
		}
	}

	//free IPP array
	ippsFree(ipsrc);
	ippsFree(ipdst);
	ippsFree(ipdst2);
	//free RF
	//vector<vector<vector<vector<short>>>>().swap(RF);
	raw.freeRF();

	string fname = "./ele.dat";
	ofstream fout(fname, ios_base::out);
	for (int i = 0; i < 4 * sample; ++i)
		fout << i << " " << elere[0][0][0][i] << "\n";
	fout.close();


	/* interpolation */
	cout << "interpolating...\n";

	vector<vector<vector<float>>> RFre(frame, vector<vector<float>>(line, vector<float>(sample, 0)));
	vector<vector<vector<float>>> RFim(frame, vector<vector<float>>(line, vector<float>(sample, 0)));

	//calculate delay
	const float c0 = 1540.0;
	int point, add;
	float eledep; //in-bound(um)
	float decimal; //decimal part in sampling point of round trip distance
	vector<float> xi(ch, 0); // x-coordinate of each element
	for (int i = 0; i < ch; ++i)
		xi[i] = 0.2 * (47.5 - i) * 1e+3; //um
	vector<float> theta(line, 0);
	for (int i = 0; i < line; ++i) //beam angle
		theta[i] = max_angle * ((line - 1) / 2 - i) * (M_PI / 180.0);
	vector<float> cendep(sample, 0);
	for (int i = 0; i < sample; ++i)
		cendep[i] = i * (c0 / (2 * frq_s)); //out-bound(um)

	//addition
	for (int i = 0; i < frame; ++i){
		for (int j = 0; j < line; ++j){
			for (int k = 0; k < sample; ++k){
				add = 0;
				for (int l = 0; l < ch; ++l){
					eledep = sqrt(pow(xi[l], 2) + pow(cendep[k], 2) - 2 * xi[l] * cendep[k] * sin(theta[j]));
					point = static_cast<int>(((cendep[k] + eledep) / 2) / (c0 / (8 * frq_s)));
					decimal = ((cendep[k] + eledep) / 2) / (c0 / (8 * frq_s)) - point;
					if (point < 4 * sample - 1){
						RFre[i][j][k] += (elere[i][j][l][point] + (elere[i][j][l][point + 1] - elere[i][j][l][point]) * decimal);
						RFim[i][j][k] += (eleim[i][j][l][point] + (eleim[i][j][l][point + 1] - eleim[i][j][l][point]) * decimal);
						//RFre[i][j][k] += elere[i][j][l][point];
						//RFim[i][j][k] += eleim[i][j][l][point];
						++add;
					}
				}
				if (add != 0){
					RFre[i][j][k] /= add;
					RFim[i][j][k] /= add;
				}
				else{
					RFre[i][j][k] = 0.0;
					RFim[i][j][k] = 0.0;
				}
			}
		}
	}

	//free eledata
	vector<vector<vector<vector<float>>>>().swap(elere);
	vector<vector<vector<vector<float>>>>().swap(eleim);

	vector<vector<vector<float>>> env(frame, vector<vector<float>>(line, vector<float>(sample, 0)));
	for (int i = 0; i < frame; ++i)
		for (int j = 0; j < line; ++j)
			for (int k = 0; k < sample; ++k)
				env[i][j][k] = sqrt(pow(RFre[i][j][k], 2) + pow(RFim[i][j][k], 2));

	/*fout.open("./env.dat", ios_base::out);
	for (int i = 0; i < line; ++i){
		for (int j = 0; j < sample; ++j)
			fout << i << " " << j << " " << env[0][i][j] << "\n";
		fout << "\n";
	}
	fout.close();*/

	//Bsector(env[0], max_angle);
	//BSector2(env[0], max_angle, frq_s);
	//Bsector3(env[0], max_angle);
	cairo(env[0], max_angle);
	
	return 0;
}

