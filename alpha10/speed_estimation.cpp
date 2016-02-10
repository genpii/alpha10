#include "stdafx.h"
#include "header.h"

using namespace std;

est_ss::est_ss(int w, int h, int beam)
{
	grid_w = w;
	grid_h = h;
	bpg = beam;
}

void est_ss::set_parameter(double dangle, double frq_s)
{
	angle = dangle;
	fs = frq_s;

	line = grid_w * bpg;
	ch = 96;

	elex = vector<double>(ch, 0);
	for (int i = 0; i < ch; ++i){
		elex[i] = 0.2 * (47.5 - i) * 1e+3;
	}
	theta = vector<double>(line, 0);
	for (int i = 0; i < line; ++i){
		theta[i] = angle * ((line - 1) / 2 - i) * (M_PI / 180.0);
	}
	phi = vector<double>(grid_w - 1, 0);
	for (int i = 0; i < grid_w - 1; ++i){
		phi[i] = ((double)(grid_w / 2) - i - 1) * angle * bpg * (M_PI / 180.0);
	}

}

void est_ss::loadRF(const vector<vector<vector<short>>>& RF0)
{
	int line_size = RF0.size();
	sample = RF0[0][0].size();
	++sample;
	
	int start = (line_size - line) / 2;
	RF = vector<vector<vector<short>>>(line, vector<vector<short>>(ch, vector<short>(sample, 0)));

	for (int i = 0; i < line; ++i){
		RF[i] = RF0[start + i];
	}

}

int est_ss::calc_delay(double depth)
{
	//reference : 95-th element
	double dep = depth * 1e+3; //mm -> um

	delay = vector<vector<double>>(line, vector<double>(ch, 0));
	scatdep = vector<double>(line, 0);
	
	int point_ref;
	double time_ref;
	double dep_ref;
	
	

	ofstream fout("corrcheck.dat", ios_base::out);

	/*Ipp64f std1, std2;
	IppStatus status;
	IppEnum Norm = (IppEnum)(ippAlgAuto | ippsNormNone);
	Ipp8u *pbuffer;
	int bufsize = 0;
	int lowlag;
	const int src1Len = 64;
	const int src2Len = sample;
	const int dstLen = 128;
	Ipp64f *pSrc1 = ippsMalloc_64f(src1Len);
	Ipp64f *pSrc2 = ippsMalloc_64f(src2Len);
	Ipp64f *pDst = ippsMalloc_64f(dstLen);
	Ipp64f maxval;
	int maxidx;*/

	//for (int i = 0; i < line; ++i){
	//	point of depth at 95-th element
	//	dep_ref = sqrt(pow(dep, 2) + pow(elex[ch - 1], 2) - 2 * dep * elex[ch - 1] * sin(theta[i]));
	//	point_ref = static_cast<int>(fs * (dep + dep_ref) / c0);

	//	vector<short>::iterator it = RF[i][ch - 1].begin();
	//	it += point_ref;
	//	auto it2 = max_element(it, RF[i][ch - 1].end());
	//	point_ref = distance(RF[i][ch - 1].begin(), it2);

	//	scatdep[i] = (pow(c0 * point_ref / fs, 2) - pow(elex[ch - 1], 2)) / 2 / (c0 * point_ref / fs - elex[ch - 1] * sin(theta[i])); //um

	//	if (point_ref + src1Len / 2  > sample){
	//		cout << "border error!\n";
	//		return 1;
	//	}
	//	time_ref = point_ref / fs;
	//	
	//	ippsZero_64f(pSrc1, src1Len);
	//	for (int j = 0; j < src1Len; ++j){
	//		pSrc1[j] = RF[i][ch - 1][point_ref + j - src1Len / 2];
	//	}
	//	ippsNorm_L2_64f(pSrc1, src1Len, &std1);

	//	lowlag = point_ref - src1Len;
	//	status = ippsCrossCorrNormGetBufferSize(src1Len, src2Len, dstLen, lowlag, ipp64f, Norm, &bufsize);
	//	pbuffer = ippsMalloc_8u(bufsize);
	//	

	//	for (int j = 0; j < ch; ++j){
	//		ippsZero_64f(pSrc2, sample);
	//		for (int k = 0; k < sample - 1; ++k){
	//			pSrc2[k] = RF[i][j][k];
	//		}

	//		status = ippsCrossCorrNorm_64f(pSrc1, src1Len, pSrc2, src2Len, pDst, dstLen, lowlag, Norm, pbuffer);
	//		
	//		for (int k = 0; k < dstLen; ++k){
	//			ippsNorm_L2_64f(pSrc2 + lowlag, src1Len, &std2);
	//			pDst[k] /= (std1 * std2);
	//		}
	//		
	//		status = ippsMaxIndx_64f(pDst, dstLen, &maxval, &maxidx);

	//		delay[i][j] = time_ref + (maxidx - dstLen / 2) / fs;
	//		if (i == line / 2){
	//			fout << j << " " << delay[i][j] << "\n";
	//		}
	//	}

	Ipp64f std1, std2;
	IppStatus status;
	IppEnum Norm = (IppEnum)(ippAlgAuto | ippsNormNone);
	Ipp8u *pbuffer;
	int bufsize = 0;
	int lowlag;
	const int src1Len = 128;
	const int src2Len = sample;
	const int dstLen = 128;
	Ipp32f *pSrc1 = ippsMalloc_32f(src1Len);
	Ipp32fc *pSrc1c = ippsMalloc_32fc(src1Len);
	Ipp32f *pSrc2 = ippsMalloc_32f(src2Len);
	Ipp32fc *pSrc2c = ippsMalloc_32fc(src2Len);
	Ipp32fc *pDst = ippsMalloc_32fc(dstLen);
	Ipp32f *pDstamp = ippsMalloc_32f(dstLen);
	Ipp32f maxval;
	int maxidx;
	IppsHilbertSpec_32f32fc *hspec;

	for (int i = 0; i < line; ++i){
		//point of depth at 95-th element
		dep_ref = sqrt(pow(dep, 2) + pow(elex[ch - 1], 2) - 2 * dep * elex[ch - 1] * sin(theta[i]));
		point_ref = static_cast<int>(fs * (dep + dep_ref) / c0);

		vector<short>::iterator it = RF[i][ch - 1].begin();
		it += point_ref;
		auto it2 = max_element(it, RF[i][ch - 1].end());
		point_ref = distance(RF[i][ch - 1].begin(), it2);

		scatdep[i] = (pow(c0 * point_ref / fs, 2) - pow(elex[ch - 1], 2)) / 2 / (c0 * point_ref / fs - elex[ch - 1] * sin(theta[i])); //um

		if (point_ref + src1Len / 2  > sample){
			cout << "border error!\n";
			return 1;
		}
		time_ref = point_ref / fs;

		ippsZero_32f(pSrc1, src1Len);
		for (int j = 0; j < src1Len; ++j){
			pSrc1[j] = RF[i][ch - 1][point_ref + j - src1Len / 2];
		}
		ippsHilbertInitAlloc_32f32fc(&hspec, src1Len, ippAlgHintFast);
		ippsHilbert_32f32fc(pSrc1, pSrc1c, hspec);
		ippsNorm_L2_32fc64f(pSrc1c, src1Len, &std1);

		lowlag = point_ref - src1Len;
		status = ippsCrossCorrNormGetBufferSize(src1Len, src2Len, dstLen, lowlag, ipp32fc, Norm, &bufsize);
		pbuffer = ippsMalloc_8u(bufsize);


		for (int j = 0; j < ch; ++j){
			ippsZero_32f(pSrc2, sample);
			for (int k = 0; k < sample - 1; ++k){
				pSrc2[k] = RF[i][j][k];
			}
			ippsHilbertInitAlloc_32f32fc(&hspec, src2Len, ippAlgHintFast);
			ippsHilbert_32f32fc(pSrc2, pSrc2c, hspec);
			status = ippsCrossCorrNorm_32fc(pSrc1c, src1Len, pSrc2c, src2Len, pDst, dstLen, lowlag, Norm, pbuffer);

			for (int k = 0; k < dstLen; ++k){
				ippsNorm_L2_32fc64f(pSrc2c + lowlag, src1Len, &std2);
				pDst[k].re /= (Ipp32f)(std1 * std2);
				pDst[k].im /= (Ipp32f)(std1 * std2);
			}
			ippsMagnitude_32fc(pDst, pDstamp, dstLen);
			status = ippsMaxIndx_32f(pDstamp, dstLen, &maxval, &maxidx);

			delay[i][j] = time_ref + (maxidx - dstLen / 2) / fs;
			if (i == line / 2){
				fout << j << " " << delay[i][j] << "\n";
			}
		}


		ippsFree(pbuffer);
	}

	fout.close();

	auto minmax = minmax_element(scatdep.begin(), scatdep.end());
	cout << "min depth :" << *minmax.first << "\n";
	cout << "max depth :" << *minmax.second << "\n";

	return 0;
}

int est_ss::calc_path()
{
	double eledep;
	//int N = grid_h * grid_w;
	path = vector<vector<vector<vector<double>>>>(line, vector<vector<vector<double>>>(ch, vector<vector<double>>(grid_w, vector<double>(grid_h, 0))));
	path_near = vector<vector<double>>(line, vector<double>(ch, 15000.0));
	
	vector<double> b(ch, 0);
	for (int i = 0; i < ch; ++i){
		b[i] = 0.2 * (47.5 - i) * 1e+3;
	}
	vector<vector<double>> a(line, vector<double>(ch, 0));
	for (int i = 0; i < line; ++i){
		for (int j = 0; j < ch; ++j){
			a[i][j] = scatdep[i] * cos(theta[i]) / (scatdep[i] * sin(theta[i]) - b[j]);
		}
	}

	vector<double> r(grid_h, 0);
	for (int i = 0; i < grid_h; ++i){
		r[i] = 35000.0 - i * 10000.0;
	}

	for (int i = 0; i < line; ++i){
		double scatdep_x = scatdep[i] * sin(theta[i]);
		double scatdep_y = scatdep[i] * cos(theta[i]);
		pair<double, pair<double, bool>> scat = make_pair(scatdep_y, make_pair(scatdep_x, false));

		for (int j = 0; j < ch; ++j){
			
			//out-bound (full path: scatdep[i]) (near field is initialized)
			for (int k = 0; k < grid_h - 1; ++k){
				path[i][j][i / bpg][k] = 10000.0;
			}
			path[i][j][i / bpg][grid_h - 1] = scatdep[i] - 35000.0;

			double ab = a[i][j] * b[j];
			double sqa = sqrt(1 + pow(a[i][j], 2));
			double alpha;
			double root;
			double axytemp;
			//in-bound (full path: eledep)
			eledep = sqrt(pow(scatdep[i] * sin(theta[i]) - b[j], 2) + pow(scatdep[i] * cos(theta[i]), 2));
			bool flag = (scatdep_x > elex[j]); //slope
			// second.first : x, first : y, bool : second.second
			vector<pair<double, pair<double, bool>>> axial(grid_h);
			for (int k = 0; k < grid_h; ++k){
				/*axxtemp = (pow(a[i][j], 2) * b[j] + sqrt((1 + pow(a[i][j], 2)) * pow(r[k], 2) - pow(a[i][j] * b[j], 2))) / (1 + pow(a[i][j], 2));
				if ((axxtemp > b[j] && a[i][j] > 0.0) || (axxtemp <= b[j] && a[i][j] <= 0.0)){
					axial[k].second.first = axxtemp;
				}
				else{
					axial[k].second.first = (pow(a[i][j], 2) * b[j] - sqrt((1 + pow(a[i][j], 2)) * pow(r[k], 2) - pow(a[i][j] * b[j], 2))) / (1 + pow(a[i][j], 2));
				}
				axial[k].first = a[i][j] * (axial[k].second.first - b[j]);*/
				//axial[k].first = a[i][j] * (sqrt((1 + pow(a[i][j], 2) * pow(r[k], 2) - pow(a[i][j] * b[j], 2))) - b[j]) / (1 + pow(a[i][j], 2));
				//root = a[i][j] * (sqrt(pow(r[k] * sqa, 2) - pow(ab, 2)) - b[j]) / r[k] / sqa / sqa;
				//if (pow(r[k], 2) >= pow(ab, 2) || b[j] >= 0.0){
				//	//axial[k].first = root * r[k];

				//}
				//else{
				//	alpha = asin(root);
				//	if (a[i][j] > 0.0){
				//		alpha = M_PI - alpha;
				//	}
				//	else{
				//		alpha = - M_PI - alpha;
				//	}
				//	axial[k].first = r[k] * sin(alpha);
				//}
				double check = sqrt((1 + pow(a[i][j], 2)) * pow(r[k], 2) - pow(ab, 2));
				axytemp = -a[i][j] * (b[j] - sqrt((1 + pow(a[i][j], 2)) * pow(r[k], 2) - pow(ab, 2))) / (1 + pow(a[i][j], 2));
				if (axytemp < 0.0){
					axytemp = -a[i][j] * (b[j] + sqrt((1 + pow(a[i][j], 2)) * pow(r[k], 2) - pow(ab, 2))) / (1 + pow(a[i][j], 2));
				}
				axial[k].first = axytemp;
				axial[k].second.first = axial[k].first / a[i][j] + b[j];
				axial[k].second.second = false; //axial:false
			}
			vector<pair<double, pair<double, bool>>> lateral(grid_w - 1);
			for (int k = 0; k < grid_w - 1; ++k){
				lateral[k].second.first = a[i][j] * b[j] / (a[i][j] - tan(phi[k]));
				lateral[k].first = a[i][j] * (lateral[k].second.first - b[j]);
				lateral[k].second.second = true; //lateral:true
			}
			sort(lateral.begin(), lateral.end()); // y-sort
			reverse(lateral.begin(), lateral.end()); // reverse
			while (axial.back().first < lateral.front().first){ //interrupt
				axial.insert(axial.end() - 2, lateral[0]);
				sort(axial.begin(), axial.end());
				lateral.erase(lateral.begin());
			}
			axial.insert(axial.begin(), scat);

			double dsum = 0.0;
			double dtmp = 0.0;
			int a_size = axial.size();
			int pos_ini = i / bpg;
			int pos_axial = 0;
			for (int k = 0; k < a_size -1; ++k){
				dtmp = sqrt(pow(axial[k].second.first - axial[k + 1].second.first, 2) + pow(axial[k].first - axial[k + 1].first, 2));
				path[i][j][pos_ini][grid_h - 1 - pos_axial] += dtmp;
				dsum += dtmp;
				if (axial[k + 1].second.second){
					if (flag){
						++pos_ini;
					}
					else{
						--pos_ini;
					}
					--pos_axial;
				}
				++pos_axial;
			}
			path_near[i][j] += (eledep - dsum);


		}
	}

	return 0;
}

int est_ss::del_brokenelement()
{
	int size = b_ele.size();

	for (int i = 0; i < line; ++i){
		for (int j = size - 1; j >= 0; --j){
			delay[i].erase(delay[i].begin() + b_ele[j]);
			path[i].erase(path[i].begin() + b_ele[j]);
			path_near[i].erase(path_near[i].begin() + b_ele[j]);
		}
	}

	ch -= size;

	/*for (int i = 0; i < line; ++i){
		for (int j = 1; j <= 4; ++j){
			for (int k = 0; k < 19; ++k){
				delay[i].erase(delay[i].end() - 1 - j);
				path[i].erase(path[i].end() - 1 - j);
				path_near[i].erase(path_near[i].end() - 1 - j);
			}
		}
	}
	ch = 5;*/

	return 0;
}

void est_ss::SVD()
{
	int matrix_layout = LAPACK_ROW_MAJOR;
	int m = line * ch;
	int n = grid_h * grid_w + 1;
	int lda = n;


	double *A, *d, *e, *tauq, *taup;
	A = new double[m * m];
	int count = 0;
	for (int i = 0; i < line; ++i){
		for (int j = 0; j < ch; ++j){
			A[count] = path_near[i][j];
			++count;
			for (int k = 0; k < grid_w; ++k){
				for (int l = 0; l < grid_h; ++l){
					A[count] = path[i][j][k][l];
					++count;
				}
			}
		}
	}

	d = new double[n];
	e = new double[n - 1];
	taup = new double[n];
	tauq = new double[n];
	for (int i = 0; i < n; ++i){
		d[i] = 0.0;
		taup[i] = 0.0;
		tauq[i] = 0.0;
		if (i != n - 1)
			e[i] = 0.0;
	}
	
	LAPACKE_dgebrd(matrix_layout, (lapack_int)m, (lapack_int)n, A, (lapack_int)lda, d, e, tauq, taup);
	double *Acpy1 = new double[m * n];
	for (int i = 0; i < m * n; ++i){
		Acpy1[i] = A[i];
	}

	LAPACKE_dorgbr(matrix_layout, 'Q', m, m, n, A, lda, tauq);

	LAPACKE_dorgbr(matrix_layout, 'P', n, n, m, Acpy1, lda, taup);



	LAPACKE_dbdsqr(matrix_layout, 'U', n, n, m, 0, d, e, Acpy1, n, A, m, 0, 0);

	int q = 0;
}

int est_ss::write_mat()
{
	//int count = 0;
	ofstream fout("path.dat", ios_base::out);
	fout.precision(10);
	for (int i = 0; i < line; ++i){
		for (int j = 0; j < ch; ++j){
			fout << path_near[i][j] << " ";
			for (int k = 0; k < grid_w; ++k){
				for (int l = 0; l < grid_h; ++l){
					fout << path[i][j][k][l];
					if (k != grid_w - 1 || l != grid_h - 1){
						fout << " ";
					}
				}
			}
			fout << "\n";
		}
	}
	fout.close();

	fout.open("delay.dat", ios_base::out);

	for (int i = 0; i < line; ++i){
		for (int j = 0; j < ch; ++j){
			fout << delay[i][j] << "\n";
		}
	}
	fout.close();

	return 0;
}

void est_ss::sim_RF()
{
	double r = 45000.0;
	double p = 200.0;
	double l = 9500.0;
	vector<vector<double>> time(line, vector<double>(ch, 0));
	for (int i = 0; i < line; ++i){
		for (int j = 0; j < ch; ++j){
			time[i][j] = (r + sqrt(pow(r, 2) + pow(elex[j], 2) - 2 * r * elex[j] * sin(theta[i]))) / c0;
		}
	}
	int pulse_wid = 3 * fs / 3.5;
	short *pulse = new short[pulse_wid];
	for (int i = 0; i < pulse_wid; ++i){
		pulse[i] = (short)(1e+4 * sin((double)(i * (fs / 3.5) / 2 / M_PI)) * (0.5 - 0.5*cos(2 * M_PI * i / pulse_wid)));
	}

	RF = vector<vector<vector<short>>>(line, vector<vector<short>>(ch, vector<short>(sample, 0)));

	int pulse_pnt;
	for (int i = 0; i < line; ++i){
		for (int j = 0; j < ch; ++j){
			for (int k = 0; k < pulse_wid; ++k){
				pulse_pnt = (int)(time[i][j] * fs);
				RF[i][j][k + pulse_pnt] = pulse[k];
			}
		}
	}
	int a;
}

void est_ss::sim_delay_setspeed()
{
	sp_g = vector<vector<double>>(grid_w, vector<double>(grid_h, 1540.0));
	sp_n = 1540.0;
	sp_n += 10.0;
	sp_g[0][0] -= 20.0; sp_g[0][1] += 10.0; sp_g[0][2] -= 30.0;
	sp_g[1][0] += 30.0; sp_g[1][1] -= 10.0; sp_g[1][2] += 20.0;
	sp_g[2][0] -= 30.0; sp_g[2][1] += 0.0; sp_g[2][2] += 10.0;
	sp_g[3][0] += 20.0; sp_g[3][1] -= 20.0; sp_g[3][2] -= 10.0;
	sp_g[4][0] += 0.0; sp_g[4][1] += 30.0; sp_g[4][2] -= 20.0;
	
	delay = vector<vector<double>>(line, vector<double>(ch, 0));
	double r = 45000.0;
	scatdep = vector<double>(line, r);
}

void est_ss::sim_delay_calcdelay()
{
	double dtmp;
	
	for (int i = 0; i < line; ++i){
		for (int j = 0; j < ch; ++j){
			dtmp = path_near[i][j] / sp_n;
			for (int k = 0; k < grid_w; ++k){
				for (int l = 0; l < grid_h; ++l){
					dtmp += path[i][j][k][l] / sp_g[k][l];
				}
			}
			delay[i][j] = dtmp;
		}
	}


}