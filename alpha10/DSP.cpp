#include "stdafx.h"

using namespace std;

void writev(vector<float>& v, string& str)
{
	ofstream fout(str, ios_base::out);
	for (auto it = v.begin(); it != v.end(); ++it){
		fout << distance(v.begin(), it) << " " << *it << "\n";
	}
	fout.close();
}

void writev(vector<float>& v, float& sc, string& str)
{
	ofstream fout(str, ios_base::out);
	for (auto it = v.begin(); it != v.end(); ++it){
		fout << distance(v.begin(), it) / sc << " " << *it << "\n";
	}
	fout.close();
}

vector<float> NormCC(const vector<float>& x1, const vector<float>& x2, int lag, int offset)
{
	IppStatus status;
	IppEnum NormA = (IppEnum)(ippAlgAuto | ippsNormA);
	int bufsize = 0;
	Ipp8u *pbuffer;
	const int src1Len = x1.size();
	const int src2Len = x2.size();
	const int dstLen = lag;
	Ipp32f *pSrc1 = ippsMalloc_32f(src1Len);
	Ipp32f *pSrc2 = ippsMalloc_32f(src2Len);
	Ipp32f *pDst = ippsMalloc_32f(dstLen);

	status = ippsCrossCorrNormGetBufferSize(src1Len, src2Len, dstLen, offset, ipp32f, NormA, &bufsize);

	pbuffer = ippsMalloc_8u(bufsize);
	Ipp32f stdev1, stdev2;
	ippsStdDev_32f(pSrc1, src1Len, &stdev1, ippAlgHintFast);
	ippsStdDev_32f(pSrc2, src2Len, &stdev2, ippAlgHintFast);

	status = ippsCrossCorrNorm_32f(pSrc1, src1Len, pSrc2, src2Len, pDst, dstLen, offset, NormA, pbuffer);
	
	vector<float> dst(dstLen, 0);
	for (int i = 0; i < dstLen; ++i){
		if (stdev1 * stdev2 >= 1e-9)
			dst[i] = pDst[i] / (stdev1 * stdev2);
		else
			dst[i] = 0.0;
	}
	ippsFree(pbuffer);
	return dst;
}

vector<float> pwrspe(const vector<float>& v, const int& order)
{
	int length = 1;
	for (int i = 0; i < order; ++i)
		length *= 2;
	vector<float> dst(length, 0);
	
	Ipp8u *sbuf, *ibuf, *wbuf;
	int s_s, s_i, s_w;
	IppsFFTSpec_C_32fc *sp = 0;
	Ipp32fc *src = ippsMalloc_32fc(length);
	ippsFFTGetSize_C_32fc(order, IPP_FFT_NODIV_BY_ANY, ippAlgHintNone, &s_s, &s_i, &s_w);
	sbuf = ippsMalloc_8u(s_s);
	ibuf = ippsMalloc_8u(s_i);
	wbuf = ippsMalloc_8u(s_w);
	ippsFFTInit_C_32fc(&sp, order, IPP_FFT_NODIV_BY_ANY, ippAlgHintNone, sbuf, ibuf);
	for (int i = 0; i < v.size(); ++i)
		src[i].re = v[i];
	for (int i = 0; i < length - v.size(); ++i)
		src[i + v.size()].re = 0.0;
	ippsFFTFwd_CToC_32fc(src, src, sp, wbuf);
	
	for (int i = 0; i < length; ++i){
		dst[i] = 10.0 * log10(pow(src[i].re, 2) + pow(src[i].im, 2));
	}
	
	return dst;
}