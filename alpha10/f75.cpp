#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>
#include <windows.h>
#include "calc.h"
#include "graphics.h"

#define  Movie    0
#define  MaxStr   0.05
#define  cStrain  1
#define  Timo     0
#define  DCor     0
#define  DispTr   1
#define  Compliance 0

#define  itrace   1
#define  ipause   0
#define  ihosei   1
#define  aux_flg  1
#define  SPAv     1

#define  Sens     0.20 /* Sens=15 (SSD-6500) */
#define  Pos      -32463

#define  izoom_B  2 /* scale factor for displaying B-mode (default: 2) */
#define  izoom_M  1 /* scale factor for displaying M-mode (default: 2) */

#define  mpp      0.1e-3 /* [m/pixel] in B(M)-mode when izoom_B(M)=1 */
/* default:0.1e-3 */
#define  nwv      4      /* nwv*wavelength=width of correlation window */
#define  nDFT    20      // averaging length of quadrature detection
#define  nfrq     4      // frequency of quadrature detection

#define  nfram0   1024
#define  nposi0   64
#define  ndep0    3000

#define  maxfrm   70

#define  B_gain   60
#define  rad      6.283185307179586   // 2*pi
#define  c0       1540.0        // in vivo

char           cdat;
unsigned char  ucdat;
unsigned short usdat;
unsigned long  uldat, tcnt;
float          fdat;

int            i, j, k, l, m, irtn, icnt, icnt2, icnt3, nnn, idmy1, idmy2, isbt, iflg, ic;
int            ii1, ifrm1, ifrm2, ifrm, iiposi, numdepth, ndft, mar_dft, ispc, nspc;
int            ip_s, ip_d, iage, icfile, nfhead, islct00, iseq, iseq1, iseq2;
int            ibx = 0, iby = 0, nx0 = 0, ny0, ipeaksum, iPRF, layer_th, izB, anabt, imesh;
int            ip1, ip2, ip3, ip4, dsh, Ecnt, iave, nslp;
static int     peak0[ndep0];
static int     ipeakflg[nposi0] = {};
static int     ipeak11[nposi0], ipeak22[nposi0];
static int     ipeak33[nposi0], ipeak44[nposi0];
int            nframe, nframe_b, nfr_low, nposi, ncount, ndmy, nwb, wwB;
int            nbeat, nbeat1, nbeat2, itrk1, itrk2, ipeak, ipeak2;
int            iadd, iaddp, dadd, iaddmin, addmin, iaux, iPCGmax, npb;
int            Bx_ini, By_ini, Mx_ini, My_ini, size;
int            ifs, ilevel, heartrate, icgain;
int            lp2s, lp2q, ix0, iy0, MpF;
static int     ihof[nposi0][nfram0], track0[maxfrm];
static int     ndep[nposi0];
short          idata, idata1, idata2;
unsigned long  wid;
float          rx0, ry0, ry0_i, deltad, beat, tbeat, period, elp;
float          PRF, PRF2, f0, f0i, fs, rmaxz, ecgmax, ecgbias, pcgmax, pcgmin, pcgbias;
float          rryy, rmaxdh, rcorr, rwb, vest00, aveV;
float          Emean, Evar;
float          maxdh, mindh, rrdh, hhh, hhh_n, rho, elt0, p_s, p_d;
float          rr1, rr2, rr3, rr4, rr5, mean1, mean2, var1, var2, var12, aaa, ddmax;
float          *sgn = NULL, *dmy = NULL;
static float   vel[ndep0], vel_p[ndep0], dh[ndep0][nfram0];
static float   disp[ndep0][nfram0], dd[nfram0];
static float   wary[500], han_w[nDFT + 1];
unsigned char  rf1, rf2;
FILE           *fin, *fout;
int            *raw_dat;

unsigned short pmax, pmin, emax, emin, plmax, plmin;
int     idltx, idltx2, idcnt, icntecg, icntecg2;
float   rr, ri, rr0, alpmin, absbeta, alpha, thetaf;
static float   xpeak0[ndep0], ecg[5 * nfram0], pcg[5 * nfram0];
static float   qdr_p[ndep0][3], qdi_p[ndep0][3];
static float   qdr_f[ndep0][3], qdi_f[ndep0][3];
fcomplex cdmy, cdmyp, cbeta, ctemp, cdmyp2, cdmyph, cdmy2, cdmyh;
fcomplex crr0, crr1, crr12, czeros, cvar12, zr1, zr2;

COLORREF gscolor;

char        fname[128], fall[256], name[64], buf[512];
static char *fhead[] = { "C:/Users/user/Desktop/etc/Lab/vs/research/data/" };

int trackingpoint_top(), trackingpoint_bot();

int main(void)
{
	/* *************************************************************************
		Create Window
		********************************************************************** */
	wwB = 1000;
	wid = pageb(100, 0, 900, wwB, "B-mode");
	start();
	Bx_ini = 100; By_ini = 100;
	/**************************************************************************
	Input Serial No.
	********************************************************************** */
	//printf("input beginning and end of #serial\n");
	//scanf("%d %d", &iseq1);
	iseq1 = 101;
	printf("iseq1 = %d \n", iseq1);
	printf("*******************************\n");
	printf("*       start analysing       *\n");
	printf("*******************************\n");
	/* *************************************************************************
	Reading Parameters from File
	********************************************************************** */
start:;
	if ((fin = fopen("./ssd6500.coe", "r")) == NULL){
		fprintf(stderr, "cannot open dpman1M.coe\n");
		exit(1);
	}

	icfile = 0;
	while (icfile == 0){
		fscanf(fin, "%d %s %d %d %s %d %f %d %d %d %d",
			&islct00, name, &iseq, &nfhead, fname, &iage, &rryy, &icnt, &ip_s, &ip_d,
			&anabt);
		if (iseq == iseq1) icfile = 1;
		if (islct00 == 9){ fclose(fin); exit(0); }
	}
	fclose(fin);
	rmaxdh = 1e-6*icnt;

	printf("islct00 = %d\n", islct00);
	printf("Name of subject = %s\n", name);
	printf("Serial number = %d\n", iseq);
	printf("Age of subject = %d years old\n", iage);
	printf("Maximal velocity range = %f m/s\n", rryy);
	printf("Systolic & Diastolic blood pressure = %d, %d mmHg\n", ip_s, ip_d);
	/* *************************************************************************
	Read Raw Data from File
	********************************************************************** */
	sprintf(fall, "%s%s.rfa", fhead[nfhead], fname);
	printf("File name = %s\n", fall);
	if ((fin = fopen(fall, "rb")) == NULL){
		printf("cannot open file: %s\n", fall);
		exit(1);
	}

	fseek(fin, 480L, SEEK_SET);
	//fread(buf,sizeof(char),480,fin);
	fread(&fdat, sizeof(float), 1, fin);   /* 送信周波数 [MHz] */
	printf("f0 = %f Hz\n", f0 = 1e6*fdat);
	fread(&fdat, sizeof(float), 1, fin);   /* Burst波数 */
	fread(&usdat, sizeof(unsigned short), 1, fin);
	npb = usdat;
	fread(&fdat, sizeof(float), 1, fin);   /* PRF [kHz] */
	printf("PRF = %f Hz\n", 1e3*fdat);
	fread(&fdat, sizeof(float), 1, fin);   /* 受信サンプリング周波数 [MHz] */
	printf("fs = %f Hz\n", fs = 1e6*fdat);
	fread(&fdat, sizeof(float), 1, fin);   /* ビーム間隔 [mm] */
	printf("beam interval = %f mm\n", elp = fdat);
	fread(&usdat, sizeof(unsigned short), 1, fin);   /* ビーム数 */
	printf("number of beams = %d\n", nposi = usdat);
	fread(&usdat, sizeof(unsigned short), 1, fin);   /* ROI上端深さ(point)=0 */
	printf("record offset = %d\n", usdat);
	fread(&usdat, sizeof(unsigned short), 1, fin);   /* ビーム数当たりサンプル数 */
	printf("number of samples per line = %d\n", ncount = usdat);
	fread(&usdat, sizeof(unsigned short), 1, fin);   /* 開始ビーム番号 */
	fread(&usdat, sizeof(unsigned short), 1, fin);   /* 終了ビーム番号 */
	fread(&usdat, sizeof(unsigned short), 1, fin);   /* フォーカス段数 */
	for (i = 0; i<7; ++i) fread(&fdat, sizeof(float), 1, fin);   /* フォーカス位置設定 */
	for (i = 0; i<4; ++i)
		fread(&usdat, sizeof(unsigned short), 1, fin);   /* フォーカス音圧設定 */
	fread(&fdat, sizeof(float), 1, fin);   /* ステア角 [degree] */
	printf("steering angle: %f\n", fdat);
	fread(&fdat, sizeof(float), 1, fin);   /* フレームレート [Hz] */
	printf("frame rate = %f Hz\n", PRF2 = fdat);
	fread(&fdat, sizeof(float), 1, fin);   /* スケールファクター */
	fread(&uldat, sizeof(unsigned long), 1, fin);   /* フレーム数 */
	printf("number of frames = %d\n", nframe = uldat);
	fread(&usdat, sizeof(unsigned short), 1, fin);   /* データ情報 */
	fread(&fdat, sizeof(float), 1, fin);   /* 受信周波数 [MHz] */
	//printf("receiving frequency: %f MHz\n",fdat);
	for (i = 0; i<14; ++i) fread(&cdat, sizeof(char), 1, fin);   /* 未使用 */
	for (i = 0; i<4; ++i) fread(&usdat, sizeof(unsigned short), 1, fin);   /* 諸設定 */
	for (i = 0; i<16; ++i) fread(&cdat, sizeof(char), 1, fin);   /* Sweepファイル名 */
	fread(&usdat, sizeof(unsigned short), 1, fin);   /* bit resolution */
	fread(&usdat, sizeof(unsigned short), 1, fin);   /* 波形種別 */
	fread(&usdat, sizeof(unsigned short), 1, fin);   /* 波形関数 */
	fread(&fdat, sizeof(float), 1, fin);   /* 振幅補正 */
	fread(&usdat, sizeof(unsigned short), 1, fin);   /* データ数 */
	fread(&fdat, sizeof(float), 1, fin);   /* 比帯域 */
	//fread(&fdat,sizeof(float),1,fin);   /* Flat比率 */
	//fread(&fdat,sizeof(float),1,fin);   /* 送信波形 台 */
	for (i = 0; i<8; ++i) fread(&cdat, sizeof(char), 1, fin);  /* スペア */
	fread(&fdat, sizeof(float), 1, fin);   /* 送信サンプル周波数 [MHz] */
	fread(&usdat, sizeof(unsigned short), 1, fin);   /* Ex-PHDのon/off */
	fread(&usdat, sizeof(unsigned short), 1, fin);   /* B FE-Gain上限 */
	fread(&usdat, sizeof(unsigned short), 1, fin);   /* B FE-Gainの下限 */
	fread(&usdat, sizeof(unsigned short), 1, fin);   /* B-Gain */
	fread(&usdat, sizeof(unsigned short), 1, fin);   /* 送信回数 */
	printf("number of transmissions per line: %d\n", usdat);
	fread(&usdat, sizeof(unsigned short), 1, fin);   /* アライメント */
	fread(&usdat, sizeof(unsigned short), 1, fin);   /* プレーン識別番号 */
	fread(&usdat, sizeof(unsigned short), 1, fin);   /* 空間コンパウンド識別番号 */
	printf("spatial compound: %d\n", usdat);
	fread(&fdat, sizeof(float), 1, fin);   /* 音響パワー */
	fread(&usdat, sizeof(unsigned short), 1, fin);   /* Interm Frame */
	fread(&fdat, sizeof(float), 1, fin);   /* Interm Inverval */
	fread(&usdat, sizeof(unsigned short), 1, fin);   /* エラスト識別番号 */
	fread(&usdat, sizeof(unsigned short), 1, fin);   /* RFビット位置 */
	for (i = 0; i<4; ++i) fread(&cdat, sizeof(char), 1, fin);
	//fread(&uldat,sizeof(unsigned long),1,fin);
	/*printf("number of data = %d bytes\n",uldat);*/
	/*printf("%d bytes\n",4*4+4*ncount*nposi*nframe);*/


	fseek(fin, 996L, SEEK_SET); //932

	// loading RF data
	tcnt = (4 + ncount)*nposi*nframe;
	raw_dat = (int *)malloc((size_t)(4 * tcnt));
	fread(raw_dat, sizeof(int), tcnt, fin);

	printf("raw size = %f\n", sizeof(raw_dat) / sizeof(int));


	// loading physio data (ECG)
	fread(&idcnt, sizeof(int), 1, fin);
	icnt = (idcnt & 0x3FFFC00) >> 6;
	printf("icnt = %d\n", icnt);
	for (i = 0; i<3; ++i) fread(&idcnt, sizeof(int), 1, fin);
	for (i = 0; i<icnt; ++i){
		fread(&emin, sizeof(unsigned short), 1, fin);
		fread(&emax, sizeof(unsigned short), 1, fin);
		fread(&pmin, sizeof(unsigned short), 1, fin);
		fread(&pmax, sizeof(unsigned short), 1, fin);
		fread(&plmin, sizeof(unsigned short), 1, fin);
		fread(&plmax, sizeof(unsigned short), 1, fin);
		for (j = 0; j<18; ++j) fread(&idcnt, sizeof(unsigned short), 1, fin);
		idcnt = PRF2*i / 1000.;
		ecg[i] = emax;
		pcg[i] = plmax;

	}

	fclose(fin);

	fin = fopen("./ecg.dat", "w");
	for (i = 0; i<icnt; ++i) fprintf(fin, "%f\n", ecg[i]);
	fclose(fin);

	/*fin=fopen("./pcg.dat","w");
	for(i=0;i<icnt;++i) fprintf(fin,"%f\n",pcg[i]);
	fclose(fin);*/

	for (i = 0; i<nframe; ++i){
		for (j = 0; j<nposi; ++j){
			ihof[j][i] = (4 + ncount)*j + (4 + ncount)*nposi*i;
		}
	}

	nfr_low = nframe;
	printf("interval of beams = %f mm\n", elp);
	dadd = 4 * 2 * elp / 0.15;

	/* Post ECG */
	for (i = 0; i<nframe; ++i){
		MpF = (int)(1000 / PRF2*i);
		ecgmax = ecg[MpF];
		ecg[i] = ecgmax;
	}

	fin = fopen("./ecg2.dat", "w");
	for (i = 0; i<nframe; ++i) fprintf(fin, "%f\n", ecg[i]);
	fclose(fin);
	/* *************************************************************************
	Setting Parameters (Temporal)
	********************************************************************** */
	if (aux_flg == 0) iaux = 5;
	if (aux_flg == 1) iaux = 6;
	f0i = f0;
	ifs = (int)(fs / 1e6);
	nspc = 1;
	if (SPAv == 1) nspc = 1;
	/* *************************************************************************
	Setting Parameters
	********************************************************************** */
	heartrate = 100;
	p_s = (float)ip_s;
	p_d = (float)ip_d;
	icgain = 4000000;
	idltx = 10;      /* nwv wavelength (default:4x) */
	idltx2 = 10;
	deltad = c0 / (fs*2.0);
	lp2s = (int)(2.0 / f0*fs);
	lp2q = (int)(2.0 / f0*fs); /* cut off frequency in quadrature demodulation */

	//f0=5.5*1e+6;
	/* *************************************************************************
	Memory Allocation
	********************************************************************** */
	if (ncount >= nframe){
		dmy = (float *)malloc((size_t)(8 * (ncount + 1)));
		sgn = (float *)malloc((size_t)(8 * (ncount + 1)));
	}
	if (nframe>ncount){
		dmy = (float *)malloc((size_t)(8 * (nframe + 1)));
		sgn = (float *)malloc((size_t)(8 * (nframe + 1)));
	}
	/* *************************************************************************
	Read ECG, PCG, AUX
	********************************************************************** */
	ecgbias = 0.0; icnt = 0;
	for (i = 0; i<nframe / 2; ++i){
		ecgbias = ecgbias + ecg[i];
		++icnt;
	}
	for (i = 0; i<nframe; ++i) ecg[i] = ecg[i] - ecgbias / (float)icnt;
	ecgmax = 0.0;
	for (i = 0; i<nframe / 2; ++i){
		if (fabs(ecg[i])>ecgmax) ecgmax = fabs(ecg[i]);
	}

	pcgbias = 0.0; icnt = 0;
	for (i = 0; i<nframe / 2; ++i){
		pcgbias = pcgbias + pcg[i];
		++icnt;
	}
	for (i = 0; i<nframe; ++i) pcg[i] = pcg[i] - pcgbias / (float)icnt;
	pcgmax = 0.0;
	for (i = 0; i<nframe / 2; ++i){
		if (fabs(pcg[i])>pcgmax) pcgmax = fabs(pcg[i]);
	}
	/* ****************************************************************************
	Detect Timing of R-wave
	************************************************************************* */
	nbeat = 0; beat = 0.0;
	ilevel = 0.6*ecgmax;
	period = 1.0 / (heartrate / 60.0);

	for (j = 0; j<maxfrm; ++j){
		if (j == 0) ii1 = 0;
		if (j >= 1) ii1 = track0[j - 1] + (int)(PRF2*period);
		if (ii1>nframe) break;
		for (i = ii1 + 1; i<nframe; ++i){
			if (ecg[i] >= ilevel){
				track0[j] = i;
				nbeat = nbeat + 1;
				break;
			}
			if (i == nframe - 1) break;
		}
		if (track0[j]<0){
			nbeat = nbeat - 1;
			break;
		}
		if (track0[j] == 0) break;
		for (i = 0; i <= j - 1; ++i){
			beat = beat + (float)(track0[i + 1] - track0[i]);
		}
		beat = beat / (float)(j + 1);
	}
	tbeat = beat / fs * 1000;

	/*for(i=0;i<nbeat;++i){
	rr=ecg[track0[i]];
	for(j>track0[i]-(int)(0.01*PRF2);j<track0[i]+(int)(0.01*PRF2);++j){
	if(j>=0 && j<nframe && fabs(ecg[j])>rr){
	rr=fabs(ecg[j]);
	icnt=j;
	}
	}
	track0[i]=j;
	}*/

	for (i = 0; i<nbeat; ++i){
		printf("track0(%d) = %d\n", i, track0[i]);
	}
	//nbeat=2;
	//track0[0]=0; track0[1]=nframe;

	/* *************************************************************************
	Setting Parameters
	********************************************************************** */
	ifrm1 = 0;
	ifrm2 = nframe;
	if (nbeat == 0){
		ifrm1 = 0;
		ifrm2 = nframe;
	}
	if (nbeat>0){
		ifrm1 = track0[0] - (int)(0.2*PRF2);
		if (ifrm1<0) ifrm1 = 0;
		ifrm2 = nframe;
	}
	printf("ifrm1, ifrm2 = %d, %d\n", ifrm1, ifrm2);
	if (ifrm1<0 || ifrm2>nframe || ifrm1>ifrm2) exit(0);
	/* *************************************************************************
	Display B-mode Image
	********************************************************************** */
	rmaxz = 0.0;
	for (i = 0;i < nposi;++i){
		for (j = 0; j < ncount; ++j) {
			dmy[j] = pow((float)raw_dat[ihof[i][track0[0]] + j + 4], 2);
		}
		wint2s(wary, lp2s);
		ndmy = ncount;
		lpfs(wary, dmy, sgn, lp2s, ndmy);
		for (j = 0; j<ncount; ++j){
			sgn[j] = sqrt(sgn[j]);
			if (fabs(sgn[j])>rmaxz) rmaxz = fabs(sgn[j]);
		}
	}

	rwb = dadd / 8.;
	nwb = (int)rwb;
	/* 直交検波 */
	for (i = 0; i<nposi; ++i){
		for (j = 0; j<ncount; ++j){
			dmy[j] = pow((float)raw_dat[ihof[i][track0[0]] + j + 1], 2);
		}
		wint2s(wary, lp2s);
		ndmy = ncount;
		lpfs(wary, dmy, sgn, lp2s, ndmy);
		for (j = 0; j<ncount; ++j) sgn[j] = sqrt(sgn[j]);

		for (j = 0; j<ncount; ++j){
			ic = 49 + (int)(100.*(20.*log10((sgn[j] + 1e-4) / rmaxz) + B_gain) / B_gain);
			if (ic<49) ic = 49;
			if (ic>148) ic = 148;
			ic = (ic - 49) * 255 / 99;
			gscolor = gscol256(ic, ic, ic);
			for (k = 0;k < (int)(izoom_B / (mpp / deltad) + 1);++k) {
				gsline(Bx_ini + (int)(izoom_B*rwb*i),
					By_ini + (int)(izoom_B*j*deltad / mpp) + k,
					Bx_ini + (int)(izoom_B*rwb*(i + 1)),
					By_ini + (int)(izoom_B*j*deltad / mpp) + k,
					gscolor);
			}
		}
	}

	gscolor = gscol256(0, 0, 0);
	gsline(Bx_ini - 15, By_ini, Bx_ini - 15,
		By_ini + (int)(izoom_B*ncount*deltad / mpp), gscolor);
	for (i = 0;i<ncount;i = i + (int)(5 * mpp / deltad)) {
		gsline(Bx_ini - 15, By_ini + (int)(izoom_B*i*deltad / mpp),
			Bx_ini - 10, By_ini + (int)(izoom_B*i*deltad / mpp), gscolor);
	}
	for (i = 0;i<ncount / 5;i = i + (int)(5 * mpp / deltad)) {
		gsline(Bx_ini - 15, By_ini + (int)(5.*izoom_B*i*deltad / mpp),
			Bx_ini - 5, By_ini + (int)(5.*izoom_B*i*deltad / mpp), gscolor);
	}
	gsline(Bx_ini, By_ini - 15, Bx_ini + (int)(izoom_B*rwb*nposi), By_ini - 15, gscolor);
	for (i = 0;i<(int)(rwb*nposi / 10.) + 1;++i) {
		gsline(Bx_ini + izoom_B * 10 * i, By_ini - 15,
			Bx_ini + izoom_B * 10 * i, By_ini - 5, gscolor);
	}
	sprintf(buf, "%4d %56s", iseq, fall);
	ptext(10, By_ini - 65, strlen(buf), buf);
	sprintf(buf, "age:%2d BP: %3d/%3d", iage, ip_s, ip_d);
	ptext(10, By_ini - 50, strlen(buf), buf);
	sprintf(buf, "mkdirhier ./images/%s", fname);
	//system(buf);
	/*sprintf(buf,"import -window 0x%x -crop %dx%d+%d+%d -silent ./images/%s/Bmode.gif",
	wid,(int)(izoom_B*rwb*nposi),
	(int)(izoom_B*ncount*deltad/mpp),
	Bx_ini,wwB-By_ini+60,fname);*/
	sprintf(buf, "import -window 0x%x -silent ./images/%s/Bmode.gif", wid, fname);
	//system(buf);
	/* *************************************************************************
	Setting Tracking Position by Tracing B-mode (itrace=1)
	********************************************************************** */
	if (itrace == 1){
		gscolor = gscol256(0,0,255);
		for (i = 0; i<nposi; ++i) ipeakflg[i] = 0;
		cursorpos(0, 0);
		trackingpoint_top();
		cursorpos(0, 0);
		for (i = 0; i<nposi0; ++i) ipeakflg[i] = 0;
		trackingpoint_bot();

		///* Save Peak Position */
		//fin = fopen("./tmp/peaks", "w");
		//for (i = 0; i<nposi; ++i)
		//	fprintf(fin, "%d %d %d %d\n", ipeak11[i], ipeak22[i], ipeak33[i], ipeak44[i]);
		//fclose(fin);
	}
	if (itrace == 0){
		sprintf(buf, "./tmp/peaks.%d", iseq1);
		fin = fopen(buf, "r");
		for (i = 0; i<nposi; ++i){
			fscanf(fin, "%d %d %d %d", &ip1, &ip2, &ip3, &ip4);
			ipeak11[i] = ip1; ipeak22[i] = ip2; ipeak33[i] = ip3; ipeak44[i] = ip4;
		}
		fclose(fin);
	}
	if (DispTr == 1){
		gscolor = gscol256(255,0,0);
		for (i = 0; i<nposi; ++i){
			gscolor = gscol256(0,0,0);
			for (k = 0; k<(int)(izoom_B / (mpp / deltad) + 1); ++k)
				gsline(Bx_ini + (int)(izoom_B*rwb*i),
				By_ini - (int)(izoom_B*ipeak33[i] * deltad / mpp) - k,
				Bx_ini + (int)(izoom_B*rwb*(i + 1)),
				By_ini - (int)(izoom_B*ipeak33[i] * deltad / mpp) - k, gscolor);
		}
	}

	for (i = 0; i<nposi; ++i){
		if (ipeak44[i] - ipeak33[i]<4 * idltx2) ipeak44[i] = ipeak33[i] + 4 * idltx2;
	}
	Sleep(1);
	gscolor = gscol256(0,0,0);
	if (ipause == 1){
		printf("Input 0 if ready to go next.\n");
		scanf("%d", &iflg);
	}
	/* *************************************************************************
	Clear Window
	********************************************************************** */
	eraseg();
	/* *************************************************************************
	Estimation of Displacements
	********************************************************************** */
	ifrm1 = track0[anabt]; ifrm2 = track0[anabt + 2]; dsh = 1;

	for (i = 0; i<nDFT; ++i) han_w[i] = 0.5 - 0.5*cos(rad*i / (float)nDFT);
	sprintf(buf, "mkdirhier ./images/%s/Mmode", fname);
	//system(buf);
	//sprintf(buf, "./ddmax_%s.dat", fname);
	sprintf(buf, "./ddmax.dat", fname);
	fin = fopen(buf, "w");
	for (iiposi = 0; iiposi<nposi; ++iiposi){
		ddmax = 0.0;
		eraseg();
		sprintf(buf, "%4d %56s", iseq, fall);
		ptext(100, wwB + 10, strlen(buf), buf);
		sprintf(buf, "age:%2d, BP: %3d/%3d, vel. max: %5.1f mm/s, dd max: %6.1f um, beam position: %2d",
			iage, ip_s, ip_d, 1e3*rryy, 1e6*rmaxdh, iiposi + 1);
		ptext(100, wwB + 22, strlen(buf), buf);
		if (aux_flg == 0)
			sprintf(buf, "max. & bias of ECG: %7.1f, %7.1f, max. PCG: %7.1f",
			ecgmax, ecgbias / PRF2, pcgmax);
		if (aux_flg == 1)
			sprintf(buf, "max & bias of ECG: %7.1f, %7.1f, max. PULSE: %7.1f",
			ecgmax, ecgbias / PRF2, pcgmax);
		ptext(100, wwB + 34, strlen(buf), buf);
		/* *************************************************************************
		Setting Tracking Position and Number of Tracking Points
		********************************************************************** */
		ndep[iiposi] = 2;
		/* *************************************************************************
		Display M-mode Image
		********************************************************************** */
		lp2s = (int)(2.0 / f0*fs);
		Mx_ini = 50; My_ini = 800;
		ry0 = nframe / 200.;
		iy0 = (int)ry0;
		if (iy0<1) iy0 = 1;
		for (i = 0; i<200; ++i){
			for (j = 0; j<ncount; ++j)
				dmy[j] = pow((float)raw_dat[ihof[iiposi][(int)(ry0*i)] + j], 2);
			wint2s(wary, lp2s);
			ndmy = ncount;
			lpfs(wary, dmy, sgn, lp2s, ndmy);
			for (j = 0; j<ncount; ++j)
				sgn[j] = sqrt(sgn[j]);
			for (j = 0; j<ncount; j = j + (int)(mpp / deltad / izoom_M)){
				ic = 49 + (int)(100.*((20.*log10(sgn[j] / rmaxz) + B_gain) / B_gain));
				if (ic<49) ic = 49;
				if (ic>148) ic = 148;
				ic = (ic - 49) * 255 / 99;
				gscolor = gscol256(ic, ic, ic);
				for (k = -1; k <= 1; ++k)
					gsline(Mx_ini + (int)(izoom_M*j*deltad / mpp), My_ini - 3 * i - k,
					Mx_ini + (int)(izoom_M*j*deltad / mpp), My_ini - 3 * (i + 1) - k, gscolor);
			}
		}
		gscolor = gscol256(0, 0, 0);
		gsline(Mx_ini, My_ini + 15,
			Mx_ini + (int)(izoom_M*ncount*deltad / mpp), My_ini + 15, gscolor);
		for (i = 0; i < ncount; i = i + (int)(5 * mpp / deltad)) {
			gsline(Mx_ini + (int)(izoom_M*i*deltad / mpp), My_ini + 15,
				Mx_ini + (int)(izoom_M*i*deltad / mpp), My_ini + 10, gscolor);
		}
		for (i = 0; i < ncount / 5; i = i + (int)(5 * mpp / deltad)) {
			gsline(Mx_ini + (int)(5.*izoom_M*i*deltad / mpp), My_ini + 15,
				Mx_ini + (int)(5.*izoom_M*i*deltad / mpp), My_ini + 5, gscolor);
		}
		gsline(Mx_ini + (int)(izoom_M*ncount*deltad / mpp) + 10, My_ini,
			Mx_ini + (int)(izoom_M*ncount*deltad / mpp) + 10, My_ini - 600, gscolor);
		for (i = 0; i < (int)(2.*nframe / PRF2) + 1; ++i) {
			gsline(Mx_ini + (int)(izoom_M*ncount*deltad / mpp) + 10,
				My_ini - (int)(300.*PRF2 / nframe*i),
				Mx_ini + (int)(izoom_M*ncount*deltad / mpp) + 5,
				My_ini - (int)(300.*PRF2 / nframe*i), gscolor);
		}
			
		/* *************************************************************************
		Frame for Displaying Signals
		********************************************************************** */
		for (i = 0; i<4; ++i){
			gsrect(Mx_ini + (int)(izoom_M*ncount*deltad / mpp) + 20 + 100 * i,
				My_ini + 15,
				Mx_ini + (int)(izoom_M*ncount*deltad / mpp) + 110 + 100 * i,
				My_ini - 615, RGB(255, 255, 255), RGB(0, 0, 0));
		}
		for (i = 0; i<4; ++i){
			for (j = 0; j<2; ++j){
				for (k = 0; k<3; ++k){
					gsline(Mx_ini + (int)(izoom_M*ncount*deltad / mpp) + 30 + 35 * k + 100 * i,
						My_ini + 15 - 622 * j,
						Mx_ini + (int)(izoom_M*ncount*deltad / mpp) + 30 + 35 * k + 100 * i,
						My_ini + 15 - 8 - 622 * j, gscolor);
				}
			}
		}
		for (i = 0; i<4; ++i){
			gsline_d(Mx_ini + (int)(izoom_M*ncount*deltad / mpp) + 20 + 45 + 100 * i, My_ini,
				Mx_ini + (int)(izoom_M*ncount*deltad / mpp) + 20 + 45 + 100 * i,
				My_ini - 600, gscolor);
		}

		for (i = 0; i<4; ++i){
			for (j = 0; j<2; ++j){
				for (k = 0; k<(int)(2.*(ifrm2 - ifrm1) / PRF2) + 1; ++k){
					gsline(Mx_ini + (int)(izoom_M*ncount*deltad / mpp) + 20 + 100 * i + 82 * j,
						My_ini - (int)(300.*PRF2 / (ifrm2 - ifrm1)*k),
						Mx_ini + (int)(izoom_M*ncount*deltad / mpp) + 20 + 8 + 100 * i + 82 * j,
						My_ini - (int)(300.*PRF2 / (ifrm2 - ifrm1)*k), gscolor);
				}
			}
		}
		/* *************************************************************************
		Display ECG
		********************************************************************** */
		ry0 = 600. / nframe; rx0 = 35. / ecgmax;
		for (i = 0; i<nframe - 2; ++i){
			nx0 = (int)(izoom_M*ncount*deltad / mpp) + 20 + 45;
			gsline(Mx_ini + nx0 - (int)(rx0*ecg[i]),
				My_ini - (int)(ry0*i),
				Mx_ini + nx0 - (int)(rx0*ecg[i + 1]),
				My_ini - (int)(ry0*(i + 1)), gscolor);
		}
		/* *************************************************************************
		Display Timing of R-wave
		********************************************************************** */
		if (nbeat != 0){
			nbeat2 = nbeat;
			for (i = 0; i<nbeat; ++i){
				if (track0[i]>0){
					nbeat1 = i;
					break;
				}
			}
			for (i = 0; i<nbeat; ++i){
				if (track0[i]>nframe){
					nbeat2 = i - 1;
					break;
				}
			}
			if (nbeat1 == 0) nbeat1 = 0;
			if (nbeat2 == 0) nbeat2 = 0;
			gscolor = gscol256(255, 0, 0);
			nx0 = (int)(izoom_M*ncount*deltad / mpp);
			for (i = nbeat1; i<nbeat2; ++i){
				gsline_d(Mx_ini, My_ini - (int)(ry0*track0[i]),
					Mx_ini + nx0, My_ini - (int)(ry0*track0[i]), gscolor);
			}
			nx0 = (int)(izoom_M*ncount*deltad / mpp) + 30;
			for (i = 0; i<4; ++i){
				for (j = nbeat1; j<nbeat2; ++j){
					gsline_d(Mx_ini + nx0 + 100 * i, My_ini - (int)(ry0*track0[j]),
						Mx_ini + nx0 + 70 + 100 * i, My_ini - (int)(ry0*track0[j]), gscolor);
				}
			}
			gscolor = gscol256(0, 0, 0);
		}

		/* *************************************************************************
		Display PCG
		********************************************************************** */
		/*ry0=600./nframe; rx0=35./pcgmax;
		for(i=0;i<nframe-1;++i){
		nx0=(int)(izoom_M*ncount*deltad/mpp)+20+145;
		gsline(Mx_ini+nx0-(int)(rx0*pcg[i]),
		My_ini+(int)(ry0*i),
		Mx_ini+nx0-(int)(rx0*pcg[i+1]),
		My_ini+(int)(ry0*(i+1)));
		}*/

		/* *************************************************************************
		Start of Loop for Frame for Displacement Estimation
		********************************************************************** */
		for (idcnt = 0; idcnt<ndep[iiposi]; ++idcnt){
			disp[idcnt][ifrm1] = 0.0;
			dd[ifrm1] = 0.0;
		}
		ip1 = 0; ip2 = ncount;


		for (ifrm = ifrm1; ifrm<ifrm2; ++ifrm){

			/* Quadrature Demodulation */
			czeros = cmplx(0.0, 0.0);
			for (i = -1; i<2; ++i){
				if (iiposi + i>-1 && iiposi + i<nposi){
					/* For Previous Frame */
					if (ifrm == ifrm1){
						for (j = ip1; j<ip2; ++j){
							qdr_p[j][i + 1] = 0.0; qdi_p[j][i + 1] = 0.0; icnt = 0;
							for (k = -nDFT / 2; k < nDFT / 2; ++k) {
								if (j + k >= 0 && j + k<ncount){
									qdr_p[j][i + 1] = qdr_p[j][i + 1] + han_w[k + nDFT / 2] *
										cos(nfrq*rad*(k + nDFT / 2) / (float)nDFT)*
										(float)raw_dat[ihof[iiposi + i][ifrm] + j + k];
									qdi_p[j][i + 1] = qdi_p[j][i + 1] - han_w[k + nDFT / 2] *
										sin(nfrq*rad*(k + nDFT / 2) / (float)nDFT)*
										(float)raw_dat[ihof[iiposi + i][ifrm] + j + k];
									++icnt;
								}
							}
							qdr_p[j][i + 1] = qdr_p[j][i + 1] / (float)icnt;
							qdi_p[j][i + 1] = qdi_p[j][i + 1] / (float)icnt;
						}
					}
					/* For Post-Frame */
					for (j = ip1; j<ip2; ++j){
						qdr_f[j][i + 1] = 0.0; qdi_f[j][i + 1] = 0.0; icnt = 0;
						for (k = -nDFT / 2; k<nDFT / 2; ++k){
							if (j + k >= 0 && j + k<ncount){
								qdr_f[j][i + 1] = qdr_f[j][i + 1] + han_w[k + nDFT / 2] *
									cos(nfrq*rad*(k + nDFT / 2) / (float)nDFT)*
									(float)raw_dat[ihof[iiposi + i][ifrm + 1] + j + k];
								qdi_f[j][i + 1] = qdi_f[j][i + 1] - han_w[k + nDFT / 2] *
									sin(nfrq*rad*(k + nDFT / 2) / (float)nDFT)*
									(float)raw_dat[ihof[iiposi + i][ifrm + 1] + j + k];
								++icnt;
							}
						}
						qdr_f[j][i + 1] = qdr_f[j][i + 1] / (float)icnt;
						qdi_f[j][i + 1] = qdi_f[j][i + 1] / (float)icnt;
					}
				}
			}
			/* Start of Loop of Depth for Displacement Estimation */
			for (idcnt = 0; idcnt<ndep[iiposi]; ++idcnt){
				if (ifrm == ifrm1){
					if (idcnt == 0) peak0[idcnt] = ipeak11[iiposi];
					if (idcnt == 1) peak0[idcnt] = ipeak22[iiposi];
				}
				xpeak0[idcnt] = (float)peak0[idcnt];
				// velocity estimation
				crr0 = czeros; crr1 = czeros;
				for (ispc = -nspc; ispc<nspc + 1; ++ispc){
					if (iiposi + ispc >= 0 && iiposi + ispc<nposi){
						for (i = 0; i <= 2 * idltx; ++i){
							if (i + peak0[idcnt] >= 0 && i + peak0[idcnt]<ncount){
								rr = qdr_p[i + peak0[idcnt]][ispc + 1];
								ri = qdi_p[i + peak0[idcnt]][ispc + 1];
								cdmyp = cmplx(rr, ri);
								rr = qdr_f[i + peak0[idcnt]][ispc + 1];
								ri = qdi_f[i + peak0[idcnt]][ispc + 1];
								cdmy = cmplx(rr, ri);
								ctemp = conjg(cdmyp);
								crr0 = cadd(crr0, cmul(cdmy, ctemp));
							}
						}
					}
				}
				vel_p[idcnt] = vel[idcnt];
				if (ifrm == ifrm1) vel_p[idcnt] = 0.0;
				vel[idcnt] = -0.5*c0*atan2(crr0.i, crr0.r) / rad / f0;

				disp[idcnt][ifrm + 1] = disp[idcnt][ifrm] + vel[idcnt];
				if (idcnt == 0){
					xpeak0[idcnt] = ipeak11[iiposi] + disp[idcnt][ifrm + 1] / deltad;
				}
				if (idcnt == 1){
					xpeak0[idcnt] = ipeak22[iiposi] + disp[idcnt][ifrm + 1] / deltad;
				}
				peak0[idcnt] = (int)xpeak0[idcnt];
			}
			dd[ifrm + 1] = disp[1][ifrm + 1] - disp[0][ifrm + 1];

			ry0 = 600. / nframe;
			for (idcnt = 0; idcnt<ndep[iiposi]; ++idcnt){
				/* Display Tracked Position */
				gscolor = gscol256(255, 0, 0);
				if (ifrm >= ifrm1 && ifrm<ifrm2)
					gsline(Mx_ini + (int)(izoom_M*peak0[idcnt] * deltad / mpp),
					My_ini - (int)(ry0*ifrm),
					Mx_ini + (int)(izoom_M*peak0[idcnt] * deltad / mpp),
					My_ini - (int)(ry0*(ifrm + 1)), gscolor);
				gscolor = gscol256(255, 255, 255);
				gsline(Mx_ini + 1 + (int)(izoom_M*peak0[idcnt] * deltad / mpp),
					My_ini - (int)(ry0*ifrm),
					Mx_ini + 1 + (int)(izoom_M*peak0[idcnt] * deltad / mpp),
					My_ini - (int)(ry0*(ifrm + 1)), gscolor);
				gscolor = gscol256(0, 0, 0);
				/* Display Velocity */
				rx0 = 35. / rryy;
				if (idcnt == 0) nx0 = (int)(izoom_M*ncount*deltad / mpp) + 20 + 145;
				if (idcnt == 1) nx0 = (int)(izoom_M*ncount*deltad / mpp) + 20 + 245;

				if (ifrm > ifrm1 && ifrm < ifrm2 - 1) {
					gsline(Mx_ini + nx0 + (int)(rx0*vel_p[idcnt] * PRF2),
						My_ini - (int)(ry0*ifrm),
						Mx_ini + nx0 + (int)(rx0*vel[idcnt] * PRF2),
						My_ini - (int)(ry0*(ifrm + 1)), gscolor);
				}
			}

			/* Display Change in Diameter */
			rx0 = 35. / rmaxdh;
			nx0 = (int)(izoom_M*ncount*deltad / mpp) + 20 + 345;
			if (ifrm >= ifrm1 && ifrm<ifrm2 - 1){
				gsline(Mx_ini + nx0 - (int)(rx0*dd[ifrm]),
					My_ini - (int)(ry0*ifrm),
					Mx_ini + nx0 - (int)(rx0*dd[ifrm + 1]),
					My_ini - (int)(ry0*(ifrm + 1)), gscolor);
			}


			for (i = 0; i<3; ++i){
				for (j = 0; j<ncount; ++j){
					qdr_p[j][i] = qdr_f[j][i]; qdi_p[j][i] = qdi_f[j][i];
				}
			}
			if (ddmax>dd[ifrm + 1])
				ddmax = dd[ifrm + 1];

		} // end of loop of frame
		sprintf(buf, "import -window 0x%x -silent ./images/%s/Mmode/pos%d.gif", wid, fname, iiposi);
		//system(buf);

		fprintf(fin, "%d %f\n", iiposi, ddmax*1e+6);

	} // end of loop of beam position
	fclose(fin);
	/* *************************************************************************
	Terminate Program
	********************************************************************** */
	free(raw_dat); free(sgn); free(dmy);
	end();
	printf("INPUT 0 for termination of program.\n");
	scanf("%d", &iflg);
	printf("*******************************\n");
	printf("*     normally terminated     *\n");
	printf("*******************************\n");
}

// サブ関数
int trackingpoint_top() {

	MSG msg;

	while (GetMessage(&msg, NULL, 0, 0)>0) {
		DispatchMessage(&msg);
		ibx = x_coordinate;
		iby = y_coordinate;
		if (iby <= By_ini + (int)(izoom_B*ncount*deltad / mpp) && iby >= By_ini) {
			for (i = 0; i<nposi; ++i) {
				nx0 = Bx_ini + (int)(izoom_B*rwb*i);
				if (ibx >= nx0 && ibx<nx0 + (int)(izoom_B*rwb) && ipeakflg[i] == 0) {
					ipeak11[i] = (int)((iby - By_ini)*mpp / deltad / izoom_B) - (int)(mpp / deltad);
					gsline(nx0, iby, nx0 + (int)(izoom_B*rwb), iby, gscolor);
					ipeakflg[i] = 1;
				}
			}
		}
		ipeaksum = 0;
		for (i = 0; i<nposi; ++i) ipeaksum = ipeaksum + ipeakflg[i];
		if (ipeaksum == nposi) break;
	}
	return 0;
}

int trackingpoint_bot() {

	MSG msg;

	while (GetMessage(&msg, NULL, 0, 0)>0) {
		DispatchMessage(&msg);
		ibx = x_coordinate;
		iby = y_coordinate;
		if (iby <= By_ini + (int)(izoom_B*ncount*deltad / mpp) && iby >= By_ini) {
			for (i = 0; i<nposi; ++i) {
				nx0 = Bx_ini + (int)(izoom_B*rwb*i);
				if (ibx >= nx0 && ibx<nx0 + (int)(izoom_B*rwb) && ipeakflg[i] == 0) {
					ipeak22[i] = (int)((iby - By_ini)*mpp / deltad / izoom_B) - (int)(mpp / deltad);
					gsline(nx0, iby, nx0 + (int)(izoom_B*rwb), iby, gscolor);
					ipeakflg[i] = 1;
				}
			}
		}
		ipeaksum = 0;
		for (i = 0; i<nposi; ++i) ipeaksum = ipeaksum + ipeakflg[i];
		if (ipeaksum == nposi) break;
	}
	return 0;
}