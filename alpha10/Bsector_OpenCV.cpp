// Bsector.cpp

#include "stdafx.h"
#include <opencv2/opencv.hpp>
#include <opencv2/opencv_lib.hpp>

using namespace std;
using namespace cv;

//const vector<vector<float>>& env, float dangle
void Bsector(const vector<vector<float>>& env, float dangle)
{
	////Mat src2 = Mat::zeros(320, 320, CV_16SC3);
	//Mat src2 = imread("./itsuki2.jpg");
	//resize(src2, src2, Size(0, 0), 0.1, 0.1);
	//Mat dst2;
	//cvtColor(src2, dst2, CV_BGR2GRAY);
	//threshold(dst2, src2, 120, 255, THRESH_BINARY);
	//imshow("sample", src2);
	//imwrite("dst.jpg", src2);
	//waitKey(0);
	
	// calculate max value for normalization
	int line = static_cast<int>(distance(env.begin(), env.end()));
	int sample = static_cast<int>(distance(env[0].begin(), env[0].end()));
	vector<float> maxvec(line, 0);
	vector<int>::iterator it;
	int tmp;
	for (int i = 0; i < line; ++i){
		maxvec[i] = *max_element(env[i].begin(), env[i].end());
	}
	float max = *max_element(maxvec.begin(), maxvec.end());

	//draw set
	Mat src = Mat(sample / 4 + 100, sample / 2, CV_8UC1, Scalar(255, 255, 255));
	Mat dst = Mat::zeros(500, 800, CV_8UC1);
	Point cen(sample / 4, 10);
	Size2d size;
	double angle = 0.0;
	double sangle = 90.0 - static_cast<double>((line * dangle) / 2.0);
	//float dangle = 2.0;
	double eangle;
	Scalar color;
	int c;

	double gain = 60.0;

	//int line = 45;
	//sample = 399;
	

	/*vector<vector<float>> env2(line, vector<float>(sample, 0));
	for (int i = 0; i < line; ++i){
		for (int j = 0; j < sample; ++j){
			if (i != 0 && i != line - 1){
				env2[i][j] = 0.5 * env[i][j] + 0.25 * (env[i - 1][j] + env[i + 1][j]);
			}
			else env2[i][j] = env[i][j];
		}
	}*/
	//draw
	//for (int i = 0; i < line; ++i){
	//	for (int j = 0; j < sample; ++j){
	//		//size.width = 400 - 400 * (j / sample);
	//		//size.height = 400 - 400 * (j / sample);
	//		//Size size2(400 * (1 - (j / sample)), 400 * (1 - (j / sample)));
	//		//c = static_cast<int>(100 * (log10(env[i][j] / max) + gain) / gain);
	//		//if (c > 255) c = 255;
	//		//if (c < 0) c = 0;

	//		c = j;
	//		color = (c,c,c);
	//		ellipse(src, cen, Size(500 * (1 - (j / sample)), 500 * (1 - (j / sample))), angle, sangle, eangle, color, 1, 8, 0);
	//	}
	//	angle += eangle;
	//}
	//size = Size(400, 400);
	for (int i = 0; i < line; ++i){
		eangle = sangle + static_cast<double>(dangle);
		for (int j = 0; j < sample / 4; ++j){
			size = Size2d(sample / 4 - j, sample / 4 - j);

			c = static_cast<int>(256 * (20 * log10(env[i][sample - j * 4] / max) + gain) / gain);
			if (c > 255) c = 255;
			if (c < 0) c = 0;
			color = Scalar(c, c, c);
			ellipse(src, cen, size, angle, sangle, eangle, color, 1, 4);
		}
		sangle = eangle;
	}
	imshow("B-mode", src);
	imwrite("./ca.png", src);
	waitKey(0);
	
}
