#include "stdafx.h"
#include "graphics.h"


using namespace std;


void Bsector3(const vector<vector<float>>& env, float dangle)
{
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
	int wwB = 800;
	unsigned long wid = pageb(0, 0, 1400, wwB, "B-mode");
	start();
	int Bx_ini = 500;
	int By_ini = 50;

	double angle = 90.0 - line * dangle / 2;
	double sangle = 0.0;
	//float dangle = 2.0;
	double eangle = static_cast<double>(dangle);
	double gain = 50.0;
	int thick = 2;
	COLORREF color;

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
	float r, theta1, theta2;
	int ic;
	for (int i = 0; i < line; ++i){
		theta1 = dangle * (line - 1) + dangle * i;
		theta2 = theta1 + dangle;

		for (int j = 0; j < sample / 4; ++j){
			ic = static_cast<int>(256 * (10.*log10(env[i][j] / max) + gain) / gain);
			if (ic < 0) ic = 0;
			if (ic > 255) ic = 255;
			color = gscol256(ic, ic, ic);
			
			gsarc(Bx_ini, By_ini, j, theta1, theta2, thick, color);
		}

	}
	end();
	char aaa;
	cin >> aaa;

}
