#include "stdafx.h"
#include "header.h"

using namespace std;

GN::GN(const vector<float> &x_ini, const vector<float> &y_ini)
{
	m = x_ini.size();
	x = ippsMalloc_32f(m);
	y = ippsMalloc_32f(m);
	for (int i = 0; i < m; ++i){
		x[i] = x_ini[i];
		y[i] = y_ini[i];
	}
}

void GN::solve(const float &r_ini, const float &c_ini)
{
	
}