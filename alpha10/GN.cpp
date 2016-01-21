#include "stdafx.h"
#include "header.h"

using namespace std;

GN::GN(const vector<double> &x_ini, const vector<double> &y_ini)
{
	m = x_ini.size();
	x = ippsMalloc_32f(m);
	y = ippsMalloc_32f(m);
	for (int i = 0; i < m; ++i){
		x[i] = x_ini[i];
		y[i] = y_ini[i];
	}
}

vector<double> GN::solve(const double &r_ini, const double &c_ini)
{
	double p = 0.2e3;
	double l = 9.5e3;
	double sn = 0.15646446504; //sin(9deg)
	double S = 0.0;

	double r = r_ini;
	double c = c_ini;
	res = vector<double>(m, 0);

	vector<double> derr(m, 0);
	vector<double> derc(m, 0);
	double Sderr = 0.0;
	double Sderc = 0.0;
	double det = 0.0;
	double crosssum = 0.0;
	double tmpr = 0.0;
	double tmpc = 0.0;
	int cnt = 0;

	for (int i = 0; i < m; ++i){
		res[i] = y[i] - (r + sqrt(p*p*x[i] * x[i]+2*(r*sn-l)*p*x[i]+r*r+l*l-2*r*l*sn))/c;
		S += res[i] * res[i];

		derr[i] = -(1 + (r + (p*x[i] - l)*sn) / sqrt(p*p*x[i] * x[i] + 2 * (r*sn - l)*p*x[i] + r*r + l*l - 2 * r*l*sn)) / c;
		Sderr += derr[i] * derr[i];

		derc[i] = (r + sqrt(p*p*x[i] * x[i] + 2 * (r*sn - l)*p*x[i] + r*r + l*l - 2 * r*l*sn)) / (c * c);
		Sderc += derc[i] * derc[i];

		crosssum += derr[i] * derc[i];
	}
	det = Sderr * Sderc - crosssum * crosssum;

	while (S > 1e-6){
		++cnt;
		tmpr = 0.0;
		tmpc = 0.0;
		S = 0.0;


		for (int i = 0; i < m; ++i){
			tmpr += res[i]*(derr[i] * Sderc - derc[i] * crosssum);
			tmpc += res[i]*(-derr[i] * crosssum + derc[i] * Sderr);
		}
		r -= tmpr / det*0.1;
		c -= tmpc / det*0.1;

		cout << "(r,c)=(" << r << "," << c << ")\n";

		Sderr = 0.0;
		Sderc = 0.0;
		crosssum = 0.0;

		for (int i = 0; i < m; ++i){
			res[i] = y[i] - (r + sqrt(p*p*x[i] * x[i] + 2 * (r*sn - l)*p*x[i] + r*r + l*l - 2 * r*l*sn)) / c;
			S += res[i] * res[i];

			derr[i] = -(1 + (r + (p*x[i] - l)*sn) / sqrt(p*p*x[i] * x[i] + 2 * (r*sn - l)*p*x[i] + r*r + l*l - 2 * r*l*sn)) / c;
			Sderr += derr[i] * derr[i];

			derc[i] = (r + sqrt(p*p*x[i] * x[i] + 2 * (r*sn - l)*p*x[i] + r*r + l*l - 2 * r*l*sn)) / (c * c);
			Sderc += derc[i] * derc[i];

			crosssum += derr[i] * derc[i];
		}
		det = Sderr * Sderc - crosssum * crosssum;
	}

	cout << "iteration: " << cnt << "\n";
	vector<double> ans = { r, c };
	cout << r << " " << c << "\n";
	return ans;
}