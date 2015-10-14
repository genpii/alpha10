#pragma once

#include <vector>

class RF {
	int in_i;
	float in_f;
	double in_d;
	unsigned short in_s;
	char *in_c;
public:
	void fileopen();
	vector<int> raw;
	

};