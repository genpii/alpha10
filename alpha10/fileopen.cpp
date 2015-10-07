// fileopen.cpp

#include "stdafx.h"
#include <fstream>
#include "fileopen.h"

using namespace std;

int fileopen()
{
	ifstream ifs("./RF/sample20151005_2.crf", ios_base::binary);
	if (!ifs){
		cout << "could't open\n";
		return 1;
	}
	
	
}