// physio.cpp
// extract physio data from DICOM

#include "stdafx.h"
#include "header.h"

using namespace std;
file::file(string filename)
{
	cout << "open file=" << filename << "\n";
	fin.open(filename, ios_base::in | ios_base::binary);
	if (!fin){
		cout << "couldn't load file.\n";
	}
}

file::~file()
{
	cout << "delete object.\n";
}

void file::open(string filename)
{
	fin.open(filename, ios_base::in | ios_base::binary);
	if (!fin){
		cout << "couldn't load file.\n";
	}
}

void file::start()
{
	fin.clear();
	fin.seekg(0, ios_base::beg);
}

void file::warp(int pos)
{
	fin.seekg(pos, ios_base::beg);
}

void file::go(int pos)
{
	fin.seekg(pos, ios_base::cur);
}

physio::physio(string filename) : file(filename)
{
}

physio::~physio()
{
}

int physio::extract(int offset)
{
	vector<unsigned char> check;
	char tmp;
	while (!fin.eof()){
		fin.read((char*)&tmp, sizeof(unsigned char));
		check.push_back(tmp);
	}
	fin.close();

	/*searching key*/
	list<unsigned char> key = { 0xFF, 0x53, 0x21, 0x10, 0x4F, 0x42, 0x00, 0x00 };
	vector<unsigned char>::iterator it = find_end(check.begin(), check.end(), key.begin(), key.end());
	int dist = distance(check.begin(), it);
	if (it == check.end()){
		cout << "not found key\n";
		return 1;
	}
	else{
		cout << "key is at check[" << dist << "]\n";
	}
	int size = (check.size() - dist - 12 - 1) / 32;
	vector<unsigned char>().swap(check);

	/*extract physio data*/
	open(fn);
	fin.seekg(dist + 12, ios_base::beg);

	short tmp2 = 0;

	ECG.reserve(size);
	PCG_min.reserve(size);
	PCG_max.reserve(size);
	for (int i = 0; i < size; ++i){
		fin.seekg(8, ios_base::cur);
		fin.read((char*)&tmp2, sizeof(short));
		ECG.push_back(tmp2 & 0x03FF);
		fin.read((char*)&tmp2, sizeof(short));
		PCG_min.push_back(tmp2 & 0x03FF);
		fin.read((char*)&tmp2, sizeof(short));
		PCG_max.push_back(tmp2 & 0x03FF);
		fin.seekg(18, ios_base::cur);
	}

	ECG.erase(ECG.begin(), ECG.begin() + offset);
	PCG_min.erase(PCG_min.begin(), PCG_min.begin() + offset);
	PCG_max.erase(PCG_max.begin(), PCG_max.begin() + offset);

	return 0;
}



int physio()
{
	extern vector<short> ECG;
	extern vector<short> PCG_min;
	extern vector<short> PCG_max;

	string fn = "D:/RFdata/study/20160104/2/2";
	ifstream fin(fn, ios_base::in | ios_base::binary);
	if (!fin){
		cout << "couldn't read physio file!\n";
		return 1;
	}
	cout << "success to open '" << fn << "'\n";
	
	/*load data for check*/
	cout << "loading physio data...\n";
	vector<unsigned char> check;
	char tmp;
	while (!fin.eof()){
		fin.read((char*)&tmp, sizeof(unsigned char));
		check.push_back(tmp);
	}
	fin.close();

	/*searching key*/
	list<unsigned char> key = {0xFF, 0x53, 0x21, 0x10, 0x4F, 0x42, 0x00, 0x00};
	vector<unsigned char>::iterator it = find_end(check.begin(), check.end(), key.begin(), key.end());
	int dist = distance(check.begin(), it);
	if (it == check.end()){
		cout << "not found key\n";
		return 1;
	}
	else{
		cout << "key is at check[" << dist << "]\n";
	}
	int size = (check.size() - dist - 12 - 1) / 32;
	vector<unsigned char>().swap(check);

	/*extract physio data*/
	fin.open(fn, ios_base::in | ios_base::binary);
	fin.seekg(dist + 12, ios_base::beg);

	short tmp2 = 0;
	
	ECG.reserve(size);
	PCG_min.reserve(size);
	PCG_max.reserve(size);
	for (int i = 0; i < size; ++i){
		fin.seekg(8, ios_base::cur);
		fin.read((char*)&tmp2, sizeof(short));
		ECG.push_back(tmp2 & 0x03FF);
		fin.read((char*)&tmp2, sizeof(short));
		PCG_min.push_back(tmp2 & 0x03FF);
		fin.read((char*)&tmp2, sizeof(short));
		PCG_max.push_back(tmp2 & 0x03FF);
		fin.seekg(18, ios_base::cur);
	}

	return 0;
}

