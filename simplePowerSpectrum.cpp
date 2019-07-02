#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <cmath>

using namespace std;
int const BINSIZE = 512;

int main(int argc, char *argv[]) {
	if(argc < 3) {
		cout << "Error: Need Input and Output Filenames\n";
		getchar();
	}
	else {
		string inputName(argv[1]);
		string outputName(argv[2]);
		cout << "File To Read: " << inputName << "\n";
		cout << "File To Write: " << outputName << "\n";
		getchar();
		ifstream inputFile (inputName);
		if(inputFile.is_open()) {
			//cout << "Able to read " << inputName << "\n";
			//getchar();
			int realGrid [3][3];
			int subPos1, subPos2;
			int count = 0;
			double x, y, z;
			string line, sub;
			while(getline(inputFile,line) && count < 3) {
				x = 0;
				y = 0;
				z = 0;
				subPos1 = 0;
				subPos2 = line.find(" ", 0);
				if(subPos2 >= 0) {
					sub = line.substr(subPos1, subPos2 + 1);
					x = stod(sub);
				}
				subPos1 = subPos2 + 1;
				subPos2 = line.find(" ", subPos1);
				if(subPos2 >= 0) {
					sub = line.substr(subPos1, subPos2 + 1);
					y = stod(sub);
				}
				subPos1 = subPos2 + 1;
				subPos2 = line.find(" ", subPos1);
				if(subPos2 >= 0) {
					sub = line.substr(subPos1, subPos2 + 1);
					z = stod(sub);
				}
				cout << x << "," << y << "," << z << "\n";
				count++;
			}
			inputFile.close();
		}
		else {
			cout << "Unable to read " << inputName << "\n";
			getchar();
		}
	}
	return 0;
}
