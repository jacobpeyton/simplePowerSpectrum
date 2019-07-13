#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <cmath>
#include <vector>
#include <vector_types.h>
#include <limits>

using namespace std;

int const TIMESRAN = 10;
double const ALPHA = 0.1;

vector<double> fftFrequencies(int N, double L) {
    vector<double> k(N);
    double dk = 2.0*M_PI/L;
    for (int i = 0; i <= N/2; ++i)
        k[i] = i*dk;
    for (int i = N/2 + 1; i < N; ++i)
        k[i] = (i - N)*dk;
    return k;
}

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
		cout << "Press Key to Continue\n";
		getchar();
		ifstream inputFile (inputName);
		if(inputFile.is_open()) {
			//cout << "Able to read " << inputName << "\n";
			//getchar();
			//Bin Galaxies to 512x512x512 bin grid
			int4 N = {512, 512, 512, 0};
			N.w = N.x*N.y*N.z;
			double3 L = {1024.0, 1024.0, 1024.0};
			double3 Delta_r = {L.x/N.x, L.y/N.y, L.z/N.z};
			vector<int> realGrid(N.w);

			int subPos1, subPos2;
			int count = 0;
			int x, y, z, index;
			string line, sub;
			while(getline(inputFile,line)) {
				x = 0;
				y = 0;
				z = 0;
				subPos1 = 0;
				subPos2 = line.find(" ", 0);
				if(subPos2 >= 0) {
					sub = line.substr(subPos1, subPos2 + 1);
					x = int (stod(sub)/Delta_r.x);
;				}
				subPos1 = subPos2 + 1;
				subPos2 = line.find(" ", subPos1);
				if(subPos2 >= 0) {
					sub = line.substr(subPos1, subPos2 + 1);
					y = int (stod(sub)/Delta_r.y);
				}
				subPos1 = subPos2 + 1;
				subPos2 = line.find(" ", subPos1);
				if(subPos2 >= 0) {
					sub = line.substr(subPos1, subPos2 + 1);
					z = int (stod(sub)/Delta_r.z);
				}
				
        			index = z + N.z*(y + N.y*x);
				if (index < N.w) {
				    realGrid[index] += 1.0;
				}
				++count;
			}
			inputFile.close();

			//Generate Random Points in second grid
			vector<int> randGrid(N.w);
			int ranCount = 0;
			double nbar = 1.0 * count / (L.x * L.y * L.z);
			double n = TIMESRAN * nbar * Delta_r.x * Delta_r.y * Delta_r.z;

			random_device seeder; // System entropy randoms used to seed pseudorandom generator
			mt19937_64 gen(seeder()); // Randomly seeded Mersenne twister
			poisson_distribution<int> pDist(n); // Poisson distribution with mean set to average number of points per cell

			for (int i = 0; i < N.x; ++i) {
				for(int j = 0; j < N.y; ++j) {
					for(int k = 0 k < N.z; ++k) {
						int numRans = pDist(gen); // Poisson sample to get number of randoms for the cell
						
        					index = z + N.z*(y + N.y*x);
						if (index < N.w) {
						    randGrid[index] += numRans;
						}
						ranCount += numRans;
					}
				}
			}

			//Take weighted difference
			vector<double> diffGrid(N.w);

			for (int i = 0; i < N.x; ++i) {
				for(int j = 0; j < N.y; ++j) {
					for(int k = 0 k < N.z; ++k) {
        					index = z + N.z*(y + N.y*x);
						if (index < N.w) {
						    diffGrid[index] = realGrid[index] - alpha * randGrid[index];
						}
					}
				}
			}

			//Feed 1D Representation of grid through fourier tranform, complex array is generated
			vector<fftw_complex> dk(N.x*N.y*(N.z/2 + 1));

			fftw_init_threads();
			fftw_import_wisdom_from_filename("fftwWisdom.dat");
			fftw_plan_with_nthreads(omp_get_max_threads());
			fftw_plan dr2dk = fftw_plan_dft_r2c_3d(N.x, N.y, N.z, diffGrid.data(), dk.data(), FFTW_MEASURE);
			fftw_export_wisdom_to_filename("fftwWisdom.dat");

			fftw_execute(dr2dk);

			//Calculate Frequences for each bin, and bin frequencies
			//Calculate power spectrum from frequencies
		}
		else {
			cout << "Unable to read " << inputName << "\n";
			getchar();
		}
	}
	return 0;
}
