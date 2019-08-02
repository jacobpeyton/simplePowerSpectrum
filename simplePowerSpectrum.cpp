#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <cmath>
#include <vector>
#include <vector_types.h>
#include <limits>
#include <fftw3.h>
#include <omp.h>

int const TIMESRAN = 10;

double const k_min = 0.01;
double const k_max = 0.2;
int const N_k = 19;

void binPoints(std::ifstream &readFile, std::vector<double> &grid, int4 N, double3 L, int &count) {
	double3 Delta_r = {L.x/N.x, L.y/N.y, L.z/N.z};
	int subPos1, subPos2;
	std::string line, sub;
	while(getline(readFile,line)) {
		int x = 0;
		int y = 0;
		int z = 0;
		subPos1 = 0;
		subPos2 = line.find(" ", 0);
		if(subPos2 >= 0) {
			sub = line.substr(subPos1, subPos2 + 1);
			x = int (stod(sub)/Delta_r.x);
		}
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

		int index = z + N.z*(y + N.y*x);
		if (index < N.w) {
			grid[index] += 1.0;
		}
		++count;
	}
	readFile.close();
}

void generateAndBinPoints(std::ofstream &writeFile, std::vector<double> &grid, int4 N, double3 L, int &count, double n) {
	double3 Delta_r = {L.x/N.x, L.y/N.y, L.z/N.z};
	std::random_device seeder; // System entropy randoms used to seed pseudorandom generator
	std::mt19937_64 gen(seeder()); // Randomly seeded Mersenne twister
	std::poisson_distribution<int> pDist(n); // Poisson distribution with mean set to average number of points per cell
	std::uniform_real_distribution<double> uDist(0.0, 1.0); // Uniform distribution to place points within cell

	for (int i = 0; i < N.x; ++i) {
		double r_x = i*Delta_r.x; // Cell offset in x direction
		for(int j = 0; j < N.y; ++j) {
			double r_y = j*Delta_r.y; // Cell offset in y direction
			for(int k = 0; k < N.z; ++k) {
				double r_z = k*Delta_r.z; // Cell offset in z direction
				int numRans = pDist(gen); // Poisson sample to get number of randoms for the cell

				int index = k + N.z*(j + N.y*i);
				for (int ran = 0; ran < numRans; ++ran) {
					writeFile << r_x + uDist(gen)*Delta_r.x << " " << r_y + uDist(gen)*Delta_r.y << " " << r_z + uDist(gen)*Delta_r.z << " \n";
				}
				if (index < N.w) {
					grid[index] += numRans;
				}
				count += numRans;
			}
		}
	}
	writeFile.close();
}

std::vector<double> fftFrequencies(int N, double L) {
    std::vector<double> k(N);
    double dk = 2.0*M_PI/L;
    for (int i = 0; i <= N/2; ++i)
        k[i] = i*dk;
    for (int i = N/2 + 1; i < N; ++i)
        k[i] = (i - N)*dk;
    return k;
}

int main(int argc, char *argv[]) {
	using namespace std;

	if(argc < 3) {
		cout << "Error: Need Real Input and Random Input Filenames\n";
		getchar();
	}
	else {
		string realInputName(argv[1]);
		string randInputName(argv[2]);
		cout << "Real Data To Read: " << realInputName << "\n";
		cout << "Random Data To Read/Write: " << randInputName << "\n";
		cout << "Press Key to Continue\n";
		getchar();
		
		//Prepare Arrays and Values
		int4 N = {512, 512, 512, 0};
		N.w = N.x*N.y*N.z;
		double3 L = {1024.0, 1024.0, 1024.0};
		double3 Delta_r = {L.x/N.x, L.y/N.y, L.z/N.z};

		int count = 0;
		vector<double> realGrid(N.w);
		
		vector<fftw_complex> dk(N.x*N.y*(N.z/2 + 1));
		fftw_init_threads();
		fftw_import_wisdom_from_filename("fftwWisdom.dat");
		fftw_plan_with_nthreads(omp_get_max_threads());
		fftw_plan dr2dk = fftw_plan_dft_r2c_3d(N.x, N.y, N.z, realGrid.data(), dk.data(), FFTW_MEASURE);
		fftw_export_wisdom_to_filename("fftwWisdom.dat");
		
		
		ifstream inputFile (realInputName);
		if(inputFile.is_open()) {
			//Bin Galaxies to 512x512x512 bin grid
			binPoints(inputFile, realGrid, N, L, count);
		}
		else {
			cout << "Unable to read " << realInputName << "\n";
			getchar();
			return 0;
		}
				
		//Generate Random Points in second grid
		vector<double> randGrid(N.w);
		int ranCount = 0;
		double nbar = 1.0 * count / (L.x * L.y * L.z);

		//Trying to read from rand file.
		ifstream randFileRead (randInputName);
		if(randFileRead.is_open()) {
			binPoints(randFileRead, randGrid, N, L, ranCount);
		}
		else {
			cout << "Unable to read " << randInputName << ", Generating Random Points\n";
			//Trying to write to rand file.
			ofstream randFileWrite (randInputName);
			if(randFileWrite.is_open()){
				generateAndBinPoints(randFileWrite, randGrid, N, L, ranCount, TIMESRAN * nbar * Delta_r.x * Delta_r.y * Delta_r.z);
			}
			else {
				cout << "Unable to write to " << realInputName << "\n";
				getchar();
				return 0;
			}
		}
		
		cout << "Read " << count << " real points and read/generated " << ranCount << " random points\n";
		
		double alpha = (1.0 * count) / ranCount;

		//Take weighted difference
		#pragma omp parallel for
		for (int i = 0; i < N.w; ++i) {
			realGrid[i] -= alpha * randGrid[i];
		}

		//Feed 1D Representation of grid through fourier tranform, complex array is generated

		fftw_execute(dr2dk);

		//Calculate Frequences for each bin, bin frequencies, calculate power spectrum from frequencies
		vector<double> k_x = fftFrequencies(N.x, L.x);
		vector<double> k_y = fftFrequencies(N.y, L.y);
		vector<double> k_z = fftFrequencies(N.z, L.z);

		vector<double> Pk(N_k, 0.0);
		vector<int> Nk(N_k, 0);
		double Pshot = count + alpha * alpha * ranCount;

		double Delta_k = (k_max - k_min)/N_k;
		for (int i = 0; i < N.x; ++i) {
			for (int j = 0; j < N.y; ++j) {
			for (int k = 0; k <= N.z/2; ++k) {
				int index = k + (N.z/2 + 1)*(j + N.y*i);
				double k_mag = sqrt(k_x[i]*k_x[i] + k_y[j]*k_y[j] + k_z[k]*k_z[k]);
				int bin = (k_mag - k_min)/Delta_k;
				if(bin < N_k && bin >= 0) {
					Pk[bin] += (dk[index][0]*dk[index][0] + dk[index][1]*dk[index][1]) - Pshot;
					++Nk[bin];
				}
			}
			}
		}

		for (int i = 0; i < N_k; ++i) {
			Pk[i] = Pk[i] / (Nk[i] * nbar * count);
			cout << k_min + (i + 0.5)*Delta_k << "," << Pk[i] << "\n";
		}
	}
	return 0;
}
