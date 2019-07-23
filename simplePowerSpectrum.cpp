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
double const ALPHA = 0.1;

double const k_min = 0.01;
double const k_max = 0.2;
int const N_k = 19;

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
			//Prepare Arrays and Values
			int4 N = {512, 512, 512, 0};
			N.w = N.x*N.y*N.z;
			double3 L = {1024.0, 1024.0, 1024.0};
			double3 Delta_r = {L.x/N.x, L.y/N.y, L.z/N.z};
			vector<double> realGrid(N.w);

			vector<fftw_complex> dk(N.x*N.y*(N.z/2 + 1));
			fftw_init_threads();
			fftw_import_wisdom_from_filename("fftwWisdom.dat");
			fftw_plan_with_nthreads(omp_get_max_threads());
			fftw_plan dr2dk = fftw_plan_dft_r2c_3d(N.x, N.y, N.z, realGrid.data(), dk.data(), FFTW_MEASURE);
			fftw_export_wisdom_to_filename("fftwWisdom.dat");

			//Bin Galaxies to 512x512x512 bin grid
			int subPos1, subPos2;
			int count = 0;
			string line, sub;
			while(getline(inputFile,line)) {
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
					for(int k = 0; k < N.z; ++k) {
						int numRans = pDist(gen); // Poisson sample to get number of randoms for the cell
						
        					int index = k + N.z*(j + N.y*i);
						if (index < N.w) {
						    randGrid[index] += numRans;
						}
						ranCount += numRans;
					}
				}
			}

			//Take weighted difference
			#pragma omp parallel for
			for (int i = 0; i < N.w; ++i) {
				realGrid[i] -= ALPHA * randGrid[i];
			}

			//Feed 1D Representation of grid through fourier tranform, complex array is generated
			
			fftw_execute(dr2dk);

			//Calculate Frequences for each bin, bin frequencies, calculate power spectrum from frequencies
			vector<double> k_x = fftFrequencies(N.x, L.x);
			vector<double> k_y = fftFrequencies(N.y, L.y);
			vector<double> k_z = fftFrequencies(N.z, L.z);

			vector<double> Pk(N_k, 0.0);
			vector<int> Nk(N_k, 0);
			double Pshot = count + ALPHA * ALPHA * ranCount;

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
				Pk[i] = Pk[i] / (Nk[i] * L.x * L.y * L.z);
				cout << Pk[i] << ",";
			}

		}
		else {
			cout << "Unable to read " << inputName << "\n";
			getchar();
		}
	}
	return 0;
}
