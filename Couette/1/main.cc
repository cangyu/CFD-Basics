#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstddef>
#include <cassert>

using namespace std;

const size_t WIDTH = 16;
const size_t DIGITS = 6;

const double L = 0.5*0.3048; // m
const double D = 0.01*0.3048; // m
const double Ue = 1.0 * 0.3048; // m/s
const double rho = 1.225055; // Kg/m3
const double Re = 63.6;

const int Nx = 21, Ny = 11;
const double dx = L / (Nx - 1), dy = D / (Ny - 1);
vector<double> x(Nx, 0.0), y(Ny, 0.0);

const double dt = 0.001;
double t = 0.0;
int iter_cnt = 0;

ofstream result;

vector<double> u(N, 0.0);

void write_result(void)
{
	result << iter_cnt << '\t' << t << endl;
	for (int i = 0; i < N; ++i)
	{
		result << setw(WIDTH) << setprecision(DIGITS) << u[i];
		result << endl;
	}
}

void init(void)
{
	// Coordinate
	for (int j = 1; j < N; ++j)
		y[j] = y[j - 1] + dy;

	// Velocity field
	u[N - 1] = 1.0;

	// Header
	result.open("flow.txt");
	result << N << endl;
	for (int i = 0; i < N; ++i)
		result << y[i] << '\t';
	result << endl;
	write_result();
}

void finalize(void)
{
	result.close();
}

void loop(void)
{

}

bool check_convergence(void)
{

}

void solve(void)
{
	bool converged = false;
	while (!converged)
	{
		++iter_cnt;
		cout << "Iter" << iter_cnt << ":" << endl;
		loop();
		write_result();
		converged = check_convergence();
	}
}

int main(int argc, char *argv[])
{
	init();
	solve();
	finalize();
	return 0;
}
