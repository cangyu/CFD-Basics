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

class Array2D
{
private:
	vector<double> data;
	size_t Nx, Ny;

public:
	Array2D(size_t nx, size_t ny, double val = 0.0) :
		Nx(nx),
		Ny(ny),
		data(nx*ny, val)
	{
		// Empty body
	}

	// 0-based indexing
	double &at(int i, int j)
	{
		int idx = i + Nx * j;
		return data[idx];
	}

	// 1-based indexing
	double &operator()(int i, int j)
	{
		return at(i - 1, j - 1);
	}
};

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
const int MAX_ITER_NUM = 2000;

ofstream result;

Array2D p(Nx, Ny), u(Nx + 1, Ny), v(Nx + 2, Ny + 1);
auto u_star = u, u_prime = u;
auto v_star = v, v_prime = v;
auto p_star = p, p_prime = p;

void write_result(void)
{
	result << iter_cnt << ' ' << t;
	
	// Pressure
	for(int j = 1; j <= Ny; ++j)
	{
		result << endl;
		for(int i = 1; i <= Nx; ++i)
			result << setw(WIDTH) << setprecision(DIGITS) << p(i, j);
	}

	result << endl;
}

void init(void)
{
	// Coordinate
	for (int i = 1; i < Nx; ++i)
		x[i] = x[i - 1] + dx;
	for (int j = 1; j < Ny; ++j)
		y[j] = y[j - 1] + dy;

	// Velocity field
	for(int i = 1; i <= Nx; ++i)
		u(i, Ny) = Ue;

	v_star(15, 5) = 0.5*0.3048; // Initial peak

	// Header
	result.open("flow.txt");
	result << Nx << '\t' << Ny << endl;
	for (int i = 0; i < Nx; ++i)
		result << x[i] << '\t';
	result << endl;
	for (int j = 0; j < Ny; ++j)
		result << y[j] << '\t';
	result << endl;
	write_result();
}

void finalize(void)
{
	result.close();
}

void loop(void)
{
	// ru_star
	for(int j = 2; j <= Ny-1; ++j)
		for (int i = 2; i <= Nx; ++i)
		{

		}
}

bool check_convergence(void)
{
	return iter_cnt > MAX_ITER_NUM;
}

void solve(void)
{
	bool converged = false;
	while (!converged)
	{
		++iter_cnt;
		cout << "Iter" << iter_cnt << ":" << endl;
		loop();
		t += dt;
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
