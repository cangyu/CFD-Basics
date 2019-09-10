#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

class Array2D
{
private:
	vector<double> m_data;
	size_t m_Nx, m_Ny;

public:
	Array2D(size_t nx, size_t ny, double val = 0.0) : m_Nx(nx), m_Ny(ny), m_data(nx*ny, val) {}

	// 0-based indexing
	double &at(int i, int j)
	{
		int idx = i + m_Nx * j;
		return m_data[idx];
	}

	double at(int i, int j) const
	{
		int idx = i + m_Nx * j;
		return m_data[idx];
	}

	// 1-based indexing
	double &operator()(int i, int j)
	{
		return at(i - 1, j - 1);
	}

	double operator()(int i, int j) const
	{
		return at(i - 1, j - 1);
	}
};

const double LHORI = 1e-5; // m
const double Re = 1000.0;
const double Pr = 0.71;
const double Ma = 4.0;
const double R = 287.0;
const size_t IMIN = 1, IMAX = 70;
const size_t JMIN = 1, JMAX = 70;
const double dx = LHORI / (IMAX - 1);
const double delta = 5 * LHORI / std::sqrt(Re);
const double LVERT = 5 * delta;
const double dy = LVERT / (JMAX - 1);
const double a = 340.0; // m/s
vector<double> x(IMAX, 0.0), y(JMAX, 0.0);
const double CFL = 0.7;

const double u_inf = a * Ma;
const double v_inf = 0.0;
const double p_inf = 101325.0;
const double T_inf = 215.0, Tw = 300.0;

Array2D u(IMAX, JMAX, 0.0);
Array2D v(IMAX, JMAX, 0.0);
Array2D p(IMAX, JMAX, 0.0);
Array2D T(IMAX, JMAX, Tw);

void init()
{
	for (size_t i = 1; i < IMAX; ++i)
		x[i] = x[i - 1] + dx;

	for (size_t j = 1; j < JMAX; ++j)
		y[j] = y[j - 1] + dy;
}

void solve()
{

}

int main(int argc, char *argv[])
{
	init();

	solve();

	return 0;
}
