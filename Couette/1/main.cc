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

	double at(int i, int j) const
	{
		int idx = i + Nx * j;
		return data[idx];
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

	void fill(double x)
	{
		std::fill(data.begin(), data.end(), x);
	}

	double max() const
	{
		return *max_element(data.begin(), data.end());
	}

	double min() const
	{
		return *min_element(data.begin(), data.end());
	}
};

const size_t WIDTH = 12;
const size_t DIGITS = 7;

const double L = 0.5; // m
const double D = 0.01; // m
const double Ue = 1.0; // m/s
const double rho = 0.002377; // Kg/m3
const double Re = 63.6;
const double mu = 3.737e-7; // Kg/m/s

const int Nx = 21, Ny = 11;
const double dx = L / (Nx - 1), dy = D / (Ny - 1);
const double dx2 = 2 * dx, dy2 = 2 * dy;
const double dxdx = dx * dx, dydy = dy * dy;
vector<double> x(Nx, 0.0), y(Ny, 0.0);

const double dt = 0.001;
double t = 0.0;
int iter_cnt = 0;
const int MAX_ITER_NUM = 500;

const double a = 2 * (dt / dxdx + dt / dydy);
const double b = -dt / dxdx;
const double c = -dt / dydy;

const double alpha_p = 0.1;

ofstream result;

Array2D p(Nx, Ny, 0.0), p_star(Nx, Ny, 0.0), p_prime(Nx, Ny, 0.0);
Array2D u(Nx + 1, Ny, 0.0), u_star(Nx + 1, Ny, 0.0), u_prime(Nx + 1, Ny, 0.0);
Array2D v(Nx + 2, Ny + 1, 0.0), v_star(Nx + 2, Ny + 1, 0.0), v_prime(Nx + 2, Ny + 1, 0.0);

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
	cout << "mu=" << mu << endl;

	// Coordinates
	for (int i = 1; i < Nx; ++i)
		x[i] = L * i / (Nx - 1);
	for (int j = 1; j < Ny; ++j)
		y[j] = D * j / (Ny - 1);

	// Velocity field
	for (int i = 1; i <= Nx + 1; ++i)
		u_star(i, Ny) = u(i, Ny) = Ue;

	v(15, 5) = v_star(15, 5) = 0.5; // Initial peak to ensure 2D flow structure

	// Header
	result.open("flow.txt");
	result << Nx << '\t' << Ny << endl;
	for (int i = 0; i < Nx; ++i)
		result << ' ' << x[i];
	result << endl;
	for (int j = 0; j < Ny; ++j)
		result << ' ' << y[j];
	result << endl;
	write_result();
}

void finalize(void)
{
	result.close();
}

void JacobiMethod()
{
	for (int r = 0; r < 200; ++r)
	{
		auto pp = p_prime;
		for (int i = 2; i < Nx; ++i)
			for (int j = 2; j < Ny; ++j)
			{
				double d = (rho*u_star(i + 1, j) - rho * u_star(i, j)) / dx + (rho*v_star(i + 1, j + 1) - rho * v_star(i + 1, j)) / dy;
				p_prime(i, j) = -(b * pp(i + 1, j) + b * pp(i - 1, j) + c * pp(i, j + 1) + c * pp(i, j - 1) + d) / a;
			}

		//cout << "\tInternal loop:" << r + 1 << ": Max(p\')=" << p_prime.max() << " Min(p\')=" << p_prime.min() << endl;
	}
}

void GaussSeidelMethod()
{
	for (int r = 0; r < 200; ++r)
	{
		for (int i = 2; i < Nx; ++i)
			for (int j = 2; j < Ny; ++j)
			{
				double d = (rho*u_star(i + 1, j) - rho * u_star(i, j)) / dx + (rho*v_star(i + 1, j + 1) - rho * v_star(i + 1, j)) / dy;
				p_prime(i, j) = -(b * p_prime(i + 1, j) + b * p_prime(i - 1, j) + c * p_prime(i, j + 1) + c * p_prime(i, j - 1) + d) / a;
			}
	}
}

void SIMPLE(void)
{
	// u_star at inner points
	for(int j = 2; j <= Ny-1; ++j)
		for (int i = 2; i <= Nx; ++i)
		{
			double v_bar1 = 0.5*(v(i, j + 1) + v(i + 1, j + 1));
			double v_bar2 = 0.5*(v(i, j) + v(i + 1, j));
			
			double t11 = rho * pow(u(i + 1, j), 2) - rho * pow(u(i - 1, j), 2);
			double t12 = rho * u(i, j + 1)*v_bar1 - rho * u(i, j - 1)*v_bar2;
			double t21 = u(i + 1, j) - 2 * u(i, j) + u(i - 1, j);
			double t22 = u(i, j + 1) - 2 * u(i, j) + u(i, j - 1);
			double A_star = -(t11 / dx2 + t12 / dy2) + mu * (t21 / dxdx + t22 / dydy);
			
			double loc_rhou_star = rho * u_star(i, j) + A_star * dt - dt / dx * (p_star(i, j) - p_star(i - 1, j));
			u_star(i, j) = loc_rhou_star / rho;
		}

	// u_star at virtual points, transparent
	for (int j = 1; j <= Ny; ++j)
	{
		u_star(1, j) = u_star(2, j); // Inlet
		u_star(Nx + 1, j) = u_star(Nx, j); // Outlet
	}

	// v_star at inner points
	for (int i = 2; i <= Nx + 1; ++i)
		for (int j = 2; j <= Ny; ++j)
		{
			double u_bar1 = 0.5 *(u(i, j - 1) + u(i, j));
			double u_bar2 = 0.5 *(u(i - 1, j - 1) + u(i - 1, j));

			double t11 = rho * v(i + 1, j) * u_bar1 - rho * v(i - 1, j) * u_bar2;
			double t12 = rho * pow(v(i, j + 1), 2) - rho * pow(v(i, j - 1), 2);
			double t21 = v(i + 1, j) - 2 * v(i, j) + v(i - 1, j);
			double t22 = v(i, j + 1) - 2 * v(i, j) + v(i, j - 1);
			double B_star = -(t11 / dx2 + t12 / dy2) + mu * (t21 / dxdx + t22 / dydy);

			double loc_rhov_star = rho * v_star(i, j) + B_star * dt - dt / dy * (p_star(i - 1, j) - p_star(i - 1, j - 1));
			v_star(i, j) = loc_rhov_star / rho;
		}	

	// v_star at virtual points
	for (int j = 1; j <= Ny + 1; ++j)
	{
		v_star(Nx + 2, j) = v_star(Nx + 1, j); // Outlet
	}

	// Solve p_prime at inner points
	JacobiMethod();

	// Correct p
	for (int i = 2; i < Nx; ++i)
		for (int j = 2; j < Ny; ++j)
		{
			// Update with relaxation
			p(i, j) = p_star(i, j) + alpha_p * p_prime(i, j);
			
			// For next round
			p_star(i, j) = p(i, j);
		}

	// Correct u
	for (int j = 2; j <= Ny - 1; ++j)
		for (int i = 1; i <= Nx + 1; ++i)
			u(i, j) = u_star(i, j);

	// Correct v
	for (int i = 2; i <= Nx + 2; ++i)
		for (int j = 2; j <= Ny; ++j)
			v(i, j) = v_star(i, j);
}

bool check_convergence(void)
{
	cout << "Max(u)=" << u.max() << " Min(u)=" << u.min() << endl;
	cout << "Max(v)=" << v.max() << " Min(v)=" << v.min() << endl;
	cout << "Max(p)=" << p.max() << " Min(p)=" << p.min() << endl;

	return iter_cnt > MAX_ITER_NUM;
}

void solve(void)
{
	bool converged = false;
	while (!converged)
	{
		++iter_cnt;
		cout << "Iter" << iter_cnt << ":" << endl;
		SIMPLE();
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
