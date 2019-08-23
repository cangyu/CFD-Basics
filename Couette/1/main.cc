#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstddef>
#include <cassert>
#include <limits>
#include <Eigen/Sparse>

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

const size_t WIDTH = 16;
const size_t DIGITS = 7;

const double L = 0.5; // m
const double D = 0.01; // m
const double Ue = 1.0; // m/s
const double rho = 1.225; // Kg/m3
const double mu = 3.737e-5; // Kg/m/s

const int Nx = 21, Ny = 11;
const double dx = L / (Nx - 1), dy = D / (Ny - 1);
const double dx2 = 2 * dx, dy2 = 2 * dy;
const double dxdx = dx * dx, dydy = dy * dy;
vector<double> x(Nx, 0.0), y(Ny, 0.0);

const double dt = 0.001;
double t = 0.0;
int iter_cnt = 0;
const int MAX_ITER_NUM = 1001;

const double a = 2 * (dt / dxdx + dt / dydy);
const double b = -dt / dxdx;
const double c = -dt / dydy;
double d_min = numeric_limits<double>::max(), d_max = numeric_limits<double>::min();

const double alpha_p = 0.1;

Array2D p(Nx, Ny, 0.0), p_star(Nx, Ny, 0.0), p_prime(Nx, Ny, 0.0);
Array2D u(Nx + 1, Ny, 0.0), u_star(Nx + 1, Ny, 0.0), u_prime(Nx + 1, Ny, 0.0);
Array2D v(Nx + 2, Ny + 1, 0.0), v_star(Nx + 2, Ny + 1, 0.0), v_prime(Nx + 2, Ny + 1, 0.0);

void output(void)
{
	Array2D u_interp(Nx, Ny, 0.0);
	for (int i = 1; i <= Nx; ++i)
		u_interp(i, 1) = 0.0; // Bottom
	for (int j = 2; j <= Ny - 1; ++j)
		for (int i = 1; i <= Nx; ++i)
			u_interp(i, j) = (u(i, j) + u(i + 1, j)) / 2; // Inner
	for (int i = 1; i <= Nx; ++i)
		u_interp(i, Ny) = Ue; // Top

	Array2D v_interp(Nx, Ny, 0.0);
	for (int i = 1; i <= Nx; ++i)
		v_interp(i, 1) = 0.0; // Bottom
	for (int j = 2; j <= Ny - 1; ++j)
	{
		v_interp(1, j) = 0.0; // Left
		for (int i = 3; i <= Nx + 1; ++i)
			v_interp(i - 1, j) = (v(i, j) + v(i, j + 1)) / 2; // Inner and Right
	}
	for (int i = 1; i <= Nx; ++i)
		v_interp(i, Ny) = 0.0; // Top

	// Create Tecplot data file.
	ofstream result("flow" + to_string(iter_cnt) + ".dat");
	if (!result)
		throw("Failed to create datafile!");

	// Header
	result << "TITLE = \"t=" << t << "\"" << endl;
	result << "VARIABLES = \"X\", \"Y\", \"P\", \"U\", \"V\"" << endl;
	result << "ZONE I=" << Nx << ", J=" << Ny << ", F=POINT" << endl;

	// Flowfield data
	for(int j=1; j<=Ny; ++j)
		for (int i = 1; i <= Nx; ++i)
		{
			result << setw(WIDTH) << setprecision(DIGITS) << x[i - 1];
			result << setw(WIDTH) << setprecision(DIGITS) << y[j - 1];
			result << setw(WIDTH) << setprecision(DIGITS) << p(i, j);
			result << setw(WIDTH) << setprecision(DIGITS) << u_interp(i, j);
			result << setw(WIDTH) << setprecision(DIGITS) << v_interp(i, j);
			result << endl;
		}

	// Finalize
	result.close();
}

void JacobiMethod() // Only suitable for Dirichlet B.C. and will blow up with Neumann B.C.
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
	}
}

void GaussSeidelMethod() // Only suitable for Dirichlet B.C. and will blow up with Neumann B.C.
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

void ImplicitMethod() // Can handle both Dirichlet and Neumann B.C.
{
	typedef Eigen::SparseMatrix<double> SpMat;
	typedef Eigen::Triplet<double> T;

	const int m = Nx * Ny;
	vector<T> coef;
	Eigen::VectorXd rhs(m);
	SpMat A(m, m);

	// Calculating coefficients
	for (int i = 0; i < Nx; ++i)
		for (int j = 0; j < Ny; ++j)
		{
			const int id = j * Nx + i;
			const int id_w = id - 1;
			const int id_e = id + 1;
			const int id_n = id + Nx;
			const int id_s = id - Nx;

			if (i == 0 || i == Nx - 1) // Inlet and Outlet
			{
				coef.push_back(T(id, id, 1.0));
				rhs(id) = 0.0;
			}
			else if (j == 0) // Bottom
			{
				coef.push_back(T(id, id, 1.0));
				coef.push_back(T(id, id_n, -1.0));
				rhs(id) = 0.0;
			}
			else if (j == Ny - 1) // Top
			{
				coef.push_back(T(id, id, 1.0));
				coef.push_back(T(id, id_s, -1.0));
				rhs(id) = 0.0;
			}
			else // Inner
			{
				// Use 0-based interface
				const double d = (rho*u_star.at(i + 1, j) - rho * u_star.at(i, j)) / dx + (rho*v_star.at(i + 1, j + 1) - rho * v_star.at(i + 1, j)) / dy;

				coef.push_back(T(id, id, a));
				coef.push_back(T(id, id_w, b));
				coef.push_back(T(id, id_e, b));
				coef.push_back(T(id, id_n, c));
				coef.push_back(T(id, id_s, c));
				rhs(id) = -d;
			}
		}

	// Construct sparse matrix
	A.setFromTriplets(coef.begin(), coef.end());

	// Solve the linear system: Ax = rhs
	Eigen::SimplicialCholesky<SpMat> chol(A);
	Eigen::VectorXd x = chol.solve(rhs);

	// Update p_prime at inner
	for (int i = 1; i < Nx - 1; ++i)
		for (int j = 1; j < Ny - 1; ++j)
		{
			const int id = j * Nx + i;
			p_prime.at(i, j) = x(id);
		}

	// Enforce B.C. of p_prime: zero at inlet and outlet, zero-gradient at top and bottom.
	for (int j = 1; j <= Ny; ++j)
	{
		p_prime(1, j) = 0.0;
		p_prime(Nx, j) = 0.0;
	}
	for (int i = 2; i <= Nx - 1; ++i)
	{
		p_prime(i, 1) = p_prime(i, 2);
		p_prime(i, Ny) = p_prime(i, Ny - 1);
	}
}

void SIMPLE(void)
{
	// u_star at inner points
	for (int j = 2; j <= Ny - 1; ++j)
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

	// v_star at inner points
	for (int i = 3; i <= Nx + 1; ++i)
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

	// Solve p_prime
	// JacobiMethod();
	// GaussSeidelMethod();
	ImplicitMethod();

	// Correct p
	for (int i = 1; i <= Nx; ++i)
		for (int j = 1; j <= Ny; ++j)
		{
			// Update with relaxation
			p(i, j) = p_star(i, j) + alpha_p * p_prime(i, j);

			// For next round
			p_star(i, j) = p(i, j);
		}

	// Correct u at inner nodes
	for (int j = 2; j <= Ny - 1; ++j)
		for (int i = 2; i <= Nx; ++i)
			u(i, j) = u_star(i, j);

	// Linear extrapolation of u at virtual nodes
	for (int j = 2; j <= Ny - 1; ++j)
	{
		u(1, j) = 2 * u(2, j) - u(3, j);
		u(Nx + 1, j) = 2 * u(Nx, j) - u(Nx - 1, j);
	}

	// Correct v at inner nodes
	for (int i = 3; i <= Nx + 1; ++i)
		for (int j = 2; j <= Ny; ++j)
			v(i, j) = v_star(i, j);

	// Linear extrapolation of v at both top and bottom virtual nodes
	// No-Penetration at both top and bottom
	for (int i = 2; i <= Nx + 1; ++i)
	{
		v(i, 1) = -v(i, 2);
		v(i, Ny + 1) = -v(i, Ny);
	}

	// Linear extrapolation of v at right virtual nodes
	for (int j = 2; j <= Ny; ++j)
		v(Nx + 2, j) = 2 * v(Nx + 1, j) - v(Nx, j);
}

bool check_convergence(void)
{
	// Statistics of u
	double u_max = numeric_limits<double>::min();
	double u_min = numeric_limits<double>::max();
	for (int i = 2; i <= Nx; ++i)
		for (int j = 1; j <= Ny; ++j)
		{
			u_max = max(u_max, u(i, j));
			u_min = min(u_min, u(i, j));
		}
	cout << "Max(u)=" << u_max << " Min(u)=" << u_min << endl;

	// Statistics of v
	double v_max = numeric_limits<double>::min();
	double v_min = numeric_limits<double>::max();
	for (int i = 2; i <= Nx + 1; ++i)
		for (int j = 2; j <= Ny; ++j)
		{
			v_max = max(v_max, v(i, j));
			v_min = min(v_min, v(i, j));
		}
	cout << "Max(v)=" << v_max << " Min(v)=" << v_min << endl;

	// Statistics of p
	double p_max = numeric_limits<double>::min();
	double p_min = numeric_limits<double>::max();
	for (int i = 1; i <= Nx; ++i)
		for (int j = 1; j <= Ny; ++j)
		{
			p_max = max(p_max, p(i, j));
			p_min = min(p_min, p(i, j));
		}
	cout << "Max(p)=" << p_max << " Min(p)=" << p_min << endl;

	return iter_cnt > MAX_ITER_NUM;
}

int main(int argc, char *argv[])
{
	cout << "mu=" << mu << endl;
	cout << "dt=" << dt << endl;

	// Init
	for (int i = 1; i < Nx; ++i)
		x[i] = L * i / (Nx - 1); // X-Coordinates
	for (int j = 1; j < Ny; ++j)
		y[j] = D * j / (Ny - 1); // Y-Coordinates

	for (int i = 1; i <= Nx + 1; ++i)
		u_star(i, Ny) = u(i, Ny) = Ue; // U at top
	v(15, 5) = v_star(15, 5) = 0.5; // Initial peak to ensure 2D flow structure

	// Solve
	output(); // I.C.
	bool converged = false;
	while (!converged)
	{
		++iter_cnt;
		cout << "Iter" << iter_cnt << ":" << endl;
		SIMPLE();
		t += dt;
		output();
		converged = check_convergence();
	}

	return 0;
}
