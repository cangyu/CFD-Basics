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
	Array2D(size_t nx, size_t ny, double val = 0.0) :
		m_data(nx*ny, val)
	{
		m_Nx = nx;
		m_Ny = ny;
	}

	// 0-based indexing
	double& at(int i, int j)
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
	double& operator()(int i, int j)
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
const double a = 340.28; // m/s
const double Ma = 4.0;
const double R = 287.0;
const double Cv = 718.0;

const size_t IMIN = 1, IMAX = 70;
const size_t JMIN = 1, JMAX = 70;
const double dx = LHORI / (IMAX - 1);
const double delta = 5 * LHORI / std::sqrt(Re);
const double LVERT = 5 * delta;
const double dy = LVERT / (JMAX - 1);
vector<double> x(IMAX, 0.0), y(JMAX, 0.0);

const double CFL = 0.7;
double dt = 1e-5;
double t = 0.0;

size_t iter = 0;
const size_t MAX_ITER = 2000;

const double u_inf = Ma * a;
const double v_inf = 0.0;
const double p_inf = 101325.0;
const double T_inf = 215.0;
const double Tw = 1.0 * T_inf;

// Primitive variables
Array2D rho(IMAX, JMAX, 0.0);
Array2D u(IMAX, JMAX, 0.0);
Array2D v(IMAX, JMAX, 0.0);
Array2D p(IMAX, JMAX, 0.0);
Array2D T(IMAX, JMAX, Tw);
Array2D e(IMAX, JMAX, 0.0); // Internal energy per unit mass

// Conservative variables
Array2D U1(IMAX, JMAX, 0.0); // rho
Array2D U2(IMAX, JMAX, 0.0); // rho u
Array2D U3(IMAX, JMAX, 0.0); // rho v
Array2D U5(IMAX, JMAX, 0.0); // rho(e+V^2 / 2)


double Sutherland(double T)
{
	static const double mu0 = 1.7894e-5; // Kg/(m*s)
	static const double T0 = 288.16; // K

	return mu0 * std::pow(T / T0, 1.5) * (T0 + 110.0) / (T + 110.0);
}

void init()
{
	/********************************** Grid **********************************/
	for (size_t i = 1; i < IMAX; ++i)
		x[i] = x[i - 1] + dx;

	for (size_t j = 1; j < JMAX; ++j)
		y[j] = y[j - 1] + dy;

	/********************************** I.C. **********************************/
	for (int j = 2; j <= JMAX; ++j)
		for (int i = 2; i <= IMAX; ++i)
		{
			u(i, j) = u_inf;
			v(i, j) = v_inf;
			p(i, j) = p_inf;
			T(i, j) = T_inf;
		}

	/********************************** B.C. **********************************/
	// Front tip
	u(1, 1) = 0.0;
	v(1, 1) = 0.0;
	p(1, 1) = p_inf;
	T(1, 1) = T_inf;

	// Inlet
	for (int j = 2; j <= JMAX; ++j)
	{
		u(1, j) = u_inf;
		v(1, j) = 0.0;
		p(1, j) = p_inf;
		T(1, j) = T_inf;
	}

	// Top(Far)
	for (int i = 2; i <= IMAX; ++i)
	{
		u(i, JMAX) = u_inf;
		v(i, JMAX) = 0.0;
		p(i, JMAX) = p_inf;
		T(i, JMAX) = T_inf;
	}

	// Bottom
	for (int i = 1; i <= IMAX; ++i)
	{
		u(i, 1) = 0.0;
		v(i, 1) = 0.0;
		p(i, 1) = 2 * p(i, 2) - p(i, 3);
		T(i, 1) = Tw;
	}

	// Outlet
	for (int j = 2; j <= JMAX - 1; ++j)
	{
		u(IMAX, j) = 2 * u(IMAX - 1, j) - u(IMAX - 2, j);
		v(IMAX, j) = 2 * v(IMAX - 1, j) - v(IMAX - 2, j);
		p(IMAX, j) = 2 * p(IMAX - 1, j) - p(IMAX - 2, j);
		T(IMAX, j) = 2 * T(IMAX - 1, j) - T(IMAX - 2, j);
	}

	/************************** Derived Variables *****************************/
	for (int j = 1; j <= JMAX; ++j)
		for (int i = 1; i <= IMAX; ++i)
		{
			rho(i, j) = p(i, j) / (R*T(i, j));
			e(i, j) = Cv * T(i, j);
		}

	/************************ Conservative Variables **************************/
	for (int j = 1; j <= JMAX; ++j)
		for (int i = 1; i <= IMAX; ++i)
		{
			U1(i, j) = rho(i, j);
			U2(i, j) = rho(i, j)*u(i, j);
			U3(i, j) = rho(i, j)*v(i, j);
			const double K = 0.5*(pow(u(i, j), 2) + pow(v(i, j), 2));
			U5(i, j) = rho(i, j)*(e(i, j) + K);
		}
}

void MacCormack()
{
	/***************************** Forward Difference *************************/
	Array2D dU1dt(IMAX, JMAX, 0.0);
	for (int j = 2; j < JMAX; ++j)
		for (int i = 2; i < IMAX; ++i)
		{
			const double dEdx = (rho(i + 1, j)*u(i + 1, j) - rho(i, j)*u(i, j)) / dx;
			const double dFdy = (rho(i, j + 1)*v(i, j + 1) - rho(i, j)*v(i, j)) / dy;
			dU1dt(i, j) -= (dEdx + dFdy);
		}
	Array2D dU2dt(IMAX, JMAX, 0.0);
	for (int j = 2; j < JMAX; ++j)
		for (int i = 2; i < IMAX; ++i)
		{

		}
	Array2D dU3dt(IMAX, JMAX, 0.0);
	for (int j = 2; j < JMAX; ++j)
		for (int i = 2; i < IMAX; ++i)
		{

		}
	Array2D dU5dt(IMAX, JMAX, 0.0);
	for (int j = 2; j < JMAX; ++j)
		for (int i = 2; i < IMAX; ++i)
		{

		}

}

void write_tecplot(size_t n)
{

}

void write_user(size_t n)
{

}

void output()
{
	write_tecplot(iter);
	write_user(iter);
}

bool check_convergence()
{
	return iter > MAX_ITER;
}

// Explicit Time-Marching
void solve()
{
	bool converged = false;
	while (!converged)
	{
		++iter;
		t += dt;
		MacCormack();
		output();
		converged = check_convergence();
	}
}

int main(int argc, char* argv[])
{
	init();
	output();
	solve();
	return 0;
}
