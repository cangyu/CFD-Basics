#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

const int Nx = 101, Ny = 81;
const double Lx = 1.0, Ly = 0.8;
const double xL = -Lx / 2, xR = Lx / 2;
const double yL = -Ly / 2, yR = Ly / 2;
const double dx = Lx / (Nx - 1), dy = Ly / (Ny - 1);
const double dx2 = pow(dx, 2), dy2 = pow(dy, 2);

const double rho = 0.996873e3; // kg/m^3
const double g = 9.80665; // m/s^2
const double p0 = 101325.0; // Pa
const double u0 = 1.0; // m/s
const double mu = 1e-5; // Pa*s

double dt = 1e-4;

vector<double> x(Nx, xL), y(Ny, yL); // m

vector<vector<double>> p(Nx, vector<double>(Ny, p0)); // Pa
vector<vector<double>> u(Nx - 1, vector<double>(Ny, 0.0)); // m/s
vector<vector<double>> v(Nx, vector<double>(Ny - 1, 0.0)); // m/s

vector<vector<double>> p_star(Nx, vector<double>(Ny, p0)); // Pa
vector<vector<double>> u_star(Nx - 1, vector<double>(Ny, 0.0)); // m/s
vector<vector<double>> v_star(Nx, vector<double>(Ny - 1, 0.0)); // m/s

vector<vector<double>> p_prime(Nx, vector<double>(Ny, p0)); // Pa
vector<vector<double>> u_prime(Nx - 1, vector<double>(Ny, 0.0)); // m/s
vector<vector<double>> v_prime(Nx, vector<double>(Ny - 1, 0.0)); // m/s

double relaxation(double a, double b, double alpha)
{
	return (1 - alpha) *a + alpha * b;
}

void setup_grid(void)
{
	for (int i = 1; i < Nx; ++i)
		x[i] = x[i - 1] + dx;

	for (int j = 1; j < Ny; ++j)
		y[j] = y[j - 1] + dy;
}

void initialize_flowfield(void)
{
	const double dp = rho * g * dy;

	for (int j = Ny - 2; j >= 0; --j)
		for (int i = 0; i < Nx; ++i)
			p[i][j] = p[i][j + 1] + dp;

	for (int i = 0; i < Nx - 1; ++i)
		u[i][Ny - 1] = u0;
}

void solve(void)
{
	bool ok = true;
	while (!ok)
	{
		// Init star value from previous calculation
		for (int i = 0; i < Nx; ++i)
			for (int j = 0; j < Ny; ++j)
				p_star[i][j] = p[i][j];

		for (int i = 0; i < Nx-1; ++i)
			for (int j = 0; j < Ny; ++j)
				u_star[i][j] = p[i][j];

		for (int i = 0; i < Nx; ++i)
			for (int j = 0; j < Ny-1; ++j)
				v_star[i][j] = p[i][j];

		// Calculate new star value
		for(int j = 1; j < Ny-1; ++j)
			for (int i = 1; i < Nx-2; ++i)
			{
				const double v_a = relaxation(v[i][j], v[i+1][j], 0.5);
				const double v_b = relaxation(v[i][j-1], v[i+1][j-1], 0.5);
				const double dru2dx = 0.5*(rho * pow(u[i+1][j], 2) - rho * pow(u[i-1][j], 2)) / dx;
				const double druvdy = 0.5*(rho * u[i][j+1] * v_a - rho * u[i][j-1] * v_b) / dy;
				const double dduddx = (u[i+1][j] - 2 * u[i][j] + u[i-1][j]) / dx2;
				const double dduddy = (u[i][j+1] - 2 * u[i][j] + u[i][j-1]) / dy2;
				const double A_star = -(dru2dx + druvdy) + mu * (dduddx + dduddy);
				u_star[i][j] += dt * (A_star - (p[i+1][j] - p[i][j]) / dx) / rho;
			}

		for(int i = 1; i < Nx-1; ++i)
			for(int j = 1; j < Ny-2; ++j)
			{
				const double u_c = relaxation(u[i-1][j], u[i-1][j+1], 0.5);
				const double u_d = relaxation(u[i][j], u[i][j+1], 0.5);
				const double drvudx = 0.5*(rho * v[i+1][j] * u_d - rho * v[i-1][j] * u_c) / dx;
				const double drv2dy = 0.5*(rho * pow(v[i][j+1], 2) - rho * pow(v[i][j-1], 2)) / dy;
				const double ddvddx = (v[i+1][j] - 2 * v[i][j] + v[i-1][j]) / dx2;
				const double ddvddy = (v[i][j+1] - 2 * v[i][j] + v[i][j-1]) / dy2;
				const double B_star = -(drvudx + drv2dy) + mu * (ddvddx + ddvddy);
				v_star[i][j] += dt * (B_star - (p[i][j+1] - p[i][j]) / dy) / rho;
			}

		ok = true;
	}
}

void output(void)
{
	ofstream flow("Cavity.dat");

	flow << "TITLE = \"2D Lid-Driven Cavity Flow\"" << endl;
	flow << "VARIABLES = \"X\", \"Y\", \"DENSITY\", \"U\", \"V\", \"P\"" << endl;
	flow << "ZONE  I = " << Nx << ", J = " << Ny << ", 	F=POINT" << endl;
	
	for (int j = 0; j < Ny; ++j)
		for (int i = 0; i < Nx; ++i)
		{
			flow << setw(16) << scientific << setprecision(6) << x[i];
			flow << setw(16) << scientific << setprecision(6) << y[j];
			flow << setw(16) << scientific << setprecision(6) << rho;

			if(i == 0)
				flow << setw(16) << scientific << setprecision(6) << u[i][j];
			else if(i == Nx-1)
				flow << setw(16) << scientific << setprecision(6) << u[i-1][j];
			else
				flow << setw(16) << scientific << setprecision(6) << relaxation(u[i-1][j], u[i][j], 0.5);

			if (j == 0)
				flow << setw(16) << scientific << setprecision(6) << v[i][j];
			else if (j == Ny - 1)
				flow << setw(16) << scientific << setprecision(6) << v[i][j - 1];
			else
				flow << setw(16) << scientific << setprecision(6) << relaxation(v[i][j - 1], v[i][j], 0.5);
			
			flow << setw(16) << scientific << setprecision(6) << p[i][j];
			flow << endl;
		}

	flow.close();
}

int main(int argc, char *argv[])
{
	setup_grid();
	initialize_flowfield();
	solve();
	output();

	return 0;
}
