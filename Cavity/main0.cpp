#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <Eigen/Sparse>

using namespace std;
using namespace Eigen;

const int Nx = 6, Ny = 5;
const double Lx = 0.1, Ly = 0.08; // m
const double xL = -Lx / 2, xR = Lx / 2;
const double yL = -Ly / 2, yR = Ly / 2;
const double dx = Lx / (Nx - 1), dy = Ly / (Ny - 1);
const double dx2 = pow(dx, 2), dy2 = pow(dy, 2);

const double Re = 200.0;
const double rho = 1.0; // kg/m^3
const double p0 = 101325.0; // Pa
const double u0 = 0.1; // m/s
const double mu = rho * u0 * max(Lx, Ly) / Re; // Pa*s

const double alpha_p = 0.8;

const double dt = 0.1;
const double coef_a = 2 * dt * (1.0 / dx2 + 1.0 / dy2);
const double coef_b = -dt / dx2;
const double coef_c = -dt / dy2;

vector<double> x(Nx, xL), y(Ny, yL); // m

vector<vector<double>> p(Nx, vector<double>(Ny, p0)); // Pa
vector<vector<double>> u(Nx - 1, vector<double>(Ny, 0.0)); // m/s
vector<vector<double>> v(Nx, vector<double>(Ny - 1, 0.0)); // m/s

vector<vector<double>> p_star(Nx, vector<double>(Ny, p0)); // Pa
vector<vector<double>> u_star(Nx - 1, vector<double>(Ny, 0.0)); // m/s
vector<vector<double>> v_star(Nx, vector<double>(Ny - 1, 0.0)); // m/s

vector<vector<double>> p_prime(Nx, vector<double>(Ny, 0.0)); // Pa
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
	cout << "Re = " << Re << endl;
	cout << "rho = " << rho << endl;
	cout << "u0 = " << u0 << endl;
	cout << "mu = " << mu << endl;

	for (int i = 0; i < Nx - 1; ++i)
		u[i][Ny - 1] = u0;
}

void output(const string &fn, const vector<vector<double>> &u, const vector<vector<double>> &v, const vector<vector<double>> &p)
{
	ofstream flow(fn);

	flow << "TITLE = \"2D Lid-Driven Cavity Flow\"" << endl;
	flow << "VARIABLES = \"X\", \"Y\", \"DENSITY\", \"U\", \"V\", \"P\"" << endl;
	flow << "ZONE  I = " << Nx << ", J = " << Ny << ", 	F=POINT" << endl;

	for (int j = 0; j < Ny; ++j)
		for (int i = 0; i < Nx; ++i)
		{
			flow << setw(16) << scientific << setprecision(6) << x[i];
			flow << setw(16) << scientific << setprecision(6) << y[j];
			flow << setw(16) << scientific << setprecision(6) << rho;

			if (i == 0 || i == Nx - 1)
			{
				if (j == Ny-1)
					flow << setw(16) << scientific << setprecision(6) << u0;
				else
					flow << setw(16) << scientific << setprecision(6) << 0.0; // u at 2 vertical boundary is 0
			}
			else
				flow << setw(16) << scientific << setprecision(6) << relaxation(u[i - 1][j], u[i][j], 0.5);

			if (j == 0 || j == Ny - 1)
				flow << setw(16) << scientific << setprecision(6) << 0.0; // v at 2 horizontal boundary is 0
			else
				flow << setw(16) << scientific << setprecision(6) << relaxation(v[i][j - 1], v[i][j], 0.5);

			flow << setw(16) << scientific << setprecision(6) << p[i][j];
			flow << endl;
		}

	flow.close();
}

void solve(void)
{
	bool ok = false;
	int iter = 0;
	while (!ok)
	{
		output("Cavity.dat", u, v, p);

		cout << "Iter " << ++iter << ":" << endl;
		cout << "\tInit star value from previous calculation..." << endl;
		for (int i = 0; i < Nx; ++i)
			for (int j = 0; j < Ny; ++j)
				p_star[i][j] = p[i][j];

		for (int i = 0; i < Nx - 1; ++i)
			for (int j = 0; j < Ny; ++j)
				u_star[i][j] = u[i][j];

		for (int i = 0; i < Nx; ++i)
			for (int j = 0; j < Ny - 1; ++j)
				v_star[i][j] = v[i][j];

		cout << "\tCalculate new star value..." << endl;
		double v_a, v_b, dru2dx, druvdy, dduddx, dduddy, A_star, dpdx;
		for (int j = 1; j < Ny - 1; ++j)
			for (int i = 0; i < Nx - 1; ++i)
			{
				v_a = relaxation(v[i][j], v[i + 1][j], 0.5);
				v_b = relaxation(v[i][j - 1], v[i + 1][j - 1], 0.5);
				if (i == 0)
				{
					dru2dx = (rho * pow(u[i + 1][j], 2) + 3 * rho * pow(u[i][j], 2)) / (3 * dx);
					dduddx = (u[i + 1][j] - 3 * u[i][j]) / (0.75 * dx2);
				}
				else if (i == Nx - 2)
				{
					dru2dx = (-rho * pow(u[i - 1][j], 2) - 3 * rho * pow(u[i][j], 2)) / (3 * dx);
					dduddx = (u[i - 1][j] - 3 * u[i][j]) / (0.75 * dx2);
				}
				else
				{
					dru2dx = (rho * pow(u[i + 1][j], 2) - rho * pow(u[i - 1][j], 2)) / (2 * dx);
					dduddx = (u[i + 1][j] - 2 * u[i][j] + u[i - 1][j]) / dx2;
				}
				druvdy = 0.5*(rho * u[i][j + 1] * v_a - rho * u[i][j - 1] * v_b) / dy;
				dduddy = (u[i][j + 1] - 2 * u[i][j] + u[i][j - 1]) / dy2;
				A_star = -(dru2dx + druvdy) + mu * (dduddx + dduddy);
				dpdx = (p[i + 1][j] - p[i][j]) / dx;
				u_star[i][j] += dt * (A_star - dpdx) / rho;
			}

		double u_c, u_d, drvudx, drv2dy, ddvddx, ddvddy, B_star, dpdy;
		for (int i = 1; i < Nx - 1; ++i)
			for (int j = 0; j < Ny - 1; ++j)
			{
				u_c = relaxation(u[i - 1][j], u[i - 1][j + 1], 0.5);
				u_d = relaxation(u[i][j], u[i][j + 1], 0.5);
				if (j == 0)
				{
					drv2dy = (rho * pow(v[i][j + 1], 2) + 3 * rho * pow(v[i][j], 2)) / (3 * dy);
					ddvddy = (v[i][j + 1] - 3 * v[i][j]) / (0.75 * dy2);
				}
				else if (j == Ny - 2)
				{
					drv2dy = (-rho * pow(v[i][j - 1], 2) - 3 * rho * pow(v[i][j], 2)) / (3 * dy);
					ddvddy = (v[i][j - 1] - 3 * v[i][j]) / (0.75 * dy2);
				}
				else
				{
					drv2dy = 0.5*(rho * pow(v[i][j + 1], 2) - rho * pow(v[i][j - 1], 2)) / dy;
					ddvddy = (v[i][j + 1] - 2 * v[i][j] + v[i][j - 1]) / dy2;
				}
				drvudx = 0.5*(rho * v[i + 1][j] * u_d - rho * v[i - 1][j] * u_c) / dx;
				ddvddx = (v[i + 1][j] - 2 * v[i][j] + v[i - 1][j]) / dx2;
				B_star = -(drvudx + drv2dy) + mu * (ddvddx + ddvddy);
				dpdy = (p[i][j + 1] - p[i][j]) / dy;
				v_star[i][j] += dt * (B_star - dpdy) / rho;
			}

		output("Cavity_star.dat", u_star, v_star, p_star);

		cout << "\tSolve the Possion equation..." << endl;
		const int NumOfUnknown = (Nx - 2) * (Ny - 2);

		SparseMatrix<double> A(NumOfUnknown, NumOfUnknown);
		A.setZero();
		VectorXd b(NumOfUnknown);
		b.setZero();
		vector<Triplet<double>> elem;

		for (int i = 1; i < Nx - 1; ++i)
			for (int j = 1; j < Ny - 1; ++j)
			{
				const int idx0 = (j - 1)*(Nx - 2) + (i - 1);
				const int idx1 = (j - 1)*(Nx - 2) + i;
				const int idx2 = j * (Nx - 2) + (i - 1);
				const int idx3 = (j - 1)*(Nx - 2) + (i - 2);
				const int idx4 = (j - 2)*(Nx - 2) + (i - 1);

				double v0 = 0.0;
				double v1 = 0.0;
				double v2 = 0.0;
				double v3 = 0.0;
				double v4 = 0.0;

				v0 += coef_a;
				if (i == Nx - 2)
					v0 += coef_b;
				else
					v1 += coef_b;

				if (j == Ny - 2)
					b[idx0] -= coef_c * 0.0;
				else
					v2 += coef_c;

				if (i == 1)
					v0 += coef_b;
				else
					v3 += coef_b;

				if (j == 1)
					v0 += coef_c;
				else
					v4 += coef_c;

				double drusdx = (rho*u_star[i][j] - rho * u_star[i - 1][j]) / dx;
				double drvsdy = (rho*v_star[i][j] - rho * v_star[i][j - 1]) / dy;
				b[idx0] -= (drusdx + drvsdy);

				elem.push_back(Triplet<double>(idx0, idx0, v0));
				if (v1 != 0.0)
					elem.push_back(Triplet<double>(idx0, idx1, v1));
				if (v2 != 0.0)
					elem.push_back(Triplet<double>(idx0, idx2, v2));
				if (v3 != 0.0)
					elem.push_back(Triplet<double>(idx0, idx3, v3));
				if (v4 != 0.0)
					elem.push_back(Triplet<double>(idx0, idx4, v4));
			}

		A.setFromTriplets(elem.begin(), elem.end());

		// Take Cholesky decomposition of A
		SimplicialCholesky<SparseMatrix<double>> chol(A);
		// Solve
		VectorXd x = chol.solve(b);

		/*
		ofstream fout("A.txt");
		fout << A;
		fout.close();
		fout.open("b.txt");
		fout << b;
		fout.close();
		fout.open("x.txt");
		fout << x;
		fout.close();
		*/

		int cnt = 0;
		for (int i = 1; i < Nx - 1; ++i)
			for (int j = 1; j < Ny - 1; ++j)
				p_prime[i][j] = x[cnt++];

		for (int j = 1; j < Ny - 1; ++j)
			for (int i = 0; i < Nx - 1; ++i)
				u_prime[i][j] = -dt / dx * (p_prime[i + 1][j] - p_prime[i][j]) / rho;

		for (int i = 1; i < Nx - 1; ++i)
			for (int j = 0; j < Ny - 1; ++j)
				v_prime[i][j] = -dt / dy * (p_prime[i][j + 1] - p_prime[i][j]) / rho;

		output("Cavity_prime.dat", u_prime, v_prime, p_prime);

		for (int i = 1; i < Nx - 1; ++i)
			for (int j = 1; j < Ny - 1; ++j)
				p[i][j] += alpha_p * p_prime[i][j];

		for (int j = 1; j < Ny - 1; ++j)
			for (int i = 0; i < Nx - 1; ++i)
				u[i][j] += u_prime[i][j];

		for (int i = 1; i < Nx - 1; ++i)
			for (int j = 0; j < Ny - 1; ++j)
				v[i][j] += v_prime[i][j];

		// Check convergence
		double rsd = 0.0;
		for (int p = 0; p < NumOfUnknown; ++p)
		{
			double tmp = fabs(b[p]);
			if (tmp > rsd)
				rsd = tmp;
		}
		cout << "\trsd=" << rsd << endl;
		ok = rsd < 1e-4 && iter > 100;
	}
}

int main(int argc, char *argv[])
{
	setup_grid();
	initialize_flowfield();
	solve();

	return 0;
}
