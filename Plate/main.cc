#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <limits>

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

	~Array2D() = default;

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
const double R = 287.0; // J/(Kg*K)
const double Cv = 718.0; // J/(Kg*K)
const double G0 = 1.4;
const double Cp = G0 * Cv; // J/(Kg*K)

// Grid params
const size_t IMIN = 1, IMAX = 70;
const size_t JMIN = 1, JMAX = 70;
const double dx = LHORI / (IMAX - 1);
const double dx2 = 2.0 * dx, dxdx = pow(dx, 2);
const double delta = 5 * LHORI / std::sqrt(Re);
const double LVERT = 5 * delta;
const double dy = LVERT / (JMAX - 1);
const double dy2 = 2.0 * dy, dydy = pow(dy, 2);
vector<double> x(IMAX, 0.0), y(JMAX, 0.0);

double dt = 1e-5;
double t = 0.0;

size_t iter = 0;
const size_t MAX_ITER = 2000;

const double u_inf = Ma * a;
const double v_inf = 0.0;
const double p_inf = 101325.0;
const double T_inf = 288.16;
const double Tw = 1.0 * T_inf;

// Primitive variables
Array2D rho(IMAX, JMAX, 0.0);
Array2D u(IMAX, JMAX, 0.0);
Array2D v(IMAX, JMAX, 0.0);
Array2D p(IMAX, JMAX, 0.0);
Array2D T(IMAX, JMAX, Tw);
Array2D e(IMAX, JMAX, 0.0); // Internal energy per unit mass

// Physical properties
Array2D mu(IMAX, JMAX, 0.0); // Kg/(m*s)
Array2D k(IMAX, JMAX, 0.0); // J/(s*m*K)
Array2D lambda(IMAX, JMAX, 0.0); // Kg/(m*s)

// Conservative variables
Array2D U1(IMAX, JMAX, 0.0); // rho
Array2D U2(IMAX, JMAX, 0.0); // rho u
Array2D U3(IMAX, JMAX, 0.0); // rho v
Array2D U5(IMAX, JMAX, 0.0); // rho(e+V^2 / 2)

inline double Sutherland(double T)
{
	static const double mu0 = 1.7894e-5; // Kg/(m*s)
	static const double T0 = 288.16; // K

	return mu0 * pow(T / T0, 1.5) * (T0 + 110.0) / (T + 110.0);
}

double TimeStep()
{
    static const double CFL = 0.5;
	static const double nu_coef = max(4.0/3, G0/Pr);

    double ret = numeric_limits<double>::max();

    for(size_t j = JMIN; j <= JMAX; ++j)
        for(size_t i= IMIN; i <= IMAX; ++i)
        {
            const double loc_nu = mu(i, j) / rho(i, j) * nu_coef;
            const double metric0 = 1.0 / dxdx + 1.0 / dydy;
            const double metric1 = sqrt(metric0);
            const double loc_C = sqrt(G0 * p(i, j) / rho(i, j));
            const double loc_dt = 1.0/(abs(u(i, j))/dx + abs(v(i,j)) / dy + loc_C * metric1 + 2 * loc_nu * metric0);
            ret = min(ret, loc_dt);
        }

    return CFL * ret;
}

inline bool at_boundary(size_t i, size_t j)
{
    return i==IMIN || i==IMAX || j==JMIN || j==JMAX;
}

void set_boundary_values(Array2D &rho, Array2D &u, Array2D &v, Array2D &p, Array2D &T, Array2D &e)
{
    // Front tip
    u(IMIN, JMIN) = 0.0;
    v(IMIN, JMIN) = 0.0;
    p(IMIN, JMIN) = p_inf;
    T(IMIN, JMIN) = T_inf;

    // Inlet
    for (int j = JMIN+1; j <= JMAX; ++j)
    {
        u(IMIN, j) = u_inf;
        v(IMIN, j) = 0.0;
        p(IMIN, j) = p_inf;
        T(IMIN, j) = T_inf;
    }

    // Top(Far)
    for (int i = IMIN+1; i <= IMAX; ++i)
    {
        u(i, JMAX) = u_inf;
        v(i, JMAX) = 0.0;
        p(i, JMAX) = p_inf;
        T(i, JMAX) = T_inf;
    }

    // Bottom
    for (int i = IMIN+1; i <= IMAX; ++i)
    {
        u(i, JMIN) = 0.0;
        v(i, JMIN) = 0.0;
        p(i, JMIN) = 2 * p(i, JMIN+1) - p(i, JMIN+2);
        T(i, JMIN) = Tw;
    }

    // Outlet
    for (int j = JMIN+1; j <= JMAX - 1; ++j)
    {
        u(IMAX, j) = 2 * u(IMAX - 1, j) - u(IMAX - 2, j);
        v(IMAX, j) = 2 * v(IMAX - 1, j) - v(IMAX - 2, j);
        p(IMAX, j) = 2 * p(IMAX - 1, j) - p(IMAX - 2, j);
        T(IMAX, j) = 2 * T(IMAX - 1, j) - T(IMAX - 2, j);
    }

    // Derived Variables
    for (int j = JMIN; j <= JMAX; ++j)
        for (int i = IMIN; i <= IMAX; ++i)
            if(at_boundary(i, j))
            {
                rho(i, j) = p(i, j) / (R*T(i, j));
                e(i, j) = Cv * T(i, j);
            }
}

void update_physical_properties(Array2D &T, Array2D &mu, Array2D &k, Array2D &lambda)
{
    for (size_t j = JMIN; j <= JMAX; ++j)
        for (size_t i = IMIN; i <= IMAX; ++i)
        {
            mu(i, j) = Sutherland(T(i, j));
            k(i, j) = mu(i, j) * Cp / Pr;
            lambda(i, j) = -2.0 / 3 * mu(i, j); // Follow Stokes's hypothesis
        }
}

void init()
{
	/********************************** Grid **********************************/
	for (size_t i = 1; i < IMAX; ++i)
		x[i] = x[i - 1] + dx;

	for (size_t j = 1; j < JMAX; ++j)
		y[j] = y[j - 1] + dy;

	/********************************** I.C. **********************************/
	// Inner
	for (int j = JMIN+1; j <= JMAX-1; ++j)
		for (int i = IMIN+1; i <= IMAX-1; ++i)
		{
			u(i, j) = u_inf;
			v(i, j) = v_inf;
			p(i, j) = p_inf;
			T(i, j) = T_inf;
            rho(i, j) = p(i, j) / (R*T(i, j));
            e(i, j) = Cv * T(i, j);
		}

	/********************************** B.C. **********************************/
    set_boundary_values(rho, u, v, p, T, e);

	/************************ Conservative Variables **************************/
	for (int j = JMIN; j <= JMAX; ++j)
		for (int i = IMIN; i <= IMAX; ++i)
		{
			U1(i, j) = rho(i, j);
			U2(i, j) = rho(i, j)*u(i, j);
			U3(i, j) = rho(i, j)*v(i, j);
			const double K = 0.5*(pow(u(i, j), 2) + pow(v(i, j), 2));
			U5(i, j) = rho(i, j)*(e(i, j) + K);
		}

    /************************** Physical Properties ***************************/
    update_physical_properties(T, mu, k, lambda);
}

void MacCormack()
{
	/***************************** Forward Difference *************************/
	cout << "\tForward Difference...\n";
	// X-derivatives
	Array2D E1(IMAX, JMAX, 0.0);
	Array2D E2(IMAX, JMAX, 0.0);
	Array2D E3(IMAX, JMAX, 0.0);
	Array2D E5(IMAX, JMAX, 0.0);
	for (size_t j = JMIN; j <= JMAX; ++j)
		for (size_t i = IMIN; i <= IMAX; ++i)
		{
			// Backward difference within E for x-derivatives
			double dudx = 0.0;
			if (i == 1)
				dudx = (u(i + 1, j) - u(i, j)) / dx;
			else
				dudx = (u(i, j) - u(i - 1, j)) / dx;

			double dvdx = 0.0;
			if (i == 1)
				dvdx = (v(i + 1, j) - v(i, j)) / dx;
			else
				dvdx = (v(i, j) - v(i - 1, j)) / dx;

			double dTdx = 0.0;
			if (i == 1)
				dTdx = (T(i + 1, j) - T(i, j)) / dx;
			else
				dTdx = (T(i, j) - T(i - 1, j)) / dx;

			// Central difference within E for y-derivatives
			double dudy = 0.0;
			if (j == 1)
				dudy = (u(i, j + 1) - u(i, j)) / dy;
			else if (j == JMAX)
				dudy = (u(i, j) - u(i, j - 1)) / dy;
			else
				dudy = (u(i, j + 1) - u(i, j - 1)) / dy2;

			double dvdy = 0.0;
			if (j == 1)
				dvdy = (v(i, j + 1) - v(i, j)) / dy;
			else if (j == JMAX)
				dvdy = (v(i, j) - v(i, j - 1)) / dy;
			else
				dvdy = (v(i, j + 1) - v(i, j - 1)) / dy2;

			// Shear stress
			const double divergence = dudx + dvdy;
			const double tau_xx = lambda(i, j) * divergence + 2 * mu(i, j) * dudx;
			const double tau_xy = mu(i, j) * (dudy + dvdx);

			// Heat conduction
			const double q_x = -k(i, j) * dTdx;

			// Elements of E
			E1(i, j) = rho(i, j)*u(i, j);
			E2(i, j) = rho(i, j)*pow(u(i, j), 2) + p(i, j) - tau_xx;
			E3(i, j) = rho(i, j)*u(i, j)*v(i, j) - tau_xy;
			E5(i, j) = (U5(i, j) + p(i, j)) * u(i, j) - u(i, j) * tau_xx - v(i, j) * tau_xy + q_x;
		}

	// Y-derivatives
	Array2D F1(IMAX, JMAX, 0.0);
	Array2D F2(IMAX, JMAX, 0.0);
	Array2D F3(IMAX, JMAX, 0.0);
	Array2D F5(IMAX, JMAX, 0.0);
	for (size_t j = JMIN; j <= JMAX; ++j)
		for (size_t i = IMIN; i <= IMAX; ++i)
		{
			// Central difference within F for x-derivatives
			double dudx = 0.0;
			if (i == 1)
				dudx = (u(i + 1, j) - u(i, j)) / dx;
			else if (i == IMAX)
				dudx = (u(i, j) - u(i - 1, j)) / dx;
			else
				dudx = (u(i + 1, j) - u(i - 1, j)) / dx2;

			double dvdx = 0.0;
			if (i == 1)
				dvdx = (v(i + 1, j) - v(i, j)) / dx;
			else if (i == IMAX)
				dvdx = (v(i, j) - v(i - 1, j)) / dx;
			else
				dvdx = (v(i + 1, j) - v(i - 1, j)) / dx2;

			// Backward difference within F for y-derivatives
			double dudy = 0.0;
			if (j == 1)
				dudy = (u(i, j + 1) - u(i, j)) / dy;
			else
				dudy = (u(i, j) - u(i, j - 1)) / dy;

			double dvdy = 0.0;
			if (j == 1)
				dvdy = (v(i, j + 1) - v(i, j)) / dy;
			else
				dvdy = (v(i, j) - v(i, j - 1)) / dy;

			double dTdy = 0.0;
			if (j == 1)
				dTdy = (T(i, j + 1) - T(i, j)) / dy;
			else
				dTdy = (T(i, j) - T(i, j - 1)) / dy;

			// Shear stress
			const double divergence = dudx + dvdy;
			const double tau_xy = mu(i, j) * (dudy + dvdx);
			const double tau_yy = lambda(i, j) * divergence + 2 * mu(i, j) * dvdy;

			// Heat conduction
			const double q_y = -k(i, j) * dTdy;

			// Elements of F
			F1(i, j) = rho(i, j)*v(i, j);
			F2(i, j) = rho(i, j)*u(i, j)*v(i, j) - tau_xy;
			F3(i, j) = rho(i, j)*pow(v(i, j), 2) + p(i, j) - tau_yy;
			F5(i, j) = (U5(i, j) + p(i, j)) * v(i, j) - u(i, j) * tau_xy - v(i, j) * tau_yy + q_y;
		}

	// Temporal derivatives
	Array2D dU1dt(IMAX, JMAX, 0.0);
	Array2D dU2dt(IMAX, JMAX, 0.0);
	Array2D dU3dt(IMAX, JMAX, 0.0);
	Array2D dU5dt(IMAX, JMAX, 0.0);
	for (int j = JMIN+1; j <= JMAX-1; ++j)
		for (int i = IMIN+1; i <= IMAX-1; ++i)
		{
			const double dE1dx = (E1(i + 1, j) - E1(i, j)) / dx;
			const double dF1dy = (F1(i, j + 1) - F1(i, j)) / dy;
			dU1dt(i, j) = -(dE1dx + dF1dy);

			const double dE2dx = (E2(i + 1, j) - E2(i, j)) / dx;
			const double dF2dy = (F2(i, j + 1) - F2(i, j)) / dy;
			dU2dt(i, j) = -(dE2dx + dF2dy);

			const double dE3dx = (E3(i + 1, j) - E3(i, j)) / dx;
			const double dF3dy = (F3(i, j + 1) - F3(i, j)) / dy;
			dU3dt(i, j) = -(dE3dx + dF3dy);

			const double dE5dx = (E5(i + 1, j) - E5(i, j)) / dx;
			const double dF5dy = (F5(i, j + 1) - F5(i, j)) / dy;
			dU5dt(i, j) = -(dE5dx + dF5dy);
		}

	/******************************* Prediction *******************************/
	cout << "\tPrediction...\n";
	Array2D U1_bar(IMAX, JMAX, 0.0);
	Array2D U2_bar(IMAX, JMAX, 0.0);
	Array2D U3_bar(IMAX, JMAX, 0.0);
	Array2D U5_bar(IMAX, JMAX, 0.0);
	for (int j = JMIN+1; j <= JMAX-1; ++j)
		for (int i = IMIN+1; i <= IMAX-1; ++i)
		{
			U1_bar(i, j) = U1(i, j) + dU1dt(i, j) * dt;
			U2_bar(i, j) = U2(i, j) + dU2dt(i, j) * dt;
			U3_bar(i, j) = U3(i, j) + dU3dt(i, j) * dt;
			U5_bar(i, j) = U5(i, j) + dU5dt(i, j) * dt;
		}

    Array2D rho_bar(IMAX, JMAX, 0.0);
    Array2D u_bar(IMAX, JMAX, 0.0);
    Array2D v_bar(IMAX, JMAX, 0.0);
    Array2D p_bar(IMAX, JMAX, 0.0);
    Array2D T_bar(IMAX, JMAX, Tw);
    Array2D e_bar(IMAX, JMAX, 0.0);

    // Update values at inner
    for (int j = JMIN+1; j <= JMAX-1; ++j)
        for (int i = IMIN+1; i <= IMAX-1; ++i)
        {
            rho_bar(i, j) = U1_bar(i, j);
            u_bar(i, j) = U2_bar(i, j) / U1_bar(i, j);
            v_bar(i, j) = U3_bar(i, j) / U1_bar(i, j);
            const double K_bar = 0.5*(pow(u_bar(i, j), 2) + pow(v_bar(i, j), 2));
            e_bar(i, j) = U5_bar(i, j) / U1_bar(i, j) - K_bar;
            T_bar(i, j) = e_bar(i, j) / Cv;
            p_bar(i, j) = rho_bar(i, j) * R * T_bar(i, j);
        }

    // Update values at boundary
    set_boundary_values(rho_bar, u_bar, v_bar, p_bar, T_bar, e_bar);

    /********************************* Checking *******************************/
    cout << "\tChecking predicted primitive values..." << endl;
    for(size_t j=JMIN; j<=JMAX; ++j)
        for(size_t i=IMIN; i<=IMAX; ++i)
        {
			string msg("("+to_string(i)+", "+to_string(j)+"):");
			bool flag = false;

			if(rho_bar(i, j) < 0)
			{
				if(!flag)
					flag = true;
				string rho_msg(" rho_bar="+to_string(rho_bar(i,j)));
				msg += rho_msg;
			}
			if(T_bar(i,j) < 0)
			{
				if(!flag)
					flag = true;
				string T_msg(" T_bar="+to_string(T_bar(i,j)));
				msg += T_msg;
            }
			if(p_bar(i, j) < 0)
			{
				if(!flag)
					flag = true;
				string p_msg(" p_bar="+to_string(p_bar(i,j)));
				msg += p_msg;
			}
            
			if(flag)
				throw runtime_error(msg);
        }

	/*************************** Backward Difference **************************/
	cout << "\tBackward Difference..." << endl;
	Array2D mu_bar(IMAX, JMAX, 0.0);
	Array2D k_bar(IMAX, JMAX, 0.0);
	Array2D lambda_bar(IMAX, JMAX, 0.0);
	update_physical_properties(T_bar, mu_bar, k_bar, lambda_bar);

	Array2D E1_bar(IMAX, JMAX, 0.0);
	Array2D E2_bar(IMAX, JMAX, 0.0);
	Array2D E3_bar(IMAX, JMAX, 0.0);
	Array2D E5_bar(IMAX, JMAX, 0.0);
	for (size_t j = JMIN; j <= JMAX; ++j)
		for (size_t i = IMIN; i <= IMAX; ++i)
		{
			// Forward difference within E for x-derivatives
			double dudx = 0.0;
			if (i == IMAX)
				dudx = (u_bar(i, j) - u_bar(i - 1, j)) / dx;
			else
				dudx = (u_bar(i + 1, j) - u_bar(i, j)) / dx;

			double dvdx = 0.0;
			if (i == IMAX)
				dvdx = (v_bar(i, j) - v_bar(i - 1, j)) / dx;
			else
				dvdx = (v_bar(i + 1, j) - v_bar(i, j)) / dx;

			double dTdx = 0.0;
			if (i == IMAX)
				dTdx = (T_bar(i, j) - T_bar(i - 1, j)) / dx;
			else
				dTdx = (T_bar(i + 1, j) - T_bar(i, j)) / dx;

			// Central difference within E for y-derivatives
			double dudy = 0.0;
			if (j == JMIN)
				dudy = (u_bar(i, j + 1) - u_bar(i, j)) / dy;
			else if (j == JMAX)
				dudy = (u_bar(i, j) - u_bar(i, j - 1)) / dy;
			else
				dudy = (u_bar(i, j + 1) - u_bar(i, j - 1)) / dy2;

			double dvdy = 0.0;
			if (j == JMIN)
				dvdy = (v_bar(i, j + 1) - v_bar(i, j)) / dy;
			else if (j == JMAX)
				dvdy = (v_bar(i, j) - v_bar(i, j - 1)) / dy;
			else
				dvdy = (v_bar(i, j + 1) - v_bar(i, j - 1)) / dy2;

			// Shear stress
			const double divergence = dudx + dvdy;
			const double tau_xx = lambda_bar(i, j) * divergence + 2 * mu_bar(i, j) * dudx;
			const double tau_xy = mu_bar(i, j) * (dudy + dvdx);

			// Heat conduction
			const double q_x = -k_bar(i, j) * dTdx;

			// Elements of E
			E1_bar(i, j) = rho_bar(i, j)*u_bar(i, j);
			E2_bar(i, j) = rho_bar(i, j)*pow(u_bar(i, j), 2) + p_bar(i, j) - tau_xx;
			E3_bar(i, j) = rho_bar(i, j)*u_bar(i, j)*v_bar(i, j) - tau_xy;
			E5_bar(i, j) = (U5_bar(i, j) + p_bar(i, j)) * u_bar(i, j) - u_bar(i, j) * tau_xx - v_bar(i, j) * tau_xy + q_x;
		}

	Array2D F1_bar(IMAX, JMAX, 0.0);
	Array2D F2_bar(IMAX, JMAX, 0.0);
	Array2D F3_bar(IMAX, JMAX, 0.0);
	Array2D F5_bar(IMAX, JMAX, 0.0);
	for (size_t j = JMIN; j <= JMAX; ++j)
		for (size_t i = IMIN; i <= IMAX; ++i)
		{
			// Central difference within F for x-derivatives
			double dudx = 0.0;
			if (i == IMIN)
				dudx = (u_bar(i + 1, j) - u_bar(i, j)) / dx;
			else if (i == IMAX)
				dudx = (u_bar(i, j) - u_bar(i - 1, j)) / dx;
			else
				dudx = (u_bar(i + 1, j) - u_bar(i - 1, j)) / dx2;

			double dvdx = 0.0;
			if (i == IMIN)
				dvdx = (v_bar(i + 1, j) - v_bar(i, j)) / dx;
			else if (i == IMAX)
				dvdx = (v_bar(i, j) - v_bar(i - 1, j)) / dx;
			else
				dvdx = (v_bar(i + 1, j) - v_bar(i - 1, j)) / dx2;

			// Forward difference within F for y-derivatives
			double dudy = 0.0;
			if (j == JMAX)
				dudy = (u_bar(i, j) - u_bar(i, j - 1)) / dy;
			else
				dudy = (u_bar(i, j + 1) - u_bar(i, j)) / dy;

			double dvdy = 0.0;
			if (j == JMAX)
				dvdy = (v_bar(i, j) - v_bar(i, j - 1)) / dy;
			else
				dvdy = (v_bar(i, j + 1) - v_bar(i, j)) / dy;

			double dTdy = 0.0;
			if (j == JMAX)
				dTdy = (T_bar(i, j) - T_bar(i, j - 1)) / dy;
			else
				dTdy = (T_bar(i, j + 1) - T_bar(i, j)) / dy;

			// Shear stress
			const double divergence = dudx + dvdy;
			const double tau_xy = mu_bar(i, j) * (dudy + dvdx);
			const double tau_yy = lambda_bar(i, j) * divergence + 2 * mu_bar(i, j) * dvdy;

			// Heat conduction
			const double q_y = -k_bar(i, j) * dTdy;

			// Elements of F
			F1_bar(i, j) = rho_bar(i, j)*v_bar(i, j);
			F2_bar(i, j) = rho_bar(i, j)*u_bar(i, j)*v_bar(i, j) - tau_xy;
			F3_bar(i, j) = rho_bar(i, j)*pow(v_bar(i, j), 2) + p_bar(i, j) - tau_yy;
			F5_bar(i, j) = (U5_bar(i, j) + p_bar(i, j)) * v_bar(i, j) - u_bar(i, j) * tau_xy - v_bar(i, j) * tau_yy + q_y;
		}

	Array2D dU1dt_bar(IMAX, JMAX, 0.0);
	Array2D dU2dt_bar(IMAX, JMAX, 0.0);
	Array2D dU3dt_bar(IMAX, JMAX, 0.0);
	Array2D dU5dt_bar(IMAX, JMAX, 0.0);
	for (int j = JMIN+1; j <= JMAX-1; ++j)
		for (int i = IMIN+1; i <= IMAX-1; ++i)
		{
			const double dE1dx = (E1_bar(i, j) - E1_bar(i - 1, j)) / dx;
			const double dF1dy = (F1_bar(i, j) - F1_bar(i, j - 1)) / dy;
			dU1dt_bar(i, j) = -(dE1dx + dF1dy);

			const double dE2dx = (E2_bar(i, j) - E2_bar(i - 1, j)) / dx;
			const double dF2dy = (F2_bar(i, j) - F2_bar(i, j - 1)) / dy;
			dU2dt_bar(i, j) = -(dE2dx + dF2dy);

			const double dE3dx = (E3_bar(i, j) - E3_bar(i - 1, j)) / dx;
			const double dF3dy = (F3_bar(i, j) - F3_bar(i, j - 1)) / dy;
			dU3dt_bar(i, j) = -(dE3dx + dF3dy);

			const double dE5dx = (E5_bar(i, j) - E5_bar(i - 1, j)) / dx;
			const double dF5dy = (F5_bar(i, j) - F5_bar(i, j - 1)) / dy;
			dU5dt_bar(i, j) = -(dE5dx + dF5dy);
		}

	/********************************* Average ********************************/
	cout << "\tAveraging derivatives..." << endl;
	Array2D dU1dt_av(IMAX, JMAX, 0.0);
	Array2D dU2dt_av(IMAX, JMAX, 0.0);
	Array2D dU3dt_av(IMAX, JMAX, 0.0);
	Array2D dU5dt_av(IMAX, JMAX, 0.0);
	for (int j = JMIN+1; j <= JMAX-1; ++j)
		for (int i = IMIN+1; i <= IMAX-1; ++i)
		{
			dU1dt_av(i, j) = 0.5*(dU1dt(i, j) + dU1dt_bar(i, j));
			dU2dt_av(i, j) = 0.5*(dU2dt(i, j) + dU2dt_bar(i, j));
			dU3dt_av(i, j) = 0.5*(dU3dt(i, j) + dU3dt_bar(i, j));
			dU5dt_av(i, j) = 0.5*(dU5dt(i, j) + dU5dt_bar(i, j));
		}

	/********************************* Update *********************************/
	cout << "\tUpdating both primitive and conservative values..." << endl;
	// Conservative values at inner
	for (int j = JMIN+1; j <= JMAX-1; ++j)
		for (int i = IMIN+1; i <= IMAX-1; ++i)
		{
			U1(i, j) += dU1dt_av(i, j) * dt;
			U2(i, j) += dU2dt_av(i, j) * dt;
			U3(i, j) += dU3dt_av(i, j) * dt;
			U5(i, j) += dU5dt_av(i, j) * dt;
		}

	// Primitive values at inner
	for (int j = JMIN+1; j <= JMAX-1; ++j)
		for (int i = IMIN+1; i <= IMAX-1; ++i)
		{
			rho(i, j) = U1(i, j);
			u(i, j) = U2(i, j) / U1(i, j);
			v(i, j) = U3(i, j) / U1(i, j);
			const double K = 0.5*(pow(u(i, j), 2) + pow(v(i, j), 2));
			e(i, j) = U5(i, j) / U1(i, j) - K;
			T(i, j) = e(i, j) / Cv;
			p(i, j) = rho(i, j) * R * T(i, j);
		}

	// Primitive values at boundary
	set_boundary_values(rho, u, v, p, T, e);

    /********************************* Checking *******************************/
    cout << "\tChecking updated primitive values..." << endl;
    for(size_t j=JMIN; j<=JMAX; ++j)
        for(size_t i=IMIN; i<=IMAX; ++i)
        {
			string msg("("+to_string(i)+", "+to_string(j)+"):");
			bool flag = false;

			if(rho(i, j) < 0)
			{
				if(!flag)
					flag = true;
				string rho_msg(" rho="+to_string(rho(i,j)));
				msg += rho_msg;
			}
			if(T(i,j) < 0)
			{
				if(!flag)
					flag = true;
				string T_msg(" T="+to_string(T(i,j)));
				msg += T_msg;
            }
			if(p(i, j) < 0)
			{
				if(!flag)
					flag = true;
				string p_msg(" p="+to_string(p(i,j)));
				msg += p_msg;
			}
            
			if(flag)
				throw runtime_error(msg);
        }

    // Physical properties
    update_physical_properties(T, mu, k, lambda);
}

void write_tecplot(size_t n)
{
	// Output format params
	static const size_t WIDTH = 16;
	static const size_t DIGITS = 7;

	// Create Tecplot data file.
	ofstream result("flow" + to_string(n) + ".dat");
	if (!result)
		throw runtime_error("Failed to create data file!");

	// Header
	result << R"(TITLE = "t=)" << t << R"(")" << endl;
	result << R"(VARIABLES = "X", "Y", "rho", "U", "V", "P", "T", "U1", "U2", "U3", "U5")" << endl;
	result << "ZONE I=" << IMAX << ", J=" << JMAX << ", F=POINT" << endl;

	// Flow-field data
	for (int j = 1; j <= JMAX; ++j)
		for (int i = 1; i <= IMAX; ++i)
		{
			result << setw(WIDTH) << setprecision(DIGITS) << x[i - 1];
			result << setw(WIDTH) << setprecision(DIGITS) << y[j - 1];
			result << setw(WIDTH) << setprecision(DIGITS) << rho(i, j);
			result << setw(WIDTH) << setprecision(DIGITS) << u(i, j);
			result << setw(WIDTH) << setprecision(DIGITS) << v(i, j);
			result << setw(WIDTH) << setprecision(DIGITS) << p(i, j);
            result << setw(WIDTH) << setprecision(DIGITS) << T(i, j);
            result << setw(WIDTH) << setprecision(DIGITS) << U1(i, j);
            result << setw(WIDTH) << setprecision(DIGITS) << U2(i, j);
            result << setw(WIDTH) << setprecision(DIGITS) << U3(i, j);
            result << setw(WIDTH) << setprecision(DIGITS) << U5(i, j);
			result << endl;
		}

	// Finalize
	result.close();
}

void write_user(size_t n)
{
	// Output format params
	static const size_t WIDTH = 16;
	static const size_t DIGITS = 7;

	// Output file name
	static const string fn("history.txt");

	ofstream fout;
	if (n == 0)
	{
		fout.open(fn, ios::out);
		if (!fout)
			throw runtime_error("Failed to open history file.");

		for (int i = 0; i < IMAX; ++i)
			fout << setw(WIDTH) << setprecision(DIGITS) << x[i];
		fout << endl;
		for (int j = 0; j < JMAX; ++j)
			fout << setw(WIDTH) << setprecision(DIGITS) << y[j];
		fout << endl;
	}
	else
	{
		fout.open(fn, ios::app);
		if (!fout)
			throw runtime_error("Failed to open history file.");
	}

	fout.close();
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
        cout << "Iter" << iter << ":" << endl;

        dt = TimeStep();
        if(dt < 0)
            throw runtime_error("dt="+to_string(dt)+"s");
        cout << "\tt=" << t <<"s, dt=" << dt << "s" << endl;

		MacCormack();

        t += dt;
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
