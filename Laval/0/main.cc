#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstddef>
#include <ctime>

using namespace std;

const size_t WIDTH = 16;
const size_t DIGITS = 6;
const size_t OUTPUT_GAP = 50;
const size_t MAX_STEP = 2000;

const double G0 = 1.4;
const double G1 = 1 / G0;
const double G2 = G0 - 1;
const double R = 8.314e3; // J/(mol*K)

const double P0 = 10 * 101325.0; // Pa
const double T0 = 300.0; // K
const double rho0 = P0 / (R*T0);
const double a0 = sqrt(G0 * R * T0);
const double A_star = 1.0;

double CFL = 0.5;
size_t iter_cnt = 0;
double t = 0.0;

const int N = 31;
const double dx = 0.1;
vector<double> A(N, 0.0), x(N, 0.0);
vector<double> a(N), rho(N), T(N), V(N), P(N), Ma(N);
vector<double> rho_bar(N), T_bar(N), V_bar(N);
vector<double> drhodt(N, 0.0), dVdt(N, 0.0), dTdt(N, 0.0);
vector<double> drhodt_bar(N, 0.0), dVdt_bar(N, 0.0), dTdt_bar(N, 0.0);
vector<double> drhodt_av(N, 0.0), dVdt_av(N, 0.0), dTdt_av(N, 0.0);

void init()
{
    // Grid
    for(int i = 1; i < N; ++i)
        x[i] = x[i-1] + dx;

    // Cross-section area
    for(int i = 0; i < N; ++i)
        A[i] = A_star + 2.2*pow(x[i]-1.5, 2);

    // I.C.
    for(int i = 0; i < N; ++i)
    {
        rho[i] = 1.0 - 0.3146 * x[i];
        T[i] = 1.0 - 0.2314 * x[i];
        V[i] = (0.1 + 1.09*x[i]) * sqrt(T[i]);
        P[i] = rho[i] * T[i];
        Ma[i] = V[i] / sqrt(T[i]);
    }
}

void loop()
{
    // Local Sound Speed
    for(int i = 0; i < N; ++i)
        a[i] = sqrt(T[i]);

    // Determine time-step
    double dt = 1.0;
    for(int i = 0; i < N; ++i)
    {
        double loc_dt = CFL * dx / (a[i] + abs(V[i]));
        if(loc_dt < dt)
            dt = loc_dt;
    }
    cout << "\tdt=" << dt << endl;

    // Forward derivative
    for(int i = 1; i < N-1; ++i)
    {
        const double drho = rho[i+1]-rho[i];
        const double dV = V[i+1]-V[i];
        const double dT = T[i+1] - T[i];
        const double dlnA = log(A[i+1])-log(A[i]);

        drhodt[i] = -rho[i]*dV/dx - rho[i]*V[i]*dlnA/dx - V[i]*drho/dx;
        dVdt[i] = -V[i]*dV/dx - G1*(dT/dx + T[i]/rho[i]*drho/dx);
        dTdt[i] = -V[i]*dT/dx - G2*T[i]*(dV/dx+V[i]*dlnA/dx);
    }

    // Prediction
    // Do NOT miss the bar value at left boundary, which will be used for backward derivatives!!!
    for(int i = 0; i < N; ++i)
    {
        rho_bar[i] = rho[i] + drhodt[i] * dt;
        V_bar[i] = V[i] + dVdt[i] * dt;
        T_bar[i] = T[i] + dTdt[i] * dt;
    }

    // Backward derivative
    for(int i = 1; i < N-1; ++i)
    {
        const double drho = rho_bar[i]-rho_bar[i-1];
        const double dV = V_bar[i]-V_bar[i-1];
        const double dT = T_bar[i] - T_bar[i-1];
        const double dlnA = log(A[i])-log(A[i-1]);

        drhodt_bar[i] = -rho_bar[i]*dV/dx - rho_bar[i]*V_bar[i]*dlnA/dx - V_bar[i]*drho/dx;
        dVdt_bar[i] = -V_bar[i]*dV/dx - G1*(dT/dx + T_bar[i]/rho_bar[i]*drho/dx);
        dTdt_bar[i] = -V_bar[i]*dT/dx - G2*T_bar[i]*(dV/dx+V_bar[i]*dlnA/dx);
    }

    // Correction
    for(int i = 1; i < N-1; ++i)
    {
        drhodt_av[i] = 0.5 * (drhodt[i] + drhodt_bar[i]);
        dVdt_av[i] = 0.5 * (dVdt[i] + dVdt_bar[i]);
        dTdt_av[i] = 0.5 * (dTdt[i] + dTdt_bar[i]);
    }
    for(int i = 1; i < N-1; ++i)
    {
        rho[i] += drhodt_av[i] * dt;
        V[i] += dVdt_av[i] * dt;
        T[i] += dTdt_av[i] * dt;
    }

    // B.C. at inlet
    // Subsonic inlet, one determined from interior
    // Eigenvalue: u-a < 0, u > 0, u+a > 0
    rho[0] = 1.0;
    V[0] = 2 * V[1] - V[2]; 
    T[0] = 1.0;

    // B.C. at outlet
    // Supersonic outlet, all determined from interior
    // Eigenvalue: u-a > 0, u > 0, u+a > 0
    rho[N-1] = 2 * rho[N-2] - rho[N-3];
    V[N-1] = 2 * V[N-2] - V[N-3];
    T[N-1] = 2 * T[N-2] - T[N-3];

    // Mach Number
    for(int i = 0; i < N; ++i)
    {
        Ma[i] = V[i] / sqrt(T[i]);
        P[i] = rho[i] * T[i];
    }
	
	t += dt;
}

double fder1(double fl, double fr)
{
    return 0.5 * (fr - fl) / dx;
}

bool check_convergence()
{
    double rms = 0.0;
    for(int i = 1; i < N-1; ++i)
    {
        double mdot_L = rho[i-1] * V[i-1] * A[i-1];
        double mdot_R = rho[i+1] * V[i+1] * A[i+1];
        rms += pow(fder1(mdot_L, mdot_R), 2);
    }
    rms = sqrt(rms / (N-1));

    cout << "\tRMS=" << rms << endl;

    return rms < 1e-3 || iter_cnt > MAX_STEP;
}

void output()
{
    if(iter_cnt % OUTPUT_GAP != 0)
        return;
    
    stringstream ss;
    ss << "iter" << iter_cnt << "_t=" << t << ".dat";
    ofstream fout(ss.str());
	fout << "Variables = ";
    fout << "x/L";
    fout << setw(WIDTH) << "A/A_star";
    fout << setw(WIDTH) << "rho/rho0";
    fout << setw(WIDTH) << "V/a0";
    fout << setw(WIDTH) << "T/T0";
    fout << setw(WIDTH) << "p/p0";
    fout << setw(WIDTH) << "Ma" << endl;
	fout << "Zone I=" << N << ", F = point" << endl;
    for(int i = 0; i < N; ++i)
    {
        fout << setw(WIDTH) << setprecision(DIGITS) << x[i];
        fout << setw(WIDTH) << setprecision(DIGITS) << A[i];
        fout << setw(WIDTH) << setprecision(DIGITS) << rho[i];
        fout << setw(WIDTH) << setprecision(DIGITS) << V[i];
        fout << setw(WIDTH) << setprecision(DIGITS) << T[i];
        fout << setw(WIDTH) << setprecision(DIGITS) << P[i];
        fout << setw(WIDTH) << setprecision(DIGITS) << Ma[i] << endl;
    }

    fout.close();
}

void solve()
{
    bool ok = false;
    while(!ok)
    {
        cout << "Iter" << ++iter_cnt << ":\n";
        loop();
        output();
        ok = check_convergence();
    }
}

int main(int argc, char *argv[])
{
    init();
    output();
    solve();
    return 0;
}
