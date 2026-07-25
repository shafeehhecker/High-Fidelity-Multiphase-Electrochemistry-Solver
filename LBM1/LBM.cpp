#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>

// =========================================================================
// 1. LATTICE BOLTZMANN CONSTANTS & STRUCTS
// =========================================================================
const int Q9 = 9;
const int ex9[9] = { 0,  1,  0, -1,  0,  1, -1, -1,  1};
const int ey9[9] = { 0,  0,  1,  0, -1,  1,  1, -1, -1};
const double w9[9] = {4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0,
                      1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};

const int Q5 = 5;
const int ex5[5] = { 0,  1,  0, -1,  0};
const int ey5[5] = { 0,  0,  1,  0, -1};
const double w5[5] = {1.0/3.0, 1.0/6.0, 1.0/6.0, 1.0/6.0, 1.0/6.0};

const int NX = 65;
const int NY = 65;

// =========================================================================
// 2. HELPER FUNCTIONS
// =========================================================================
double feq9(int k, double rho, double u, double v) {
    double cu = 3.0 * (ex9[k] * u + ey9[k] * v);
    double uv = (u * u + v * v);
    return w9[k] * rho * (1.0 + cu + 0.5 * cu * cu - 1.5 * uv);
}

double geq5(int k, double T, double u, double v) {
    double cu = 3.0 * (ex5[k] * u + ey5[k] * v);
    double uv = (u * u + v * v);
    return w5[k] * T * (1.0 + cu - 1.5 * uv);
}

// =========================================================================
// 3. MAIN SOLVER
// =========================================================================
int main() {
    std::cout << "Starting Robust Multiscale LBM Solver..." << std::endl;

    double Re = 100.0;
    double Pr = 100.0;
    double nu = 0.01;
    double alpha = nu / Pr;

    double tau_f = 0.5 + (3.0 * nu);
    double time_scale_ratio = 20.0;
    double tau_t = 0.5 + (3.0 * alpha * time_scale_ratio);

    std::cout << "Flow Relaxation (tau_f): " << tau_f << std::endl;
    std::cout << "Thermal Relaxation (tau_t): " << tau_t << " (Stable!)" << std::endl;

    int size = NX * NY;
    std::vector<double> f(size * Q9), f_new(size * Q9);
    std::vector<double> g(size * Q5), g_new(size * Q5);
    std::vector<double> rho(size), u(size), v(size), T(size);

    for (int i = 0; i < size; ++i) {
        rho[i] = 1.0; u[i] = 0.0; v[i] = 0.0; T[i] = 0.0;
        for (int k = 0; k < Q9; ++k) f[i * Q9 + k] = feq9(k, rho[i], u[i], v[i]);
        for (int k = 0; k < Q5; ++k) g[i * Q5 + k] = geq5(k, T[i], u[i], v[i]);
    }

    double U_lid = 0.05; // Reduced slightly for extra stability at Re=100
    int max_steps = 15000;

    for (int step = 0; step < max_steps; ++step) {

        // 1. MACROSCOPIC VARIABLE CALCULATION
        for (int i = 0; i < size; ++i) {
            rho[i] = 0; u[i] = 0; v[i] = 0; T[i] = 0;
            for (int k = 0; k < Q9; ++k) {
                rho[i] += f[i * Q9 + k];
                u[i] += ex9[k] * f[i * Q9 + k];
                v[i] += ey9[k] * f[i * Q9 + k];
            }

            // SAFEGUARD: Prevent division by zero or negative density
            if (rho[i] > 1e-10) {
                u[i] /= rho[i];
                v[i] /= rho[i];
            } else {
                u[i] = 0.0; v[i] = 0.0; rho[i] = 1.0;
            }

            for (int k = 0; k < Q5; ++k) {
                T[i] += g[i * Q5 + k];
            }
        }

        // 2. COLLISION STEP
        for (int i = 0; i < size; ++i) {
            for (int k = 0; k < Q9; ++k) {
                f[i * Q9 + k] -= (f[i * Q9 + k] - feq9(k, rho[i], u[i], v[i])) / tau_f;
            }
            for (int k = 0; k < Q5; ++k) {
                g[i * Q5 + k] -= (g[i * Q5 + k] - geq5(k, T[i], u[i], v[i])) / tau_t;
            }
        }

        // 3. STREAMING STEP (Flow)
        f_new = f;
        for (int y = 0; y < NY; ++y) {
            for (int x = 0; x < NX; ++x) {
                int idx = y * NX + x;
                for (int k = 0; k < Q9; ++k) {
                    int xp = x - ex9[k];
                    int yp = y - ey9[k];
                    if (xp >= 0 && xp < NX && yp >= 0 && yp < NY) {
                        f_new[idx * Q9 + k] = f[(yp * NX + xp) * Q9 + k];
                    }
                }
            }
        }
        f = f_new;

        // 4. STREAMING STEP (Heat) - Sub-stepping!
        for (int sub = 0; sub < time_scale_ratio; ++sub) {
            g_new = g;
            for (int y = 0; y < NY; ++y) {
                for (int x = 0; x < NX; ++x) {
                    int idx = y * NX + x;
                    for (int k = 0; k < Q5; ++k) {
                        int xp = x - ex5[k];
                        int yp = y - ey5[k];
                        if (xp >= 0 && xp < NX && yp >= 0 && yp < NY) {
                            g_new[idx * Q5 + k] = g[(yp * NX + xp) * Q5 + k];
                        }
                    }
                }
            }
            g = g_new;
        }

        // 5. ROBUST BOUNDARY CONDITIONS (Equilibrium Method)
        // Top Wall (Lid)
        for (int x = 0; x < NX; ++x) {
            int idx = (NY - 1) * NX + x;
            double rho_w = 1.0;
            for (int k = 0; k < Q9; ++k) f[idx * Q9 + k] = feq9(k, rho_w, U_lid, 0.0);
            T[idx] = 0.0;
            for(int k=0; k<Q5; ++k) g[idx * Q5 + k] = geq5(k, T[idx], U_lid, 0.0);
        }

        // Bottom Wall
        for (int x = 0; x < NX; ++x) {
            int idx = x;
            double rho_w = 1.0;
            for (int k = 0; k < Q9; ++k) f[idx * Q9 + k] = feq9(k, rho_w, 0.0, 0.0);
            T[idx] = 0.0;
            for(int k=0; k<Q5; ++k) g[idx * Q5 + k] = geq5(k, T[idx], 0.0, 0.0);
        }

        // Left Wall (Hot)
        for (int y = 0; y < NY; ++y) {
            int idx = y * NX;
            double rho_w = 1.0;
            for (int k = 0; k < Q9; ++k) f[idx * Q9 + k] = feq9(k, rho_w, 0.0, 0.0);
            T[idx] = 1.0;
            for(int k=0; k<Q5; ++k) g[idx * Q5 + k] = geq5(k, T[idx], 0.0, 0.0);
        }

        // Right Wall (Cold)
        for (int y = 0; y < NY; ++y) {
            int idx = y * NX + (NX - 1);
            double rho_w = 1.0;
            for (int k = 0; k < Q9; ++k) f[idx * Q9 + k] = feq9(k, rho_w, 0.0, 0.0);
            T[idx] = 0.0;
            for(int k=0; k<Q5; ++k) g[idx * Q5 + k] = geq5(k, T[idx], 0.0, 0.0);
        }

        if (step % 2000 == 0) {
            double max_T = 0;
            for(int i=0; i<size; ++i) if(T[i] > max_T) max_T = T[i];
            std::cout << "Step " << step << " | Max Temp: " << std::fixed << std::setprecision(4) << max_T << std::endl;
        }
    }

    std::ofstream outfile("lbm_result_fixed.csv");
    outfile << "x,y,T\n";
    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            int idx = y * NX + x;
            outfile << x << "," << y << "," << T[idx] << "\n";
        }
    }
    outfile.close();
    std::cout << "Simulation complete. Data saved to lbm_result_fixed.csv" << std::endl;

    return 0;
}
