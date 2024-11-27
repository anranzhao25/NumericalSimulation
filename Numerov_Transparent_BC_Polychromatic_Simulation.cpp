// Numerov_Transparent_BC_Simulation.cpp : This file contains the 'main' function. Program execution begins and ends there.


#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <string>
#include <complex>
#include <cmath>
#include <fstream>
#include "Header.h"
#include <boost/math/special_functions/legendre.hpp>
#include <sys/stat.h>
#include <sys/types.h>

using namespace std;
using namespace std::complex_literals;
using namespace boost::math;

//Physical Parameters
double hbar = 1.054571817e-34; // Reduced Planck's Constant
double mn = 1.67492749804e-27; // Neutron Mass in kg

// Mirror Parameter
double R; // Radius of Curvature of Mirror in m
double Theta; // Angular Size in Degrees
double a; //Roughness Parameter in meters
double U0;  //Re Part of Optical Potential of MgF2 
double U1;  //Im Part of Optical Potential of MgF2 

//Units
double A = 1e-10; //Angstrom in m
double nm = 1e-9; //nm in m
double um = 1e-6; //um in m
double mm = 1e-3; //mm in m
double rad = 2 * M_PI / 360; //Degrees to Radians
double deg = 1 / rad; //Radians to Degrees
double eV = 1.602176634e-19; //eV in Joules
double neV = eV * 1e-9; //eV in Joules
double barn = 1e-28; //Barn in m^2

vector<complex<double>> gaussian_wavepacket(const vector<double>& x, double x0, double sigma0, double p0) {
    double A = pow(2 * M_PI * pow(sigma0, 2), -0.25);
    vector<complex<double>> psi(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        if (i == 0 || i == 1 || i == x.size() - 2 || i == x.size() - 1) {
            psi[i] = 0;
            //psi[i] = A * exp(1i * p0 * x[i] - pow((x[i] - x0) / (2 * sigma0), 2));
        }
        else {
            psi[i] = A * exp(1i * p0 * x[i] - pow((x[i] - x0) / (2 * sigma0), 2));
        }
    }
    return psi;
}

vector<double> potential_step(const vector<double>& x, double x0, double V0) {
    vector<double> V(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        if (x[i] < x0) {
            V[i] = V0;
        }
        else {
            V[i] = 0.0;
        }
    }
    return V;
}

vector<double> diffuse_potential_step(const vector<double>& x, double x0, double V0, double a) {
    vector<double> V(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        V[i] = V0 / (1 + exp(x[i] / a));
    }
    return V;
}

vector<double> linear_potential(const vector<double>& x, double F0) {
    vector<double> V(x.size());
    for (size_t j = 0; j < x.size(); ++j) {
        V[j] = F0 * x[j];
    }
    return V;
}

vector<complex<double>> triangular_potential(const vector<double>& x, double x0, complex<double> V0, double a, double F0) {
    vector<complex<double>> V(x.size());
    
    if (abs(a) > 1e-6) {
        for (size_t i = 0; i < x.size(); ++i) {
            V[i] = V0 / (1.0 + exp(x[i] / a)) + F0 * x[i];
        }
    }

    else {
        for (size_t i = 0; i < x.size(); ++i) {
            if (x[i] <= x0) {
                V[i] = V0 + F0 * x[i];
            }
            else {
                V[i] = 0.0 + F0 * x[i];
            }
        }
    }
    return V;
}

vector<double> CP_potential(const vector<double>& x, double x0, double C1, double C2) {
    vector<double> V(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        V[i] = -C1 / pow(x[i] - x0, 4);
    }
    return V;
}



vector<double> CP_triangular_potential(const vector<double>& x, double x0, double C1, double C2, double F0) {
    vector<double> V(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        V[i] = -C1 / pow(x[i] - x0, 4) + F0 * x[i];
    }
    return V;
}


vector<double> rectangular_potential_barrier(const vector<double>& x, double x1, double x2, double V0) {
    vector<double> V(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        if (x[i] < x1) {
            V[i] = 0;
        }
        else if (x[i] <= x2 && x[i] >= x1) {
            V[i] = V0;
        }
        else
            V[i] = 0;
    }
    return V;
}



double Norm(double dx, vector<complex<double>>& psi) {
    // Trapezoidal Integration
    double N;
    N = 0;
    for (size_t i = 0; i < psi.size(); ++i) {
        if (i == 0 || i == psi.size() - 1) {
            N += pow(abs(psi[i]), 2) * dx / 2.0;
        }
        else {
            N += pow(abs(psi[i]), 2) * dx;
        }
    }
    return N;
}

complex<double> calculate_mu(const complex<double>& a) {
    return (1.0 - pow(abs(a), 2)) / abs(1.0 - pow(a, 2));
}

complex<double> calculate_q_boundary(int index, const vector<complex<double>>& psi_history, const complex<double>& d, const complex<double>& a, const complex<double>& alpha, int n, double dx, double dt) {

    double lamb = 2 * pow(dx, 2) / dt;
    complex<double> c = 1.0 - 1.0i * lamb / (6.0 * d);
    double phi = arg((pow(a, 2) - 1.0) / c);
    complex<double> mu = calculate_mu(a);
    complex<double> l_nk;

    complex<double> q = conj(d) * (conj(a) - alpha) * psi_history[n];

    for (int k = 1; k <= n; ++k) {
        if ((n - k + 1) - 2 < 0) {
            l_nk = exp(complex<double>(0, -(n - k + 1) * phi)) * (legendre_p((n - k + 1), mu.real()) - 0 * legendre_p((n - k + 1) - 2, mu.real())) / (2.0 * (n - k + 1) - 1.0);
        }
        else {
            l_nk = exp(complex<double>(0, -(n - k + 1) * phi)) * (legendre_p((n - k + 1), mu.real()) - legendre_p((n - k + 1) - 2, mu.real())) / (2.0 * (n - k + 1) - 1.0);
        }

        if (n == 0) {
            cout << "l_nk(n = " << n << ") = " << l_nk << "\n";
        }

        q += d * (a - alpha) * l_nk * psi_history[k];
    }

    return q;
}

vector<complex<double>> crank_nicholson_step(const vector<complex<double>>& psi, double dx, int N_x, double dt, int t, vector<complex<double>>& psi_left, vector<complex<double>>& psi_right, const vector<complex<double>>& V,
    vector<complex<double>> a, vector<complex<double>> d, vector<complex<double>> e, vector<complex<double>> g) {
    vector<complex<double>> q(N_x, 0.0), w(N_x, 0.0), f(N_x, 0.0), psi_next(N_x, 0.0);
    // Find q_j given e_j 
    for (int j = 0; j < N_x; ++j) {
        f[j] = 4i * psi[j] / dt;
        if (j == 0) {
            q[j] = calculate_q_boundary(j, psi_left, d[j], a[j], e[j], t, dx, dt);
        }
        else if (j == N_x - 1) {
            q[j] = calculate_q_boundary(j, psi_right, d[j], a[j], e[j], t, dx, dt);
        }
        else {
            q[j] = q[j - 1] / e[j - 1] + pow(dx, 2) * f[j] / d[j];
        }
    }
    // Find w_j given w_J
    for (int j = N_x - 1; j >= 0; --j) {
        if (j == N_x - 1) {
            w[j] = (q[j - 1] + q[j] * e[j - 1]) / (1.0 - e[j] * e[j - 1]);
        }
        else {
            w[j] = (w[j + 1] - q[j]) / e[j];
        }
        psi_next[j] = psi[j] * (1i * pow(dx, 2) / (3.0 * dt) / d[j] - 1.0) + w[j] / d[j];
    }

    return psi_next;
}

vector<complex<double>> simulate_crank_nicholson(vector<complex<double>> psi, const vector<double>& x, double dx, int N_x, double dt, int N_t, const vector<complex<double>>& V) {
    vector<complex<double>> psi_left;
    vector<complex<double>> psi_right;
    vector<complex<double>> e(N_x, 0.0), d(N_x, 0.0), g(N_x, 0.0), a(N_x, 0.0);

    for (int j = 0; j < N_x; ++j) {
        g[j] = V[j] - 2i / dt;
        d[j] = 1.0 - pow(dx, 2) * g[j] / 12.0;
        a[j] = 1.0 + pow(dx, 2) * g[j] / (2.0 * d[j]);
        //  Calculate e vector once, independent of time 
        if (j == 0) {
            e[j] = a[j] + sqrt(pow(a[j], 2) - 1.0);
            if (abs(e[j]) < 1.0) {
                e[j] = a[j] - sqrt(pow(a[j], 2) - 1.0);
            }
        }
        else if (j == N_x - 1) {
            e[j] = a[j] + sqrt(pow(a[j], 2) - 1.0);
            if (abs(e[j]) < 1.0) {
                e[j] = a[j] - sqrt(pow(a[j], 2) - 1.0);
            }
        }
        else {
            e[j] = 2.0 + pow(dx, 2) * g[j] / d[j] - 1.0 / e[j - 1];
            //cout << "e[" << j << "] = " << e[j] << "\n";
        }
    }

    psi_left.push_back(psi[0]);
    psi_right.push_back(psi[N_x - 1]);

    cout << 0.0 << " % " << "\n";
    for (int i = 0; i < N_t; ++i) {
        if ((i + 1) % int(N_t / 10) == 0) {
            double per = round(100 * (double(i + 1) / double(N_t)));
            cout << per << " % " << "\n";
        }
        psi_left.push_back(psi[0]);
        psi_right.push_back(psi[N_x - 1]);
        psi = crank_nicholson_step(psi, dx, N_x, dt, i, psi_left, psi_right, V, a, d, e, g);
    }

    return psi;
}


void simulate_interference(double dx, double dt, double pos, double sigma, double k_perp, complex<double> U, double a, double R, double Theta, const vector<double>& lamb, string data_directory, string par_name, string data_name) {
   
    vector<complex<double>> psi_left;
    vector<complex<double>> psi_right;
    vector<complex<double>> V;
    vector<complex<double>> psi0;
    vector<complex<double>> psi;
    double k;
    double l0;
    double e0;
    double t0;
    double T;

    double xi;
    double xf;
    int N_x;

    double ti;
    double tf;
    int N_t;

    int N_lamb = size(lamb);
    
    /////////////////////////////////////////////////
    //Create Parameter File
    ////////////////////////////////////////////////

    string par_filename = data_directory + par_name;
    ofstream par_file(par_filename);

    par_file << "/* Parameters used in simulating the neutron whispering gallery: \n";
    par_file << "/* dx, dt used in simulations: \n";
    par_file << dx << "\t" << dt << "\n";
    par_file << "/* Woods-Saxon Triangular Potential: U/(1 + exp(x/a)) + mv^2/R x \n";
    par_file << "/* Mirror Parameters: Optical Potential ReU (neV), ImU (neV), Radius R (mm), Angular Size Theta (Degrees), Roughness a0 (A) \n";
    par_file << U0 << "\t" << U1 << "\t" << R << "\t" << Theta * deg << "\t" << a << "\n";
    par_file << "/* Gaussian Initial Wavefunction: x0 (nm), sigma (nm), k perpindicular (1/nm) \n";
    par_file << pos / nm << "\t" << sigma / nm << "\t" << k_perp * nm << "\n";
    par_file << "/* Wavelenghth vector: Lambda (A) \n";
    
    par_file << lamb[0];
    for (int n = 1; n < N_lamb; ++n) {
        par_file << "\t" << lamb[n];
    }
    par_file << "\n";

    par_file << "/* Position matrix: x[lambda_i] (dimensionless) \n";
    par_file << "/* The ith x row corresponds to the ith lambda \n";
    par_file << "/* Each x row has a different size, maintaining dx spacing \n";

    /////////////////////////////////////////////////
    //Create Data File
    ////////////////////////////////////////////////

    string data_filename = data_directory + data_name;
    ofstream data_file(data_filename);

    data_file << "/* Neutron whispering gallery simulation for a polychromatic beam\n";
    data_file << "/* Psi[Theta,lambda]\n";
    data_file << "/* The wave function was propagated to Theta using the parameters described in the param file\n";
    data_file << "/* The ith Psi row corresponds to the ith lambda \n";
    data_file << "/* Each Psi row has a different size, maintaining dx spacing \n";

    /////////////////////////////////////////////////
    //Calculate Polychromatic Interference Pattern
    ////////////////////////////////////////////////
    for (int n = 0; n < N_lamb; ++n) {
        cout << "###############" << " Lambda = " << lamb[n] << " A " << "###############" << "\n";

        /////////////////////////////////////////////////
        //Define Wavelength Dependent Parameter
        ////////////////////////////////////////////////
        k = 2 * M_PI / lamb[n] / A;
        // Scales and Dimensionless Quantities
        l0 = pow(R / pow(k, 2.0) / 2, 1 / 3.0); //WG Length Scale in m
        //l0 = pow(27, 1/3.0);
        e0 = pow(hbar, 2) / (2 * mn * pow(l0, 2)); // WG Energy Scale in J
        t0 = hbar / e0; // WG Time Scale in s
        T = R * Theta / (hbar * k / mn); // Time to propagate along the mirrro

        complex<double> u = U / e0; //dimensionless Optical Potential

        xi = -10;
        xf = u.real() * 2;
        N_x = int((xf - xi) / dx) + 1;

        ti = 0;
        tf = T / t0;
        N_t = int((tf - ti) / dt) + 1;

        vector<double> x(N_x);
        for (int i = 0; i < N_x; ++i) {
            x[i] = xi + i * dx;
        }

        V = triangular_potential(x, 0.0, U / e0, a/l0, 1);
        psi0 = gaussian_wavepacket(x, pos / l0, sigma / l0, k_perp * l0);
        psi = simulate_crank_nicholson(psi0, x, dx, N_x, dt, N_t, V);

        for (int j = 1; j < N_x; ++j) {
            par_file << "\t" << x[j];
            data_file << "\t" << psi[j].real() << "\t" << psi[j].imag();
        }

        par_file << "\n";
        data_file << "\n";
    }

    par_file.close();
    data_file.close();
}


int main() {
    // Simulation Resolution: 
    double dx = 0.05;
    double dt = 0.05;

    // Wavelength Vector:
    int N_Lamb = 10;
    double Lambi = 2; 
    double Lambf = 10;
    double dLamb = (Lambf - Lambi) / (N_Lamb - 1);


    //Initial Wave Function Parameters (Gaussian)
    double pos = 0.0 * nm; // Peak Position, prefactor in nm
    double sig = 1000 * nm; // Packed Width, prefactor in nm
    double k_perp = 0 / nm; // Wave Vector Perpendicular to Surface, prefactor in 1/nm

    //Beam Spectra
    vector<double> Lamb(N_Lamb);
    for (int i = 0; i < N_Lamb; ++i) {
        Lamb[i] = Lambi + i * dLamb;
    }

    //Mirror Properties
    R = 30.008 * mm; // Mirror Radius, prefactor in mm
    a = 10 * A; //Mirror Roughness, prefactor in Angstrom
    Theta = 39.742 * rad; // Mirror Angular Size, prefactor in degrees
    U0 = 2.117850587471434e-26;
    U1 = -2.901520570435523e-32;
    complex<double> U = U0 + U1*1i;

    // Simuation Location:
    string folder = "/Polychromatic_Rough_Step";
    string base_directory = "/Users/anranzhao/Downloads/NumSimData";
    string data_directory = base_directory + folder;

    // Create the directory if it doesn't exist
    mkdir(base_directory.c_str(), 0777);       // Create base directory
    mkdir(data_directory.c_str(), 0777);      // Create subdirectory

    string par_file_name = "/sim_par.dat";
    string sim_file_name = "/sim_data.dat";

    simulate_interference(dx, dt, pos, sig, k_perp, U, a, R, Theta, Lamb, data_directory, par_file_name, sim_file_name);
    
    cout << "FINISHED" << "\n";
    return 0;
}
