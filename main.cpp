#include <iostream>
#include <cmath>
#include <string>
#include <tgmath.h>
#include <fstream>
#include <vector>
#include <gsl/gsl_sf_dawson.h>
//#include <boost/chrono.hpp>
//#include <Eigen/Dense>


using namespace std;
//using namespace Eigen;

double pi = 3.14159;
double hbar = 6.626 * pow(10, -34) / (2 * pi);
//double epsilon_r = 11.68;
double e = 1.602 * pow(10.0, -19);
//double k_const = 8.9876 * pow(10.0, 9);
double m = 9.11 * pow(10, -31);
double m_perp = .198 * m;
double m_z = .92 * m;
double omega0 = .008 * e / hbar; //8 meV
double E_VS = .33*e/1000; // valley splitting in meV, for 21 use .75 meV, for 31 use .33 meV
double g = 1.998;
double r = 1.1 * pow(10.0, -9); //dipole size, 1.1 nm
double rho = 2200; //silicon mass density 2200 kg/m^3 in SiO2 (Peihao's paper)
double v_t = 3750; //3750, 5420
double v_l = 5900; //5900, 9330
double T = .15;
double log_e = 2.718;
double z0 = .11*pow(10.0,-9);

double Rashba = 45;
double Dressel = 0;

double xi_d = 5.0*e; //dilation deformation 5 eV
double xi_u = 8.77*e; //uniaxial sheer deformation 8.77 eV



double Zeeman_Energy(double B) {
    double ub = e * hbar / (2 * m);
    //cout << g * ub * B / hbar << endl;
    return g * ub * B; 
}

double sin_gamma_sq(double B) {
    double e3 = (E_VS - Zeeman_Energy(B));
    double delta23 = r * m_perp * (E_VS) * (Dressel + Rashba) / (sqrt(2) * hbar);
    double e3_tilde = sqrt(e3 * e3 + delta23 * delta23);
    double gamma = atan(delta23 / e3_tilde);
    return (sqrt(e3*e3 + delta23*delta23 ) + e3) / (2 * sqrt(e3*e3+ delta23*delta23));

}

double cos_gamma_sq(double B) {
    double e3 = (E_VS - Zeeman_Energy(B));
    double delta23 = r * m_perp * (E_VS) * (Dressel + Rashba) / (sqrt(2) * hbar);
    double e3_tilde = sqrt(e3 * e3 + delta23 * delta23);
    double gamma = atan(delta23 / e3_tilde);
    return (sqrt(e3 * e3 + delta23 * delta23) - e3) / (2 * sqrt(e3 * e3 + delta23 * delta23));

}

double Single_Electron_Energy(double B, int spin, int valley) {
    //double omega_c = e * B / m_perp;
    //double Omega = sqrt(.25 * omega_c * omega_c + omega0 * omega0);
    //double E = hbar * Omega * 1000 / e;
    double E = valley * E_VS + spin * Zeeman_Energy(B);
    return E;
}
double lambdasq(double B) {
    return hbar / sqrt(m_perp * m_perp * omega0 * omega0 + e * e * B * B / 4.0);
}
// order: spin1, valley1 are for lower energy state
          //spin2, valley2 are for upper energy state
double pure_transition(double B, int spin1, int valley1, int spin2, int valley2) {
    double Energy_diff = Single_Electron_Energy(B, spin2, valley2) - Single_Electron_Energy(B, spin1, valley1);
    double pre = pow(Energy_diff, 5)*r*r / (8.0 * pi * pi * rho * pow(hbar, 6)); //hbar in meV: 6.582 * pow(10, -13)
    double omega_z = Zeeman_Energy(B) / hbar;

    double pre_l = pre / pow(v_l, 7);
    double integrals_l = ((2.0/3.0) * xi_d * xi_d + (4.0 / 15.0) * xi_d * xi_u + (2.0 / 35.0) * xi_u * xi_u);
    double integrals_l_z = ((2.0 / 3.0) * xi_d * xi_d + (4.0 / 5.0) * xi_d * xi_u + (2.0 / 7.0) * xi_u * xi_u);
    double rate_l = pre_l * (integrals_l);
   
    double pre_t = pre / pow(v_t, 7);
    double integrals_t = (xi_u * xi_u * (8.0 / 105.0));
    double integrals_t_z = (xi_u * xi_u * (4.0 / 35.0));
    double rate_t = pre_t * (integrals_t);

    return (rate_l + rate_t) / tanh(hbar * omega_z / (2 * 1.38 * pow(10, -23) * T));

}

double erfi(double x) {
    return (2.0 / sqrt(pi)) * exp(x * x) * gsl_sf_dawson(x);
}

double intravalley_SO(double B) {
    double omega_z = Zeeman_Energy(B) / hbar;
    double pre = pow(Dressel + Rashba, 2.0) * pow(omega_z, 7.0) * pow(omega0, -4.0) / (4.0 * pi * rho * hbar);
    double pre_l = pre / pow(v_l, 7.0);
    double integrals_l = ((2.0 / 3.0) * xi_d * xi_d + (4.0 / 15.0) * xi_d * xi_u + (2.0 / 35.0) * xi_u * xi_u);
    //double integrals_l_z = ((2.0 / 3.0) * xi_d * xi_d + (4.0 / 5.0) * xi_d * xi_u + (2.0 / 7.0) * xi_u * xi_u);
    double rate_l = pre_l * (integrals_l);

    double pre_t = pre / pow(v_t, 7);
    double integrals_t = (xi_u * xi_u * (8.0 / 105.0));
    //double integrals_t_z = (xi_u * xi_u * (4.0 / 35.0));
    double rate_t = pre_t * (integrals_t);
    double expl = exp(-lambdasq(B) * omega_z * omega_z / (2 * v_l * v_l));
    double expt = exp(-lambdasq(B) * omega_z * omega_z / (2 * v_t * v_t));
    return (rate_l*expl + rate_t*expt) / tanh(hbar * omega_z / (2 * 1.38 * pow(10, -23) * T));
}   

double A_int(int n, double delta, double v) {
    double q = delta / (hbar * v);
    double long a = .25 * (r * r - z0 * z0) * q * q;
    if (n == 0) {
        return 0.5*sqrt(pi) * erfi(sqrt(a)) / sqrt(a);
        //return -sqrt(pi) * erfi(sqrt(a)) / sqrt(a);
    }
    if (n == 2) {
        return 0.5*( - 1.0 * sqrt(pi) * erfi(sqrt(a)) / (2 * pow(a, 3.0 / 2.0)) + exp(a) / a);
        //return sqrt(pi) * erfi(sqrt(a)) / (2 * pow(a, 3.0 / 2.0)) - exp(a) / a;
    }
    if (n == 4) {
        return 0.5*(3.0 * sqrt(pi) * erfi(sqrt(a))/(4*pow(a, 5.0/2.0)) + exp(a)*(2*a-3.0)/(2*a*a));
        //return -3.0 * sqrt(pi) * erfi(sqrt(a)) / (4 * pow(a, 5.0 / 2.0)) - exp(a) * (2 * a - 3.0) / (2 * a * a);
    }
    if (n == 6) {
        return 0.5*((4.0 * a * a - 10.0 * a + 15.0) * exp(a) / (4 * pow(a, 3.0)) - (15.0 * sqrt(pi) * erfi(sqrt(a))) / (8 * pow(a, 7.0/2.0)));
        //return -(4.0 * a * a - 10.0 * a + 15.0) * exp(a) / (4.0 * pow(a, 3.0)) + (15.0 * sqrt(pi) * erfi(sqrt(a))) / (8 * pow(a, 7.0 / 2.0));
    }
    else return 0;
}
double bottleneck_SO(double B) {
    double delta = Zeeman_Energy(B);
    double omega_z = delta / (hbar);
    double pre = pow(Dressel + Rashba, 2) * pow(omega_z, 7) * pow(omega0, -4) / (4.0 * pi * rho * hbar);
    double pre_l = pre*exp(-lambdasq(B) * omega_z * omega_z / (2 * v_l * v_l)) / pow(v_l, 7.0);
    double A_term_l = xi_d * xi_d * (A_int(0, delta, v_l) - A_int(2, delta, v_l))
        + 2.0 * xi_u * xi_d * (A_int(2, delta, v_l) - A_int(4, delta, v_l))
        + xi_u * xi_u * (A_int(4, delta, v_l) - A_int(6, delta, v_l));
    double pre_t = pre*exp(-lambdasq(B) * omega_z * omega_z / (2 * v_t * v_t)) / pow(v_t, 7);
    double A_term_t = xi_u * xi_u * (A_int(2, delta, v_t) - 2.0 * A_int(4, delta, v_t) + A_int(6, delta, v_t));
    //cout << exp(-lambdasq(B) * omega_z * omega_z / (2 * v_l * v_l)) << endl;
    return (pre_l*A_term_l + pre_t*A_term_t) / tanh(hbar * omega_z / (2 * 1.38 * pow(10, -23) * T));
}

int main(){
    double data_points = 120;
    double beginning = 1.0;
    double ending = 10.0;
    int spin1 = -1.0;
    int valley1 = -1.0;
    int spin2 = 1.0;
    int valley2 = -1.0;
    int spin3 = -1.0;
    int valley3 = 1.0;

    ofstream myfile;
    myfile.open("relaxation.txt");
    for (double i = beginning*data_points; i <= ending*data_points; i++) {
        double B = i/data_points;
        myfile << B << " ";
        if (E_VS >= Zeeman_Energy(B)) {
            myfile << pure_transition(B, spin1, valley1, spin2, valley2) * cos_gamma_sq(B) + bottleneck_SO(B) << " ";
            myfile << pure_transition(B, spin1, valley1, spin2, valley2) * cos_gamma_sq(B) << " ";
            myfile << bottleneck_SO(B) << endl;
            
        }
        else {
            myfile << pure_transition(B, spin1, valley1, spin3, valley3) * sin_gamma_sq(B) + bottleneck_SO(B) << " ";
            myfile << pure_transition(B, spin1, valley1, spin3, valley3) * sin_gamma_sq(B) << " ";
            myfile << bottleneck_SO(B) << endl;
     
        }
    }
    myfile.close();
    return 0;
}