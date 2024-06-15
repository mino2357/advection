//
// 1-dim Navier-Stokes solver.
//

#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

// mesh with equal intervals.
constexpr double Lx = 3.0;
constexpr unsigned int Nx = 100;
constexpr double dx = Lx / (Nx - 1);
// time step, and simulation time.
constexpr double dt = 0.1 * (dx / 350.0);
constexpr double time_end = 0.02;

class Fluid {
    private:
        double start_x = -1.0;
        // parameters
        double T = 300.0;
        double R = 8.31;
        double mu = 1.83e-5;
        double time = 0.0;
        std::vector<double> q; // q = rho * u.
        std::vector<double> rho;
    public:
        Fluid();
        void set_initial_conditions();
        void print_q();
        void print_u();
        void print_rho();
        void print_p();
        void integrate_by_Euler();
        void integrate_to_end_time_by_Euler();
};

Fluid::Fluid(){
    q.resize(Nx);
    rho.resize(Nx);
}

void Fluid::set_initial_conditions(){
    for (unsigned int i = 0; i < Nx; i++){
        double x = start_x + i * Lx / Nx;
        if (x < 0.5){
            rho[i] = 40.643802647412755716004813; // 101325.0 / (R * T) = 40.6438026474...
        } else {
            rho[i] = 0.9 * 40.643802647412755716004813;
        }
        q[i] = 0.0;
    }
}

void Fluid::print_q(){
    for (unsigned int i = 0; i < Nx; i++){
        double x = start_x + i * dx;
        std::cout << std::setprecision(14) << x << " " << q[i] << std::endl;
    }
    std::cout << std::endl;
}

void Fluid::print_u(){
    for (unsigned int i = 0; i < Nx; i++){
        double x = start_x + i * dx;
        std::cout << std::setprecision(14) << x << " " << q[i] / rho[i] << std::endl;
    }
    std::cout << std::endl;
}

void Fluid::print_rho(){
    for (unsigned int i = 0; i < Nx; i++){
        double x = start_x + i * dx;
        std::cout << std::setprecision(14) << x << " " << rho[i] << std::endl;
    }
    std::cout << std::endl;
}

void Fluid::print_p(){
    for (unsigned int i = 0; i < Nx; i++){
        double x = start_x + i * dx;
        std::cout << std::setprecision(14) << x << " " << rho[i] * R * T << std::endl;
    }
    std::cout << std::endl;
}

void Fluid::integrate_by_Euler(){
    std::vector<double> q_new; q_new.resize(Nx);
    std::vector<double> rho_new; rho_new.resize(Nx);
    // navier-stokes equation.
    // advection term.
    for (unsigned int i = 1; i < Nx-1; i++){
        if (q[i] > 0.0){
            q_new[i] = q[i] + dt * ((q[i+1] * q[i+1] / rho[i+1]) - (q[i] * q[i] / rho[i])) / dx;
        } else {
            q_new[i] = q[i] + dt * ((q[i] * q[i] / rho[i]) - (q[i-1] * q[i-1] / rho[i-1])) / dx;
        }
    }
    // pressure gradient term and viscous term.
    for (unsigned int i = 1; i < Nx-1; i++){
        double u_e = q[i-1] / rho[i-1];
        double u_c = q[i] / rho[i];
        double u_w = q[i+1] / rho[i+1];
        q_new[i] = q[i] + dt * (- R * T * (rho[i+1] - rho[i-1]) / (2.0 * dx) + mu * (u_e - 2.0 * u_c + u_w) / (dx * dx));
    }
    // q boundary condition. No Slip.
    q_new[0] = 0.0;
    q_new[Nx-1] = 0.0;
    // continuity equation.
    for(unsigned int i = 1; i < Nx-1; i++){
        rho_new[i] = rho[i] + dt * (- (q[i+1] - q[i-1]) / (2.0 * dx));
    }
    // rho boundary condition.
    rho_new[0] = rho[0];
    rho_new[Nx-1] = rho[Nx-1];
    // update.
    q = q_new;
    rho = rho_new;
    time += dt;
}

void Fluid::integrate_to_end_time_by_Euler(){
    while (time < time_end){
        integrate_by_Euler();
    }
}

int main(){
    Fluid fluid;
    fluid.set_initial_conditions();
    fluid.integrate_to_end_time_by_Euler();
    fluid.print_rho();
}
