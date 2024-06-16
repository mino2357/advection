//
// 1-dim Navier-Stokes solver.
//

#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

// mesh with equal intervals.
constexpr double Lx = 3.0;
constexpr unsigned int Nx = 129;
constexpr double dx = Lx / (Nx - 1);
// time step, and simulation time.
constexpr double dt = 0.1 * (dx / 360.0);
constexpr double time_end = 0.2;

constexpr bool rho_const = true;

class Fluid {
    private:
        double start_x = -1.0;
        // parameters
        double T = 300.0;
        double R = 8.31;
        double p0 = 101325.0;
        double mu = 0.0; //1.83e-5;
        double time = 0.0;
        std::vector<double> u;
        std::vector<double> rho;
        double exact_solution(double end_time, double x);
    public:
        Fluid();
        void set_initial_conditions();
        void print_u();
        void print_rho();
        void integrate_by_Euler();
        void integrate_to_end_time_by_Euler();
        void print_exact_solution();
        void print_exact_solution_2_times();
};

Fluid::Fluid(){
    u.resize(Nx);
    rho.resize(Nx);
}

void Fluid::set_initial_conditions(){
    // set rho.
    for (unsigned int i = 0; i < Nx; i++){
        //double x = start_x + i * Lx / Nx;
        rho[i] = p0 / (R * T);// + 0.001 * exp(- 100.0 * (x - 1.0) * (x - 1.0));
    }
    // set q.
    for (unsigned int i = 0; i < Nx; i++){
        double x = start_x + i * Lx / Nx;
        if (x < 0){
            u[i] = 1.0;
        } else if (x >= 0 && x < 1){
            u[i] = 2.0;
        } else {
            u[i] = 0.0;
        }
    }
}

void Fluid::print_u(){
    for (unsigned int i = 0; i < Nx; i++){
        double x = start_x + i * dx;
        std::cout << std::setprecision(14) << x << " " << u[i] << std::endl;
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

void Fluid::integrate_by_Euler(){
    std::vector<double> u_new; u_new.resize(Nx);
    std::vector<double> rho_new; rho_new.resize(Nx);
    // navier-stokes equation.
    // pressure gradient term and viscous term.
    for (unsigned int i = 1; i < Nx-1; i++){
        double u_e = u[i-1];
        double u_c = u[i];
        double u_w = u[i+1];
        u_new[i] = u[i] + dt * (- R * T * (rho[i+1] - rho[i-1]) / (2.0 * dx) + mu * (u_e - 2.0 * u_c + u_w) / (dx * dx)) / rho[i];
    }
    // advection term.
    // (fully) linear upwind scheme.
    for (unsigned int i = 2; i < Nx-2; i++){
        if (u[i] >= 0.0){
            double phi_e = u[i] + 0.5 * (u[i] - u[i-1]);
            double phi_w = u[i-1] + 0.5 * (u[i-1] - u[i-2]);
            u_new[i] += - dt * u[i] * (phi_e - phi_w) / dx;
        } else {
            double phi_e = u[i+1] + 0.5 * (u[i+1] - u[i-1]);
            double phi_w = u[i] + 0.5 * (u[i] - u[i-2]);
            u_new[i] += - dt * u[i] * (phi_e - phi_w) / dx;
        }
    }
    // upwind scheme for boundary.
    if (u[1] >= 0.0){
        u_new[1] = u_new[1] + dt * u[0] * (u[1] - u[0]) / dx;
    } else {
        u_new[1] = u_new[1] + dt * u[2] * (u[2] - u[1]) / dx;
    }
    if (u[Nx-2] >= 0.0){
        u_new[Nx-2] = u_new[Nx-2] + dt * u[Nx-3] * (u[Nx-2] - u[Nx-3]) / dx;
    } else {
        u_new[Nx-2] = u_new[Nx-2] + dt * u[Nx-1] * (u[Nx-1] - u[Nx-2]) / dx;
    }
    // u boundary condition. Neumann 0.
    u_new[0] = u_new[1];
    u_new[Nx-1] = u_new[Nx-2];
    // continuity equation.
    if(rho_const){
        for(unsigned int i = 1; i < Nx-1; i++){
            rho_new[i] = rho[i];
        }
    } else {
        for(unsigned int i = 1; i < Nx-1; i++){
            if (u[i] >= 0.0){
                rho_new[i] = rho[i] - dt * (rho[i] * u[i] - rho[i-1] * u[i-1]) / dx;
            } else {
                rho_new[i] = rho[i] - dt * (rho[i+1] * u[i+1] - rho[i] * u[i]) / dx;
            }
        }
    }
    // rho boundary condition. Neumann 0.
    rho_new[0] = rho[0];
    rho_new[Nx-1] = rho[Nx-1];
    // update.
    u = u_new;
    rho = rho_new;
    time += dt;
}

void Fluid::integrate_to_end_time_by_Euler(){
    while (time < time_end){
        integrate_by_Euler();
    }
}

// exact solution. ref: https://math.stackexchange.com/questions/602613/entropy-solution-of-the-burgers-equation
double Fluid::exact_solution(double end_time, double x) {
	if (x <= end_time) {
		return 1.0;
	}
	if (end_time < x && x <= 2.0*end_time) {
		return x / end_time;
	}
	if (2.0*end_time < x && x <= end_time+1.0) {
		return 2.0;
	}
	if (end_time+1.0 <= x) {
		return 0.0;
	}
	return 0.0;
}

void Fluid::print_exact_solution() {
    for (unsigned int i = 0; i < Nx; i++){
        double x = start_x + i * dx;
        std::cout << std::setprecision(14) << x << " " << exact_solution(time, x) << std::endl;
    }
    std::cout << std::endl;
}

void Fluid::print_exact_solution_2_times() {
    for (unsigned int i = 0; i < Nx; i++){
        double x = start_x + i * dx;
        std::cout << std::setprecision(14) << x << " " << exact_solution(2.0*time, x) << std::endl;
    }
    std::cout << std::endl;
}

int main(){
    Fluid fluid;
    fluid.set_initial_conditions();
    fluid.integrate_to_end_time_by_Euler();
    fluid.print_u();
    //fluid.print_rho();
    fluid.print_exact_solution();
    //fluid.print_exact_solution_2_times();
}
