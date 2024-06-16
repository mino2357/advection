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
constexpr double time_end = 0.002;

constexpr bool rho_const = false;

class Fluid {
    private:
        double start_x = -1.0;
        // parameters
        double T = 300.0;
        double R = 8.31;
        double p0 = 101325.0;
        double u0 = 0.0;
        double mu = 0.0; //1.83e-5;
        double time = 0.0;
        std::vector<double> q; // q = rho * u.
        std::vector<double> rho;
        double exact_solution(double end_time, double x);
    public:
        Fluid();
        void set_initial_conditions();
        void print_q();
        void print_u();
        void print_rho();
        void print_p();
        void integrate_by_Euler();
        void integrate_to_end_time_by_Euler();
        void print_exact_solution();
        void print_exact_solution_2_times();
};

Fluid::Fluid(){
    q.resize(Nx);
    rho.resize(Nx);
}

void Fluid::set_initial_conditions(){
    // set rho.
    for (unsigned int i = 0; i < Nx; i++){
        double x = start_x + i * Lx / Nx;
        rho[i] = p0 / (R * T) + 0.001 * exp(- 100.0 * (x - 1.0) * (x - 1.0));
    }
    // set q.
    for (unsigned int i = 0; i < Nx; i++){
        double x = start_x + i * Lx / Nx;
        if (x < 0){
            q[i] = rho[i] * u0;
        } else if (x >= 0 && x < 1){
            q[i] = rho[i] * u0;
        } else {
            q[i] = rho[i] * u0;
        }
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
        if (q[i] >= 0.0){
            q_new[i] = q[i] - dt * ((q[i] * q[i] / rho[i]) - (q[i-1] * q[i-1] / rho[i-1])) / dx;
        } else {
            q_new[i] = q[i] - dt * ((q[i+1] * q[i+1] / rho[i+1]) - (q[i] * q[i] / rho[i])) / dx;
        }
    }
    // pressure gradient term and viscous term.
    for (unsigned int i = 1; i < Nx-1; i++){
        double u_e = q[i-1] / rho[i-1];
        double u_c = q[i] / rho[i];
        double u_w = q[i+1] / rho[i+1];
        q_new[i] += dt * (- R * T * (rho[i+1] - rho[i-1]) / (2.0 * dx) + mu * (u_e - 2.0 * u_c + u_w) / (dx * dx));
    }
    // q boundary condition. Neumann 0.
    q_new[0] = q_new[1];
    q_new[Nx-1] = q_new[Nx-2];
    // continuity equation.
    if(rho_const){
        for(unsigned int i = 1; i < Nx-1; i++){
            rho_new[i] = rho[i];
        }
    } else {
        for(unsigned int i = 1; i < Nx-1; i++){
            if (q[i] >= 0.0){
                rho_new[i] = rho[i] - dt * (q[i]/rho[i] * (rho[i] - rho[i-1]) / dx + rho[i] * (q[i]/rho[i] - q[i-1]/rho[i-1]) / dx);
            } else {
                rho_new[i] = rho[i] - dt * (q[i]/rho[i] * (rho[i+1] - rho[i]) / dx + rho[i] * (q[i+1]/rho[i+1] - q[i]/rho[i]) / dx);
            }
        }
    }
    // rho boundary condition. Neumann 0.
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
    //fluid.print_u();
    fluid.print_rho();
    //fluid.print_exact_solution();
    //fluid.print_exact_solution_2_times();
}
