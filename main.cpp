//
// advection
//

#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

constexpr double Lx = 3.0;
constexpr unsigned int Nx = 300;
constexpr double dx = Lx / (Nx - 1);
constexpr double Courant = 0.5;
constexpr double velocity = 2.0;
constexpr double dt = Courant * dx / velocity;
constexpr double boundary_left = 0.0;
constexpr double boundary_right = 0.0;
constexpr double time_end = 0.2;

void initial_sin(std::vector<double> &s){
    for(unsigned int i = 0; i < s.size(); i++){
        auto x = i * dx - 1.0;
        s[i] = std::sin(2.0 * M_PI * x);
    }
}

void initial_cos(std::vector<double> &s){
    for(unsigned int i = 0; i < s.size(); i++){
        auto x = i * dx - 1.0;
        s[i] = 2.0 * M_PI * std::cos(2.0 * M_PI * x);
    }
}

void initial_zero(std::vector<double> &s){
    for(unsigned int i = 0; i < s.size(); i++){
        s[i] = 0.0;
    }
}

void advection(std::vector<double> &s){
    std::vector<double> s_new(s.size());
    s_new[0] = boundary_left;
    s_new[s.size()-1] = boundary_right;
    for(unsigned int i = 1; i < s.size()-1; i++){
        s_new[i] = s[i] - dt * (s[i] * (s[i+1] - s[i-1]) / (2.0 * dx) - std::abs(s[i]) * (s[i+1] - 2.0 * s[i] + s[i-1]) / (2.0 * dx));
    }
    s = s_new;
}

void advection_2(std::vector<double> &x, std::vector<double> &y){
    std::vector<double> x_t_new(y.size());
    std::vector<double> y_t_new(x.size());
    y_t_new[1]          = - y[1]          * (y[2]           - boundary_left) / (2.0 * dx) + std::abs(y[1])          * (y[2]           - 2.0 * y[1]          + boundary_left) / (2.0 * dx);
    y_t_new[y.size()-2] = - y[y.size()-2] * (boundary_right - y[y.size()-3]) / (2.0 * dx) + std::abs(y[y.size()-2]) * (boundary_right - 2.0 * y[y.size()-2] + y[y.size()-3]) / (2.0 * dx);
    for(unsigned int i = 2; i < x.size()-2; i++){
        x_t_new[i] = - x[i] * y[i];
        y_t_new[i] = - y[i] * y[i]
                     - x[i] * (y[i+1] - y[i-1]) / (2.0 * dx) + std::abs(x[i]) * (y[i+1] - 2.0 * y[i] + y[i-1]) / (2.0 * dx);
    }
    for(unsigned int i = 1; i < x.size()-1; i++){
        x[i] = x[i] + dt * x_t_new[i];
        y[i] = y[i] + dt * y_t_new[i];
    }
}

int main(){
    auto x = std::vector<double>(Nx);
    auto y = std::vector<double>(Nx);
    auto z = std::vector<double>(Nx);
    initial_sin(x);
    initial_cos(y);
    initial_zero(z);
    auto time = 0.0;
    for(unsigned int i = 0; time < time_end; i++){
        time = i * dt;
        advection_2(x, y);
        advection(z);
    }

    for(unsigned int i = 0; i < x.size(); i++){
        auto pos = i * dx - 1.0;
        std::cout << std::setprecision(15) << pos << " " << x[i] << " " << z[i] << std::endl; // exact_solution(pos, time_end) << " " << y[i] << std::endl;
    }
}
