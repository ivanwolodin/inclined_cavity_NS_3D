
#pragma once

#include "constants.h"
#include "utilities.h"

const int VELOCITY_PART = 1;

using std::to_string;

void calculate_Nx_pressure() {
    for (int j = 1; j < Ny - 1; j++) {
        p[Nx - 1][j] = pressure_part * ((p_previous[Nx - 2][j] + p_previous[1][j]) / (dx * dx) +
                                        (p_previous[Nx - 1][j - 1] + p_previous[Nx - 1][j + 1]) /
                                        (dy * dy)
        ) -
                       VELOCITY_PART * (pressure_part / dt) * (
                               (u_auxilliary[1][j] - u_auxilliary[Nx - 2][j]) / (2 * dx) +
                               (v_auxilliary[Nx - 1][j + 1] - v_auxilliary[Nx - 1][j - 1]) / (2 * dy)
                       );
    }
}

void apply_pressure_bc() {

    // periodic Nx-1
    for (int j = 0; j <= Ny - 1; j++) {
        p[0][j] = p[Nx - 1][j];
    }

    // solid
    for (int i = 0; i <= Nx - 1; i++) {
        p[i][Ny - 1] = 0.0;
        p[i][0] = 0.0;
    }
}

void calculate_pressure(int mainStep) {
    int step = 1;
    while (true) {
//        cout << "Poisson" << step << endl;
        p_previous = copy_previous(p_previous, p);
        for (int i = 1; i < Nx - 1; i++) {
            for (int j = 1; j < Ny - 1; j++) {

//			cout<<i<<" " <<j<<" "<<k<<endl; 
                p[i][j] = p_previous[i][j] +
                          pressure_time * (
                                  (p_previous[i - 1][j] - 2 * p_previous[i][j] + p_previous[i + 1][j]) /
                                  (dx * dx) +
                                  (p_previous[i][j - 1] - 2 * p_previous[i][j] + p_previous[i][j + 1]) /
                                  (dy * dy)
                          ) -
                          VELOCITY_PART * (pressure_time / dt) * (
                                  (u_auxilliary[i + 1][j] - u_auxilliary[i - 1][j]) / (2 * dx) +
                                  (v_auxilliary[i][j + 1] - v_auxilliary[i][j - 1]) / (2 * dy)
                          );
            }
        }

        calculate_Nx_pressure();
        apply_pressure_bc();
        if (check_convergence(p, p_previous, pressure_calculation_precision)) {
//            cout << "Poisson calculated" << step << endl;
            //tooutput_to_file(" p_final.txt", to_string(step), p);
            output_poisson_solution_dynamics(mainStep, step);
            return;
        }

//        if (step % 1000 == 0) {
//            output_to_file(" p.txt", to_string(step), p);
//        }

        step++;
    }
}

