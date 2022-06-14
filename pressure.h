
#pragma once

#include "constants.h"
#include "utilities.h"

const int VELOCITY_PART = 1;

using std::to_string;

void calculate_Nx_pressure() {
    for (int j = 1; j < Ny - 1; j++) {
        for (int k = 1; k < Nz - 1; k++) {
            p[Nx - 1][j][k] = pressure_part * ((p_previous[Nx - 2][j][k] + p_previous[1][j][k]) / (dx * dx) +
                                               (p_previous[Nx - 1][j - 1][k] + p_previous[Nx - 1][j + 1][k]) /
                                               (dy * dy) +
                                               (p_previous[Nx - 1][j][k - 1] + p_previous[Nx - 1][j][k + 1]) /
                                               (dz * dz)) -
                              VELOCITY_PART * (pressure_part / dt) * (
                                      (u_auxilliary[1][j][k] - u_auxilliary[Nx - 2][j][k]) / (2 * dx) +
                                      (v_auxilliary[Nx - 1][j + 1][k] - v_auxilliary[Nx - 1][j - 1][k]) / (2 * dy) +
                                      (w_auxilliary[Nx - 1][j][k + 1] - w_auxilliary[Nx - 1][j][k - 1]) / (2 * dz)
                              );
        }
    }
}


void calculate_Ny_pressure() {
    for (int i = 1; i < Nx - 1; i++) {
        for (int k = 1; k < Nz - 1; k++) {
            p[i][Ny - 1][k] =
                    pressure_part * ((p_previous[i - 1][Ny - 1][k] + p_previous[i + 1][Ny - 1][k]) / (dx * dx) +
                                     (p_previous[i][Ny - 2][k] + p_previous[i][1][k]) / (dy * dy) +
                                     (p_previous[i][Ny - 1][k - 1] + p_previous[i][Ny - 1][k + 1]) / (dz * dz)) -
                    VELOCITY_PART * (pressure_part / dt) * (
                            (u_auxilliary[i + 1][Ny - 1][k] - u_auxilliary[i - 1][Ny - 1][k]) / (2 * dx) +
                            (v_auxilliary[i][1][k] - v_auxilliary[i][Ny - 2][k]) / (2 * dy) +
                            (w_auxilliary[i][Ny - 1][k + 1] - w_auxilliary[i][Ny - 1][k - 1]) / (2 * dz)
                    );
        }
    }
}


void calculate_Ny_Nx_edge_pressure() {
    for (int k = 1; k < Nz - 1; k++) {
        p[Nx - 1][Ny - 1][k] = pressure_part * ((p_previous[Nx - 2][Ny - 1][k] + p_previous[1][Ny - 1][k]) / (dx * dx) +
                                                (p_previous[Nx - 1][Ny - 2][k] + p_previous[Nx - 1][1][k]) / (dy * dy) +
                                                (p_previous[Nx - 1][Ny - 1][k - 1] +
                                                 p_previous[Nx - 1][Ny - 1][k + 1]) / (dz * dz)) -
                               VELOCITY_PART * (pressure_part / dt) * (
                                       (u_auxilliary[1][Ny - 1][k] - u_auxilliary[Nx - 2][Ny - 1][k]) / (2 * dx) +
                                       (v_auxilliary[Nx - 1][1][k] - v_auxilliary[Nx - 1][Ny - 2][k]) / (2 * dy) +
                                       (w_auxilliary[Nx - 1][Ny - 1][k + 1] - w_auxilliary[Nx - 1][Ny - 1][k - 1]) /
                                       (2 * dz)
                               );

    }

    // spread to all edges
    for (int k = 1; k < Nz - 1; k++) {
        p[0][0][k] = p[Nx - 1][Ny - 1][k];
        p[Nx - 1][0][k] = p[Nx - 1][Ny - 1][k];
        p[0][Ny - 1][k] = p[Nx - 1][Ny - 1][k];
    }

}


void apply_pressure_bc() {

    // upper wall
    for (int i = 0; i <= Nx - 1; i++) {
        for (int j = 0; j <= Ny - 1; j++) {
            p[i][j][Nz - 1] = 0;
            p[i][j][0] = 0;
        }
    }

    // periodic Nx-1
    for (int j = 0; j <= Ny - 1; j++) {
        for (int k = 0; k <= Nz - 1; k++) {
            p[0][j][k] = p[Nx - 1][j][k];
        }
    }

    // periodic Ny-1
    for (int i = 0; i <= Nx - 1; i++) {
        for (int k = 0; k <= Nz - 1; k++) {
            p[i][0][k] = p[i][Ny - 1][k];
        }
    }

//    for (int i = 0; i < Nx; i++) {
//        for (int k = 0; k < Nz; k++) {
//            p[i][0][k] = 0.0;
//            p[i][Ny - 1][k] = 0.0;
//        }
//    }
//
//    // x walls
//    for (int j = 0; j < Ny; j++) {
//        for (int k = 0; k < Nz; k++) {
//            p[0][j][k] = 0.0;
//            p[Nx - 1][j][k] = 0.0;
//        }
//    }
}

void calculate_pressure(int mainStep) {
    int step = 1;
    while (true) {
//        cout << "Poisson" << step << endl;
        p_previous = copy_previous(p_previous, p);
        for (int i = 1; i < Nx - 1; i++) {
            for (int j = 1; j < Ny - 1; j++) {
                for (int k = 1; k < Nz - 1; k++) {
//			cout<<i<<" " <<j<<" "<<k<<endl; 
                    p[i][j][k] = p_previous[i][j][k] +
                              pressure_time * (
                                      (p_previous[i - 1][j][k] - 2 * p_previous[i][j][k] + p_previous[i + 1][j][k]) /
                                      (dx * dx) +
                                      (p_previous[i][j - 1][k] - 2 * p_previous[i][j][k] + p_previous[i][j + 1][k]) /
                                      (dy * dy) +
                                      (p_previous[i][j][k - 1] - 2 * p_previous[i][j][k] + p_previous[i][j][k + 1]) /
                                      (dz* dz)
                              ) -
                              VELOCITY_PART * (pressure_time / dt) * (
                                         (u_auxilliary[i + 1][j][k] - u_auxilliary[i - 1][j][k]) / (2 * dx) +
                                         (v_auxilliary[i][j + 1][k] - v_auxilliary[i][j - 1][k]) / (2 * dy) +
                                         (w_auxilliary[i][j][k + 1] - w_auxilliary[i][j][k - 1]) / (2 * dz)
                                 );

                }
            }
        }

        calculate_Nx_pressure();
        calculate_Ny_pressure();
        calculate_Ny_Nx_edge_pressure();
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

