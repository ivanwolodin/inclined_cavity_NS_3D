
#pragma once

#include "constants.h"
#include "utilities.h"
#include <omp.h>

const int VELOCITY_PART = 1;

using std::to_string;


void pressure_boundary() {

    // left and right walls
    #pragma omp parallel for
    for (int j = 1; j <= Ny - 1; j++) {
        p[0][j] = p[Nx - 1][j];
        p[Nx][j] = p[1][j];
    }

    // upper and bottom walls
    #pragma omp parallel for
    for (int i = 1; i <= Nx - 1; i++) {
        p[i][0] = -p[i][1];
        p[i][Ny] = -p[i][Ny - 1];
    }

    p[0][0] = p[0][Ny] = p[Nx][0] = p[Nx][Ny] = 0; // ??
}


void calculate_pressure(int mainStep) {
    int step = 1;
    while (true) {
//        if (debug == 1) {
//            cout << "Poisson" << step << endl;
//        }
//        cout << "Poisson: " << step << endl;
//        cout<<pressure_pseudo_time;
        p_previous = copy_previous(p_previous, p, Nx + 1, Ny + 1);
//cout<<"p1"<<endl;
        #pragma omp parallel for
        for (int i = 1; i <= Nx - 1; i++) {
            for (int j = 1; j <= Ny - 1; j++) {
//cout<<i<<" "<<j<<endl;
                p[i][j] = p_previous[i][j] +
                          pressure_pseudo_time * (
                                  (p_previous[i - 1][j] - 2 * p_previous[i][j] + p_previous[i + 1][j]) / (dx * dx) +
                                  (p_previous[i][j - 1] - 2 * p_previous[i][j] + p_previous[i][j + 1]) * dydy2[j]
                          ) -
                          VELOCITY_PART * (pressure_pseudo_time / dt) * (
                                  (u_auxilliary[i][j] - u_auxilliary[i - 1][j]) / (dx) +
                                  (v_auxilliary[i][j] - v_auxilliary[i][j - 1]) * dy1[j]
                          );

            }
        }
//cout<<"p2"<<endl;
        pressure_boundary();
//cout<<"p3"<<endl;
        if (check_convergence(p, p_previous, pressure_calculation_precision, Nx, Ny)) {
            output_poisson_solution_dynamics(mainStep, step);
//            cout << "Poisson calculated" << step << endl<< endl<< endl;
//            output_all_fields(step);
//            output_to_file(" p_final.txt", to_string(step), p);
            return;
        }

//        output_all_fields(step);
//        if (step % 1000 == 0) {
//            output_to_file(" p.txt", to_string(step), p);
//        }

        step++;
    }
}

