#pragma once

#include "utilities.h"
#include "pressure.h"
#include "constants.h"

const int NONLINEAR = 1;


void calculate_Nx_auxiliary_velocity() {
    for (int j = 1; j < Ny - 1; j++) {

        // u
        u_auxilliary[Nx - 1][j] = u[Nx - 1][j] + (dt / Re) * (
                (u[1][j] - 2 * u[Nx - 1][j] + u[Nx - 2][j]) / (dx * dx) +
                (u[Nx - 1][j + 1] - 2 * u[Nx - 1][j] + u[Nx - 1][j - 1]) / (dy * dy)
        ) -
                                  NONLINEAR * (
                                          (dt / dx) *
                                          (u[1][j] * u[1][j] - u[Nx - 1][j] * u[Nx - 1][j]) +

                                          (dt / dy) * (u[Nx - 1][j + 1] * v[Nx - 1][j + 1] -
                                                       u[Nx - 1][j] * v[Nx - 1][j])
                                  ) + gravity;

        // v
        v_auxilliary[Nx - 1][j] = v[Nx - 1][j] + (dt / Re) * (
                (v[1][j] - 2 * v[Nx - 1][j] + v[Nx - 2][j]) / (dx * dx) +
                (v[Nx - 1][j + 1] - 2 * v[Nx - 1][j] + v[Nx - 1][j - 1]) / (dy * dy)
        ) -
                                  NONLINEAR * (
                                          (dt / dx) *
                                          (u[1][j][k] * v[1][j][k] - u[Nx - 1][j][k] * v[Nx - 1][j][k]) +

                                          (dt / dy) * (v[Nx - 1][j + 1][k] * v[Nx - 1][j + 1][k] -
                                                       v[Nx - 1][j][k] * v[Nx - 1][j][k])
                                  );
    }

    for (int j = 1; j < Ny - 1; j++) {
        u_auxilliary[0][j] = u_auxilliary[Nx - 1][j];
        v_auxilliary[0][j] = v_auxilliary[Nx - 1][j];
    }

}

void calculate_Ny_auxiliary_velocity() {
    for (int i = 1; i < Nx - 1; i++) {

        v_auxilliary[i][Ny - 1] = - (dt / dy) * p[i][Ny - 2];
        v_auxilliary[i][0] = - (dt / dy) * p[i][1];

        u_auxilliary[i][0] = (dt / dx) * p[i][1];
        u_auxilliary[i][Ny - 1] = u_auxilliary[i][Ny - 2] + (dt / dx) * p[i][Ny - 3] -  2* (dt / dx) * p[i][Ny - 2];

    }
}

void calculate_auxilliary_velocity() {
    for (int i = 1; i < Nx - 1; i++) {
        for (int j = 1; j < Ny - 1; j++) {
            // u
            u_auxilliary[i][j] = u[i][j] + (dt / Re) * (
                    (u[i + 1][j] - 2 * u[i][j] + u[i - 1][j]) / (dx * dx) +
                    (u[i][j + 1] - 2 * u[i][j] + u[i][j - 1]) / (dy * dy)
            ) -
                                 NONLINEAR * (

                                         (dt / dx) *
                                         (u[i + 1][j] * u[i + 1][j] - u[i][j] * u[i][j]) +

                                         (dt / dy) *
                                         (u[i][j + 1] * v[i][j + 1] - u[i][j] * v[i][j])

                                 ) + gravity;
            // v
            v_auxilliary[i][j] = v[i][j] + (dt / Re) * (
                    (v[i + 1][j] - 2 * v[i][j] + v[i - 1][j]) / (dx * dx) +
                    (v[i][j + 1] - 2 * v[i][j] + v[i][j - 1]) / (dy * dy)
            ) -
                                 NONLINEAR * (

                                         (dt / dx) *
                                         (u[i + 1][j] * v[i + 1][j] - u[i][j] * v[i][j]) +

                                         (dt / dy) *
                                         (v[i][j + 1] * v[i][j + 1] - v[i][j] * v[i][j])
                                 );

        }
    }
    calculate_Nx_auxiliary_velocity();
    calculate_Ny_auxiliary_velocity();
}


void calculate_velocity() {

    for (int i = 1; i < Nx - 1; i++) {
        for (int j = 1; j < Ny - 1; j++) {
            u[i][j] = u_auxilliary[i][j] - (dt / (dx)) * (p[i + 1][j] - p[i][j]);
            v[i][j] = v_auxilliary[i][j] - (dt / (dy)) * (p[i][j + 1] - p[i][j]);
        }
    }

    for (int j = 1; j < Ny - 1; j++) {
        u[Nx - 1][j] = u_auxilliary[Nx - 1][j] - (dt / (dx)) * (p[1][j] - p[Nx - 1][j]);
        v[Nx - 1][j] = v_auxilliary[Nx - 1][j] - (dt / (dy)) * (p[1][j + 1] - p[Nx - 1][j]);
    }

    for (int j = 1; j < Ny - 1; j++) {
        u[0][j] = u[Nx - 1][j];
        v[0][j] = u[Nx - 1][j];
    }

    for (int i = 1; i < Nx - 1; i++) {
        v[i][0] = 0;
        v[i][Ny - 1] = 0;
        u[i][0] = 0;
        u[i][Ny - 1] = u[i][Ny - 2];
    }
}
