#pragma once

#include "utilities.h"
#include "pressure.h"
#include "constants.h"

const int NONLINEAR = 1;


void calculate_Nx_auxiliary_velocity() {
    for (int j = 1; j < Ny - 1; j++) {
        for (int k = 1; k < Nz - 1; k++) {

            // u
            u_auxilliary[Nx - 1][j][k] = u[Nx - 1][j][k] + (dt / Re) * (
                    (u[1][j][k] - 2 * u[Nx - 1][j][k] + u[Nx - 2][j][k]) / (dx * dx) +
                    (u[Nx - 1][j + 1][k] - 2 * u[Nx - 1][j][k] + u[Nx - 1][j - 1][k]) / (dy * dy) +
                    (u[Nx - 1][j][k + 1] - 2 * u[Nx - 1][j][k] + u[Nx - 1][j][k - 1]) / (dz * dz)
            ) -
                                         NONLINEAR * (
                                                 (dt / dx) *
                                                 (u[1][j][k] * u[1][j][k] - u[Nx - 1][j][k] * u[Nx - 1][j][k]) +

                                                 (dt / dy) * (u[Nx - 1][j + 1][k] * v[Nx - 1][j + 1][k] -
                                                              u[Nx - 1][j][k] * v[Nx - 1][j][k]) +

                                                 (dt / dz) * (u[Nx - 1][j][k + 1] * w[Nx - 1][j][k + 1] -
                                                              u[Nx - 1][j][k] * w[Nx - 1][j][k])
                                         );


            // v
            v_auxilliary[Nx - 1][j][k] = v[Nx - 1][j][k] + (dt / Re) * (
                    (v[1][j][k] - 2 * v[Nx - 1][j][k] + v[Nx - 2][j][k]) / (dx * dx) +
                    (v[Nx - 1][j + 1][k] - 2 * v[Nx - 1][j][k] + v[Nx - 1][j - 1][k]) / (dy * dy) +
                    (v[Nx - 1][j][k + 1] - 2 * v[Nx - 1][j][k] + v[Nx - 1][j][k - 1]) / (dz * dz)
            ) -
                                         NONLINEAR * (
                                                 (dt / dx) *
                                                 (u[1][j][k] * v[1][j][k] - u[Nx - 1][j][k] * v[Nx - 1][j][k]) +

                                                 (dt / dy) * (v[Nx - 1][j + 1][k] * v[Nx - 1][j + 1][k] -
                                                              v[Nx - 1][j][k] * v[Nx - 1][j][k]) +

                                                 (dt / dz) * (w[Nx - 1][j][k + 1] * v[Nx - 1][j][k + 1] -
                                                              w[Nx - 1][j][k] * v[Nx - 1][j][k])
                                         ) + gravity;


            // w
            w_auxilliary[Nx - 1][j][k] = w[Nx - 1][j][k] + (dt / Re) * (
                    (w[1][j][k] - 2 * w[Nx - 1][j][k] + w[Nx - 2][j][k]) / (dx * dx) +
                    (w[Nx - 1][j + 1][k] - 2 * w[Nx - 1][j][k] + w[Nx - 1][j - 1][k]) / (dy * dy) +
                    (w[Nx - 1][j][k + 1] - 2 * w[Nx - 1][j][k] + w[Nx - 1][j][k - 1]) / (dz * dz)
            ) -
                                         NONLINEAR * (
                                                 (dt / dx) *
                                                 (u[1][j][k] * w[1][j][k] - u[Nx - 1][j][k] * w[Nx - 1][j][k]) +

                                                 (dt / dy) * (v[Nx - 1][j + 1][k] * w[Nx - 1][j + 1][k] -
                                                              v[Nx - 1][j][k] * w[Nx - 1][j][k]) +

                                                 (dt / dz) * (w[Nx - 1][j][k + 1] * w[Nx - 1][j][k + 1] -
                                                              w[Nx - 1][j][k] * w[Nx - 1][j][k])
                                         );
        }
    }
}


void calculate_Ny_auxiliary_velocity() {
    for (int i = 1; i < Nx - 1; i++) {
        for (int k = 1; k < Nz - 1; k++) {
            // u
            u_auxilliary[i][Ny - 1][k] = u[i][Ny - 1][k] + (dt / Re) * (
                    (u[i + 1][Ny - 1][k] - 2 * u[i][Ny - 1][k] + u[i - 1][Ny - 1][k]) / (dx * dx) +
                    (u[i][1][k] - 2 * u[i][Ny - 1][k] + u[i][Ny - 2][k]) / (dy * dy) +
                    (u[i][Ny - 1][k + 1] - 2 * u[i][Ny - 1][k] + u[i][Ny - 1][k - 1]) / (dz * dz)
            ) -
                                         NONLINEAR * (
                                                 (dt / dx) * (u[i + 1][Ny - 1][k] * u[i + 1][Ny - 1][k] -
                                                              u[i][Ny - 1][k] * u[i][Ny - 1][k]) +

                                                 (dt / dy) *
                                                 (u[i][1][k] * v[i][1][k] - u[i][Ny - 1][k] * v[i][Ny - 1][k]) +

                                                 (dt / dz) * (u[i][Ny - 1][k + 1] * w[i][Ny - 1][k + 1] -
                                                              u[i][Ny - 1][k] * w[i][Ny - 1][k])
                                         );


            // v
            v_auxilliary[i][Ny - 1][k] = v[i][Ny - 1][k] + (dt / Re) * (
                    (v[i + 1][Ny - 1][k] - 2 * v[i][Ny - 1][k] + v[i - 1][Ny - 1][k]) / (dx * dx) +
                    (v[i][1][k] - 2 * v[i][Ny - 1][k] + v[i][Ny - 2][k]) / (dy * dy) +
                    (v[i][Ny - 1][k + 1] - 2 * v[i][Ny - 1][k] + v[i][Ny - 1][k - 1]) / (dz * dz)
            ) -
                                         NONLINEAR * (
                                                 (dt / dx) * (u[i + 1][Ny - 1][k] * v[i + 1][Ny - 1][k] -
                                                              u[i][Ny - 1][k] * v[i][Ny - 1][k]) +

                                                 (dt / dy) *
                                                 (v[i][1][k] * v[i][1][k] - v[i][Ny - 1][k] * v[i][Ny - 1][k]) +

                                                 (dt / dz) * (w[i][Ny - 1][k + 1] * v[i][Ny - 1][k + 1] -
                                                              w[i][Ny - 1][k] * v[i][Ny - 1][k])
                                         ) + gravity;


            // w
            w_auxilliary[i][Ny - 1][k] = w[i][Ny - 1][k] + (dt / Re) * (
                    (w[i + 1][Ny - 1][k] - 2 * w[i][Ny - 1][k] + w[i - 1][Ny - 1][k]) / (dx * dx) +
                    (w[i][1][k] - 2 * w[i][Ny - 1][k] + w[i][Ny - 2][k]) / (dy * dy) +
                    (w[i][Ny - 1][k + 1] - 2 * w[i][Ny - 1][k] + w[i][Ny - 1][k - 1]) / (dz * dz)
            ) -
                                         NONLINEAR * (
                                                 (dt / dx) * (u[i + 1][Ny - 1][k] * w[i + 1][Ny - 1][k] -
                                                              u[i][Ny - 1][k] * w[i][Ny - 1][k]) +

                                                 (dt / dy) *
                                                 (v[i][1][k] * w[i][1][k] - v[i][Ny - 1][k] * w[i][Ny - 1][k]) +

                                                 (dt / dz) * (w[i][Ny - 1][k + 1] * w[i][Ny - 1][k + 1] -
                                                              w[i][Ny - 1][k] * w[i][Ny - 1][k])
                                         );
        }
    }
}


void calculate_Ny_Nx_edge_auxiliary_velocity() {
    for (int k = 1; k < Nz - 1; k++) {

        // u
        u_auxilliary[Nx - 1][Ny - 1][k] = u[Nx - 1][Ny - 1][k] + (dt / Re) * (
                (u[1][Ny - 1][k] - 2 * u[Nx - 1][Ny - 1][k] + u[Nx - 2][Ny - 1][k]) / (dx * dx) +
                (u[Nx - 1][1][k] - 2 * u[Nx - 1][Ny - 1][k] + u[Nx - 1][Ny - 2][k]) / (dy * dy) +
                (u[Nx - 1][Ny - 1][k + 1] - 2 * u[Nx - 1][Ny - 1][k] + u[Nx - 1][Ny - 1][k - 1]) / (dz * dz)
        ) -
                                          NONLINEAR * (
                                                  (dt / dx) * (u[1][Ny - 1][k] * u[1][Ny - 1][k] -
                                                               u[Nx - 1][Ny - 1][k] * u[Nx - 1][Ny - 1][k]) +

                                                  (dt / dy) *
                                                  (u[Nx - 1][1][k] * v[Nx - 1][1][k] -
                                                   u[Nx - 1][Ny - 1][k] * v[Nx - 1][Ny - 1][k]) +

                                                  (dt / dz) *
                                                  (u[Nx - 1][Ny - 1][k + 1] * w[Nx - 1][Ny - 1][k + 1] -
                                                   u[Nx - 1][Ny - 1][k] * w[Nx - 1][Ny - 1][k])
                                          );


        // v
        v_auxilliary[Nx - 1][Ny - 1][k] = v[Nx - 1][Ny - 1][k] + (dt / Re) * (
                (v[1][Ny - 1][k] - 2 * v[Nx - 1][Ny - 1][k] + v[Nx - 2][Ny - 1][k]) / (dx * dx) +
                (v[Nx - 1][1][k] - 2 * v[Nx - 1][Ny - 1][k] + v[Nx - 1][Ny - 2][k]) / (dy * dy) +
                (v[Nx - 1][Ny - 1][k + 1] - 2 * v[Nx - 1][Ny - 1][k] + v[Nx - 1][Ny - 1][k - 1]) / (dz * dz)
        ) -
                                          NONLINEAR * (
                                                  (dt / dx) * (u[1][Ny - 1][k] * v[1][Ny - 1][k] -
                                                               u[Nx - 1][Ny - 1][k] * v[Nx - 1][Ny - 1][k]) +

                                                  (dt / dy) *
                                                  (v[Nx - 1][1][k] * v[Nx - 1][1][k] -
                                                   v[Nx - 1][Ny - 1][k] * v[Nx - 1][Ny - 1][k]) +

                                                  (dt / dz) *
                                                  (w[Nx - 1][Ny - 1][k + 1] * v[Nx - 1][Ny - 1][k + 1] -
                                                   w[Nx - 1][Ny - 1][k] * v[Nx - 1][Ny - 1][k])
                                          ) + gravity;


        // w
        w_auxilliary[Nx - 1][Ny - 1][k] = w[Nx - 1][Ny - 1][k] + (dt / Re) * (
                (w[1][Ny - 1][k] - 2 * w[Nx - 1][Ny - 1][k] + w[Nx - 2][Ny - 1][k]) / (dx * dx) +
                (w[Nx - 1][1][k] - 2 * w[Nx - 1][Ny - 1][k] + w[Nx - 1][Ny - 2][k]) / (dy * dy) +
                (w[Nx - 1][Ny - 1][k + 1] - 2 * w[Nx - 1][Ny - 1][k] + w[Nx - 1][Ny - 1][k - 1]) / (dz * dz)
        ) -
                                          NONLINEAR * (
                                                  (dt / dx) * (u[1][Ny - 1][k] * w[1][Ny - 1][k] -
                                                               u[Nx - 1][Ny - 1][k] * w[Nx - 1][Ny - 1][k]) +

                                                  (dt / dy) *
                                                  (v[Nx - 1][1][k] * w[Nx - 1][1][k] -
                                                   v[Nx - 1][Ny - 1][k] * w[Nx - 1][Ny - 1][k]) +

                                                  (dt / dz) *
                                                  (w[Nx - 1][Ny - 1][k + 1] * w[Nx - 1][Ny - 1][k + 1] -
                                                   w[Nx - 1][Ny - 1][k] * w[Nx - 1][Ny - 1][k])
                                          );
    }

    // spread to all edges
    for (int k = 1; k < Nz - 1; k++) {
        u_auxilliary[0][0][k] = u_auxilliary[Nx - 1][Ny - 1][k];
        u_auxilliary[Nx - 1][0][k] = u_auxilliary[Nx - 1][Ny - 1][k];
        u_auxilliary[0][Ny - 1][k] = u_auxilliary[Nx - 1][Ny - 1][k];

        v_auxilliary[0][0][k] = v_auxilliary[Nx - 1][Ny - 1][k];
        v_auxilliary[Nx - 1][0][k] = v_auxilliary[Nx - 1][Ny - 1][k];
        v_auxilliary[0][Ny - 1][k] = v_auxilliary[Nx - 1][Ny - 1][k];

        w_auxilliary[0][0][k] = w_auxilliary[Nx - 1][Ny - 1][k];
        w_auxilliary[Nx - 1][0][k] = w_auxilliary[Nx - 1][Ny - 1][k];
        w_auxilliary[0][Ny - 1][k] = w_auxilliary[Nx - 1][Ny - 1][k];
    }

}


void apply_auxilliary_velocity_bc() {

    // upper wall
    for (int i = 0; i <= Nx - 1; i++) {
        for (int j = 0; j <= Ny - 1; j++) {
            u_auxilliary[i][j][0] = 0;
            v_auxilliary[i][j][0] = 0;
            w_auxilliary[i][j][0] = (p[i][j][1] * dt) / dz;

            // u_auxilliary[i][j][Nz - 1] = 0;
            // v_auxilliary[i][j][Nz - 1] = 0;
            w_auxilliary[i][j][Nz - 1] = -(p[i][j][Nz - 2] * dt) / dz;
        }
    }

    for (int i = 0; i < Nx - 1; i++) {
        for (int j = 0; j <= Ny - 1; j++) {
            u_auxilliary[i][j][Nz - 1] =
                    u_auxilliary[i][j][Nz - 2] - (dt / dx) * (p[i + 1][j][Nz - 2] - p[i][j][Nz - 2]);
        }
    }

    for (int j = 0; j <= Ny - 1; j++) {
        u_auxilliary[Nx - 1][j][Nz - 1] =
                u_auxilliary[Nx - 1][j][Nz - 2] - (dt / dx) * (p[1][j][Nz - 2] - p[Nx - 1][j][Nz - 2]);
    }

    for (int i = 0; i <= Nx - 1; i++) {
        for (int j = 0; j < Ny - 1; j++) {
            v_auxilliary[i][j][Nz - 1] =
                    v_auxilliary[i][j][Nz - 2] - (dt / dy) * (p[i][j + 1][Nz - 2] - p[i][j][Nz - 2]);
        }
    }

    for (int i = 0; i <= Nx - 1; i++) {
        v_auxilliary[i][Ny - 1][Nz - 1] =
                v_auxilliary[i][Ny - 1][Nz - 2] - (dt / dy) * (p[i][1][Nz - 2] - p[i][Ny - 1][Nz - 2]);
    }


    // periodic Nx-1
    for (int j = 0; j <= Ny - 1; j++) {
        for (int k = 1; k < Nz - 1; k++) {
            u_auxilliary[0][j][k] = u_auxilliary[Nx - 1][j][k];
            v_auxilliary[0][j][k] = v_auxilliary[Nx - 1][j][k];
            w_auxilliary[0][j][k] = w_auxilliary[Nx - 1][j][k];
        }
    }

    // periodic Ny-1
    for (int i = 0; i <= Nx - 1; i++) {
        for (int k = 1; k < Nz - 1; k++) {
            u_auxilliary[i][0][k] = u_auxilliary[i][Ny - 1][k];
            v_auxilliary[i][0][k] = v_auxilliary[i][Ny - 1][k];
            w_auxilliary[i][0][k] = w_auxilliary[i][Ny - 1][k];
        }
    }
}


void calculate_auxilliary_velocity() {
    for (int i = 1; i < Nx - 1; i++) {
        for (int j = 1; j < Ny - 1; j++) {
            for (int k = 1; k < Nz - 1; k++) {
                // u
                u_auxilliary[i][j][k] = u[i][j][k] + (dt / Re) * (
                        (u[i + 1][j][k] - 2 * u[i][j][k] + u[i - 1][j][k]) / (dx * dx) +
                        (u[i][j + 1][k] - 2 * u[i][j][k] + u[i][j - 1][k]) / (dy * dy) +
                        (u[i][j][k + 1] - 2 * u[i][j][k] + u[i][j][k - 1]) / (dz * dz)
                ) -
                                        NONLINEAR * (

                                                (dt / dx) *
                                                (u[i + 1][j][k] * u[i + 1][j][k] - u[i][j][k] * u[i][j][k]) +

                                                (dt / dy) *
                                                (u[i][j + 1][k] * v[i][j + 1][k] - u[i][j][k] * v[i][j][k]) +

                                                (dt / dz) *
                                                (u[i][j][k + 1] * w[i][j][k + 1] - u[i][j][k] * w[i][j][k])

                                        );


                // v
                v_auxilliary[i][j][k] = v[i][j][k] + (dt / Re) * (
                        (v[i + 1][j][k] - 2 * v[i][j][k] + v[i - 1][j][k]) / (dx * dx) +
                        (v[i][j + 1][k] - 2 * v[i][j][k] + v[i][j - 1][k]) / (dy * dy) +
                        (v[i][j][k + 1] - 2 * v[i][j][k] + v[i][j][k - 1]) / (dz * dz)
                ) -
                                        NONLINEAR * (

                                                (dt / dx) *
                                                (u[i + 1][j][k] * v[i + 1][j][k] - u[i][j][k] * v[i][j][k]) +

                                                (dt / dy) *
                                                (v[i][j + 1][k] * v[i][j + 1][k] - v[i][j][k] * v[i][j][k]) +

                                                (dt / dz) *
                                                (w[i][j][k + 1] * v[i][j][k + 1] - w[i][j][k] * v[i][j][k])
                                        ) + gravity;


                // w
                w_auxilliary[i][j][k] = w[i][j][k] + (dt / Re) * (
                        (w[i + 1][j][k] - 2 * w[i][j][k] + w[i - 1][j][k]) / (dx * dx) +
                        (w[i][j + 1][k] - 2 * w[i][j][k] + w[i][j - 1][k]) / (dy * dy) +
                        (w[i][j][k + 1] - 2 * w[i][j][k] + w[i][j][k - 1]) / (dz * dz)
                ) -
                                        NONLINEAR * (

                                                (dt / dx) *
                                                (u[i + 1][j][k] * w[i + 1][j][k] - u[i][j][k] * w[i][j][k]) +

                                                (dt / dy) *
                                                (v[i][j + 1][k] * w[i][j + 1][k] - v[i][j][k] * w[i][j][k]) +

                                                (dt / dz) *
                                                (w[i][j][k + 1] * w[i][j][k + 1] - w[i][j][k] * w[i][j][k])
                                        );
            }
        }
    }
    calculate_Nx_auxiliary_velocity();
    calculate_Ny_auxiliary_velocity();
    calculate_Ny_Nx_edge_auxiliary_velocity();
    apply_auxilliary_velocity_bc();
}


void calculate_velocity() {

    for (int i = 1; i < Nx - 1; i++) {
        for (int j = 1; j < Ny - 1; j++) {
            for (int k = 1; k < Nz - 1; k++) {
                u[i][j][k] = u_auxilliary[i][j][k] - (dt / (dx)) * (p[i + 1][j][k] - p[i][j][k]);
                v[i][j][k] = v_auxilliary[i][j][k] - (dt / (dy)) * (p[i][j + 1][k] - p[i][j][k]);
                w[i][j][k] = w_auxilliary[i][j][k] - (dt / (dz)) * (p[i][j][k + 1] - p[i][j][k]);
            }
        }
    }

    for (int i = 0; i <= Nx - 1; i++) {
        for (int j = 0; j <= Ny - 1; j++) {

            u[i][j][Nz - 1] = u_auxilliary[i][j][Nz - 1];
            u[i][j][0] = 0;

            v[i][j][Nz - 1] = v_auxilliary[i][j][Nz - 1];
            v[i][j][0] = 0;

            w[i][j][Nz - 1] = 0;  // w_auxilliary[i][j][Nz - 1] - (dt / (dz)) * (p[i][j][1] - p[i][j][Nz - 1]);
            w[i][j][0] = 0;       // w_auxilliary[i][j][0] - (dt / (dz)) * (p[i][j][1] - p[i][j][0]);

        }
    }

    // Nx - 1
    for (int j = 1; j < Ny - 1; j++) {
        for (int k = 1; k < Nz - 1; k++) {
            u[Nx - 1][j][k] = u_auxilliary[Nx - 1][j][k] - (dt / (dx)) * (p[1][j][k] - p[Nx - 1][j][k]);
            v[Nx - 1][j][k] =
                    v_auxilliary[Nx - 1][j][k] - (dt / (dy)) * (p[Nx - 1][j + 1][k] - p[Nx - 1][j][k]);
            w[Nx - 1][j][k] =
                    w_auxilliary[Nx - 1][j][k] - (dt / (dz)) * (p[Nx - 1][j][k + 1] - p[Nx - 1][j][k]);
        }
    }

    // Ny - 1
    for (int i = 1; i < Nx - 1; i++) {
        for (int k = 1; k < Nz - 1; k++) {
            u[i][Ny - 1][k] =
                    u_auxilliary[i][Ny - 1][k] - (dt / (dx)) * (p[i + 1][Ny - 1][k] - p[i][Ny - 1][k]);
            v[i][Ny - 1][k] = v_auxilliary[i][Ny - 1][k] - (dt / (dy)) * (p[i][1][k] - p[i][Ny - 1][k]);
            w[i][Ny - 1][k] =
                    w_auxilliary[i][Ny - 1][k] - (dt / (dz)) * (p[i][Ny - 1][k + 1] - p[i][Ny - 1][k]);
        }
    }

    // edge
    for (int k = 1; k < Nz - 1; k++) {
        u[Nx - 1][Ny - 1][k] =
                u_auxilliary[Nx - 1][Ny - 1][k] - (dt / (dx)) * (p[1][Ny - 1][k] - p[Nx - 1][Ny - 1][k]);
        v[Nx - 1][Ny - 1][k] =
                v_auxilliary[Nx - 1][Ny - 1][k] - (dt / (dy)) * (p[Nx - 1][1][k] - p[Nx - 1][Ny - 1][k]);
        w[Nx - 1][Ny - 1][k] =
                w_auxilliary[Nx - 1][Ny - 1][k] -
                (dt / (dz)) * (p[Nx - 1][Ny - 1][k + 1] - p[Nx - 1][Ny - 1][k]);
    }

    // spread to all edges
    for (int k = 1; k < Nz - 1; k++) {

        u[0][0][k] = u[Nx - 1][Ny - 1][k];
        u[Nx - 1][0][k] = u[Nx - 1][Ny - 1][k];
        u[0][Ny - 1][k] = u[Nx - 1][Ny - 1][k];

        v[0][0][k] = v[Nx - 1][Ny - 1][k];
        v[Nx - 1][0][k] = v[Nx - 1][Ny - 1][k];
        v[0][Ny - 1][k] = v[Nx - 1][Ny - 1][k];

        w[0][0][k] = w[Nx - 1][Ny - 1][k];
        w[Nx - 1][0][k] = w[Nx - 1][Ny - 1][k];
        w[0][Ny - 1][k] = w[Nx - 1][Ny - 1][k];

    }

    // periodic 0 <--> Nx-1
    for (int j = 0; j <= Ny - 1; j++) {
        for (int k = 0; k <= Nz - 1; k++) {
            u[0][j][k] = u[Nx - 1][j][k];
            v[0][j][k] = v[Nx - 1][j][k];
            w[0][j][k] = w[Nx - 1][j][k];
        }
    }

    // periodic 0 <--> Ny-1
    for (int i = 0; i <= Nx - 1; i++) {
        for (int k = 0; k <= Nz - 1; k++) {
            u[i][0][k] = u[i][Ny - 1][k];
            v[i][0][k] = v[i][Ny - 1][k];
            w[i][0][k] = w[i][Ny - 1][k];
        }
    }

}
