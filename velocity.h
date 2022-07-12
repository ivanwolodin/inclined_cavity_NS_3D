#pragma once

#include "utilities.h"
#include "pressure.h"
#include "constants.h"

const int NONLINEAR = 1;
using std::cout;

void velocity_boundary() {
//cout<<"bound"<<endl;
    // left and right walls
    for (int j = 1; j <= Ny - 1; j++) {
        u[0][j] = u[Nx - 1][j];
        u[Nx + 1][j] = u[2][j];
    }

    // upper and bottom walls
    for (int i = 1; i < Nx + 1; i++) {
        u[i][0] = -u[i][1];
        //u[i][Ny] = -u[i][Ny - 1];
        u[i][Ny] = u_auxilliary[i - 1][Ny];
        u[i][Ny - 1] = u_auxilliary[i - 1][Ny - 1];
    }
    u[0][0] = u[Nx + 1][0] = u[0][Ny] = u[Nx + 1][Ny] = 0; // ?
    // ========================================
    u[1][Ny] = u_auxilliary[0][Ny-1];
    // upper and bottom walls
    for (int i = 0; i <= Nx - 1; i++) {
        v[i][0] = 0;
        v[i][Ny - 1] = 0;
    }

    // left and right walls
    for (int j = 1; j < Ny - 1; j++) {
        v[0][j] = v[Nx - 1][j];
        v[Nx][j] = v[1][j];
    }
    v[0][0] = v[0][Ny - 1] = v[Nx][0] = v[Nx][Ny - 1]; // ?
    // ========================================
}


void auxiliary_velocity_boundary() {
//cout<<"aux"<<endl;
    // left and right walls
//cout<<"1"<<endl;
    for (int j = 1; j <= Ny - 1; j++) {

        u_auxilliary[0][j] =
                u[1][j] +
                (dt / Re) * (
                        (u[2][j] - 2 * u[1][j] + u[0][j]) / (dx * dx) +
                        (u[1][j + 1] - 2 * u[1][j] + u[1][j - 1]) / (dy * dy)
                ) -
                NONLINEAR * (
                        (dt / dx) *
                        (
                                0.5 * (u[2][j] + u[1][j]) *
                                0.5 * (u[2][j] + u[1][j]) -

                                0.5 * (u[1][j] + u[0][j]) *
                                0.5 * (u[1][j] + u[0][j])
                        ) +

                        (dt / dy) *
                        (
                                0.5 * (u[1][j + 1] + u[1][j]) *
                                0.5 * (v[2][j] + v[1][j]) -

                                0.5 * (u[1][j] + u[1][j - 1]) *
                                0.5 * (v[2][j - 1] + v[1][j - 1])
                        )
                ) + gravity;
    }
//cout<<"2"<<endl;
    for (int j = 1; j <= Ny - 1; j++) {
        u_auxilliary[Nx - 1][j] = u_auxilliary[0][j];
    }

//cout<<"3"<<endl;
    // bottom and upper walls
    for (int i = 1; i < Nx - 1; i++) {
        u_auxilliary[i][0] = u[i + 1][0] + (dt / (dx)) * (p[i + 1][0] - p[i][0]);
        //u_auxilliary[i][Ny] = u[i][Ny] + (dt / (dx)) * (p[i + 1][Ny] - p[i][Ny]);
        u_auxilliary[i][Ny] = u_auxilliary[i][Ny - 1];
    } 
    u_auxilliary[0][0] = u_auxilliary[Nx - 1][0] = u_auxilliary[0][Ny] = u_auxilliary[Nx - 1][Ny] = 0;

    // ==============================================================================================================
//cout<<"4"<<endl;
    // bottom and upper walls
    for (int i = 1; i <= Nx - 1; i++) {
        v_auxilliary[i][0] = 0 + (dt / (dy)) * (p[i][1] - p[i][0]);
        v_auxilliary[i][Ny - 1] = 0 + (dt / (dy)) * (p[i][Ny] - p[i][Ny - 1]);
    }
//cout<<"5"<<endl;
    // left and right walls
    for (int j = 1; j < Ny - 1; j++) {
        v_auxilliary[0][j] = v_auxilliary[Nx - 1][j];
        v_auxilliary[Nx][j] = v_auxilliary[1][j];

    }

    v_auxilliary[0][0] = v_auxilliary[Nx - 1][0] = v_auxilliary[0][Ny - 1] = v_auxilliary[Nx - 1][Ny - 1] = 0; // ??
    // ==============================================================================================================

}


void calculate_auxiliary_velocity() {
//cout<<"aux_vel"<<endl;
    for (int i = 1; i < Nx - 1; i++) {
        for (int j = 1; j <= Ny - 1; j++) {
            u_auxilliary[i][j] =
                    u[i + 1][j] +
                    (dt / Re) * (
                            (u[i + 2][j] - 2 * u[i + 1][j] + u[i][j]) / (dx * dx) +
                            (u[i + 1][j + 1] - 2 * u[i + 1][j] + u[i + 1][j - 1]) / (dy * dy)
                    ) -
                    NONLINEAR * (
                            (dt / dx) *
                            (
                                    0.5 * (u[i + 1][j] + u[i + 2][j]) *
                                    0.5 * (u[i + 1][j] + u[i + 2][j]) -

                                    0.5 * (u[i][j] + u[i + 1][j]) *
                                    0.5 * (u[i][j] + u[i + 1][j])
                            ) +

                            (dt / dy) *
                            (
                                    0.5 * (u[i + 1][j + 1] + u[i + 1][j]) *
                                    0.5 * (v[i + 1][j] + v[i][j]) -

                                    0.5 * (u[i + 1][j] + u[i + 1][j - 1]) *
                                    0.5 * (v[i + 1][j - 1] + v[i][j - 1])
                            )
                    ) + gravity;
            // ==================================================================================================
        }
    }

    for (int i = 1; i <= Nx - 1; i++) {
        for (int j = 1; j < Ny - 1; j++) {

            v_auxilliary[i][j] =
                    v[i][j] +
                    (dt / Re) * (
                            (v[i + 1][j] - 2 * v[i][j] + v[i - 1][j]) / (dx * dx) +
                            (v[i][j + 1] - 2 * v[i][j] + v[i][j - 1]) / (dy * dy)
                    ) -
                    NONLINEAR * (
                            (dt / dx) *
                            (
                                    0.5 * (u[i][j + 1] + u[i][j]) *
                                    0.5 * (v[i + 1][j] + v[i][j]) -

                                    0.5 * (u[i - 1][j + 1] + u[i - 1][j]) *
                                    0.5 * (v[i][j] + v[i - 1][j])
                            ) +

                            (dt / dy) *
                            (
                                    0.5 * (v[i][j + 1] + v[i][j]) *
                                    0.5 * (v[i][j + 1] + v[i][j]) -

                                    0.5 * (v[i][j] + v[i][j - 1]) *
                                    0.5 * (v[i][j] + v[i][j - 1])
                            )
                    );
            // ==================================================================================================
        }
    }

    auxiliary_velocity_boundary();

}


void calculate_velocity() {
//cout<<"true_vel"<<endl;
    // since u[][][], v[][][] and w[][][] have different size, we need to handle it differently
    for (int i = 1; i <= Nx; i++) {
        for (int j = 1; j < Ny - 1; j++) {
            u[i][j] = u_auxilliary[i - 1][j] - (dt / (dx)) * (p[i][j] - p[i - 1][j]);
        }
    }

    for (int i = 1; i <= Nx - 1; i++) {
        for (int j = 1; j < Ny - 1; j++) {
            v[i][j] = v_auxilliary[i][j] - (dt / (dy)) * (p[i][j + 1] - p[i][j]);
        }
    }

    velocity_boundary();
}
