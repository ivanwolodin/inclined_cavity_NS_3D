#pragma once

#include "constants.h"
#include "iostream"
#include <fstream>
#include <sstream>

using std::cout;
using std::endl;
using std::string;
using std::ofstream;
using std::to_string;
using std::ifstream;
using std:: stringstream;

void zero_values() {
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < Nz; k++) {

                u_auxilliary[i][j][k] = 0.0;
                v_auxilliary[i][j][k] = 0.0;
                w_auxilliary[i][j][k] = 0.0;

                u_previous[i][j][k] = 0.0;
                v_previous[i][j][k] = 0.0;
                w_previous[i][j][k] = 0.0;

                u[i][j][k] = 0.0;
                v[i][j][k] = 0.0;
                w[i][j][k] = 0.0;
                p[i][j][k] = 0.0;

            }
        }
    }
}

void read_data_into_array(vector<vector<vector<double> > > &array, string fileName){
    ifstream inputFile(fileName);        // Input file stream object
    int i, j, k;
    double value;
    while (inputFile >> i >> j >>k >> value){
        array[i][j][k] = value;
    }
}


void initial_distribution() {
    zero_values();

    read_data_into_array(u, "u.txt");
    read_data_into_array(v, "v.txt");
    read_data_into_array(w, "w.txt");
    read_data_into_array(p, "p.txt");

//    float a = 1 / (((float) Nz / 2) * ((float) Nz / 2) - ((float) Nz - 1) * ((float) Nz / 2));
//    float b = -a * ((float) Nz - 1);
//
//    for (int j = 0; j < Ny; j++) {
//        for (int i = 0; i < Nx; i++) {
//            for (int k = 0; k < Nz; k++) {
////                p[i][j][k] = ((1 / (1 - (float) Nz)) * (float) k + 1) * 1;
////                u[i][j][k] = ((1 / (1 - (float) Nz)) * (float) k + 1) * 1;
////                v[i][j][k] = ((1 / (1 - (float) Nz)) * (float) k + 1) * 1;
////                w[i][j][k] = ((1 / (1 - (float) Nz)) * (float) k + 1) * 1;
//                // ((1/ / (1 - (float) Nz)) * (float) k + 1) * 0.1;
////                a * k * k + b * k;
//            }
//        }
//    }
    //u[Nx/2][Ny/2][Nz/2] = 0.99;

    // read from file


}


bool check_convergence(
        vector<vector<vector<double> > > & previous_vector,
        vector<vector<vector<double> > > & current_vector,
const float &precision
) {
    // check convergence
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < Nz; k++) {
                if (abs(current_vector[i][j][k] - previous_vector[i][j][k]) > precision) {
                    return false;
                }
            }
        }
    }
    return true;
}


vector<vector<vector<double> > > copy_previous(
        vector<vector<vector<double> > > & previous_vector,
        vector<vector<vector<double> > > & current_vector
) {

    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < Nz; k++) {
                previous_vector[i][j][k] = current_vector[i][j][k];
            }
        }
    }

    return previous_vector;
}


void output_to_file(const string &fname, const string &step, vector<vector<vector<double> > > &current_vector) {
    time_t t = time(0);   // get time now
    struct tm *now = localtime(&t);

    char buffer[80];
    strftime(buffer, 80, "%d-%m %H.%M", now);

    ofstream myfile;
    myfile.open(string("res/") + buffer + " step=" + step + " " + fname);

    if (myfile.is_open()) {
        myfile << "Re=" << Re << endl;
        myfile << "x-plain x=3\n";
        for (int k = 0; k < Nz; k++) {
            for (int j = 0; j < Ny; j++) {
                myfile << current_vector[2][j][k] << " ";
            }
            myfile << endl;
        }

        myfile << endl;
        myfile << endl;
        myfile << "x-plain x=Nx/2\n";
        for (int k = 0; k < Nz; k++) {
            for (int j = 0; j < Ny; j++) {
                myfile << current_vector[(int) Nx / 2][j][k] << " ";
            }
            myfile << endl;
        }

        myfile << endl;
        myfile << endl;
        myfile << "x-plain x=Nx-3\n";
        for (int k = 0; k < Nz; k++) {
            for (int j = 0; j < Ny; j++) {
                myfile << current_vector[(int) Nx - 3][j][k] << " ";
            }
            myfile << endl;
        }

        myfile << endl;
        myfile << endl;
        myfile << "y-plain y=3\n";
        for (int k = 0; k < Nz; k++) {
            for (int i = 0; i < Nx; i++) {
                myfile << current_vector[i][2][k] << " ";
            }
            myfile << endl;
        }

        myfile << endl;
        myfile << endl;
        myfile << "y-plain y=Ny/2\n";
        for (int k = 0; k < Nz; k++) {
            for (int i = 0; i < Nx; i++) {
                myfile << current_vector[i][(int) Ny / 2][k] << " ";
            }
            myfile << endl;
        }

        myfile << endl;
        myfile << endl;
        myfile << "y-plain y=Ny-2\n";
        for (int k = 0; k < Nz; k++) {
            for (int i = 0; i < Nx; i++) {
                myfile << current_vector[i][(int) Ny - 2][k] << " ";
            }
            myfile << endl;
        }


        myfile << endl;
        myfile << endl;
        myfile << "z-plain z=2\n";
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                myfile << current_vector[i][j][3] << " ";
            }
            myfile << endl;
        }

        myfile << endl;
        myfile << endl;
        myfile << "z-plain z=Nz-2\n";
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                myfile << current_vector[i][j][(int) Nz - 2] << " ";
            }
            myfile << endl;
        }


        myfile.close();
    } else cout << "Unable to open file";
}


void
output_line(const string &fname, const string &step, vector<vector<vector<double> > > &current_vector, int xx, int yy) {
    time_t t = time(0);   // get time now
    struct tm *now = localtime(&t);

    char buffer[80];
    strftime(buffer, 80, "%d-%m %H.%M", now);

    ofstream myfile;
    myfile.open(string("res/") + buffer + " step=" + step + " " + fname);

    if (myfile.is_open()) {

        myfile << "center of Nx-Ny planes\n";
        for (int k = 0; k < Nz; k++) {
            myfile << current_vector[xx][yy][k] << " ";
            myfile << endl;
        }

        myfile.close();
    } else cout << "Unable to open file";
}

void output_plane(const string &fname, const string &step, vector<vector<vector<double> > > &current_vector) {
    time_t t = time(0);   // get time now
    struct tm *now = localtime(&t);

    char buffer[80];
    strftime(buffer, 80, "%d-%m %H.%M", now);

    ofstream myfile;
    myfile.open(string("res/") + buffer + " step=" + step + " " + fname);

    if (myfile.is_open()) {

        myfile << "plane\n";
        for (int i = 0; i < Nx; i++) {
            myfile << "i=";
            myfile << i << endl;
            for (int j = 0; j < Ny; j++) {
                myfile << current_vector[j][i][(int) Nz - 3] << " ";
                myfile << endl;
            }
            myfile << endl;
            myfile << endl;
        }

        myfile.close();
    } else cout << "Unable to open file";
}


void save_previous_velocity() {
    u_previous = copy_previous(u_previous, u);
    v_previous = copy_previous(v_previous, v);
    w_previous = copy_previous(w_previous, w);
}


void save_results(int &step) {
    output_to_file(" u.txt", to_string(step), u);
    output_to_file(" v.txt", to_string(step), v);
    output_to_file(" w.txt", to_string(step), w);
    output_to_file(" p.txt", to_string(step), p);
}


void output_boundary_plane(const string &fname, const string &step, vector<vector<vector<double> > > &current_vector) {
    time_t t = time(0);   // get time now
    struct tm *now = localtime(&t);

    char buffer[80];
    strftime(buffer, 80, "%d-%m %H.%M", now);

    ofstream myfile;
    myfile.open(string("res/") + buffer + " step=" + step + " " + fname);

    if (myfile.is_open()) {

        myfile << "y-plain x=0\n";
        for (int k = 0; k < Nz; k++) {
            for (int j = 0; j < Ny; j++) {
                myfile << current_vector[0][j][k] << " ";
            }
            myfile << endl;
        }

        myfile << endl;
        myfile << endl;
        myfile << "y-plain x=Nx-1\n";
        for (int k = 0; k < Nz; k++) {
            for (int j = 0; j < Ny; j++) {
                myfile << current_vector[Nx - 1][j][k] << " ";
            }
            myfile << endl;
        }

        myfile << endl;
        myfile << endl;
        myfile << "x-plain y=0\n";
        for (int k = 0; k < Nz; k++) {
            for (int i = 0; i < Nx; i++) {
                myfile << current_vector[i][0][k] << " ";
            }
            myfile << endl;
        }

        myfile << endl;
        myfile << endl;
        myfile << "x-plain y=Ny-1\n";
        for (int k = 0; k < Nz; k++) {
            for (int i = 0; i < Nx; i++) {
                myfile << current_vector[i][Ny - 1][k] << " ";
            }
            myfile << endl;
        }

        myfile.close();
    } else cout << "Unable to open file";
}


void output_whole_field(const string &fname, const string &step, vector<vector<vector<double> > > &current_vector) {
    time_t t = time(0);   // get time now
    struct tm *now = localtime(&t);

    char buffer[80];
    strftime(buffer, 80, "%d-%m %H.%M", now);

    ofstream myfile;
    myfile.open(string("res/") + buffer + " step=" + step + " " + fname);

    if (myfile.is_open()) {

        myfile << "whole_field\n";
        myfile << "Re=" << Re << endl;

        for (int i = 0; i < Nx; i++) {
            //myfile << "i="<<i<<endl;
            for (int j = 0; j < Ny; j++) {
                myfile << "i=" << i << "  " << "j=" << j << endl;
                for (int k = 0; k < Nz; k++) {
                    myfile << current_vector[i][j][k] << " " << endl;
                }
                myfile << "_________________________" << endl;
            }
        }

        myfile.close();
    } else cout << "Unable to open file";
}


bool check_velocity_convergence() {
    if (check_convergence(u, u_previous, velocity_calculation_precision) and
        check_convergence(v, v_previous, velocity_calculation_precision) and
        check_convergence(w, w_previous, velocity_calculation_precision)) {
        return true;
    } else {
        return false;
    }
}


void output_all_fields(int &step) {
    output_whole_field(" whole_u.txt", to_string(step), u);
    output_whole_field(" whole_v.txt", to_string(step), v);
    output_whole_field(" whole_w.txt", to_string(step), w);
    output_whole_field(" whole_p.txt", to_string(step), p);
}


bool check_stability() {
    const float stability_condition_one = 0.25 * (mod_u + mod_v + mod_w) * (mod_u + mod_v + mod_w) * dt * Re;
    const float stability_condition_two = dt / (Re * dx * dx);

    if (stability_condition_one <= 1 and stability_condition_two <= 0.25) {
        return true;
    } else {
        return false;
    }
}

void resize_array(vector<vector<vector<double> > > &array, int X, int Y, int Z){
    array.resize(X);
    for (int i = 0; i < X; i++)
    {
        array[i].resize(Y);
        for (int j = 0; j < Y; j++)
        {
            array[i][j].resize(Z);
        }
    }
}
void resize_arrays(){
    resize_array(p, Nx, Ny, Nz);
    resize_array(p_previous, Nx, Ny, Nz);

    resize_array(u, Nx, Ny, Nz);
    resize_array(u_previous, Nx, Ny, Nz);
    resize_array(u_auxilliary, Nx, Ny, Nz);

    resize_array(v, Nx, Ny, Nz);
    resize_array(v_previous, Nx, Ny, Nz);
    resize_array(v_auxilliary, Nx, Ny, Nz);

    resize_array(w, Nx, Ny, Nz);
    resize_array(w_previous, Nx, Ny, Nz);
    resize_array(w_auxilliary, Nx, Ny, Nz);

}
bool read_parameters() {
    ifstream in("params.txt");
    if (!in.is_open()) {
        std::cerr << "Error: Could not open file.\n";
        return false;
    }
    string beforeEqual = "";
    string afterEqual = "";
    while (!in.eof()) {
        getline(in, beforeEqual, '='); //grtting string upto =
        getline(in, afterEqual, '\n'); //getting string after =

        if (beforeEqual == "Re") {
            Re = std::stof(afterEqual);
        }

        if (beforeEqual == "Nx") {
            Nx = std::stoi(afterEqual);
            dx = (float) 4 / Nx;
        }

        if (beforeEqual == "Ny") {
            Ny = std::stoi(afterEqual);
            dy = (float) 4 / Ny;
        }

        if (beforeEqual == "Nz") {
            Nz = std::stoi(afterEqual);
            dz = (float) 1 / Nz;
        }

        if (beforeEqual == "dt") {
            dt = std::stof(afterEqual);
        }

        if (beforeEqual == "gravity") {
            gravity = std::stof(afterEqual);
        }

        if (beforeEqual == "pressure_calculation_precision") {
            pressure_calculation_precision = std::stof(afterEqual);
        }

        if (beforeEqual == "velocity_calculation_precision") {
            velocity_calculation_precision = std::stof(afterEqual);
        }

        if (beforeEqual == "pressure_time") {
            pressure_time = std::stof(afterEqual);
        }

        pressure_part = (dx * dx * dy * dy * dz * dz) / (2 * (dy * dy * dz * dz + dx * dx * dz * dz + dx * dx * dy * dy));
	resize_arrays();
    }
    return true;
}

bool dump_simulation_parameters() {
    time_t t = time(0);   // get time now
    struct tm *now = localtime(&t);

    char buffer[80];
    strftime(buffer, 80, "%d-%m %H.%M", now);

    ofstream myfile;
    myfile.open(string("./") + buffer + " " + "simulation_params.txt");

    if (myfile.is_open()) {

        myfile << "Re=" << Re << endl;
//        cout << "Re=" << Re << endl;

        myfile << "Nx=" << Nx << endl;
//        cout << "Nx=" << Nx << endl;

        myfile << "Ny=" << Ny << endl;
//        cout << "Ny=" << Ny << endl;

        myfile << "Nz=" << Nz << endl;
//        cout << "Nz=" << Nz << endl;

        myfile << "dt=" << dt << endl;
//        cout << "dt=" << dt << endl;

        myfile << "gravity=" << gravity << endl;
//        cout << "gravity=" << gravity << endl;

        myfile << "dx=" << dx << endl;
        myfile << "dy=" << dy << endl;
        myfile << "dz=" << dz << endl;
        myfile << "pressure_part=" << pressure_part << endl;
        myfile << "velocity_calculation_precision=" << velocity_calculation_precision << endl;
        myfile << "pressure_calculation_precision=" << pressure_calculation_precision << endl;

        myfile.close();
        return true;
    } else {
        std::cerr << "Error: Could not open file.\n";
        return false;
    }
}

void output_poisson_solution_dynamics(int mainStep, int poissonStep) {

ofstream myfile;
myfile.open("poisson_dynamics.txt", std::ios_base::app);

if (myfile.is_open()) {
    myfile <<mainStep<< ": " << poissonStep << endl;
    myfile.close();
} else cout << "Unable to open file";

}
