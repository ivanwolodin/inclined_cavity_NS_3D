#pragma once

#include "constants.h"
#include "iostream"
#include <fstream>

using std::cout;
using std::endl;
using std::string;
using std::ofstream;
using std::to_string;
using std::ifstream;
using std:: stringstream;

void read_data_into_array(vector<vector<double> > &array, string fileName){
    ifstream inputFile(fileName);        // Input file stream object
    int i, k;
    double value;

    while (inputFile >> i >>k >> value){
        array[i][k] = value;
        if(i == Nx or k == Ny){
            continue;
        }
    }
}

void zero_values() {
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
                u_auxilliary[i][j] = 0.0;
                v_auxilliary[i][j] = 0.0;

                u_previous[i][j] = 0.0;
                v_previous[i][j] = 0.0;

                u[i][j] = 0.0;
                v[i][j] = 0.0;

                p[i][j] = 0.0;
        }
    }
}


void initial_distribution() {
    zero_values();
    read_data_into_array(u, "u.txt");
    read_data_into_array(v, "v.txt");
    read_data_into_array(p, "p.txt");
}


bool check_convergence(
        vector<vector<double> > & previous_vector,
        vector<vector<double> > & current_vector,
const float &precision
) {
    // check convergence
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            if (abs(current_vector[i][j] - previous_vector[i][j]) > precision) {
                return false;
            }
        }
    }
    return true;
}


vector<vector<double> > copy_previous(
        vector<vector<double> > & previous_vector,
        vector<vector<double> > & current_vector
) {
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            previous_vector[i][j] = current_vector[i][j];
        }
    }
    return previous_vector;
}

void save_previous_velocity() {
    u_previous = copy_previous(u_previous, u);
    v_previous = copy_previous(v_previous, v);
}


void output_whole_field(const string &fname, const string &step, vector<vector<double> > &current_vector) {
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
            myfile << "i=" << i << endl;
            for (int j = 0; j < Ny; j++) {
                myfile << current_vector[i][j] << " " << endl;
            }
            myfile << "_________________________" << endl;
        }
        myfile.close();

    } else cout << "Unable to open file";
}


bool check_velocity_convergence() {
    if (check_convergence(u, u_previous, velocity_calculation_precision) and
        check_convergence(v, v_previous, velocity_calculation_precision)
        ) {
        return true;
    } else {
        return false;
    }
}


void output_all_fields(int &step) {
    output_whole_field(" whole_u.txt", to_string(step), u);
    output_whole_field(" whole_v.txt", to_string(step), v);
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

void resize_array(vector<vector<double> > &array, int X, int Y){
    array.resize(X);
    for (int i = 0; i < X; i++)
    {
        array[i].resize(Y);
    }
}
void resize_arrays(){
    resize_array(p, Nx, Ny);
    resize_array(p_previous, Nx, Ny);

    resize_array(u, Nx, Ny);
    resize_array(u_previous, Nx, Ny);
    resize_array(u_auxilliary, Nx, Ny);

    resize_array(v, Nx, Ny);
    resize_array(v_previous, Nx, Ny);
    resize_array(v_auxilliary, Nx, Ny);

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

        pressure_part = (dx * dx * dy * dy) / (2 * (dy * dy + dx * dx + dx * dx * dy * dy));
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

//        myfile << "Nz=" << Nz << endl;
//        cout << "Nz=" << Nz << endl;

        myfile << "dt=" << dt << endl;
//        cout << "dt=" << dt << endl;

        myfile << "gravity=" << gravity << endl;
//        cout << "gravity=" << gravity << endl;

        myfile << "dx=" << dx << endl;
        myfile << "dy=" << dy << endl;
//        myfile << "dz=" << dz << endl;
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
