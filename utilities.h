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
/*
void read_data_into_array(vector<vector<double> >  &array, string fileName, int X, int Y){
    ifstream inputFile(fileName);        // Input file stream object
    int i, j;
    double value;

    vector<vector<double>  > temporary;
    temporary.resize(X);

    for (int i = 0; i < X; i++) {
        temporary[i].resize(Y);
    }

    while (inputFile >> i >>j >> value){
        temporary[i][j] = value;
    }

    for (int i = 0; i < X; i++) {
        for (int j = 0; j < Y; j++) {
            if (j >= 200 and j <= 500) {
                array[i][j] = temporary[i][199];
            }
            if (j > 500) {
                array[i][j] = temporary[i][j - 300];
            } else {
                array[i][j] = temporary[i][j];
            }
        }
    }
}
 */

void read_data_into_array(vector<vector<double> >  &array, string fileName){
    ifstream inputFile(fileName);        // Input file stream object
    int i, k;
    double value;

    while (inputFile >> i >>k >> value){
        array[i][k] = value;
    }
}

void zero_values() {
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
//            for (int k = 0; k < Nz; k++) {

            u_auxilliary[i][j] = 0.0;
            v_auxilliary[i][j] = 0.0;

            u_previous[i][j] = 0.0;
            v_previous[i][j] = 0.0;

            u[i][j] = 0.0;
            v[i][j] = 0.0;
            p[i][j] = 0.0;

//            }
        }
    }
}


void initial_distribution() {
    zero_values();
    read_data_into_array(u, "u.txt");
    read_data_into_array(v, "v.txt");
    read_data_into_array(p, "p.txt");

//    float a = 1 / (((float) Nz / 2) * ((float) Nz / 2) - ((float) Nz - 1) * ((float) Nz / 2));
//    float b = -a * ((float) Nz - 1);
//
//    for (int j = 0; j < Ny; j++) {
//        for (int i = 0; i < Nx; i++) {
////            for (int k = 0; k < Nz; k++) {
//            p[i][j] = 0;//((1 / (1 - (float) Ny)) * (float) j + 1) * 1;
//            u[i][j] = 0;//((1 / (1 - (float) Ny)) * (float) j + 1) * 1;
//            v[i][j] = 0;//((1 / (1 - (float) Ny)) * (float) j + 1) * 1;
////                w[i][j] = ((1 / (1 - (float) Ny)) * (float) j + 1) * 1;
//            // ((1/ / (1 - (float) Nz)) * (float) k + 1) * 0.1;
////                a * k * k + b * k;
////            }
//        }
//    }
////    p[Nx/2][Ny/2][Nz/2] = 0.01;
////    u[Nx/2][Ny/2][Nz/2] = 0.99;
}


bool check_convergence(
        vector <vector<double> > &previous_vector,
        vector <vector<double> > &current_vector,
        const float &precision,
        int X, int Y
) {
    // check convergence
    for (int i = 1; i < X; i++) {
        for (int j = 1; j < Y; j++) {
            if (abs(current_vector[i][j] - previous_vector[i][j]) > precision) {
//                    cout<<"False in"<< i<<" "<< j <<" "<< k<<" "<<abs(current_vector[i][j][k] - previous_vector[i][j][k]) <<endl;
                return false;
            }
        }
    }
    return true;
}


vector <vector<double> > copy_previous(
        vector <vector<double> > &previous_vector,
        vector <vector<double> > &current_vector,
        int X, int Y
) {

    for (int i = 0; i < X; i++) {
        for (int j = 0; j < Y; j++) {
            previous_vector[i][j] = current_vector[i][j];
        }
    }

    return previous_vector;
}


void save_previous_velocity() {
    u_previous = copy_previous(u_previous, u, Nx + 2, Ny + 1);
    v_previous = copy_previous(v_previous, v, Nx + 1, Ny);
//    w_previous = copy_previous(w_previous, w, Nx + 1, Ny + 1, Nz);
}


void output_whole_field(
        const string &fname,
        const string &step,
        vector <vector<double> > &current_vector,
        int X,
        int Y
) {
    time_t t = time(0);   // get time now
    struct tm *now = localtime(&t);

    char buffer[80];
    strftime(buffer, 80, "%d-%m %H.%M", now);

    ofstream myfile;
    myfile.open(string("res/") + buffer + " step=" + step + " " + fname);

    if (myfile.is_open()) {

        myfile << "whole_field\n";
        myfile << "Re=" << Re << endl;

        for (int i = 0; i < X; i++) {
            myfile << "i=" << i << "  "<< endl;
            for (int j = 0; j < Y; j++) {
                myfile << current_vector[i][j] << " " << endl;
            }
             myfile << "_________________________" << endl;
        }

        myfile.close();
    } else cout << "Unable to open file";
}


bool check_velocity_convergence() {
    if (check_convergence(u, u_previous, velocity_calculation_precision, Nx + 2, Ny + 1) and
        check_convergence(v, v_previous, velocity_calculation_precision, Nx + 1, Ny)
//        and
//        check_convergence(w, w_previous, velocity_calculation_precision, Nx + 1, Ny + 1)
            ) {
        return true;
    } else {
        return false;
    }
}


void output_all_fields(int &step) {
    output_whole_field(
            " whole_u.txt",
            to_string(step),
            u,
            Nx,
            Ny + 1
    );
    output_whole_field(
            " whole_v.txt",
            to_string(step),
            v,
            Nx + 1,
            Ny
    );

    output_whole_field(
            " whole_p.txt",
            to_string(step),
            p,
            Nx + 1,
            Ny + 1
    );
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


void resize_array(vector <vector<double> > &array, int X, int Y) {
    array.resize(X);
    for (int i = 0; i < X; i++) {
        array[i].resize(Y);
//        for (int j = 0; j < Y; j++)
//        {
//            array[i][j].resize(Z);
//        }
    }
}


void resize_arrays() {
    resize_array(p, Nx + 1, Ny + 1);
    resize_array(p_previous, Nx + 1, Ny + 1);

    resize_array(u, Nx + 2, Ny + 1);
    resize_array(u_previous, Nx + 2, Ny + 1);
    resize_array(u_auxilliary, Nx, Ny + 1);

    resize_array(v, Nx + 1, Ny);
    resize_array(v_previous, Nx + 1, Ny);
    resize_array(v_auxilliary, Nx + 1, Ny);

//    resize_array(w, Nx + 1, Ny + 1, Nz);
//    resize_array(w_previous, Nx + 1, Ny + 1, Nz);



//    resize_array(w_auxilliary, Nx, Ny, Nz);
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
            dx = (float) 5 / Nx;
        }

        if (beforeEqual == "Ny") {
            Ny = std::stoi(afterEqual);
//            dy = (float) 1 / Ny;
            dy.resize(Ny + 1);
            for (int j = 0; j < Ny + 1; j++) {
                if (j <= 70) {
                    dy[j] = 0.002;
                }
                if (j > 70 && j <= 250) {
                    dy[j] = 0.002;
                }
                if (j > 250) {
                    dy[j] = 0.002;
                }
            }
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

        if (beforeEqual == "debug") {
            debug = std::stoi(afterEqual);
        }
        if (beforeEqual == "velocity_calculation_precision") {
            velocity_calculation_precision = std::stof(afterEqual);
        }
        if (beforeEqual == "pressure_calculation_precision") {
            pressure_calculation_precision = std::stof(afterEqual);
        }
        if (beforeEqual == "pressure_pseudo_time") {
            pressure_pseudo_time = std::stof(afterEqual);
        }

    }
//    pressure_part =
//            (dx * dx * dy * dy * dz * dz) / (2 * (dy * dy * dz * dz + dx * dx * dz * dz + dx * dx * dy * dy));
    resize_arrays();

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

        myfile << "Nx=" << Nx << endl;

        myfile << "Ny=" << Ny << endl;

        myfile << "Nz=" << Nz << endl;

        myfile << "dt=" << dt << endl;

        myfile << "gravity=" << gravity << endl;

        myfile << "dx=" << dx << endl;
//        myfile << "dy=" << dy << endl;

        for (int j = 0; j < Ny + 1; j++) {
            myfile << "dy["<<j<<"]= " << dy[j] << endl;
        }

        myfile << "dz=" << dz << endl;
        myfile << "poisson_precision=" << pressure_calculation_precision << endl;
        myfile << "velocity_precision=" << velocity_calculation_precision << endl;


        myfile.close();

        if (debug == 1) {
            cout << "Re=" << Re << endl;
            cout << "Nx=" << Nx << endl;
            cout << "Ny=" << Ny << endl;
            cout << "Nz=" << Nz << endl;
            cout << "dt=" << dt << endl;
            cout << "dx=" << dx << endl;
//            cout << "dy=" << dy << endl;
            cout << "dz=" << dz << endl;
            cout << "gravity=" << gravity << endl;
//            cout << "pressure_part=" << pressure_part << endl;
        }

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
        myfile << mainStep << ": " << poissonStep << endl;
        myfile.close();
    } else cout << "Unable to open file";

}

void output_averaged_veloicty(int mainStep) {

    time_t t = time(0);   // get time now
    struct tm *now = localtime(&t);

    char buffer[80];
    strftime(buffer, 80, "%d-%m %H.%M", now);

    ofstream myfile;
    myfile.open(string("res/") + buffer + " step=" + to_string(mainStep) + " " + "averaged_velocity");

    if (myfile.is_open()) {

       
        myfile << "mainStep=" << mainStep << endl;
        
//        for (int i = 0; i < X; i++) {
//            myfile << "i=" << i << "  "<< endl;
//            for (int j = 0; j < Y; j++) {
//                myfile << current_vector[i][j] << " " << endl;
//            }
//            myfile << "_________________________" << endl;
//        }
        float sum = 0;

        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx + 1; i++) {
                sum = sum + v[i][j];
            }
//            sum = (sum / Nx );
        }
        myfile <<"summ = "<< sum << endl;
        myfile.close();
    } else cout << "Unable to open file";

}
