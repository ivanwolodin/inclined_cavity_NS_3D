#pragma once

#include <vector>

using std::vector;


int Nx = 30;
int Ny = 30;
int Nz = 1;
const int k = 1;

float Re = 2;

float dx = (float) 1 / Nx;
float dy = (float) 1 / Ny;
float dz = 1; //(float) 1 / Nz;

float dt = 0.001;
float pressure_pseudo_time = 0.0001;

//float pressure_part =
//        (dx * dx * dy * dy * dz * dz) / (2 * (dy * dy * dz * dz + dx * dx * dz * dz + dx * dx * dy * dy));

float gravity = 0.00000;

float velocity_calculation_precision = 0.00001;
float pressure_calculation_precision = 0.000001;

vector<vector<double>  > u;
vector<vector<double>  > u_auxilliary;

vector<vector<double> > v;
vector<vector<double> > v_auxilliary;

vector<vector<double> >  p;
vector<vector<double> >  p_previous;

// for convergence check
vector<vector<double> > u_previous;
vector<vector<double> > v_previous;


const float mod_u = 1.0;
const float mod_v = 1.0;
const float mod_w = 1.0;

int debug = 1; // 1 stands for True, 0 - for False