#pragma once

#include <vector>

using std::vector;


int Nx = 60;
int Ny = 60;
//int Nz = 30;


float Re = 3000;

float dx = (float) 2 / Nx;
float dy = (float) 2 / Ny;
//float dz = (float) 1 / Nz;

float dt = 0.001;

float pressure_part =
        (dx * dx * dy * dy) / (2 * (dy * dy + dx * dx + dx * dx * dy * dy));

float gravity = 0.000001;

float velocity_calculation_precision = 0.0000000001;
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

float pressure_time = 0.001;
