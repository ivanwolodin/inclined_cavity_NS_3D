#pragma once

#include <vector>

using std::vector;
typedef double floatT;

int Nx = 30;
int Ny = 30;
int Nz = 1;
const int k = 1;

floatT Re = 2;

floatT dx = (floatT) 1 / Nx;
//float dy = (float) 1 / Ny;

vector< floatT> dy;
vector< floatT> dy1;
vector< floatT> dydy2;
vector< floatT> dtdy2;

floatT dz = 1; //(float) 1 / Nz;

floatT dt = 0.001;
floatT pressure_pseudo_time = 0.0001;

//float pressure_part =
//        (dx * dx * dy * dy * dz * dz) / (2 * (dy * dy * dz * dz + dx * dx * dz * dz + dx * dx * dy * dy));

floatT gravity = 0.00000;

floatT velocity_calculation_precision = 0.00001;
floatT pressure_calculation_precision = 0.000001;

vector<vector<floatT>  > u;
vector<vector<floatT>  > u_auxilliary;

vector<vector<floatT> > v;
vector<vector<floatT> > v_auxilliary;

vector<vector<floatT> >  p;
vector<vector<floatT> >  p_previous;

// for convergence check
vector<vector<floatT> > u_previous;
vector<vector<floatT> > v_previous;


const floatT mod_u = 1.0;
const floatT mod_v = 1.0;
const floatT mod_w = 1.0;

int debug = 1; // 1 stands for True, 0 - for False