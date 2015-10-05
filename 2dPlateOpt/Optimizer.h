#pragma once

#include <time.h>
#include <vector>
#include <iostream>
#include <fstream>
#include "Solver.h"
#include "plate_var_types.h"
#include "asa_user.h"

#define SOLVER_E1 102970000000.0
#define SOLVER_E2 7550000000.0
#define SOLVER_BY0 1.0 
#define SOLVER_AP 0.1524 * 2.0
#define SOLVER_BP 0.1524

#define STRESS_P0 30000.0
#define STRESS_TAU 0.01
#define STRESS_RAD_FACTOR 8.0

#define CHAR_TIME 0.03
#define J0_SCALE 1.e8

#define GRAD_SIZE 2

using std::vector;
using std::ofstream;
using std::cout;
using std::endl;

void optimizeASA( const vector<N_PRES>& params );

N_PRES calcValTaus( N_PRES* x, long n );
N_PRES calcValGradTaus( N_PRES* g, N_PRES* x, long n );

//HagerASA-compatible wrappers around funcs to compute 1st order optimization info
N_PRES calcValASAcurSin( asa_objective* asa );
void calcGradASAcurSin( asa_objective* asa );
N_PRES calcValGradASAcurSin( asa_objective* asa );